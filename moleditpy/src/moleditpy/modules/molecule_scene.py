import traceback

from PyQt6.QtWidgets import (
    QApplication, QGraphicsScene, QGraphicsItem,
    QGraphicsLineItem
)

from PyQt6.QtGui import (
    QPen, QCursor
)


from PyQt6.QtCore import (
    Qt, QPointF, QRectF, QLineF
)
import math

try:
    from .template_preview_item import TemplatePreviewItem
    from .atom_item import AtomItem
    from .bond_item import BondItem
except Exception:
    from modules.template_preview_item import TemplatePreviewItem
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem

try:
    from .constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE
except Exception:
    from modules.constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE

# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    from . import sip_isdeleted_safe
except Exception:
    from modules import sip_isdeleted_safe


class MoleculeScene(QGraphicsScene):
    def clear_template_preview(self):
        """テンプレートプレビュー用のゴースト線を全て消す"""
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(item, '_is_template_preview', False):
                try:
                    # If SIP reports the wrapper as deleted, skip it. Otherwise
                    # ensure it is still in a scene before attempting removal.
                    if sip_isdeleted_safe(item):
                        continue
                    sc = None
                    try:
                        sc = item.scene() if hasattr(item, 'scene') else None
                    except Exception:
                        sc = None
                    if sc is None:
                        continue
                    try:
                        self.removeItem(item)
                    except Exception:
                        # Best-effort: ignore removal errors to avoid crashes during teardown
                        pass
                except Exception:
                    # Non-fatal: continue with other items
                    continue
        self.template_context = {}
        if hasattr(self, 'template_preview'):
            self.template_preview.hide()

    def __init__(self, data, window):
        super().__init__()
        self.data, self.window = data, window
        self.mode, self.current_atom_symbol = 'select', 'C'
        self.bond_order, self.bond_stereo = 1, 0
        self.start_atom, self.temp_line, self.start_pos = None, None, None; self.press_pos = None
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.hovered_item = None
        
        self.key_to_symbol_map = {
            Qt.Key.Key_C: 'C', Qt.Key.Key_N: 'N', Qt.Key.Key_O: 'O', Qt.Key.Key_S: 'S',
            Qt.Key.Key_F: 'F', Qt.Key.Key_B: 'B', Qt.Key.Key_I: 'I', Qt.Key.Key_H: 'H',
            Qt.Key.Key_P: 'P',
        }
        self.key_to_symbol_map_shift = { Qt.Key.Key_C: 'Cl', Qt.Key.Key_B: 'Br', Qt.Key.Key_S: 'Si',}

        self.key_to_bond_mode_map = {
            Qt.Key.Key_1: 'bond_1_0',
            Qt.Key.Key_2: 'bond_2_0',
            Qt.Key.Key_3: 'bond_3_0',
            Qt.Key.Key_W: 'bond_1_1',
            Qt.Key.Key_D: 'bond_1_2',
        }
        self.reinitialize_items()

    def reinitialize_items(self):
        self.template_preview = TemplatePreviewItem(); self.addItem(self.template_preview)
        self.template_preview.hide(); self.template_preview_points = []; self.template_context = {}
        # Hold strong references to deleted wrappers for the lifetime of the scene
        # to avoid SIP/C++ finalization causing segfaults when Python still
        # briefly touches those objects elsewhere in the app. Items collected
        # here are hidden and never accessed again by normal code paths.
        self._deleted_items = []
        # Ensure we purge any held deleted-wrapper references when the
        # application is shutting down. Connecting here is safe even if
        # multiple scenes exist; the slot is defensive and idempotent.
        try:
            app = QApplication.instance()
            if app is not None:
                try:
                    app.aboutToQuit.connect(self.purge_deleted_items)
                except Exception:
                    # If connecting fails for any reason, continue without
                    # the connection — at worst holders will be freed by
                    # process teardown.
                    pass
        except Exception:
            pass

    def clear_all_problem_flags(self):
        """全ての AtomItem の has_problem フラグをリセットし、再描画する"""
        needs_update = False
        for atom_data in self.data.atoms.values():
            item = atom_data.get('item')
            # hasattr は安全性のためのチェック
            if item and hasattr(item, 'has_problem') and item.has_problem: 
                item.has_problem = False
                item.update()
                needs_update = True
        return needs_update

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        
        # 削除されたオブジェクトを安全にチェックして初期位置を記録
        self.initial_positions_in_event = {}
        for item in self.items():
            if isinstance(item, AtomItem):
                try:
                    self.initial_positions_in_event[item] = item.pos()
                except RuntimeError:
                    # オブジェクトが削除されている場合はスキップ
                    continue

        if not self.window.is_2d_editable:
            return

        if event.button() == Qt.MouseButton.RightButton:
            item = self.itemAt(event.scenePos(), self.views()[0].transform())
            if not isinstance(item, (AtomItem, BondItem)):
                return # 対象外のものをクリックした場合は何もしない
            data_changed = False
            # If the user has a rectangular multi-selection and the clicked item
            # is part of that selection, delete all selected items (atoms/bonds).
            try:
                selected_items = [it for it in self.selectedItems() if isinstance(it, (AtomItem, BondItem))]
            except Exception:
                selected_items = []

            if len(selected_items) > 1 and item in selected_items and not self.mode.startswith(('template', 'charge', 'radical')):
                # Delete the entire rectangular selection
                data_changed = self.delete_items(set(selected_items))
                if data_changed:
                    self.window.push_undo_state()
                self.press_pos = None
                event.accept()
                return
            # --- E/Zモード専用処理 ---
            if self.mode == 'bond_2_5':
                if isinstance(item, BondItem):
                    try:
                        # E/Zラベルを消す（ノーマルに戻す）
                        if item.stereo in [3, 4]:
                            item.set_stereo(0)
                            # データモデルも更新
                            for (id1, id2), bdata in self.data.bonds.items():
                                if bdata.get('item') is item:
                                    bdata['stereo'] = 0
                                    break
                            self.window.push_undo_state()
                            data_changed = False  # ここでundo済みなので以降で積まない
                    except Exception as e:
                        print(f"Error clearing E/Z label: {e}")
                        
                        traceback.print_exc()
                        if hasattr(self.window, 'statusBar'):
                            self.window.statusBar().showMessage(f"Error clearing E/Z label: {e}", 5000)
                # AtomItemは何もしない
            # --- 通常の処理 ---
            elif isinstance(item, AtomItem):
                # ラジカルモードの場合、ラジカルを0にする
                if self.mode == 'radical' and item.radical != 0:
                    item.prepareGeometryChange()
                    item.radical = 0
                    self.data.atoms[item.atom_id]['radical'] = 0
                    item.update_style()
                    data_changed = True
                # 電荷モードの場合、電荷を0にする
                elif self.mode in ['charge_plus', 'charge_minus'] and item.charge != 0:
                    item.prepareGeometryChange()
                    item.charge = 0
                    self.data.atoms[item.atom_id]['charge'] = 0
                    item.update_style()
                    data_changed = True
                # 上記以外のモード（テンプレート、電荷、ラジカルを除く）では原子を削除
                elif not self.mode.startswith(('template', 'charge', 'radical')):
                    data_changed = self.delete_items({item})
            elif isinstance(item, BondItem):
                # テンプレート、電荷、ラジカルモード以外で結合を削除
                if not self.mode.startswith(('template', 'charge', 'radical')):
                    data_changed = self.delete_items({item})

            if data_changed:
                self.window.push_undo_state()
            self.press_pos = None
            event.accept()
            return # 右クリック処理を完了し、左クリックの処理へ進ませない

        if self.mode.startswith('template'):
            self.clearSelection() # テンプレートモードでは選択処理を一切行わず、クリック位置の記録のみ行う
            return

        # Z,Eモードの時は選択処理を行わないようにする
        if self.mode in ['bond_2_5']:
            self.clearSelection()
            event.accept()
            return

        if getattr(self, "mode", "") != "select":
            self.clearSelection()
            event.accept()

        item = self.itemAt(self.press_pos, self.views()[0].transform())

        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != 'select':
                self.clearSelection()
                self.temp_line = QGraphicsLineItem(QLineF(self.start_atom.pos(), self.press_pos))
                self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
                self.addItem(self.temp_line)
            else:
                super().mousePressEvent(event)
        elif item is None and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos)); self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine)); self.addItem(self.temp_line)
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.window.is_2d_editable:
            return 

        if self.mode.startswith('template'):
            self.update_template_preview(event.scenePos())
        
        if not self.mouse_moved_since_press and self.press_pos:
            if (event.scenePos() - self.press_pos).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True
        
        if self.temp_line and not self.mode.startswith('template'):
            start_point = self.start_atom.pos() if self.start_atom else self.start_pos
            if not start_point:
                super().mouseMoveEvent(event)
                return

            current_pos = event.scenePos()
            end_point = current_pos

            target_atom = None
            for item in self.items(current_pos):
                if isinstance(item, AtomItem):
                    target_atom = item
                    break
            
            is_valid_snap_target = (
                target_atom is not None and
                (self.start_atom is None or target_atom is not self.start_atom)
            )

            if is_valid_snap_target:
                end_point = target_atom.pos()
            
            self.temp_line.setLine(QLineF(start_point, end_point))
        else: 
            # テンプレートモードであっても、ホバーイベントはここで伝播する
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if not self.window.is_2d_editable:
            return 

        end_pos = event.scenePos()
        is_click = self.press_pos and (end_pos - self.press_pos).manhattanLength() < QApplication.startDragDistance()

        if self.temp_line:
            try:
                if not sip_isdeleted_safe(self.temp_line):
                    try:
                        if getattr(self.temp_line, 'scene', None) and self.temp_line.scene():
                            self.removeItem(self.temp_line)
                    except Exception:
                        pass
            except Exception:
                try:
                    self.removeItem(self.temp_line)
                except Exception:
                    pass
            finally:
                self.temp_line = None

        if self.mode.startswith('template') and is_click:
            if self.template_context and self.template_context.get('points'):
                context = self.template_context
                # Check if this is a user template
                if self.mode.startswith('template_user'):
                    self.add_user_template_fragment(context)
                else:
                    self.add_molecule_fragment(context['points'], context['bonds_info'], existing_items=context.get('items', []))
                self.data_changed_in_event = True
                # イベント処理をここで完了させ、下のアイテムが選択されるのを防ぐ
                self.start_atom=None; self.start_pos = None; self.press_pos = None
                if self.data_changed_in_event: self.window.push_undo_state()
                return

        released_item = self.itemAt(end_pos, self.views()[0].transform())

        # 1. 特殊モード（ラジカル/電荷）の処理
        if (self.mode == 'radical') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            atom.prepareGeometryChange()
            # ラジカルの状態をトグル (0 -> 1 -> 2 -> 0)
            atom.radical = (atom.radical + 1) % 3 
            self.data.atoms[atom.atom_id]['radical'] = atom.radical
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom=None; self.start_pos = None; self.press_pos = None
            if self.data_changed_in_event: self.window.push_undo_state()
            return
        elif (self.mode == 'charge_plus' or self.mode == 'charge_minus') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            atom.prepareGeometryChange()
            delta = 1 if self.mode == 'charge_plus' else -1
            atom.charge += delta
            self.data.atoms[atom.atom_id]['charge'] = atom.charge
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom=None; self.start_pos = None; self.press_pos = None
            if self.data_changed_in_event: self.window.push_undo_state()
            return

        elif self.mode.startswith('bond') and is_click and isinstance(released_item, BondItem):
            b = released_item 
            if self.mode == 'bond_2_5':
                try:
                    if b.order == 2:
                        current_stereo = b.stereo
                        if current_stereo not in [3, 4]:
                            new_stereo = 3  # None -> Z
                        elif current_stereo == 3:
                            new_stereo = 4  # Z -> E
                        else:  # current_stereo == 4
                            new_stereo = 0  # E -> None
                        self.update_bond_stereo(b, new_stereo)
                        self.window.push_undo_state()  # ここでUndo stackに積む
                except Exception as e:
                    print(f"Error in E/Z stereo toggle: {e}")
                    
                    traceback.print_exc()
                    if hasattr(self.window, 'statusBar'):
                        self.window.statusBar().showMessage(f"Error changing E/Z stereochemistry: {e}", 5000)
                return # この後の処理は行わない
            elif self.bond_stereo != 0 and b.order == self.bond_order and b.stereo == self.bond_stereo:
                # 方向性を反転させる
                old_id1, old_id2 = b.atom1.atom_id, b.atom2.atom_id
                # 1. 古い方向の結合をデータから削除
                self.data.remove_bond(old_id1, old_id2)
                # 2. 逆方向で結合をデータに再追加
                new_key, _ = self.data.add_bond(old_id2, old_id1, self.bond_order, self.bond_stereo)
                # 3. BondItemの原子参照を入れ替え、新しいデータと関連付ける
                b.atom1, b.atom2 = b.atom2, b.atom1
                self.data.bonds[new_key]['item'] = b
                # 4. 見た目を更新
                b.update_position()
            else:
                # 既存の結合を一度削除
                self.data.remove_bond(b.atom1.atom_id, b.atom2.atom_id)
                # BondItemが記憶している方向(b.atom1 -> b.atom2)で、新しい結合様式を再作成
                # これにより、修正済みのadd_bondが呼ばれ、正しい方向で保存される
                new_key, _ = self.data.add_bond(b.atom1.atom_id, b.atom2.atom_id, self.bond_order, self.bond_stereo)
                # BondItemの見た目とデータ参照を更新
                b.prepareGeometryChange()
                b.order = self.bond_order
                b.stereo = self.bond_stereo
                self.data.bonds[new_key]['item'] = b
                b.update()
            self.clearSelection()
            self.data_changed_in_event = True
        # 3. 新規原子・結合の作成処理 (atom_* モード および すべての bond_* モードで許可)
        elif self.start_atom and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            line = QLineF(self.start_atom.pos(), end_pos); end_item = self.itemAt(end_pos, self.views()[0].transform())
            # 使用する結合様式を決定
            # atomモードの場合は bond_order/stereo を None にして create_bond にデフォルト値(1, 0)を適用
            # bond_* モードの場合は現在の設定 (self.bond_order/stereo) を使用
            order_to_use = self.bond_order if self.mode.startswith('bond') else None
            stereo_to_use = self.bond_stereo if self.mode.startswith('bond') else None
            if is_click:
                # 短いクリック: 既存原子のシンボル更新 (atomモードのみ)
                if self.mode.startswith('atom') and self.start_atom.symbol != self.current_atom_symbol:
                    self.start_atom.symbol=self.current_atom_symbol; self.data.atoms[self.start_atom.atom_id]['symbol']=self.current_atom_symbol; self.start_atom.update_style()
                    self.data_changed_in_event = True
            else:
                # ドラッグ: 新規結合または既存原子への結合
                if isinstance(end_item, AtomItem) and self.start_atom!=end_item: 
                    self.create_bond(self.start_atom, end_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                else:
                    new_id = self.create_atom(self.current_atom_symbol, end_pos); new_item = self.data.atoms[new_id]['item']
                    self.create_bond(self.start_atom, new_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                self.data_changed_in_event = True
        # 4. 空白領域からの新規作成処理 (atom_* モード および すべての bond_* モードで許可)
        elif self.start_pos and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            line = QLineF(self.start_pos, end_pos)
            # 使用する結合様式を決定
            order_to_use = self.bond_order if self.mode.startswith('bond') else None
            stereo_to_use = self.bond_stereo if self.mode.startswith('bond') else None
            if line.length() < 10:
                self.create_atom(self.current_atom_symbol, end_pos); self.data_changed_in_event = True
            else:
                end_item = self.itemAt(end_pos, self.views()[0].transform())
                if isinstance(end_item, AtomItem):
                    start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                    start_item = self.data.atoms[start_id]['item']
                    self.create_bond(start_item, end_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                else:
                    start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                    end_id = self.create_atom(self.current_atom_symbol, end_pos)
                    self.create_bond(
                        self.data.atoms[start_id]['item'], 
                        self.data.atoms[end_id]['item'], 
                        bond_order=order_to_use, 
                        bond_stereo=stereo_to_use
                    )
                self.data_changed_in_event = True 
        # 5. それ以外の処理 (Selectモードなど)
        else: super().mouseReleaseEvent(event)

        # 削除されたオブジェクトを安全にチェック
        moved_atoms = []
        for item, old_pos in self.initial_positions_in_event.items():
            try:
                # オブジェクトが有効で、シーンに存在し、位置が変更されているかチェック
                if item.scene() and item.pos() != old_pos:
                    moved_atoms.append(item)
            except RuntimeError:
                # オブジェクトが削除されている場合はスキップ
                continue
        if moved_atoms:
            self.data_changed_in_event = True
            bonds_to_update = set()
            for atom in moved_atoms:
                try:
                    self.data.atoms[atom.atom_id]['pos'] = atom.pos()
                    bonds_to_update.update(atom.bonds)
                except RuntimeError:
                    # オブジェクトが削除されている場合はスキップ
                    continue
            for bond in bonds_to_update: bond.update_position()
            # 原子移動後に測定ラベルの位置を更新
            self.window.update_2d_measurement_labels()
            if self.views(): self.views()[0].viewport().update()
        self.start_atom=None; self.start_pos = None; self.press_pos = None; self.temp_line = None
        self.template_context = {}
        # Clear user template data when switching modes
        if hasattr(self, 'user_template_data'):
            self.user_template_data = None
        if self.data_changed_in_event: self.window.push_undo_state()

    def mouseDoubleClickEvent(self, event):
        """ダブルクリックイベントを処理する"""
        item = self.itemAt(event.scenePos(), self.views()[0].transform())

        if self.mode in ['charge_plus', 'charge_minus', 'radical'] and isinstance(item, AtomItem):
            if self.mode == 'radical':
                item.prepareGeometryChange()
                item.radical = (item.radical + 1) % 3
                self.data.atoms[item.atom_id]['radical'] = item.radical
                item.update_style()
            else:
                item.prepareGeometryChange()
                delta = 1 if self.mode == 'charge_plus' else -1
                item.charge += delta
                self.data.atoms[item.atom_id]['charge'] = item.charge
                item.update_style()

            self.window.push_undo_state()

            event.accept()
            return
        
        # Select-mode: double-click should select the clicked atom/bond and
        # only the atoms/bonds connected to it (the connected component).
        if self.mode == 'select' and isinstance(item, (AtomItem, BondItem)):
            try:
                start_atoms = set()
                if isinstance(item, AtomItem):
                    start_atoms.add(item)
                else:
                    # BondItem: start from both ends if available
                    a1 = getattr(item, 'atom1', None)
                    a2 = getattr(item, 'atom2', None)
                    if a1 is not None:
                        start_atoms.add(a1)
                    if a2 is not None:
                        start_atoms.add(a2)

                # BFS/DFS over atoms via bond references (defensive checks)
                atoms_to_visit = list(start_atoms)
                connected_atoms = set()
                connected_bonds = set()

                while atoms_to_visit:
                    a = atoms_to_visit.pop()
                    if a is None:
                        continue
                    if a in connected_atoms:
                        continue
                    connected_atoms.add(a)
                    # iterate bonds attached to atom
                    for b in getattr(a, 'bonds', []) or []:
                        if b is None:
                            continue
                        connected_bonds.add(b)
                        # find the other atom at the bond
                        other = None
                        try:
                            if getattr(b, 'atom1', None) is a:
                                other = getattr(b, 'atom2', None)
                            else:
                                other = getattr(b, 'atom1', None)
                        except Exception:
                            other = None
                        if other is not None and other not in connected_atoms:
                            atoms_to_visit.append(other)

                # Apply selection: clear previous and select only these
                try:
                    self.clearSelection()
                except Exception:
                    pass

                for a in connected_atoms:
                    try:
                        a.setSelected(True)
                    except Exception:
                        try:
                            # fallback: set selected attribute if exists
                            setattr(a, 'selected', True)
                        except Exception:
                            pass
                for b in connected_bonds:
                    try:
                        b.setSelected(True)
                    except Exception:
                        try:
                            setattr(b, 'selected', True)
                        except Exception:
                            pass

                event.accept()
                return
            except Exception:
                # On any unexpected error, fall back to default handling
                pass

        elif self.mode in ['bond_2_5']:
                event.accept()
                return

        super().mouseDoubleClickEvent(event)

    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self.data.add_atom(symbol, pos, charge=charge, radical=radical)
        atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
        self.data.atoms[atom_id]['item'] = atom_item; self.addItem(atom_item); return atom_id


    def create_bond(self, start_atom, end_atom, bond_order=None, bond_stereo=None):
        try:
            if start_atom is None or end_atom is None:
                print("Error: Cannot create bond with None atoms")
                return
                
            exist_b = self.find_bond_between(start_atom, end_atom)
            if exist_b:
                return

            # 引数で次数が指定されていればそれを使用し、なければ現在のモードの値を使用する
            order_to_use = self.bond_order if bond_order is None else bond_order
            stereo_to_use = self.bond_stereo if bond_stereo is None else bond_stereo

            key, status = self.data.add_bond(start_atom.atom_id, end_atom.atom_id, order_to_use, stereo_to_use)
            if status == 'created':
                bond_item = BondItem(start_atom, end_atom, order_to_use, stereo_to_use)
                self.data.bonds[key]['item'] = bond_item
                if hasattr(start_atom, 'bonds'):
                    start_atom.bonds.append(bond_item)
                if hasattr(end_atom, 'bonds'):
                    end_atom.bonds.append(bond_item)
                self.addItem(bond_item)
            
            if hasattr(start_atom, 'update_style'):
                start_atom.update_style()
            if hasattr(end_atom, 'update_style'):
                end_atom.update_style()
                
        except Exception as e:
            print(f"Error creating bond: {e}")
            
            traceback.print_exc()

    def add_molecule_fragment(self, points, bonds_info, existing_items=None, symbol='C'):
        """
        add_molecule_fragment の最終確定版。
        - 既存の結合次数を変更しないポリシーを徹底（最重要）。
        - ベンゼン環テンプレートは、フューズされる既存結合の次数に基づき、
          「新規に作られる二重結合が2本になるように」回転を決定するロジックを適用（条件分岐あり）。
        """
    
        num_points = len(points)
        atom_items = [None] * num_points

        is_benzene_template = (num_points == 6 and any(o == 2 for _, _, o in bonds_info))

    
        def coords(p):
            if hasattr(p, 'x') and hasattr(p, 'y'):
                return (p.x(), p.y())
            try:
                return (p[0], p[1])
            except Exception:
                raise ValueError("point has no x/y")
    
        def dist_pts(a, b):
            ax, ay = coords(a); bx, by = coords(b)
            return math.hypot(ax - bx, ay - by)
    
        # --- 1) 既にクリックされた existing_items をテンプレート頂点にマップ ---
        existing_items = existing_items or []
        used_indices = set()
        ref_lengths = [dist_pts(points[i], points[j]) for i, j, _ in bonds_info if i < num_points and j < num_points]
        avg_len = (sum(ref_lengths) / len(ref_lengths)) if ref_lengths else 20.0
        map_threshold = max(0.5 * avg_len, 8.0)
    
        for ex_item in existing_items:
            try:
                ex_pos = ex_item.pos()
                best_idx, best_d = -1, float('inf')
                for i, p in enumerate(points):
                    if i in used_indices: continue
                    d = dist_pts(p, ex_pos)
                    if best_d is None or d < best_d:
                        best_d, best_idx = d, i
                if best_idx != -1 and best_d <= max(map_threshold, 1.5 * avg_len):
                    atom_items[best_idx] = ex_item
                    used_indices.add(best_idx)
            except Exception:
                pass
    
        # --- 2) シーン内既存原子を self.data.atoms から列挙してマップ ---
        mapped_atoms = {it for it in atom_items if it is not None}
        for i, p in enumerate(points):
            if atom_items[i] is not None: continue
            
            nearby = None
            best_d = float('inf')
            
            for atom_data in self.data.atoms.values():
                a_item = atom_data.get('item')
                if not a_item or a_item in mapped_atoms: continue
                try:
                    d = dist_pts(p, a_item.pos())
                except Exception:
                    continue
                if d < best_d:
                    best_d, nearby = d, a_item

            if nearby and best_d <= map_threshold:
                atom_items[i] = nearby
                mapped_atoms.add(nearby)
    
        # --- 3) 足りない頂点は新規作成　---
        for i, p in enumerate(points):
            if atom_items[i] is None:
                atom_id = self.create_atom(symbol, p)
                atom_items[i] = self.data.atoms[atom_id]['item']
    
        # --- 4) テンプレートのボンド配列を決定（ベンゼン回転合わせの処理） ---
        template_bonds_to_use = list(bonds_info)
        is_6ring = (num_points == 6 and len(bonds_info) == 6)
        template_has_double = any(o == 2 for (_, _, o) in bonds_info)
    
        if is_6ring and template_has_double:
            existing_orders = {} # key: bonds_infoのインデックス, value: 既存の結合次数
            for k, (i_idx, j_idx, _) in enumerate(bonds_info):
                if i_idx < len(atom_items) and j_idx < len(atom_items):
                    a, b = atom_items[i_idx], atom_items[j_idx]
                    if a is None or b is None: continue
                    eb = self.find_bond_between(a, b)
                    if eb:
                        existing_orders[k] = getattr(eb, 'order', 1) 

            if existing_orders:
                orig_orders = [o for (_, _, o) in bonds_info]
                best_rot = 0
                max_score = -999 # スコアは「適合度」を意味する

                # --- フューズされた辺の数による条件分岐 ---
                if len(existing_orders) >= 2:
                    # 2辺以上フューズ: 単純に既存の辺の次数とテンプレートの辺の次数が一致するものを最優先する
                    # (この場合、新しい環を交互配置にするのは難しく、単に既存の構造を壊さないことを優先)
                    for rot in range(num_points):
                        current_score = sum(100 for k, exist_order in existing_orders.items() 
                                            if orig_orders[(k + rot) % num_points] == exist_order)
                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot

                elif len(existing_orders) == 1:
                    # 1辺フューズ: 既存の辺を維持しつつ、その両隣で「反転一致」を達成し、新しい環を交互配置にする
                    
                    # フューズされた辺のインデックスと次数を取得
                    k_fuse = next(iter(existing_orders.keys()))
                    exist_order = existing_orders[k_fuse]
                    
                    # 目標: フューズされた辺の両隣（k-1とk+1）に来るテンプレートの次数が、既存の辺の次数と逆であること
                    # k_adj_1 -> (k_fuse - 1) % 6
                    # k_adj_2 -> (k_fuse + 1) % 6
                    
                    for rot in range(num_points):
                        current_score = 0
                        rotated_template_order = orig_orders[(k_fuse + rot) % num_points]

                        # 1. まず、フューズされた辺自体が次数を反転させられる位置にあるかチェック（必須ではないが、回転を絞る）
                        if (exist_order == 1 and rotated_template_order == 2) or \
                           (exist_order == 2 and rotated_template_order == 1):
                            current_score += 100 # 大幅ボーナス: 理想的な回転

                        # 2. 次に、両隣の辺の次数をチェック（交互配置維持の主目的）
                        # 既存辺の両隣は、新規に作成されるため、テンプレートの次数でボンドが作成されます。
                        # ここで、テンプレートの次数が既存辺の次数と逆になる回転を選ぶ必要があります。
                        
                        # テンプレートの辺は、回転後のk_fuseの両隣（m_adj1, m_adj2）
                        m_adj1 = (k_fuse - 1 + rot) % num_points 
                        m_adj2 = (k_fuse + 1 + rot) % num_points
                        
                        neighbor_order_1 = orig_orders[m_adj1]
                        neighbor_order_2 = orig_orders[m_adj2]

                        # 既存が単結合(1)の場合、両隣は二重結合(2)であってほしい
                        if exist_order == 1:
                            if neighbor_order_1 == 2: current_score += 50
                            if neighbor_order_2 == 2: current_score += 50
                        
                        # 既存が二重結合(2)の場合、両隣は単結合(1)であってほしい
                        elif exist_order == 2:
                            if neighbor_order_1 == 1: current_score += 50
                            if neighbor_order_2 == 1: current_score += 50
                            
                        # 3. タイブレーク: その他の既存結合（フューズ辺ではない）との次数一致度も加味
                        for k, e_order in existing_orders.items():
                             if k != k_fuse:
                                r_t_order = orig_orders[(k + rot) % num_points]
                                if r_t_order == e_order: current_score += 10 # 既存構造維持のボーナス
                        
                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot
                
                # 最終的な回転を反映
                new_tb = []
                for m in range(num_points):
                    i_idx, j_idx, _ = bonds_info[m]
                    new_order = orig_orders[(m + best_rot) % num_points]
                    new_tb.append((i_idx, j_idx, new_order))
                template_bonds_to_use = new_tb
    
        # --- 5) ボンド作成／更新---
        for id1_idx, id2_idx, order in template_bonds_to_use:
            if id1_idx < len(atom_items) and id2_idx < len(atom_items):
                a_item, b_item = atom_items[id1_idx], atom_items[id2_idx]
                if not a_item or not b_item or a_item is b_item: continue

                id1, id2 = a_item.atom_id, b_item.atom_id
                if id1 > id2: id1, id2 = id2, id1

                exist_b = self.find_bond_between(a_item, b_item)

                if exist_b:
                    # デフォルトでは既存の結合を維持する
                    should_overwrite = False

                    # 条件1: ベンゼン環テンプレートであること
                    # 条件2: 接続先が単結合であること
                    if is_benzene_template and exist_b.order == 1:

                        # 条件3: 接続先の単結合が共役系の一部ではないこと
                        # (つまり、両端の原子が他に二重結合を持たないこと)
                        atom1 = exist_b.atom1
                        atom2 = exist_b.atom2

                        # atom1が他に二重結合を持つかチェック
                        atom1_has_other_double_bond = any(b.order == 2 for b in atom1.bonds if b is not exist_b)

                        # atom2が他に二重結合を持つかチェック
                        atom2_has_other_double_bond = any(b.order == 2 for b in atom2.bonds if b is not exist_b)

                        # 両方の原子が他に二重結合を持たない「孤立した単結合」の場合のみ上書きフラグを立てる
                        if not atom1_has_other_double_bond and not atom2_has_other_double_bond:
                            should_overwrite = True

                    if should_overwrite:
                        # 上書き条件が全て満たされた場合にのみ、結合次数を更新
                        exist_b.order = order
                        exist_b.stereo = 0
                        self.data.bonds[(id1, id2)]['order'] = order
                        self.data.bonds[(id1, id2)]['stereo'] = 0
                        exist_b.update()
                    else:
                        # 上書き条件を満たさない場合は、既存の結合を維持する
                        continue
                else:
                    # 新規ボンド作成
                    self.create_bond(a_item, b_item, bond_order=order, bond_stereo=0)
        
        # --- 6) 表示更新　---
        for at in atom_items:
            try:
                if at: at.update_style() 
            except Exception:
                pass
    
        return atom_items


    def update_template_preview(self, pos):
        mode_parts = self.mode.split('_')
        
        # Check if this is a user template
        if len(mode_parts) >= 3 and mode_parts[1] == 'user':
            self.update_user_template_preview(pos)
            return
        
        is_aromatic = False
        if mode_parts[1] == 'benzene':
            n = 6
            is_aromatic = True
        else:
            try: n = int(mode_parts[1])
            except ValueError: return

        items_under = self.items(pos)  # top-most first
        item = None
        for it in items_under:
            if isinstance(it, (AtomItem, BondItem)):
                item = it
                break

        points, bonds_info = [], []
        l = DEFAULT_BOND_LENGTH
        self.template_context = {}


        if isinstance(item, AtomItem):
            p0 = item.pos()
            continuous_angle = math.atan2(pos.y() - p0.y(), pos.x() - p0.x())
            snap_angle_rad = math.radians(15)
            snapped_angle = round(continuous_angle / snap_angle_rad) * snap_angle_rad
            p1 = p0 + QPointF(l * math.cos(snapped_angle), l * math.sin(snapped_angle))
            points = self._calculate_polygon_from_edge(p0, p1, n)
            self.template_context['items'] = [item]

        elif isinstance(item, BondItem):
            # 結合にスナップ
            p0, p1 = item.atom1.pos(), item.atom2.pos()
            points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=pos, use_existing_length=True)
            self.template_context['items'] = [item.atom1, item.atom2]

        else:
            angle_step = 2 * math.pi / n
            start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2
            points = [
                pos + QPointF(l * math.cos(start_angle + i * angle_step), l * math.sin(start_angle + i * angle_step))
                for i in range(n)
            ]

        if points:
            if is_aromatic:
                bonds_info = [(i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)]
            else:
                bonds_info = [(i, (i + 1) % n, 1) for i in range(n)]

            self.template_context['points'] = points
            self.template_context['bonds_info'] = bonds_info

            self.template_preview.set_geometry(points, is_aromatic)

            self.template_preview.show()
            if self.views():
                self.views()[0].viewport().update()
        else:
            self.template_preview.hide()
            if self.views():
                self.views()[0].viewport().update()

    def _calculate_polygon_from_edge(self, p0, p1, n, cursor_pos=None, use_existing_length=False):
        if n < 3: return []
        v_edge = p1 - p0
        edge_length = math.sqrt(v_edge.x()**2 + v_edge.y()**2)
        if edge_length == 0: return []
        
        target_length = edge_length if use_existing_length else DEFAULT_BOND_LENGTH
        
        v_edge = (v_edge / edge_length) * target_length
        
        if not use_existing_length:
             p1 = p0 + v_edge

        points = [p0, p1]
        
        interior_angle = (n - 2) * math.pi / n
        rotation_angle = math.pi - interior_angle
        
        if cursor_pos:
            # Note: v_edgeは正規化済みだが、方向は同じなので判定には問題ない
            v_cursor = cursor_pos - p0
            cross_product_z = (p1 - p0).x() * v_cursor.y() - (p1 - p0).y() * v_cursor.x()
            if cross_product_z < 0:
                rotation_angle = -rotation_angle

        cos_a, sin_a = math.cos(rotation_angle), math.sin(rotation_angle)
        
        current_p, current_v = p1, v_edge
        for _ in range(n - 2):
            new_vx = current_v.x() * cos_a - current_v.y() * sin_a
            new_vy = current_v.x() * sin_a + current_v.y() * cos_a
            current_v = QPointF(new_vx, new_vy)
            current_p = current_p + current_v
            points.append(current_p)
        return points

    def delete_items(self, items_to_delete):
        """指定されたアイテムセット（原子・結合）を安全な順序で削除する修正版"""
        # Hardened deletion: perform data-model removals first, then scene removals,
        # and always defensively check attributes to avoid accessing partially-deleted objects.
        if not items_to_delete:
            return False

        # First sanitize the incoming collection: only keep live, expected QGraphics wrappers
        try:
            sanitized = set()
            for it in items_to_delete:
                try:
                    if it is None:
                        continue
                    # Skip SIP-deleted wrappers early to avoid native crashes
                    if sip_isdeleted_safe(it):
                        continue
                    # Only accept AtomItem/BondItem or other QGraphicsItem subclasses
                    if isinstance(it, (AtomItem, BondItem, QGraphicsItem)):
                        sanitized.add(it)
                except Exception:
                    # If isinstance or sip check raises, skip this entry
                    continue
            items_to_delete = sanitized
        except Exception:
            # If sanitization fails, fall back to original input and proceed defensively
            pass

        try:
            atoms_to_delete = {item for item in items_to_delete if isinstance(item, AtomItem)}
            bonds_to_delete = {item for item in items_to_delete if isinstance(item, BondItem)}

            # Include bonds attached to atoms being deleted
            for atom in list(atoms_to_delete):
                try:
                    if hasattr(atom, 'bonds') and atom.bonds:
                        for b in list(atom.bonds):
                            bonds_to_delete.add(b)
                except Exception:
                    # If accessing bonds raises (item partially deleted), skip
                    continue

            # Determine atoms that will remain but whose bond lists must be updated
            atoms_to_update = set()
            for bond in list(bonds_to_delete):
                try:
                    a1 = getattr(bond, 'atom1', None)
                    a2 = getattr(bond, 'atom2', None)
                    if a1 and a1 not in atoms_to_delete:
                        atoms_to_update.add(a1)
                    if a2 and a2 not in atoms_to_delete:
                        atoms_to_update.add(a2)
                except Exception:
                    continue

            # 1) Update surviving atoms' bond lists to remove references to bonds_to_delete
            #    (Important: remove BondItem references so atoms properly reflect
            #     that they have no remaining bonds and update visibility accordingly.)
            for atom in list(atoms_to_update):
                try:
                    if sip_isdeleted_safe(atom):
                        continue
                    # Defensive: if the atom has a bonds list, filter out bonds being deleted
                    if hasattr(atom, 'bonds') and atom.bonds:
                        try:
                            # Replace in-place to preserve any other references.
                            # Avoid touching SIP-deleted bond wrappers: build a set
                            # of live bonds-to-delete and also prune any SIP-deleted
                            # entries that may exist in atom.bonds.
                            live_btd = {b for b in bonds_to_delete if not sip_isdeleted_safe(b)}

                            # First, remove any SIP-deleted bond wrappers from atom.bonds
                            atom.bonds[:] = [b for b in atom.bonds if not sip_isdeleted_safe(b)]

                            # Then remove bonds which are in the live_btd set
                            if live_btd:
                                atom.bonds[:] = [b for b in atom.bonds if b not in live_btd]
                        except Exception:
                            # Fall back to iterative removal if list comprehension fails
                            try:
                                live_btd = [b for b in list(bonds_to_delete) if not sip_isdeleted_safe(b)]
                                for b in live_btd:
                                    if b in atom.bonds:
                                        atom.bonds.remove(b)
                            except Exception:
                                pass

                    # After pruning bond references, update visual style so carbons without
                    # bonds become visible again.
                    if hasattr(atom, 'update_style'):
                        atom.update_style()
                except Exception:
                    continue

            # 2) Remove bonds/atoms from the data model first (so other code reading the model
            #    doesn't encounter stale entries while we are removing graphics)
            for bond in list(bonds_to_delete):
                try:
                    a1 = getattr(bond, 'atom1', None)
                    a2 = getattr(bond, 'atom2', None)
                    if a1 and a2 and hasattr(self, 'data'):
                        try:
                            self.data.remove_bond(a1.atom_id, a2.atom_id)
                        except Exception:
                            # try reverse order if remove_bond expects ordered tuple
                            try:
                                self.data.remove_bond(a2.atom_id, a1.atom_id)
                            except Exception:
                                pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    if hasattr(atom, 'atom_id') and hasattr(self, 'data'):
                        try:
                            self.data.remove_atom(atom.atom_id)
                        except Exception:
                            pass
                except Exception:
                    continue

            # Invalidate any pending implicit-hydrogen UI updates because the
            # underlying data model changed. This prevents a scheduled
            # update_implicit_hydrogens closure from touching atoms/bonds that
            # were just removed. Do a single increment rather than one per-atom.
            try:
                self._ih_update_counter += 1
            except Exception:
                try:
                    self._ih_update_counter = 0
                except Exception:
                    pass

            # 3) Remove graphic items from the scene (bonds first)
            # To avoid calling into methods on wrappers that may refer to
            # already-deleted C++ objects (which can cause a native crash when
            # SIP is not available), take a snapshot of the current scene's
            # items and use membership tests instead of calling item.scene().
            try:
                current_scene_items = set(self.items())
            except Exception:
                # If for any reason items() fails, fall back to an empty set
                current_scene_items = set()

            for bond in list(bonds_to_delete):
                try:
                    # If the SIP wrapper is already deleted, skip it.
                    if sip_isdeleted_safe(bond):
                        continue
                    # Only attempt to remove the bond if it is present in the
                    # scene snapshot. This avoids calling bond.scene() which
                    # may invoke C++ on a deleted object.
                    if bond in current_scene_items:
                        try:
                            self.removeItem(bond)
                        except Exception:
                            pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    # Skip if wrapper is reported deleted by SIP
                    if sip_isdeleted_safe(atom):
                        continue
                    if atom in current_scene_items:
                        try:
                            self.removeItem(atom)
                        except Exception:
                            pass
                except Exception:
                    continue

            # 4) Instead of aggressively nullling object attributes (which can
            #    lead to C++/SIP finalization races and segfaults), keep a
            #    strong reference to the deleted wrappers for the lifetime of
            #    the scene. This prevents their underlying SIP wrappers from
            #    being finalized while other code may still touch them.
            try:
                if not hasattr(self, '_deleted_items') or self._deleted_items is None:
                    self._deleted_items = []
            except Exception:
                self._deleted_items = []

            for bond in list(bonds_to_delete):
                try:
                    # Hide the graphics item if possible and stash it
                    if not sip_isdeleted_safe(bond):
                        try:
                            bond.hide()
                        except Exception:
                            pass
                        try:
                            self._deleted_items.append(bond)
                        except Exception:
                            # Swallow any error while stashing
                            pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    if not sip_isdeleted_safe(atom):
                        try:
                            atom.hide()
                        except Exception:
                            pass
                        try:
                            self._deleted_items.append(atom)
                        except Exception:
                            pass
                except Exception:
                    continue

            # 5) Final visual updates for surviving atoms
            for atom in list(atoms_to_update):
                try:
                    if hasattr(atom, 'update_style'):
                        atom.update_style()
                except Exception:
                    continue

            return True

        except Exception as e:
            # Keep the application alive on unexpected errors
            print(f"Error during delete_items operation: {e}")
            
            traceback.print_exc()
            return False
    def purge_deleted_items(self):
        """Purge and release any held deleted-wrapper references.

        This is intended to be invoked on application shutdown to allow
        the process to release references to SIP/C++ wrappers that were
        kept around to avoid finalization races during normal runtime.
        The method is defensive: it tolerates partially-deleted wrappers
        and any SIP unavailability.
        """
        try:
            if not hasattr(self, '_deleted_items') or not self._deleted_items:
                return

            # Iterate a copy since we will clear the list.
            for obj in list(self._deleted_items):
                try:
                    # If the wrapper is still alive, attempt to hide it so
                    # the graphics subsystem isn't holding on to resources.
                    if not sip_isdeleted_safe(obj):
                        try:
                            obj.hide()
                        except Exception:
                            pass

                    # Try to clear container attributes that may hold refs
                    # to other scene objects (bonds, etc.) to help GC.
                    try:
                        if hasattr(obj, 'bonds') and getattr(obj, 'bonds') is not None:
                            try:
                                obj.bonds.clear()
                            except Exception:
                                # Try assignment fallback
                                try:
                                    obj.bonds = []
                                except Exception:
                                    pass
                    except Exception:
                        pass

                except Exception:
                    # Continue purging remaining items even if one fails.
                    continue

            # Finally, drop our references.
            try:
                self._deleted_items.clear()
            except Exception:
                try:
                    self._deleted_items = []
                except Exception:
                    pass

        except Exception as e:
            # Never raise during shutdown
            try:
                print(f"Error purging deleted items: {e}")
            except Exception:
                pass
    
    def add_user_template_fragment(self, context):
        """ユーザーテンプレートフラグメントを配置"""
        points = context.get('points', [])
        bonds_info = context.get('bonds_info', [])
        atoms_data = context.get('atoms_data', [])
        attachment_atom = context.get('attachment_atom')
        
        if not points or not atoms_data:
            return
        
        # Create atoms
        atom_id_map = {}  # template id -> scene atom id
        
        for i, (pos, atom_data) in enumerate(zip(points, atoms_data)):
            # Skip first atom if attaching to existing atom
            if i == 0 and attachment_atom:
                atom_id_map[atom_data['id']] = attachment_atom.atom_id
                continue
            
            symbol = atom_data.get('symbol', 'C')
            charge = atom_data.get('charge', 0)
            radical = atom_data.get('radical', 0)
            
            atom_id = self.data.add_atom(symbol, pos, charge, radical)
            atom_id_map[atom_data['id']] = atom_id
            
            # Create visual atom item
            atom_item = AtomItem(atom_id, symbol, pos, charge, radical)
            self.data.atoms[atom_id]['item'] = atom_item
            self.addItem(atom_item)
        
        # Create bonds (bonds_infoは必ずidベースで扱う)
        # まずindex→id変換テーブルを作る
        index_to_id = [atom_data.get('id', i) for i, atom_data in enumerate(atoms_data)]
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                # bonds_infoの0,1番目がindexならidに変換
                atom1_idx = bond_info[0]
                atom2_idx = bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                stereo = bond_info[3] if len(bond_info) > 3 else 0

                # index→id変換（すでにidならそのまま）
                if isinstance(atom1_idx, int) and atom1_idx < len(index_to_id):
                    template_atom1_id = index_to_id[atom1_idx]
                else:
                    template_atom1_id = atom1_idx
                if isinstance(atom2_idx, int) and atom2_idx < len(index_to_id):
                    template_atom2_id = index_to_id[atom2_idx]
                else:
                    template_atom2_id = atom2_idx

                atom1_id = atom_id_map.get(template_atom1_id)
                atom2_id = atom_id_map.get(template_atom2_id)

                if atom1_id is not None and atom2_id is not None:
                    # Skip if bond already exists
                    existing_bond = None
                    if (atom1_id, atom2_id) in self.data.bonds:
                        existing_bond = (atom1_id, atom2_id)
                    elif (atom2_id, atom1_id) in self.data.bonds:
                        existing_bond = (atom2_id, atom1_id)

                    if not existing_bond:
                        bond_key, _ = self.data.add_bond(atom1_id, atom2_id, order, stereo)
                        # Create visual bond item
                        atom1_item = self.data.atoms[atom1_id]['item']
                        atom2_item = self.data.atoms[atom2_id]['item']
                        if atom1_item and atom2_item:
                            bond_item = BondItem(atom1_item, atom2_item, order, stereo)
                            self.data.bonds[bond_key]['item'] = bond_item
                            self.addItem(bond_item)
                            atom1_item.bonds.append(bond_item)
                            atom2_item.bonds.append(bond_item)
        
        # Update atom visuals
        for atom_id in atom_id_map.values():
            if atom_id in self.data.atoms and self.data.atoms[atom_id]['item']:
                self.data.atoms[atom_id]['item'].update_style()
    
    def update_user_template_preview(self, pos):
        """ユーザーテンプレートのプレビューを更新"""
        # Robust user template preview: do not access self.data.atoms for preview-only atoms
        if not hasattr(self, 'user_template_data') or not self.user_template_data:
            return

        template_data = self.user_template_data
        atoms = template_data.get('atoms', [])
        bonds = template_data.get('bonds', [])

        if not atoms:
            return

        # Find attachment point (first atom or clicked item)
        items_under = self.items(pos)
        attachment_atom = None
        for item in items_under:
            if isinstance(item, AtomItem):
                attachment_atom = item
                break

        # Calculate template positions
        points = []
        # Find template bounds for centering
        if atoms:
            min_x = min(atom['x'] for atom in atoms)
            max_x = max(atom['x'] for atom in atoms)
            min_y = min(atom['y'] for atom in atoms)
            max_y = max(atom['y'] for atom in atoms)
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
        # Position template
        if attachment_atom:
            # Attach to existing atom
            attach_pos = attachment_atom.pos()
            offset_x = attach_pos.x() - atoms[0]['x']
            offset_y = attach_pos.y() - atoms[0]['y']
        else:
            # Center at cursor position
            offset_x = pos.x() - center_x
            offset_y = pos.y() - center_y
        # Calculate atom positions
        for atom in atoms:
            new_pos = QPointF(atom['x'] + offset_x, atom['y'] + offset_y)
            points.append(new_pos)
        # Create atom ID to index mapping (for preview only)
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get('id', i)
            atom_id_to_index[atom_id] = i
        # bonds_info をテンプレートの bonds から生成
        bonds_info = []
        for bond in bonds:
            atom1_idx = atom_id_to_index.get(bond['atom1'])
            atom2_idx = atom_id_to_index.get(bond['atom2'])
            if atom1_idx is not None and atom2_idx is not None:
                order = bond.get('order', 1)
                stereo = bond.get('stereo', 0)
                bonds_info.append((atom1_idx, atom2_idx, order, stereo))
        # プレビュー用: points, bonds_info から線を描画
        # 設置用 context を保存
        self.template_context = {
            'points': points,
            'bonds_info': bonds_info,
            'atoms_data': atoms,
            'attachment_atom': attachment_atom,
        }
        # 既存のプレビューアイテムを一旦クリア
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(item, '_is_template_preview', False):
                self.removeItem(item)

        # Draw preview lines only using calculated points (do not access self.data.atoms)
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                i, j = bond_info[0], bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                # stereo = bond_info[3] if len(bond_info) > 3 else 0
                if i < len(points) and j < len(points):
                    line = QGraphicsLineItem(QLineF(points[i], points[j]))
                    pen = QPen(Qt.black, 2 if order == 2 else 1)
                    line.setPen(pen)
                    line._is_template_preview = True  # フラグで区別
                    self.addItem(line)
        # Never access self.data.atoms here for preview-only atoms

    def leaveEvent(self, event):
        self.template_preview.hide(); super().leaveEvent(event)

    def set_hovered_item(self, item):
        """BondItemから呼ばれ、ホバー中のアイテムを記録する"""
        self.hovered_item = item

    def keyPressEvent(self, event):
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())
        key = event.key()
        modifiers = event.modifiers()
        
        if not self.window.is_2d_editable:
            return    


        if key == Qt.Key.Key_4:
            # --- 動作1: カーソルが原子/結合上にある場合 (ワンショットでテンプレート配置) ---
            if isinstance(item_at_cursor, (AtomItem, BondItem)):
                
                # ベンゼンテンプレートのパラメータを設定
                n, is_aromatic = 6, True
                points, bonds_info, existing_items = [], [], []
                
                # update_template_preview と同様のロジックで配置情報を計算
                if isinstance(item_at_cursor, AtomItem):
                    p0 = item_at_cursor.pos()
                    l = DEFAULT_BOND_LENGTH
                    direction = QLineF(p0, cursor_pos).unitVector()
                    p1 = p0 + direction.p2() * l if direction.length() > 0 else p0 + QPointF(l, 0)
                    points = self._calculate_polygon_from_edge(p0, p1, n)
                    existing_items = [item_at_cursor]

                elif isinstance(item_at_cursor, BondItem):
                    p0, p1 = item_at_cursor.atom1.pos(), item_at_cursor.atom2.pos()
                    points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=cursor_pos, use_existing_length=True)
                    existing_items = [item_at_cursor.atom1, item_at_cursor.atom2]
                
                if points:
                    bonds_info = [(i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)]
                    
                    # 計算した情報を使って、その場にフラグメントを追加
                    self.add_molecule_fragment(points, bonds_info, existing_items=existing_items)
                    self.window.push_undo_state()

            # --- 動作2: カーソルが空白領域にある場合 (モード切替) ---
            else:
                self.window.set_mode_and_update_toolbar('template_benzene')

            event.accept()
            return

        # --- 0a. ラジカルの変更 (.) ---
        if key == Qt.Key.Key_Period:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                for atom in target_atoms:
                    # ラジカルの状態をトグル (0 -> 1 -> 2 -> 0)
                    atom.prepareGeometryChange()
                    atom.radical = (atom.radical + 1) % 3
                    self.data.atoms[atom.atom_id]['radical'] = atom.radical
                    atom.update_style()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 0b. 電荷の変更 (+/-キー) ---
        if key == Qt.Key.Key_Plus or key == Qt.Key.Key_Minus:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                delta = 1 if key == Qt.Key.Key_Plus else -1
                for atom in target_atoms:
                    atom.prepareGeometryChange()
                    atom.charge += delta
                    self.data.atoms[atom.atom_id]['charge'] = atom.charge
                    atom.update_style()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 1. Atomに対する操作 (元素記号の変更) ---
        if isinstance(item_at_cursor, AtomItem):
            new_symbol = None
            if modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_symbol_map:
                new_symbol = self.key_to_symbol_map[key]
            elif modifiers == Qt.KeyboardModifier.ShiftModifier and key in self.key_to_symbol_map_shift:
                new_symbol = self.key_to_symbol_map_shift[key]

            if new_symbol and item_at_cursor.symbol != new_symbol:
                item_at_cursor.prepareGeometryChange()
                
                item_at_cursor.symbol = new_symbol
                self.data.atoms[item_at_cursor.atom_id]['symbol'] = new_symbol
                item_at_cursor.update_style()


                atoms_to_update = {item_at_cursor}
                for bond in item_at_cursor.bonds:
                    bond.update()
                    other_atom = bond.atom1 if bond.atom2 is item_at_cursor else bond.atom2
                    atoms_to_update.add(other_atom)

                for atom in atoms_to_update:
                    atom.update_style()

                self.window.push_undo_state()
                event.accept()
                return

        # --- 2. Bondに対する操作 (次数・立体化学の変更) ---
        target_bonds = []
        if isinstance(item_at_cursor, BondItem):
            target_bonds = [item_at_cursor]
        else:
            target_bonds = [it for it in self.selectedItems() if isinstance(it, BondItem)]

        if target_bonds:
            any_bond_changed = False
            for bond in target_bonds:
                # 1. 結合の向きを考慮して、データ辞書内の現在のキーを正しく特定する
                id1, id2 = bond.atom1.atom_id, bond.atom2.atom_id
                current_key = None
                if (id1, id2) in self.data.bonds:
                    current_key = (id1, id2)
                elif (id2, id1) in self.data.bonds:
                    current_key = (id2, id1)
                
                if not current_key: continue

                # 2. 変更前の状態を保存
                old_order, old_stereo = bond.order, bond.stereo

                # 3. キー入力に応じてBondItemのプロパティを変更
                if key == Qt.Key.Key_W:
                    if bond.stereo == 1:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1; bond.stereo = 1

                elif key == Qt.Key.Key_D:
                    if bond.stereo == 2:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1; bond.stereo = 2

                elif key == Qt.Key.Key_1 and (bond.order != 1 or bond.stereo != 0):
                    bond.order = 1; bond.stereo = 0
                elif key == Qt.Key.Key_2 and (bond.order != 2 or bond.stereo != 0):
                    bond.order = 2; bond.stereo = 0; needs_update = True
                elif key == Qt.Key.Key_3 and bond.order != 3:
                    bond.order = 3; bond.stereo = 0; needs_update = True

                # 4. 実際に変更があった場合のみデータモデルを更新
                if old_order != bond.order or old_stereo != bond.stereo:
                    any_bond_changed = True
                    
                    # 5. 古いキーでデータを辞書から一度削除
                    bond_data = self.data.bonds.pop(current_key)
                    bond_data['order'] = bond.order
                    bond_data['stereo'] = bond.stereo

                    # 6. 変更後の種類に応じて新しいキーを決定し、再登録する
                    new_key_id1, new_key_id2 = bond.atom1.atom_id, bond.atom2.atom_id
                    if bond.stereo == 0:
                        if new_key_id1 > new_key_id2:
                            new_key_id1, new_key_id2 = new_key_id2, new_key_id1
                    
                    new_key = (new_key_id1, new_key_id2)
                    self.data.bonds[new_key] = bond_data
                    
                    bond.update()

            if any_bond_changed:
                self.window.push_undo_state()
            
            if key in [Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3, Qt.Key.Key_W, Qt.Key.Key_D]:
                event.accept()
                return

        if isinstance(self.hovered_item, BondItem) and self.hovered_item.order == 2:
            if event.key() == Qt.Key.Key_Z:
                self.update_bond_stereo(self.hovered_item, 3)  # Z-isomer
                self.window.push_undo_state()
                event.accept()
                return
            elif event.key() == Qt.Key.Key_E:
                self.update_bond_stereo(self.hovered_item, 4)  # E-isomer
                self.window.push_undo_state()
                event.accept()
                return
                    
        # --- 3. Atomに対する操作 (原子の追加 - マージされた機能) ---
        if key == Qt.Key.Key_1:
            start_atom = None
            if isinstance(item_at_cursor, AtomItem):
                start_atom = item_at_cursor
            else:
                selected_atoms = [item for item in self.selectedItems() if isinstance(item, AtomItem)]
                if len(selected_atoms) == 1:
                    start_atom = selected_atoms[0]

            if start_atom:
                start_pos = start_atom.pos()
                l = DEFAULT_BOND_LENGTH
                new_pos_offset = QPointF(0, -l) # デフォルトのオフセット (上)

                # 接続している原子のリストを取得 (H原子以外)
                neighbor_positions = []
                for bond in start_atom.bonds:
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    if other_atom.symbol != 'H': # 水素原子を無視 (四面体構造の考慮のため)
                        neighbor_positions.append(other_atom.pos())

                num_non_H_neighbors = len(neighbor_positions)
                
                if num_non_H_neighbors == 0:
                    # 結合ゼロ: デフォルト方向
                    new_pos_offset = QPointF(0, -l)
                
                elif num_non_H_neighbors == 1:
                    # 結合1本: 既存結合と約120度（または60度）の角度
                    bond = start_atom.bonds[0]
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    existing_bond_vector = start_pos - other_atom.pos()
                    
                    # 既存の結合から時計回り60度回転 (ベンゼン環のような構造にしやすい)
                    angle_rad = math.radians(60) 
                    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
                    vx, vy = existing_bond_vector.x(), existing_bond_vector.y()
                    new_vx, new_vy = vx * cos_a - vy * sin_a, vx * sin_a + vy * cos_a
                    rotated_vector = QPointF(new_vx, new_vy)
                    line = QLineF(QPointF(0, 0), rotated_vector)
                    line.setLength(l)
                    new_pos_offset = line.p2()

                elif num_non_H_neighbors == 3:

                    bond_vectors_sum = QPointF(0, 0)
                    for pos in neighbor_positions:
                        # start_pos から neighbor_pos へのベクトル
                        vec = pos - start_pos 
                        # 単位ベクトルに変換
                        line_to_other = QLineF(QPointF(0,0), vec)
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2()
                    
                    # SUM_TOLERANCE is now a module-level constant
                    if bond_vectors_sum.manhattanLength() > SUM_TOLERANCE:
                        new_direction_line = QLineF(QPointF(0,0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        new_pos_offset = QPointF(l * 0.7071, -l * 0.7071) 


                else: # 2本または4本以上の場合 (一般的な骨格の継続、または過結合)
                    bond_vectors_sum = QPointF(0, 0)
                    for bond in start_atom.bonds:
                        other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                        line_to_other = QLineF(start_pos, other_atom.pos())
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2() - line_to_other.p1()
                    
                    if bond_vectors_sum.manhattanLength() > 0.01:
                        new_direction_line = QLineF(QPointF(0,0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        # 総和がゼロの場合は、デフォルト（上）
                        new_pos_offset = QPointF(0, -l)


                # SNAP_DISTANCE is a module-level constant
                target_pos = start_pos + new_pos_offset
                
                # 近くに原子を探す
                near_atom = self.find_atom_near(target_pos, tol=SNAP_DISTANCE)
                
                if near_atom and near_atom is not start_atom:
                    # 近くに既存原子があれば結合
                    self.create_bond(start_atom, near_atom)
                else:
                    # 新規原子を作成し結合
                    new_atom_id = self.create_atom('C', target_pos)
                    new_atom_item = self.data.atoms[new_atom_id]['item']
                    self.create_bond(start_atom, new_atom_item)

                self.clearSelection()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 4. 全体に対する操作 (削除、モード切替など) ---
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            if self.temp_line:
                try:
                    if not sip_isdeleted_safe(self.temp_line):
                        try:
                            if getattr(self.temp_line, 'scene', None) and self.temp_line.scene():
                                self.removeItem(self.temp_line)
                        except Exception:
                            pass
                except Exception:
                    try:
                        self.removeItem(self.temp_line)
                    except Exception:
                        pass
                self.temp_line = None; self.start_atom = None; self.start_pos = None
                self.initial_positions_in_event = {}
                event.accept()
                return

            items_to_process = set(self.selectedItems()) 
            # カーソル下のアイテムも削除対象に加える
            if item_at_cursor and isinstance(item_at_cursor, (AtomItem, BondItem)):
                items_to_process.add(item_at_cursor)

            if self.delete_items(items_to_process):
                self.window.push_undo_state()
                self.window.statusBar().showMessage("Deleted selected items.")

            # もしデータモデル内の原子が全て無くなっていたら、シーンをクリアして初期状態に戻す
            if not self.data.atoms:
                # 1. シーン上の全グラフィックアイテムを削除する
                self.clear() 

                # 2. テンプレートプレビューなど、初期状態で必要なアイテムを再生成する
                self.reinitialize_items()
                
                # 3. 結合描画中などの一時的な状態も完全にリセットする
                self.temp_line = None
                self.start_atom = None
                self.start_pos = None
                self.initial_positions_in_event = {}
                
                # このイベントはここで処理完了とする
                event.accept()
                return
    
            # 描画の強制更新
            if self.views():
                self.views()[0].viewport().update() 
                QApplication.processEvents()
    
                event.accept()
                return
        

        if key == Qt.Key.Key_Space:
            if self.mode != 'select':
                self.window.activate_select_mode()
            else:
                self.window.select_all()
            event.accept()
            return

        # グローバルな描画モード切替
        mode_to_set = None

        # 1. 原子描画モードへの切り替え
        symbol_for_mode_change = None
        if modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_symbol_map:
            symbol_for_mode_change = self.key_to_symbol_map[key]
        elif modifiers == Qt.KeyboardModifier.ShiftModifier and key in self.key_to_symbol_map_shift:
            symbol_for_mode_change = self.key_to_symbol_map_shift[key]
        
        if symbol_for_mode_change:
            mode_to_set = f'atom_{symbol_for_mode_change}'

        # 2. 結合描画モードへの切り替え
        elif modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_bond_mode_map:
            mode_to_set = self.key_to_bond_mode_map[key]

        # モードが決定されていれば、モード変更を実行
        if mode_to_set:
            if hasattr(self.window, 'set_mode_and_update_toolbar'):
                 self.window.set_mode_and_update_toolbar(mode_to_set)
                 event.accept()
                 return
        
        # --- どの操作にも当てはまらない場合 ---
        super().keyPressEvent(event)
        
    def find_atom_near(self, pos, tol=14.0):
        # Create a small search rectangle around the position
        search_rect = QRectF(pos.x() - tol, pos.y() - tol, 2 * tol, 2 * tol)
        nearby_items = self.items(search_rect)

        for it in nearby_items:
            if isinstance(it, AtomItem):
                # Check the precise distance only for candidate items
                if QLineF(it.pos(), pos).length() <= tol:
                    return it
        return None

    def find_bond_between(self, atom1, atom2):
        for b in atom1.bonds:
            if (b.atom1 is atom1 and b.atom2 is atom2) or \
               (b.atom1 is atom2 and b.atom2 is atom1):
                return b
        return None

    def update_bond_stereo(self, bond_item, new_stereo):
        """結合の立体化学を更新する共通メソッド"""
        try:
            if bond_item is None:
                print("Error: bond_item is None in update_bond_stereo")
                return
                
            if bond_item.order != 2 or bond_item.stereo == new_stereo:
                return

            if not hasattr(bond_item, 'atom1') or not hasattr(bond_item, 'atom2'):
                print("Error: bond_item missing atom references")
                return
                
            if bond_item.atom1 is None or bond_item.atom2 is None:
                print("Error: bond_item has None atom references")
                return
                
            if not hasattr(bond_item.atom1, 'atom_id') or not hasattr(bond_item.atom2, 'atom_id'):
                print("Error: bond atoms missing atom_id")
                return

            id1, id2 = bond_item.atom1.atom_id, bond_item.atom2.atom_id

            # E/Z結合は方向性を持つため、キーは(id1, id2)のまま探す
            key_to_update = (id1, id2)
            if key_to_update not in self.data.bonds:
                # Wedge/Dashなど、逆順で登録されている可能性も考慮
                key_to_update = (id2, id1)
                if key_to_update not in self.data.bonds:
                    # Log error instead of printing to console
                    if hasattr(self.window, 'statusBar'):
                        self.window.statusBar().showMessage(f"Warning: Bond between atoms {id1} and {id2} not found in data model.", 3000)
                    print(f"Error: Bond key not found: {id1}-{id2} or {id2}-{id1}")
                    return
                    
            # Update data model
            self.data.bonds[key_to_update]['stereo'] = new_stereo
            
            # Update visual representation
            bond_item.set_stereo(new_stereo)
            
            self.data_changed_in_event = True
            
        except Exception as e:
            print(f"Error in update_bond_stereo: {e}")
            
            traceback.print_exc()
            if hasattr(self.window, 'statusBar'):
                self.window.statusBar().showMessage(f"Error updating bond stereochemistry: {e}", 5000)
