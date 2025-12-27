#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""


"""
main_window_edit_actions.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowEditActions
"""


import numpy as np
import pickle
import math
import io
import itertools
import traceback

from collections import deque

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem


# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication
)

from PyQt6.QtGui import (
    QCursor
)


from PyQt6.QtCore import (
    QPointF, QLineF, QMimeData, QByteArray, QTimer
)


# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None
    
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
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import CLIPBOARD_MIME_TYPE
    from .molecular_data import MolecularData
    from .atom_item import AtomItem
    from .bond_item import BondItem
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import CLIPBOARD_MIME_TYPE
    from modules.molecular_data import MolecularData
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem


try:
    # Import the shared SIP helper used across the package. This is
    # defined in modules/__init__.py and centralizes sip.isdeleted checks.
    from . import sip_isdeleted_safe
except Exception:
    from modules import sip_isdeleted_safe

# --- クラス定義 ---
class MainWindowEditActions(object):
    """ main_window.py から分離された機能クラス """


    def copy_selection(self):
        """選択された原子と結合をクリップボードにコピーする"""
        try:
            selected_atoms = [item for item in self.scene.selectedItems() if isinstance(item, AtomItem)]
            if not selected_atoms:
                return

            # 選択された原子のIDセットを作成
            selected_atom_ids = {atom.atom_id for atom in selected_atoms}
            
            # 選択された原子の幾何学的中心を計算
            center = QPointF(
                sum(atom.pos().x() for atom in selected_atoms) / len(selected_atoms),
                sum(atom.pos().y() for atom in selected_atoms) / len(selected_atoms)
            )
            
            # コピー対象の原子データをリストに格納（位置は中心からの相対座標）
            # 同時に、元のatom_idから新しいインデックス(0, 1, 2...)へのマッピングを作成
            atom_id_to_idx_map = {}
            fragment_atoms = []
            for i, atom in enumerate(selected_atoms):
                atom_id_to_idx_map[atom.atom_id] = i
                fragment_atoms.append({
                    'symbol': atom.symbol,
                    'rel_pos': atom.pos() - center,
                    'charge': atom.charge,
                    'radical': atom.radical,
                })
                
            # 選択された原子同士を結ぶ結合のみをリストに格納
            fragment_bonds = []
            for (id1, id2), bond_data in self.data.bonds.items():
                if id1 in selected_atom_ids and id2 in selected_atom_ids:
                    fragment_bonds.append({
                        'idx1': atom_id_to_idx_map[id1],
                        'idx2': atom_id_to_idx_map[id2],
                        'order': bond_data['order'],
                        'stereo': bond_data.get('stereo', 0),  # E/Z立体化学情報も保存
                    })

            # pickleを使ってデータをバイト配列にシリアライズ
            data_to_pickle = {'atoms': fragment_atoms, 'bonds': fragment_bonds}
            byte_array = QByteArray()
            buffer = io.BytesIO()
            pickle.dump(data_to_pickle, buffer)
            byte_array.append(buffer.getvalue())

            # カスタムMIMEタイプでクリップボードに設定
            mime_data = QMimeData()
            mime_data.setData(CLIPBOARD_MIME_TYPE, byte_array)
            QApplication.clipboard().setMimeData(mime_data)
            self.statusBar().showMessage(f"Copied {len(fragment_atoms)} atoms and {len(fragment_bonds)} bonds.")
            
        except Exception as e:
            print(f"Error during copy operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during copy operation: {e}")



    def cut_selection(self):
        """選択されたアイテムを切り取り（コピーしてから削除）"""
        try:
            selected_items = self.scene.selectedItems()
            if not selected_items:
                return
            
            # 最初にコピー処理を実行
            self.copy_selection()
            
            if self.scene.delete_items(set(selected_items)):
                self.push_undo_state()
                self.statusBar().showMessage("Cut selection.", 2000)
                
        except Exception as e:
            print(f"Error during cut operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during cut operation: {e}")



    def paste_from_clipboard(self):
        """クリップボードから分子フラグメントを貼り付け"""
        try:
            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            if not mime_data.hasFormat(CLIPBOARD_MIME_TYPE):
                return

            byte_array = mime_data.data(CLIPBOARD_MIME_TYPE)
            buffer = io.BytesIO(byte_array)
            try:
                fragment_data = pickle.load(buffer)
            except pickle.UnpicklingError:
                self.statusBar().showMessage("Error: Invalid clipboard data format")
                return
            
            paste_center_pos = self.view_2d.mapToScene(self.view_2d.mapFromGlobal(QCursor.pos()))
            self.scene.clearSelection()

            new_atoms = []
            for atom_data in fragment_data['atoms']:
                pos = paste_center_pos + atom_data['rel_pos']
                new_id = self.scene.create_atom(
                    atom_data['symbol'], pos,
                    charge=atom_data.get('charge', 0),
                    radical=atom_data.get('radical', 0)
                )
                new_item = self.data.atoms[new_id]['item']
                new_atoms.append(new_item)
                new_item.setSelected(True)

            for bond_data in fragment_data['bonds']:
                atom1 = new_atoms[bond_data['idx1']]
                atom2 = new_atoms[bond_data['idx2']]
                self.scene.create_bond(
                    atom1, atom2,
                    bond_order=bond_data.get('order', 1),
                    bond_stereo=bond_data.get('stereo', 0)  # E/Z立体化学情報も復元
                )
            
            self.push_undo_state()
            self.statusBar().showMessage(f"Pasted {len(fragment_data['atoms'])} atoms and {len(fragment_data['bonds'])} bonds.", 2000)
            
        except Exception as e:
            print(f"Error during paste operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during paste operation: {e}")
        self.statusBar().showMessage(f"Pasted {len(new_atoms)} atoms.", 2000)
        self.activate_select_mode()



    def remove_hydrogen_atoms(self):
        """2Dビューで水素原子とその結合を削除する"""
        try:
            # Collect hydrogen atom items robustly (store atom_id -> item)
            hydrogen_map = {}

            # Iterate over a snapshot of atoms to avoid "dictionary changed size"
            for atom_id, atom_data in list(self.data.atoms.items()):
                try:
                    if atom_data.get('symbol') != 'H':
                        continue
                    item = atom_data.get('item')
                    # Only collect live AtomItem wrappers
                    if item is None:
                        continue
                    if sip_isdeleted_safe(item):
                        continue
                    if not isinstance(item, AtomItem):
                        continue
                    # Prefer storing by original atom id to detect actual removals later
                    hydrogen_map[atom_id] = item
                except Exception:
                    # Ignore problematic entries and continue scanning
                    continue

            if not hydrogen_map:
                self.statusBar().showMessage("No hydrogen atoms found to remove.", 2000)
                return

            # To avoid blocking the UI or causing large, monolithic deletions that may
            # trigger internal re-entrancy issues, delete in batches and process UI events
            items = list(hydrogen_map.values())
            total = len(items)
            batch_size = 200  # tuned conservative batch size
            deleted_any = False

            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                batch = set()
                # Filter out items that are already deleted or invalid just before deletion
                for it in items[start:end]:
                    try:
                        if it is None:
                            continue
                        if sip_isdeleted_safe(it):
                            continue
                        if not isinstance(it, AtomItem):
                            continue
                        batch.add(it)
                    except Exception:
                        continue

                if not batch:
                    # Nothing valid to delete in this batch
                    continue

                try:
                    # scene.delete_items is expected to handle bond cleanup; call it per-batch
                    success = False
                    try:
                        success = bool(self.scene.delete_items(batch))
                    except Exception:
                        # If scene.delete_items raises for a batch, attempt a safe per-item fallback
                        success = False

                    if not success:
                        # Fallback: try deleting items one-by-one to isolate problematic items
                        for it in list(batch):
                            try:
                                # Use scene.delete_items for single-item as well
                                ok = bool(self.scene.delete_items({it}))
                                if ok:
                                    deleted_any = True
                            except Exception:
                                # If single deletion also fails, skip that item
                                continue
                    else:
                        deleted_any = True

                except Exception:
                    # Continue with next batch on unexpected errors
                    continue

                # Allow the GUI to process events between batches to remain responsive
                try:
                    QApplication.processEvents()
                except Exception:
                    pass

            # Determine how many hydrogens actually were removed by re-scanning data
            remaining_h = 0
            try:
                for _, atom_data in list(self.data.atoms.items()):
                    try:
                        if atom_data.get('symbol') == 'H':
                            remaining_h += 1
                    except Exception:
                        continue
            except Exception:
                remaining_h = 0

            removed_count = max(0, len(hydrogen_map) - remaining_h)

            if removed_count > 0:
                # Only push a single undo state once for the whole operation
                try:
                    self.push_undo_state()
                except Exception:
                    # Do not allow undo stack problems to crash the app
                    pass
                self.statusBar().showMessage(f"Removed {removed_count} hydrogen atoms.", 2000)
            else:
                # If nothing removed but we attempted, show an informative message
                if deleted_any:
                    # Deleted something but couldn't determine count reliably
                    self.statusBar().showMessage("Removed hydrogen atoms (count unknown).", 2000)
                else:
                    self.statusBar().showMessage("Failed to remove hydrogen atoms or none found.")

        except Exception as e:
            # Capture and log unexpected errors but don't let them crash the UI
            print(f"Error during hydrogen removal: {e}")
            traceback.print_exc()
            try:
                self.statusBar().showMessage(f"Error removing hydrogen atoms: {e}")
            except Exception:
                pass



    def add_hydrogen_atoms(self):
        """RDKitで各原子の暗黙の水素数を調べ、その数だけ明示的な水素原子と単結合を作成する（2Dビュー）。

        実装上の仮定:
        - `self.data.to_rdkit_mol()` は各RDKit原子に `_original_atom_id` プロパティを設定している。
        - 原子の2D座標は `self.data.atoms[orig_id]['item'].pos()` で得られる。
        - 新しい原子は `self.scene.create_atom(symbol, pos, ...)` で追加し、
          結合は `self.scene.create_bond(atom_item, hydrogen_item, bond_order=1)` で作成する。
        """
        try:
            

            mol = self.data.to_rdkit_mol(use_2d_stereo=False)
            if not mol or mol.GetNumAtoms() == 0:
                self.statusBar().showMessage("No molecule available to compute hydrogens.", 2000)
                return

            added_count = 0
            added_items = []

            # すべてのRDKit原子について暗黙水素数を確認
            for idx in range(mol.GetNumAtoms()):
                rd_atom = mol.GetAtomWithIdx(idx)
                try:
                    orig_id = rd_atom.GetIntProp("_original_atom_id")
                except Exception:
                    # 元のエディタ側のIDがない場合はスキップ
                    continue

                if orig_id not in self.data.atoms:
                    continue

                # 暗黙水素数を優先して取得。存在しない場合は総水素数 - 明示水素数を使用
                implicit_h = int(rd_atom.GetNumImplicitHs()) if hasattr(rd_atom, 'GetNumImplicitHs') else 0
                if implicit_h is None or implicit_h < 0:
                    implicit_h = 0
                if implicit_h == 0:
                    # フォールバック
                    try:
                        total_h = int(rd_atom.GetTotalNumHs())
                        explicit_h = int(rd_atom.GetNumExplicitHs()) if hasattr(rd_atom, 'GetNumExplicitHs') else 0
                        implicit_h = max(0, total_h - explicit_h)
                    except Exception:
                        implicit_h = 0

                if implicit_h <= 0:
                    continue

                parent_item = self.data.atoms[orig_id]['item']
                parent_pos = parent_item.pos()

                # 周囲の近接原子の方向を取得して、水素を邪魔しないように角度を決定
                neighbor_angles = []
                try:
                    for (a1, a2), bdata in self.data.bonds.items():
                        # 対象原子に結合している近傍の原子角度を収集する。
                        # ただし既存の水素は配置に影響させない（すでにあるHで埋めない）。
                        try:
                            if a1 == orig_id and a2 in self.data.atoms:
                                neigh = self.data.atoms[a2]
                                if neigh.get('symbol') == 'H':
                                    continue
                                if neigh.get('item') is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get('item')):
                                    continue
                                vec = neigh['item'].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                            elif a2 == orig_id and a1 in self.data.atoms:
                                neigh = self.data.atoms[a1]
                                if neigh.get('symbol') == 'H':
                                    continue
                                if neigh.get('item') is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get('item')):
                                    continue
                                vec = neigh['item'].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                        except Exception:
                            # 個々の近傍読み取りの問題は無視して続行
                            continue
                except Exception:
                    neighbor_angles = []

                # 画面上の適当な結合長（ピクセル）を使用
                bond_length = 75

                # ヘルパー: 指定インデックスの水素に使うbond_stereoを決定
                def _choose_stereo(i):
                    # 0: plain, 1: wedge, 2: dash, 3: plain, 4+: all plain
                    if i == 0:
                        return 0
                    if i == 1:
                        return 1
                    if i == 2:
                        return 2
                    return 0  #4th+ hydrogens are all plain

                # 角度配置を改善: 既存の結合角度の最大ギャップを見つけ、
                # そこに水素を均等配置する。既存結合が無ければ全周に均等配置。
                target_angles = []
                try:
                    if not neighbor_angles:
                        # 既存結合が無い -> 全円周に均等配置
                        for h_idx in range(implicit_h):
                            angle = (2.0 * math.pi * h_idx) / implicit_h
                            target_angles.append(angle)
                    else:
                        # 正規化してソート
                        angs = [((a + 2.0 * math.pi) if a < 0 else a) for a in neighbor_angles]
                        angs = sorted(angs)
                        # ギャップを計算（循環含む）
                        gaps = []  # list of (gap_size, start_angle, end_angle)
                        for i in range(len(angs)):
                            a1 = angs[i]
                            a2 = angs[(i + 1) % len(angs)]
                            if i == len(angs) - 1:
                                # wrap-around gap
                                gap = (a2 + 2.0 * math.pi) - a1
                                start = a1
                                end = a2 + 2.0 * math.pi
                            else:
                                gap = a2 - a1
                                start = a1
                                end = a2
                            gaps.append((gap, start, end))

                        # 最大ギャップを選ぶ
                        gaps.sort(key=lambda x: x[0], reverse=True)
                        max_gap, gstart, gend = gaps[0]
                        # もし最大ギャップが小さい（つまり周りに均等に原子がある）でも
                        # そのギャップ内に均等配置することで既存結合と重ならないようにする
                        # ギャップ内に implicit_h 個を等間隔で配置（分割数 = implicit_h + 1）
                        for i in range(implicit_h):
                            seg = max_gap / (implicit_h + 1)
                            angle = gstart + (i + 1) * seg
                            # 折り返しを戻して 0..2pi に正規化
                            angle = angle % (2.0 * math.pi)
                            target_angles.append(angle)
                except Exception:
                    # フォールバック: 単純な等間隔配置
                    for h_idx in range(implicit_h):
                        angle = (2.0 * math.pi * h_idx) / implicit_h
                        target_angles.append(angle)

                # 角度から位置を計算して原子と結合を追加
                for h_idx, angle in enumerate(target_angles):
                    dx = bond_length * math.cos(angle)
                    dy = bond_length * math.sin(angle)
                    pos = QPointF(parent_pos.x() + dx, parent_pos.y() + dy)

                    # 新しい水素原子を作成
                    try:
                        new_id = self.scene.create_atom('H', pos)
                        new_item = self.data.atoms[new_id]['item']
                        # bond_stereo を指定（最初は plain=0, 次に wedge/dash）
                        stereo = _choose_stereo(h_idx)
                        self.scene.create_bond(parent_item, new_item, bond_order=1, bond_stereo=stereo)
                        added_items.append(new_item)
                        added_count += 1
                    except Exception as e:
                        # 個々の追加失敗はログに残して続行
                        print(f"Failed to add H for atom {orig_id}: {e}")

            if added_count > 0:
                self.push_undo_state()
                self.statusBar().showMessage(f"Added {added_count} hydrogen atoms.", 2000)
                # 選択を有効化して追加した原子を選択状態にする
                try:
                    self.scene.clearSelection()
                    for it in added_items:
                        it.setSelected(True)
                except Exception:
                    pass
            else:
                self.statusBar().showMessage("No implicit hydrogens found to add.", 2000)

        except Exception as e:
            print(f"Error during hydrogen addition: {e}")
            traceback.print_exc()
            self.statusBar().showMessage(f"Error adding hydrogen atoms: {e}")



    def update_edit_menu_actions(self):
        """選択状態やクリップボードの状態に応じて編集メニューを更新"""
        try:
            has_selection = len(self.scene.selectedItems()) > 0
            self.cut_action.setEnabled(has_selection)
            self.copy_action.setEnabled(has_selection)
            
            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            self.paste_action.setEnabled(mime_data is not None and mime_data.hasFormat(CLIPBOARD_MIME_TYPE))
        except RuntimeError:
            pass




    def select_all(self):
        for item in self.scene.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.setSelected(True)



    def clear_all(self):
        # 未保存の変更があるかチェック
        if not self.check_unsaved_changes():
            return  # ユーザーがキャンセルした場合は何もしない

        self.restore_ui_for_editing()

        # データが存在しない場合は何もしない
        if not self.data.atoms and self.current_mol is None:
            return
        
        # 3Dモードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化
        
        # 3D原子選択をクリア
        self.clear_3d_selection()
        
        self.dragged_atom_info = None
            
        # 2Dエディタをクリアする（Undoスタックにはプッシュしない）
        self.clear_2d_editor(push_to_undo=False)
        
        # 3Dモデルをクリアする
        self.current_mol = None
        self.plotter.clear()
        self.constraints_3d = []
        
        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        # Undo/Redoスタックをリセットする
        self.reset_undo_stack()
        
        # ファイル状態をリセット（新規ファイル状態に）
        self.has_unsaved_changes = False
        self.current_file_path = None
        self.update_window_title()
        
        # 2Dビューのズームをリセット
        self.reset_zoom()
        
        # シーンとビューの明示的な更新
        self.scene.update()
        if self.view_2d:
            self.view_2d.viewport().update()

        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        # 3Dプロッターの再描画
        self.plotter.render()
        
        # メニューテキストと状態を更新（分子がクリアされたので通常の表示に戻す）
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()
        
        # アプリケーションのイベントループを強制的に処理し、画面の再描画を確実に行う
        QApplication.processEvents()
        
        self.statusBar().showMessage("Cleared all data.")
        


    def clear_2d_editor(self, push_to_undo=True):
        self.data = MolecularData()
        self.scene.data = self.data
        self.scene.clear()
        self.scene.reinitialize_items()
        self.is_xyz_derived = False  # 2Dエディタをクリアする際にXYZ由来フラグもリセット
        
        # 測定ラベルもクリア
        self.clear_2d_measurement_labels()
        
        # Clear 3D data and disable 3D-related menus
        self.current_mol = None
        self.plotter.clear()
        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        if push_to_undo:
            self.push_undo_state()



    def update_implicit_hydrogens(self):
        """現在の2D構造に基づいて各原子の暗黙の水素数を計算し、AtomItemに反映する"""
        # Quick guards: nothing to do if no atoms or no QApplication
        if not self.data.atoms:
            return

        # If called from non-GUI thread, schedule the heavy RDKit work here but
        # always perform UI mutations on the main thread via QTimer.singleShot.
        try:
            # Bump a local token to identify this request. The closure we
            # schedule below will capture `my_token` and will only apply UI
            # changes if the token still matches the most recent global
            # counter. This avoids applying stale updates after deletions or
            # teardown.
            try:
                self._ih_update_counter += 1
            except Exception:
                self._ih_update_counter = getattr(self, '_ih_update_counter', 0) or 1
            my_token = self._ih_update_counter

            mol = None
            try:
                mol = self.data.to_rdkit_mol()
            except Exception:
                mol = None

            # Build a mapping of original_id -> hydrogen count without touching Qt items
            h_count_map = {}

            if mol is None:
                # Invalid/unsanitizable structure: reset all counts to 0
                for atom_id in list(self.data.atoms.keys()):
                    h_count_map[atom_id] = 0
            else:
                for atom in mol.GetAtoms():
                    try:
                        if not atom.HasProp("_original_atom_id"):
                            continue
                        original_id = atom.GetIntProp("_original_atom_id")

                        # Robust retrieval of H counts: prefer implicit, fallback to total or 0
                        try:
                            h_count = int(atom.GetNumImplicitHs())
                        except Exception:
                            try:
                                h_count = int(atom.GetTotalNumHs())
                            except Exception:
                                h_count = 0

                        h_count_map[int(original_id)] = h_count
                    except Exception:
                        # Skip problematic RDKit atoms
                        continue

            # Compute a per-atom problem map (original_id -> bool) so the
            # UI closure can safely set AtomItem.has_problem on the main thread.
            problem_map = {}
            try:
                if mol is not None:
                    try:
                        problems = Chem.DetectChemistryProblems(mol)
                    except Exception:
                        problems = None

                    if problems:
                        for prob in problems:
                            try:
                                atom_idx = prob.GetAtomIdx()
                                rd_atom = mol.GetAtomWithIdx(atom_idx)
                                if rd_atom and rd_atom.HasProp("_original_atom_id"):
                                    orig = int(rd_atom.GetIntProp("_original_atom_id"))
                                    problem_map[orig] = True
                            except Exception:
                                continue
                else:
                    # Fallback: use a lightweight valence heuristic similar to
                    # check_chemistry_problems_fallback() so we still flag atoms
                    # when RDKit conversion wasn't possible.
                    for atom_id, atom_data in self.data.atoms.items():
                        try:
                            symbol = atom_data.get('symbol')
                            charge = atom_data.get('charge', 0)
                            bond_count = 0
                            for (id1, id2), bond_data in self.data.bonds.items():
                                if id1 == atom_id or id2 == atom_id:
                                    bond_count += bond_data.get('order', 1)

                            is_problematic = False
                            if symbol == 'C' and bond_count > 4:
                                is_problematic = True
                            elif symbol == 'N' and bond_count > 3 and charge == 0:
                                is_problematic = True
                            elif symbol == 'O' and bond_count > 2 and charge == 0:
                                is_problematic = True
                            elif symbol == 'H' and bond_count > 1:
                                is_problematic = True
                            elif symbol in ['F', 'Cl', 'Br', 'I'] and bond_count > 1 and charge == 0:
                                is_problematic = True

                            if is_problematic:
                                problem_map[atom_id] = True
                        except Exception:
                            continue
            except Exception:
                # If any unexpected error occurs while building the map, fall back
                # to an empty map so we don't accidentally crash the UI.
                problem_map = {}

            # Schedule UI updates on the main thread to avoid calling Qt methods from
            # background threads or during teardown (which can crash the C++ layer).
            def _apply_ui_updates():
                # If the global counter changed since this closure was
                # created, bail out — the update is stale.
                try:
                    if my_token != getattr(self, '_ih_update_counter', None):
                        return
                except Exception:
                    # If anything goes wrong checking the token, be conservative
                    # and skip the update to avoid touching possibly-damaged
                    # Qt wrappers.
                    return

                # Work on a shallow copy/snapshot of the data.atoms mapping so
                # that concurrent mutations won't raise KeyError during
                # iteration. We still defensively check each item below.
                try:
                    atoms_snapshot = dict(self.data.atoms)
                except Exception:
                    atoms_snapshot = {}
                # Prefer the module-level SIP helper to avoid repeated imports
                # and centralize exception handling. Use the safe wrapper
                # `sip_isdeleted_safe` provided by the package which already
                # handles the optional presence of sip.isdeleted.
                is_deleted_func = sip_isdeleted_safe

                items_to_update = []
                for atom_id, atom_data in atoms_snapshot.items():
                    try:
                        item = atom_data.get('item')
                        if not item:
                            continue

                        # If sip.isdeleted is available, skip deleted C++ wrappers
                        try:
                            if is_deleted_func and is_deleted_func(item):
                                continue
                        except Exception:
                            # If sip check itself fails, continue with other lightweight guards
                            pass

                        # If the item is no longer in a scene, skip updating it to avoid
                        # touching partially-deleted objects during scene teardown.
                        try:
                            sc = item.scene() if hasattr(item, 'scene') else None
                            if sc is None:
                                continue
                        except Exception:
                            # Accessing scene() might fail for a damaged object; skip it
                            continue

                        # Desired new count (default to 0 if not computed)
                        new_count = h_count_map.get(atom_id, 0)

                        current = getattr(item, 'implicit_h_count', None)
                        current_prob = getattr(item, 'has_problem', False)
                        desired_prob = problem_map.get(atom_id, False)

                        # If neither the implicit-H count nor the problem flag
                        # changed, skip this item.
                        if current == new_count and current_prob == desired_prob:
                            continue

                        # Only prepare a geometry change if the implicit H count
                        # changes (this may affect the item's bounding rect).
                        need_geometry = (current != new_count)
                        try:
                            if need_geometry and hasattr(item, 'prepareGeometryChange'):
                                try:
                                    item.prepareGeometryChange()
                                except Exception:
                                    pass

                            # Apply implicit hydrogen count (guarded)
                            try:
                                item.implicit_h_count = new_count
                            except Exception:
                                # If setting the count fails, continue but still
                                # attempt to set the problem flag below.
                                pass

                            # Apply problem flag (visual red-outline)
                            try:
                                item.has_problem = bool(desired_prob)
                            except Exception:
                                pass

                            # Ensure the item is updated in the scene so paint() runs
                            # when either geometry or problem-flag changed.
                            items_to_update.append(item)
                        except Exception:
                            # Non-fatal: skip problematic items
                            continue

                    except Exception:
                        continue

                # Trigger updates once for unique items; wrap in try/except to avoid crashes
                # Trigger updates once for unique items; dedupe by object id so
                # we don't attempt to hash QGraphicsItem wrappers which may
                # behave oddly when partially deleted.
                seen = set()
                for it in items_to_update:
                    try:
                        if it is None:
                            continue
                        oid = id(it)
                        if oid in seen:
                            continue
                        seen.add(oid)
                        if hasattr(it, 'update'):
                            try:
                                it.update()
                            except Exception:
                                # ignore update errors for robustness
                                pass
                    except Exception:
                        # Ignore any unexpected errors when touching the item
                        continue

            # Always schedule on main thread asynchronously
            try:
                QTimer.singleShot(0, _apply_ui_updates)
            except Exception:
                # Fallback: try to call directly (best-effort)
                try:
                    _apply_ui_updates()
                except Exception:
                    pass

        except Exception:
            # Make sure update failures never crash the application
            pass




    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")
        
        # 最初に既存の化学的問題フラグをクリア
        self.scene.clear_all_problem_flags()
        
        # 2Dエディタに原子が存在しない場合
        if not self.data.atoms:
            self.statusBar().showMessage("Error: No atoms to optimize.")
            return
        
        mol = self.data.to_rdkit_mol()
        if mol is None or mol.GetNumAtoms() == 0:
            # RDKit変換が失敗した場合は化学的問題をチェック
            self.check_chemistry_problems_fallback()
            return

        try:
            # 安定版：原子IDとRDKit座標の確実なマッピング
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())
            new_positions_map = {}
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            for rdkit_atom in mol.GetAtoms():
                original_id = rdkit_atom.GetIntProp("_original_atom_id")
                new_positions_map[original_id] = conf.GetAtomPosition(rdkit_atom.GetIdx())

            if not new_positions_map:
                self.statusBar().showMessage("Optimization failed to generate coordinates."); return

            target_atom_items = [self.data.atoms[atom_id]['item'] for atom_id in new_positions_map.keys() if atom_id in self.data.atoms and 'item' in self.data.atoms[atom_id]]
            if not target_atom_items:
                self.statusBar().showMessage("Error: Atom items not found for optimized atoms."); return

            # 元の図形の中心を維持
            #original_center_x = sum(item.pos().x() for item in target_atom_items) / len(target_atom_items)
            #original_center_y = sum(item.pos().y() for item in target_atom_items) / len(target_atom_items)

            positions = list(new_positions_map.values())
            rdkit_cx = sum(p.x for p in positions) / len(positions)
            rdkit_cy = sum(p.y for p in positions) / len(positions)

            SCALE = 50.0

            # 新しい座標を適用
            for atom_id, rdkit_pos in new_positions_map.items():
                if atom_id in self.data.atoms:
                    item = self.data.atoms[atom_id]['item']
                    sx = ((rdkit_pos.x - rdkit_cx) * SCALE) + view_center.x()
                    sy = (-(rdkit_pos.y - rdkit_cy) * SCALE) + view_center.y()
                    new_scene_pos = QPointF(sx, sy)
                    item.setPos(new_scene_pos)
                    self.data.atoms[atom_id]['pos'] = new_scene_pos

            # 最終的な座標に基づき、全ての結合表示を一度に更新
            # Guard against partially-deleted Qt wrappers: skip items that
            # SIP reports as deleted or which are no longer in a scene.
            for bond_data in self.data.bonds.values():
                item = bond_data.get('item') if bond_data else None
                if not item:
                    continue
                try:
                    # If SIP is available, skip wrappers whose C++ object is gone
                    if sip_isdeleted_safe(item):
                        continue
                except Exception:
                    # If the sip check fails, continue with other lightweight guards
                    pass
                try:
                    sc = None
                    try:
                        sc = item.scene() if hasattr(item, 'scene') else None
                    except Exception:
                        sc = None
                    if sc is None:
                        continue
                    try:
                        item.update_position()
                    except Exception:
                        # Best-effort: skip any bond items that raise when updating
                        continue
                except Exception:
                    continue

            # 重なり解消ロジックを実行
            self. resolve_overlapping_groups()
            
            # 測定ラベルの位置を更新
            self.update_2d_measurement_labels()
            
            # シーン全体の再描画を要求
            self.scene.update()

            self.statusBar().showMessage("2D structure optimization successful.")
            self.push_undo_state()

        except Exception as e:
            self.statusBar().showMessage(f"Error during 2D optimization: {e}")
        finally:
            self.view_2d.setFocus()



    def resolve_overlapping_groups(self):
        """
        誤差範囲で完全に重なっている原子のグループを検出し、
        IDが大きい方のフラグメントを左下に平行移動して解消する。
        """

        # --- パラメータ設定 ---
        # 重なっているとみなす距離の閾値。構造に合わせて調整してください。
        OVERLAP_THRESHOLD = 0.5  
        # 左下へ移動させる距離。
        MOVE_DISTANCE = 20

        # self.data.atoms.values() から item を安全に取得
        all_atom_items = [
            data['item'] for data in self.data.atoms.values() 
            if data and 'item' in data
        ]

        if len(all_atom_items) < 2:
            return

        # --- ステップ1: 重なっている原子ペアを全てリストアップ ---
        overlapping_pairs = []
        for item1, item2 in itertools.combinations(all_atom_items, 2):
            # 結合で直接結ばれているペアは重なりと見なさない
            if self.scene.find_bond_between(item1, item2):
                continue

            dist = QLineF(item1.pos(), item2.pos()).length()
            if dist < OVERLAP_THRESHOLD:
                overlapping_pairs.append((item1, item2))

        if not overlapping_pairs:
            self.statusBar().showMessage("No overlapping atoms found.", 2000)
            return

        # --- ステップ2: Union-Findアルゴリズムで重なりグループを構築 ---
        # 各原子がどのグループに属するかを管理する
        parent = {item.atom_id: item.atom_id for item in all_atom_items}

        def find_set(atom_id):
            # atom_idが属するグループの代表（ルート）を見つける
            if parent[atom_id] == atom_id:
                return atom_id
            parent[atom_id] = find_set(parent[atom_id])  # 経路圧縮による最適化
            return parent[atom_id]

        def unite_sets(id1, id2):
            # 2つの原子が属するグループを統合する
            root1 = find_set(id1)
            root2 = find_set(id2)
            if root1 != root2:
                parent[root2] = root1

        for item1, item2 in overlapping_pairs:
            unite_sets(item1.atom_id, item2.atom_id)

        # --- ステップ3: グループごとに移動計画を立てる ---
        # 同じ代表を持つ原子でグループを辞書にまとめる
        groups_by_root = {}
        for item in all_atom_items:
            root_id = find_set(item.atom_id)
            if root_id not in groups_by_root:
                groups_by_root[root_id] = []
            groups_by_root[root_id].append(item.atom_id)

        move_operations = []
        processed_roots = set()

        for root_id, group_atom_ids in groups_by_root.items():
            # 処理済みのグループや、メンバーが1つしかないグループはスキップ
            if root_id in processed_roots or len(group_atom_ids) < 2:
                continue
            processed_roots.add(root_id)

            # 3a: グループを、結合に基づいたフラグメントに分割する (BFSを使用)
            fragments = []
            visited_in_group = set()
            group_atom_ids_set = set(group_atom_ids)

            for atom_id in group_atom_ids:
                if atom_id not in visited_in_group:
                    current_fragment = set()
                    q = deque([atom_id])
                    visited_in_group.add(atom_id)
                    current_fragment.add(atom_id)

                    while q:
                        current_id = q.popleft()
                        # 隣接リスト self.adjacency_list があれば、ここでの探索が高速になります
                        for neighbor_id in self.data.adjacency_list.get(current_id, []):
                            if neighbor_id in group_atom_ids_set and neighbor_id not in visited_in_group:
                                visited_in_group.add(neighbor_id)
                                current_fragment.add(neighbor_id)
                                q.append(neighbor_id)
                    fragments.append(current_fragment)

            if len(fragments) < 2:
                continue  # 複数のフラグメントが重なっていない場合

            # 3b: 移動するフラグメントを決定する
            # このグループの重なりの原因となった代表ペアを一つ探す
            rep_item1, rep_item2 = None, None
            for i1, i2 in overlapping_pairs:
                if find_set(i1.atom_id) == root_id:
                    rep_item1, rep_item2 = i1, i2
                    break

            if not rep_item1: continue

            # 代表ペアがそれぞれどのフラグメントに属するかを見つける
            frag1 = next((f for f in fragments if rep_item1.atom_id in f), None)
            frag2 = next((f for f in fragments if rep_item2.atom_id in f), None)

            # 同一フラグメント内の重なりなどはスキップ
            if not frag1 or not frag2 or frag1 == frag2:
                continue

            # 仕様: IDが大きい方の原子が含まれるフラグメントを動かす
            if rep_item1.atom_id > rep_item2.atom_id:
                ids_to_move = frag1
            else:
                ids_to_move = frag2

            # 3c: 移動計画を作成
            translation_vector = QPointF(-MOVE_DISTANCE, MOVE_DISTANCE)  # 左下方向へのベクトル
            move_operations.append((ids_to_move, translation_vector))

        # --- ステップ4: 計画された移動を一度に実行 ---
        if not move_operations:
            self.statusBar().showMessage("No actionable overlaps found.", 2000)
            return

        for group_ids, vector in move_operations:
            for atom_id in group_ids:
                item = self.data.atoms[atom_id]['item']
                new_pos = item.pos() + vector
                item.setPos(new_pos)
                self.data.atoms[atom_id]['pos'] = new_pos

        # --- ステップ5: 表示と状態を更新 ---
        for bond_data in self.data.bonds.values():
            item = bond_data.get('item') if bond_data else None
            if not item:
                continue
            try:
                if sip_isdeleted_safe(item):
                    continue
            except Exception:
                pass
            try:
                sc = None
                try:
                    sc = item.scene() if hasattr(item, 'scene') else None
                except Exception:
                    sc = None
                if sc is None:
                    continue
                try:
                    item.update_position()
                except Exception:
                    continue
            except Exception:
                continue
        
        # 重なり解消後に測定ラベルの位置を更新
        self.update_2d_measurement_labels()
        
        self.scene.update()
        self.push_undo_state()
        self.statusBar().showMessage("Resolved overlapping groups.", 2000)




    def adjust_molecule_positions_to_avoid_collisions(self, mol, frags):
        """
        複数分子の位置を調整して、衝突を回避する（バウンディングボックス最適化版）
        """
        if len(frags) <= 1:
            return
        
        conf = mol.GetConformer()
        pt = Chem.GetPeriodicTable()
        
        # --- 1. 各フラグメントの情報（原子インデックス、VDW半径）を事前計算 ---
        frag_info = []
        for frag_indices in frags:
            positions = []
            vdw_radii = []
            for idx in frag_indices:
                pos = conf.GetAtomPosition(idx)
                positions.append(np.array([pos.x, pos.y, pos.z]))
                
                atom = mol.GetAtomWithIdx(idx)
                # GetRvdw() はファンデルワールス半径を返す
                vdw_radii.append(pt.GetRvdw(atom.GetAtomicNum()))

            positions_np = np.array(positions)
            vdw_radii_np = np.array(vdw_radii)
            
            # このフラグメントで最大のVDW半径を計算（ボックスのマージンとして使用）
            max_vdw = np.max(vdw_radii_np) if len(vdw_radii_np) > 0 else 0.0

            frag_info.append({
                'indices': frag_indices,
                'centroid': np.mean(positions_np, axis=0),
                'positions_np': positions_np, # Numpy配列として保持
                'vdw_radii_np': vdw_radii_np,  # Numpy配列として保持
                'max_vdw_radius': max_vdw,
                'bbox_min': np.zeros(3), # 後で計算
                'bbox_max': np.zeros(3)  # 後で計算
            })
        
        # --- 2. 衝突判定のパラメータ ---
        collision_scale = 1.2  # VDW半径の120%
        max_iterations = 100
        moved = True
        iteration = 0
        
        while moved and iteration < max_iterations:
            moved = False
            iteration += 1
            
            # --- 3. フラグメントのバウンディングボックスを毎イテレーション更新 ---
            for i in range(len(frag_info)):
                # 現在の座標からボックスを再計算
                current_positions = []
                for idx in frag_info[i]['indices']:
                    pos = conf.GetAtomPosition(idx)
                    current_positions.append([pos.x, pos.y, pos.z])
                
                positions_np = np.array(current_positions)
                frag_info[i]['positions_np'] = positions_np # 座標情報を更新
                
                # VDW半径とスケールを考慮したマージンを計算
                # (最大VDW半径 * スケール) をマージンとして使う
                margin = frag_info[i]['max_vdw_radius'] * collision_scale
                
                frag_info[i]['bbox_min'] = np.min(positions_np, axis=0) - margin
                frag_info[i]['bbox_max'] = np.max(positions_np, axis=0) + margin

            # --- 4. 衝突判定ループ ---
            for i in range(len(frag_info)):
                for j in range(i + 1, len(frag_info)):
                    frag_i = frag_info[i]
                    frag_j = frag_info[j]
                    
                    # === バウンディングボックス判定 ===
                    # 2つのボックスが重なっているかチェック (AABB交差判定)
                    # X, Y, Zの各軸で重なりをチェック
                    overlap_x = (frag_i['bbox_min'][0] <= frag_j['bbox_max'][0] and frag_i['bbox_max'][0] >= frag_j['bbox_min'][0])
                    overlap_y = (frag_i['bbox_min'][1] <= frag_j['bbox_max'][1] and frag_i['bbox_max'][1] >= frag_j['bbox_min'][1])
                    overlap_z = (frag_i['bbox_min'][2] <= frag_j['bbox_max'][2] and frag_i['bbox_max'][2] >= frag_j['bbox_min'][2])
                    
                    # ボックスがX, Y, Zのいずれかの軸で離れている場合、原子間の詳細なチェックをスキップ
                    if not (overlap_x and overlap_y and overlap_z):
                        continue
                    # =================================

                    # ボックスが重なっている場合のみ、高コストな原子間の総当たりチェックを実行
                    total_push_vector = np.zeros(3)
                    collision_count = 0
                    
                    # 事前計算したNumpy配列を使用
                    positions_i = frag_i['positions_np']
                    positions_j = frag_j['positions_np']
                    vdw_i_all = frag_i['vdw_radii_np']
                    vdw_j_all = frag_j['vdw_radii_np']

                    for k, idx_i in enumerate(frag_i['indices']):
                        pos_i = positions_i[k]
                        vdw_i = vdw_i_all[k]
                        
                        for l, idx_j in enumerate(frag_j['indices']):
                            pos_j = positions_j[l]
                            vdw_j = vdw_j_all[l]
                            
                            distance_vec = pos_i - pos_j
                            distance_sq = np.dot(distance_vec, distance_vec) # 平方根を避けて高速化
                            
                            min_distance = (vdw_i + vdw_j) * collision_scale
                            min_distance_sq = min_distance * min_distance
                            
                            if distance_sq < min_distance_sq and distance_sq > 0.0001:
                                distance = np.sqrt(distance_sq)
                                push_direction = distance_vec / distance
                                push_magnitude = (min_distance - distance) / 2 # 押し出し量は半分ずつ
                                total_push_vector += push_direction * push_magnitude
                                collision_count += 1
                    
                    if collision_count > 0:
                        # 平均的な押し出しベクトルを適用
                        avg_push_vector = total_push_vector / collision_count
                        
                        # Conformerの座標を更新
                        for idx in frag_i['indices']:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos + avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())
                        
                        for idx in frag_j['indices']:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos - avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())
                        
                        moved = True
                        # (この移動により、このイテレーションで使う frag_info の座標キャッシュが古くなりますが、
                        #  次のイテレーションの最初でボックスと共に再計算されるため問題ありません)



    def _apply_chem_check_and_set_flags(self, mol, source_desc=None):
        """Central helper to apply chemical sanitization (or skip it) and set
        chem_check_tried / chem_check_failed flags consistently.

        When sanitization fails, a warning is shown and the Optimize 3D button
        is disabled. If the user setting 'skip_chemistry_checks' is True, no
        sanitization is attempted and both flags remain False.
        """
        try:
            self.chem_check_tried = False
            self.chem_check_failed = False
        except Exception:
            # Ensure attributes exist even if called very early
            self.chem_check_tried = False
            self.chem_check_failed = False

        if self.settings.get('skip_chemistry_checks', False):
            # User asked to skip chemistry checks entirely
            return

        try:
            Chem.SanitizeMol(mol)
            self.chem_check_tried = True
            self.chem_check_failed = False
        except Exception:
            # Mark that we tried sanitization and it failed
            self.chem_check_tried = True
            self.chem_check_failed = True
            try:
                desc = f" ({source_desc})" if source_desc else ''
                self.statusBar().showMessage(f"Molecule sanitization failed{desc}; file may be malformed.")
            except Exception:
                pass
            # Disable 3D optimization UI to prevent running on invalid molecules
            if hasattr(self, 'optimize_3d_button'):
                try:
                    self.optimize_3d_button.setEnabled(False)
                except Exception:
                    pass
        


    def _clear_xyz_flags(self, mol=None):
        """Clear XYZ-derived markers from a molecule (or current_mol) and
        reset UI flags accordingly.

        This is a best-effort cleanup to remove properties like
        _xyz_skip_checks and _xyz_atom_data that may have been attached when
        an XYZ file was previously loaded. After clearing molecule-level
        markers, the UI flag self.is_xyz_derived is set to False and the
        Optimize 3D button is re-evaluated (enabled unless chem_check_failed
        is True).
        """
        target = mol if mol is not None else getattr(self, 'current_mol', None)
        try:
            if target is not None:
                # Remove RDKit property if present
                try:
                    if hasattr(target, 'HasProp') and target.HasProp('_xyz_skip_checks'):
                        try:
                            target.ClearProp('_xyz_skip_checks')
                        except Exception:
                            try:
                                target.SetIntProp('_xyz_skip_checks', 0)
                            except Exception:
                                pass
                except Exception:
                    pass

                # Remove attribute-style markers if present
                try:
                    if hasattr(target, '_xyz_skip_checks'):
                        try:
                            delattr(target, '_xyz_skip_checks')
                        except Exception:
                            try:
                                del target._xyz_skip_checks
                            except Exception:
                                try:
                                    target._xyz_skip_checks = False
                                except Exception:
                                    pass
                except Exception:
                    pass

                try:
                    if hasattr(target, '_xyz_atom_data'):
                        try:
                            delattr(target, '_xyz_atom_data')
                        except Exception:
                            try:
                                del target._xyz_atom_data
                            except Exception:
                                try:
                                    target._xyz_atom_data = None
                                except Exception:
                                    pass
                except Exception:
                    pass

        except Exception:
            # best-effort only
            pass

        # Reset UI flags
        try:
            self.is_xyz_derived = False
        except Exception:
            pass

        # Enable Optimize 3D unless sanitization failed
        try:
            if hasattr(self, 'optimize_3d_button'):
                if getattr(self, 'chem_check_failed', False):
                    try:
                        self.optimize_3d_button.setEnabled(False)
                    except Exception:
                        pass
                else:
                    try:
                        self.optimize_3d_button.setEnabled(True)
                    except Exception:
                        pass
        except Exception:
            pass
            
