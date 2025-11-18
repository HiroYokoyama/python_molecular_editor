#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_window_app_state.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowAppState
"""


import numpy as np
import copy
import os
import base64


# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QMessageBox
)



from PyQt6.QtCore import (
    Qt, QPointF, QDateTime
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
    from .constants import VERSION
    from .atom_item import AtomItem
    from .bond_item import BondItem
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import VERSION
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem


# --- クラス定義 ---
class MainWindowAppState(object):
    """ main_window.py から分離された機能クラス """

    def __init__(self):
        """
        クラスの初期化
        BoundFeature経由で呼ばれるため、'self' には MainWindow インスタンスが渡されます。
        """
        self.DEBUG_UNDO = False



    def get_current_state(self):
        atoms = {atom_id: {'symbol': data['symbol'],
                           'pos': (data['item'].pos().x(), data['item'].pos().y()),
                           'charge': data.get('charge', 0),
                           'radical': data.get('radical', 0)} 
                 for atom_id, data in self.data.atoms.items()}
        bonds = {key: {'order': data['order'], 'stereo': data.get('stereo', 0)} for key, data in self.data.bonds.items()}
        state = {'atoms': atoms, 'bonds': bonds, '_next_atom_id': self.data._next_atom_id}

        state['version'] = VERSION 
        
        if self.current_mol: state['mol_3d'] = self.current_mol.ToBinary()

        state['is_3d_viewer_mode'] = not self.is_2d_editable

        json_safe_constraints = []
        try:
            for const in self.constraints_3d:
                # (Type, (Idx...), Value, Force) -> [Type, [Idx...], Value, Force]
                if len(const) == 4:
                    json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                else:
                    # 後方互換性: 3要素の場合はデフォルトForceを追加
                    json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
        except Exception:
            pass # 失敗したら空リスト
        state['constraints_3d'] = json_safe_constraints
            
        return state



    def set_state_from_data(self, state_data):
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)
        
        loaded_data = copy.deepcopy(state_data)

        # ファイルのバージョンを取得（存在しない場合は '0.0.0' とする）
        file_version_str = loaded_data.get('version', '0.0.0')

        try:
            app_version_parts = tuple(map(int, VERSION.split('.')))
            file_version_parts = tuple(map(int, file_version_str.split('.')))

            # ファイルのバージョンがアプリケーションのバージョンより新しい場合に警告
            if file_version_parts > app_version_parts:
                QMessageBox.warning(
                    self,
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).\n\n"
                    f"Your current version is {VERSION}.\n\n"
                    "Some features may not load or work correctly."
                )
        except (ValueError, AttributeError):
            pass

        raw_atoms = loaded_data.get('atoms', {})
        raw_bonds = loaded_data.get('bonds', {})

        # 制約データの復元 (pmeraw)
        try:
            loaded_constraints = loaded_data.get("constraints_3d", [])
            # pmerawもJSON互換形式 [Type, [Idx...], Value, Force] で保存されている想定
            self.constraints_3d = []
            for const in loaded_constraints:
                if isinstance(const, list):
                    if len(const) == 4:
                        # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                        self.constraints_3d.append((const[0], tuple(const[1]), const[2], const[3]))
                    elif len(const) == 3:
                        # 後方互換性: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                        self.constraints_3d.append((const[0], tuple(const[1]), const[2], 1.0e5))
        except Exception:
            self.constraints_3d = [] # 読み込み失敗時はリセット

        for atom_id, data in raw_atoms.items():
            pos = QPointF(data['pos'][0], data['pos'][1])
            charge = data.get('charge', 0)
            radical = data.get('radical', 0)  # <-- ラジカル情報を取得
            # AtomItem生成時にradicalを渡す
            atom_item = AtomItem(atom_id, data['symbol'], pos, charge=charge, radical=radical)
            # self.data.atomsにもradical情報を格納する
            self.data.atoms[atom_id] = {'symbol': data['symbol'], 'pos': pos, 'item': atom_item, 'charge': charge, 'radical': radical}
            self.scene.addItem(atom_item)
        
        self.data._next_atom_id = loaded_data.get('_next_atom_id', max(self.data.atoms.keys()) + 1 if self.data.atoms else 0)

        for key_tuple, data in raw_bonds.items():
            id1, id2 = key_tuple
            if id1 in self.data.atoms and id2 in self.data.atoms:
                atom1_item = self.data.atoms[id1]['item']; atom2_item = self.data.atoms[id2]['item']
                bond_item = BondItem(atom1_item, atom2_item, data.get('order', 1), data.get('stereo', 0))
                self.data.bonds[key_tuple] = {'order': data.get('order', 1), 'stereo': data.get('stereo', 0), 'item': bond_item}
                atom1_item.bonds.append(bond_item); atom2_item.bonds.append(bond_item)
                self.scene.addItem(bond_item)

        for atom_data in self.data.atoms.values():
            if atom_data['item']: atom_data['item'].update_style()
        self.scene.update()

        if 'mol_3d' in loaded_data and loaded_data['mol_3d'] is not None:
            try:
                self.current_mol = Chem.Mol(loaded_data['mol_3d'])
                # デバッグ：3D構造が有効かチェック
                if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                    self.draw_molecule_3d(self.current_mol)
                    self.plotter.reset_camera()
                    # 3D関連機能を統一的に有効化
                    self._enable_3d_features(True)
                    
                    # 3D原子情報ホバー表示を再設定
                    self.setup_3d_hover()
                else:
                    # 無効な3D構造の場合
                    self.current_mol = None
                    self.plotter.clear()
                    # 3D関連機能を統一的に無効化
                    self._enable_3d_features(False)
            except Exception as e:
                self.statusBar().showMessage(f"Could not load 3D model from project: {e}")
                self.current_mol = None
                # 3D関連機能を統一的に無効化
                self._enable_3d_features(False)
        else:
            self.current_mol = None; self.plotter.clear(); self.analysis_action.setEnabled(False)
            self.optimize_3d_button.setEnabled(False)
            # 3D関連機能を統一的に無効化
            self._enable_3d_features(False)

        self.update_implicit_hydrogens()
        self.update_chiral_labels()

        if loaded_data.get('is_3d_viewer_mode', False):
            self._enter_3d_viewer_ui_mode()
            self.statusBar().showMessage("Project loaded in 3D Viewer Mode.")
        else:
            self.restore_ui_for_editing()
            # 3D分子がある場合は、2Dエディタモードでも3D編集機能を有効化
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                self._enable_3d_edit_actions(True)
        
        # undo/redo後に測定ラベルの位置を更新
        self.update_2d_measurement_labels()
        



    def push_undo_state(self):
        if self._is_restoring_state:
            return
            
        current_state_for_comparison = {
            'atoms': {k: (v['symbol'], v['item'].pos().x(), v['item'].pos().y(), v.get('charge', 0), v.get('radical', 0)) for k, v in self.data.atoms.items()},
            'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in self.data.bonds.items()},
            '_next_atom_id': self.data._next_atom_id,
            'mol_3d': self.current_mol.ToBinary() if self.current_mol else None
        }
        
        last_state_for_comparison = None
        if self.undo_stack:
            last_state = self.undo_stack[-1]
            last_atoms = last_state.get('atoms', {})
            last_bonds = last_state.get('bonds', {})
            last_state_for_comparison = {
                'atoms': {k: (v['symbol'], v['pos'][0], v['pos'][1], v.get('charge', 0), v.get('radical', 0)) for k, v in last_atoms.items()},
                'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in last_bonds.items()},
                '_next_atom_id': last_state.get('_next_atom_id'),
                'mol_3d': last_state.get('mol_3d', None)
            }

        if not last_state_for_comparison or current_state_for_comparison != last_state_for_comparison:
            # Deepcopy state to ensure saved states are immutable and not affected
            # by later modifications to objects referenced from the state.
            state = copy.deepcopy(self.get_current_state())
            self.undo_stack.append(state)
            if getattr(self, 'DEBUG_UNDO', False):
                try:
                    print(f"DEBUG_UNDO: push_undo_state -> new stack size: {len(self.undo_stack)}")
                except Exception:
                    pass
            self.redo_stack.clear()
            # 初期化完了後のみ変更があったことを記録
            if self.initialization_complete:
                self.has_unsaved_changes = True
                self.update_window_title()
        
        self.update_implicit_hydrogens()
        self.update_realtime_info()
        self.update_undo_redo_actions()



    def update_window_title(self):
        """ウィンドウタイトルを更新（保存状態を反映）"""
        base_title = f"MoleditPy Ver. {VERSION}"
        if self.current_file_path:
            filename = os.path.basename(self.current_file_path)
            title = f"{filename} - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        else:
            # Untitledファイルとして扱う
            title = f"Untitled - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        self.setWindowTitle(title)



    def check_unsaved_changes(self):
        """未保存の変更があるかチェックし、警告ダイアログを表示"""
        if not self.has_unsaved_changes:
            return True  # 保存済みまたは変更なし
        
        if not self.data.atoms and self.current_mol is None:
            return True  # 空のドキュメント
        
        reply = QMessageBox.question(
            self,
            "Unsaved Changes",
            "You have unsaved changes. Do you want to save them?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            # 拡張子がPMEPRJでなければ「名前を付けて保存」
            file_path = self.current_file_path
            if not file_path or not file_path.lower().endswith('.pmeprj'):
                self.save_project_as()
            else:
                self.save_project()
            return not self.has_unsaved_changes  # 保存に成功した場合のみTrueを返す
        elif reply == QMessageBox.StandardButton.No:
            return True  # 保存せずに続行
        else:
            return False  # キャンセル



    def reset_undo_stack(self):
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.push_undo_state()
        if getattr(self, 'DEBUG_UNDO', False):
            try:
                print(f"DEBUG_UNDO: reset_undo_stack -> undo={len(self.undo_stack)} redo={len(self.redo_stack)}")
            except Exception:
                pass



    def undo(self):
        if len(self.undo_stack) > 1:
            self.redo_stack.append(self.undo_stack.pop())
            state = self.undo_stack[-1]
            self._is_restoring_state = True
            try:
                self.set_state_from_data(state)
            finally:
                self._is_restoring_state = False

            
            # Undo後に3D構造の状態に基づいてメニューを再評価
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D構造がある場合は3D編集機能を有効化
                self._enable_3d_edit_actions(True)
            else:
                # 3D構造がない場合は3D編集機能を無効化
                self._enable_3d_edit_actions(False)
                    
        if getattr(self, 'DEBUG_UNDO', False):
            try:
                print(f"DEBUG_UNDO: undo -> undo_stack size: {len(self.undo_stack)}, redo_stack size: {len(self.redo_stack)}")
            except Exception:
                pass
        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus() 



    def redo(self):
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self._is_restoring_state = True
            try:
                self.set_state_from_data(state)
            finally:
                self._is_restoring_state = False
            
            # Redo後に3D構造の状態に基づいてメニューを再評価
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D構造がある場合は3D編集機能を有効化
                self._enable_3d_edit_actions(True)
            else:
                # 3D構造がない場合は3D編集機能を無効化
                self._enable_3d_edit_actions(False)
                    
        if getattr(self, 'DEBUG_UNDO', False):
            try:
                print(f"DEBUG_UNDO: redo -> undo_stack size: {len(self.undo_stack)}, redo_stack size: {len(self.redo_stack)}")
            except Exception:
                pass
        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus() 
        


    def update_undo_redo_actions(self):
        self.undo_action.setEnabled(len(self.undo_stack) > 1)
        self.redo_action.setEnabled(len(self.redo_stack) > 0)



    def update_realtime_info(self):
        """ステータスバーの右側に現在の分子情報を表示する"""
        if not self.data.atoms:
            self.formula_label.setText("")  # 原子がなければ右側のラベルをクリア
            return

        try:
            mol = self.data.to_rdkit_mol()
            if mol:
                # 水素原子を明示的に追加した分子オブジェクトを生成
                mol_with_hs = Chem.AddHs(mol)
                mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                # 水素を含む分子オブジェクトから原子数を取得
                num_atoms = mol_with_hs.GetNumAtoms()
                # 右側のラベルのテキストを更新
                self.formula_label.setText(f"Formula: {mol_formula}   |   Atoms: {num_atoms}")
        except Exception:
            # 計算に失敗してもアプリは継続
            self.formula_label.setText("Invalid structure")



    def create_json_data(self):
        """現在の状態をPMEJSON形式のデータに変換"""
        # 基本的なメタデータ
        json_data = {
            "format": "PME Project",
            "version": "1.0",
            "application": "MoleditPy",
            "application_version": VERSION,
            "created": str(QDateTime.currentDateTime().toString(Qt.DateFormat.ISODate)),
            "is_3d_viewer_mode": not self.is_2d_editable
        }
        
        # 2D構造データ
        if self.data.atoms:
            atoms_2d = []
            for atom_id, data in self.data.atoms.items():
                pos = data['item'].pos()
                atom_data = {
                    "id": atom_id,
                    "symbol": data['symbol'],
                    "x": pos.x(),
                    "y": pos.y(),
                    "charge": data.get('charge', 0),
                    "radical": data.get('radical', 0)
                }
                atoms_2d.append(atom_data)
            
            bonds_2d = []
            for (atom1_id, atom2_id), bond_data in self.data.bonds.items():
                bond_info = {
                    "atom1": atom1_id,
                    "atom2": atom2_id,
                    "order": bond_data['order'],
                    "stereo": bond_data.get('stereo', 0)
                }
                bonds_2d.append(bond_info)
            
            json_data["2d_structure"] = {
                "atoms": atoms_2d,
                "bonds": bonds_2d,
                "next_atom_id": self.data._next_atom_id
            }
        
        # 3D分子データ
        if self.current_mol and self.current_mol.GetNumConformers() > 0:
            try:
                # MOLデータをBase64エンコードで保存（バイナリデータの安全な保存）
                mol_binary = self.current_mol.ToBinary()
                mol_base64 = base64.b64encode(mol_binary).decode('ascii')
                
                # 3D座標を抽出
                atoms_3d = []
                if self.current_mol.GetNumConformers() > 0:
                    conf = self.current_mol.GetConformer()
                    for i in range(self.current_mol.GetNumAtoms()):
                        atom = self.current_mol.GetAtomWithIdx(i)
                        pos = conf.GetAtomPosition(i)

                        # Try to preserve original editor atom ID (if present) so it can be
                        # restored when loading PMEPRJ files. RDKit atom properties may
                        # contain _original_atom_id when the molecule was created from
                        # the editor's 2D structure.
                        original_id = None
                        try:
                            if atom.HasProp("_original_atom_id"):
                                original_id = atom.GetIntProp("_original_atom_id")
                        except Exception:
                            original_id = None

                        atom_3d = {
                            "index": i,
                            "symbol": atom.GetSymbol(),
                            "atomic_number": atom.GetAtomicNum(),
                            "x": pos.x,
                            "y": pos.y,
                            "z": pos.z,
                            "formal_charge": atom.GetFormalCharge(),
                            "num_explicit_hs": atom.GetNumExplicitHs(),
                            "num_implicit_hs": atom.GetNumImplicitHs(),
                            # include original editor atom id when available for round-trip
                            "original_id": original_id
                        }
                        atoms_3d.append(atom_3d)
                
                # 結合情報を抽出
                bonds_3d = []
                for bond in self.current_mol.GetBonds():
                    bond_3d = {
                        "atom1": bond.GetBeginAtomIdx(),
                        "atom2": bond.GetEndAtomIdx(),
                        "order": int(bond.GetBondType()),
                        "is_aromatic": bond.GetIsAromatic(),
                        "stereo": int(bond.GetStereo())
                    }
                    bonds_3d.append(bond_3d)
                
                # constraints_3dをJSON互換形式に変換
                json_safe_constraints = []
                try:
                    for const in self.constraints_3d:
                        if len(const) == 4:
                            json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                        else:
                            json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
                except Exception:
                    json_safe_constraints = []
                
                json_data["3d_structure"] = {
                    "mol_binary_base64": mol_base64,
                    "atoms": atoms_3d,
                    "bonds": bonds_3d,
                    "num_conformers": self.current_mol.GetNumConformers(),
                    "constraints_3d": json_safe_constraints
                }
                
                # 分子の基本情報
                json_data["molecular_info"] = {
                    "num_atoms": self.current_mol.GetNumAtoms(),
                    "num_bonds": self.current_mol.GetNumBonds(),
                    "molecular_weight": Descriptors.MolWt(self.current_mol),
                    "formula": rdMolDescriptors.CalcMolFormula(self.current_mol)
                }
                
                # SMILESとInChI（可能であれば）
                try:
                    json_data["identifiers"] = {
                        "smiles": Chem.MolToSmiles(self.current_mol),
                        "canonical_smiles": Chem.MolToSmiles(self.current_mol, canonical=True)
                    }
                    
                    # InChI生成を試行
                    try:
                        inchi = Chem.MolToInchi(self.current_mol)
                        inchi_key = Chem.MolToInchiKey(self.current_mol)
                        json_data["identifiers"]["inchi"] = inchi
                        json_data["identifiers"]["inchi_key"] = inchi_key
                    except:
                        pass  # InChI生成に失敗した場合は無視
                        
                except Exception as e:
                    print(f"Warning: Could not generate molecular identifiers: {e}")
                    
            except Exception as e:
                print(f"Warning: Could not process 3D molecular data: {e}")
        else:
            # 3D情報がない場合の記録
            json_data["3d_structure"] = None
            json_data["note"] = "No 3D structure available. Generate 3D coordinates first."

        # Record the last-successful optimization method (if any)
        # This is a convenience field so saved projects remember which
        # optimizer variant was last used (e.g. "MMFF94s", "MMFF94", "UFF").
        try:
            json_data["last_successful_optimization_method"] = getattr(self, 'last_successful_optimization_method', None)
        except Exception:
            json_data["last_successful_optimization_method"] = None
        
        return json_data



    def load_from_json_data(self, json_data):
        """JSONデータから状態を復元"""
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)
        self._enable_3d_edit_actions(False)
        self._enable_3d_features(False)

        # 3Dビューアーモードの設定
        is_3d_mode = json_data.get("is_3d_viewer_mode", False)
        # Restore last successful optimization method if present in file
        try:
            self.last_successful_optimization_method = json_data.get("last_successful_optimization_method", None)
        except Exception:
            self.last_successful_optimization_method = None


        # 2D構造データの復元
        if "2d_structure" in json_data:
            structure_2d = json_data["2d_structure"]
            atoms_2d = structure_2d.get("atoms", [])
            bonds_2d = structure_2d.get("bonds", [])

            # 原子の復元
            for atom_data in atoms_2d:
                atom_id = atom_data["id"]
                symbol = atom_data["symbol"]
                pos = QPointF(atom_data["x"], atom_data["y"])
                charge = atom_data.get("charge", 0)
                radical = atom_data.get("radical", 0)

                atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
                self.data.atoms[atom_id] = {
                    'symbol': symbol,
                    'pos': pos,
                    'item': atom_item,
                    'charge': charge,
                    'radical': radical
                }
                self.scene.addItem(atom_item)

            # next_atom_idの復元
            self.data._next_atom_id = structure_2d.get(
                "next_atom_id",
                max([atom["id"] for atom in atoms_2d]) + 1 if atoms_2d else 0
            )

            # 結合の復元
            for bond_data in bonds_2d:
                atom1_id = bond_data["atom1"]
                atom2_id = bond_data["atom2"]

                if atom1_id in self.data.atoms and atom2_id in self.data.atoms:
                    atom1_item = self.data.atoms[atom1_id]['item']
                    atom2_item = self.data.atoms[atom2_id]['item']

                    bond_order = bond_data["order"]
                    stereo = bond_data.get("stereo", 0)

                    bond_item = BondItem(atom1_item, atom2_item, bond_order, stereo=stereo)
                    # 原子の結合リストに追加（重要：炭素原子の可視性判定で使用）
                    atom1_item.bonds.append(bond_item)
                    atom2_item.bonds.append(bond_item)

                    self.data.bonds[(atom1_id, atom2_id)] = {
                        'order': bond_order,
                        'item': bond_item,
                        'stereo': stereo
                    }
                    self.scene.addItem(bond_item)

            # --- ここで全AtomItemのスタイルを更新（炭素原子の可視性を正しく反映） ---
            for atom in self.data.atoms.values():
                atom['item'].update_style()
        # 3D構造データの復元
        if "3d_structure" in json_data:
            structure_3d = json_data["3d_structure"]

            # 制約データの復元 (JSONはタプルをリストとして保存するので、タプルに再変換)
            try:
                loaded_constraints = structure_3d.get("constraints_3d", [])
                self.constraints_3d = []
                for const in loaded_constraints:
                    if isinstance(const, list):
                        if len(const) == 4:
                            # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                            self.constraints_3d.append((const[0], tuple(const[1]), const[2], const[3]))
                        elif len(const) == 3:
                            # 後方互換性: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                            self.constraints_3d.append((const[0], tuple(const[1]), const[2], 1.0e5))
            except Exception:
                self.constraints_3d = [] # 読み込み失敗時はリセット

            try:
                # バイナリデータの復元
                mol_base64 = structure_3d.get("mol_binary_base64")
                if mol_base64:
                    mol_binary = base64.b64decode(mol_base64.encode('ascii'))
                    self.current_mol = Chem.Mol(mol_binary)
                    if self.current_mol:
                        # 3D座標の設定
                        if self.current_mol.GetNumConformers() > 0:
                            conf = self.current_mol.GetConformer()
                            atoms_3d = structure_3d.get("atoms", [])
                            self.atom_positions_3d = np.zeros((len(atoms_3d), 3))
                            for atom_data in atoms_3d:
                                idx = atom_data["index"]
                                if idx < len(self.atom_positions_3d):
                                    self.atom_positions_3d[idx] = [
                                        atom_data["x"], 
                                        atom_data["y"], 
                                        atom_data["z"]
                                    ]
                                # Restore original editor atom id into RDKit atom property
                                try:
                                    original_id = atom_data.get("original_id", None)
                                    if original_id is not None and idx < self.current_mol.GetNumAtoms():
                                        rd_atom = self.current_mol.GetAtomWithIdx(idx)
                                        # set as int prop so other code expecting _original_atom_id works
                                        rd_atom.SetIntProp("_original_atom_id", int(original_id))
                                except Exception:
                                    pass
                            # Build mapping from original 2D atom IDs to RDKit indices so
                            # 3D picks can be synchronized back to 2D AtomItems.
                            try:
                                self.create_atom_id_mapping()
                                # update menu and UI states that depend on original IDs
                                try:
                                    self.update_atom_id_menu_text()
                                    self.update_atom_id_menu_state()
                                except Exception:
                                    pass
                            except Exception:
                                # non-fatal if mapping creation fails
                                pass

                        # 3D分子があれば必ず3D表示
                        self.draw_molecule_3d(self.current_mol)
                        # ViewerモードならUIも切り替え
                        if is_3d_mode:
                            self._enter_3d_viewer_ui_mode()
                        else:
                            self.is_2d_editable = True
                        self.plotter.reset_camera()

                        # 成功的に3D分子が復元されたので、3D関連UIを有効にする
                        try:
                            self._enable_3d_edit_actions(True)
                            self._enable_3d_features(True)
                        except Exception:
                            pass
                            
            except Exception as e:
                print(f"Warning: Could not restore 3D molecular data: {e}")
                self.current_mol = None

