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
main_window_view_loaders.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowViewLoaders
"""


import os
import traceback


# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QFileDialog
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
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import VERSION

# --- クラス定義 ---
class MainWindowViewLoaders(object):
    """ main_window.py から分離された機能クラス """

    def load_xyz_for_3d_viewing(self, file_path=None):
        """XYZファイルを読み込んで3Dビューアで表示する"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Load 3D XYZ (View Only)", "", "XYZ Files (*.xyz);;All Files (*)")
            if not file_path:
                return

        try:
            mol = self.load_xyz_file(file_path)
            if mol is None:
                raise ValueError("Failed to create molecule from XYZ file.")
            if mol.GetNumConformers() == 0:
                raise ValueError("XYZ file has no 3D coordinates.")

            # 2Dエディタをクリア
            self.clear_2d_editor(push_to_undo=False)
            
            # 3D構造をセットして描画
            # Set the molecule. If bonds were determined (mol has bonds),
            # treat this the same as loading a MOL file: clear the XYZ-derived
            # flag and enable 3D optimization. Only mark as XYZ-derived and
            # disable 3D optimization when the molecule has no bond information.
            self.current_mol = mol

            # XYZファイル読み込み時はマッピングをクリア（2D構造がないため）
            self.atom_id_to_rdkit_idx_map = {}

            # If the loader marked the molecule as produced under skip_chemistry_checks,
            # always treat it as XYZ-derived and disable optimization. Otherwise
            # fall back to the existing behavior based on bond presence.
            skip_flag = False
            try:
                # Prefer RDKit int prop
                skip_flag = bool(self.current_mol.GetIntProp("_xyz_skip_checks"))
            except Exception:
                try:
                    skip_flag = bool(getattr(self.current_mol, '_xyz_skip_checks', False))
                except Exception:
                    skip_flag = False

            if skip_flag:
                self.is_xyz_derived = True
                if hasattr(self, 'optimize_3d_button'):
                    try:
                        self.optimize_3d_button.setEnabled(False)
                    except Exception:
                        pass
            else:
                try:
                    has_bonds = (self.current_mol.GetNumBonds() > 0)
                except Exception:
                    has_bonds = False

                if has_bonds:
                    self.is_xyz_derived = False
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            # Only enable optimize if the molecule is not considered XYZ-derived
                            if not getattr(self, 'is_xyz_derived', False):
                                self.optimize_3d_button.setEnabled(True)
                            else:
                                self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass
                else:
                    self.is_xyz_derived = True
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass

            self.draw_molecule_3d(self.current_mol)
            self.plotter.reset_camera()

            # UIを3Dビューアモードに設定
            self._enter_3d_viewer_ui_mode()

            # 3D関連機能を統一的に有効化
            self._enable_3d_features(True)
            
            # メニューテキストと状態を更新
            self.update_atom_id_menu_text()
            self.update_atom_id_menu_state()
            
            self.statusBar().showMessage(f"3D Viewer Mode: Loaded {os.path.basename(file_path)}")
            self.reset_undo_stack()
            # XYZファイル名をcurrent_file_pathにセットし、未保存状態はFalse
            self.current_file_path = file_path
            self.has_unsaved_changes = False
            self.update_window_title()

        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
            self.restore_ui_for_editing()
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid XYZ file: {e}")
            self.restore_ui_for_editing()
        except Exception as e:
            self.statusBar().showMessage(f"Error loading XYZ file: {e}")
            self.restore_ui_for_editing()
            
            traceback.print_exc()



    def save_3d_as_mol(self):
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        try:
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    name = os.path.splitext(base)[0]
                    default_name = f"{name}"
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save 3D MOL File", default_path, "MOL Files (*.mol);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.mol'):
                file_path += '.mol'

            mol_to_save = Chem.Mol(self.current_mol)

            if mol_to_save.HasProp("_2D"):
                mol_to_save.ClearProp("_2D")

            mol_block = Chem.MolToMolBlock(mol_to_save, includeStereo=True)
            lines = mol_block.split('\n')
            if len(lines) > 1 and 'RDKit' in lines[1]:
                lines[1] = '  MoleditPy Ver. ' + VERSION + '  3D'
            modified_mol_block = '\n'.join(lines)
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(modified_mol_block)
            self.statusBar().showMessage(f"3D data saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except UnicodeEncodeError as e:
            self.statusBar().showMessage(f"Text encoding error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving 3D MOL file: {e}")
            
            traceback.print_exc()



    def load_mol_file_for_3d_viewing(self, file_path=None):
        """MOL/SDFファイルを3Dビューアーで開く"""
        if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open MOL/SDF File", "", 
                "MOL/SDF Files (*.mol *.sdf);;All Files (*)"
            )
            if not file_path:
                return
        
        try:
            
            # Determine extension early and handle .mol specially by reading the
            # raw block and running it through fix_mol_block before parsing.
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ''

            if ext == '.sdf':
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)

            elif ext == '.mol':
                # Read the file contents and attempt to fix malformed counts lines
                with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)

                # If parsing the fixed block fails, fall back to RDKit's file reader
                # as a last resort (keeps behavior conservative).
                if mol is None:
                    try:
                        mol = Chem.MolFromMolFile(file_path, removeHs=False)
                    except Exception:
                        mol = None

                if mol is None:
                    self.statusBar().showMessage(f"Failed to load molecule from {file_path}")
                    return

            else:
                # Default: let RDKit try to read the file (most common case)
                if file_path.lower().endswith('.sdf'):
                    suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                    mol = next(suppl, None)
                else:
                    mol = Chem.MolFromMolFile(file_path, removeHs=False)

                if mol is None:
                    self.statusBar().showMessage(f"Failed to load molecule from {file_path}")
                    return

            # 3D座標がない場合は2Dから3D変換（最適化なし）
            if mol.GetNumConformers() == 0:
                self.statusBar().showMessage("No 3D coordinates found. Converting to 3D...")
                try:
                    try:
                        AllChem.EmbedMolecule(mol)
                        # 最適化は実行しない
                        # 3D変換直後にUndoスタックに積む
                        self.current_mol = mol
                        self.push_undo_state()
                    except Exception:
                        # If skipping chemistry checks, allow molecule to be displayed without 3D embedding
                        if self.settings.get('skip_chemistry_checks', False):
                            self.statusBar().showMessage("Warning: failed to generate 3D coordinates but skip_chemistry_checks is enabled; continuing.")
                            # Keep mol as-is (may lack conformer); downstream code checks for conformers
                        else:
                            raise
                except Exception:
                    self.statusBar().showMessage("Failed to generate 3D coordinates")
                    return
            
            # Clear XYZ markers on the newly loaded MOL/SDF so Optimize 3D is
            # correctly enabled when appropriate.
            try:
                self._clear_xyz_flags(mol)
            except Exception:
                pass

            # 3Dビューアーに表示
            # Centralized chemical/sanitization handling
            # Ensure the skip_chemistry_checks setting is respected and flags are set
            self._apply_chem_check_and_set_flags(mol, source_desc='MOL/SDF')

            self.current_mol = mol
            self.draw_molecule_3d(mol)
            
            # カメラをリセット
            self.plotter.reset_camera()
            
            # UIを3Dビューアーモードに設定
            self._enter_3d_viewer_ui_mode()
            
            # メニューテキストと状態を更新
            self.update_atom_id_menu_text()
            self.update_atom_id_menu_state()
            
            self.statusBar().showMessage(f"Loaded {file_path} in 3D viewer")
            
            self.reset_undo_stack()
            self.has_unsaved_changes = False  # ファイル読込直後は未変更扱い
            self.current_file_path = file_path
            self.update_window_title()
            

        except Exception as e:
            self.statusBar().showMessage(f"Error loading MOL/SDF file: {e}")
    
