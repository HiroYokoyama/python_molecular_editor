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
main_window_molecular_parsers.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowMolecularParsers
"""


import io
import os
import contextlib
import traceback
import logging


# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolTransforms, rdGeometry
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout,
    QPushButton, QDialog, QFileDialog, QLabel, QLineEdit, QInputDialog, QDialogButtonBox
)



from PyQt6.QtCore import (
    QPointF, QTimer
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
        logging.warning("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
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
class MainWindowMolecularParsers(object):
    """ main_window.py から分離された機能クラス """


    def load_mol_file(self, file_path=None):
        if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Import MOL File", "", "Chemical Files (*.mol *.sdf);;All Files (*)")
            if not file_path: 
                return

        try:
            self.dragged_atom_info = None
            # If this is a single-record .mol file, read & fix the counts line
            # before parsing. For multi-record .sdf files, keep using SDMolSupplier.
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ''
            if ext == '.mol':
                # Read file text, fix CTAB counts line if needed, then parse
                with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
                if mol is None:
                    raise ValueError("Failed to read molecule from .mol file after fixing counts line.")
            else:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
                if mol is None:
                    raise ValueError("Failed to read molecule from file.")

            Chem.Kekulize(mol)

            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)
            
            # 1. 座標がなければ2D座標を生成する
            if mol.GetNumConformers() == 0: 
                AllChem.Compute2DCoords(mol)
            
            # 2. 座標の有無にかかわらず、常に立体化学を割り当て、2D表示用にくさび結合を設定する
            # これにより、3D座標を持つMOLファイルからでも正しく2Dの立体表現が生成される
            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            conf = mol.GetConformer()

            SCALE_FACTOR = 50.0
            
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())

            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if positions:
                mol_center_x = sum(p.x for p in positions) / len(positions)
                mol_center_y = sum(p.y for p in positions) / len(positions)
            else:
                mol_center_x, mol_center_y = 0.0, 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id
                        
            for bond in mol.GetBonds():
                b_idx,e_idx=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble(); b_dir = bond.GetBondDir()
                stereo = 0
                # Check for single bond Wedge/Dash
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2
                # ADDED: Check for double bond E/Z stereochemistry
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3 # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4 # E

                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item,a2_item=self.data.atoms[a1_id]['item'],self.data.atoms[a2_id]['item']

                self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            self.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.reset_undo_stack()
            # NEWファイル扱い: ファイルパスをクリアし未保存状態はFalse（変更なければ保存警告なし）
            self.current_file_path = file_path
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid MOL file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading file: {e}")
            
            traceback.print_exc()
    


    def load_xyz_file(self, file_path):
        """XYZファイルを読み込んでRDKitのMolオブジェクトを作成する"""
        
        if not self.check_unsaved_changes():
            return  # ユーザーがキャンセルした場合は何もしない

        try:
            # We will attempt one silent load with default charge=0 (no dialog).
            # If RDKit emits chemistry warnings (for example "Explicit valence ..."),
            # prompt the user once for an overall charge and retry. Only one retry is allowed.


            # Helper: prompt for charge once when needed
            # Returns a tuple: (charge_value_or_0, accepted:bool, skip_chemistry:bool)
            def prompt_for_charge():
                try:
                    # Create a custom dialog so we can provide a "Skip chemistry" button
                    dialog = QDialog(self)
                    dialog.setWindowTitle("Import XYZ Charge")
                    layout = QVBoxLayout(dialog)

                    label = QLabel("Enter total molecular charge:")
                    line_edit = QLineEdit(dialog)
                    line_edit.setText("")

                    # Standard OK/Cancel buttons
                    btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog)

                    # Additional Skip chemistry button
                    skip_btn = QPushButton("Skip chemistry", dialog)

                    # Horizontal layout for buttons
                    hl = QHBoxLayout()
                    hl.addWidget(btn_box)
                    hl.addWidget(skip_btn)

                    layout.addWidget(label)
                    layout.addWidget(line_edit)
                    layout.addLayout(hl)

                    result = {"accepted": False, "skip": False}

                    def on_ok():
                        result["accepted"] = True
                        dialog.accept()

                    def on_cancel():
                        dialog.reject()

                    def on_skip():
                        # Mark skip and accept so caller can proceed with skip behavior
                        result["skip"] = True
                        dialog.accept()

                    try:
                        btn_box.button(QDialogButtonBox.Ok).clicked.connect(on_ok)
                        btn_box.button(QDialogButtonBox.Cancel).clicked.connect(on_cancel)
                    except Exception:
                        # Fallback if button lookup fails
                        btn_box.accepted.connect(on_ok)
                        btn_box.rejected.connect(on_cancel)

                    skip_btn.clicked.connect(on_skip)

                    # Execute dialog modally
                    if dialog.exec_() != QDialog.Accepted:
                        return None, False, False

                    if result["skip"]:
                        # User chose to skip chemistry checks; return skip flag
                        return 0, True, True

                    if not result["accepted"]:
                        return None, False, False

                    charge_text = line_edit.text()
                except Exception:
                    # On any dialog creation error, fall back to simple input dialog
                    try:
                        charge_text, ok = QInputDialog.getText(self, "Import XYZ Charge", "Enter total molecular charge:", text="0")
                    except Exception:
                        return 0, True, False
                    if not ok:
                        return None, False, False
                    try:
                        return int(str(charge_text).strip()), True, False
                    except Exception:
                        try:
                            return int(float(str(charge_text).strip())), True, False
                        except Exception:
                            return 0, True, False

                if charge_text is None:
                    return None, False, False

                try:
                    return int(str(charge_text).strip()), True, False
                except Exception:
                    try:
                        return int(float(str(charge_text).strip())), True, False
                    except Exception:
                        return 0, True, False

            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 空行とコメント行を除去（但し、先頭2行は保持）
            non_empty_lines = []
            for i, line in enumerate(lines):
                stripped = line.strip()
                if i < 2:  # 最初の2行は原子数とコメント行なので保持
                    non_empty_lines.append(stripped)
                elif stripped and not stripped.startswith('#'):  # 空行とコメント行をスキップ
                    non_empty_lines.append(stripped)
            
            if len(non_empty_lines) < 2:
                raise ValueError("XYZ file format error: too few lines")
            
            # 原子数を読み取り
            try:
                num_atoms = int(non_empty_lines[0])
            except ValueError:
                raise ValueError("XYZ file format error: invalid atom count")
            
            if num_atoms <= 0:
                raise ValueError("XYZ file format error: atom count must be positive")
            
            # コメント行（2行目）
            comment = non_empty_lines[1] if len(non_empty_lines) > 1 else ""
            
            # 原子データを読み取り
            atoms_data = []
            data_lines = non_empty_lines[2:]
            
            if len(data_lines) < num_atoms:
                raise ValueError(f"XYZ file format error: expected {num_atoms} atom lines, found {len(data_lines)}")
            
            for i, line in enumerate(data_lines[:num_atoms]):
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"XYZ file format error: invalid atom data at line {i+3}")
                
                symbol = parts[0].strip()
                
                # 元素記号の妥当性をチェック
                try:
                    # RDKitで認識される元素かどうかをチェック
                    test_atom = Chem.Atom(symbol)
                except Exception:
                    # 認識されない場合、最初の文字を大文字にして再試行
                    symbol = symbol.capitalize()
                    try:
                        test_atom = Chem.Atom(symbol)
                    except Exception:
                        # If user requested to skip chemistry checks, coerce unknown symbols to C
                        if self.settings.get('skip_chemistry_checks', False):
                            symbol = 'C'
                        else:
                            raise ValueError(f"Unrecognized element symbol: {parts[0]} at line {i+3}")
                
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError:
                    raise ValueError(f"XYZ file format error: invalid coordinates at line {i+3}")
                
                atoms_data.append((symbol, x, y, z))
            
            if len(atoms_data) == 0:
                raise ValueError("XYZ file format error: no atoms found")
            
            # RDKitのMolオブジェクトを作成
            mol = Chem.RWMol()
            
            # 原子を追加
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                # XYZファイルでの原子のUniqueID（0ベースのインデックス）を保存
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
            
            # 3D座標を設定
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)
            # If user requested to skip chemistry checks, bypass RDKit's
            # DetermineBonds/sanitization flow entirely and use only the
            # distance-based bond estimation. Treat the resulting molecule
            # as "XYZ-derived" (disable 3D optimization) and return it.
            try:
                skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
            except Exception:
                skip_checks = False

            if skip_checks:
                used_rd_determine = False
                try:
                    # Use the conservative distance-based heuristic to add bonds
                    self.estimate_bonds_from_distances(mol)
                except Exception:
                    # Non-fatal: continue even if distance-based estimation fails
                    pass

                # Finalize and return a plain Mol object
                try:
                    candidate_mol = mol.GetMol()
                except Exception:
                    try:
                        candidate_mol = Chem.Mol(mol)
                    except Exception:
                        candidate_mol = None

                if candidate_mol is None:
                    raise ValueError("Failed to create valid molecule object when skip_chemistry_checks=True")

                # Attach a default charge property
                try:
                    candidate_mol.SetIntProp("_xyz_charge", 0)
                except Exception:
                    try:
                        candidate_mol._xyz_charge = 0
                    except Exception:
                        pass

                # Mark that this molecule was produced via the skip-chemistry path
                try:
                    candidate_mol.SetIntProp("_xyz_skip_checks", 1)
                except Exception:
                    try:
                        candidate_mol._xyz_skip_checks = True
                    except Exception:
                        pass

                # Set UI flags consistently: mark as XYZ-derived and disable optimize
                try:
                    self.current_mol = candidate_mol
                    self.is_xyz_derived = True
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass
                except Exception:
                    pass

                # Store atom data for later analysis and return
                candidate_mol._xyz_atom_data = atoms_data
                return candidate_mol
            # We'll attempt silently first with charge=0 and only prompt the user
            # for a charge when the RDKit processing block fails (raises an
            # exception). If the user provides a charge, retry; allow repeated
            # prompts until the user cancels. This preserves the previous
            # fallback behaviors (skip_chemistry_checks, distance-based bond
            # estimation) and property attachments.
            used_rd_determine = False
            final_mol = None

            # First, try silently with charge=0. If that raises an exception we
            # will enter a loop prompting the user for a charge and retrying as
            # long as the user provides values. If the user cancels, return None.
            def _process_with_charge(charge_val):
                """Inner helper: attempt to build/finalize molecule with given charge.

                Returns the finalized RDKit Mol on success. May raise exceptions
                which will be propagated to the caller.
                """
                nonlocal used_rd_determine
                # Capture RDKit stderr while we run the processing to avoid
                # spamming the console. We won't treat warnings specially here;
                # only exceptions will trigger a prompt/retry. We also want to
                # distinguish failures originating from DetermineBonds so the
                # outer logic can decide whether to prompt the user repeatedly
                # for different charge values.
                buf = io.StringIO()
                determine_failed = False
                with contextlib.redirect_stderr(buf):
                    # Try DetermineBonds if available
                    try:
                        from rdkit.Chem import rdDetermineBonds
                        try:
                            try:
                                mol_candidate = Chem.RWMol(Chem.Mol(mol))
                            except Exception:
                                mol_candidate = Chem.RWMol(mol)

                            # This call may raise. If it does, mark determine_failed
                            # so the caller can prompt for a different charge.
                            rdDetermineBonds.DetermineBonds(mol_candidate, charge=charge_val)
                            mol_to_finalize = mol_candidate
                            used_rd_determine = True
                        except Exception:
                            # DetermineBonds failed for this charge value. We
                            # should allow the caller to prompt for another
                            # charge (or cancel). Mark the flag and re-raise a
                            # dedicated exception to be handled by the outer
                            # loop.
                            determine_failed = True
                            used_rd_determine = False
                            mol_to_finalize = mol
                            # Raise a sentinel exception to indicate DetermineBonds failure
                            raise RuntimeError("DetermineBondsFailed")
                    except RuntimeError:
                        # Propagate our sentinel so outer code can catch it.
                        raise
                    except Exception:
                        # rdDetermineBonds not available or import failed; use
                        # distance-based fallback below.
                        used_rd_determine = False
                        mol_to_finalize = mol

                    if not used_rd_determine:
                        # distance-based fallback
                        self.estimate_bonds_from_distances(mol_to_finalize)

                    # Finalize molecule
                    try:
                        candidate_mol = mol_to_finalize.GetMol()
                    except Exception:
                        candidate_mol = None

                    if candidate_mol is None:
                        # Try salvage path
                        try:
                            candidate_mol = mol.GetMol()
                        except Exception:
                            candidate_mol = None

                    if candidate_mol is None:
                        raise ValueError("Failed to create valid molecule object")

                    # Attach charge property if possible
                    try:
                        try:
                            candidate_mol.SetIntProp("_xyz_charge", int(charge_val))
                        except Exception:
                            try:
                                candidate_mol._xyz_charge = int(charge_val)
                            except Exception:
                                pass
                    except Exception:
                        pass

                    # Preserve whether the user requested skip_chemistry_checks
                    try:
                        if bool(self.settings.get('skip_chemistry_checks', False)):
                            try:
                                candidate_mol.SetIntProp("_xyz_skip_checks", 1)
                            except Exception:
                                try:
                                    candidate_mol._xyz_skip_checks = True
                                except Exception:
                                    pass
                    except Exception:
                        pass

                    # Run chemistry checks which may emit warnings to stderr
                    self._apply_chem_check_and_set_flags(candidate_mol, source_desc='XYZ')

                    # Accept the candidate
                    return candidate_mol

            # Decide whether to silently try charge=0 first, or prompt user first.
            always_ask = bool(self.settings.get('always_ask_charge', False))

            try:
                if not always_ask:
                    # Silent first attempt (existing behavior)
                    try:
                        final_mol = _process_with_charge(0)
                    except RuntimeError:
                        # DetermineBonds explicitly failed for charge=0. In this
                        # situation, repeatedly prompt the user for charges until
                        # DetermineBonds succeeds or the user cancels.
                        while True:
                            charge_val, ok, skip_flag = prompt_for_charge()
                            if not ok:
                                # user cancelled the prompt -> abort
                                return None
                            if skip_flag:
                                # User selected Skip chemistry: attempt distance-based salvage
                                try:
                                    self.estimate_bonds_from_distances(mol)
                                except Exception:
                                    pass
                                salvaged = None
                                try:
                                    salvaged = mol.GetMol()
                                except Exception:
                                    salvaged = None

                                if salvaged is not None:
                                    try:
                                        salvaged.SetIntProp("_xyz_skip_checks", 1)
                                    except Exception:
                                        try:
                                            salvaged._xyz_skip_checks = True
                                        except Exception:
                                            pass
                                    final_mol = salvaged
                                    break
                                else:
                                    # Could not salvage; abort
                                    try:
                                        self.statusBar().showMessage("Skip chemistry selected but failed to create salvaged molecule.")
                                    except Exception:
                                        pass
                                    return None

                            try:
                                final_mol = _process_with_charge(charge_val)
                                # success -> break out of prompt loop
                                break
                            except RuntimeError:
                                # DetermineBonds still failing for this charge -> loop again
                                try:
                                    self.statusBar().showMessage("DetermineBonds failed for that charge; please try a different total charge or cancel.")
                                except Exception:
                                    pass
                                continue
                            except Exception as e_prompt:
                                # Some other failure occurred after DetermineBonds or in
                                # finalization. If skip_chemistry_checks is enabled we
                                # try the salvaged mol once; otherwise prompt again.
                                try:
                                    skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
                                except Exception:
                                    skip_checks = False

                                salvaged = None
                                try:
                                    salvaged = mol.GetMol()
                                except Exception:
                                    salvaged = None

                                if skip_checks and salvaged is not None:
                                    final_mol = salvaged
                                    # mark salvaged molecule as produced under skip_checks
                                    try:
                                        final_mol.SetIntProp("_xyz_skip_checks", 1)
                                    except Exception:
                                        try:
                                            final_mol._xyz_skip_checks = True
                                        except Exception:
                                            pass
                                    break
                                else:
                                    try:
                                        self.statusBar().showMessage(f"Retry failed: {e_prompt}")
                                    except Exception:
                                        pass
                                    # Continue prompting
                                    continue
                else:
                    # User has requested to always be asked for charge — prompt before any silent try
                    while True:
                        charge_val, ok, skip_flag = prompt_for_charge()
                        if not ok:
                            # user cancelled the prompt -> abort
                            return None
                        if skip_flag:
                            # User selected Skip chemistry: attempt distance-based salvage
                            try:
                                self.estimate_bonds_from_distances(mol)
                            except Exception:
                                pass
                            salvaged = None
                            try:
                                salvaged = mol.GetMol()
                            except Exception:
                                salvaged = None
    
                            if salvaged is not None:
                                try:
                                    salvaged.SetIntProp("_xyz_skip_checks", 1)
                                except Exception:
                                    try:
                                        salvaged._xyz_skip_checks = True
                                    except Exception:
                                        pass
                                final_mol = salvaged
                                break
                            else:
                                try:
                                    self.statusBar().showMessage("Skip chemistry selected but failed to create salvaged molecule.")
                                except Exception:
                                    pass
                                return None
    
                        try:
                            final_mol = _process_with_charge(charge_val)
                            # success -> break out of prompt loop
                            break
                        except RuntimeError:
                            # DetermineBonds still failing for this charge -> loop again
                            try:
                                self.statusBar().showMessage("DetermineBonds failed for that charge; please try a different total charge or cancel.")
                            except Exception:
                                pass
                            continue
                        except Exception as e_prompt:
                            try:
                                skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
                            except Exception:
                                skip_checks = False
    
                            salvaged = None
                            try:
                                salvaged = mol.GetMol()
                            except Exception:
                                salvaged = None
    
                            if skip_checks and salvaged is not None:
                                final_mol = salvaged
                                try:
                                    final_mol.SetIntProp("_xyz_skip_checks", 1)
                                except Exception:
                                    try:
                                        final_mol._xyz_skip_checks = True
                                    except Exception:
                                        pass
                                break
                            else:
                                try:
                                    self.statusBar().showMessage(f"Retry failed: {e_prompt}")
                                except Exception:
                                    pass
                                continue
                
            except Exception:
                # If the silent attempt failed for reasons other than
                # DetermineBonds failing (e.g., finalization errors), fall
                # back to salvaging or prompting depending on settings.
                salvaged = None
                try:
                    salvaged = mol.GetMol()
                except Exception:
                    salvaged = None

                try:
                    skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
                except Exception:
                    skip_checks = False

                if skip_checks and salvaged is not None:
                    final_mol = salvaged
                else:
                    # Repeatedly prompt until the user cancels or processing
                    # succeeds.
                    while True:
                        charge_val, ok, skip_flag = prompt_for_charge()
                        if not ok:
                            # user cancelled the prompt -> abort
                            return None
                        if skip_flag:
                            # User selected Skip chemistry: attempt distance-based salvage
                            try:
                                self.estimate_bonds_from_distances(mol)
                            except Exception:
                                pass
                            salvaged = None
                            try:
                                salvaged = mol.GetMol()
                            except Exception:
                                salvaged = None

                            if salvaged is not None:
                                try:
                                    salvaged.SetIntProp("_xyz_skip_checks", 1)
                                except Exception:
                                    try:
                                        salvaged._xyz_skip_checks = True
                                    except Exception:
                                        pass
                                final_mol = salvaged
                                break
                            else:
                                try:
                                    self.statusBar().showMessage("Skip chemistry selected but failed to create salvaged molecule.")
                                except Exception:
                                    pass
                                return None

                        try:
                            final_mol = _process_with_charge(charge_val)
                            # success -> break out of prompt loop
                            break
                        except RuntimeError:
                            # DetermineBonds failed for this charge -> let the
                            # user try another
                            try:
                                self.statusBar().showMessage("DetermineBonds failed for that charge; please try a different total charge or cancel.")
                            except Exception:
                                pass
                            continue
                        except Exception as e_prompt:
                            try:
                                self.statusBar().showMessage(f"Retry failed: {e_prompt}")
                            except Exception:
                                pass
                            continue

            # If we have a finalized molecule, apply the same UI flags and return
            if final_mol is not None:
                mol = final_mol
                try:
                    self.current_mol = mol

                    self.is_xyz_derived = not used_rd_determine
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            has_bonds = mol.GetNumBonds() > 0
                            # Respect the XYZ-derived flag: if the molecule is XYZ-derived,
                            # keep Optimize disabled regardless of bond detection.
                            if getattr(self, 'is_xyz_derived', False):
                                self.optimize_3d_button.setEnabled(False)
                            else:
                                self.optimize_3d_button.setEnabled(bool(has_bonds))
                        except Exception:
                            pass
                except Exception:
                    pass

                # Store original atom data for analysis
                mol._xyz_atom_data = atoms_data
                return mol
            
            # 元のXYZ原子データを分子オブジェクトに保存（分析用）
            mol._xyz_atom_data = atoms_data
            
            return mol
            
        except (OSError, IOError) as e:
            raise ValueError(f"File I/O error: {e}")
        except Exception as e:
            if "XYZ file format error" in str(e) or "Unrecognized element" in str(e):
                raise e
            else:
                raise ValueError(f"Error parsing XYZ file: {e}")



    def estimate_bonds_from_distances(self, mol):
        """原子間距離に基づいて結合を推定する"""
        
        # 一般的な共有結合半径（Ångström）- より正確な値
        covalent_radii = {
            'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76,
            'N': 0.75, 'O': 0.73, 'F': 0.71, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41,
            'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
            'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.14, 'Kr': 1.16,
            'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
            'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.33, 'Xe': 1.40
        }
        
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        # 追加された結合をトラッキング
        bonds_added = []
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)
                
                # 原子間距離を計算
                distance = rdMolTransforms.GetBondLength(conf, i, j)
                
                # 期待される結合距離を計算
                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()
                
                radius_i = covalent_radii.get(symbol_i, 1.0)  # デフォルト半径
                radius_j = covalent_radii.get(symbol_j, 1.0)
                
                expected_bond_length = radius_i + radius_j
                
                # 結合タイプによる許容範囲を調整
                # 水素結合は通常の共有結合より短い
                if symbol_i == 'H' or symbol_j == 'H':
                    tolerance_factor = 1.2  # 水素は結合が短くなりがち
                else:
                    tolerance_factor = 1.3  # 他の原子は少し余裕を持たせる
                
                max_bond_length = expected_bond_length * tolerance_factor
                min_bond_length = expected_bond_length * 0.5  # 最小距離も設定
                
                # 距離が期待値の範囲内なら結合を追加
                if min_bond_length <= distance <= max_bond_length:
                    try:
                        mol.AddBond(i, j, Chem.BondType.SINGLE)
                        bonds_added.append((i, j, distance))
                    except Exception:
                        # 既に結合が存在する場合はスキップ
                        pass
        
        # デバッグ情報（オプション）
        # Added bonds based on distance analysis
        
        return len(bonds_added)



    def save_as_mol(self):
        try:
            mol_block = self.data.to_mol_block()
            if not mol_block: 
                self.statusBar().showMessage("Error: No 2D data to save.") 
                return
                
            lines = mol_block.split('\n')
            if len(lines) > 1 and 'RDKit' in lines[1]:
                lines[1] = '  MoleditPy Ver. ' + VERSION + '  2D'
            modified_mol_block = '\n'.join(lines)
            
            # default filename: based on current_file_path, append -2d for 2D mol
            default_name = "untitled-2d"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    name = os.path.splitext(base)[0]
                    default_name = f"{name}-2d"
            except Exception:
                default_name = "untitled-2d"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save 2D MOL File", default_path, "MOL Files (*.mol);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.mol'): 
                file_path += '.mol'
                
            with open(file_path, 'w', encoding='utf-8') as f: 
                f.write(modified_mol_block)
            self.statusBar().showMessage(f"2D data saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except UnicodeEncodeError as e:
            self.statusBar().showMessage(f"Text encoding error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving file: {e}")
            
            traceback.print_exc()
            


    def save_as_xyz(self):
        if not self.current_mol: self.statusBar().showMessage("Error: Please generate a 3D structure first."); return
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

        file_path,_=QFileDialog.getSaveFileName(self,"Save 3D XYZ File",default_path,"XYZ Files (*.xyz);;All Files (*)")
        if file_path:
            if not file_path.lower().endswith('.xyz'): file_path += '.xyz'
            try:
                conf=self.current_mol.GetConformer(); num_atoms=self.current_mol.GetNumAtoms()
                xyz_lines=[str(num_atoms)]
                # 電荷と多重度を計算
                try:
                    charge = Chem.GetFormalCharge(self.current_mol)
                except Exception:
                    charge = 0 # 取得失敗時は0
                
                try:
                    # 全原子のラジカル電子の合計を取得
                    num_radicals = Descriptors.NumRadicalElectrons(self.current_mol)
                    # スピン多重度を計算 (M = N + 1, N=ラジカル電子数)
                    multiplicity = num_radicals + 1
                except Exception:
                    multiplicity = 1 # 取得失敗時は 1 (singlet)

                smiles=Chem.MolToSmiles(Chem.RemoveHs(self.current_mol))
                xyz_lines.append(f"chrg = {charge}  mult = {multiplicity} | Generated by MoleditPy Ver. {VERSION}")
                for i in range(num_atoms):
                    pos=conf.GetAtomPosition(i); symbol=self.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
                with open(file_path,'w') as f: f.write("\n".join(xyz_lines) + "\n")
                self.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")


    def fix_mol_counts_line(self, line: str) -> str:
        """
        Check and fix the CTAB counts line in a MOL file.
        If the line already contains 'V3000' or 'V2000' it is left unchanged.
        Otherwise the line is treated as V2000 and the proper 39-character
        format (33 chars of counts + ' V2000') is returned.
        """
        # If already V3000 or V2000, leave as-is
        if 'V3000' in line or 'V2000' in line:
            return line

        # Prepare prefix (first 33 characters for the 11 * I3 fields)
        prefix = line.rstrip().ljust(33)[0:33]
        version_str = ' V2000'
        return prefix + version_str

    def fix_mol_block(self, mol_block: str) -> str:
        """
        Given an entire MOL block as a string, ensure the 4th line (CTAB counts
        line) is valid. If the file has fewer than 4 lines, return as-is.
        """
        lines = mol_block.splitlines()
        if len(lines) < 4:
            # Not a valid MOL block — return unchanged
            return mol_block

        counts_line = lines[3]
        fixed_counts_line = self.fix_mol_counts_line(counts_line)
        lines[3] = fixed_counts_line
        return "\n".join(lines)

