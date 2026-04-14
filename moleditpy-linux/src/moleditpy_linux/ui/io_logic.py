#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import os
import json
import pickle
from typing import Any, Optional

from PyQt6.QtCore import QPointF, QTimer
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QDialogButtonBox,
    QPushButton,
    QFileDialog,
    QMessageBox,
)

from rdkit import Chem
from rdkit.Chem import AllChem, rdGeometry, rdMolTransforms, Descriptors

try:
    # package relative imports
    from ..utils.constants import COVALENT_RADII, VERSION
except (ImportError, ValueError, AttributeError):
    # Fallback to absolute imports
    try:
        from moleditpy_linux.utils.constants import COVALENT_RADII, VERSION
    except ImportError:
        COVALENT_RADII = {}
        VERSION = "unknown"


class IOManager:
    """Independent manager for IO actions (Load/Save), ported from Mixins."""

    def __init__(self, host: Any = None) -> None:
        if host is not None:
            self.host = host

    def _get_default_basename(self) -> str:
        """Helper to get a default filename base from the current file path."""
        try:
            if (
                hasattr(self.host.init_manager, "current_file_path")
                and self.host.init_manager.current_file_path
            ):
                base = os.path.basename(self.host.init_manager.current_file_path)
                name = os.path.splitext(base)[0]
                if name:
                    return name
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass
        return "untitled"

    def _get_default_path(self, suffix: str = "") -> str:
        """Get the full default path (dir + basename + suffix) based on the current file."""
        basename = self._get_default_basename() + suffix
        try:
            cur_path = self.host.init_manager.current_file_path
            if cur_path:
                return os.path.join(os.path.dirname(cur_path), basename)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass
        return basename

    def fix_mol_counts_line(self, line: str) -> str:
        if "V3000" in line or "V2000" in line:
            return line
        prefix = line.rstrip().ljust(33)[0:33]
        return prefix + " V2000"

    def fix_mol_block(self, mol_block: str) -> str:
        lines = mol_block.splitlines()
        if len(lines) < 4:
            return mol_block
        lines[3] = self.fix_mol_counts_line(lines[3])
        return "\n".join(lines)

    def load_xyz_file(self, file_path: str) -> Optional[Any]:
        """Load XYZ file and create RDKit Mol with charge prompt and bond determination."""
        if not self.host.state_manager.check_unsaved_changes():
            return None

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                raw_lines = f.readlines()

            lines = [ln.strip() for ln in raw_lines if not ln.strip().startswith("#")]
            while lines and not lines[0]:
                lines.pop(0)

            if len(lines) < 2:
                raise ValueError("XYZ file format error: too few lines")

            num_atoms = int(lines[0])
            if num_atoms == 0:
                raise ValueError("XYZ file has zero atoms")

            atoms_data = []
            for i, line in enumerate(lines[2 : 2 + num_atoms]):
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"Invalid atom data at line {i + 3}")
                symbol = parts[0].capitalize()
                try:
                    Chem.Atom(symbol)
                except (RuntimeError, ValueError):
                    settings = getattr(self.host.init_manager, "settings", {}) or {}
                    if settings.get("skip_chemistry_checks", False):
                        symbol = "C"
                    else:
                        raise ValueError(f"Unrecognized element symbol: {parts[0]}")
                atoms_data.append(
                    (symbol, float(parts[1]), float(parts[2]), float(parts[3]))
                )

            if not atoms_data:
                raise ValueError("No valid atoms found in XYZ file")

            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)

            settings = getattr(self.host.init_manager, "settings", {}) or {}
            skip_checks = bool(settings.get("skip_chemistry_checks", False))

            def _set_prop(m, key, val):
                try:
                    if isinstance(val, int):
                        m.SetIntProp(key, val)
                    elif isinstance(val, float):
                        m.SetDoubleProp(key, val)
                except (RuntimeError, TypeError, ValueError):
                    pass

            def _process(charge_val, use_rd_determine=True):
                if use_rd_determine:
                    try:
                        from rdkit.Chem import rdDetermineBonds

                        mol_copy = Chem.RWMol(mol)
                        rdDetermineBonds.DetermineBonds(mol_copy, charge=charge_val)
                        candidate = mol_copy.GetMol()
                        _set_prop(candidate, "_xyz_charge", charge_val)
                        return candidate
                    except (RuntimeError, ValueError, TypeError) as e:
                        raise e
                else:
                    self.estimate_bonds_from_distances(mol)
                    candidate = mol.GetMol()
                    _set_prop(candidate, "_xyz_charge", charge_val)
                    return candidate

            if skip_checks:
                final_mol = _process(0, use_rd_determine=False)
                _set_prop(final_mol, "_xyz_skip_checks", 1)
            else:
                final_mol = None
                settings = getattr(self.host.init_manager, "settings", {}) or {}
                # First try with charge 0 (per user's 'first try with 0 then ask' requirement)
                # but only if "Always ask" is not explicitly enabled in settings.
                if not settings.get("always_ask_charge", False):
                    try:
                        final_mol = _process(0, use_rd_determine=True)
                    except (RuntimeError, ValueError, TypeError):
                        final_mol = None

                # If still no final_mol (because always_ask is True, or charge 0 failed)
                if final_mol is None:
                    while True:
                        prompt_fn = getattr(self, "prompt_for_charge", None)
                        if callable(prompt_fn):
                            result = prompt_fn()
                            if isinstance(result, tuple) and len(result) == 3:
                                charge_val, ok, skip_flag = result
                            else:
                                charge_val, ok, skip_flag = 0, True, False
                        else:
                            charge_val, ok, skip_flag = 0, True, False

                        if not ok:
                            return None
                        if skip_flag:
                            final_mol = _process(0, use_rd_determine=False)
                            _set_prop(final_mol, "_xyz_skip_checks", 1)
                            break
                        try:
                            final_mol = _process(charge_val, use_rd_determine=True)
                            break
                        except (RuntimeError, ValueError, TypeError) as e:
                            if (
                                hasattr(self.host, "statusBar")
                                and self.host.statusBar()
                            ):
                                self.host.statusBar().showMessage(
                                    f"Chemistry failed for charge {charge_val}: {e}. Try a different charge or skip."
                                )
                            if not callable(prompt_fn):
                                raise e

                if final_mol:
                    final_mol._xyz_atom_data = atoms_data
            return final_mol

        except (RuntimeError, TypeError, ValueError, UnicodeDecodeError) as e:
            self.host.statusBar().showMessage(f"Error parsing XYZ file: {e}")
            return None

    def prompt_for_charge(self):
        """Show dialog to prompt user for molecular charge when loading XYZ files."""
        dialog = QDialog(self.host)
        dialog.setWindowTitle("Import XYZ Charge")
        layout = QVBoxLayout(dialog)
        line_edit = QLineEdit(dialog)
        line_edit.setText("0")
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        skip_btn = QPushButton("Skip chemistry", dialog)
        hl = QHBoxLayout()
        hl.addWidget(btn_box)
        hl.addWidget(skip_btn)
        layout.addWidget(QLabel("Enter total molecular charge:"))
        layout.addWidget(line_edit)
        layout.addLayout(hl)
        result = {"accepted": False, "skip": False}
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)
        skip_btn.clicked.connect(
            lambda: (result.update({"skip": True}), dialog.accept())
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None, False, False
        if result["skip"]:
            return 0, True, True
        try:
            val = int(float(line_edit.text().strip() or "0"))
            return val, True, False
        except ValueError:
            return 0, True, False

    def estimate_bonds_from_distances(self, mol) -> int:
        """Estimate bonds based on interatomic distances using covalent radii."""
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        bonds_added = []

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)
                distance = rdMolTransforms.GetBondLength(conf, i, j)
                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()
                radius_i = COVALENT_RADII.get(symbol_i, 1.0)
                radius_j = COVALENT_RADII.get(symbol_j, 1.0)
                expected = radius_i + radius_j
                tolerance = 1.2 if (symbol_i == "H" or symbol_j == "H") else 1.3
                if expected * 0.5 <= distance <= expected * tolerance:
                    if mol.GetBondBetweenAtoms(i, j) is None:
                        if (
                            symbol_i == "H" and mol.GetAtomWithIdx(i).GetDegree() >= 1
                        ) or (
                            symbol_j == "H" and mol.GetAtomWithIdx(j).GetDegree() >= 1
                        ):
                            continue
                        try:
                            mol.AddBond(i, j, Chem.BondType.SINGLE)
                            bonds_added.append((i, j, distance))
                        except (RuntimeError, ValueError, TypeError):
                            pass

        return len(bonds_added)

    def save_project(self) -> None:
        """Save (Ctrl+S) - Defaults to PMEPRJ format."""
        if (
            not self.host.state_manager.data.atoms
            and not self.host.view_3d_manager.current_mol
        ):
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return

        native_exts = [".pmeprj", ".pmeraw"]
        if self.host.init_manager.current_file_path and any(
            self.host.init_manager.current_file_path.lower().endswith(ext)
            for ext in native_exts
        ):
            try:
                if self.host.init_manager.current_file_path.lower().endswith(".pmeraw"):
                    save_data = self.host.state_manager.get_current_state()
                    with open(self.host.init_manager.current_file_path, "wb") as f:
                        pickle.dump(save_data, f)
                else:
                    json_data = self.host.state_manager.create_json_data()
                    with open(
                        self.host.init_manager.current_file_path, "w", encoding="utf-8"
                    ) as f:
                        json.dump(json_data, f, indent=2, ensure_ascii=False)

                self.host.state_manager.has_unsaved_changes = False
                self.host.state_manager.update_window_title()
                self.host.statusBar().showMessage(
                    f"Project saved to {self.host.init_manager.current_file_path}"
                )
            except (OSError, IOError) as e:
                self.host.statusBar().showMessage(f"File I/O error: {e}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving project: {e}")
        else:
            self.save_project_as()

    def save_project_as(self) -> None:
        """Save As (Ctrl+Shift+S) — always saves in PMEPRJ format."""
        import copy

        if (
            not self.host.state_manager.data.atoms
            and not self.host.view_3d_manager.current_mol
        ):
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return

        default_path = self._get_default_path()

        file_path, _ = QFileDialog.getSaveFileName(
            self.host,
            "Save Project As",
            default_path,
            "PME Project Files (*.pmeprj);;All Files (*)",
        )
        if not file_path:
            return

        if not file_path.lower().endswith(".pmeprj"):
            file_path += ".pmeprj"

        try:
            json_data = self.host.state_manager.create_json_data()
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)

            self.host.state_manager.has_unsaved_changes = False
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.update_window_title()
            try:
                self.host.state_manager._saved_state = copy.deepcopy(
                    self.host.state_manager.get_current_state()
                )
            except Exception:
                pass
            self.host.statusBar().showMessage(f"Project saved to {file_path}")
        except (OSError, IOError) as e:
            self.host.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Error saving project: {e}")

    def open_project(self) -> None:
        """Open an existing project file (legacy name — delegates to open_project_file)."""
        self.open_project_file()

    def open_project_file(self, file_path: Optional[str] = None) -> None:
        """Open project file (.pmeprj or .pmeraw)."""
        if not self.host.state_manager.check_unsaved_changes():
            return None

        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Open Project File",
                default_dir,
                "PME Project Files (*.pmeprj);;PME Raw Files (*.pmeraw);;All Files (*)",
            )
            if not file_path:
                return

        if file_path.lower().endswith(".pmeprj"):
            self.load_json_data(file_path)
        elif file_path.lower().endswith(".pmeraw"):
            self.load_raw_data(file_path)
        else:
            try:
                self.load_json_data(file_path)
            except Exception:
                try:
                    self.load_raw_data(file_path)
                except Exception:
                    self.host.statusBar().showMessage(
                        "Error: Unable to determine file format or file corrupted."
                    )

    def save_as_json(self) -> None:
        """Save as PME Project (JSON) format."""
        if (
            not self.host.state_manager.data.atoms
            and not self.host.view_3d_manager.current_mol
        ):
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return

        default_path = self._get_default_path()

        try:
            file_path, _ = QFileDialog.getSaveFileName(
                self.host,
                "Save as PME Project",
                default_path,
                "PME Project Files (*.pmeprj);;All Files (*)",
            )
            if not file_path:
                return

            if not file_path.lower().endswith(".pmeprj"):
                file_path += ".pmeprj"

            json_data = self.host.state_manager.create_json_data()
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(json_data, f, indent=2, ensure_ascii=False)

            self.host.state_manager.has_unsaved_changes = False
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.update_window_title()
            self.host.statusBar().showMessage(f"PME Project saved to {file_path}")
        except (OSError, IOError) as e:
            self.host.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Error saving PME Project file: {e}")

    def load_json_data(self, file_path: Optional[str] = None) -> None:
        """Load PME Project (.pmeprj) file."""
        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Open PME Project File",
                default_dir,
                "PME Project Files (*.pmeprj);;All Files (*)",
            )
            if not file_path:
                return

        if not self.host.edit_actions_manager.clear_all(skip_check=True):
            return

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                json_data = json.load(f)

            if json_data.get("format") != "PME Project":
                QMessageBox.warning(
                    self.host,
                    "Invalid Format",
                    "This file is not a valid PME Project format.",
                )
                return

            file_version = json_data.get("version", "1.0")
            if file_version != "1.0":
                QMessageBox.information(
                    self.host,
                    "Version Notice",
                    f"This file was created with PME Project version {file_version}.\n"
                    "Loading will be attempted but some features may not work correctly.",
                )

            self.host.ui_manager.restore_ui_for_editing()
            self.host.state_manager.load_from_json_data(json_data)
            self.host.state_manager.reset_undo_stack()
            self.host.state_manager.has_unsaved_changes = False
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.update_window_title()
            self.host.statusBar().showMessage(f"PME Project loaded from {file_path}")

            # Reset camera/zoom after drawing
            QTimer.singleShot(
                50, lambda: self.host.view_3d_manager.plotter.view_isometric()
            )
            QTimer.singleShot(100, lambda: self.host.view_3d_manager.plotter.render())
            QTimer.singleShot(100, self.host.view_3d_manager.fit_to_view)

        except FileNotFoundError:
            self.host.statusBar().showMessage(f"File not found: {file_path}")
        except json.JSONDecodeError as e:
            self.host.statusBar().showMessage(f"Invalid JSON format: {e}")
        except (OSError, IOError) as e:
            self.host.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e:
            self.host.statusBar().showMessage(
                f"Data corruption in PME Project file: {e}"
            )

    def save_raw_data(self) -> None:
        """Save as PME Raw (pickle) format."""
        import copy

        if (
            not self.host.state_manager.data.atoms
            and not self.host.view_3d_manager.current_mol
        ):
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return

        default_path = self._get_default_path()

        try:
            file_path, _ = QFileDialog.getSaveFileName(
                self.host,
                "Save Project File",
                default_path,
                "Project Files (*.pmeraw);;All Files (*)",
            )
            if not file_path:
                return

            if not file_path.lower().endswith(".pmeraw"):
                file_path += ".pmeraw"

            save_data = self.host.state_manager.get_current_state()
            with open(file_path, "wb") as f:
                pickle.dump(save_data, f)

            self.host.state_manager.has_unsaved_changes = False
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.update_window_title()
            try:
                self.host.state_manager._saved_state = copy.deepcopy(
                    self.host.state_manager.get_current_state()
                )
            except Exception:
                pass
            self.host.statusBar().showMessage(f"Project saved to {file_path}")
        except (OSError, IOError) as e:
            self.host.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Export error: {e}")

    def load_mol_file(self, file_path: Optional[str] = None) -> None:
        """Import a MOL/SDF file and add its contents to the 2D editor."""
        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Import MOL File",
                default_dir,
                "Chemical Files (*.mol *.sdf);;All Files (*)",
            )
            if not file_path:
                return

        if not os.path.exists(file_path):
            self.host.statusBar().showMessage(f"File not found: {file_path}")
            return

        try:
            if file_path.lower().endswith(".mol"):
                with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
            else:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)

            if mol is None:
                raise ValueError("Failed to read molecule from file.")

            Chem.Kekulize(mol)

            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            SCALE_FACTOR = 50.0
            existing_atoms = self.host.state_manager.data.atoms
            if existing_atoms:
                max_x = max(
                    v["pos"].x() if hasattr(v["pos"], "x") else v["pos"][0]
                    for v in existing_atoms.values()
                )
                avg_y = sum(
                    v["pos"].y() if hasattr(v["pos"], "y") else v["pos"][1]
                    for v in existing_atoms.values()
                ) / len(existing_atoms)
                place_center = QPointF(max_x + 80.0, avg_y)
            else:
                place_center = self.host.init_manager.view_2d.mapToScene(
                    self.host.init_manager.view_2d.viewport().rect().center()
                )

            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = (
                sum(p.x for p in positions) / len(positions) if positions else 0.0
            )
            mol_center_y = (
                sum(p.y for p in positions) / len(positions) if positions else 0.0
            )

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                scene_x = ((pos.x - mol_center_x) * SCALE_FACTOR) + place_center.x()
                scene_y = (-(pos.y - mol_center_y) * SCALE_FACTOR) + place_center.y()
                atom_id = self.host.init_manager.scene.create_atom(
                    atom.GetSymbol(),
                    QPointF(scene_x, scene_y),
                    charge=atom.GetFormalCharge(),
                )
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                stereo = 0
                b_dir = bond.GetBondDir()
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4
                self.host.init_manager.scene.create_bond(
                    self.host.state_manager.data.atoms[
                        rdkit_idx_to_my_id[bond.GetBeginAtomIdx()]
                    ]["item"],
                    self.host.state_manager.data.atoms[
                        rdkit_idx_to_my_id[bond.GetEndAtomIdx()]
                    ]["item"],
                    bond_order=int(bond.GetBondTypeAsDouble()),
                    bond_stereo=stereo,
                )

            self.host.statusBar().showMessage(f"Successfully imported {file_path}")
            self.host.init_manager.scene.update_all_items()
            self.host.edit_actions_manager.push_undo_state()
            QTimer.singleShot(100, self.host.view_3d_manager.fit_to_view)
        except Exception as e:
            self.host.statusBar().showMessage(f"Error loading file: {e}")

    def save_as_mol(self) -> None:
        """Save current 2D structure as MOL file."""
        try:
            mol_block = self.host.state_manager.data.to_mol_block()
            if not mol_block:
                self.host.statusBar().showMessage("Error: No 2D data to save.")
                return
            lines = mol_block.split("\n")
            if len(lines) > 1 and "RDKit" in lines[1]:
                lines[1] = f"  MoleditPy Ver. {VERSION}  2D"

            default_path = self._get_default_path(suffix="-2d")

            file_path, _ = QFileDialog.getSaveFileName(
                self.host,
                "Save 2D MOL File",
                default_path,
                "MOL Files (*.mol);;All Files (*)",
            )
            if file_path:
                if not file_path.lower().endswith(".mol"):
                    file_path += ".mol"
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
                self.host.statusBar().showMessage(f"2D data saved to {file_path}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Error saving MOL: {e}")

    def load_xyz_for_3d_viewing(self, file_path: Optional[str] = None) -> None:
        """Load XYZ file and display in 3D viewer."""
        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Load 3D XYZ (View Only)",
                default_dir,
                "XYZ Files (*.xyz);;All Files (*)",
            )
            if not file_path:
                return

        try:
            mol = self.load_xyz_file(file_path)
            if mol is None:
                raise ValueError("Failed to create molecule from XYZ file.")

            self.host.edit_actions_manager.clear_all(skip_check=True)
            self.host.view_3d_manager.current_mol = mol
            self.host.view_3d_manager.atom_id_to_rdkit_idx_map = {}

            # Determine is_xyz_derived: True only when bond estimation was skipped
            import contextlib

            skip_flag = False
            with contextlib.suppress(AttributeError, KeyError, RuntimeError, TypeError):
                if mol.HasProp("_xyz_skip_checks"):
                    skip_flag = bool(mol.GetIntProp("_xyz_skip_checks"))
            self.host.is_xyz_derived = skip_flag or (mol.GetNumBonds() == 0)

            if hasattr(self.host.view_3d_manager, "draw_molecule_3d"):
                self.host.view_3d_manager.draw_molecule_3d(mol)
            else:
                logging.error(
                    "REPORT ERROR: Missing attribute 'draw_molecule_3d' on object"
                )

            # Reset camera/zoom after drawing
            QTimer.singleShot(
                50, lambda: self.host.view_3d_manager.plotter.view_isometric()
            )
            QTimer.singleShot(100, lambda: self.host.view_3d_manager.plotter.render())

            if hasattr(self.host.ui_manager, "_enter_3d_viewer_ui_mode"):
                self.host.ui_manager._enter_3d_viewer_ui_mode()
            else:
                logging.error(
                    "REPORT ERROR: Missing attribute '_enter_3d_viewer_ui_mode' on object"
                )

            if hasattr(self.host.ui_manager, "_enable_3d_features"):
                self.host.ui_manager._enable_3d_features(True)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute '_enable_3d_features' on object"
                )

            if hasattr(self.host.view_3d_manager, "update_atom_id_menu_text"):
                self.host.view_3d_manager.update_atom_id_menu_text()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'update_atom_id_menu_text' on object"
                )
            if hasattr(self.host.view_3d_manager, "update_atom_id_menu_state"):
                self.host.view_3d_manager.update_atom_id_menu_state()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'update_atom_id_menu_state' on object"
                )

            if hasattr(self.host, "statusBar") and self.host.statusBar():
                self.host.statusBar().showMessage(
                    f"3D Viewer Mode: Loaded {os.path.basename(file_path)}"
                )
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.has_unsaved_changes = False
            self.host.state_manager.update_window_title()
        except Exception as e:
            if hasattr(self.host, "statusBar") and self.host.statusBar():
                self.host.statusBar().showMessage(f"XYZ Load failed: {e}")

    def load_mol_file_for_3d_viewing(self, file_path: Optional[str] = None) -> None:
        """Open MOL/SDF file in 3D viewer."""
        if not self.host.state_manager.check_unsaved_changes():
            return
        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Open MOL/SDF File",
                default_dir,
                "MOL/SDF Files (*.mol *.sdf);;All Files (*)",
            )
            if not file_path:
                return
        try:
            if file_path.lower().endswith(".sdf"):
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
            else:
                with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)

            if mol is None:
                raise ValueError("Failed to load molecule.")
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol)

            self.host.edit_actions_manager.clear_all(skip_check=True)
            self.host.view_3d_manager.current_mol = mol
            self.host.view_3d_manager.atom_id_to_rdkit_idx_map = {}
            self.host.is_xyz_derived = False
            self.host.view_3d_manager.draw_molecule_3d(mol)

            # Reset camera/zoom after drawing
            QTimer.singleShot(
                50, lambda: self.host.view_3d_manager.plotter.view_isometric()
            )
            QTimer.singleShot(100, lambda: self.host.view_3d_manager.plotter.render())

            self.host.ui_manager._enter_3d_viewer_ui_mode()
            if hasattr(self.host.ui_manager, "_enable_3d_features"):
                self.host.ui_manager._enable_3d_features(True)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute '_enable_3d_features' on object"
                )
            if hasattr(self.host.view_3d_manager, "update_atom_id_menu_text"):
                self.host.view_3d_manager.update_atom_id_menu_text()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'update_atom_id_menu_text' on object"
                )
            if hasattr(self.host.view_3d_manager, "update_atom_id_menu_state"):
                self.host.view_3d_manager.update_atom_id_menu_state()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'update_atom_id_menu_state' on object"
                )
            self.host.statusBar().showMessage(f"Loaded {file_path} in 3D viewer")
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.has_unsaved_changes = False
            self.host.state_manager.update_window_title()
        except Exception as e:
            self.host.statusBar().showMessage(f"3D MOL Load failed: {e}")

    def save_3d_as_mol(self) -> None:
        """Save current 3D structure as MOL file."""
        if not self.host.view_3d_manager.current_mol:
            self.host.statusBar().showMessage("Error: No 3D structure to save.")
            return
        default_path = self._get_default_path()
        file_path, _ = QFileDialog.getSaveFileName(
            self.host,
            "Save 3D MOL File",
            default_path,
            "MOL Files (*.mol);;All Files (*)",
        )
        if file_path:
            try:
                mol_block = Chem.MolToMolBlock(
                    self.host.view_3d_manager.current_mol, includeStereo=True
                )
                lines = mol_block.split("\n")
                if len(lines) > 1 and "RDKit" in lines[1]:
                    lines[1] = f"  MoleditPy Ver. {VERSION}  3D"
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
                self.host.statusBar().showMessage(f"3D data saved to {file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving 3D MOL: {e}")

    def save_as_xyz(self) -> None:
        """Save current 3D structure as XYZ file."""
        if not self.host.view_3d_manager.current_mol:
            self.host.statusBar().showMessage(
                "Error: Please generate a 3D structure first."
            )
            return

        default_path = self._get_default_path()
        file_path, _ = QFileDialog.getSaveFileName(
            self.host,
            "Save 3D XYZ File",
            default_path,
            "XYZ Files (*.xyz);;All Files (*)",
        )
        if file_path:
            if not file_path.lower().endswith(".xyz"):
                file_path += ".xyz"
            try:
                conf = self.host.view_3d_manager.current_mol.GetConformer()
                num_atoms = self.host.view_3d_manager.current_mol.GetNumAtoms()
                xyz_lines = [str(num_atoms)]

                charge = 0
                if self.host.view_3d_manager.current_mol.HasProp("_xyz_charge"):
                    charge = self.host.view_3d_manager.current_mol.GetIntProp(
                        "_xyz_charge"
                    )
                else:
                    try:
                        charge = Chem.GetFormalCharge(
                            self.host.view_3d_manager.current_mol
                        )
                    except Exception as e:
                        logging.warning("Could not compute formal charge: %s", e)

                multiplicity = 1
                try:
                    multiplicity = (
                        Descriptors.NumRadicalElectrons(
                            self.host.view_3d_manager.current_mol
                        )
                        + 1
                    )
                except Exception as e:
                    logging.warning("Could not compute multiplicity: %s", e)

                xyz_lines.append(
                    f"chrg = {charge}  mult = {multiplicity} | Generated by MoleditPy Ver. {VERSION}"
                )
                for i in range(num_atoms):
                    pos = conf.GetAtomPosition(i)
                    symbol = self.host.view_3d_manager.current_mol.GetAtomWithIdx(
                        i
                    ).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(xyz_lines) + "\n")
                self.host.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving XYZ: {e}")

    def load_raw_data(self, file_path: Optional[str] = None) -> None:
        """Open a .pmeraw pickle project file."""
        if not file_path:
            default_dir = (
                os.path.dirname(self.host.init_manager.current_file_path)
                if self.host.init_manager.current_file_path
                else ""
            )
            file_path, _ = QFileDialog.getOpenFileName(
                self.host,
                "Open Project File",
                default_dir,
                "Project Files (*.pmeraw);;All Files (*)",
            )
            if not file_path:
                return

        if not self.host.edit_actions_manager.clear_all(skip_check=True):
            return

        try:
            with open(file_path, "rb") as f:
                loaded_data = pickle.load(f)
            self.host.ui_manager.restore_ui_for_editing()
            self.host.state_manager.set_state_from_data(loaded_data)
            self.host.state_manager.reset_undo_stack()
            self.host.state_manager.has_unsaved_changes = False
            self.host.init_manager.current_file_path = file_path
            self.host.state_manager.update_window_title()
            self.host.statusBar().showMessage(f"Project loaded from {file_path}")

            # Reset camera/zoom after drawing
            QTimer.singleShot(
                50, lambda: self.host.view_3d_manager.plotter.view_isometric()
            )
            QTimer.singleShot(100, lambda: self.host.view_3d_manager.plotter.render())

        except FileNotFoundError:
            self.host.statusBar().showMessage(f"File not found: {file_path}")
        except (OSError, IOError) as e:
            self.host.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.UnpicklingError as e:
            self.host.statusBar().showMessage(f"Invalid project file format: {e}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Error loading project file: {e}")

    def _set_mol_prop(self, mol: Any, prop_name: str, value: Any) -> None:
        """Set an RDKit molecule property safely."""
        try:
            if isinstance(value, int):
                mol.SetIntProp(prop_name, value)
            elif isinstance(value, float):
                mol.SetDoubleProp(prop_name, value)
            else:
                mol.SetProp(prop_name, str(value))
        except Exception:
            pass

    def _get_mol_prop(self, mol: Any, prop_name: str, default: Any = None) -> Any:
        """Get an RDKit molecule property safely, trying int/float/str in order."""
        try:
            if not mol.HasProp(prop_name):
                return default
            for getter in [mol.GetIntProp, mol.GetDoubleProp, mol.GetProp]:
                try:
                    return getter(prop_name)
                except Exception:
                    continue
            return default
        except Exception:
            return default


def _set_mol_prop_safe(mol: Any, key: str, val: Any) -> None:
    """Module-level helper: set an int or float property on an RDKit mol silently."""
    import contextlib

    with contextlib.suppress(RuntimeError, TypeError, ValueError, AttributeError):
        if isinstance(val, int):
            mol.SetIntProp(key, val)
        elif isinstance(val, float):
            mol.SetDoubleProp(key, val)


# Backward-compat aliases
IOManager = IOManager
