#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import contextlib
import os

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdGeometry, rdMolTransforms
from ..utils.constants import VERSION

# PyQt6 Modules
from PyQt6.QtCore import QPointF, QTimer
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from ..utils.constants import (
        VERSION,
        COVALENT_RADII,
        BOND_ESTIMATION_MAX_DIST,
        BOND_ESTIMATION_MIN_DIST,
        BOND_ESTIMATION_TOLERANCE,
    )
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.utils.constants import (
        VERSION,
        COVALENT_RADII,
    )


def _set_mol_prop_safe(mol, key, val):
    """Set an integer property on an RDKit mol, silently ignoring failures."""
    # Suppress potential RDKit-specific crashes when setting properties on invalid/malformed atoms
    with contextlib.suppress(RuntimeError, TypeError, ValueError):
        mol.SetIntProp(key, int(val))


# --- Class Definition ---
class MainWindowMolecularParsers:
    def load_mol_file(self, file_path=None):
        if not self.check_unsaved_changes():
            return

        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Import MOL File",
                "",
                "Chemical Files (*.mol *.sdf);;All Files (*)",
            )
            if not file_path:
                return

        if not os.path.exists(file_path):
            self.statusBar().showMessage(f"File not found: {file_path}")
            return

        try:
            self.dragged_atom_info = None
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ""

            if ext == ".mol":
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
            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)

            # Generate 2D coords if missing
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            SCALE_FACTOR = 50.0
            view_center = self.view_2d.mapToScene(
                self.view_2d.viewport().rect().center()
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
                charge = atom.GetFormalCharge()

                scene_x = ((pos.x - mol_center_x) * SCALE_FACTOR) + view_center.x()
                scene_y = (-(pos.y - mol_center_y) * SCALE_FACTOR) + view_center.y()

                atom_id = self.scene.create_atom(
                    atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge
                )
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2

                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4

                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item, a2_item = (
                    self.data.atoms[a1_id]["item"],
                    self.data.atoms[a2_id]["item"],
                )
                self.scene.create_bond(
                    a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo
                )

            self.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.reset_undo_stack()
            self.current_file_path = file_path
            self.has_unsaved_changes = False
            self.update_window_title()
            # Request scene update and ring re-analysis
            self.scene.update_all_items()
            QTimer.singleShot(0, self.fit_to_view)

        except (RuntimeError, ValueError, TypeError, UnicodeDecodeError) as e:
            # File loading or coordinate validation error reported to user via status bar.
            # We catch RuntimeError/ValueError for RDKit parsing and UnicodeDecodeError for encoding issues.
            import logging

            logging.error(f"Error loading MOL/SDF file: {e}", exc_info=True)
            self.statusBar().showMessage(f"Error loading file: {e}")

    def _set_mol_prop(self, mol, prop_name, value):
        """Internal helper to set molecule properties safely."""
        try:
            if isinstance(value, int):
                mol.SetIntProp(prop_name, value)
            elif isinstance(value, float):
                mol.SetDoubleProp(prop_name, value)
            else:
                mol.SetProp(prop_name, str(value))
        except (
            AttributeError,
            RuntimeError,
            TypeError,
            ValueError,
            Chem.rdchem.MolPropException,
        ):
            # Suppress errors if RDKit property setter fails (e.g., invalid type or concurrent mol access).
            # These properties are metadata and should not block core molecule loading logic.
            pass  # Metadata failure is non-critical

    def _get_mol_prop(self, mol, prop_name, default=None):
        """Internal helper to get molecule properties safely."""
        try:
            if not mol.HasProp(prop_name):
                return default
            # Try specific types first, then string fallback
            for getter in [mol.GetIntProp, mol.GetDoubleProp, mol.GetProp]:
                try:
                    return getter(prop_name)
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    # If one getter fails, try the next one in the sequence.
                    continue
            return default
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Property retrieval failure returns the default value.
            return default

    def prompt_for_charge(self):
        """Helper dialog to prompt for charge or skip chemistry."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Import XYZ Charge")
        layout = QVBoxLayout(dialog)
        line_edit = QLineEdit(dialog)
        line_edit.setText("0")
        btn_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog
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

        if dialog.exec() != QDialog.Accepted:
            return None, False, False

        if result["skip"]:
            return 0, True, True

        try:
            val = int(float(line_edit.text().strip() or "0"))
            return val, True, False
        except ValueError:
            return 0, True, False

    def load_xyz_file(self, file_path):
        """Load XYZ file and create RDKit Mol."""
        if not self.check_unsaved_changes():
            return

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                raw_lines = f.readlines()

            # Filter out comment lines ONLY if they start with #
            # (Standard XYZ doesn't use # but some extended versions do)
            lines = [ln.strip() for ln in raw_lines if not ln.strip().startswith("#")]
            # Remove any leading empty lines to find the atom count
            while lines and not lines[0]:
                lines.pop(0)

            if len(lines) < 2:
                raise ValueError("XYZ file format error: too few lines")

            num_atoms = int(lines[0])
            atoms_data = []
            # In XYZ, lines[0] is count, lines[1] is title, lines[2...] are atoms
            for i, line in enumerate(lines[2 : 2 + num_atoms]):
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"Invalid atom data at line {i + 3}")
                symbol = parts[0].capitalize()
                try:
                    Chem.Atom(symbol)
                except (RuntimeError, ValueError):
                    # Catch RDKit atom creation failures due to unrecognized symbols.
                    if self.settings.get("skip_chemistry_checks", False):
                        symbol = "C"
                    else:
                        raise ValueError(f"Unrecognized element symbol: {parts[0]}")
                atoms_data.append(
                    (symbol, float(parts[1]), float(parts[2]), float(parts[3]))
                )

            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)

            skip_checks = bool(self.settings.get("skip_chemistry_checks", False))

            final_mol = None
            used_rd_determine = False

            def _process(charge_val, use_rd_determine=True):
                nonlocal used_rd_determine
                if use_rd_determine:
                    try:
                        from rdkit.Chem import rdDetermineBonds

                        mol_copy = Chem.RWMol(mol)
                        rdDetermineBonds.DetermineBonds(mol_copy, charge=charge_val)
                        candidate = mol_copy.GetMol()
                        used_rd_determine = True
                        self._apply_chem_check_and_set_flags(
                            candidate, source_desc="XYZ"
                        )
                        if getattr(self, "chem_check_failed", False):
                            raise ValueError("Sanitization failed")
                        _set_mol_prop_safe(candidate, "_xyz_charge", charge_val)
                        return candidate
                    except (RuntimeError, ValueError, TypeError) as e:
                        # Catch RDKit bond determination failures (e.g., valence overflow).
                        used_rd_determine = False
                        raise e
                else:
                    self.estimate_bonds_from_distances(mol)
                    candidate = mol.GetMol()
                    used_rd_determine = False
                    self._apply_chem_check_and_set_flags(
                        candidate, source_desc="XYZ (Skipped)", force_skip=True
                    )
                    _set_mol_prop_safe(candidate, "_xyz_charge", charge_val)
                    return candidate

            if skip_checks:
                # Skip chemistry: just estimate bonds, no dialog
                final_mol = _process(0, use_rd_determine=False)
                _set_mol_prop_safe(final_mol, "_xyz_skip_checks", 1)
            else:
                # Always ask user for charge (dialog defaults to 0)
                while True:
                    charge_val, ok, skip_flag = self.prompt_for_charge()
                    if not ok:
                        return None
                    if skip_flag:
                        final_mol = _process(0, use_rd_determine=False)
                        _set_mol_prop_safe(final_mol, "_xyz_skip_checks", 1)
                        break
                    try:
                        final_mol = _process(charge_val, use_rd_determine=True)
                        break
                    except (RuntimeError, ValueError, TypeError) as e:
                        # Catch chemistry failures and prompt the user to adjust charge or skip.
                        self.statusBar().showMessage(
                            f"Chemistry failed for charge {charge_val}: {e}. Try a different charge or skip."
                        )
            return final_mol

        except (RuntimeError, TypeError, ValueError, UnicodeDecodeError) as e:
            # XYZ parsing error reported to user via status bar.
            # We catch RuntimeError/ValueError for parsing and UnicodeDecodeError for file encoding.
            import logging

            logging.error(f"Error parsing XYZ file: {e}", exc_info=True)
            self.statusBar().showMessage(f"Error parsing XYZ file: {e}")
            return None

    def estimate_bonds_from_distances(self, mol):
        """Estimate bonds based on interatomic distances."""

        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()

        # Track added bonds
        bonds_added = []

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)

                distance = rdMolTransforms.GetBondLength(conf, i, j)

                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()

                radius_i = COVALENT_RADII.get(symbol_i, 1.0)  # Default
                radius_j = COVALENT_RADII.get(symbol_j, 1.0)

                expected_bond_length = radius_i + radius_j

                if symbol_i == "H" or symbol_j == "H":
                    tolerance_factor = 1.2  # H bonds often shorter
                else:
                    tolerance_factor = 1.3  # Margin for others

                max_bond_length = expected_bond_length * tolerance_factor
                min_bond_length = expected_bond_length * 0.5

                if min_bond_length <= distance <= max_bond_length:
                    if mol.GetBondBetweenAtoms(i, j) is None:
                        # Limit Hydrogen to only one bond
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
                            # Best-effort: ignore bond adding errors for invalid distances
                            # or if the molecule state becomes inconsistent during loop.
                            pass

        # Debug information (optional)
        # Added bonds based on distance analysis

        return len(bonds_added)

    def save_as_mol(self):
        try:
            mol_block = self.data.to_mol_block()
            if not mol_block:
                self.statusBar().showMessage("Error: No 2D data to save.")
                return

            lines = mol_block.split("\n")
            if len(lines) > 1 and "RDKit" in lines[1]:
                lines[1] = "  MoleditPy Ver. " + VERSION + "  2D"
            modified_mol_block = "\n".join(lines)

            # default filename: based on current_file_path, append -2d for 2D mol
            default_name = "untitled-2d"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    name = os.path.splitext(base)[0]
                    default_name = f"{name}-2d"
            except (AttributeError, RuntimeError, ValueError, TypeError):
                default_name = "untitled-2d"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(
                        os.path.dirname(self.current_file_path), default_name
                    )
            except (AttributeError, RuntimeError, ValueError, TypeError):
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self,
                "Save 2D MOL File",
                default_path,
                "MOL Files (*.mol);;All Files (*)",
            )
            if not file_path:
                return

            if not file_path.lower().endswith(".mol"):
                file_path += ".mol"

            with open(file_path, "w", encoding="utf-8") as f:
                f.write(modified_mol_block)
            self.statusBar().showMessage(f"2D data saved to {file_path}")

        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except UnicodeEncodeError as e:
            self.statusBar().showMessage(f"Text encoding error: {e}")
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            import logging

            logging.error(f"Error saving MOL file: {e}", exc_info=True)
            self.statusBar().showMessage(f"Error saving file: {e}")

    def save_as_xyz(self):
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        default_name = "untitled"
        if self.current_file_path:
            with contextlib.suppress(
                AttributeError, RuntimeError, ValueError, TypeError
            ):
                default_name = os.path.splitext(
                    os.path.basename(self.current_file_path)
                )[0]

        default_path = default_name
        if self.current_file_path:
            with contextlib.suppress(
                AttributeError, RuntimeError, ValueError, TypeError
            ):
                default_path = os.path.join(
                    os.path.dirname(self.current_file_path), default_name
                )

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save 3D XYZ File", default_path, "XYZ Files (*.xyz);;All Files (*)"
        )
        if not file_path:
            return
        if not file_path.lower().endswith(".xyz"):
            file_path += ".xyz"

        try:
            conf = self.current_mol.GetConformer()
            num_atoms = self.current_mol.GetNumAtoms()
            xyz_lines = [str(num_atoms)]

            charge = 0
            if hasattr(self.current_mol, "HasProp") and self.current_mol.HasProp(
                "_xyz_charge"
            ):
                charge = self.current_mol.GetIntProp("_xyz_charge")
            else:
                # Suppress potential errors during formal charge calculation for unvalidated RDKit mools
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    charge = Chem.GetFormalCharge(self.current_mol)

            multiplicity = 1
            with contextlib.suppress(
                AttributeError, RuntimeError, TypeError, ValueError
            ):
                multiplicity = Descriptors.NumRadicalElectrons(self.current_mol) + 1

            xyz_lines.append(
                f"chrg = {charge}  mult = {multiplicity} | Generated by MoleditPy Ver. {VERSION}"
            )
            for i in range(num_atoms):
                pos = conf.GetAtomPosition(i)
                symbol = self.current_mol.GetAtomWithIdx(i).GetSymbol()
                xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

            with open(file_path, "w", encoding="utf-8") as f:
                f.write("\n".join(xyz_lines) + "\n")
            self.statusBar().showMessage(f"Successfully saved to {file_path}")
        except (IOError, OSError, RuntimeError, ValueError, TypeError) as e:
            # Catch I/O errors and RDKit serialization failures (RuntimeError/ValueError).
            import logging

            logging.error(f"Error saving XYZ file: {e}", exc_info=True)
            self.statusBar().showMessage(f"Error saving file: {e}")

    def fix_mol_counts_line(self, line: str) -> str:
        """
        Check and fix the CTAB counts line in a MOL file.
        If the line already contains 'V3000' or 'V2000' it is left unchanged.
        Otherwise the line is treated as V2000 and the proper 39-character
        format (33 chars of counts + ' V2000') is returned.
        """
        # If already V3000 or V2000, leave as-is
        if "V3000" in line or "V2000" in line:
            return line

        # Prepare prefix (first 33 characters for the 11 * I3 fields)
        prefix = line.rstrip().ljust(33)[0:33]
        version_str = " V2000"
        return prefix + version_str

    def fix_mol_block(self, mol_block: str) -> str:
        """Ensure CTAB counts line (4th) is valid."""
        lines = mol_block.splitlines()
        if len(lines) < 4:
            # Not a valid MOL block — return unchanged
            return mol_block

        counts_line = lines[3]
        fixed_counts_line = self.fix_mol_counts_line(counts_line)
        lines[3] = fixed_counts_line
        return "\n".join(lines)
