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
Module separated from MainWindow (main_window.py)
Functional class: MainWindowViewLoaders
"""

import contextlib
import logging
import os

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem

# PyQt6 Modules
from PyQt6.QtWidgets import QFileDialog

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except (ImportError, AttributeError, TypeError):
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import VERSION
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.utils.constants import VERSION


# --- Class Definition ---
class MainWindowViewLoaders:
    """Mixin class separated from main_window.py"""

    def load_xyz_for_3d_viewing(self, file_path=None):
        """Load XYZ file and display in 3D viewer"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Load 3D XYZ (View Only)", "", "XYZ Files (*.xyz);;All Files (*)"
            )
            if not file_path:
                return

        try:
            mol = self.load_xyz_file(file_path)
            if mol is None:
                raise ValueError("Failed to create molecule from XYZ file.")
            if mol.GetNumConformers() == 0:
                raise ValueError("XYZ file has no 3D coordinates.")

            # Clear 2D editor
            if hasattr(self, "clear_2d_editor"):
                self.clear_2d_editor(push_to_undo=False)

            self.current_mol = mol
            self.atom_id_to_rdkit_idx_map = {}

            # Check skip flag
            skip_flag = False
            # Suppress potential errors if the coordinate skip flag is missing or RDKit property system fails
            with contextlib.suppress(AttributeError, KeyError, RuntimeError, TypeError):
                if mol.HasProp("_xyz_skip_checks"):
                    skip_flag = bool(mol.GetIntProp("_xyz_skip_checks"))
                else:
                    skip_flag = bool(getattr(mol, "_xyz_skip_checks", False))

            # Update optimization state
            self.is_xyz_derived = skip_flag or (mol.GetNumBonds() == 0)
            opt_btn = getattr(self, "optimize_3d_button", None)
            if opt_btn:
                opt_btn.setEnabled(not self.is_xyz_derived)

            # Update visualization and UI mode
            if hasattr(self, "draw_molecule_3d"):
                self.draw_molecule_3d(self.current_mol)

            if hasattr(self, "plotter"):
                self.plotter.reset_camera()

            if hasattr(self, "_enter_3d_viewer_ui_mode"):
                self._enter_3d_viewer_ui_mode()

            if hasattr(self, "_enable_3d_features"):
                self._enable_3d_features(True)

            if hasattr(self, "update_atom_id_menu_text"):
                self.update_atom_id_menu_text()
            if hasattr(self, "update_atom_id_menu_state"):
                self.update_atom_id_menu_state()

            self.statusBar().showMessage(
                f"3D Viewer Mode: Loaded {os.path.basename(file_path)}"
            )
            if hasattr(self, "reset_undo_stack"):
                self.reset_undo_stack()

            self.current_file_path = file_path
            self.has_unsaved_changes = False
            if hasattr(self, "update_window_title"):
                self.update_window_title()

        except (
            FileNotFoundError,
            ValueError,
            RuntimeError,
            TypeError,
            AttributeError,
        ) as e:
            # File loading or coordinate validation error reported to user via status bar.
            # We catch AttributeError/TypeError here to handle malformed RDKit molecule objects gracefully.
            self.statusBar().showMessage(f"Load failed: {e}")
            if hasattr(self, "restore_ui_for_editing"):
                self.restore_ui_for_editing()

    def save_3d_as_mol(self):
        """Save the current 3D structure as a MOL file."""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        # Determine default save path
        curr_path = getattr(self, "current_file_path", None)
        default_dir = os.path.dirname(curr_path) if curr_path else ""
        default_name = (
            os.path.splitext(os.path.basename(curr_path))[0]
            if curr_path
            else "untitled"
        )
        default_save_path = os.path.join(default_dir, default_name)

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save 3D MOL File",
            default_save_path,
            "MOL Files (*.mol);;All Files (*)",
        )
        if not file_path:
            return

        if not file_path.lower().endswith(".mol"):
            file_path += ".mol"

        try:
            mol_to_save = Chem.Mol(self.current_mol)
            if mol_to_save.HasProp("_2D"):
                mol_to_save.ClearProp("_2D")

            mol_block = Chem.MolToMolBlock(mol_to_save, includeStereo=True)
            lines = mol_block.split("\n")
            if len(lines) > 1 and "RDKit" in lines[1]:
                lines[1] = f"  MoleditPy Ver. {VERSION}  3D"

            with open(file_path, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
            self.statusBar().showMessage(f"3D data saved to {file_path}")

        except (IOError, OSError, RuntimeError, ValueError, TypeError) as e:
            # Catch I/O errors and RDKit serialization failures (RuntimeError/ValueError).
            # We also catch TypeError in case of unexpected nil objects during teardown or property access.
            msg = f"Error saving 3D MOL: {str(e)}"
            self.statusBar().showMessage(msg)

    def load_mol_file_for_3d_viewing(self, file_path=None):
        """Open MOL/SDF file in 3D viewer"""
        if not self.check_unsaved_changes():
            return  # Do nothing if user cancels
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Open MOL/SDF File",
                "",
                "MOL/SDF Files (*.mol *.sdf);;All Files (*)",
            )
            if not file_path:
                return
        try:
            # Determine extension
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ""

            # Load molecule
            mol = None
            if ext == ".sdf":
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
            elif ext == ".mol":
                with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
                if mol is None:
                    mol = Chem.MolFromMolFile(file_path, removeHs=False)
            else:
                mol = Chem.MolFromMolFile(file_path, removeHs=False)

            if mol is None:
                self.statusBar().showMessage(
                    f"Failed to load molecule from {file_path}"
                )
                return

            # Convert 2D to 3D if needed
            if mol.GetNumConformers() == 0:
                self.statusBar().showMessage(
                    "No 3D coordinates found. Converting to 3D..."
                )
                try:
                    AllChem.EmbedMolecule(mol)
                    self.current_mol = mol
                    if hasattr(self, "push_undo_state"):
                        self.push_undo_state()
                except (RuntimeError, ValueError, TypeError):
                    if not getattr(self, "settings", {}).get(
                        "skip_chemistry_checks", False
                    ):
                        self.statusBar().showMessage(
                            "Failed to generate 3D coordinates"
                        )
                        return
                    self.statusBar().showMessage(
                        "Warning: failed to generate 3D coordinates (continuing due to settings)"
                    )

            # Finalize load
            try:
                if hasattr(self, "_clear_xyz_flags"):
                    self._clear_xyz_flags(mol)
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                # Suppress non-critical errors if the molecule lacks XYZ source flags to clear.
                logging.debug(
                    f"Failed to clear XYZ flags (likely not an XYZ source): {e}"
                )
            if hasattr(self, "_apply_chem_check_and_set_flags"):
                self._apply_chem_check_and_set_flags(mol, source_desc="MOL/SDF")

            self.current_mol = mol
            if hasattr(self, "draw_molecule_3d"):
                self.draw_molecule_3d(mol)
            if hasattr(self, "plotter"):
                self.plotter.reset_camera()
            if hasattr(self, "_enter_3d_viewer_ui_mode"):
                self._enter_3d_viewer_ui_mode()

            if hasattr(self, "update_atom_id_menu_text"):
                self.update_atom_id_menu_text()
            if hasattr(self, "update_atom_id_menu_state"):
                self.update_atom_id_menu_state()

            self.statusBar().showMessage(f"Loaded {file_path} in 3D viewer")
            if hasattr(self, "reset_undo_stack"):
                self.reset_undo_stack()

            self.has_unsaved_changes = False
            self.current_file_path = file_path
            if hasattr(self, "update_window_title"):
                self.update_window_title()

        except (
            RuntimeError,
            ValueError,
            TypeError,
            ImportError,
            UnicodeDecodeError,
        ) as e:
            # Catch RDKit parsing errors (RuntimeError), UI data mismatches (ValueError/TypeError),
            # and file encoding issues during manual fix_mol_block processing.
            self.statusBar().showMessage(f"Error processing MOL/SDF file: {e}")
