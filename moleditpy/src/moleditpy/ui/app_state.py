#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

import base64
import re
import binascii
import copy
import logging
import os
import contextlib
from typing import Any, Dict, Optional, Tuple


import numpy as np

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# PyQt6 Modules
from PyQt6.QtCore import QDateTime, Qt
from PyQt6.QtWidgets import QMessageBox

from ..utils.constants import VERSION
from ..core.molecular_data import MolecularData


from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .main_window import MainWindow


# --- Class Definition ---
def _serialize_constraints(constraints: list) -> list:
    """Convert internal constraint tuples to JSON-serializable lists."""
    result = []
    for const in constraints:
        if len(const) == 4:
            result.append([const[0], list(const[1]), const[2], const[3]])
        else:
            result.append([const[0], list(const[1]), const[2], 1.0e5])
    return result


def _deserialize_constraints(raw: list) -> list:
    """Convert JSON-loaded constraint lists back to internal tuples."""
    result = []
    for const in raw:
        if not isinstance(const, (list, tuple)):
            continue
        try:
            if len(const) == 4:
                result.append((const[0], tuple(const[1]), const[2], const[3]))
            elif len(const) == 3:
                result.append((const[0], tuple(const[1]), const[2], 1.0e5))
        except (TypeError, ValueError, IndexError) as e:
            logging.debug(f"Failed to parse constraint {const}: {e}")
    return result


class StateManager:
    def __init__(self, host: MainWindow) -> None:
        self.host = host
        self.data: MolecularData  # Dynamically assigned in main_window_init.py
        self.has_unsaved_changes = False
        self._preserved_plugin_data: Dict[str, Any] = {}
        self.dragged_atom_info: Optional[Dict[str, Any]] = None
        self.saved_state: Optional[Dict[str, Any]] = None

    def get_current_state(self) -> Dict[str, Any]:
        atoms = {
            atom_id: {
                "symbol": data["symbol"],
                "pos": data["pos"],  # Already a tuple (x, y)
                "charge": data.get("charge", 0),
                "radical": data.get("radical", 0),
            }
            for atom_id, data in self.data.atoms.items()
        }
        bonds = {
            key: {"order": data["order"], "stereo": data.get("stereo", 0)}
            for key, data in self.data.bonds.items()
        }
        state = {
            "atoms": atoms,
            "bonds": bonds,
            "_next_atom_id": self.data.next_atom_id,
        }

        state["version"] = VERSION

        if self.host.current_mol:
            state["mol_3d"] = self.host.current_mol.ToBinary()
            # RDKit binary serialization does not preserve custom properties like _original_atom_id.
            # We store them separately to ensure we can restore the 2D-3D link after undo/redo.
            mol_3d_atom_ids = []
            for i in range(self.host.current_mol.GetNumAtoms()):
                atom = self.host.current_mol.GetAtomWithIdx(i)
                if atom and atom.HasProp("_original_atom_id"):
                    try:
                        mol_3d_atom_ids.append(atom.GetIntProp("_original_atom_id"))
                    except (RuntimeError, ValueError, TypeError):
                        mol_3d_atom_ids.append(None)
                else:
                    mol_3d_atom_ids.append(None)
            state["mol_3d_atom_ids"] = mol_3d_atom_ids

        state["is_3d_viewer_mode"] = not self.host.ui_manager.is_2d_editable

        state["constraints_3d"] = _serialize_constraints(
            self.host.edit_3d_manager.constraints_3d
        )

        return state

    def set_state_from_data(self, state_data: Dict[str, Any]) -> None:
        self.dragged_atom_info = None
        self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)

        loaded_data = copy.deepcopy(state_data)

        # Get file version (default '0.0.0')
        file_version_str = loaded_data.get("version", "0.0.0")

        def parse_v(v_str: str) -> Tuple[int, ...]:
            if not v_str or not isinstance(v_str, str):
                return (0, 0, 0)
            try:
                parts = []
                for p in v_str.split("."):
                    m = re.match(r"(\d+)", p)
                    parts.append(int(m.group(1)) if m else 0)
                return tuple(parts)
            except (ValueError, AttributeError, IndexError):
                return (0, 0, 0)

        app_version_parts = parse_v(VERSION)
        file_version_parts = parse_v(file_version_str)

        # Warn if file version is newer than app version
        if file_version_parts > app_version_parts:
            self.host.warning_message_box(
                "Version Mismatch",
                f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).",
            )

        raw_atoms = loaded_data.get("atoms", {})
        raw_bonds = loaded_data.get("bonds", {})

        loaded_constraints = loaded_data.get("constraints_3d", [])
        self.host.set_constraints_3d(
            _deserialize_constraints(loaded_constraints)
            if isinstance(loaded_constraints, list)
            else []
        )

        self.host.scene.restore_atoms_and_bonds(raw_atoms, raw_bonds)
        self.data.next_atom_id = loaded_data.get(
            "_next_atom_id",
            max(self.data.atoms.keys()) + 1 if self.data.atoms else 0,
        )
        mol_3d_data = loaded_data.get("mol_3d")
        if mol_3d_data is not None:
            try:
                self.host.set_current_molecule(Chem.Mol(mol_3d_data))
                # Debug: check if 3D structure is valid
                if (
                    self.host.current_mol
                    and self.host.current_mol.GetNumAtoms() > 0
                ):
                    # Restore _original_atom_id if present in saved state
                    mol_3d_atom_ids = loaded_data.get("mol_3d_atom_ids")
                    if (
                        mol_3d_atom_ids
                        and len(mol_3d_atom_ids)
                        == self.host.current_mol.GetNumAtoms()
                    ):
                        for i, aid in enumerate(mol_3d_atom_ids):
                            if aid is not None:
                                rd_atom = self.host.current_mol.GetAtomWithIdx(
                                    i
                                )
                                if rd_atom:
                                    try:
                                        rd_atom.SetIntProp(
                                            "_original_atom_id", int(aid)
                                        )
                                    except (RuntimeError, ValueError, TypeError):
                                        # Safe defensive fallback catching RuntimeError, ValueError, TypeError
                                        pass

                    # Sync 2D atoms with 3D actors
                    if self.host.compute_manager:
                        try:
                            self.host.compute_manager.create_atom_id_mapping()
                            self.host.view_3d_manager.update_atom_id_menu_text()
                            self.host.view_3d_manager.update_atom_id_menu_state()
                        except (RuntimeError, AttributeError) as e:
                            logging.debug(
                                f"Partial failure during ID mapping restoration: {e}"
                            )

                    # draw_molecule_3d will use restored IDs
                    self.host.view_3d_manager.draw_molecule_3d(
                        self.host.current_mol
                    )
                    if self.host.plotter:
                        self.host.plotter.reset_camera()

                    self.host.ui_manager.enable_3d_features(True)
                    self.host.view_3d_manager.setup_3d_hover()
                else:
                    self.host.clear_3d_view()
                    self.host.ui_manager.enable_3d_features(False)
            except (RuntimeError, ValueError, TypeError) as e:
                logging.error(f"Could not load 3D model from state data: {e}")
                self.host.update_status_message(f"Error loading 3D model: {e}", 5000)
                self.host.set_current_molecule(None)
                self.host.ui_manager.enable_3d_features(False)

        else:
            self.host.clear_3d_view()
            self.host.init_manager.analysis_action.setEnabled(False)
            self.host.init_manager.optimize_3d_button.setEnabled(False)
            # Disable 3D features
            self.host.ui_manager.enable_3d_features(False)

        self.host.edit_actions_manager.update_implicit_hydrogens()
        self.host.view_3d_manager.update_chiral_labels()

        if loaded_data.get("is_3d_viewer_mode", False):
            self.host.ui_manager.enter_3d_viewer_mode()
            self.host.statusBar().showMessage("Project loaded in 3D Viewer Mode.")
        else:
            self.host.ui_manager.restore_ui_for_editing()
            # Enable 3D edit features even in 2D editor mode if 3D molecule exists
            if (
                self.host.current_mol
                and self.host.current_mol.GetNumAtoms() > 0
            ):
                self.host.ui_manager.enable_3d_edit_actions(True)

        # Update labels after undo/redo
        self.host.edit_3d_manager.update_2d_measurement_labels()

    def push_undo_state(self) -> None:
        self.host.edit_actions_manager.push_undo_state()

    def update_window_title(self) -> None:
        """Update window title to reflect save state."""
        base_title = f"MoleditPy Ver. {VERSION}"
        if self.host.init_manager.current_file_path:
            filename = os.path.basename(self.host.init_manager.current_file_path)
            title = f"{filename} - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        else:
            # Handle as Untitled
            title = f"Untitled - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        self.host.setWindowTitle(title)

    def check_unsaved_changes(self) -> bool:
        """Check for unsaved changes and show warning."""
        if not self.has_unsaved_changes:
            return True  # Saved or no changes

        if (
            not self.data.atoms
            and self.host.current_mol is None
        ):
            return True  # Empty document

        reply = QMessageBox.question(
            self.host,
            "Unsaved Changes",
            "You have unsaved changes. Do you want to save them?",
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )

        if reply == QMessageBox.StandardButton.Yes:
            # 'Save As' if not PMEPRJ
            file_path = self.host.init_manager.current_file_path
            if not file_path or not file_path.lower().endswith(".pmeprj"):
                self.host.io_manager.save_project_as()
            else:
                self.host.io_manager.save_project()
            return (
                not self.has_unsaved_changes
            )  # Return True only if save was successful
        elif reply == QMessageBox.StandardButton.No:
            return True  # Continue without saving
        else:
            return False  # Cancel

    def reset_undo_stack(self) -> None:
        self.host.reset_undo_redo_stacks()

    def update_realtime_info(self) -> None:
        """Show molecular info in status bar."""
        if not self.data.atoms:
            self.host.update_formula_label("")
            return

        if self.data:
            try:
                mol = self.data.to_rdkit_mol()
                if mol:
                    mol_with_hs = Chem.AddHs(mol)
                    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                    num_atoms = mol_with_hs.GetNumAtoms()
                    self.host.update_formula_label(
                        f"Formula: {mol_formula}   |   Atoms: {num_atoms}"
                    )
            except (RuntimeError, TypeError, ValueError) as e:
                logging.debug(
                    f"Molecular info update suppressed for unstable structure: {e}"
                )
                self.host.update_formula_label("Invalid structure")

    def create_json_data(self) -> Dict[str, Any]:
        """Convert current state to PMEJSON."""
        # Metadata
        json_data: Dict[str, Any] = {
            "format": "PME Project",
            "version": "1.0",
            "application": "MoleditPy",
            "application_version": VERSION,
            "created": str(QDateTime.currentDateTime().toString(Qt.DateFormat.ISODate)),
            "is_3d_viewer_mode": not self.host.ui_manager.is_2d_editable,
        }

        # 2D data
        if self.data.atoms:
            atoms_2d = []
            for atom_id, data in self.data.atoms.items():
                pos = data["pos"]
                atom_data = {
                    "id": atom_id,
                    "symbol": data["symbol"],
                    "x": pos[0],
                    "y": pos[1],
                    "charge": data.get("charge", 0),
                    "radical": data.get("radical", 0),
                }
                atoms_2d.append(atom_data)

            bonds_2d = []
            for (
                atom1_id,
                atom2_id,
            ), bond_data in self.data.bonds.items():
                bond_info = {
                    "atom1": atom1_id,
                    "atom2": atom2_id,
                    "order": bond_data["order"],
                    "stereo": bond_data.get("stereo", 0),
                }
                bonds_2d.append(bond_info)

            json_data["2d_structure"] = {
                "atoms": atoms_2d,
                "bonds": bonds_2d,
                "next_atom_id": self.data.next_atom_id,
            }

        # 3D data
        if (
            self.host.current_mol
            and self.host.current_mol.GetNumConformers() > 0
        ):
            try:
                # Save as Base64 encoded binary
                mol_binary = self.host.current_mol.ToBinary()
                mol_base64 = base64.b64encode(mol_binary).decode("ascii")

                # Extract 3D coordinates
                atoms_3d = []
                if self.host.current_mol.GetNumConformers() > 0:
                    conf = self.host.current_mol.GetConformer()
                    for i in range(self.host.current_mol.GetNumAtoms()):
                        atom = self.host.current_mol.GetAtomWithIdx(i)
                        pos = conf.GetAtomPosition(i)

                        # Try to preserve original editor atom ID (if present) so it can be
                        # restored when loading PMEPRJ files. RDKit atom properties may
                        # contain _original_atom_id when the molecule was created from
                        # the editor's 2D structure.
                        original_id = None
                        try:
                            if atom.HasProp("_original_atom_id"):
                                original_id = atom.GetIntProp("_original_atom_id")
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            # Skip property access if atom data is inconsistent
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
                            "original_id": original_id,
                        }
                        atoms_3d.append(atom_3d)

                # Extract bonds
                bonds_3d = []
                for bond in self.host.current_mol.GetBonds():
                    bond_3d = {
                        "atom1": bond.GetBeginAtomIdx(),
                        "atom2": bond.GetEndAtomIdx(),
                        "order": int(bond.GetBondType()),
                        "is_aromatic": bond.GetIsAromatic(),
                        "stereo": int(bond.GetStereo()),
                    }
                    bonds_3d.append(bond_3d)

                try:
                    json_safe_constraints = _serialize_constraints(
                        self.host.edit_3d_manager.constraints_3d
                    )
                except (AttributeError, TypeError, KeyError, ValueError):
                    json_safe_constraints = []

                json_data["3d_structure"] = {
                    "mol_binary_base64": mol_base64,
                    "atoms": atoms_3d,
                    "bonds": bonds_3d,
                    "num_conformers": self.host.current_mol.GetNumConformers(),
                    "constraints_3d": json_safe_constraints,
                }

                # Molecular info
                json_data["molecular_info"] = {
                    "num_atoms": self.host.current_mol.GetNumAtoms(),
                    "num_bonds": self.host.current_mol.GetNumBonds(),
                    "molecular_weight": Descriptors.MolWt(
                        self.host.current_mol
                    ),
                    "formula": rdMolDescriptors.CalcMolFormula(
                        self.host.current_mol
                    ),
                }

                # Identifiers (SMILES/InChI)
                try:
                    json_data["identifiers"] = {
                        "smiles": Chem.MolToSmiles(
                            self.host.current_mol
                        ),
                        "canonical_smiles": Chem.MolToSmiles(
                            self.host.current_mol, canonical=True
                        ),
                    }

                    # Attempt InChI generation
                    try:
                        inchi = Chem.MolToInchi(self.host.current_mol)
                        inchi_key = Chem.MolToInchiKey(
                            self.host.current_mol
                        )
                        json_data["identifiers"]["inchi"] = inchi
                        json_data["identifiers"]["inchi_key"] = inchi_key
                    except (AttributeError, RuntimeError, TypeError, ValueError):
                        # Suppress InChI generation errors during project save if RDKit lacks InChI support
                        # Safe defensive fallback catching AttributeError, RuntimeError, TypeError, ValueError
                        pass

                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    logging.warning("Could not generate molecular identifiers: %s", e)

            except (AttributeError, RuntimeError, ValueError) as e:
                logging.warning("Could not process 3D molecular data: %s", e)
        else:
            # Record if no 3D data
            json_data["3d_structure"] = None
            json_data["note"] = (
                "No 3D structure available. Generate 3D coordinates first."
            )

        # Record the last-successful optimization method (if any)
        try:
            json_data["last_successful_optimization_method"] = getattr(
                self.host.compute_manager, "last_successful_optimization_method", None
            )
        except (AttributeError, RuntimeError, TypeError):
            json_data["last_successful_optimization_method"] = None

        # Plugin State Persistence (Phase 3)
        # Start with preserved data from missing plugins
        plugin_data = (
            self._preserved_plugin_data.copy() if self._preserved_plugin_data else {}
        )

        pm = getattr(self.host, "plugin_manager", None)
        if pm and pm.save_handlers:
            for name, callback in pm.save_handlers.items():
                if callable(callback):
                    try:
                        p_state = callback()
                        plugin_data[name] = p_state
                    except (RuntimeError, TypeError, ValueError, AttributeError) as e:
                        logging.error(f"Error saving state for plugin {name}: {e}")

        if plugin_data:
            json_data["plugins"] = plugin_data

        return json_data

    def load_from_json_data(self, json_data: Dict[str, Any]) -> None:
        """Restore state from JSON."""
        self.dragged_atom_info = None
        self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)
        self.host.ui_manager.enable_3d_edit_actions(False)
        self.host.ui_manager.enable_3d_features(False)

        # 3D viewer mode
        is_3d_mode = json_data.get("is_3d_viewer_mode", False)
        # Restore last successful optimization method if present in file
        try:
            method = json_data.get("last_successful_optimization_method", None)
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.host.set_last_successful_optimization_method(method)
        except (AttributeError, RuntimeError, TypeError):
            # Safe defensive fallback catching AttributeError, RuntimeError, TypeError
            pass

        # Plugin State Restoration (Phase 3)
        self._preserved_plugin_data = {}  # Reset preserved data on new load
        pm = getattr(self.host, "plugin_manager", None)
        if "plugins" in json_data:
            plugin_data = json_data["plugins"]
            if isinstance(plugin_data, dict):
                for name, p_state in plugin_data.items():
                    load_hand = None
                    if pm:
                        load_hand = pm.load_handlers.get(name)

                    if load_hand and callable(load_hand):
                        try:
                            load_hand(p_state)
                        except (RuntimeError, ValueError, TypeError) as e:
                            logging.error(f"Error loading state for plugin {name}: {e}")
                    else:
                        # No handler found (plugin disabled or missing)
                        # Preserve data so it's not lost on next save
                        self._preserved_plugin_data[name] = p_state

        # Restore 2D data
        if "2d_structure" in json_data:
            structure_2d = json_data["2d_structure"]
            atoms_2d = structure_2d.get("atoms", [])
            bonds_2d = structure_2d.get("bonds", [])

            self.host.scene.restore_atoms_and_bonds_from_json(atoms_2d, bonds_2d)
            self.data.next_atom_id = structure_2d.get(
                "next_atom_id",
                max([atom["id"] for atom in atoms_2d]) + 1 if atoms_2d else 0,
            )
        # Restore 3D data
        structure_3d = json_data.get("3d_structure")
        if isinstance(structure_3d, dict):
            loaded_constraints = structure_3d.get("constraints_3d", [])
            self.host.set_constraints_3d(
                _deserialize_constraints(loaded_constraints)
                if isinstance(loaded_constraints, list)
                else []
            )
            try:
                # Restore binary data
                mol_base64 = structure_3d.get("mol_binary_base64")
                if mol_base64:
                    mol_binary = base64.b64decode(mol_base64.encode("ascii"))
                    self.host.set_current_molecule(Chem.Mol(mol_binary))
                    if self.host.current_mol:
                        # Set 3D coordinates
                        if self.host.current_mol.GetNumConformers() > 0:
                            atoms_3d = structure_3d.get("atoms", [])
                            # Ensure numpy array size matches atoms in file
                            num_atoms_file = len(atoms_3d)
                            if num_atoms_file > 0:
                                positions_3d = np.zeros((num_atoms_file, 3))
                                for atom_data in atoms_3d:
                                    idx = atom_data.get("index", -1)
                                    if 0 <= idx < num_atoms_file:
                                        positions_3d[idx] = [
                                            atom_data.get("x", 0.0),
                                            atom_data.get("y", 0.0),
                                            atom_data.get("z", 0.0),
                                        ]
                                        # Restore original editor atom id into RDKit atom property
                                        original_id = atom_data.get("original_id")
                                        if original_id is not None:
                                            try:
                                                rd_atom = self.host.current_mol.GetAtomWithIdx(
                                                    idx
                                                )
                                                if rd_atom:
                                                    rd_atom.SetIntProp(
                                                        "_original_atom_id",
                                                        int(original_id),
                                                    )
                                            except (
                                                RuntimeError,
                                                ValueError,
                                                TypeError,
                                                IndexError,
                                            ):  # [RDKIT GUARD] RDKit atom property assignment may fail on invalid ids; skip silently.
                                                pass
                                self.host.set_3d_atom_positions(positions_3d)

                            # Build mapping
                            try:
                                self.host.compute_manager.create_atom_id_mapping()
                                self.host.view_3d_manager.update_atom_id_menu_text()
                                self.host.view_3d_manager.update_atom_id_menu_state()
                            except (RuntimeError, TypeError, AttributeError):
                                # Safe defensive fallback catching RuntimeError, TypeError, AttributeError
                                pass

                        # Always show 3D if 3D molecule exists
                        self.host.view_3d_manager.draw_molecule_3d(
                            self.host.current_mol
                        )

                        # Switch UI in Viewer mode
                        if is_3d_mode:
                            self.host.ui_manager.enter_3d_viewer_mode()
                        else:
                            self.host.set_is_2d_editable(True)

                        if self.host.plotter:
                            self.host.plotter.reset_camera()

                        # Enable 3D-related UI
                        try:
                            self.host.ui_manager.enable_3d_edit_actions(True)
                            self.host.ui_manager.enable_3d_features(True)
                        except (RuntimeError, TypeError, AttributeError):
                            # Safe defensive fallback catching RuntimeError, TypeError, AttributeError
                            pass
            except (RuntimeError, ValueError, TypeError, binascii.Error) as e:
                logging.error(f"Could not restore 3D molecular data: {e}")
                self.host.set_current_molecule(None)
