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
# Module separated from MainWindow (main_window.py)
# Functional class: MainWindowAppState
"""

import base64
import copy
import logging
import os
import numpy as np

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# PyQt6 Modules
from PyQt6.QtCore import QDateTime, QPointF, Qt
from PyQt6.QtWidgets import QMessageBox

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .atom_item import AtomItem
    from .bond_item import BondItem
    from .constants import VERSION
except ImportError:
    # Fallback to absolute imports for script-style execution
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem
    from modules.constants import VERSION


# --- Class Definition ---
class MainWindowAppState:
    """Mixin class separated from main_window.py"""

    _cls = None

    def __init__(self):
        """Initialize class. 'self' is MainWindow instance."""
        pass

    def get_current_state(self):
        atoms = {
            atom_id: {
                "symbol": data["symbol"],
                "pos": (data["item"].pos().x(), data["item"].pos().y()),
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
            "_next_atom_id": self.data._next_atom_id,
        }

        state["version"] = VERSION

        if self.current_mol:
            state["mol_3d"] = self.current_mol.ToBinary()
            # RDKit binary serialization does not preserve custom properties like _original_atom_id.
            # We store them separately to ensure we can restore the 2D-3D link after undo/redo.
            mol_3d_atom_ids = []
            for i in range(self.current_mol.GetNumAtoms()):
                atom = self.current_mol.GetAtomWithIdx(i)
                if atom and atom.HasProp("_original_atom_id"):
                    try:
                        mol_3d_atom_ids.append(atom.GetIntProp("_original_atom_id"))
                    except (RuntimeError, ValueError, TypeError):
                        mol_3d_atom_ids.append(None)
                else:
                    mol_3d_atom_ids.append(None)
            state["mol_3d_atom_ids"] = mol_3d_atom_ids

        state["is_3d_viewer_mode"] = not self.is_2d_editable

        json_safe_constraints = []
        constraints = getattr(self, "constraints_3d", [])
        for const in constraints:
            # (Type, (Idx...), Value, Force) -> [Type, [Idx...], Value, Force]
            if len(const) == 4:
                json_safe_constraints.append(
                    [const[0], list(const[1]), const[2], const[3]]
                )
            else:
                # Backward compatibility
                json_safe_constraints.append(
                    [const[0], list(const[1]), const[2], 1.0e5]
                )
        state["constraints_3d"] = json_safe_constraints

        return state

    def set_state_from_data(self, state_data):
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)

        loaded_data = copy.deepcopy(state_data)

        # Get file version (default '0.0.0')
        file_version_str = loaded_data.get("version", "0.0.0")

        def parse_v(v_str):
            if not v_str or not isinstance(v_str, str):
                return (0, 0, 0)
            try:
                parts = []
                import re

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
            if hasattr(self, "warning_message_box"):  # Use helper if available
                self.warning_message_box(
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).",
                )
            else:
                QMessageBox.warning(
                    self,
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).\n\n"
                    f"Your current version is {VERSION}.\n\n"
                    "Some features may not load or work correctly.",
                )

        raw_atoms = loaded_data.get("atoms", {})
        raw_bonds = loaded_data.get("bonds", {})

        # Restore constraints
        loaded_constraints = loaded_data.get("constraints_3d", [])
        if loaded_constraints and isinstance(loaded_constraints, list):
            self.constraints_3d = []
            for const in loaded_constraints:
                if isinstance(const, (list, tuple)):
                    try:
                        if len(const) == 4:
                            self.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], const[3])
                            )
                        elif len(const) == 3:
                            self.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], 1.0e5)
                            )
                    except (TypeError, ValueError, IndexError) as e:
                        logging.debug(f"Failed to parse constraint {const}: {e}")
        else:
            self.constraints_3d = []

        for atom_id, data in raw_atoms.items():
            pos = QPointF(data["pos"][0], data["pos"][1])
            charge = data.get("charge", 0)
            radical = data.get("radical", 0)
            # Pass radical to AtomItem
            atom_item = AtomItem(
                atom_id, data["symbol"], pos, charge=charge, radical=radical
            )
            # Store radical in data
            self.data.atoms[atom_id] = {
                "symbol": data["symbol"],
                "pos": pos,
                "item": atom_item,
                "charge": charge,
                "radical": radical,
            }
            self.scene.addItem(atom_item)

        self.data._next_atom_id = loaded_data.get(
            "_next_atom_id", max(self.data.atoms.keys()) + 1 if self.data.atoms else 0
        )

        for key_tuple, data in raw_bonds.items():
            id1, id2 = key_tuple
            if id1 in self.data.atoms and id2 in self.data.atoms:
                atom1_item = self.data.atoms[id1]["item"]
                atom2_item = self.data.atoms[id2]["item"]
                bond_item = BondItem(
                    atom1_item, atom2_item, data.get("order", 1), data.get("stereo", 0)
                )
                self.data.bonds[key_tuple] = {
                    "order": data.get("order", 1),
                    "stereo": data.get("stereo", 0),
                    "item": bond_item,
                }
                atom1_item.bonds.append(bond_item)
                atom2_item.bonds.append(bond_item)
                self.scene.addItem(bond_item)

        for atom_data in self.data.atoms.values():
            if atom_data["item"]:
                atom_data["item"].update_style()
        self.scene.update_all_items()
        mol_3d_data = loaded_data.get("mol_3d")
        if mol_3d_data is not None:
            try:
                self.current_mol = Chem.Mol(mol_3d_data)
                # Debug: check if 3D structure is valid
                if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                    # Restore _original_atom_id if present in saved state
                    mol_3d_atom_ids = loaded_data.get("mol_3d_atom_ids")
                    if (
                        mol_3d_atom_ids
                        and len(mol_3d_atom_ids) == self.current_mol.GetNumAtoms()
                    ):
                        for i, aid in enumerate(mol_3d_atom_ids):
                            if aid is not None:
                                rd_atom = self.current_mol.GetAtomWithIdx(i)
                                if rd_atom:
                                    try:
                                        rd_atom.SetIntProp(
                                            "_original_atom_id", int(aid)
                                        )
                                    except (RuntimeError, ValueError, TypeError):
                                        pass

                    # Sync 2D atoms with 3D actors
                    if hasattr(self, "create_atom_id_mapping"):
                        try:
                            self.create_atom_id_mapping()
                            if hasattr(self, "update_atom_id_menu_text"):
                                self.update_atom_id_menu_text()
                            if hasattr(self, "update_atom_id_menu_state"):
                                self.update_atom_id_menu_state()
                        except Exception as e:
                            logging.debug(
                                f"Partial failure during ID mapping restoration: {e}"
                            )

                    # draw_molecule_3d will use restored IDs
                    if hasattr(self, "draw_molecule_3d"):
                        self.draw_molecule_3d(self.current_mol)
                    if (
                        hasattr(self, "plotter")
                        and self.plotter
                        and hasattr(self.plotter, "reset_camera")
                    ):
                        self.plotter.reset_camera()

                    self._enable_3d_features(True)
                    if hasattr(self, "setup_3d_hover"):
                        self.setup_3d_hover()
                else:
                    self.current_mol = None
                    if hasattr(self, "plotter") and self.plotter:
                        self.plotter.clear()
                    self._enable_3d_features(False)
            except (RuntimeError, ValueError, TypeError) as e:
                logging.error(f"Could not load 3D model from state data: {e}")
                if hasattr(self, "statusBar") and self.statusBar():
                    self.statusBar().showMessage(f"Error loading 3D model: {e}", 5000)
                self.current_mol = None
                self._enable_3d_features(False)

        else:
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)
            self.optimize_3d_button.setEnabled(False)
            # Disable 3D features
            self._enable_3d_features(False)

        self.update_implicit_hydrogens()
        self.update_chiral_labels()

        if loaded_data.get("is_3d_viewer_mode", False):
            self._enter_3d_viewer_ui_mode()
            self.statusBar().showMessage("Project loaded in 3D Viewer Mode.")
        else:
            self.restore_ui_for_editing()
            # Enable 3D edit features even in 2D editor mode if 3D molecule exists
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                self._enable_3d_edit_actions(True)

        # Update labels after undo/redo
        self.update_2d_measurement_labels()

    def push_undo_state(self):
        if self._is_restoring_state:
            return

        current_state_for_comparison = {
            "atoms": {
                k: (
                    v["symbol"],
                    v["item"].pos().x(),
                    v["item"].pos().y(),
                    v.get("charge", 0),
                    v.get("radical", 0),
                )
                for k, v in self.data.atoms.items()
            },
            "bonds": {
                k: (v["order"], v.get("stereo", 0)) for k, v in self.data.bonds.items()
            },
            "_next_atom_id": self.data._next_atom_id,
            "mol_3d": self.current_mol.ToBinary() if self.current_mol else None,
            "mol_3d_atom_ids": [
                (
                    a.GetIntProp("_original_atom_id")
                    if a and a.HasProp("_original_atom_id")
                    else None
                )
                for a in self.current_mol.GetAtoms()
            ]
            if self.current_mol
            else None,
        }

        last_state_for_comparison = None
        if self.undo_stack:
            last_state = self.undo_stack[-1]
            last_atoms = last_state.get("atoms", {})
            last_bonds = last_state.get("bonds", {})
            last_state_for_comparison = {
                "atoms": {
                    k: (
                        v["symbol"],
                        v["pos"][0],
                        v["pos"][1],
                        v.get("charge", 0),
                        v.get("radical", 0),
                    )
                    for k, v in last_atoms.items()
                },
                "bonds": {
                    k: (v["order"], v.get("stereo", 0)) for k, v in last_bonds.items()
                },
                "_next_atom_id": last_state.get("_next_atom_id"),
                "mol_3d": last_state.get("mol_3d", None),
                "mol_3d_atom_ids": last_state.get("mol_3d_atom_ids", None),
            }

        if (
            not last_state_for_comparison
            or current_state_for_comparison != last_state_for_comparison
        ):
            # Deepcopy state to ensure saved states are immutable and not affected
            # by later modifications to objects referenced from the state.
            state = copy.deepcopy(self.get_current_state())
            self.undo_stack.append(state)

            self.redo_stack.clear()
            # Record changes after initialization
            if self.initialization_complete:
                self.has_unsaved_changes = True
                self.update_window_title()

        self.update_implicit_hydrogens()
        self.update_realtime_info()
        self.update_undo_redo_actions()

    def update_window_title(self):
        """Update window title to reflect save state."""
        base_title = f"MoleditPy Ver. {VERSION}"
        if self.current_file_path:
            filename = os.path.basename(self.current_file_path)
            title = f"{filename} - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        else:
            # Handle as Untitled
            title = f"Untitled - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        self.setWindowTitle(title)

    def check_unsaved_changes(self):
        """Check for unsaved changes and show warning."""
        if not self.has_unsaved_changes:
            return True  # Saved or no changes

        if not self.data.atoms and self.current_mol is None:
            return True  # Empty document

        reply = QMessageBox.question(
            self,
            "Unsaved Changes",
            "You have unsaved changes. Do you want to save them?",
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )

        if reply == QMessageBox.StandardButton.Yes:
            # 'Save As' if not PMEPRJ
            file_path = self.current_file_path
            if not file_path or not file_path.lower().endswith(".pmeprj"):
                self.save_project_as()
            else:
                self.save_project()
            return (
                not self.has_unsaved_changes
            )  # Return True only if save was successful
        elif reply == QMessageBox.StandardButton.No:
            return True  # Continue without saving
        else:
            return False  # Cancel

    def reset_undo_stack(self):
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.push_undo_state()

    def undo(self):
        if len(self.undo_stack) > 1:
            self.redo_stack.append(self.undo_stack.pop())
            state = self.undo_stack[-1]
            self._is_restoring_state = True
            try:
                self.set_state_from_data(state)
            finally:
                self._is_restoring_state = False

            # Re-evaluate menu states based on 3D structure after Undo
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D structure exists: enable 3D edit features
                self._enable_3d_edit_actions(True)
            else:
                # No 3D structure: disable 3D edit features
                self._enable_3d_edit_actions(False)

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

            # Re-evaluate menu states based on 3D structure after Redo
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D structure exists: enable 3D edit features
                self._enable_3d_edit_actions(True)
            else:
                # No 3D structure: disable 3D edit features
                self._enable_3d_edit_actions(False)

        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus()

    def update_undo_redo_actions(self):
        self.undo_action.setEnabled(len(self.undo_stack) > 1)
        self.redo_action.setEnabled(len(self.redo_stack) > 0)

    def update_realtime_info(self):
        """Show molecular info in status bar."""
        if not self.data.atoms:
            self.formula_label.setText("")  # Clear label if no atoms
            return

        if hasattr(self, "data") and self.data and hasattr(self.data, "to_rdkit_mol"):
            try:
                mol = self.data.to_rdkit_mol()
                if mol:
                    # Generate mol with explicit Hs
                    mol_with_hs = Chem.AddHs(mol)
                    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                    # Get atom count with Hs
                    num_atoms = mol_with_hs.GetNumAtoms()
                    # Update label text
                    if hasattr(self, "formula_label") and self.formula_label:
                        self.formula_label.setText(
                            f"Formula: {mol_formula}   |   Atoms: {num_atoms}"
                        )
            except (RuntimeError, TypeError, ValueError) as e:
                logging.debug(
                    f"Molecular info update suppressed for unstable structure: {e}"
                )
                if hasattr(self, "formula_label") and self.formula_label:
                    self.formula_label.setText("Invalid structure")

    def create_json_data(self):
        """Convert current state to PMEJSON."""
        # Metadata
        json_data = {
            "format": "PME Project",
            "version": "1.0",
            "application": "MoleditPy",
            "application_version": VERSION,
            "created": str(QDateTime.currentDateTime().toString(Qt.DateFormat.ISODate)),
            "is_3d_viewer_mode": not self.is_2d_editable,
        }

        # 2D data
        if self.data.atoms:
            atoms_2d = []
            for atom_id, data in self.data.atoms.items():
                pos = data["item"].pos()
                atom_data = {
                    "id": atom_id,
                    "symbol": data["symbol"],
                    "x": pos.x(),
                    "y": pos.y(),
                    "charge": data.get("charge", 0),
                    "radical": data.get("radical", 0),
                }
                atoms_2d.append(atom_data)

            bonds_2d = []
            for (atom1_id, atom2_id), bond_data in self.data.bonds.items():
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
                "next_atom_id": self.data._next_atom_id,
            }

        # 3D data
        if self.current_mol and self.current_mol.GetNumConformers() > 0:
            try:
                # Save as Base64 encoded binary
                mol_binary = self.current_mol.ToBinary()
                mol_base64 = base64.b64encode(mol_binary).decode("ascii")

                # Extract 3D coordinates
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
                for bond in self.current_mol.GetBonds():
                    bond_3d = {
                        "atom1": bond.GetBeginAtomIdx(),
                        "atom2": bond.GetEndAtomIdx(),
                        "order": int(bond.GetBondType()),
                        "is_aromatic": bond.GetIsAromatic(),
                        "stereo": int(bond.GetStereo()),
                    }
                    bonds_3d.append(bond_3d)

                # Convert constraints to JSON compatible format
                json_safe_constraints = []
                try:
                    for const in self.constraints_3d:
                        if len(const) == 4:
                            json_safe_constraints.append(
                                [const[0], list(const[1]), const[2], const[3]]
                            )
                        else:
                            json_safe_constraints.append(
                                [const[0], list(const[1]), const[2], 1.0e5]
                            )
                except (AttributeError, TypeError, KeyError, ValueError):
                    # Reset constraints if serialization fails
                    json_safe_constraints = []

                json_data["3d_structure"] = {
                    "mol_binary_base64": mol_base64,
                    "atoms": atoms_3d,
                    "bonds": bonds_3d,
                    "num_conformers": self.current_mol.GetNumConformers(),
                    "constraints_3d": json_safe_constraints,
                }

                # Molecular info
                json_data["molecular_info"] = {
                    "num_atoms": self.current_mol.GetNumAtoms(),
                    "num_bonds": self.current_mol.GetNumBonds(),
                    "molecular_weight": Descriptors.MolWt(self.current_mol),
                    "formula": rdMolDescriptors.CalcMolFormula(self.current_mol),
                }

                # Identifiers (SMILES/InChI)
                try:
                    json_data["identifiers"] = {
                        "smiles": Chem.MolToSmiles(self.current_mol),
                        "canonical_smiles": Chem.MolToSmiles(
                            self.current_mol, canonical=True
                        ),
                    }

                    # Attempt InChI generation
                    try:
                        inchi = Chem.MolToInchi(self.current_mol)
                        inchi_key = Chem.MolToInchiKey(self.current_mol)
                        json_data["identifiers"]["inchi"] = inchi
                        json_data["identifiers"]["inchi_key"] = inchi_key
                    except (AttributeError, RuntimeError, TypeError, ValueError):
                        # Suppress InChI generation errors during project save if RDKit lacks InChI support
                        pass

                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    print(f"Warning: Could not generate molecular identifiers: {e}")

            except (AttributeError, RuntimeError, ValueError) as e:
                print(f"Warning: Could not process 3D molecular data: {e}")
        else:
            # Record if no 3D data
            json_data["3d_structure"] = None
            json_data["note"] = (
                "No 3D structure available. Generate 3D coordinates first."
            )

        # Record the last-successful optimization method (if any)
        # This is a convenience field so saved projects remember which
        # optimizer variant was last used (e.g. "MMFF94s", "MMFF94", "UFF").
        try:
            json_data["last_successful_optimization_method"] = getattr(
                self, "last_successful_optimization_method", None
            )
        except (AttributeError, RuntimeError, TypeError):
            json_data["last_successful_optimization_method"] = None

        # Plugin State Persistence (Phase 3)
        # Start with preserved data from missing plugins
        plugin_data = (
            self._preserved_plugin_data.copy() if self._preserved_plugin_data else {}
        )

        pm = getattr(self, "plugin_manager", None)
        if pm and hasattr(pm, "save_handlers") and pm.save_handlers:
            for name, callback in pm.save_handlers.items():
                if callable(callback):
                    try:
                        p_state = callback()
                        plugin_data[name] = p_state
                    except Exception as e:
                        logging.error(f"Error saving state for plugin {name}: {e}")

        if plugin_data:
            json_data["plugins"] = plugin_data

        return json_data

    def load_from_json_data(self, json_data):
        """Restore state from JSON."""
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)
        self._enable_3d_edit_actions(False)
        self._enable_3d_features(False)

        # 3D viewer mode
        is_3d_mode = json_data.get("is_3d_viewer_mode", False)
        # Restore last successful optimization method if present in file
        try:
            self.last_successful_optimization_method = json_data.get(
                "last_successful_optimization_method", None
            )
        except (AttributeError, RuntimeError, TypeError):
            self.last_successful_optimization_method = None

        # Plugin State Restoration (Phase 3)
        self._preserved_plugin_data = {}  # Reset preserved data on new load
        pm = getattr(self, "plugin_manager", None)
        if "plugins" in json_data:
            plugin_data = json_data["plugins"]
            if isinstance(plugin_data, dict):
                for name, p_state in plugin_data.items():
                    load_hand = None
                    if pm and hasattr(pm, "load_handlers"):
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

            # Restore atoms
            for atom_data in atoms_2d:
                atom_id = atom_data["id"]
                symbol = atom_data["symbol"]
                pos = QPointF(atom_data["x"], atom_data["y"])
                charge = atom_data.get("charge", 0)
                radical = atom_data.get("radical", 0)

                atom_item = AtomItem(
                    atom_id, symbol, pos, charge=charge, radical=radical
                )
                self.data.atoms[atom_id] = {
                    "symbol": symbol,
                    "pos": pos,
                    "item": atom_item,
                    "charge": charge,
                    "radical": radical,
                }
                self.scene.addItem(atom_item)

            # Restore next_atom_id
            self.data._next_atom_id = structure_2d.get(
                "next_atom_id",
                max([atom["id"] for atom in atoms_2d]) + 1 if atoms_2d else 0,
            )

            # Restore bonds
            for bond_data in bonds_2d:
                atom1_id = bond_data["atom1"]
                atom2_id = bond_data["atom2"]

                if atom1_id in self.data.atoms and atom2_id in self.data.atoms:
                    atom1_item = self.data.atoms[atom1_id]["item"]
                    atom2_item = self.data.atoms[atom2_id]["item"]

                    bond_order = bond_data["order"]
                    stereo = bond_data.get("stereo", 0)

                    bond_item = BondItem(
                        atom1_item, atom2_item, bond_order, stereo=stereo
                    )
                    # Add to bond list (used for C visibility)
                    atom1_item.bonds.append(bond_item)
                    atom2_item.bonds.append(bond_item)

                    self.data.bonds[(atom1_id, atom2_id)] = {
                        "order": bond_order,
                        "item": bond_item,
                        "stereo": stereo,
                    }
                    self.scene.addItem(bond_item)

            # Update all AtomItem styles
            for atom in self.data.atoms.values():
                atom["item"].update_style()
            self.scene.update_all_items()
        # Restore 3D data
        structure_3d = json_data.get("3d_structure")
        if isinstance(structure_3d, dict):
            # Restore constraints
            loaded_constraints = structure_3d.get("constraints_3d", [])
            if loaded_constraints and isinstance(loaded_constraints, list):
                self.constraints_3d = []
                for const in loaded_constraints:
                    if isinstance(const, (list, tuple)):
                        try:
                            if len(const) == 4:
                                self.constraints_3d.append(
                                    (const[0], tuple(const[1]), const[2], const[3])
                                )
                            elif len(const) == 3:
                                self.constraints_3d.append(
                                    (const[0], tuple(const[1]), const[2], 1.0e5)
                                )
                        except (TypeError, ValueError, IndexError):
                            pass
            else:
                self.constraints_3d = []
            try:
                # Restore binary data
                mol_base64 = structure_3d.get("mol_binary_base64")
                if mol_base64:
                    mol_binary = base64.b64decode(mol_base64.encode("ascii"))
                    self.current_mol = Chem.Mol(mol_binary)
                    if self.current_mol:
                        # Set 3D coordinates
                        if self.current_mol.GetNumConformers() > 0:
                            conf = self.current_mol.GetConformer()
                            atoms_3d = structure_3d.get("atoms", [])
                            # Ensure numpy array size matches atoms in file
                            num_atoms_file = len(atoms_3d)
                            if num_atoms_file > 0:
                                self.atom_positions_3d = np.zeros((num_atoms_file, 3))
                                for atom_data in atoms_3d:
                                    idx = atom_data.get("index", -1)
                                    if 0 <= idx < num_atoms_file:
                                        self.atom_positions_3d[idx] = [
                                            atom_data.get("x", 0.0),
                                            atom_data.get("y", 0.0),
                                            atom_data.get("z", 0.0),
                                        ]
                                        # Restore original editor atom id into RDKit atom property
                                        original_id = atom_data.get("original_id")
                                        if original_id is not None:
                                            try:
                                                rd_atom = (
                                                    self.current_mol.GetAtomWithIdx(idx)
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
                                            ):
                                                pass

                            # Build mapping
                            if hasattr(self, "create_atom_id_mapping"):
                                try:
                                    self.create_atom_id_mapping()
                                    if hasattr(self, "update_atom_id_menu_text"):
                                        self.update_atom_id_menu_text()
                                    if hasattr(self, "update_atom_id_menu_state"):
                                        self.update_atom_id_menu_state()
                                except (RuntimeError, TypeError, AttributeError):
                                    pass

                        # Always show 3D if 3D molecule exists
                        if hasattr(self, "draw_molecule_3d"):
                            self.draw_molecule_3d(self.current_mol)

                        # Switch UI in Viewer mode
                        if is_3d_mode and hasattr(self, "_enter_3d_viewer_ui_mode"):
                            self._enter_3d_viewer_ui_mode()
                        else:
                            self.is_2d_editable = True

                        if (
                            hasattr(self, "plotter")
                            and self.plotter
                            and hasattr(self.plotter, "reset_camera")
                        ):
                            self.plotter.reset_camera()

                        # Enable 3D-related UI
                        try:
                            self._enable_3d_edit_actions(True)
                            self._enable_3d_features(True)
                        except (RuntimeError, TypeError, AttributeError):
                            pass
            except (RuntimeError, ValueError, TypeError, base64.binascii.Error) as e:
                logging.error(f"Could not restore 3D molecular data: {e}")
                self.current_mol = None


MainWindowAppState._cls = MainWindowAppState
