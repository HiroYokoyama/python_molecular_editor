#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
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
    from ..utils.constants import VERSION
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.atom_item import AtomItem
    from moleditpy.ui.bond_item import BondItem
    from moleditpy.utils.constants import VERSION


# --- Class Definition ---
class StateManager:
    _cls = None

    def __init__(self, host):
        self.host = host
        """Initialize class. 'self' is MainWindow instance."""
        pass

    def get_current_state(self):
        atoms = {
            atom_id: {
                "symbol": data["symbol"],
                "pos": data["pos"],  # Already a tuple (x, y)
                "charge": data.get("charge", 0),
                "radical": data.get("radical", 0),
            }
            for atom_id, data in self.host.data.atoms.items()
        }
        bonds = {
            key: {"order": data["order"], "stereo": data.get("stereo", 0)}
            for key, data in self.host.data.bonds.items()
        }
        state = {
            "atoms": atoms,
            "bonds": bonds,
            "_next_atom_id": self.host.data._next_atom_id,
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

        state["is_3d_viewer_mode"] = not self.host.is_2d_editable

        json_safe_constraints = []
        constraints = getattr(self.host, "constraints_3d", [])
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
        self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)

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
            if hasattr(self.host, "warning_message_box"):
                self.host.warning_message_box(
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).",
                )
            else:
                QMessageBox.warning(
                    self.host,
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
            self.host.constraints_3d = []
            for const in loaded_constraints:
                if isinstance(const, (list, tuple)):
                    try:
                        if len(const) == 4:
                            self.host.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], const[3])
                            )
                        elif len(const) == 3:
                            self.host.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], 1.0e5)
                            )
                    except (TypeError, ValueError, IndexError) as e:
                        logging.debug(f"Failed to parse constraint {const}: {e}")
        else:
            self.host.constraints_3d = []

        for atom_id, data in raw_atoms.items():
            raw_pos = tuple(data["pos"])
            pos_q = QPointF(raw_pos[0], raw_pos[1])
            charge = data.get("charge", 0)
            radical = data.get("radical", 0)
            # Pass QPointF to AtomItem for UI positioning
            atom_item = AtomItem(
                atom_id, data["symbol"], pos_q, charge=charge, radical=radical
            )
            # Store raw tuple in data
            self.host.data.atoms[atom_id] = {
                "symbol": data["symbol"],
                "pos": raw_pos,
                "item": atom_item,
                "charge": charge,
                "radical": radical,
            }
            self.host.scene.addItem(atom_item)

        self.host.data._next_atom_id = loaded_data.get(
            "_next_atom_id", max(self.host.data.atoms.keys()) + 1 if self.host.data.atoms else 0
        )

        for key_tuple, data in raw_bonds.items():
            id1, id2 = key_tuple
            if id1 in self.host.data.atoms and id2 in self.host.data.atoms:
                atom1_item = self.host.data.atoms[id1]["item"]
                atom2_item = self.host.data.atoms[id2]["item"]
                bond_item = BondItem(
                    atom1_item, atom2_item, data.get("order", 1), data.get("stereo", 0)
                )
                self.host.data.bonds[key_tuple] = {
                    "order": data.get("order", 1),
                    "stereo": data.get("stereo", 0),
                    "item": bond_item,
                }
                atom1_item.bonds.append(bond_item)
                atom2_item.bonds.append(bond_item)
                self.host.scene.addItem(bond_item)

        for atom_data in self.host.data.atoms.values():
            if atom_data["item"]:
                atom_data["item"].update_style()
        self.host.scene.update_all_items()
        mol_3d_data = loaded_data.get("mol_3d")
        if mol_3d_data is not None:
            try:
                self.host.current_mol = Chem.Mol(mol_3d_data)
                # Debug: check if 3D structure is valid
                if self.host.current_mol and self.host.current_mol.GetNumAtoms() > 0:
                    # Restore _original_atom_id if present in saved state
                    mol_3d_atom_ids = loaded_data.get("mol_3d_atom_ids")
                    if (
                        mol_3d_atom_ids
                        and len(mol_3d_atom_ids) == self.host.current_mol.GetNumAtoms()
                    ):
                        for i, aid in enumerate(mol_3d_atom_ids):
                            if aid is not None:
                                rd_atom = self.host.current_mol.GetAtomWithIdx(i)
                                if rd_atom:
                                    try:
                                        rd_atom.SetIntProp(
                                            "_original_atom_id", int(aid)
                                        )
                                    except (RuntimeError, ValueError, TypeError):
                                        pass

                    # Sync 2D atoms with 3D actors
                    if self.host.compute_manager:
                        try:
                            self.host.compute_manager.create_atom_id_mapping()
                            self.host.view_3d_manager.update_atom_id_menu_text()
                            self.host.view_3d_manager.update_atom_id_menu_state()
                        except Exception as e:
                            logging.debug(
                                f"Partial failure during ID mapping restoration: {e}"
                            )

                    # draw_molecule_3d will use restored IDs
                    self.host.view_3d_manager.draw_molecule_3d(self.host.current_mol)
                    if (
                        hasattr(self.host, "plotter")
                        and self.host.plotter
                        and hasattr(self.host.plotter, "reset_camera")
                    ):
                        self.host.plotter.reset_camera()

                    self.host.ui_manager._enable_3d_features(True)
                    self.host.view_3d_manager.setup_3d_hover()
                else:
                    self.host.current_mol = None
                    if hasattr(self.host, "plotter") and self.host.plotter:
                        self.host.plotter.clear()
                    self.host.ui_manager._enable_3d_features(False)
            except (RuntimeError, ValueError, TypeError) as e:
                logging.error(f"Could not load 3D model from state data: {e}")
                if hasattr(self.host, "statusBar") and self.host.statusBar():
                    self.host.statusBar().showMessage(f"Error loading 3D model: {e}", 5000)
                self.host.current_mol = None
                self.host.ui_manager._enable_3d_features(False)

        else:
            self.host.current_mol = None
            self.host.plotter.clear()
            self.host.analysis_action.setEnabled(False)
            self.host.optimize_3d_button.setEnabled(False)
            # Disable 3D features
            self.host.ui_manager._enable_3d_features(False)

        self.host.edit_actions_manager.update_implicit_hydrogens()
        self.host.view_3d_manager.update_chiral_labels()

        if loaded_data.get("is_3d_viewer_mode", False):
            self.host.ui_manager._enter_3d_viewer_ui_mode()
            self.host.statusBar().showMessage("Project loaded in 3D Viewer Mode.")
        else:
            self.host.ui_manager.restore_ui_for_editing()
            # Enable 3D edit features even in 2D editor mode if 3D molecule exists
            if self.host.current_mol and self.host.current_mol.GetNumAtoms() > 0:
                self.host.ui_manager._enable_3d_edit_actions(True)

        # Update labels after undo/redo
        self.host.edit_3d_manager.update_2d_measurement_labels()

    def push_undo_state(self):
        self.host.edit_actions_manager.push_undo_state()

    def update_window_title(self):
        """Update window title to reflect save state."""
        base_title = f"MoleditPy Ver. {VERSION}"
        if self.host.current_file_path:
            filename = os.path.basename(self.host.current_file_path)
            title = f"{filename} - {base_title}"
            if self.host.has_unsaved_changes:
                title = f"*{title}"
        else:
            # Handle as Untitled
            title = f"Untitled - {base_title}"
            if self.host.has_unsaved_changes:
                title = f"*{title}"
        self.host.setWindowTitle(title)

    def check_unsaved_changes(self):
        """Check for unsaved changes and show warning."""
        if not self.host.has_unsaved_changes:
            return True  # Saved or no changes

        if not self.host.data.atoms and self.host.current_mol is None:
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
            file_path = self.host.current_file_path
            if not file_path or not file_path.lower().endswith(".pmeprj"):
                self.save_project_as()
            else:
                self.save_project()
            return (
                not self.host.has_unsaved_changes
            )  # Return True only if save was successful
        elif reply == QMessageBox.StandardButton.No:
            return True  # Continue without saving
        else:
            return False  # Cancel

    def reset_undo_stack(self):
        self.host.edit_actions_manager.undo_stack.clear()
        self.host.edit_actions_manager.redo_stack.clear()
        self.push_undo_state()

    def undo(self):
        self.host.edit_actions_manager.undo()

    def redo(self):
        self.host.edit_actions_manager.redo()

    def update_undo_redo_actions(self):
        self.host.edit_actions_manager.update_undo_redo_actions()

    def update_realtime_info(self):
        """Show molecular info in status bar."""
        if not self.host.data.atoms:
            self.host.formula_label.setText("")  # Clear label if no atoms
            return

        if hasattr(self.host, "data") and self.host.data and hasattr(self.host.data, "to_rdkit_mol"):
            try:
                mol = self.host.data.to_rdkit_mol()
                if mol:
                    # Generate mol with explicit Hs
                    mol_with_hs = Chem.AddHs(mol)
                    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                    # Get atom count with Hs
                    num_atoms = mol_with_hs.GetNumAtoms()
                    # Update label text
                    if hasattr(self.host, "formula_label") and self.host.formula_label:
                        self.host.formula_label.setText(
                            f"Formula: {mol_formula}   |   Atoms: {num_atoms}"
                        )
            except (RuntimeError, TypeError, ValueError) as e:
                logging.debug(
                    f"Molecular info update suppressed for unstable structure: {e}"
                )
                if hasattr(self.host, "formula_label") and self.host.formula_label:
                    self.host.formula_label.setText("Invalid structure")

    def create_json_data(self):
        """Convert current state to PMEJSON."""
        # Metadata
        json_data = {
            "format": "PME Project",
            "version": "1.0",
            "application": "MoleditPy",
            "application_version": VERSION,
            "created": str(QDateTime.currentDateTime().toString(Qt.DateFormat.ISODate)),
            "is_3d_viewer_mode": not self.host.is_2d_editable,
        }

        # 2D data
        if self.host.data.atoms:
            atoms_2d = []
            for atom_id, data in self.host.data.atoms.items():
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
            for (atom1_id, atom2_id), bond_data in self.host.data.bonds.items():
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
                "next_atom_id": self.host.data._next_atom_id,
            }

        # 3D data
        if self.host.current_mol and self.host.current_mol.GetNumConformers() > 0:
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

                # Convert constraints to JSON compatible format
                json_safe_constraints = []
                try:
                    for const in self.host.constraints_3d:
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
                    "num_conformers": self.host.current_mol.GetNumConformers(),
                    "constraints_3d": json_safe_constraints,
                }

                # Molecular info
                json_data["molecular_info"] = {
                    "num_atoms": self.host.current_mol.GetNumAtoms(),
                    "num_bonds": self.host.current_mol.GetNumBonds(),
                    "molecular_weight": Descriptors.MolWt(self.host.current_mol),
                    "formula": rdMolDescriptors.CalcMolFormula(self.host.current_mol),
                }

                # Identifiers (SMILES/InChI)
                try:
                    json_data["identifiers"] = {
                        "smiles": Chem.MolToSmiles(self.host.current_mol),
                        "canonical_smiles": Chem.MolToSmiles(
                            self.host.current_mol, canonical=True
                        ),
                    }

                    # Attempt InChI generation
                    try:
                        inchi = Chem.MolToInchi(self.host.current_mol)
                        inchi_key = Chem.MolToInchiKey(self.host.current_mol)
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
        self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)
        self.host.ui_manager._enable_3d_edit_actions(False)
        self.host.ui_manager._enable_3d_features(False)

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
                raw_pos = (float(atom_data["x"]), float(atom_data["y"]))
                pos_q = QPointF(raw_pos[0], raw_pos[1])
                charge = atom_data.get("charge", 0)
                radical = atom_data.get("radical", 0)

                atom_item = AtomItem(
                    atom_id, symbol, pos_q, charge=charge, radical=radical
                )
                self.host.data.atoms[atom_id] = {
                    "symbol": symbol,
                    "pos": raw_pos,
                    "item": atom_item,
                    "charge": charge,
                    "radical": radical,
                }
                self.host.scene.addItem(atom_item)

            # Restore next_atom_id
            self.host.data._next_atom_id = structure_2d.get(
                "next_atom_id",
                max([atom["id"] for atom in atoms_2d]) + 1 if atoms_2d else 0,
            )

            # Restore bonds
            for bond_data in bonds_2d:
                atom1_id = bond_data["atom1"]
                atom2_id = bond_data["atom2"]

                if atom1_id in self.host.data.atoms and atom2_id in self.host.data.atoms:
                    atom1_item = self.host.data.atoms[atom1_id]["item"]
                    atom2_item = self.host.data.atoms[atom2_id]["item"]

                    bond_order = bond_data["order"]
                    stereo = bond_data.get("stereo", 0)

                    bond_item = BondItem(
                        atom1_item, atom2_item, bond_order, stereo=stereo
                    )
                    # Add to bond list (used for C visibility)
                    atom1_item.bonds.append(bond_item)
                    atom2_item.bonds.append(bond_item)

                    self.host.data.bonds[(atom1_id, atom2_id)] = {
                        "order": bond_order,
                        "item": bond_item,
                        "stereo": stereo,
                    }
                    self.host.scene.addItem(bond_item)

            # Update all AtomItem styles
            for atom in self.host.data.atoms.values():
                atom["item"].update_style()
            self.host.scene.update_all_items()
        # Restore 3D data
        structure_3d = json_data.get("3d_structure")
        if isinstance(structure_3d, dict):
            # Restore constraints
            loaded_constraints = structure_3d.get("constraints_3d", [])
            if loaded_constraints and isinstance(loaded_constraints, list):
                self.host.constraints_3d = []
                for const in loaded_constraints:
                    if isinstance(const, (list, tuple)):
                        try:
                            if len(const) == 4:
                                self.host.constraints_3d.append(
                                    (const[0], tuple(const[1]), const[2], const[3])
                                )
                            elif len(const) == 3:
                                self.host.constraints_3d.append(
                                    (const[0], tuple(const[1]), const[2], 1.0e5)
                                )
                        except (TypeError, ValueError, IndexError):
                            pass
            else:
                self.host.constraints_3d = []
            try:
                # Restore binary data
                mol_base64 = structure_3d.get("mol_binary_base64")
                if mol_base64:
                    mol_binary = base64.b64decode(mol_base64.encode("ascii"))
                    self.host.current_mol = Chem.Mol(mol_binary)
                    if self.host.current_mol:
                        # Set 3D coordinates
                        if self.host.current_mol.GetNumConformers() > 0:
                            conf = self.host.current_mol.GetConformer()
                            atoms_3d = structure_3d.get("atoms", [])
                            # Ensure numpy array size matches atoms in file
                            num_atoms_file = len(atoms_3d)
                            if num_atoms_file > 0:
                                self.host.view_3d_manager.atom_positions_3d = np.zeros((num_atoms_file, 3))
                                for atom_data in atoms_3d:
                                    idx = atom_data.get("index", -1)
                                    if 0 <= idx < num_atoms_file:
                                        self.host.view_3d_manager.atom_positions_3d[idx] = [
                                            atom_data.get("x", 0.0),
                                            atom_data.get("y", 0.0),
                                            atom_data.get("z", 0.0),
                                        ]
                                        # Restore original editor atom id into RDKit atom property
                                        original_id = atom_data.get("original_id")
                                        if original_id is not None:
                                            try:
                                                rd_atom = (
                                                    self.host.current_mol.GetAtomWithIdx(idx)
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
                            if hasattr(self.host.compute_manager, "create_atom_id_mapping"):
                                try:
                                    self.host.compute_manager.create_atom_id_mapping()
                                    if hasattr(self.host.view_3d_manager, "update_atom_id_menu_text"):
                                        self.host.view_3d_manager.update_atom_id_menu_text()
                                    if hasattr(self.host.view_3d_manager, "update_atom_id_menu_state"):
                                        self.host.view_3d_manager.update_atom_id_menu_state()
                                except (RuntimeError, TypeError, AttributeError):
                                    pass

                        # Always show 3D if 3D molecule exists
                        if hasattr(self.host.view_3d_manager, "draw_molecule_3d"):
                            self.host.view_3d_manager.draw_molecule_3d(self.host.current_mol)

                        # Switch UI in Viewer mode
                        if is_3d_mode and hasattr(self.host.ui_manager, "_enter_3d_viewer_ui_mode"):
                            self.host.ui_manager._enter_3d_viewer_ui_mode()
                        else:
                            self.host.is_2d_editable = True

                        if (
                            hasattr(self.host, "plotter")
                            and self.host.plotter
                            and hasattr(self.host.plotter, "reset_camera")
                        ):
                            self.host.plotter.reset_camera()

                        # Enable 3D-related UI
                        try:
                            self.host.ui_manager._enable_3d_edit_actions(True)
                            self.host.ui_manager._enable_3d_features(True)
                        except (RuntimeError, TypeError, AttributeError):
                            pass
            except (RuntimeError, ValueError, TypeError, base64.binascii.Error) as e:
                logging.error(f"Could not restore 3D molecular data: {e}")
                self.host.current_mol = None


StateManager._cls = StateManager


# Backward-compat aliases
MainWindowAppState = StateManager
