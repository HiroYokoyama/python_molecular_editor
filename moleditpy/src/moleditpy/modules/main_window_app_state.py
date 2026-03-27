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
                try:
                    if atom.HasProp("_original_atom_id"):
                        mol_3d_atom_ids.append(atom.GetIntProp("_original_atom_id"))
                    else:
                        mol_3d_atom_ids.append(None)
                except (AttributeError, RuntimeError, TypeError):
                    # Suppress minor property access errors if RDKit object state is inconsistent
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
            try:
                parts = []
                for p in v_str.split("."):
                    # Extract numeric part from start of string (e.g., '0a1' -> 0)
                    import re
                    m = re.match(r"(\d+)", p)
                    if m:
                        parts.append(int(m.group(1)))
                    else:
                        parts.append(0)
                return tuple(parts)
            except (ValueError, AttributeError):
                return (0, 0, 0)

        try:
            app_version_parts = parse_v(VERSION)
            file_version_parts = parse_v(file_version_str)

            # Warn if file version is newer than app version
            if file_version_parts > app_version_parts:
                QMessageBox.warning(
                    self,
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).\n\n"
                    f"Your current version is {VERSION}.\n\n"
                    "Some features may not load or work correctly.",
                )
        except (ValueError, AttributeError, IndexError):
            # Suppress non-critical version parsing errors if file metadata is malformed
            pass

        raw_atoms = loaded_data.get("atoms", {})
        raw_bonds = loaded_data.get("bonds", {})

        # Restore constraints
        try:
            loaded_constraints = loaded_data.get("constraints_3d", [])
            # Assumes [Type, [Idx...], Value, Force] format
            self.constraints_3d = []
            for const in loaded_constraints:
                if isinstance(const, list):
                    if len(const) == 4:
                        # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                        self.constraints_3d.append(
                            (const[0], tuple(const[1]), const[2], const[3])
                        )
                    elif len(const) == 3:
                        # Backward compatibility: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                        self.constraints_3d.append(
                            (const[0], tuple(const[1]), const[2], 1.0e5)
                        )
        except (AttributeError, RuntimeError, TypeError, ValueError, IndexError):
            # Reset constraints on load failure to prevent inconsistent state
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

        if "mol_3d" in loaded_data and loaded_data["mol_3d"] is not None:
            try:
                self.current_mol = Chem.Mol(loaded_data["mol_3d"])
                # Debug: check if 3D structure is valid
                if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                    # Restore _original_atom_id if present in saved state
                    if "mol_3d_atom_ids" in loaded_data:
                        atom_ids = loaded_data["mol_3d_atom_ids"]
                        if len(atom_ids) == self.current_mol.GetNumAtoms():
                            for i, aid in enumerate(atom_ids):
                                if aid is not None:
                                    if self.current_mol and i < self.current_mol.GetNumAtoms():
                                        try:
                                            self.current_mol.GetAtomWithIdx(i).SetIntProp(
                                                "_original_atom_id", int(aid)
                                            )
                                        except (AttributeError, RuntimeError, TypeError, ValueError):  
                                            # Skip property assignment on failure if RDKit object is partially invalid
                                            pass
                    # Sync 2D atoms with 3D actors
                    try:
                        self.create_atom_id_mapping()
                        self.update_atom_id_menu_text()
                        self.update_atom_id_menu_state()
                    except (AttributeError, RuntimeError, TypeError, ValueError):  
                        # Suppress UI sync noise during state restoration if widgets are being recreated
                        pass
                    # draw_molecule_3d will use restored IDs
                    self.draw_molecule_3d(self.current_mol)
                    self.plotter.reset_camera()
                    # Enable 3D features
                    self._enable_3d_features(True)

                    # Setup 3D hover
                    self.setup_3d_hover()
                else:
                    # Handle invalid 3D structure
                    self.current_mol = None
                    self.plotter.clear()
                    # Disable 3D features
                    self._enable_3d_features(False)
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                self.statusBar().showMessage(
                    f"Could not load 3D model from project: {e}"
                )
                self.current_mol = None
                # Disable 3D features
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
                    if a.HasProp("_original_atom_id")
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
            return not self.has_unsaved_changes  # Return True only if save was successful
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

        try:
            mol = self.data.to_rdkit_mol()
            if mol:
                # Generate mol with explicit Hs
                mol_with_hs = Chem.AddHs(mol)
                mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                # Get atom count with Hs
                num_atoms = mol_with_hs.GetNumAtoms()
                # Update label text
                self.formula_label.setText(
                    f"Formula: {mol_formula}   |   Atoms: {num_atoms}"
                )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Suppress molecular info update errors if structure is partially invalid during real-time updates
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

        if self.plugin_manager and self.plugin_manager.save_handlers:
            for name, callback in self.plugin_manager.save_handlers.items():
                try:
                    p_state = callback()
                    # Ensure serializable? Use primitive types ideally.
                    plugin_data[name] = p_state
                except (AttributeError, RuntimeError, TypeError, KeyError) as e:
                    print(f"Error saving state for plugin {name}: {e}")

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
        if "plugins" in json_data:
            plugin_data = json_data["plugins"]
            for name, p_state in plugin_data.items():
                if self.plugin_manager and name in self.plugin_manager.load_handlers:
                    try:
                        self.plugin_manager.load_handlers[name](p_state)
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        # Log but suppress individual plugin load errors to prevent app-wide crash
                        print(f"Error loading state for plugin {name}: {e}")
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
            try:
                loaded_constraints = structure_3d.get("constraints_3d", [])
                self.constraints_3d = []
                for const in loaded_constraints:
                    if isinstance(const, list):
                        if len(const) == 4:
                            # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                            self.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], const[3])
                            )
                        elif len(const) == 3:
                            # Backward compatibility: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                            self.constraints_3d.append(
                                (const[0], tuple(const[1]), const[2], 1.0e5)
                            )
            except (AttributeError, RuntimeError, TypeError, ValueError, IndexError):
                # Reset constraints if restoration fails to prevent inconsistent state
                self.constraints_3d = []  # Reset on load failure

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
                            self.atom_positions_3d = np.zeros((len(atoms_3d), 3))
                            for atom_data in atoms_3d:
                                idx = atom_data["index"]
                                if idx < len(self.atom_positions_3d):
                                    self.atom_positions_3d[idx] = [
                                        atom_data["x"],
                                        atom_data["y"],
                                        atom_data["z"],
                                    ]
                                # Restore original editor atom id into RDKit atom property
                                try:
                                    original_id = atom_data.get("original_id", None)
                                    if (
                                        original_id is not None
                                        and idx < self.current_mol.GetNumAtoms()
                                    ):
                                        rd_atom = self.current_mol.GetAtomWithIdx(idx)
                                        # set as int prop so other code expecting _original_atom_id works
                                        rd_atom.SetIntProp(
                                            "_original_atom_id", int(original_id)
                                        )
                                except (AttributeError, RuntimeError, TypeError, ValueError):  
                                    # Skip UI update if widgets or mapping is partially invalid during restoration
                                    pass
                            # Build mapping from original 2D atom IDs to RDKit indices so
                            # 3D picks can be synchronized back to 2D AtomItems.
                            try:
                                self.create_atom_id_mapping()
                                # update menu and UI states that depend on original IDs
                                try:
                                    self.update_atom_id_menu_text()
                                    self.update_atom_id_menu_state()
                                except (AttributeError, RuntimeError, TypeError):  
                                    # Suppress traceback
                                    pass
                            except (AttributeError, RuntimeError, TypeError):
                                # non-fatal if mapping creation fails
                                pass

                        # Always show 3D if 3D molecule exists
                        self.draw_molecule_3d(self.current_mol)
                        # Switch UI in Viewer mode
                        if is_3d_mode:
                            self._enter_3d_viewer_ui_mode()
                        else:
                            self.is_2d_editable = True
                        self.plotter.reset_camera()

                        # Enable 3D-related UI as 3D molecule was successfully restored
                        try:
                            self._enable_3d_edit_actions(True)
                            self._enable_3d_features(True)
                        except (AttributeError, RuntimeError, TypeError, ValueError):  
                            # Suppress 3D UI activation noise if features or widgets are unavailable
                            pass
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Warning: Could not restore 3D molecular data: {e}")
                self.current_mol = None

MainWindowAppState._cls = MainWindowAppState
