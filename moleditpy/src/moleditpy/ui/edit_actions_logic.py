#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy  EA Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations
import contextlib
import io
import copy
import logging
import math
import pickle
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np

try:
    from .mol_geometry import is_problematic_valence, identify_valence_problems
except ImportError:
    from moleditpy.core.mol_geometry import identify_valence_problems

# RDKit imports (explicit to satisfy flake8 and used features)
from PyQt6.QtCore import QByteArray, QMimeData, QPointF, Qt, QTimer
from PyQt6.QtGui import QCursor

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

class Rotate2DDialog(QDialog):
    def __init__(self, parent: Optional[QWidget] = None, initial_angle: float = 0) -> None:
        super().__init__(parent)
        self.setWindowTitle("Rotate 2D")
        self.setFixedWidth(300)

        layout = QVBoxLayout(self)

        # Angle input
        input_layout = QHBoxLayout()
        input_layout.addWidget(QLabel("Angle (degrees):"))
        self.angle_spin = QSpinBox()
        self.angle_spin.setRange(-360, 360)
        self.angle_spin.setValue(initial_angle)
        input_layout.addWidget(self.angle_spin)
        layout.addLayout(input_layout)

        # Slider
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(-180, 180)
        self.slider.setValue(initial_angle)
        self.slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.slider.setTickInterval(15)
        layout.addWidget(self.slider)

        # Sync slider and spinbox
        self.angle_spin.valueChanged.connect(self.slider.setValue)
        self.slider.valueChanged.connect(self.angle_spin.setValue)

        # Buttons
        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("Rotate")
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

    def get_angle(self) -> float:
        return self.angle_spin.value()

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
    from .constants import CLIPBOARD_MIME_TYPE
    from .molecular_data import MolecularData
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.atom_item import AtomItem
    from moleditpy.ui.bond_item import BondItem
    from moleditpy.utils.constants import CLIPBOARD_MIME_TYPE
    from moleditpy.core.molecular_data import MolecularData

try:
    # Import the shared SIP helper used across the package. This is
    # defined in modules/__init__.py and centralizes sip.isdeleted checks.
    from . import sip_isdeleted_safe
except ImportError:
    from moleditpy.utils import sip_isdeleted_safe

# --- Class Definition ---
class EditActionsManager:
    """Independent manager for molecular editing actions, ported from MainWindowEditActions mixin."""

    _cls = None

    def __init__(self, host: Any) -> None:
        self.host = host
        # State variables previously held by mixin
        self.last_rotation_angle: float = 0
        self.undo_stack: List[Dict[str, Any]] = []
        self.redo_stack: List[Dict[str, Any]] = []

    def push_undo_state(self) -> None:
        """Saves current molecular state to undo stack for history tracking."""
        if getattr(self.host, "_is_restoring_state", False):
            return

        current_state_for_comparison = {
            "atoms": {
                k: (
                    v["symbol"],
                    v["pos"].x() if hasattr(v["pos"], "x") else v["pos"][0],
                    v["pos"].y() if hasattr(v["pos"], "y") else v["pos"][1],
                    v.get("charge", 0),
                    v.get("radical", 0),
                )
                for k, v in self.host.data.atoms.items()
            },
            "bonds": {
                k: (v["order"], v.get("stereo", 0)) for k, v in self.host.data.bonds.items()
            },
            "_next_atom_id": self.host.data._next_atom_id,
            "mol_3d": self.host.current_mol.ToBinary() if self.host.current_mol else None,
            "mol_3d_atom_ids": [
                (
                    a.GetIntProp("_original_atom_id")
                    if a and a.HasProp("_original_atom_id")
                    else None
                )
                for a in self.host.current_mol.GetAtoms()
            ]
            if self.host.current_mol
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
                        v["pos"].x() if hasattr(v["pos"], "x") else v["pos"][0],
                        v["pos"].y() if hasattr(v["pos"], "y") else v["pos"][1],
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
            # Deepcopy state via pickling or explicit get_current_state call
            state = copy.deepcopy(self.host.state_manager.get_current_state())
            self.undo_stack.append(state)

            self.redo_stack.clear()
            # Record changes after initialization
            if getattr(self.host, "initialization_complete", False):
                self.host.has_unsaved_changes = True
                self.host.state_manager.update_window_title()

        self.update_implicit_hydrogens()
        if hasattr(self.host.state_manager, "update_realtime_info"):
            self.host.state_manager.update_realtime_info()
        self.update_undo_redo_actions()

    def undo(self) -> None:
        """Revert to the previous molecular state."""
        if len(self.undo_stack) > 1:
            self.redo_stack.append(self.undo_stack.pop())
            state = self.undo_stack[-1]
            self.host._is_restoring_state = True
            try:
                self.host.state_manager.set_state_from_data(state)
            finally:
                self.host._is_restoring_state = False

            # Re-evaluate menu states based on 3D structure after Undo
            if self.host.current_mol and self.host.current_mol.GetNumAtoms() > 0:
                self.host.ui_manager._enable_3d_edit_actions(True)
            else:
                self.host.ui_manager._enable_3d_edit_actions(False)

        self.update_undo_redo_actions()
        if hasattr(self.host.state_manager, "update_realtime_info"):
            self.host.state_manager.update_realtime_info()
        if hasattr(self.host, "view_2d") and self.host.view_2d:
            self.host.view_2d.setFocus()

    def redo(self) -> None:
        """Restore the previously undone molecular state."""
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self.host._is_restoring_state = True
            try:
                self.host.state_manager.set_state_from_data(state)
            finally:
                self.host._is_restoring_state = False

            # Re-evaluate menu states based on 3D structure after Redo
            if self.host.current_mol and self.host.current_mol.GetNumAtoms() > 0:
                self.host.ui_manager._enable_3d_edit_actions(True)
            else:
                self.host.ui_manager._enable_3d_edit_actions(False)

        self.update_undo_redo_actions()
        if hasattr(self.host.state_manager, "update_realtime_info"):
            self.host.state_manager.update_realtime_info()
        if hasattr(self.host, "view_2d") and self.host.view_2d:
            self.host.view_2d.setFocus()

    def update_undo_redo_actions(self) -> None:
        """Enable or disable Undo/Redo UI actions based on stack counts."""
        if hasattr(self.host, "undo_action"):
            self.host.undo_action.setEnabled(len(self.undo_stack) > 1)
        if hasattr(self.host, "redo_action"):
            self.host.redo_action.setEnabled(len(self.redo_stack) > 0)

    def copy_selection(self) -> None:
        """Copy selected atoms and bonds to clipboard"""
        try:
            selected_atoms = [
                item
                for item in self.host.scene.selectedItems()
                if isinstance(item, AtomItem)
            ]
            if not selected_atoms:
                return

            # Create set of selected atom IDs
            selected_atom_ids = {atom.atom_id for atom in selected_atoms}
            # Calculate geometric center of selected atoms
            center = QPointF(
                sum(atom.pos().x() for atom in selected_atoms) / len(selected_atoms),
                sum(atom.pos().y() for atom in selected_atoms) / len(selected_atoms),
            )

            # Store atom data with relative positions and map IDs to new indices
            atom_id_to_idx_map = {}
            fragment_atoms = []
            for i, atom in enumerate(selected_atoms):
                atom_id_to_idx_map[atom.atom_id] = i
                fragment_atoms.append(
                    {
                        "symbol": atom.symbol,
                        "rel_pos": atom.pos() - center,
                        "charge": atom.charge,
                        "radical": atom.radical,
                    }
                )

            # Store bonds between selected atoms
            fragment_bonds = []
            for (id1, id2), bond_data in self.host.data.bonds.items():
                if id1 in selected_atom_ids and id2 in selected_atom_ids:
                    fragment_bonds.append(
                        {
                            "idx1": atom_id_to_idx_map[id1],
                            "idx2": atom_id_to_idx_map[id2],
                            "order": bond_data["order"],
                            "stereo": bond_data.get(
                                "stereo", 0
                            ),  # Also save E/Z stereochemistry info
                        }
                    )

            # Serialize data to byte array using pickle
            data_to_pickle = {"atoms": fragment_atoms, "bonds": fragment_bonds}
            byte_array = QByteArray()
            buffer = io.BytesIO()
            pickle.dump(data_to_pickle, buffer)
            byte_array.append(buffer.getvalue())

            # Set clipboard with custom MIME type
            mime_data = QMimeData()
            mime_data.setData(CLIPBOARD_MIME_TYPE, byte_array)
            QApplication.clipboard().setMimeData(mime_data)
            self.host.statusBar().showMessage(
                f"Copied {len(fragment_atoms)} atoms and {len(fragment_bonds)} bonds."
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during copy operation: {e}")
            self.host.statusBar().showMessage(f"Error during copy operation: {e}")

    def cut_selection(self) -> None:
        """Cut selected items (copy then delete)"""
        try:
            selected_items = self.host.scene.selectedItems()
            if not selected_items:
                return

            # Execute copy process first
            self.copy_selection()

            if self.host.scene.delete_items(set(selected_items)):
                self.host.state_manager.push_undo_state()
                self.host.statusBar().showMessage("Cut selection.", 2000)

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during cut operation: {e}")
            self.host.statusBar().showMessage(f"Error during cut operation: {e}")

    def paste_from_clipboard(self) -> None:
        """Paste molecular fragment from clipboard"""
        try:
            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            if not mime_data.hasFormat(CLIPBOARD_MIME_TYPE):
                return

            byte_array = mime_data.data(CLIPBOARD_MIME_TYPE)
            buffer = io.BytesIO(byte_array)
            try:
                fragment_data = pickle.load(buffer)
            except pickle.UnpicklingError:
                self.host.statusBar().showMessage(
                    "Error: Invalid clipboard data format"
                )
                return

            paste_center_pos = self.host.view_2d.mapToScene(
                self.host.view_2d.mapFromGlobal(QCursor.pos())
            )
            self.host.scene.clearSelection()

            new_atoms = []
            for atom_data in fragment_data["atoms"]:
                pos = paste_center_pos + atom_data["rel_pos"]
                new_id = self.host.scene.create_atom(
                    atom_data["symbol"],
                    pos,
                    charge=atom_data.get("charge", 0),
                    radical=atom_data.get("radical", 0),
                )
                new_item = self.host.data.atoms[new_id]["item"]
                new_atoms.append(new_item)
                new_item.setSelected(True)

            for bond_data in fragment_data["bonds"]:
                atom1 = new_atoms[bond_data["idx1"]]
                atom2 = new_atoms[bond_data["idx2"]]
                self.host.scene.create_bond(
                    atom1,
                    atom2,
                    bond_order=bond_data.get("order", 1),
                    bond_stereo=bond_data.get(
                        "stereo", 0
                    ),  # Restore E/Z stereochemistry information as well
                )

            self.host.state_manager.push_undo_state()
            self.host.statusBar().showMessage(
                f"Pasted {len(fragment_data['atoms'])} atoms and {len(fragment_data['bonds'])} bonds.",
                2000,
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during paste operation: {e}")
            self.host.statusBar().showMessage(f"Error during paste operation: {e}")
        self.activate_select_mode()
        self.host.scene.update_all_items()

    def remove_hydrogen_atoms(self) -> None:
        """Delete hydrogen atoms and their bonds in 2D view"""
        try:
            # Collect hydrogen atom items robustly (store atom_id -> item)
            hydrogen_map = {}

            # Iterate over a snapshot of atoms to avoid "dictionary changed size"
            for atom_id, atom_data in list(self.host.data.atoms.items()):
                try:
                    if atom_data.get("symbol") != "H":
                        continue
                    item = atom_data.get("item")
                    # Only collect live AtomItem wrappers
                    if item is None:
                        continue
                    if sip_isdeleted_safe(item):
                        continue
                    if not isinstance(item, AtomItem):
                        continue
                    # Prefer storing by original atom id to detect actual removals later
                    hydrogen_map[atom_id] = item
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    # Ignore problematic entries and continue scanning
                    continue

            if not hydrogen_map:
                self.host.statusBar().showMessage(
                    "No hydrogen atoms found to remove.", 2000
                )
                return

            # To avoid blocking the UI or causing large, monolithic deletions that may
            # trigger internal re-entrancy issues, delete in batches and process UI events
            items = list(hydrogen_map.values())
            total = len(items)
            batch_size = 200  # tuned conservative batch size
            deleted_any = False

            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                batch = set()
                # Filter out items that are already deleted or invalid just before deletion
                for it in items[start:end]:
                    try:
                        if it is None:
                            continue
                        if sip_isdeleted_safe(it):
                            continue
                        if not isinstance(it, AtomItem):
                            continue
                        batch.add(it)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        continue

                if not batch:
                    # Nothing valid to delete in this batch
                    continue

                try:
                    # scene.delete_items is expected to handle bond cleanup; call it per-batch
                    success = False
                    try:
                        success = bool(self.host.scene.delete_items(batch))
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        # If scene.delete_items raises for a batch, attempt a safe per-item fallback
                        success = False

                    if not success:
                        # Fallback: try deleting items one-by-one to isolate problematic items
                        for it in list(batch):
                            try:
                                # Use scene.delete_items for single-item as well
                                ok = bool(self.host.scene.delete_items({it}))
                                if ok:
                                    deleted_any = True
                            except (
                                AttributeError,
                                RuntimeError,
                                ValueError,
                                TypeError,
                            ):
                                # If single deletion also fails, skip that item
                                continue
                    else:
                        deleted_any = True

                except (AttributeError, RuntimeError, ValueError, TypeError):
                    # Continue with next batch on unexpected errors
                    continue

                try:
                    QApplication.processEvents()
                except RuntimeError:
                    # Suppress non-critical error
                    pass
            # Determine how many hydrogens actually were removed by re-scanning data
            remaining_h = 0
            try:
                for _, atom_data in list(self.host.data.atoms.items()):
                    try:
                        if atom_data.get("symbol") == "H":
                            remaining_h += 1
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        continue
            except (AttributeError, RuntimeError, ValueError, TypeError):
                remaining_h = 0

            removed_count = max(0, len(hydrogen_map) - remaining_h)

            if removed_count > 0:
                # Only push a single undo state once for the whole operation
                self.host.state_manager.push_undo_state()
                self.host.statusBar().showMessage(
                    f"Removed {removed_count} hydrogen atoms.", 2000
                )
            else:
                # If nothing removed but we attempted, show an informative message
                if deleted_any:
                    # Deleted something but couldn't determine count reliably
                    self.host.statusBar().showMessage(
                        "Removed hydrogen atoms (count unknown).", 2000
                    )
                else:
                    self.host.statusBar().showMessage(
                        "Failed to remove hydrogen atoms or none found."
                    )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during hydrogen removal: {e}")
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                # Suppress transient errors during UI status reporting.
                self.host.statusBar().showMessage(f"Error removing hydrogen atoms: {e}")

    def add_hydrogen_atoms(self) -> None:
        """Compute and add explicit hydrogens in 2D view using RDKit."""

        try:
            mol = self.host.data.to_rdkit_mol(use_2d_stereo=False)
            if not mol or mol.GetNumAtoms() == 0:
                self.host.statusBar().showMessage(
                    "No molecule available to compute hydrogens.", 2000
                )
                return

            added_count = 0
            added_items = []

            # Check implicit hydrogens for all RDKit atoms
            for idx in range(mol.GetNumAtoms()):
                rd_atom = mol.GetAtomWithIdx(idx)
                try:
                    orig_id = rd_atom.GetIntProp("_original_atom_id")
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    # Skip if original editor ID is missing or entry is already stale.
                    continue

                if orig_id not in self.host.data.atoms:
                    continue

                # Get implicit hydrogens; fallback to total - explicit
                implicit_h = (
                    int(rd_atom.GetNumImplicitHs())
                    if hasattr(rd_atom, "GetNumImplicitHs")
                    else 0
                )
                if implicit_h is None or implicit_h < 0:
                    implicit_h = 0
                if implicit_h == 0:
                    # Fallback
                    try:
                        total_h = int(rd_atom.GetTotalNumHs())
                        explicit_h = (
                            int(rd_atom.GetNumExplicitHs())
                            if hasattr(rd_atom, "GetNumExplicitHs")
                            else 0
                        )
                        implicit_h = max(0, total_h - explicit_h)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        implicit_h = 0

                if implicit_h <= 0:
                    continue

                parent_item = self.host.data.atoms[orig_id]["item"]
                parent_pos = parent_item.pos()

                # Determine angles based on neighbors to avoid collisions
                neighbor_angles = []
                try:
                    for (a1, a2), bdata in self.host.data.bonds.items():
                        # Collect neighboring atom angles (ignore H)
                        try:
                            if a1 == orig_id and a2 in self.host.data.atoms:
                                neigh = self.host.data.atoms[a2]
                                if neigh.get("symbol") == "H":
                                    continue
                                if neigh.get("item") is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get("item")):
                                    continue
                                vec = neigh["item"].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                            elif a2 == orig_id and a1 in self.host.data.atoms:
                                neigh = self.host.data.atoms[a1]
                                if neigh.get("symbol") == "H":
                                    continue
                                if neigh.get("item") is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get("item")):
                                    continue
                                vec = neigh["item"].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            continue
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    neighbor_angles = []

                # Set bond length in pixels
                bond_length = 75

                # Helper: determine bond_stereo for hydrogen
                def _choose_stereo(i):
                    # 0: plain, 1: wedge, 2: dash, 3: plain, 4+: all plain
                    if i == 0:
                        return 0
                    if i == 1:
                        return 1
                    if i == 2:
                        return 2
                    return 0  # 4th+ hydrogens are all plain

                # Improve placement: find largest gap in existing angles
                target_angles = []
                try:
                    if not neighbor_angles:
                        # No existing bonds: space evenly around circle
                        for h_idx in range(implicit_h):
                            angle = (2.0 * math.pi * h_idx) / implicit_h
                            target_angles.append(angle)
                    else:
                        # Normalize and sort
                        angs = [
                            ((a + 2.0 * math.pi) if a < 0 else a)
                            for a in neighbor_angles
                        ]
                        angs = sorted(angs)
                        # Calculate gaps (including wrap-around)
                        gaps = []  # list of (gap_size, start_angle, end_angle)
                        for i in range(len(angs)):
                            a1 = angs[i]
                            a2 = angs[(i + 1) % len(angs)]
                            if i == len(angs) - 1:
                                # wrap-around gap
                                gap = (a2 + 2.0 * math.pi) - a1
                                start = a1
                                end = a2 + 2.0 * math.pi
                            else:
                                gap = a2 - a1
                                start = a1
                                end = a2
                            gaps.append((gap, start, end))

                        # Select largest gap and space hydrogens evenly
                        gaps.sort(key=lambda x: x[0], reverse=True)
                        max_gap, gstart, gend = gaps[0]
                        for i in range(implicit_h):
                            seg = max_gap / (implicit_h + 1)
                            angle = gstart + (i + 1) * seg
                            # Normalize back to 0..2pi
                            angle = angle % (2.0 * math.pi)
                            target_angles.append(angle)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    # Fallback: simple even spacing
                    for h_idx in range(implicit_h):
                        angle = (2.0 * math.pi * h_idx) / implicit_h
                        target_angles.append(angle)

                # Add atom/bond at calculated position
                for h_idx, angle in enumerate(target_angles):
                    dx = bond_length * math.cos(angle)
                    dy = bond_length * math.sin(angle)
                    pos = QPointF(parent_pos.x() + dx, parent_pos.y() + dy)

                    # Create new hydrogen atom
                    try:
                        new_id = self.host.scene.create_atom("H", pos)
                        new_item = self.host.data.atoms[new_id]["item"]
                        # Set bond_stereo (plain, wedge, dash)
                        stereo = _choose_stereo(h_idx)
                        self.host.scene.create_bond(
                            parent_item, new_item, bond_order=1, bond_stereo=stereo
                        )
                        added_items.append(new_item)
                        added_count += 1
                    except (AttributeError, RuntimeError, ValueError) as e:
                        print(f"Failed to add H for atom {orig_id}: {e}")

            if added_count > 0:
                self.host.scene.update_all_items()
                self.host.state_manager.push_undo_state()
                self.host.statusBar().showMessage(
                    f"Added {added_count} hydrogen atoms.", 2000
                )
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    # Suppress selection errors if the scene is being cleared or items are invalid.
                    self.host.scene.clearSelection()
                    for it in added_items:
                        it.setSelected(True)
            else:
                self.host.statusBar().showMessage(
                    "No implicit hydrogens found to add.", 2000
                )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during hydrogen addition: {e}")
            self.host.statusBar().showMessage(f"Error adding hydrogen atoms: {e}")

    def update_edit_menu_actions(self) -> None:
        """Update edit menu based on selection and clipboard"""
        try:
            has_selection = len(self.host.scene.selectedItems()) > 0
            self.host.cut_action.setEnabled(has_selection)
            self.host.copy_action.setEnabled(has_selection)

            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            self.host.paste_action.setEnabled(
                mime_data is not None and mime_data.hasFormat(CLIPBOARD_MIME_TYPE)
            )
        except RuntimeError:
            # Suppress non-critical error
            pass

    def open_rotate_2d_dialog(self) -> None:
        """Open 2D rotation dialog"""
        # Initialize last_rotation_angle if not present
        if not hasattr(self, "last_rotation_angle"):
            self.last_rotation_angle = 0

        dialog = Rotate2DDialog(self, initial_angle=self.last_rotation_angle)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            angle = dialog.get_angle()
            self.last_rotation_angle = angle  # Remember for next time
            self.rotate_molecule_2d(angle)

    def rotate_molecule_2d(self, angle_degrees: float) -> None:
        """Rotate 2D molecule (selection or entire)"""
        try:
            # Determine target atoms
            selected_items = self.host.scene.selectedItems()
            target_atoms = [
                item for item in selected_items if isinstance(item, AtomItem)
            ]

            # If no selection, rotate everything
            if not target_atoms:
                target_atoms = [
                    data["item"]
                    for data in self.host.data.atoms.values()
                    if data.get("item") and not sip_isdeleted_safe(data["item"])
                ]

            if not target_atoms:
                self.host.statusBar().showMessage("No atoms to rotate.")
                return

            # Calculate Center
            xs = [atom.pos().x() for atom in target_atoms]
            ys = [atom.pos().y() for atom in target_atoms]
            if not xs:
                return

            center_x = sum(xs) / len(xs)
            center_y = sum(ys) / len(ys)

            # Map for core rotation
            points_map = {
                atom.atom_id: (atom.pos().x(), atom.pos().y()) for atom in target_atoms
            }

            from moleditpy.core.mol_geometry import rotate_2d_points

            new_positions = rotate_2d_points(
                points_map, center_x, center_y, angle_degrees
            )

            for atom in target_atoms:
                if atom.atom_id in new_positions:
                    nx, ny = new_positions[atom.atom_id]
                    new_pos = QPointF(nx, ny)
                    atom.setPos(new_pos)
                    self.host.data.set_atom_pos(atom.atom_id, new_pos)

            # Update bonds
            self.host.scene.update_connected_bonds(target_atoms)
            self.host.scene.update_all_items()

            self.host.state_manager.push_undo_state()
            self.host.statusBar().showMessage(
                f"Rotated {len(target_atoms)} atoms by {angle_degrees} degrees."
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error rotating molecule: {e}")
            self.host.statusBar().showMessage(f"Error rotating: {e}")

    def select_all(self) -> None:
        for item in self.host.scene.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.setSelected(True)

    def clear_all(self) -> bool:
        # Check for unsaved changes
        if not self.check_unsaved_changes():
            # Cancel if requested
            return False

        self.host.ui_manager.restore_ui_for_editing()

        # Reset 3D mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        if self.host.edit_3d_manager.is_3d_edit_mode:
            self.host.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)

        # Clear 3D selection
        self.host.edit_3d_manager.clear_3d_selection()

        self.dragged_atom_info = None

        # Clear 2D editor (no undo push)
        self.clear_2d_editor(push_to_undo=False)

        # Clear 3D model
        self.host.current_mol = None
        self.host.plotter.clear()
        self.host.constraints_3d = []

        # Disable 3D features
        self.host.ui_manager._enable_3d_features(False)

        # Reset undo/redo stack
        self.host.state_manager.reset_undo_stack()

        # Reset file state
        self.host.has_unsaved_changes = False
        self.host.current_file_path = None
        self.host.state_manager.update_window_title()

        # Reset 2D zoom
        self.reset_zoom()

        # Update scene and view
        self.host.scene.update()
        if self.host.view_2d:
            self.host.view_2d.viewport().update()

        # Disable 3D features
        self.host.ui_manager._enable_3d_features(False)

        # Redraw 3D plotter
        self.host.plotter.render()

        # Update menu text and state
        self.host.view_3d_manager.update_atom_id_menu_text()
        self.host.view_3d_manager.update_atom_id_menu_state()

        # Force UI event processing
        QApplication.processEvents()

        # Call plugin document reset handlers
        if hasattr(self.host, "plugin_manager") and self.host.plugin_manager:
            self.host.plugin_manager.invoke_document_reset_handlers()

        self.host.statusBar().showMessage("Cleared all data.")
        return True

    def clear_2d_editor(self, push_to_undo: bool = True) -> None:
        # Clear 2D editor (no undo push)
        self.host.data = MolecularData()
        self.host.scene.data = self.host.data
        self.host.scene.clear()
        self.host.scene.reinitialize_items()
        self.host.is_xyz_derived = False
        # self.host.current_mol is now cleared via self.host.current_mol = None if needed,
        # but usually it's handled in clear_all.

        # Also clear measurement labels
        self.clear_2d_measurement_labels()

        # Clear 3D data and disable 3D-related menus
        self.host.current_mol = None
        self.host.plotter.clear()
        # Disable 3D features
        self.host.ui_manager._enable_3d_features(False)

        if push_to_undo:
            self.host.state_manager.push_undo_state()

    def _compute_h_counts(self, mol: Any) -> Dict[int, int]:
        """Build a mapping of original_id -> hydrogen count without touching Qt items."""
        h_count_map = {}
        if mol is None:
            # Invalid/unsanitizable structure: reset all counts to 0
            for atom_id in list(self.host.data.atoms.keys()):
                h_count_map[atom_id] = 0
            return h_count_map

        for atom in mol.GetAtoms():
            try:
                if not atom.HasProp("_original_atom_id"):
                    continue
                original_id = atom.GetIntProp("_original_atom_id")

                # Robust retrieval of H counts: prefer implicit, fallback to total or 0
                try:
                    h_count = int(atom.GetNumImplicitHs())
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    try:
                        h_count = int(atom.GetTotalNumHs())
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        h_count = 0

                h_count_map[int(original_id)] = h_count
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Skip problematic RDKit atoms
                continue
        return h_count_map

    def _detect_chemistry_problems(self, mol: Optional[Chem.Mol]) -> Dict[int, bool]:
        """Compute a per-atom problem map (original_id -> bool)."""
        problem_map = {}
        try:
            if mol is not None:
                try:
                    problems = Chem.DetectChemistryProblems(mol)
                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    logging.warning(f"RDKit DetectChemistryProblems failed: {e}")
                    problems = None

                if problems:
                    for prob in problems:
                        try:
                            atom_idx = prob.GetAtomIdx()
                            rd_atom = mol.GetAtomWithIdx(atom_idx)
                            if rd_atom and rd_atom.HasProp("_original_atom_id"):
                                orig = int(rd_atom.GetIntProp("_original_atom_id"))
                                problem_map[orig] = True
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            continue
            else:
                # Fallback: use the pure-logic valence heuristic from mol_geometry
                for atom_id in identify_valence_problems(
                    self.host.data.atoms, self.host.data.bonds
                ):
                    problem_map[atom_id] = True
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.error(f"Error during chemistry problem detection: {e}")

        return problem_map

    def _apply_ui_h_counts(self, h_count_map: Dict[int, int], problem_map: Dict[int, bool], my_token: int) -> None:
        """Apply the computed H counts and problem flags to UI items on the main thread."""
        # If the global counter changed since this closure was
        # created, bail out  Ethe update is stale.
        try:
            if my_token != getattr(self.host, "_ih_update_counter", None):
                return
        except (AttributeError, RuntimeError, ValueError, TypeError):
            return

        atoms_snapshot = (
            dict(self.host.data.atoms)
            if (hasattr(self.host, "data") and hasattr(self.host.data, "atoms"))
            else {}
        )
        is_deleted_func = sip_isdeleted_safe

        items_to_update = []
        for atom_id, atom_data in atoms_snapshot.items():
            try:
                item = atom_data.get("item")
                if not item:
                    continue

                # Suppress potential errors if the item is already destroyed by SIP during iteration
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    if is_deleted_func and is_deleted_func(item):
                        continue

                # Check if the item is no longer in a scene: skip updating it to avoid
                # touching partially-deleted objects during scene teardown.
                sc = item.scene() if hasattr(item, "scene") else None
                if sc is None:
                    continue

                # Desired new count (default to 0 if not computed)
                new_count = h_count_map.get(atom_id, 0)

                current = getattr(item, "implicit_h_count", None)
                current_prob = getattr(item, "has_problem", False)
                desired_prob = problem_map.get(atom_id, False)

                # If neither the implicit-H count nor the problem flag
                # changed, skip this item.
                if current == new_count and current_prob == desired_prob:
                    continue

                # Only prepare a geometry change if the implicit H count
                # changes (this may affect the item's bounding rect).
                need_geometry = current != new_count
                if need_geometry and hasattr(item, "prepareGeometryChange"):
                    item.prepareGeometryChange()
                item.implicit_h_count = new_count
                item.has_problem = bool(desired_prob)
                # Ensure the item is updated in the scene so paint() runs
                # when either geometry or problem-flag changed.
                items_to_update.append(item)

            except (AttributeError, RuntimeError, ValueError, TypeError):
                continue

        # Trigger updates once for unique items; wrap in try/except to avoid crashes
        seen = set()
        for it in items_to_update:
            try:
                if it is None:
                    continue
                oid = id(it)
                if oid in seen:
                    continue
                seen.add(oid)
                if hasattr(it, "update"):
                    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                        # Suppress transient errors during item update.
                        it.update()
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Ignore any unexpected errors when touching the item
                continue

    def update_implicit_hydrogens(self):
        """Update implicit hydrogen counts on AtomItems."""
        # Quick guards: nothing to do if no atoms or no QApplication
        if not self.host.data.atoms:
            return

        try:
            try:
                self.host._ih_update_counter += 1
            except (AttributeError, RuntimeError, ValueError, TypeError):
                self.host._ih_update_counter = getattr(self.host, "_ih_update_counter", 0) or 1
            my_token = self.host._ih_update_counter

            mol = None
            try:
                mol = self.host.data.to_rdkit_mol()
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                logging.debug(f"to_rdkit_mol failed during H-update: {e}")
                mol = None

            h_count_map = self._compute_h_counts(mol)
            problem_map = self._detect_chemistry_problems(mol)

            def _ui_closure():
                self._apply_ui_h_counts(h_count_map, problem_map, my_token)

            try:
                QTimer.singleShot(0, _ui_closure)
            except (RuntimeError, TypeError):
                # Fallback if QTimer fails (e.g. during app shutdown)
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    _ui_closure()
        except (AttributeError, RuntimeError, TypeError, ValueError) as e:
            logging.exception(f"Unexpected error in update_implicit_hydrogens: {e}")

    def clean_up_2d_structure(self):
        self.host.statusBar().showMessage("Optimizing 2D structure...")

        # Clear existing problem flags
        self.host.scene.clear_all_problem_flags()

        # Case: no atoms in 2D editor
        if not self.host.data.atoms:
            self.host.statusBar().showMessage("Error: No atoms to optimize.")
            return

        mol = self.host.data.to_rdkit_mol()
        if mol is None or mol.GetNumAtoms() == 0:
            # If RDKit conversion fails, check for chemistry problems
            if hasattr(self.host, "compute_manager") and hasattr(self.host.compute_manager, "check_chemistry_problems_fallback"):
                self.host.compute_manager.check_chemistry_problems_fallback()
            return

        try:
            from moleditpy.core.mol_geometry import optimize_2d_coords

            new_positions = optimize_2d_coords(mol)

            if not new_positions:
                self.host.statusBar().showMessage(
                    "Optimization failed to generate coordinates."
                )
                return

            # Centering logic: identify 2D view center and RDKit molecule centroid
            view_center = self.host.view_2d.mapToScene(
                self.host.view_2d.viewport().rect().center()
            )

            coords = list(new_positions.values())
            rdkit_cx = sum(p[0] for p in coords) / len(coords)
            rdkit_cy = sum(p[1] for p in coords) / len(coords)

            SCALE = 50.0  # 1.0 A = 50 pixels (matches 0.02 A/pixel constant)

            # Apply new coordinates with centering and scaling
            for atom_id, (nx, ny) in new_positions.items():
                if atom_id in self.host.data.atoms:
                    # Centered scaling: (coord - rdkit_center) * scale + scene_view_center
                    sx = ((nx - rdkit_cx) * SCALE) + view_center.x()
                    sy = (-(ny - rdkit_cy) * SCALE) + view_center.y()
                    new_pos = QPointF(sx, sy)

                    item = self.host.data.atoms[atom_id]["item"]
                    if item:
                        item.setPos(new_pos)
                    # Cache back to model
                    self.host.data.set_atom_pos(atom_id, new_pos)

            # Update all bond positions
            for bond_data in self.host.data.bonds.values():
                item = bond_data.get("item") if bond_data else None
                if not item:
                    continue
                if sip_isdeleted_safe(item):
                    continue
                if hasattr(item, "scene") and item.scene():
                    # Suppress potential errors if the item is already destroyed during coordinate adjustment
                    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                        item.update_position()

            # Run overlap resolution
            self.resolve_overlapping_groups()

            # Update measurement labels
            if hasattr(self.host.edit_3d_manager, "update_2d_measurement_labels"):
                self.host.edit_3d_manager.update_2d_measurement_labels()

            # Request scene update and ring re-analysis
            self.host.scene.update_all_items()

            self.host.statusBar().showMessage("2D structure optimization successful.")
            self.host.state_manager.push_undo_state()

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(f"Error during 2D optimization: {e}")
        finally:
            self.host.view_2d.setFocus()

    def redraw_molecule_3d(self):
        """Manually trigger redraw of the 3D molecule."""
        if hasattr(self.host, "current_mol") and self.host.current_mol:
            self.host.view_3d_manager.draw_molecule_3d(self.host.current_mol)
            self.host.statusBar().showMessage("Redraw complete.", 2000)
        else:
            self.host.statusBar().showMessage("No 3D molecule to redraw.")

    def resolve_overlapping_groups(self):
        """Detect and resolve overlapping atom groups."""

        # --- Parameters ---
        # Distance threshold for overlap
        OVERLAP_THRESHOLD = 0.5
        # Translation distance (bottom-left)
        MOVE_DISTANCE = 20

        # Safely retrieve item from self.host.data.atoms.values()
        all_atom_items = [
            data["item"]
            for data in self.host.data.atoms.values()
            if data and "item" in data
        ]

        if len(all_atom_items) < 2:
            return  # Step 1-3: Handled by core logic
        positions_map = {aid: data["pos"] for aid, data in self.host.data.atoms.items()}

        from moleditpy.core.mol_geometry import resolve_2d_overlaps

        def has_bond_check(id1, id2):
            item1 = self.host.data.atoms[id1]["item"]
            item2 = self.host.data.atoms[id2]["item"]
            return self.host.scene.find_bond_between(item1, item2) is not None

        move_operations = resolve_2d_overlaps(
            set(self.host.data.atoms.keys()),
            positions_map,
            self.host.data.adjacency_list,
            overlap_threshold=OVERLAP_THRESHOLD,
            move_distance=MOVE_DISTANCE,
            has_bond_check_func=has_bond_check,
        )

        if not move_operations:
            self.host.statusBar().showMessage("No overlapping atoms found.", 2000)
            return

        # Step 4: Execute translations
        for group_ids, (vx, vy) in move_operations:
            vector = QPointF(vx, vy)
            for atom_id in group_ids:
                item = self.host.data.atoms[atom_id]["item"]
                new_pos = item.pos() + vector
                item.setPos(new_pos)
                self.host.data.set_atom_pos(atom_id, new_pos)

        # Step 5: Update display and state
        for bond_data in self.host.data.bonds.values():
            item = bond_data.get("item") if bond_data else None
            if not item:
                continue
            try:
                if sip_isdeleted_safe(item):
                    continue
                if hasattr(item, "scene") and item.scene():
                    item.update_position()
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(f"Bond position update suppressed: {e}")

        # Update labels after resolution
        if hasattr(self.host.edit_3d_manager, "update_2d_measurement_labels"):
            self.host.edit_3d_manager.update_2d_measurement_labels()

        self.host.scene.update()
        self.host.state_manager.push_undo_state()
        self.host.statusBar().showMessage("Resolved overlapping groups.", 2000)

    def adjust_molecule_positions_to_avoid_collisions(self, mol, frags):
        """Adjust molecule positions to avoid collisions (BBox optimized)."""
        if len(frags) <= 1:
            return

        conf = mol.GetConformer()
        pt = Chem.GetPeriodicTable()

        # 1. Precompute fragment info (indices, VDW radii)
        frag_info = []
        for frag_indices in frags:
            positions = []
            vdw_radii = []
            for idx in frag_indices:
                pos = conf.GetAtomPosition(idx)
                positions.append(np.array([pos.x, pos.y, pos.z]))

                atom = mol.GetAtomWithIdx(idx)
                # GetRvdw() returns the Van der Waals radius
                try:
                    vdw_radii.append(pt.GetRvdw(atom.GetAtomicNum()))
                except RuntimeError:
                    vdw_radii.append(1.5)

            positions_np = np.array(positions)
            vdw_radii_np = np.array(vdw_radii)

            # Max VDW for box margin
            max_vdw = np.max(vdw_radii_np) if len(vdw_radii_np) > 0 else 0.0

            frag_info.append(
                {
                    "indices": frag_indices,
                    "centroid": np.mean(positions_np, axis=0),
                    "positions_np": positions_np,  # Keep as Numpy array
                    "vdw_radii_np": vdw_radii_np,  # Keep as Numpy array
                    "max_vdw_radius": max_vdw,
                    "bbox_min": np.zeros(3),  # To be calculated later
                    "bbox_max": np.zeros(3),  # To be calculated later
                }
            )

        # 2. Collision parameters
        collision_scale = 1.2  # 120% VDW
        max_iterations = 100
        moved = True
        iteration = 0

        while moved and iteration < max_iterations:
            moved = False
            iteration += 1

            # 3. Update BBoxes per iteration
            for i in range(len(frag_info)):
                # Recalculate box
                current_positions = []
                for idx in frag_info[i]["indices"]:
                    pos = conf.GetAtomPosition(idx)
                    current_positions.append([pos.x, pos.y, pos.z])

                positions_np = np.array(current_positions)
                frag_info[i]["positions_np"] = positions_np  # Update position info

                # Margin based on VDW and scale
                margin = frag_info[i]["max_vdw_radius"] * collision_scale

                frag_info[i]["bbox_min"] = np.min(positions_np, axis=0) - margin
                frag_info[i]["bbox_max"] = np.max(positions_np, axis=0) + margin

            # 4. Collision detection loop
            for i in range(len(frag_info)):
                for j in range(i + 1, len(frag_info)):
                    frag_i = frag_info[i]
                    frag_j = frag_info[j]

                    # AABB intersection check
                    # Check all axes
                    overlap_x = (
                        frag_i["bbox_min"][0] <= frag_j["bbox_max"][0]
                        and frag_i["bbox_max"][0] >= frag_j["bbox_min"][0]
                    )
                    overlap_y = (
                        frag_i["bbox_min"][1] <= frag_j["bbox_max"][1]
                        and frag_i["bbox_max"][1] >= frag_j["bbox_min"][1]
                    )
                    overlap_z = (
                        frag_i["bbox_min"][2] <= frag_j["bbox_max"][2]
                        and frag_i["bbox_max"][2] >= frag_j["bbox_min"][2]
                    )

                    if not (overlap_x and overlap_y and overlap_z):
                        # Skip if BBoxes don't overlap
                        continue

                    # Heavy per-atom check if BBoxes overlap
                    total_push_vector = np.zeros(3)
                    collision_count = 0

                    # Use precomputed Numpy arrays
                    positions_i = frag_i["positions_np"]
                    positions_j = frag_j["positions_np"]
                    vdw_i_all = frag_i["vdw_radii_np"]
                    vdw_j_all = frag_j["vdw_radii_np"]

                    for k, idx_i in enumerate(frag_i["indices"]):
                        pos_i = positions_i[k]
                        vdw_i = vdw_i_all[k]

                        for l, idx_j in enumerate(frag_j["indices"]):
                            pos_j = positions_j[l]
                            vdw_j = vdw_j_all[l]

                            distance_vec = pos_i - pos_j
                            distance_sq = np.dot(
                                distance_vec, distance_vec
                            )  # Avoid square root for speed

                            min_distance = (vdw_i + vdw_j) * collision_scale
                            min_distance_sq = min_distance * min_distance

                            if distance_sq < min_distance_sq and distance_sq > 0.0001:
                                distance = np.sqrt(distance_sq)
                                push_direction = distance_vec / distance
                                push_magnitude = (
                                    min_distance - distance
                                ) / 2  # Split push
                                total_push_vector += push_direction * push_magnitude
                                collision_count += 1

                    if collision_count > 0:
                        # Apply average push vector
                        avg_push_vector = total_push_vector / collision_count

                        # Update conformer coordinates
                        for idx in frag_i["indices"]:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos + avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())

                        for idx in frag_j["indices"]:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos - avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())

                        moved = True

    def _apply_chem_check_and_set_flags(self, mol, source_desc=None, force_skip=False):
        """Central helper to apply chemical sanitization (or skip it) and set
        chem_check_tried / chem_check_failed flags consistently.

        When sanitization fails, a warning is shown and the Optimize 3D button
        is disabled. If the user setting 'skip_chemistry_checks' is True, no
        sanitization is attempted and both flags remain False.
        """
        self.host.chem_check_tried = False
        self.host.chem_check_failed = False

        if force_skip or self.host.settings.get("skip_chemistry_checks", False):
            # User asked to skip chemistry checks entirely
            return

        try:
            Chem.SanitizeMol(mol)
            self.host.chem_check_tried = True
            self.host.chem_check_failed = False
        except (AttributeError, RuntimeError, ValueError, TypeError):
            self.host.chem_check_tried = True
            self.host.chem_check_failed = True
            desc = f" ({source_desc})" if source_desc else ""
            # Suppress potential status bar or button state errors if the window is being closed or destroyed
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.host.statusBar().showMessage(
                    f"Molecule sanitization failed{desc}; file may be malformed."
                )
            # Disable 3D optimization UI to prevent running on invalid molecules
            if hasattr(self.host, "optimize_3d_button"):
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    self.host.optimize_3d_button.setEnabled(False)

    def _clear_xyz_flags(self, mol=None):
        """Clear XYZ-derived markers from a molecule (or current_mol) and
        reset UI flags accordingly.

        This is a best-effort cleanup to remove properties like
        _xyz_skip_checks and _xyz_atom_data that may have been attached when
        an XYZ file was previously loaded. After clearing molecule-level
        markers, the UI flag self.host.is_xyz_derived is set to False and the
        Optimize 3D button is re-evaluated (enabled unless chem_check_failed
        is True).
        """
        target = mol if mol is not None else getattr(self.host, "current_mol", None)
        if target is not None:
            # Remove RDKit property _xyz_skip_checks
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                # Suppress error if HasProp or ClearProp is unavailable.
                if hasattr(target, "HasProp") and target.HasProp("_xyz_skip_checks"):
                    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                        target.ClearProp("_xyz_skip_checks")
            # Remove attribute-style markers
            target.__dict__.pop("_xyz_skip_checks", None)
            target.__dict__.pop("_xyz_atom_data", None)

        # Reset UI flag
        self.host.is_xyz_derived = False

        # Enable Optimize 3D unless sanitization failed
        if hasattr(self.host, "optimize_3d_button"):
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                # Suppress error if optimize_3d_button is partially destroyed.
                self.host.optimize_3d_button.setEnabled(
                    not getattr(self.host, "chem_check_failed", False)
                )

EditActionsManager._cls = EditActionsManager
