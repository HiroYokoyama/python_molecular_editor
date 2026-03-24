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
# Functional class: MainWindowEditActions
"""


import io
import itertools
import math
import pickle
import traceback
from collections import deque

import numpy as np

try:
    from .mol_geometry import is_problematic_valence
except ImportError:
    from modules.mol_geometry import is_problematic_valence

# RDKit imports (explicit to satisfy flake8 and used features)
from PyQt6.QtCore import QByteArray, QLineF, QMimeData, QPointF, Qt, QTimer
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
)
from rdkit import Chem
from rdkit.Chem import AllChem


class Rotate2DDialog(QDialog):
    def __init__(self, parent=None, initial_angle=0):
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

    def get_angle(self):
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
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem
    from modules.constants import CLIPBOARD_MIME_TYPE
    from modules.molecular_data import MolecularData


try:
    # Import the shared SIP helper used across the package. This is
    # defined in modules/__init__.py and centralizes sip.isdeleted checks.
    from . import sip_isdeleted_safe
except ImportError:
    from modules import sip_isdeleted_safe


# --- Class Definition ---
class MainWindowEditActions(object):
    """Functional class separated from main_window.py"""

    def copy_selection(self):
        """Copy selected atoms and bonds to clipboard"""
        try:
            selected_atoms = [
                item
                for item in self.scene.selectedItems()
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
            for (id1, id2), bond_data in self.data.bonds.items():
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
            self.statusBar().showMessage(
                f"Copied {len(fragment_atoms)} atoms and {len(fragment_bonds)} bonds."
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during copy operation: {e}")
            self.statusBar().showMessage(f"Error during copy operation: {e}")

    def cut_selection(self):
        """Cut selected items (copy then delete)"""
        try:
            selected_items = self.scene.selectedItems()
            if not selected_items:
                return

            # Execute copy process first
            self.copy_selection()

            if self.scene.delete_items(set(selected_items)):
                self.push_undo_state()
                self.statusBar().showMessage("Cut selection.", 2000)

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during cut operation: {e}")
            self.statusBar().showMessage(f"Error during cut operation: {e}")

    def paste_from_clipboard(self):
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
                self.statusBar().showMessage("Error: Invalid clipboard data format")
                return

            paste_center_pos = self.view_2d.mapToScene(
                self.view_2d.mapFromGlobal(QCursor.pos())
            )
            self.scene.clearSelection()

            new_atoms = []
            for atom_data in fragment_data["atoms"]:
                pos = paste_center_pos + atom_data["rel_pos"]
                new_id = self.scene.create_atom(
                    atom_data["symbol"],
                    pos,
                    charge=atom_data.get("charge", 0),
                    radical=atom_data.get("radical", 0),
                )
                new_item = self.data.atoms[new_id]["item"]
                new_atoms.append(new_item)
                new_item.setSelected(True)

            for bond_data in fragment_data["bonds"]:
                atom1 = new_atoms[bond_data["idx1"]]
                atom2 = new_atoms[bond_data["idx2"]]
                self.scene.create_bond(
                    atom1,
                    atom2,
                    bond_order=bond_data.get("order", 1),
                    bond_stereo=bond_data.get("stereo", 0),  # Restore E/Z stereochemistry information as well
                )

            self.push_undo_state()
            self.statusBar().showMessage(
                f"Pasted {len(fragment_data['atoms'])} atoms and {len(fragment_data['bonds'])} bonds.",
                2000,
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during paste operation: {e}")
            self.statusBar().showMessage(f"Error during paste operation: {e}")
        self.statusBar().showMessage(f"Pasted {len(new_atoms)} atoms.", 2000)
        self.activate_select_mode()

    def remove_hydrogen_atoms(self):
        """Delete hydrogen atoms and their bonds in 2D view"""
        try:
            # Collect hydrogen atom items robustly (store atom_id -> item)
            hydrogen_map = {}

            # Iterate over a snapshot of atoms to avoid "dictionary changed size"
            for atom_id, atom_data in list(self.data.atoms.items()):
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
                except (AttributeError, RuntimeError):
                    # Ignore problematic entries and continue scanning
                    continue

            if not hydrogen_map:
                self.statusBar().showMessage("No hydrogen atoms found to remove.", 2000)
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
                    except (AttributeError, RuntimeError):
                        continue

                if not batch:
                    # Nothing valid to delete in this batch
                    continue

                try:
                    # scene.delete_items is expected to handle bond cleanup; call it per-batch
                    success = False
                    try:
                        success = bool(self.scene.delete_items(batch))
                    except (AttributeError, RuntimeError):
                        # If scene.delete_items raises for a batch, attempt a safe per-item fallback
                        success = False

                    if not success:
                        # Fallback: try deleting items one-by-one to isolate problematic items
                        for it in list(batch):
                            try:
                                # Use scene.delete_items for single-item as well
                                ok = bool(self.scene.delete_items({it}))
                                if ok:
                                    deleted_any = True
                            except (AttributeError, RuntimeError):
                                # If single deletion also fails, skip that item
                                continue
                    else:
                        deleted_any = True

                except (AttributeError, RuntimeError):
                    # Continue with next batch on unexpected errors
                    continue

                # Allow the GUI to process events between batches to remain responsive
                try:
                    QApplication.processEvents()
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
            # Determine how many hydrogens actually were removed by re-scanning data
            remaining_h = 0
            try:
                for _, atom_data in list(self.data.atoms.items()):
                    try:
                        if atom_data.get("symbol") == "H":
                            remaining_h += 1
                    except (AttributeError, RuntimeError):
                        continue
            except (AttributeError, RuntimeError):
                remaining_h = 0

            removed_count = max(0, len(hydrogen_map) - remaining_h)

            if removed_count > 0:
                # Only push a single undo state once for the whole operation
                try:
                    self.push_undo_state()
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
                self.statusBar().showMessage(
                    f"Removed {removed_count} hydrogen atoms.", 2000
                )
            else:
                # If nothing removed but we attempted, show an informative message
                if deleted_any:
                    # Deleted something but couldn't determine count reliably
                    self.statusBar().showMessage(
                        "Removed hydrogen atoms (count unknown).", 2000
                    )
                else:
                    self.statusBar().showMessage(
                        "Failed to remove hydrogen atoms or none found."
                    )

        except (AttributeError, RuntimeError, ValueError) as e:
            # Capture and log unexpected errors but don't let them crash the UI
            print(f"Error during hydrogen removal: {e}")
            try:
                self.statusBar().showMessage(f"Error removing hydrogen atoms: {e}")
            except (AttributeError, RuntimeError):  # pragma: no cover
                pass

    def add_hydrogen_atoms(self):
        """Compute and add explicit hydrogens in 2D view using RDKit."""
        try:
            mol = self.data.to_rdkit_mol(use_2d_stereo=False)
            if not mol or mol.GetNumAtoms() == 0:
                self.statusBar().showMessage(
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
                except (AttributeError, RuntimeError):
                    # Skip if no original editor ID
                    continue

                if orig_id not in self.data.atoms:
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
                    except (AttributeError, RuntimeError):
                        implicit_h = 0

                if implicit_h <= 0:
                    continue

                parent_item = self.data.atoms[orig_id]["item"]
                parent_pos = parent_item.pos()

                # Determine angles based on neighbors to avoid collisions
                neighbor_angles = []
                try:
                    for (a1, a2), bdata in self.data.bonds.items():
                        # Collect neighboring atom angles (ignore H)
                        try:
                            if a1 == orig_id and a2 in self.data.atoms:
                                neigh = self.data.atoms[a2]
                                if neigh.get("symbol") == "H":
                                    continue
                                if neigh.get("item") is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get("item")):
                                    continue
                                vec = neigh["item"].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                            elif a2 == orig_id and a1 in self.data.atoms:
                                neigh = self.data.atoms[a1]
                                if neigh.get("symbol") == "H":
                                    continue
                                if neigh.get("item") is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get("item")):
                                    continue
                                vec = neigh["item"].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                        except (AttributeError, RuntimeError):
                            continue
                except (AttributeError, RuntimeError):
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
                except (AttributeError, RuntimeError):
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
                        new_id = self.scene.create_atom("H", pos)
                        new_item = self.data.atoms[new_id]["item"]
                        # Set bond_stereo (plain, wedge, dash)
                        stereo = _choose_stereo(h_idx)
                        self.scene.create_bond(
                            parent_item, new_item, bond_order=1, bond_stereo=stereo
                        )
                        added_items.append(new_item)
                        added_count += 1
                    except (AttributeError, RuntimeError, ValueError) as e:
                        print(f"Failed to add H for atom {orig_id}: {e}")

            if added_count > 0:
                self.push_undo_state()
                self.statusBar().showMessage(
                    f"Added {added_count} hydrogen atoms.", 2000
                )
                # Select newly added atoms
                try:
                    self.scene.clearSelection()
                    for it in added_items:
                        it.setSelected(True)
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
            else:
                self.statusBar().showMessage(
                    "No implicit hydrogens found to add.", 2000
                )

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during hydrogen addition: {e}")
            self.statusBar().showMessage(f"Error adding hydrogen atoms: {e}")

    def update_edit_menu_actions(self):
        """Update edit menu based on selection and clipboard"""
        try:
            has_selection = len(self.scene.selectedItems()) > 0
            self.cut_action.setEnabled(has_selection)
            self.copy_action.setEnabled(has_selection)

            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            self.paste_action.setEnabled(
                mime_data is not None and mime_data.hasFormat(CLIPBOARD_MIME_TYPE)
            )
        except RuntimeError:
            pass

    def open_rotate_2d_dialog(self):
        """Open 2D rotation dialog"""
        # Initialize last_rotation_angle if not present
        if not hasattr(self, "last_rotation_angle"):
            self.last_rotation_angle = 0

        dialog = Rotate2DDialog(self, initial_angle=self.last_rotation_angle)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            angle = dialog.get_angle()
            self.last_rotation_angle = angle  # Remember for next time
            self.rotate_molecule_2d(angle)

    def rotate_molecule_2d(self, angle_degrees):
        """Rotate 2D molecule (selection or entire)"""
        try:
            # Determine target atoms
            selected_items = self.scene.selectedItems()
            target_atoms = [
                item for item in selected_items if isinstance(item, AtomItem)
            ]

            # If no selection, rotate everything
            if not target_atoms:
                target_atoms = [
                    data["item"]
                    for data in self.data.atoms.values()
                    if data.get("item") and not sip_isdeleted_safe(data["item"])
                ]

            if not target_atoms:
                self.statusBar().showMessage("No atoms to rotate.")
                return

            # Calculate Center
            xs = [atom.pos().x() for atom in target_atoms]
            ys = [atom.pos().y() for atom in target_atoms]
            if not xs:
                return

            center_x = sum(xs) / len(xs)
            center_y = sum(ys) / len(ys)
            center = QPointF(center_x, center_y)

            rad = math.radians(angle_degrees)
            cos_a = math.cos(rad)
            sin_a = math.sin(rad)

            for atom in target_atoms:
                # Relative pos
                dx = atom.pos().x() - center_x
                dy = atom.pos().y() - center_y

                # Rotate
                new_dx = dx * cos_a - dy * sin_a
                new_dy = dx * sin_a + dy * cos_a

                new_pos = QPointF(center_x + new_dx, center_y + new_dy)
                atom.setPos(new_pos)

            # Update bonds
            self.scene.update_connected_bonds(target_atoms)

            self.push_undo_state()
            self.statusBar().showMessage(
                f"Rotated {len(target_atoms)} atoms by {angle_degrees} degrees."
            )
            self.scene.update()
            self.scene.update_all_items()

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error rotating molecule: {e}")
            self.statusBar().showMessage(f"Error rotating: {e}")

    def select_all(self):
        for item in self.scene.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.setSelected(True)

    def clear_all(self):
        # Check for unsaved changes
        if not self.check_unsaved_changes():
            # Cancel if requested
            return False

        self.restore_ui_for_editing()

        # Reset 3D mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)

        # Clear 3D selection
        self.clear_3d_selection()

        self.dragged_atom_info = None

        # Clear 2D editor (no undo push)
        self.clear_2d_editor(push_to_undo=False)

        # Clear 3D model
        self.current_mol = None
        self.plotter.clear()
        self.constraints_3d = []

        # Disable 3D features
        self._enable_3d_features(False)

        # Reset undo/redo stack
        self.reset_undo_stack()

        # Reset file state
        self.has_unsaved_changes = False
        self.current_file_path = None
        self.update_window_title()

        # Reset 2D zoom
        self.reset_zoom()

        # Update scene and view
        self.scene.update()
        if self.view_2d:
            self.view_2d.viewport().update()

        # Disable 3D features
        self._enable_3d_features(False)

        # Redraw 3D plotter
        self.plotter.render()

        # Update menu text and state
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

        # Force UI event processing
        QApplication.processEvents()

        # Call plugin document reset handlers
        if hasattr(self, "plugin_manager") and self.plugin_manager:
            self.plugin_manager.invoke_document_reset_handlers()

        self.statusBar().showMessage("Cleared all data.")
        return True

    def clear_2d_editor(self, push_to_undo=True):
        # Clear 2D editor (no undo push)
        self.data = MolecularData()
        self.scene.data = self.data
        self.scene.clear()
        self.scene.reinitialize_items()
        self.is_xyz_derived = False

        # Also clear measurement labels
        self.clear_2d_measurement_labels()

        # Clear 3D data and disable 3D-related menus
        self.current_mol = None
        self.plotter.clear()
        # Disable 3D features
        self._enable_3d_features(False)

        if push_to_undo:
            self.push_undo_state()

    def update_implicit_hydrogens(self):
        """Update implicit hydrogen counts on AtomItems."""
        # Quick guards: nothing to do if no atoms or no QApplication
        if not self.data.atoms:
            return

        # If called from non-GUI thread, schedule the heavy RDKit work here but
        # always perform UI mutations on the main thread via QTimer.singleShot.
        try:
            try:
                self._ih_update_counter += 1
            except (AttributeError, RuntimeError):
                self._ih_update_counter = getattr(self, "_ih_update_counter", 0) or 1
            my_token = self._ih_update_counter

            mol = None
            try:
                mol = self.data.to_rdkit_mol()
            except (AttributeError, RuntimeError):
                mol = None

            # Build a mapping of original_id -> hydrogen count without touching Qt items
            h_count_map = {}

            if mol is None:
                # Invalid/unsanitizable structure: reset all counts to 0
                for atom_id in list(self.data.atoms.keys()):
                    h_count_map[atom_id] = 0
            else:
                for atom in mol.GetAtoms():
                    try:
                        if not atom.HasProp("_original_atom_id"):
                            continue
                        original_id = atom.GetIntProp("_original_atom_id")

                        # Robust retrieval of H counts: prefer implicit, fallback to total or 0
                        try:
                            h_count = int(atom.GetNumImplicitHs())
                        except (AttributeError, RuntimeError):
                            try:
                                h_count = int(atom.GetTotalNumHs())
                            except (AttributeError, RuntimeError):
                                h_count = 0

                        h_count_map[int(original_id)] = h_count
                    except (AttributeError, RuntimeError):
                        # Skip problematic RDKit atoms
                        continue

            # Compute a per-atom problem map (original_id -> bool) so the
            # UI closure can safely set AtomItem.has_problem on the main thread.
            problem_map = {}
            try:
                if mol is not None:
                    try:
                        problems = Chem.DetectChemistryProblems(mol)
                    except (AttributeError, RuntimeError):
                        problems = None

                    if problems:
                        for prob in problems:
                            try:
                                atom_idx = prob.GetAtomIdx()
                                rd_atom = mol.GetAtomWithIdx(atom_idx)
                                if rd_atom and rd_atom.HasProp("_original_atom_id"):
                                    orig = int(rd_atom.GetIntProp("_original_atom_id"))
                                    problem_map[orig] = True
                            except (AttributeError, RuntimeError):
                                continue
                else:
                    # Fallback: use a lightweight valence heuristic similar to
                    # check_chemistry_problems_fallback() so we still flag atoms
                    # when RDKit conversion wasn't possible.
                    for atom_id, atom_data in self.data.atoms.items():
                        try:
                            symbol = atom_data.get("symbol")
                            charge = atom_data.get("charge", 0)
                            bond_count = 0
                            for (id1, id2), bond_data in self.data.bonds.items():
                                if id1 == atom_id or id2 == atom_id:
                                    bond_count += bond_data.get("order", 1)

                            if is_problematic_valence(symbol, bond_count, charge):
                                problem_map[atom_id] = True
                        except (AttributeError, RuntimeError):
                            continue
            except (AttributeError, RuntimeError):
                problem_map = {}

            def _apply_ui_updates():
                # If the global counter changed since this closure was
                # created, bail out — the update is stale.
                try:
                    if my_token != getattr(self, "_ih_update_counter", None):
                        return
                except (AttributeError, RuntimeError):
                    return

                try:
                    atoms_snapshot = dict(self.data.atoms)
                except (AttributeError, RuntimeError):
                    atoms_snapshot = {}
                is_deleted_func = sip_isdeleted_safe

                items_to_update = []
                for atom_id, atom_data in atoms_snapshot.items():
                    try:
                        item = atom_data.get("item")
                        if not item:
                            continue

                        # If sip.isdeleted is available, skip deleted C++ wrappers
                        try:
                            if is_deleted_func and is_deleted_func(item):
                                continue
                        except (AttributeError, RuntimeError):  # pragma: no cover
                            pass

                        # If the item is no longer in a scene, skip updating it to avoid
                        # touching partially-deleted objects during scene teardown.
                        try:
                            sc = item.scene() if hasattr(item, "scene") else None
                            if sc is None:
                                continue
                        except (AttributeError, RuntimeError):
                            # Accessing scene() might fail for a damaged object; skip it
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
                        try:
                            if need_geometry and hasattr(item, "prepareGeometryChange"):
                                try:
                                    item.prepareGeometryChange()
                                except (AttributeError, RuntimeError):  # pragma: no cover
                                    pass
                            # Apply implicit hydrogen count (guarded)
                            try:
                                item.implicit_h_count = new_count
                            except (AttributeError, RuntimeError):  # pragma: no cover
                                pass

                            # Apply problem flag (visual red-outline)
                            try:
                                item.has_problem = bool(desired_prob)
                            except (AttributeError, RuntimeError):  # pragma: no cover
                                pass
                            # Ensure the item is updated in the scene so paint() runs
                            # when either geometry or problem-flag changed.
                            items_to_update.append(item)
                        except (AttributeError, RuntimeError):
                            # Non-fatal: skip problematic items
                            continue

                    except (AttributeError, RuntimeError):
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
                            try:
                                it.update()
                            except (AttributeError, RuntimeError):  # pragma: no cover
                                pass
                    except (AttributeError, RuntimeError):
                        # Ignore any unexpected errors when touching the item
                        continue

            # Always schedule on main thread asynchronously
            try:
                QTimer.singleShot(0, _apply_ui_updates)
            except (AttributeError, RuntimeError):
                # Fallback: try to call directly (best-effort)
                try:
                    _apply_ui_updates()
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
        except (AttributeError, RuntimeError):  # pragma: no cover
            pass

    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")

        # Clear existing problem flags
        self.scene.clear_all_problem_flags()

        # Case: no atoms in 2D editor
        if not self.data.atoms:
            self.statusBar().showMessage("Error: No atoms to optimize.")
            return

        mol = self.data.to_rdkit_mol()
        if mol is None or mol.GetNumAtoms() == 0:
            # If RDKit conversion fails, check for chemistry problems
            self.check_chemistry_problems_fallback()
            return

        try:
            # Map atom IDs to RDKit coordinates
            view_center = self.view_2d.mapToScene(
                self.view_2d.viewport().rect().center()
            )
            new_positions_map = {}
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            for rdkit_atom in mol.GetAtoms():
                original_id = rdkit_atom.GetIntProp("_original_atom_id")
                new_positions_map[original_id] = conf.GetAtomPosition(
                    rdkit_atom.GetIdx()
                )

            if not new_positions_map:
                self.statusBar().showMessage(
                    "Optimization failed to generate coordinates."
                )
                return

            target_atom_items = [
                self.data.atoms[atom_id]["item"]
                for atom_id in new_positions_map.keys()
                if atom_id in self.data.atoms and "item" in self.data.atoms[atom_id]
            ]
            if not target_atom_items:
                self.statusBar().showMessage(
                    "Error: Atom items not found for optimized atoms."
                )
                return

            # Maintain original centroid
            # original_center_x = sum(item.pos().x() for item in target_atom_items) / len(target_atom_items)
            # original_center_y = sum(item.pos().y() for item in target_atom_items) / len(target_atom_items)

            positions = list(new_positions_map.values())
            rdkit_cx = sum(p.x for p in positions) / len(positions)
            rdkit_cy = sum(p.y for p in positions) / len(positions)

            SCALE = 50.0

            # Apply new coordinates
            for atom_id, rdkit_pos in new_positions_map.items():
                if atom_id in self.data.atoms:
                    item = self.data.atoms[atom_id]["item"]
                    sx = ((rdkit_pos.x - rdkit_cx) * SCALE) + view_center.x()
                    sy = (-(rdkit_pos.y - rdkit_cy) * SCALE) + view_center.y()
                    new_scene_pos = QPointF(sx, sy)
                    item.setPos(new_scene_pos)
                    self.data.atoms[atom_id]["pos"] = new_scene_pos

            # Update all bond positions
            # Guard against partially-deleted Qt wrappers: skip items that
            # SIP reports as deleted or which are no longer in a scene.
            for bond_data in self.data.bonds.values():
                item = bond_data.get("item") if bond_data else None
                if not item:
                    continue
                try:
                    # If SIP is available, skip wrappers whose C++ object is gone
                    if sip_isdeleted_safe(item):
                        continue
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
                try:
                    sc = None
                    try:
                        sc = item.scene() if hasattr(item, "scene") else None
                    except (AttributeError, RuntimeError):
                        sc = None
                    if sc is None:
                        continue
                    try:
                        item.update_position()
                    except (AttributeError, RuntimeError):
                        # Best-effort: skip any bond items that raise when updating
                        continue
                except (AttributeError, RuntimeError):
                    continue

            # Run overlap resolution
            self.resolve_overlapping_groups()

            # Update measurement labels
            self.update_2d_measurement_labels()

            # Request scene update
            self.scene.update()

            self.statusBar().showMessage("2D structure optimization successful.")
            self.push_undo_state()

        except (AttributeError, RuntimeError, ValueError) as e:
            self.statusBar().showMessage(f"Error during 2D optimization: {e}")
        finally:
            self.view_2d.setFocus()

    def redraw_molecule_3d(self):
        """Manually trigger redraw of the 3D molecule."""
        if hasattr(self, "current_mol") and self.current_mol:
            self.draw_molecule_3d(self.current_mol)
            self.statusBar().showMessage("Redraw complete.", 2000)
        else:
            self.statusBar().showMessage("No 3D molecule to redraw.")

    def resolve_overlapping_groups(self):
        """Detect and resolve overlapping atom groups."""

        # --- Parameters ---
        # Distance threshold for overlap
        OVERLAP_THRESHOLD = 0.5
        # Translation distance (bottom-left)
        MOVE_DISTANCE = 20

        # Safely retrieve item from self.data.atoms.values()
        all_atom_items = [
            data["item"] for data in self.data.atoms.values() if data and "item" in data
        ]

        if len(all_atom_items) < 2:
            return

        # Step 1: List overlapping pairs
        overlapping_pairs = []
        for item1, item2 in itertools.combinations(all_atom_items, 2):
            # Ignore directly bonded pairs
            if self.scene.find_bond_between(item1, item2):
                continue

            dist = QLineF(item1.pos(), item2.pos()).length()
            if dist < OVERLAP_THRESHOLD:
                overlapping_pairs.append((item1, item2))

        if not overlapping_pairs:
            self.statusBar().showMessage("No overlapping atoms found.", 2000)
            return

        # Step 2: Build overlap groups (Union-Find)
        # Track group membership
        parent = {item.atom_id: item.atom_id for item in all_atom_items}

        def find_set(atom_id):
            # Find group representative
            if parent[atom_id] == atom_id:
                return atom_id
            parent[atom_id] = find_set(parent[atom_id])
            return parent[atom_id]

        def unite_sets(id1, id2):
            # Unite two groups
            root1 = find_set(id1)
            root2 = find_set(id2)
            if root1 != root2:
                parent[root2] = root1

        for item1, item2 in overlapping_pairs:
            unite_sets(item1.atom_id, item2.atom_id)

        # Step 3: Plan translations per group
        # Group atoms by representative
        groups_by_root = {}
        for item in all_atom_items:
            root_id = find_set(item.atom_id)
            if root_id not in groups_by_root:
                groups_by_root[root_id] = []
            groups_by_root[root_id].append(item.atom_id)

        move_operations = []
        processed_roots = set()

        for root_id, group_atom_ids in groups_by_root.items():
            # Skip processed or single-atom groups
            if root_id in processed_roots or len(group_atom_ids) < 2:
                continue
            processed_roots.add(root_id)

            # 3a: Split group into fragments (BFS)
            fragments = []
            visited_in_group = set()
            group_atom_ids_set = set(group_atom_ids)

            for atom_id in group_atom_ids:
                if atom_id not in visited_in_group:
                    current_fragment = set()
                    q = deque([atom_id])
                    visited_in_group.add(atom_id)
                    current_fragment.add(atom_id)

                    while q:
                        current_id = q.popleft()
                        # Faster search if adjacency_list exists
                        for neighbor_id in self.data.adjacency_list.get(current_id, []):
                            if (
                                neighbor_id in group_atom_ids_set
                                and neighbor_id not in visited_in_group
                            ):
                                visited_in_group.add(neighbor_id)
                                current_fragment.add(neighbor_id)
                                q.append(neighbor_id)
                    fragments.append(current_fragment)

            if len(fragments) < 2:
                continue  # Skip if fragments don't overlap

            # 3b: Determine fragment to move
            # Find representative overlapping pair
            rep_item1, rep_item2 = None, None
            for i1, i2 in overlapping_pairs:
                if find_set(i1.atom_id) == root_id:
                    rep_item1, rep_item2 = i1, i2
                    break

            if not rep_item1:
                continue

            # Map pair to fragments
            frag1 = next((f for f in fragments if rep_item1.atom_id in f), None)
            frag2 = next((f for f in fragments if rep_item2.atom_id in f), None)

            # Skip overlaps within the same fragment
            if not frag1 or not frag2 or frag1 == frag2:
                continue

            # Move fragment with higher atom ID
            if rep_item1.atom_id > rep_item2.atom_id:
                ids_to_move = frag1
            else:
                ids_to_move = frag2

            # 3c: Plan translation
            translation_vector = QPointF(
                -MOVE_DISTANCE, MOVE_DISTANCE
            )  # Vector towards bottom-left
            move_operations.append((ids_to_move, translation_vector))

        # Step 4: Execute translations
        if not move_operations:
            self.statusBar().showMessage("No actionable overlaps found.", 2000)
            return

        for group_ids, vector in move_operations:
            for atom_id in group_ids:
                item = self.data.atoms[atom_id]["item"]
                new_pos = item.pos() + vector
                item.setPos(new_pos)
                self.data.atoms[atom_id]["pos"] = new_pos

        # Step 5: Update display and state
        for bond_data in self.data.bonds.values():
            item = bond_data.get("item") if bond_data else None
            if not item:
                continue
            try:
                if sip_isdeleted_safe(item):
                    continue
            except (AttributeError, RuntimeError):  # pragma: no cover
                pass

            try:
                sc = None
                try:
                    sc = item.scene() if hasattr(item, "scene") else None
                except (AttributeError, RuntimeError):
                    sc = None
                if sc is None:
                    continue
                try:
                    item.update_position()
                except (AttributeError, RuntimeError):
                    continue
            except (AttributeError, RuntimeError):
                continue

        # Update labels after resolution
        self.update_2d_measurement_labels()

        self.scene.update()
        self.push_undo_state()
        self.statusBar().showMessage("Resolved overlapping groups.", 2000)

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

    def _apply_chem_check_and_set_flags(self, mol, source_desc=None):
        """Central helper to apply chemical sanitization (or skip it) and set
        chem_check_tried / chem_check_failed flags consistently.

        When sanitization fails, a warning is shown and the Optimize 3D button
        is disabled. If the user setting 'skip_chemistry_checks' is True, no
        sanitization is attempted and both flags remain False.
        """
        try:
            self.chem_check_tried = False
            self.chem_check_failed = False
        except (AttributeError, RuntimeError):
            # Ensure attributes exist even if called very early
            self.chem_check_tried = False
            self.chem_check_failed = False

        if self.settings.get("skip_chemistry_checks", False):
            # User asked to skip chemistry checks entirely
            return

        try:
            Chem.SanitizeMol(mol)
            self.chem_check_tried = True
            self.chem_check_failed = False
        except (AttributeError, RuntimeError):
            # Mark that we tried sanitization and it failed
            self.chem_check_tried = True
            self.chem_check_failed = True
            try:
                desc = f" ({source_desc})" if source_desc else ""
                self.statusBar().showMessage(
                    f"Molecule sanitization failed{desc}; file may be malformed."
                )
            except (AttributeError, RuntimeError):  # pragma: no cover
                pass
            # Disable 3D optimization UI to prevent running on invalid molecules
            if hasattr(self, "optimize_3d_button"):
                try:
                    self.optimize_3d_button.setEnabled(False)
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass

    def _clear_xyz_flags(self, mol=None):
        """Clear XYZ-derived markers from a molecule (or current_mol) and
        reset UI flags accordingly.

        This is a best-effort cleanup to remove properties like
        _xyz_skip_checks and _xyz_atom_data that may have been attached when
        an XYZ file was previously loaded. After clearing molecule-level
        markers, the UI flag self.is_xyz_derived is set to False and the
        Optimize 3D button is re-evaluated (enabled unless chem_check_failed
        is True).
        """
        target = mol if mol is not None else getattr(self, "current_mol", None)
        try:
            if target is not None:
                # Remove RDKit property if present
                try:
                    if hasattr(target, "HasProp") and target.HasProp(
                        "_xyz_skip_checks"
                    ):
                        try:
                            target.ClearProp("_xyz_skip_checks")
                        except (AttributeError, RuntimeError):
                            try:
                                target.SetIntProp("_xyz_skip_checks", 0)
                            except (AttributeError, RuntimeError):  # pragma: no cover
                                pass
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
                # Remove attribute-style markers if present
                try:
                    if hasattr(target, "_xyz_skip_checks"):
                        try:
                            delattr(target, "_xyz_skip_checks")
                        except (AttributeError, RuntimeError):
                            try:
                                del target._xyz_skip_checks
                            except (AttributeError, RuntimeError):
                                try:
                                    target._xyz_skip_checks = False
                                except (AttributeError, RuntimeError):  # pragma: no cover
                                    pass
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
                try:
                    if hasattr(target, "_xyz_atom_data"):
                        try:
                            delattr(target, "_xyz_atom_data")
                        except (AttributeError, RuntimeError):
                            try:
                                del target._xyz_atom_data
                            except (AttributeError, RuntimeError):
                                try:
                                    target._xyz_atom_data = None
                                except (AttributeError, RuntimeError):  # pragma: no cover
                                    pass
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
        except (AttributeError, RuntimeError):  # pragma: no cover
            pass

        # Reset UI flags
        try:
            self.is_xyz_derived = False
        except (AttributeError, RuntimeError):  # pragma: no cover
            pass

        # Enable Optimize 3D unless sanitization failed
        try:
            if hasattr(self, "optimize_3d_button"):
                if getattr(self, "chem_check_failed", False):
                    try:
                        self.optimize_3d_button.setEnabled(False)
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        pass
                else:
                    try:
                        self.optimize_3d_button.setEnabled(True)
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        pass
        except (AttributeError, RuntimeError):  # pragma: no cover
            pass



