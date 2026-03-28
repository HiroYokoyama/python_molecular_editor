#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import numpy as np

try:
    from .mol_geometry import (
        calc_angle_deg,
        calc_distance,
        calculate_dihedral as _calculate_dihedral,
    )
except ImportError:
    from moleditpy_linux.core.mol_geometry import (
        calc_angle_deg,
        calc_distance,
        calculate_dihedral as _calculate_dihedral,
    )

# RDKit imports (explicit to satisfy flake8 and used features)
try:
    from . import sip_isdeleted_safe
except ImportError:
    from modules import sip_isdeleted_safe

# PyQt6 Modules
import pyvista as pv
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import QGraphicsTextItem

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import VDW_RADII
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy_linux.utils.constants import VDW_RADII


# --- Classes ---
class MainWindowEdit3d:
    _cls = None

    def toggle_measurement_mode(self, checked):
        """Toggle measurement mode on/off."""
        if checked:
            # Disable 3D Drag mode when measurement mode is on
            if self.is_3d_edit_mode:
                self.edit_3d_action.setChecked(False)
                self.toggle_3d_edit_mode(False)

            # Close active 3D edit dialogs
            self.close_all_3d_edit_dialogs()

        self.measurement_mode = checked

        if not checked:
            self.clear_measurement_selection()

        # Update status message
        if checked:
            self.statusBar().showMessage(
                "Measurement mode enabled. Click atoms to measure distances/angles/dihedrals."
            )
        else:
            self.statusBar().showMessage("Measurement mode disabled.")

    def close_all_3d_edit_dialogs(self):
        """Close all active 3D edit dialogs."""
        dialogs_to_close = self.active_3d_dialogs.copy()
        for dialog in dialogs_to_close:
            try:
                dialog.close()
            except (AttributeError, RuntimeError):
                # Suppress non-critical 3D edit/UI sync errors during bulk dialog teardown.
                # If a dialog is already closed or its C++ object is gone, we ignore it.
                pass

        self.active_3d_dialogs.clear()

    def handle_measurement_atom_selection(self, atom_idx):
        """Handle atom selection for measurement."""
        # Skip if already selected
        if atom_idx in self.selected_atoms_for_measurement:
            return

        self.selected_atoms_for_measurement.append(atom_idx)

        # Add atom labels
        self.add_measurement_label(atom_idx, len(self.selected_atoms_for_measurement))

        # Calculate and display results
        self.calculate_and_display_measurements()

    def add_measurement_label(self, atom_idx, label_number):
        """Add numeric labels to atoms."""
        if not self.current_mol or atom_idx >= self.current_mol.GetNumAtoms():
            return

        # Update label list
        self.measurement_labels.append((atom_idx, str(label_number)))

        # Redraw 3D measurement labels
        self.update_measurement_labels_display()

        # Update 2D measurement labels
        self.update_2d_measurement_labels()

    def update_measurement_labels_display(self):
        """Draw measurement labels in 3D (atom centers)."""
        try:
            # Remove existing labels
            self.plotter.remove_actor("measurement_labels")
        except (AttributeError, RuntimeError):
            # Suppress if the actor is already destroyed or not found.
            pass

        if not self.measurement_labels or not self.current_mol:
            return

        # Prepare label positions and text
        pts, labels = [], []
        for atom_idx, label_text in self.measurement_labels:
            if atom_idx < len(self.atom_positions_3d):
                coord = self.atom_positions_3d[atom_idx].copy()
                # Place at atom center
                pts.append(coord)
                labels.append(label_text)

        if pts and labels:
            # Use PyVista's point_labels
            self.plotter.add_point_labels(
                np.array(pts),
                labels,
                font_size=16,
                point_size=0,
                text_color="red",  # Always red for measurement
                name="measurement_labels",
                always_visible=True,
                tolerance=0.01,
                show_points=False,
            )

    def clear_measurement_selection(self):
        """Clear measurement selection."""
        self.selected_atoms_for_measurement.clear()

        # Remove 3D labels
        self.measurement_labels.clear()
        try:
            self.plotter.remove_actor("measurement_labels")
        except (AttributeError, RuntimeError):
            # Suppress if the actor is already destroyed or not found.
            pass

        # Remove 2D labels
        self.clear_2d_measurement_labels()

        # Remove result text
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
                self.measurement_text_actor = None
            except (AttributeError, RuntimeError):
                # Suppress if the actor is already destroyed or not found.
                pass

        self.plotter.render()

    def update_2d_measurement_labels(self):
        """Update 2D measurement labels."""
        # Remove existing 2D labels
        self.clear_2d_measurement_labels()

        # Create atom-to-AtomItem mapping
        if not self.current_mol or not hasattr(self, "data") or not self.data.atoms:
            return

        # Map RDKit index to 2D AtomItem
        atom_idx_to_item = {}

        # Get AtomItems from scene
        if hasattr(self, "scene"):
            for item in self.scene.items():
                if hasattr(item, "atom_id") and hasattr(
                    item, "symbol"
                ):  # Check if AtomItem
                    # Find RDKit index from atom ID
                    rdkit_idx = self.find_rdkit_atom_index(item)
                    if rdkit_idx is not None:
                        atom_idx_to_item[rdkit_idx] = item

        # Add to 2D view
        if not hasattr(self, "measurement_label_items_2d"):
            self.measurement_label_items_2d = []

        for atom_idx, label_text in self.measurement_labels:
            if atom_idx in atom_idx_to_item:
                atom_item = atom_idx_to_item[atom_idx]
                self.add_2d_measurement_label(atom_item, label_text)

    def add_2d_measurement_label(self, atom_item, label_text):
        """Add measurement label to specific AtomItem."""
        # Create label item
        label_item = QGraphicsTextItem(label_text)
        label_item.setDefaultTextColor(QColor(255, 0, 0))  # Red
        label_item.setFont(QFont("Arial", 12, QFont.Weight.Bold))

        # Set Z-value for top-most display
        label_item.setZValue(2000)  # Ensure it stays on top

        # Position near top-right of atom
        atom_pos = atom_item.pos()
        atom_rect = atom_item.boundingRect()
        label_pos = QPointF(
            atom_pos.x() + atom_rect.width() / 4 + 2,
            atom_pos.y() - atom_rect.height() / 4 - 8,
        )
        label_item.setPos(label_pos)

        # Add to scene
        self.scene.addItem(label_item)
        self.measurement_label_items_2d.append(label_item)

    def clear_2d_measurement_labels(self):
        """Remove all 2D measurement labels."""
        if hasattr(self, "measurement_label_items_2d"):
            for label_item in self.measurement_label_items_2d:
                try:
                    # Avoid touching partially-deleted wrappers
                    if sip_isdeleted_safe(label_item):
                        continue
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except (AttributeError, RuntimeError):
                        # Scene access or removal failed; skip this item.
                        pass
                except (AttributeError, RuntimeError):
                    # If sip check itself fails, fall back to best-effort removal
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        # Best-effort removal failed after sip check failed; skip.
                        continue
            self.measurement_label_items_2d.clear()

    def find_rdkit_atom_index(self, atom_item):
        """Find RDKit index from AtomItem."""
        if not self.current_mol or not atom_item:
            return None

        # Use mapping dictionary
        if (
            hasattr(self, "atom_id_to_rdkit_idx_map")
            and atom_item.atom_id in self.atom_id_to_rdkit_idx_map
        ):
            return self.atom_id_to_rdkit_idx_map[atom_item.atom_id]

        # Return None if no mapping exists
        return None

    def calculate_and_display_measurements(self):
        """Calculate and display measurement values."""
        num_selected = len(self.selected_atoms_for_measurement)
        if num_selected < 2:
            return

        measurement_text = []

        if num_selected >= 2:
            # Distance calculation
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            distance = self.calculate_distance(atom1_idx, atom2_idx)
            measurement_text.append(f"Distance 1-2: {distance:.3f} Å")

        if num_selected >= 3:
            # Angle calculation
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            atom3_idx = self.selected_atoms_for_measurement[2]
            angle = self.calculate_angle(atom1_idx, atom2_idx, atom3_idx)
            measurement_text.append(f"Angle 1-2-3: {angle:.2f}°")

        if num_selected >= 4:
            # Dihedral calculation
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            atom3_idx = self.selected_atoms_for_measurement[2]
            atom4_idx = self.selected_atoms_for_measurement[3]
            dihedral = self.calculate_dihedral(
                atom1_idx, atom2_idx, atom3_idx, atom4_idx
            )
            measurement_text.append(f"Dihedral 1-2-3-4: {dihedral:.2f}°")

        # Display results on 3D view
        self.display_measurement_text(measurement_text)

    def calculate_distance(self, atom1_idx, atom2_idx):
        """Calculate distance between two atoms."""
        return calc_distance(
            self.atom_positions_3d[atom1_idx], self.atom_positions_3d[atom2_idx]
        )

    def calculate_angle(self, atom1_idx, atom2_idx, atom3_idx):
        """Calculate angle (center is vertex)."""
        return calc_angle_deg(
            self.atom_positions_3d[atom1_idx],
            self.atom_positions_3d[atom2_idx],  # vertex
            self.atom_positions_3d[atom3_idx],
        )

    def calculate_dihedral(self, atom1_idx, atom2_idx, atom3_idx, atom4_idx):
        """Calculate dihedral angle."""
        return _calculate_dihedral(
            self.atom_positions_3d, atom1_idx, atom2_idx, atom3_idx, atom4_idx
        )

    def display_measurement_text(self, measurement_lines):
        """Display results text on 3D view."""
        # Remove existing text
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Suppress non-critical 3D edit/UI sync errors if the plotter or actor is already destroyed
                pass

        if not measurement_lines:
            self.measurement_text_actor = None
            return

        # Combine text
        text = "\n".join(measurement_lines)

        # Determine text color from background
        try:
            bg_color_hex = self.settings.get("background_color", "#919191")
            bg_qcolor = QColor(bg_color_hex)
            if bg_qcolor.isValid():
                luminance = bg_qcolor.toHsl().lightness()
                text_color = "black" if luminance > 128 else "white"
            else:
                text_color = "white"
        except (AttributeError, RuntimeError, ValueError, TypeError):
            # Fallback for determining text contrast; suppress if settings or plotter state is inconsistent.
            text_color = "white"

        # Display upper-left
        self.measurement_text_actor = self.plotter.add_text(
            text,
            position="upper_left",
            font_size=10,  # Smaller font
            color=text_color,  # Color matching background
            font="courier",  # Monospace font
            name="measurement_display",
        )

        self.plotter.render()

    def toggle_atom_selection_3d(self, atom_idx):
        """Toggle atom selection in 3D."""
        if atom_idx in self.selected_atoms_3d:
            self.selected_atoms_3d.remove(atom_idx)
        else:
            self.selected_atoms_3d.add(atom_idx)

        # Update feedback
        self.update_3d_selection_display()

    def clear_3d_selection(self):
        """Clear 3D selection."""
        self.selected_atoms_3d.clear()
        self.update_3d_selection_display()

    def update_3d_selection_display(self):
        """Update 3D selection highlight."""
        try:
            # Remove existing highlight
            self.plotter.remove_actor("selection_highlight")
        except (AttributeError, RuntimeError, ValueError, TypeError):
            # Suppress non-critical UI/rendering/measurement noise if the plotter or actor is already destroyed.
            pass

        if not self.selected_atoms_3d or not self.current_mol:
            self.plotter.render()
            return

        # Create index list
        selected_indices = list(self.selected_atoms_3d)

        # Get atom positions
        selected_positions = self.atom_positions_3d[selected_indices]

        # Highlight with slightly larger radius
        selected_radii = np.array(
            [
                VDW_RADII.get(self.current_mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3
                for i in selected_indices
            ]
        )

        # Create highlight dataset
        highlight_source = pv.PolyData(selected_positions)
        highlight_source["radii"] = selected_radii

        # Yellow semi-transparent highlight
        highlight_glyphs = highlight_source.glyph(
            scale="radii",
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
            orient=False,
        )

        self.plotter.add_mesh(
            highlight_glyphs, color="yellow", opacity=0.3, name="selection_highlight"
        )

        self.plotter.render()

    def remove_dialog_from_list(self, dialog):
        """Remove dialog from active list."""
        if dialog in self.active_3d_dialogs:
            self.active_3d_dialogs.remove(dialog)


MainWindowEdit3d._cls = MainWindowEdit3d
