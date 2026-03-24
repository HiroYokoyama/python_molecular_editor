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
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin


class AlignPlaneDialog(Dialog3DPickingMixin, QDialog):  # pragma: no cover
    def __init__(self, mol, main_window, plane, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.plane = plane
        self.selected_atoms = set()

        # Add preselected atoms
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)

        self.init_ui()

        # Add labels to preselected atoms
        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    def init_ui(self):
        plane_names = {"xy": "XY", "xz": "XZ", "yz": "YZ"}
        self.setWindowTitle(f"Align to {plane_names[self.plane]} Plane")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            f"Click atoms in the 3D view to select them for align to the {plane_names[self.plane]} plane. At least 3 atoms are required."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip(
            "Select all atoms in the molecule for alignment"
        )
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)

        button_layout.addStretch()

        self.apply_button = QPushButton("Apply align")
        self.apply_button.clicked.connect(self.apply_PlaneAlign)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()

    def enable_picking(self):
        """Enable atom selection in the 3D view."""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True

    def disable_picking(self):
        """Disable atom selection in the 3D view."""
        if hasattr(self, "picking_enabled") and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False

    def on_atom_picked(self, atom_idx):
        """Handle the event when an atom is picked in the 3D view."""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)

        # Display labels on the atoms
        self.show_atom_labels()
        self.update_display()

    def keyPressEvent(self, event):
        """Handle keyboard shortcut events."""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_PlaneAlign()
            event.accept()
        else:
            super().keyPressEvent(event)

    def clear_selection(self):
        """Clear the current atom selection."""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def select_all_atoms(self):
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, "mol") and self.mol is not None:
                try:
                    n = self.mol.GetNumAtoms()
                    # create a set of indices [0..n-1]
                    self.selected_atoms = set(range(n))
                except (AttributeError, RuntimeError):
                    # fallback to main_window data map
                    self.selected_atoms = (
                        set(self.main_window.data.atoms.keys())
                        if hasattr(self.main_window, "data")
                        else set()
                    )
            else:
                # fallback to main_window data map
                self.selected_atoms = (
                    set(self.main_window.data.atoms.keys())
                    if hasattr(self.main_window, "data")
                    else set()
                )

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except (AttributeError, RuntimeError, TypeError, KeyError) as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def update_display(self):
        """Update the UI display with current selection info."""
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText(
                "Click atoms to select for align (minimum 3 required)"
            )
            self.apply_button.setEnabled(False)
        else:
            atom_list = sorted(self.selected_atoms)
            atom_display = []
            for i, atom_idx in enumerate(atom_list):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_display.append(f"#{i + 1}: {symbol}({atom_idx})")

            self.selection_label.setText(
                f"Selected {count} atoms: {', '.join(atom_display)}"
            )
            self.apply_button.setEnabled(count >= 3)

    def show_atom_labels(self):
        """Show numeric labels for the selected atoms in the 3D view."""
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            pairs = [(idx, f"#{i + 1}") for i, idx in enumerate(sorted_atoms)]
            self.show_atom_labels_for(pairs, color="blue")
        else:
            self.clear_atom_labels()

    def apply_PlaneAlign(self):
        """Apply plane alignment (rotation-based)."""
        if len(self.selected_atoms) < 3:
            QMessageBox.warning(
                self, "Warning", "Please select at least 3 atoms for align."
            )
            return
        try:
            # Get positions of selected atoms
            selected_indices = list(self.selected_atoms)
            selected_positions = self.main_window.atom_positions_3d[
                selected_indices
            ].copy()

            # Calculate centroid
            centroid = np.mean(selected_positions, axis=0)

            # Move centroid to origin
            centered_positions = selected_positions - centroid

            # Find optimal plane using PCA
            # Calculate covariance matrix for centered coordinates
            cov_matrix = np.cov(centered_positions.T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

            # Normal vector of the plane corresponds to the smallest eigenvalue
            normal_vector = eigenvectors[:, 0]  # Normal vector of the plane corresponds to the smallest eigenvalue

            # Define target plane normal vector
            if self.plane == "xy":
                target_normal = np.array([0, 0, 1])  # Z-axis direction
            elif self.plane == "xz":
                target_normal = np.array([0, 1, 0])  # Y-axis direction
            elif self.plane == "yz":
                target_normal = np.array([1, 0, 0])  # X-axis direction
            else:
                target_normal = np.array([0, 0, 1])  # Default to Z-axis (XY plane)

            # Adjust normal vector direction (ensure it's in the target direction)
            if np.dot(normal_vector, target_normal) < 0:
                normal_vector = -normal_vector

            # Calculate rotation axis and angle
            rotation_axis = np.cross(normal_vector, target_normal)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-10:  # If rotation is necessary
                rotation_axis = rotation_axis / rotation_axis_norm
                cos_angle = np.dot(normal_vector, target_normal)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                rotation_angle = np.arccos(cos_angle)

                # Rotate the entire molecule using Rodrigues' rotation formula
                def rodrigues_rotation(v, axis, angle):
                    cos_a = np.cos(angle)
                    sin_a = np.sin(angle)
                    return (
                        v * cos_a
                        + np.cross(axis, v) * sin_a
                        + axis * np.dot(axis, v) * (1 - cos_a)
                    )

                # Rotate all atoms
                conf = self.mol.GetConformer()
                for i in range(self.mol.GetNumAtoms()):
                    current_pos = np.array(conf.GetAtomPosition(i))
                    # Rotate relative to centroid
                    centered_pos = current_pos - centroid
                    rotated_pos = rodrigues_rotation(
                        centered_pos, rotation_axis, rotation_angle
                    )
                    new_pos = rotated_pos + centroid
                    conf.SetAtomPosition(i, new_pos.tolist())
                    self.main_window.atom_positions_3d[i] = new_pos

            # Update 3D visualization
            self.main_window.draw_molecule_3d(self.mol)

            # Update chirality labels
            self.main_window.update_chiral_labels()

            # Save state for Undo
            self.main_window.push_undo_state()

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(self, "Error", f"Failed to apply align: {str(e)}")

    def closeEvent(self, event):
        """Clean up labels and picking when closed."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """Handle cancellation or manual closing."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()

    def accept(self):
        """Handle acceptance (Apply/OK)."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
