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
from PyQt6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)
from rdkit import Geometry

try:
    from .dialog_3d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from moleditpy.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin


class AlignmentDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, axis, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.axis = axis
        self.selected_atoms = set()

        # Add preselected atoms (maximum 2)
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms[:2])

        self.init_ui()

        # Add labels to preselected atoms
        if self.selected_atoms:
            for i, atom_idx in enumerate(sorted(self.selected_atoms), 1):
                self.add_selection_label(atom_idx, f"Atom {i}")
            self.update_display()

    def init_ui(self):
        axis_names = {"x": "X-axis", "y": "Y-axis", "z": "Z-axis"}
        self.setWindowTitle(f"Align to {axis_names[self.axis]}")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            f"Click atoms in the 3D view to select them for alignment to the {axis_names[self.axis]}. Exactly 2 atoms are required. The first atom will be moved to the origin, and the second atom will be positioned on the {axis_names[self.axis]}."
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

        button_layout.addStretch()

        self.apply_button = QPushButton("Apply Alignment")
        self.apply_button.clicked.connect(self.apply_alignment)
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
        """Enable atom selection in the 3D window."""
        # Use functionality from Dialog3DPickingMixin
        super().enable_picking()

    def disable_picking(self):
        """Disable atom selection in the 3D window."""
        # Use functionality from Dialog3DPickingMixin
        super().disable_picking()

    def on_atom_picked(self, atom_idx):
        """Handle atom selection event."""
        if self.main_window.current_mol is None:
            return

        if atom_idx in self.selected_atoms:
            # De-select if already selected
            self.selected_atoms.remove(atom_idx)
            self.remove_atom_label(atom_idx)
        else:
            # Maximum of 2 atoms can be selected
            if len(self.selected_atoms) < 2:
                self.selected_atoms.add(atom_idx)
                # Show label indicating selection order
                label_text = f"Atom {len(self.selected_atoms)}"
                self.add_selection_label(atom_idx, label_text)

        self.update_display()

    def update_display(self):
        """Update the UI based on current selection."""
        if len(self.selected_atoms) == 0:
            self.selection_label.setText(
                "Click atoms to select for alignment (exactly 2 required)"
            )
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 1:
            selected_list = list(self.selected_atoms)
            atom = self.mol.GetAtomWithIdx(selected_list[0])
            self.selection_label.setText(
                f"Selected 1 atom: {atom.GetSymbol()}{selected_list[0] + 1}"
            )
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 2:
            selected_list = sorted(list(self.selected_atoms))
            atom1 = self.mol.GetAtomWithIdx(selected_list[0])
            atom2 = self.mol.GetAtomWithIdx(selected_list[1])
            self.selection_label.setText(
                f"Selected 2 atoms: {atom1.GetSymbol()}{selected_list[0] + 1}, {atom2.GetSymbol()}{selected_list[1] + 1}"
            )
            self.apply_button.setEnabled(True)

    def clear_selection(self):
        """Clear the current selection and labels."""
        self.clear_selection_labels()
        self.selected_atoms.clear()
        self.update_display()

    def remove_atom_label(self, atom_idx):
        """Remove a label for a specific atom."""
        # Re-draw all labels for simplicity
        self.clear_selection_labels()
        for i, idx in enumerate(sorted(self.selected_atoms), 1):
            if idx != atom_idx:
                self.add_selection_label(idx, f"Atom {i}")

    def apply_alignment(self):
        """Apply the specific axial alignment to the molecule."""
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(
                self, "Warning", "Please select exactly 2 atoms for alignment."
            )
            return
        try:
            selected_list = sorted(list(self.selected_atoms))
            atom1_idx, atom2_idx = selected_list[0], selected_list[1]

            conf = self.mol.GetConformer()

            # Get current atom positions
            pos1 = np.array(conf.GetAtomPosition(atom1_idx))
            pos2 = np.array(conf.GetAtomPosition(atom2_idx))

            # Translate entire molecule so atom1 is at the origin
            translation = -pos1
            for i in range(self.mol.GetNumAtoms()):
                current_pos = np.array(conf.GetAtomPosition(i))
                new_pos = current_pos + translation
                conf.SetAtomPosition(i, new_pos.tolist())

            # Get new position of atom2 after translation
            pos2_translated = pos2 + translation

            # Calculate rotation to align atom2 relative to the chosen axis
            axis_vectors = {
                "x": np.array([1.0, 0.0, 0.0]),
                "y": np.array([0.0, 1.0, 0.0]),
                "z": np.array([0.0, 0.0, 1.0]),
            }
            target_axis = axis_vectors[self.axis]

            # Direction vector from origin to translated atom2
            current_vector = pos2_translated
            current_length = np.linalg.norm(current_vector)

            if current_length > 1e-10:  # If not a zero vector
                current_vector_normalized = current_vector / current_length

                # Calculate rotation axis and angle
                rotation_axis = np.cross(current_vector_normalized, target_axis)
                rotation_axis_length = np.linalg.norm(rotation_axis)

                if rotation_axis_length > 1e-10:  # Rotation required
                    rotation_axis = rotation_axis / rotation_axis_length
                    cos_angle = np.dot(current_vector_normalized, target_axis)
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    rotation_angle = np.arccos(cos_angle)

                    # Use Rodrigues' rotation formula
                    def rodrigues_rotation(v, k, theta):
                        cos_theta = np.cos(theta)
                        sin_theta = np.sin(theta)
                        return (
                            v * cos_theta
                            + np.cross(k, v) * sin_theta
                            + k * np.dot(k, v) * (1 - cos_theta)
                        )

                    # Apply rotation to all atoms
                    for i in range(self.mol.GetNumAtoms()):
                        current_pos = np.array(conf.GetAtomPosition(i))
                        rotated_pos = rodrigues_rotation(
                            current_pos, rotation_axis, rotation_angle
                        )
                        conf.SetAtomPosition(
                            i,
                            Geometry.Point3D(
                                float(rotated_pos[0]),
                                float(rotated_pos[1]),
                                float(rotated_pos[2]),
                            ),
                        )

            # Update 3D positions
            self.main_window.atom_positions_3d = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(self.mol.GetNumAtoms())]
            )

            # Update 3D visualization
            self.main_window.draw_molecule_3d(self.mol)

            # Update chirality labels
            self.main_window.update_chiral_labels()

            # Save state for Undo
            self.main_window.push_undo_state()

            QMessageBox.information(
                self, "Success", f"Alignment to {self.axis.upper()}-axis completed."
            )

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(self, "Error", f"Failed to apply alignment: {str(e)}")

    def closeEvent(self, event):
        """Clean up when the dialog is closed."""
        self.clear_selection_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """Handle dialog rejection/cancellation."""
        self.clear_selection_labels()
        self.disable_picking()
        super().reject()

    def accept(self):
        """Handle dialog acceptance (OK/Apply)."""
        self.clear_selection_labels()
        self.disable_picking()
        super().accept()
