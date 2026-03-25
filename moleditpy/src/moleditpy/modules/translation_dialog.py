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
    QCheckBox,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from .dialog3_d_picking_mixin import Dialog3DPickingMixin


class TranslationDialog(Dialog3DPickingMixin, QDialog):  
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()  # For selecting multiple atoms
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Translation")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click atoms in the 3D view to select them. The centroid of selected atoms will be moved to the target coordinates, translating the entire molecule."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Coordinate inputs
        coord_layout = QGridLayout()
        coord_layout.addWidget(QLabel("Target X:"), 0, 0)
        self.x_input = QLineEdit("0.0")
        coord_layout.addWidget(self.x_input, 0, 1)

        coord_layout.addWidget(QLabel("Target Y:"), 1, 0)
        self.y_input = QLineEdit("0.0")
        coord_layout.addWidget(self.y_input, 1, 1)

        coord_layout.addWidget(QLabel("Target Z:"), 2, 0)
        self.z_input = QLineEdit("0.0")
        coord_layout.addWidget(self.z_input, 2, 1)

        layout.addLayout(coord_layout)

        # Translation target toggle: Entire molecule (default) or Selected atoms only
        self.translate_selected_only_checkbox = QCheckBox(
            "Translate selected atoms only"
        )
        self.translate_selected_only_checkbox.setToolTip(
            "When checked, only the atoms you selected will be moved so their centroid matches the target.\n"
            "When unchecked (default), the entire molecule will be translated so the selected atoms' centroid moves to the target."
        )
        self.translate_selected_only_checkbox.setChecked(
            False
        )  # default: entire molecule
        layout.addWidget(self.translate_selected_only_checkbox)

        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip(
            "Select all atoms in the molecule for translation"
        )
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)

        button_layout.addStretch()

        self.apply_button = QPushButton("Apply Translation")
        self.apply_button.clicked.connect(self.apply_translation)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()

    def on_atom_picked(self, atom_idx):
        """Handle the event when an atom is picked in the 3D view."""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def keyPressEvent(self, event):
        """Handle keyboard shortcut events."""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_translation()
            event.accept()
        else:
            super().keyPressEvent(event)

    def update_display(self):
        """Update the UI labels and button state based on the current selection."""
        if not self.selected_atoms:
            self.selection_label.setText("No atoms selected")
            self.apply_button.setEnabled(False)
        else:
            # Check molecule validity
            if not self.mol or self.mol.GetNumConformers() == 0:
                self.selection_label.setText("Error: No valid molecule or conformer")
                self.apply_button.setEnabled(False)
                return

            try:
                conf = self.mol.GetConformer()
                # Calculate the centroid of selected atoms
                centroid = self.calculate_centroid()

                # Display information about selected atoms
                atom_info = []
                for atom_idx in sorted(self.selected_atoms):
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    atom_info.append(f"{symbol}({atom_idx})")

                self.selection_label.setText(
                    f"Selected atoms: {', '.join(atom_info)}\n"
                    f"Centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})"
                )
                self.apply_button.setEnabled(True)
            except (AttributeError, RuntimeError, ValueError) as e:
                self.selection_label.setText(f"Error accessing atom data: {str(e)}")
                self.apply_button.setEnabled(False)

        # Update the coordinate input fields when selection changes
        if self.selected_atoms:
            try:
                coords = centroid
            except NameError:
                coords = self.calculate_centroid()
            self.x_input.setText(f"{coords[0]:.4f}")
            self.y_input.setText(f"{coords[1]:.4f}")
            self.z_input.setText(f"{coords[2]:.4f}")
        else:
            self.x_input.setText("0.0")
            self.y_input.setText("0.0")
            self.z_input.setText("0.0")

    def calculate_centroid(self):
        """Calculate the geometric center (centroid) of the selected atoms."""
        if not self.selected_atoms:
            return np.array([0.0, 0.0, 0.0])

        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.selected_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])

        return np.mean(positions, axis=0)

        
    def apply_translation(self):
        """Apply the translation to either the selected atoms or the entire molecule."""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select at least one atom.")
            return

        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Warning", "No valid molecule or conformer available.")
            return

        try:
            target_x = float(self.x_input.text())
            target_y = float(self.y_input.text())
            target_z = float(self.z_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid coordinates.")
            return

        try:
            # Calculate translation
            current_centroid = self.calculate_centroid()
            target_pos = np.array([target_x, target_y, target_z])
            translation_vector = target_pos - current_centroid

            conf = self.mol.GetConformer()
            atom_positions = getattr(self.main_window, "atom_positions_3d", {})

            # Apply translation
            translate_selected = self.translate_selected_only_checkbox.isChecked()
            for i in range(self.mol.GetNumAtoms()):
                if not translate_selected or i in self.selected_atoms:
                    atom_pos = np.array(conf.GetAtomPosition(i))
                    new_pos = atom_pos + translation_vector
                    conf.SetAtomPosition(i, new_pos.tolist())
                    
                    # Update cache in main window
                    if i in atom_positions:
                        atom_positions[i] = new_pos

            # Update visualization and state
            if hasattr(self.main_window, "draw_molecule_3d"):
                self.main_window.draw_molecule_3d(self.mol)
            
            if hasattr(self.main_window, "update_chiral_labels"):
                self.main_window.update_chiral_labels()

            self.clear_selection()

            if hasattr(self.main_window, "push_undo_state"):
                self.main_window.push_undo_state()

        except (AttributeError, RuntimeError, ValueError) as e:
            pass  # Suppress errors during translation application
            QMessageBox.critical(self, "Error", f"Failed to apply translation: {str(e)}")

    def clear_selection(self):
        """Clear the current atom selection and labels."""
        self.selected_atoms.clear()
        self.clear_atom_labels()

    def select_all_atoms(self):
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            if hasattr(self, "mol") and self.mol is not None:
                self.selected_atoms = set(range(self.mol.GetNumAtoms()))
            else:
                self.selected_atoms = (
                    set(self.main_window.data.atoms.keys())
                    if hasattr(self.main_window, "data")
                    else set()
                )
            self.show_atom_labels()
            self.update_display()
        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def show_atom_labels(self):
        """Display geometric labels for the selected atoms in the 3D window."""
        # Clear existing labels
        self.clear_atom_labels()

        if not hasattr(self, "selection_labels"):
            self.selection_labels = []

        if self.selected_atoms:
            positions = []
            labels = []

            for i, atom_idx in enumerate(sorted(self.selected_atoms)):
                pos = self.main_window.atom_positions_3d[atom_idx]
                positions.append(pos)
                labels.append(f"S{i + 1}")

            # Indicate the centroid position
            if len(self.selected_atoms) > 1:
                centroid = self.calculate_centroid()
                positions.append(centroid)
                labels.append("CEN")

            # Add labels to the plotter
            if positions:
                label_actor = self.main_window.plotter.add_point_labels(
                    positions,
                    labels,
                    point_size=20,
                    font_size=12,
                    text_color="cyan",
                    always_visible=True,
                )
                # Handle case where add_point_labels returns a list
                if isinstance(label_actor, list):
                    self.selection_labels.extend(label_actor)
                else:
                    self.selection_labels.append(label_actor)

    def clear_atom_labels(self):
        """Clear atom labels and force a re-render of the 3D scene."""
        super().clear_atom_labels()
        
        # Force re-render
        if hasattr(self.main_window, "plotter"):
            try:
                self.main_window.plotter.render()
            except (RuntimeError, ValueError, TypeError):
                # Suppress non-critical error
                pass
    def closeEvent(self, event):
        """Clean up when the dialog is closed directly."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """Clean up when the dialog is rejected (e.g., Cancel clicked)."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()

    def accept(self):
        """Clean up when the dialog is accepted (e.g., OK clicked)."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
