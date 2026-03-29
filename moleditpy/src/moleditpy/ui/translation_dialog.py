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
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from rdkit import Geometry

try:
    from .base_picking_dialog import BasePickingDialog
except ImportError:
    from moleditpy.ui.base_picking_dialog import BasePickingDialog


class TranslationDialog(BasePickingDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        super().__init__(mol, main_window, parent)
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
        self.setWindowTitle("Translate Atoms")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click atoms in the 3D view to select them for translation. Specify the translation vector in Å."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Translation vector input
        vector_layout = QHBoxLayout()
        vector_layout.addWidget(QLabel("dX:"))
        self.dx_input = QLineEdit("0.0")
        vector_layout.addWidget(self.dx_input)

        vector_layout.addWidget(QLabel("dY:"))
        self.dy_input = QLineEdit("0.0")
        vector_layout.addWidget(self.dy_input)

        vector_layout.addWidget(QLabel("dZ:"))
        self.dz_input = QLineEdit("0.0")
        vector_layout.addWidget(self.dz_input)

        layout.addLayout(vector_layout)

        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule")
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

        # Display labels on the atoms
        self.show_atom_labels()
        self.update_display()

    def clear_selection(self):
        """Clear the current atom selection."""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def select_all_atoms(self):
        """Select all atoms in the current molecule."""
        try:
            if hasattr(self, "mol") and self.mol is not None:
                n = self.mol.GetNumAtoms()
                self.selected_atoms = set(range(n))
            else:
                self.selected_atoms = (
                    set(self.main_window.data.atoms.keys())
                    if hasattr(self.main_window, "data")
                    else set()
                )

            self.show_atom_labels()
            self.update_display()
        except (AttributeError, RuntimeError, TypeError, KeyError) as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def update_display(self):
        """Update the UI display with current selection info."""
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select (minimum 1 required)")
            self.apply_button.setEnabled(False)
        else:
            self.selection_label.setText(f"Selected {count} atoms")
            self.apply_button.setEnabled(True)

    def show_atom_labels(self):
        """Show numeric labels for the selected atoms."""
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            pairs = [(idx, str(i + 1)) for i, idx in enumerate(sorted_atoms)]
            self.show_atom_labels_for(pairs)
        else:
            self.clear_atom_labels()

    def apply_translation(self):
        """Apply the translation to selected atoms."""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select at least one atom.")
            return

        try:
            dx = float(self.dx_input.text())
            dy = float(self.dy_input.text())
            dz = float(self.dz_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid numbers for dx, dy, dz.")
            return

        if dx == 0 and dy == 0 and dz == 0:
            return

        conf = self.mol.GetConformer()
        translation_vec = np.array([dx, dy, dz])

        # Push Undo state
        if hasattr(self.main_window, "state_manager"):
            self.main_window.state_manager.push_undo_state()
        elif hasattr(self.main_window, "edit_actions_manager"):
            self.main_window.edit_actions_manager.push_undo_state()

        # Update positions
        for atom_idx in self.selected_atoms:
            old_pos = np.array(conf.GetAtomPosition(atom_idx))
            new_pos = old_pos + translation_vec
            conf.SetAtomPosition(
                atom_idx,
                Geometry.Point3D(float(new_pos[0]), float(new_pos[1]), float(new_pos[2])),
            )
            self.main_window.view_3d_manager.atom_positions_3d[atom_idx] = new_pos

        # Redraw
        self.main_window.view_3d_manager.draw_molecule_3d(self.mol)

        # Update labels
        self.show_atom_labels()
