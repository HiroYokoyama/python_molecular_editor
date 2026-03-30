#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import numpy as np
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)
from rdkit import Geometry

try:
    from .base_picking_dialog import BasePickingDialog
except ImportError:
    from moleditpy.ui.base_picking_dialog import BasePickingDialog


class PlanarizeDialog(BasePickingDialog):
    """Dialog to planarize a selected set of atoms by projecting them onto a best-fit plane."""

    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        super().__init__(mol, main_window, parent)
        self.selected_atoms = set()

        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)

        self.init_ui()

        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    def init_ui(self):
        self.setWindowTitle("Planarize")
        self.setModal(False)
        layout = QVBoxLayout(self)

        instruction_label = QLabel(
            "Click atoms in the 3D view to select them for planarization (minimum 3 required)."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        # Add Select All Atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip(
            "Select all atoms in the molecule for planarization"
        )
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)

        self.apply_button = QPushButton("Apply planarize")
        self.apply_button.clicked.connect(self.apply_planarize)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        button_layout.addStretch()

        layout.addLayout(button_layout)

        # enable picking
        self.picker_connection = None
        self.enable_picking()

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def clear_selection(self):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def update_display(self):
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText(
                "Click atoms to select for planarize (minimum 3 required)"
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

        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def show_atom_labels(self):
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            pairs = [(idx, f"#{i + 1}") for i, idx in enumerate(sorted_atoms)]
            self.show_atom_labels_for(pairs, color="cyan")
        else:
            self.clear_atom_labels()

    def apply_planarize(self):
        if not self.selected_atoms or len(self.selected_atoms) < 3:
            QMessageBox.warning(
                self, "Warning", "Please select at least 3 atoms for planarize."
            )
            return

        try:
            # Get positions of selected atoms
            selected_indices = list(sorted(self.selected_atoms))
            selected_positions = self.main_window.view_3d_manager.atom_positions_3d[
                selected_indices
            ].copy()

            centroid = np.mean(selected_positions, axis=0)
            centered_positions = selected_positions - centroid

            from moleditpy.core.mol_geometry import calculate_best_fit_plane_projection

            # Get normal of the least-squares plane via SVD
            u, s, vh = np.linalg.svd(centered_positions, full_matrices=False)
            normal = vh[-1]
            norm = np.linalg.norm(normal)
            if norm == 0:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "Cannot determine fit plane (degenerate positions).",
                )
                return
            normal = normal / norm

            new_positions = calculate_best_fit_plane_projection(
                centered_positions, normal, centroid
            )

            # Update molecular coordinates
            positions = self.mol.GetConformer().GetPositions()
            for i, new_pos in zip(selected_indices, new_positions):
                positions[int(i)] = new_pos
 
            # Write updated positions back using inherited helper
            self._update_molecule_geometry(positions)

            # Push Undo state AFTER modification
            self._push_undo()

            QMessageBox.information(
                self,
                "Success",
                f"Planarized {len(selected_indices)} atoms to best-fit plane.",
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to planarize: {e}")
