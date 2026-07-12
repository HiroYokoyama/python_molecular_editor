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
from typing import TYPE_CHECKING, Any, Literal, Optional, Sequence

import numpy as np

from PyQt6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

from .base_picking_dialog import BasePickingDialog, SelectionList

from ..core.mol_geometry import rodrigues_rotate

if TYPE_CHECKING:
    from .main_window import MainWindow


class AlignPlaneDialog(BasePickingDialog):
    """Dialog for aligning selected atoms to a principal plane (XY, XZ, or YZ)."""

    def __init__(
        self,
        mol: Chem.Mol,
        main_window: "MainWindow",
        plane: Literal["xy", "xz", "yz"],
        preselected_atoms: Optional[Sequence[int]] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(mol, main_window, parent)
        self.plane = plane
        self._selected_atoms = SelectionList()

        # Add preselected atoms
        if preselected_atoms:
            self._selected_atoms.update(preselected_atoms)

        self.init_ui()

        # Add labels to preselected atoms
        if self._selected_atoms:
            self.show_atom_labels()
            self.update_display()

    @property
    def selected_atoms(self) -> SelectionList:
        """Return the ordered list of selected atom indices."""
        return self._selected_atoms

    @selected_atoms.setter
    def selected_atoms(self, val: Any) -> None:
        """Replace the selection with a new SelectionList built from val."""
        self._selected_atoms = SelectionList(val)

    def init_ui(self) -> None:
        """Build and lay out all widgets for the plane-alignment dialog."""
        plane_names = {"xy": "XY", "xz": "XZ", "yz": "YZ"}
        self.setWindowTitle(f"Align to {plane_names[self.plane]} Plane")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            f"Click atoms in the 3D view to select them for align to "
            f"the {plane_names[self.plane]} plane. At least 3 atoms "
            f"are required."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Move to zero plane option (default False)
        self.move_to_zero_plane_checkbox = QCheckBox("Move the plane to the zero plane")
        self.move_to_zero_plane_checkbox.setChecked(False)
        layout.addWidget(self.move_to_zero_plane_checkbox)

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

    def on_atom_picked(self, atom_idx: int) -> None:
        """Handle the event when an atom is picked in the 3D view."""
        if atom_idx in self._selected_atoms:
            self._selected_atoms.remove(atom_idx)
        else:
            self._selected_atoms.append(atom_idx)

        self.show_atom_labels()
        self.update_display()

    def clear_selection(self) -> None:
        """Clear the current atom selection."""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def select_all_atoms(self) -> None:
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, "mol") and self.mol is not None:
                n = self.mol.GetNumAtoms()
                self.selected_atoms = set(range(n))
            else:
                # fallback to main_window data map
                self.selected_atoms = (
                    set(self.main_window.state_manager.data.atoms.keys())
                    if hasattr(self.main_window.state_manager, "data")
                    else set()
                )

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except (AttributeError, RuntimeError, TypeError, KeyError) as e:
            logging.exception("Failed to select all atoms")
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def update_display(self) -> None:
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

    def show_atom_labels(self) -> None:
        """Show numeric labels for the selected atoms in the 3D view."""
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            pairs = [(idx, f"#{i + 1}") for i, idx in enumerate(sorted_atoms)]
            self.show_atom_labels_for(pairs, color="yellow")
        else:
            self.clear_atom_labels()

    def apply_PlaneAlign(self) -> None:
        """Apply plane alignment (rotation-based)."""
        if len(self.selected_atoms) < 3:
            QMessageBox.warning(
                self,
                "Warning",
                "Please select at least 3 atoms for align.",
            )
            return
        try:
            # Get positions of selected atoms
            selected_indices = list(self.selected_atoms)
            positions = self.mol.GetConformer().GetPositions()
            selected_positions = positions[selected_indices]

            # Calculate centroid
            centroid = np.mean(selected_positions, axis=0)

            # Move centroid to origin
            centered_positions = selected_positions - centroid

            # Find optimal plane using PCA
            cov_matrix = np.cov(centered_positions.T)
            _, eigenvectors = np.linalg.eigh(cov_matrix)

            # Normal vector of the plane corresponds to the smallest eigenvalue
            normal_vector = eigenvectors[:, 0]

            # Define target plane normal vector
            if self.plane == "xy":
                target_normal = np.array([0, 0, 1])  # Z-axis direction
            elif self.plane == "xz":
                target_normal = np.array([0, 1, 0])  # Y-axis direction
            elif self.plane == "yz":
                target_normal = np.array([1, 0, 0])  # X-axis direction
            else:
                # Default to Z-axis (XY plane)
                target_normal = np.array([0, 0, 1])

            # Adjust normal vector direction
            if np.dot(normal_vector, target_normal) < 0:
                normal_vector = -normal_vector

            # Calculate rotation axis and angle
            rotation_axis = np.cross(normal_vector, target_normal)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            # Calculate new positions (rotated, centered back by default)
            conf = self.mol.GetConformer()
            new_positions = np.empty_like(positions)
            for i in range(self.mol.GetNumAtoms()):
                current_pos = np.array(conf.GetAtomPosition(i))
                centered_pos = current_pos - centroid
                if rotation_axis_norm > 1e-10:
                    rot_norm = rotation_axis_norm
                    rotation_axis_normalized = rotation_axis / rot_norm
                    cos_angle = np.dot(normal_vector, target_normal)
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    rotation_angle = np.arccos(cos_angle)
                    rotated_pos = rodrigues_rotate(
                        centered_pos,
                        rotation_axis_normalized,
                        rotation_angle,
                    )
                else:
                    rotated_pos = centered_pos
                new_pos = rotated_pos + centroid
                new_positions[i] = new_pos

            # If move_to_zero_plane is True, translate so the plane
            # of selected atoms is at zero
            if self.move_to_zero_plane_checkbox.isChecked():
                selected_new_positions = new_positions[selected_indices]
                new_centroid = np.mean(selected_new_positions, axis=0)
                translation_offset = np.zeros(3)
                if self.plane == "xy":
                    translation_offset[2] = new_centroid[2]
                elif self.plane == "xz":
                    translation_offset[1] = new_centroid[1]
                elif self.plane == "yz":
                    translation_offset[0] = new_centroid[0]
                new_positions = new_positions - translation_offset

            # Update the conformer positions array in place
            for i in range(self.mol.GetNumAtoms()):
                positions[i] = new_positions[i]

            self._update_molecule_geometry(positions)
            self.show_atom_labels()

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.exception("Failed to apply align: %s", e)
