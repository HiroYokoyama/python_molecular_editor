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
from typing import TYPE_CHECKING, Literal, Optional, Sequence

import numpy as np

from PyQt6.QtGui import QCloseEvent
from PyQt6.QtWidgets import (
    QCheckBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem, Geometry

from .dialog_3d_picking_mixin import Dialog3DPickingMixin

from .base_picking_dialog import SelectionList

from ..core.mol_geometry import rodrigues_rotate

if TYPE_CHECKING:
    from .main_window import MainWindow


class AlignmentDialog(Dialog3DPickingMixin, QDialog):
    """Dialog for aligning two selected atoms along a principal axis (X, Y, or Z)."""

    def __init__(
        self,
        mol: Chem.Mol,
        main_window: "MainWindow",
        axis: Literal["x", "y", "z"],
        preselected_atoms: Optional[Sequence[int]] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.axis = axis
        self._selected_atoms = SelectionList()

        # Add preselected atoms (maximum 2)
        if preselected_atoms:
            self._selected_atoms.update(preselected_atoms[:2])

        self.init_ui()

        # Add labels to preselected atoms
        if self._selected_atoms:
            for i, atom_idx in enumerate(self._selected_atoms, 1):
                self.add_selection_label(atom_idx, f"#{i}", color="yellow")
            self.update_display()

    @property
    def selected_atoms(self) -> SelectionList:
        """Return the ordered list of selected atom indices."""
        return self._selected_atoms

    @selected_atoms.setter
    def selected_atoms(self, val: object) -> None:
        """Replace the selection with a new SelectionList built from val."""
        self._selected_atoms = SelectionList(val)  # type: ignore[arg-type]

    def init_ui(self) -> None:
        """Build and lay out all widgets for the alignment dialog."""
        axis_names = {"x": "X-axis", "y": "Y-axis", "z": "Z-axis"}
        self.setWindowTitle(f"Align to {axis_names[self.axis]}")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            f"Click atoms in the 3D view to select them for alignment to the "
            f"{axis_names[self.axis]}. Exactly 2 atoms are required."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Move to origin option (default False)
        self.move_to_origin_checkbox = QCheckBox("Move the first atom to the origin")
        self.move_to_origin_checkbox.setChecked(False)
        layout.addWidget(self.move_to_origin_checkbox)

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

    def on_atom_picked(self, atom_idx: int) -> None:
        """Handle atom selection event."""
        if self.main_window.view_3d_manager.current_mol is None:
            return

        if atom_idx in self.selected_atoms:
            # De-select if already selected
            self.selected_atoms.remove(atom_idx)
            self.remove_atom_label(atom_idx)
        else:
            # Maximum of 2 atoms can be selected
            if len(self.selected_atoms) < 2:
                self.selected_atoms.append(atom_idx)
                # Show label indicating selection order
                label_text = f"#{len(self.selected_atoms)}"
                self.add_selection_label(atom_idx, label_text, color="yellow")

        self.update_display()

    def update_display(self) -> None:
        """Update the UI based on current selection."""
        if len(self.selected_atoms) == 0:
            self.selection_label.setText(
                "Click atoms to select for alignment (exactly 2 required)"
            )
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 1:
            atom = self.mol.GetAtomWithIdx(self.selected_atoms[0])
            self.selection_label.setText(f"Selected 1 atom: {atom.GetSymbol()}")
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 2:
            atom1 = self.mol.GetAtomWithIdx(self.selected_atoms[0])
            atom2 = self.mol.GetAtomWithIdx(self.selected_atoms[1])
            self.selection_label.setText(
                f"Selected 2 atoms: {atom1.GetSymbol()}, {atom2.GetSymbol()}"
            )
            self.apply_button.setEnabled(True)

    def clear_selection(self) -> None:
        """Clear the current selection and labels."""
        self.clear_selection_labels()
        self.selected_atoms.clear()
        self.update_display()

    def remove_atom_label(self, _atom_idx: int) -> None:
        """Remove a label for a specific atom (redraws all labels)."""
        self.clear_selection_labels()
        for i, idx in enumerate(self.selected_atoms, 1):
            self.add_selection_label(idx, f"#{i}", color="yellow")

    def apply_alignment(self) -> None:
        """Apply the specific axial alignment to the molecule."""
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(
                self, "Warning", "Please select exactly 2 atoms for alignment."
            )
            return
        try:
            atom1_idx = self.selected_atoms[0]
            atom2_idx = self.selected_atoms[1]

            conf = self.mol.GetConformer()

            # Get original atom positions
            positions = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(self.mol.GetNumAtoms())]
            )
            centroid = np.mean(positions, axis=0)

            pos1 = positions[atom1_idx]
            pos2 = positions[atom2_idx]

            # Calculate rotation to align atom1 -> atom2 relative to the chosen axis
            axis_vectors = {
                "x": np.array([1.0, 0.0, 0.0]),
                "y": np.array([0.0, 1.0, 0.0]),
                "z": np.array([0.0, 0.0, 1.0]),
            }
            target_axis = axis_vectors[self.axis]

            # Direction vector from atom1 to atom2
            current_vector = pos2 - pos1
            current_length = np.linalg.norm(current_vector)

            # Keep track of rotated positions (initially original positions)
            new_positions = np.copy(positions)

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

                    # Apply rotation to all atoms about the molecule's centroid
                    for i in range(self.mol.GetNumAtoms()):
                        rel_pos = positions[i] - centroid
                        rotated_pos = rodrigues_rotate(
                            rel_pos, rotation_axis, rotation_angle
                        )
                        new_positions[i] = rotated_pos + centroid

            # If move_to_origin is True, translate entire molecule so atom1 ends up at origin
            if (
                hasattr(self, "move_to_origin_checkbox")
                and self.move_to_origin_checkbox.isChecked()
            ):
                new_pos1 = new_positions[atom1_idx]
                new_positions = new_positions - new_pos1

            # Update conformer positions
            for i in range(self.mol.GetNumAtoms()):
                conf.SetAtomPosition(
                    i,
                    Geometry.Point3D(
                        float(new_positions[i][0]),
                        float(new_positions[i][1]),
                        float(new_positions[i][2]),
                    ),
                )

            # Update 3D positions
            self.main_window.view_3d_manager.atom_positions_3d = new_positions

            # Update 3D visualization
            self.main_window.view_3d_manager.draw_molecule_3d(self.mol)

            # Restore selection labels
            self.clear_selection_labels()
            for i, idx in enumerate(self.selected_atoms, 1):
                self.add_selection_label(idx, f"#{i}", color="yellow")

            # Update chirality labels
            self.main_window.view_3d_manager.update_chiral_labels()

            # Save state for Undo
            self.main_window.edit_actions_manager.push_undo_state()

            QMessageBox.information(
                self,
                "Success",
                f"Alignment to {self.axis.upper()}-axis completed.",
            )

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.exception("Failed to apply alignment")
            QMessageBox.critical(self, "Error", f"Failed to apply alignment: {str(e)}")

    def closeEvent(self, event: Optional[QCloseEvent]) -> None:
        """Clean up when the dialog is closed."""
        self.clear_selection_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self) -> None:
        """Handle dialog rejection/cancellation."""
        self.clear_selection_labels()
        self.disable_picking()
        super().reject()

    def accept(self) -> None:
        """Handle dialog acceptance (OK/Apply)."""
        self.clear_selection_labels()
        self.disable_picking()
        super().accept()
