#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING, Optional, Sequence

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

try:
    from .geometry_base_dialog import GeometryBaseDialog
    from ..core.mol_geometry import calc_distance, get_connected_group
except ImportError:
    from moleditpy_linux.ui.geometry_base_dialog import GeometryBaseDialog
    from moleditpy_linux.core.mol_geometry import calc_distance, get_connected_group

if TYPE_CHECKING:
    from .main_window import MainWindow


class BondLengthDialog(GeometryBaseDialog):
    def __init__(
        self,
        mol: Chem.Mol,
        main_window: "MainWindow",
        preselected_atoms: Optional[Sequence[int]] = None,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(mol, main_window, parent)
        self.atom1_idx: Optional[int] = None
        self.atom2_idx: Optional[int] = None

        # Set preselected atoms
        if preselected_atoms and len(preselected_atoms) >= 2:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]

        self.init_ui()

    def init_ui(self) -> None:
        self.setWindowTitle("Adjust Bond Length")
        self.setModal(False)

        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click two atoms in the 3D view to select a bond, then specify the new length."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Current distance display
        self.distance_label = QLabel("")
        layout.addWidget(self.distance_label)

        # New distance input
        distance_layout = QHBoxLayout()
        distance_layout.addWidget(QLabel("New distance (Å):"))
        self.distance_input = QLineEdit()
        self.distance_input.setPlaceholderText("1.54")
        self.distance_input.textChanged.connect(self.on_distance_input_changed)
        distance_layout.addWidget(self.distance_input)

        self.distance_slider = QSlider(Qt.Orientation.Horizontal)
        self.distance_slider.setMinimum(10)  # 0.1 A
        self.distance_slider.setMaximum(1000)  # 10.0 A
        self.distance_slider.setValue(154)  # 1.54 A
        self.distance_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.distance_slider.setTickInterval(100)
        self.distance_slider.setEnabled(False)

        # Connect to base class real-time handlers
        self.distance_slider.sliderPressed.connect(self.on_slider_pressed)
        self.distance_slider.sliderMoved.connect(
            lambda v: self.on_slider_moved_realtime(v, self.distance_input, 100.0)
        )
        self.distance_slider.sliderReleased.connect(self.on_slider_released)
        self.distance_slider.valueChanged.connect(
            lambda v: self.on_slider_value_changed_click(v, self.distance_input, 100.0)
        )

        layout.addLayout(distance_layout)
        layout.addWidget(self.distance_slider)

        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Movement Options:"))

        self.atom1_fix_group_radio = QRadioButton(
            "Atom 1: Fixed, Atom 2: Move connected group"
        )
        self.atom1_fix_group_radio.setChecked(True)
        group_layout.addWidget(self.atom1_fix_group_radio)

        self.atom1_fix_radio = QRadioButton("Atom 1: Fixed, Atom 2: Move atom only")
        group_layout.addWidget(self.atom1_fix_radio)

        self.both_groups_radio = QRadioButton(
            "Both groups: Move towards center equally"
        )
        group_layout.addWidget(self.both_groups_radio)

        layout.addWidget(group_box)

        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        button_layout.addStretch()

        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()

        # Update display if atoms are preselected
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()

    def on_atom_picked(self, atom_idx: int) -> None:
        """Handle atom picking event in the 3D view."""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None

        self.update_display()

    def clear_selection(self) -> None:
        """Clear the current atom selection."""
        self.atom1_idx = None
        self.atom2_idx = None
        self.clear_selection_labels()
        self.update_display()

    def show_atom_labels(self) -> None:
        """Display labels on the selected atoms."""
        selected_atoms = [self.atom1_idx, self.atom2_idx]
        labels = ["1st", "2nd"]
        pairs = [
            (idx, labels[i]) for i, idx in enumerate(selected_atoms) if idx is not None
        ]
        self.show_atom_labels_for(pairs)

    def update_display(self) -> None:
        """Update the UI display."""
        # Clear existing labels
        self.clear_selection_labels()

        if self.atom1_idx is None:
            self.selection_label.setText("No atoms selected")
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear distance input when no selection
            try:
                self.distance_input.blockSignals(True)
                self.distance_input.clear()
                self.distance_input.blockSignals(False)
                self.distance_slider.blockSignals(True)
                self.distance_slider.setValue(154)
                self.distance_slider.setEnabled(False)
                self.distance_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass

        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(
                f"First atom: {symbol1} (index {self.atom1_idx})"
            )
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # Add label
            self.add_selection_label(self.atom1_idx, "1")
            # Clear distance input while selection is incomplete
            try:
                self.distance_input.blockSignals(True)
                self.distance_input.clear()
                self.distance_input.blockSignals(False)
                self.distance_slider.blockSignals(True)
                self.distance_slider.setValue(154)
                self.distance_slider.setEnabled(False)
                self.distance_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(
                f"Bond: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx})"
            )

            # Calculate current distance
            conf = self.mol.GetConformer()
            current_distance = calc_distance(
                conf.GetAtomPosition(self.atom1_idx),
                conf.GetAtomPosition(self.atom2_idx),
            )
            self.distance_label.setText(f"Current distance: {current_distance:.3f} Å")
            self.apply_button.setEnabled(True)
            # Update the distance input box and slider
            try:
                self.distance_input.blockSignals(True)
                self.distance_input.setText(f"{current_distance:.3f}")
                self.distance_input.blockSignals(False)
                self.distance_slider.blockSignals(True)
                slider_val = int(current_distance * 100)
                slider_val = max(10, min(1000, slider_val))
                self.distance_slider.setValue(slider_val)
                self.distance_slider.setEnabled(True)
                self.distance_slider.blockSignals(False)
            except (AttributeError, RuntimeError, TypeError):
                pass

            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")

    def on_distance_input_changed(self, text: str) -> None:
        """Line edit text changed, update slider."""
        if not self.distance_input.isEnabled() or not self.apply_button.isEnabled():
            return
        self._sync_input_to_slider(text, self.distance_slider, 100.0)

    def apply_changes(self) -> None:
        """Apply the bond length changes to the molecule."""
        if not self._is_selection_complete():
            return

        try:
            new_distance = float(self.distance_input.text())
            if new_distance <= 0:
                QMessageBox.warning(self, "Invalid Input", "Distance must be positive.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return

        # Apply the update
        self.apply_geometry_update(float(new_distance))

        # Push Undo state AFTER modification
        self._push_undo()

    def _is_selection_complete(self) -> bool:
        return self.atom1_idx is not None and self.atom2_idx is not None

    def apply_geometry_update(self, new_distance: float) -> None:  # pylint: disable=arguments-renamed
        """Adjust the bond length."""
        conf = self.mol.GetConformer()

        # Use snapshot if available (slider dragging) to keep base positions stable and accurate
        # Reverting to dictionary per user request for clarity/consistency
        snapshot = self._snapshot_positions
        if snapshot is not None:
            # If snapshot is a numpy array (from base class), convert to dict or use as is
            if isinstance(snapshot, dict):
                positions = snapshot.copy()
            else:
                positions = {i: snapshot[i].copy() for i in range(len(snapshot))}
        else:
            positions = {
                i: np.array(conf.GetAtomPosition(i))
                for i in range(self.mol.GetNumAtoms())
            }

        if self.atom1_idx is None or self.atom2_idx is None:
            return
        idx1: int = self.atom1_idx
        idx2: int = self.atom2_idx
        pos1 = positions[idx1]
        pos2 = positions[idx2]

        # Direction vector from atom1 to atom2
        direction = pos2 - pos1
        current_distance = calc_distance(pos1, pos2)

        if current_distance == 0:
            return

        direction = direction / current_distance

        if self.both_groups_radio.isChecked():
            # Both groups move towards center equally
            bond_center = (pos1 + pos2) / 2
            half_distance = new_distance / 2

            # New positions for both atoms
            new_pos1 = bond_center - direction * half_distance
            new_pos2 = bond_center + direction * half_distance

            # Get both connected groups
            group1_atoms = get_connected_group(self.mol, idx1, exclude=idx2)
            group2_atoms = get_connected_group(self.mol, idx2, exclude=idx1)

            # Calculate displacements
            displacement1 = new_pos1 - pos1
            displacement2 = new_pos2 - pos2

            # Move group 1
            for atom_idx in group1_atoms:
                positions[atom_idx] += displacement1

            # Move group 2
            for atom_idx in group2_atoms:
                positions[atom_idx] += displacement2
        elif self.atom1_fix_radio.isChecked():
            # Move only the second atom
            new_pos2 = pos1 + direction * new_distance
            positions[idx2] = new_pos2
        else:
            # Move the connected group (default behavior)
            new_pos2 = pos1 + direction * new_distance
            atoms_to_move = get_connected_group(self.mol, idx2, exclude=idx1)
            displacement = new_pos2 - pos2

            for atom_idx in atoms_to_move:
                positions[atom_idx] += displacement

        # Write updated positions back using inherited helper
        self._update_molecule_geometry(positions)
