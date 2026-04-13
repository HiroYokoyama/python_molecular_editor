#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSlider,
    QVBoxLayout,
    QWidget,
    QMessageBox,
)
from PyQt6.QtCore import Qt

try:
    from .geometry_base_dialog import GeometryBaseDialog
    from ..core.mol_geometry import (
        adjust_bond_angle,
        calc_angle_deg,
        get_connected_group,
    )
except ImportError:
    from moleditpy.ui.geometry_base_dialog import GeometryBaseDialog
    from moleditpy.core.mol_geometry import (
        adjust_bond_angle,
        calc_angle_deg,
        get_connected_group,
    )


class AngleDialog(GeometryBaseDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        super().__init__(mol, main_window, parent)
        self.atom1_idx = None
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None

        # Set preselected atoms
        if preselected_atoms and len(preselected_atoms) >= 3:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # vertex
            self.atom3_idx = preselected_atoms[2]

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Adjust Angle")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click three atoms in order: first-vertex-third. The angle around the vertex atom will be adjusted."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Current angle display
        self.angle_label = QLabel("")
        layout.addWidget(self.angle_label)

        # New angle input
        angle_layout = QHBoxLayout()
        angle_layout.addWidget(QLabel("New angle (degrees):"))
        self.angle_input = QLineEdit()
        self.angle_input.setPlaceholderText("109.5")
        self.angle_input.textChanged.connect(self.on_angle_input_changed)
        angle_layout.addWidget(self.angle_input)

        self.angle_slider = QSlider(Qt.Orientation.Horizontal)
        self.angle_slider.setMinimum(-180)
        self.angle_slider.setMaximum(180)
        self.angle_slider.setValue(109)
        self.angle_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.angle_slider.setTickInterval(45)
        self.angle_slider.setEnabled(False)

        # Connect to base class real-time handlers
        self.angle_slider.sliderPressed.connect(self.on_slider_pressed)
        self.angle_slider.sliderMoved.connect(
            lambda v: self.on_slider_moved_realtime(v, self.angle_input, 1.0)
        )
        self.angle_slider.sliderReleased.connect(self.on_slider_released)
        self.angle_slider.valueChanged.connect(
            lambda v: self.on_slider_value_changed_click(v, self.angle_input, 1.0)
        )

        layout.addLayout(angle_layout)
        layout.addWidget(self.angle_slider)

        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Rotation Options:"))

        self.rotate_group_radio = QRadioButton(
            "Atom 1,2: Fixed, Atom 3: Rotate connected group"
        )
        self.rotate_group_radio.setChecked(True)
        group_layout.addWidget(self.rotate_group_radio)

        self.rotate_atom_radio = QRadioButton(
            "Atom 1,2: Fixed, Atom 3: Rotate atom only"
        )
        group_layout.addWidget(self.rotate_atom_radio)

        self.both_groups_radio = QRadioButton("Vertex fixed: Both arms rotate equally")
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

    def on_atom_picked(self, atom_idx):
        """Handle atom picking event in the 3D view."""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
            # Take a fresh snapshot immediately upon completing the triad selection
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self._snapshot_positions = None

        self.update_display()

    def clear_selection(self):
        """Clear the current atom selection."""
        self.atom1_idx = None
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None
        self._snapshot_positions = None
        self.clear_selection_labels()
        self.update_display()

    def show_atom_labels(self):
        """Display labels on the selected atoms."""
        selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx]
        labels = ["1st", "2nd (vertex)", "3rd"]
        pairs = [
            (idx, labels[i]) for i, idx in enumerate(selected_atoms) if idx is not None
        ]
        self.show_atom_labels_for(pairs)

    def update_display(self):
        """Update the UI display."""
        # Clear existing labels
        self.clear_selection_labels()

        if self.atom1_idx is None:
            self.selection_label.setText("No atoms selected")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            self._snapshot_positions = None
            # Clear angle input when no selection
            try:
                self.angle_input.blockSignals(True)
                self.angle_input.clear()
                self.angle_input.blockSignals(False)
                self.angle_slider.blockSignals(True)
                self.angle_slider.setValue(109)
                self.angle_slider.setEnabled(False)
                self.angle_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(
                f"First atom: {symbol1} (index {self.atom1_idx})"
            )
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            self._snapshot_positions = None
            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.blockSignals(True)
                self.angle_input.clear()
                self.angle_input.blockSignals(False)
                self.angle_slider.blockSignals(True)
                self.angle_slider.setValue(109)
                self.angle_slider.setEnabled(False)
                self.angle_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
        elif self.atom3_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(
                f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - ?"
            )
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            self._snapshot_positions = None
            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.blockSignals(True)
                self.angle_input.clear()
                self.angle_input.blockSignals(False)
                self.angle_slider.blockSignals(True)
                self.angle_slider.setValue(109)
                self.angle_slider.setEnabled(False)
                self.angle_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            self.selection_label.setText(
                f"Angle: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - {symbol3}({self.atom3_idx})"
            )

            # Calculate current angle
            current_angle = self.calculate_angle()
            self.angle_label.setText(f"Current angle: {current_angle:.2f}°")
            self.apply_button.setEnabled(True)
            # Update angle input box with current angle
            try:
                self.angle_input.blockSignals(True)
                self.angle_input.setText(f"{current_angle:.2f}")
                self.angle_input.blockSignals(False)
                self.angle_slider.blockSignals(True)
                slider_val = int(round(current_angle))
                slider_val = max(-180, min(180, slider_val))
                self.angle_slider.setValue(slider_val)
                self.angle_slider.setEnabled(True)
                self.angle_slider.blockSignals(False)
            except (AttributeError, RuntimeError, TypeError):
                pass

            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            self.add_selection_label(self.atom3_idx, "3")

    def calculate_angle(self):
        """Calculate the current bond angle."""
        conf = self.mol.GetConformer()
        pos1 = conf.GetAtomPosition(self.atom1_idx)
        pos2 = conf.GetAtomPosition(self.atom2_idx)  # vertex
        pos3 = conf.GetAtomPosition(self.atom3_idx)
        return calc_angle_deg(pos1, pos2, pos3)

    def on_angle_input_changed(self, text):
        """Line edit text changed, update slider."""
        if not self.angle_input.isEnabled() or not self.apply_button.isEnabled():
            return
        self._sync_input_to_slider(text, self.angle_slider, 1.0, wrap=True)

    def apply_changes(self):
        """Apply the angle changes to the molecule."""
        if not self._is_selection_complete():
            return

        try:
            raw_angle = float(self.angle_input.text())
            # Automatic Range Wrapping
            new_angle = (raw_angle + 180) % 360 - 180

            # Formally update the input to reflect wrapping
            self.angle_input.blockSignals(True)
            self.angle_input.setText(f"{new_angle:.2f}")
            self.angle_input.blockSignals(False)
            self._snapshot_positions = None
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return

        # Apply the update
        self.apply_geometry_update(new_angle)

        # Push Undo state AFTER modification
        self._push_undo()

    def _is_selection_complete(self):
        return (
            self.atom1_idx is not None
            and self.atom2_idx is not None
            and self.atom3_idx is not None
        )

    def apply_geometry_update(self, new_angle_deg):
        """Adjust the bond angle."""
        conf = self.mol.GetConformer()

        # Use snapshot if available (slider dragging) to keep the rotation axis stable
        snapshot = self._snapshot_positions
        if snapshot is not None:
            positions = snapshot.copy()
        else:
            positions = conf.GetPositions()

        idx_a = self.atom1_idx
        idx_b = self.atom2_idx  # vertex
        idx_c = self.atom3_idx

        if self.both_groups_radio.isChecked():
            # Both arms rotate equally (half angle each)
            current_angle = self.calculate_angle()
            half_delta_deg = (new_angle_deg - current_angle) / 2.0

            group1 = get_connected_group(self.mol, idx_a, exclude=idx_b)
            group3 = get_connected_group(self.mol, idx_c, exclude=idx_b)

            # Arm 1 rotates by -half relative to C-B-A
            adjust_bond_angle(
                positions,
                idx_c,
                idx_b,
                idx_a,
                current_angle + half_delta_deg,
                group1,
            )
            # Arm 3 rotates to the FINAL angle (relative to the now-moved Arm 1)
            adjust_bond_angle(
                positions,
                idx_a,
                idx_b,
                idx_c,
                new_angle_deg,
                group3,
            )
        elif self.rotate_atom_radio.isChecked():
            # Move only atom C
            adjust_bond_angle(
                positions,
                idx_a,
                idx_b,
                idx_c,
                new_angle_deg,
                {idx_c},
            )
        else:
            # Default: rotate atom C and its connected sub-structure
            atoms_to_move = get_connected_group(
                self.mol,
                idx_c,
                exclude=idx_b,
            )
            adjust_bond_angle(
                positions,
                idx_a,
                idx_b,
                idx_c,
                new_angle_deg,
                atoms_to_move,
            )

        # Write updated positions back using inherited helper
        self._update_molecule_geometry(positions)
