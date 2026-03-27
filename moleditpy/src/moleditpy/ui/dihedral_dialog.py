import logging
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
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

try:
    from .dialog_3d_picking_mixin import Dialog3DPickingMixin
    from ..core.mol_geometry import calculate_dihedral, get_connected_group
except ImportError:
    from moleditpy.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin
    from moleditpy.core.mol_geometry import calculate_dihedral, get_connected_group

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMessageBox
from rdkit import Geometry


class DihedralDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None

        # Set pre-selected atoms
        if preselected_atoms and len(preselected_atoms) >= 4:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # central bond start
            self.atom3_idx = preselected_atoms[2]  # central bond end
            self.atom4_idx = preselected_atoms[3]

        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Adjust Dihedral Angle")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click four atoms in order to define a dihedral angle. The rotation will be around the bond between the 2nd and 3rd atoms."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Current dihedral angle display
        self.dihedral_label = QLabel("")
        layout.addWidget(self.dihedral_label)

        # New dihedral angle input
        dihedral_layout = QHBoxLayout()
        dihedral_layout.addWidget(QLabel("New dihedral angle (degrees):"))
        self.dihedral_input = QLineEdit()
        self.dihedral_input.setPlaceholderText("180.0")
        self.dihedral_input.textChanged.connect(self.on_dihedral_input_changed)
        dihedral_layout.addWidget(self.dihedral_input)

        self.dihedral_slider = QSlider(Qt.Orientation.Horizontal)
        self.dihedral_slider.setMinimum(-180)
        self.dihedral_slider.setMaximum(180)
        self.dihedral_slider.setValue(180)
        self.dihedral_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.dihedral_slider.setTickInterval(45)
        self.dihedral_slider.setEnabled(False)
        self.dihedral_slider.sliderPressed.connect(self.on_slider_pressed)
        self.dihedral_slider.sliderMoved.connect(self.on_slider_moved)
        self.dihedral_slider.sliderReleased.connect(self.on_slider_released)
        self.dihedral_slider.valueChanged.connect(self.on_slider_value_changed)
        self._slider_dragging = False
        layout.addLayout(dihedral_layout)
        layout.addWidget(self.dihedral_slider)

        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Move:"))

        self.move_group_radio = QRadioButton("Atom 1,2,3: Fixed, Atom 4 group: Rotate")
        self.move_group_radio.setChecked(True)
        group_layout.addWidget(self.move_group_radio)

        self.move_atom_radio = QRadioButton(
            "Atom 1,2,3: Fixed, Atom 4: Rotate atom only"
        )
        group_layout.addWidget(self.move_atom_radio)

        self.both_groups_radio = QRadioButton(
            "Central bond fixed: Both groups rotate equally"
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

        # Connect to main window's picker for DihedralDialog
        self.picker_connection = None
        self.enable_picking()

        # Update initial display if atoms are pre-selected
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()

    def on_atom_picked(self, atom_idx):
        """Handle atom picked event."""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
        elif self.atom4_idx is None:
            self.atom4_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self.atom4_idx = None

        # Display atom labels
        self.show_atom_labels()
        self.update_display()

    def keyPressEvent(self, event):
        """Handle keyboard events."""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_changes()
            event.accept()
        else:
            super().keyPressEvent(event)

    def closeEvent(self, event):
        """Handle dialog close event."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def accept(self):
        """Handle OK action."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()

    def clear_selection(self):
        """Clear selection."""
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None
        self.clear_atom_labels()
        self.update_display()

    def show_atom_labels(self):
        """Display labels on selected atoms."""
        selected_atoms = [
            self.atom1_idx,
            self.atom2_idx,
            self.atom3_idx,
            self.atom4_idx,
        ]
        labels = ["1st", "2nd (bond start)", "3rd (bond end)", "4th"]
        pairs = [
            (idx, labels[i]) for i, idx in enumerate(selected_atoms) if idx is not None
        ]
        self.show_atom_labels_for(pairs)

    def update_display(self):
        """Update display."""
        selected_count = sum(
            x is not None
            for x in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        )

        if selected_count == 0:
            self.selection_label.setText("No atoms selected")
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input when no selection
            try:
                self.dihedral_input.blockSignals(True)
                self.dihedral_input.clear()
                self.dihedral_input.blockSignals(False)
                self.dihedral_slider.blockSignals(True)
                self.dihedral_slider.setValue(180)
                self.dihedral_slider.setEnabled(False)
                self.dihedral_slider.blockSignals(False)
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress errors during dihedral input clearing

        elif selected_count < 4:
            selected_atoms = [
                self.atom1_idx,
                self.atom2_idx,
                self.atom3_idx,
                self.atom4_idx,
            ]

            display_parts = []
            for atom_idx in selected_atoms:
                if atom_idx is not None:
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    display_parts.append(f"{symbol}({atom_idx})")
                else:
                    display_parts.append("?")

            self.selection_label.setText(" - ".join(display_parts))
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input while selection is incomplete
            try:
                self.dihedral_input.blockSignals(True)
                self.dihedral_input.clear()
                self.dihedral_input.blockSignals(False)
                self.dihedral_slider.blockSignals(True)
                self.dihedral_slider.setValue(180)
                self.dihedral_slider.setEnabled(False)
                self.dihedral_slider.blockSignals(False)
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress non-critical UI update errors
        else:
            selected_atoms = [
                self.atom1_idx,
                self.atom2_idx,
                self.atom3_idx,
                self.atom4_idx,
            ]

            display_parts = []
            for atom_idx in selected_atoms:
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                display_parts.append(f"{symbol}({atom_idx})")

            self.selection_label.setText(" - ".join(display_parts))

            # Calculate current dihedral angle
            current_dihedral = calculate_dihedral(
                self.mol.GetConformer().GetPositions(),
                self.atom1_idx,
                self.atom2_idx,
                self.atom3_idx,
                self.atom4_idx,
            )
            self.dihedral_label.setText(f"Current dihedral: {current_dihedral:.2f}°")
            self.apply_button.setEnabled(True)
            # Update dihedral input box with current dihedral
            try:
                self.dihedral_input.blockSignals(True)
                self.dihedral_input.setText(f"{current_dihedral:.2f}")
                self.dihedral_input.blockSignals(False)
                self.dihedral_slider.blockSignals(True)
                slider_val = int(round(current_dihedral))
                slider_val = max(-180, min(180, slider_val))
                self.dihedral_slider.setValue(slider_val)
                self.dihedral_slider.setEnabled(True)
                self.dihedral_slider.blockSignals(False)
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress non-critical UI update errors

    def on_dihedral_input_changed(self, text):
        """Line edit text changed, update slider."""
        if not self.dihedral_input.isEnabled() or not self.apply_button.isEnabled():
            return
        try:
            val = float(text)
            wrapped_val = (val + 180) % 360 - 180
            self.dihedral_slider.blockSignals(True)
            self.dihedral_slider.setValue(int(round(wrapped_val)))
            self.dihedral_slider.blockSignals(False)
        except ValueError as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Ignore invalid numeric input during typing

    def on_slider_pressed(self):
        """Remember the state before slider dragging starts."""
        if any(
            idx is None
            for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        ):
            return
        self._slider_dragging = True
        self.main_window.push_undo_state()

    def on_slider_moved(self, value):
        """Update geometry in real-time while dragging."""
        if any(
            idx is None
            for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        ):
            return

        self.dihedral_input.blockSignals(True)
        self.dihedral_input.setText(f"{value}")
        self.dihedral_input.blockSignals(False)

        self.adjust_dihedral(float(value))

    def on_slider_released(self):
        """Finalize slider dragging."""
        self._slider_dragging = False
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()

    def on_slider_value_changed(self, value):
        """Handle click-to-position on the slider track."""
        if self._slider_dragging:
            return
        if any(
            idx is None
            for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        ):
            return
        self.main_window.push_undo_state()
        self.dihedral_input.blockSignals(True)
        self.dihedral_input.setText(f"{value}")
        self.dihedral_input.blockSignals(False)
        self.adjust_dihedral(float(value))
        self.main_window.update_chiral_labels()

    def apply_changes(self):
        """Apply changes."""
        if any(
            idx is None
            for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        ):
            return

        try:
            raw_dihedral = float(self.dihedral_input.text())
            # Automatic Range Wrapping
            new_dihedral = (raw_dihedral + 180) % 360 - 180

            # Formally update the input to reflect wrapping
            self.dihedral_input.blockSignals(True)
            self.dihedral_input.setText(f"{new_dihedral:.2f}")
            self.dihedral_input.blockSignals(False)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return

        # Apply the dihedral angle change
        self.adjust_dihedral(new_dihedral)

        # Update chiral labels
        self.main_window.update_chiral_labels()

        # Save undo state
        self.main_window.push_undo_state()

    def adjust_dihedral(self, new_dihedral_deg):
        """Adjust dihedral angle (improved algorithm)."""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        pos4 = np.array(conf.GetAtomPosition(self.atom4_idx))

        # Current dihedral angle
        current_dihedral = calculate_dihedral(
            self.mol.GetConformer().GetPositions(),
            self.atom1_idx,
            self.atom2_idx,
            self.atom3_idx,
            self.atom4_idx,
        )

        # Calculate rotation angle needed
        rotation_angle_deg = new_dihedral_deg - current_dihedral

        # Handle angle wrapping for shortest rotation
        if rotation_angle_deg > 180:
            rotation_angle_deg -= 360
        elif rotation_angle_deg < -180:
            rotation_angle_deg += 360

        rotation_angle_rad = np.radians(rotation_angle_deg)

        # Skip if no rotation needed
        if abs(rotation_angle_rad) < 1e-6:
            return

        # Rotation axis is the bond between atom2 and atom3
        rotation_axis = pos3 - pos2
        axis_length = np.linalg.norm(rotation_axis)

        if axis_length == 0:
            return  # Atoms are at the same position

        rotation_axis = rotation_axis / axis_length

        # Rodrigues' rotation formula implementation
        def rotate_point_around_axis(point, axis_point, axis_direction, angle):
            """Rotate a point around an axis using Rodrigues' formula"""
            # Translate point so axis passes through origin
            translated_point = point - axis_point

            # Apply Rodrigues' rotation formula
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)

            rotated = (
                translated_point * cos_a
                + np.cross(axis_direction, translated_point) * sin_a
                + axis_direction
                * np.dot(axis_direction, translated_point)
                * (1 - cos_a)
            )

            # Translate back to original coordinate system
            return rotated + axis_point

        if self.both_groups_radio.isChecked():
            # Both groups rotate equally around the central bond (half angle each in opposite directions)
            half_rotation = rotation_angle_rad / 2

            # Get both connected groups
            group1_atoms = get_connected_group(
                self.mol, self.atom2_idx, exclude=self.atom3_idx
            )
            group4_atoms = get_connected_group(
                self.mol, self.atom3_idx, exclude=self.atom2_idx
            )

            # Rotate group1 (atom1 side) by -half_rotation
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(
                    current_pos, pos2, rotation_axis, -half_rotation
                )
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos

            # Rotate group4 (atom4 side) by +half_rotation
            for atom_idx in group4_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(
                    current_pos, pos2, rotation_axis, half_rotation
                )
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos

        elif self.move_group_radio.isChecked():
            # Move the connected group containing atom4
            # Find all atoms connected to atom3 (excluding atom2 side)
            atoms_to_rotate = get_connected_group(
                self.mol, self.atom3_idx, exclude=self.atom2_idx
            )

            # Rotate all atoms in the group
            for atom_idx in atoms_to_rotate:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(
                    current_pos, pos2, rotation_axis, rotation_angle_rad
                )
                conf.SetAtomPosition(
                    atom_idx,
                    Geometry.Point3D(
                        float(new_pos[0]), float(new_pos[1]), float(new_pos[2])
                    ),
                )
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        else:
            # Move only atom4
            new_pos4 = rotate_point_around_axis(
                pos4, pos2, rotation_axis, rotation_angle_rad
            )
            conf.SetAtomPosition(
                self.atom4_idx,
                Geometry.Point3D(
                    float(new_pos4[0]), float(new_pos4[1]), float(new_pos4[2])
                ),
            )
            self.main_window.atom_positions_3d[self.atom4_idx] = new_pos4

        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)

    def reject(self):
        """Handle cancel action."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
        try:
            if self.main_window.current_mol:
                self.main_window.draw_molecule_3d(self.main_window.current_mol)
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress errors during dialog teardown
