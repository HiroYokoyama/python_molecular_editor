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
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
    from .mol_geometry import adjust_bond_angle, calc_angle_deg, get_connected_group, rodrigues_rotate
except ImportError:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin
    from modules.mol_geometry import adjust_bond_angle, calc_angle_deg, get_connected_group, rodrigues_rotate

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMessageBox


class AngleDialog(Dialog3DPickingMixin, QDialog):  # pragma: no cover
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
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
        self.angle_slider.sliderPressed.connect(self.on_slider_pressed)
        self.angle_slider.sliderMoved.connect(self.on_slider_moved)
        self.angle_slider.sliderReleased.connect(self.on_slider_released)
        self.angle_slider.valueChanged.connect(self.on_slider_value_changed)
        self._slider_dragging = False
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

        # Connect to main window's picker for AngleDialog
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
            # This locks in the "original" initial geometry the user started modifying from
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self._snapshot_positions = None

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
        """Handle OK event."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()

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
            except (AttributeError, RuntimeError):
                import traceback
                traceback.print_exc()

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
            except (AttributeError, RuntimeError):
                import traceback
                traceback.print_exc()

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
            except (AttributeError, RuntimeError):
                import traceback
                traceback.print_exc()
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
            except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
                import traceback
                traceback.print_exc()

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
        try:
            val = float(text)
            wrapped_val = (val + 180) % 360 - 180
            self.angle_slider.blockSignals(True)
            self.angle_slider.setValue(int(round(wrapped_val)))
            self.angle_slider.blockSignals(False)
        except ValueError:
            import traceback
            traceback.print_exc()

    def on_slider_pressed(self):
        """Remember the state before slider dragging starts."""
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return
        self._slider_dragging = True
        self.main_window.push_undo_state()
        # Snapshot positions so the rotation axis stays stable during drag
        # Only take snapshot if one doesn't exist to preserve directional info
        if getattr(self, '_snapshot_positions', None) is None:
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()

    def on_slider_moved(self, value):
        """Update geometry in real-time while dragging."""
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return
        
        self.angle_input.blockSignals(True)
        self.angle_input.setText(f"{value}")
        self.angle_input.blockSignals(False)
        
        self.adjust_angle(float(value))

    def on_slider_released(self):
        """Finalize slider dragging."""
        self._slider_dragging = False
        # Do NOT clear snapshot here. Keep it to preserve the turning direction
        # even if they approach ±180°. It gets cleared when selection is changed.
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()

    def on_slider_value_changed(self, value):
        """Handle click-to-position on the slider track."""
        if self._slider_dragging:
            return
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return
        self.main_window.push_undo_state()
        
        # Ensure we have a snapshot for click-to-position as well to maintain direction
        if getattr(self, '_snapshot_positions', None) is None:
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()
            
        self.angle_input.blockSignals(True)
        self.angle_input.setText(f"{value}")
        self.angle_input.blockSignals(False)
        self.adjust_angle(float(value))
        self.main_window.update_chiral_labels()

    def apply_changes(self):
        """Apply the angle changes to the molecule."""
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return

        try:
            raw_angle = float(self.angle_input.text())
            # Automatic Range Wrapping
            new_angle = (raw_angle + 180) % 360 - 180
            
            # Formally update the input to reflect wrapping
            self.angle_input.blockSignals(True)
            self.angle_input.setText(f"{new_angle:.2f}")
            self.angle_input.blockSignals(False)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return

        # Save undo state
        self.main_window.push_undo_state()

        # Apply the angle change
        self.adjust_angle(new_angle)

        # Update chirality labels
        self.main_window.update_chiral_labels()

    def adjust_angle(self, new_angle_deg):
        """Adjust the bond angle (with options for group rotation).

        Uses the difference-based rotation approach via
        :func:`~mol_geometry.adjust_bond_angle` to avoid 3D
        rotational ambiguity.

        During slider dragging, positions are restored from a snapshot
        taken at press-time so that the rotation axis (cross product)
        never flips direction.
        """
        conf = self.mol.GetConformer()

        # Use snapshot if available (slider dragging) to keep the
        # rotation axis stable; otherwise use current positions.
        snapshot = getattr(self, '_snapshot_positions', None)
        if snapshot is not None:
            positions = snapshot.copy()
        else:
            positions = conf.GetPositions()  # N×3 ndarray (copy)

        idx_a = self.atom1_idx
        idx_b = self.atom2_idx  # vertex
        idx_c = self.atom3_idx

        if self.both_groups_radio.isChecked():
            # Both arms rotate equally (half angle each)
            current_angle = self.calculate_angle()
            half_delta_deg = (new_angle_deg - current_angle) / 2.0

            group1 = get_connected_group(self.mol, idx_a, exclude=idx_b)
            group3 = get_connected_group(self.mol, idx_c, exclude=idx_b)

            # Arm 1 rotates by −half (note: reversed A/C roles)
            adjust_bond_angle(
                positions, idx_c, idx_b, idx_a,
                current_angle + half_delta_deg, group1,
            )
            # Arm 3 rotates by +half
            adjust_bond_angle(
                positions, idx_a, idx_b, idx_c,
                current_angle + half_delta_deg, group3,
            )
        elif self.rotate_atom_radio.isChecked():
            # Move only atom C
            adjust_bond_angle(
                positions, idx_a, idx_b, idx_c,
                new_angle_deg, {idx_c},
            )
        else:
            # Default: rotate atom C and its connected sub-structure
            atoms_to_move = get_connected_group(
                self.mol, idx_c, exclude=idx_b,
            )
            adjust_bond_angle(
                positions, idx_a, idx_b, idx_c,
                new_angle_deg, atoms_to_move,
            )

        # Write updated positions back to the conformer and 3D cache
        for i in range(conf.GetNumAtoms()):
            conf.SetAtomPosition(i, positions[i].tolist())
            self.main_window.atom_positions_3d[i] = positions[i]

        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)

    def reject(self):
        """Handle cancellation event."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
        try:
            if self.main_window.current_mol:
                self.main_window.draw_molecule_3d(self.main_window.current_mol)
        except (AttributeError, RuntimeError):
            import traceback
            traceback.print_exc()
