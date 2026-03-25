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
    QDialog,
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

from .dialog_3d_picking_mixin import Dialog3DPickingMixin
from .mol_geometry import calc_distance, get_connected_group


class BondLengthDialog(Dialog3DPickingMixin, QDialog):  
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None

        # Set preselected atoms
        if preselected_atoms and len(preselected_atoms) >= 2:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]

        self.init_ui()

    def init_ui(self):
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
        self.distance_slider.setMinimum(10)   # 0.1 A
        self.distance_slider.setMaximum(1000) # 10.0 A
        self.distance_slider.setValue(154)    # 1.54 A
        self.distance_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.distance_slider.setTickInterval(100)
        self.distance_slider.setEnabled(False)
        self.distance_slider.sliderPressed.connect(self.on_slider_pressed)
        self.distance_slider.sliderMoved.connect(self.on_slider_moved)
        self.distance_slider.sliderReleased.connect(self.on_slider_released)
        self.distance_slider.valueChanged.connect(self.on_slider_value_changed)
        self._slider_dragging = False
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

    def on_atom_picked(self, atom_idx):
        """Handle atom picking event in the 3D view."""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None

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
        self.atom2_idx = None
        self.clear_selection_labels()
        self.update_display()

    def show_atom_labels(self):
        """Display labels on the selected atoms."""
        selected_atoms = [self.atom1_idx, self.atom2_idx]
        labels = ["1st", "2nd"]
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
                pass  # Suppress non-critical UI update errors (distance input/slider)

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
                pass  # Suppress non-critical UI update errors (distance input/slider)
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
            # Update the distance input box to show current distance
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
                pass  # Suppress errors during distance UI sync

            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")

    def on_distance_input_changed(self, text):
        """Line edit text changed, update slider."""
        if not self.distance_input.isEnabled() or not self.apply_button.isEnabled():
            return
        try:
            val = float(text)
            if 0.1 <= val <= 10.0:
                self.distance_slider.blockSignals(True)
                self.distance_slider.setValue(int(val * 100))
                self.distance_slider.blockSignals(False)
        except ValueError:
            pass  # Ignore invalid numeric input during typing

    def on_slider_pressed(self):
        """Remember the state before slider dragging starts."""
        if self.atom1_idx is None or self.atom2_idx is None:
            return
        self._slider_dragging = True
        self.main_window.push_undo_state()

    def on_slider_moved(self, value):
        """Update geometry in real-time while dragging."""
        if self.atom1_idx is None or self.atom2_idx is None:
            return
        
        new_distance = value / 100.0
        self.distance_input.blockSignals(True)
        self.distance_input.setText(f"{new_distance:.3f}")
        self.distance_input.blockSignals(False)
        
        self.adjust_bond_length(new_distance)

    def on_slider_released(self):
        """Finalize slider dragging."""
        self._slider_dragging = False
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()

    def on_slider_value_changed(self, value):
        """Handle click-to-position on the slider track."""
        if self._slider_dragging:
            return  # Already handled by on_slider_moved
        if self.atom1_idx is None or self.atom2_idx is None:
            return
        self.main_window.push_undo_state()
        new_distance = value / 100.0
        self.distance_input.blockSignals(True)
        self.distance_input.setText(f"{new_distance:.3f}")
        self.distance_input.blockSignals(False)
        self.adjust_bond_length(new_distance)
        self.main_window.update_chiral_labels()

    def apply_changes(self):
        """Apply the bond length changes to the molecule."""
        if self.atom1_idx is None or self.atom2_idx is None:
            return

        try:
            new_distance = float(self.distance_input.text())
            if new_distance <= 0:
                QMessageBox.warning(self, "Invalid Input", "Distance must be positive.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return

        # Save undo state
        self.main_window.push_undo_state()

        # Apply the bond length change
        self.adjust_bond_length(new_distance)

        # Update chirality labels
        self.main_window.update_chiral_labels()

    def adjust_bond_length(self, new_distance):
        """Adjust the bond length."""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))

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
            group1_atoms = get_connected_group(
                self.mol, self.atom1_idx, exclude=self.atom2_idx
            )
            group2_atoms = get_connected_group(
                self.mol, self.atom2_idx, exclude=self.atom1_idx
            )

            # Calculate displacements
            displacement1 = new_pos1 - pos1
            displacement2 = new_pos2 - pos2

            # Move group 1
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement1
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos

            # Move group 2
            for atom_idx in group2_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement2
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos

        elif self.atom1_fix_radio.isChecked():
            # Move only the second atom
            new_pos2 = pos1 + direction * new_distance
            conf.SetAtomPosition(self.atom2_idx, new_pos2.tolist())
            self.main_window.atom_positions_3d[self.atom2_idx] = new_pos2
        else:
            # Move the connected group (default behavior)
            new_pos2 = pos1 + direction * new_distance
            atoms_to_move = get_connected_group(
                self.mol, self.atom2_idx, exclude=self.atom1_idx
            )
            displacement = new_pos2 - pos2

            for atom_idx in atoms_to_move:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos

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
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass  # Suppress errors during dialog teardown
