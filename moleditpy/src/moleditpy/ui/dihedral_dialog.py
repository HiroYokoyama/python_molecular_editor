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
    from ..core.mol_geometry import (
        adjust_dihedral,
        calculate_dihedral,
        get_connected_group,
    )
except ImportError:
    from moleditpy.ui.geometry_base_dialog import GeometryBaseDialog
    from moleditpy.core.mol_geometry import (
        adjust_dihedral,
        calculate_dihedral,
        get_connected_group,
    )

if TYPE_CHECKING:
    from .main_window import MainWindow


class DihedralDialog(GeometryBaseDialog):
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
        self.atom3_idx: Optional[int] = None
        self.atom4_idx: Optional[int] = None

        # Set preselected atoms
        if preselected_atoms and len(preselected_atoms) >= 4:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]
            self.atom3_idx = preselected_atoms[2]
            self.atom4_idx = preselected_atoms[3]

        self.init_ui()

    def init_ui(self) -> None:
        self.setWindowTitle("Adjust Dihedral Angle")
        self.setModal(False)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click four atoms in order (1-2-3-4). The dihedral angle around the 2-3 bond will be adjusted."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Current dihedral display
        self.dihedral_label = QLabel("")
        layout.addWidget(self.dihedral_label)

        # New dihedral input
        dihedral_layout = QHBoxLayout()
        dihedral_layout.addWidget(QLabel("New dihedral (degrees):"))
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

        # Connect to base class real-time handlers
        self.dihedral_slider.sliderPressed.connect(self.on_slider_pressed)
        self.dihedral_slider.sliderMoved.connect(
            lambda v: self.on_slider_moved_realtime(v, self.dihedral_input, 1.0)
        )
        self.dihedral_slider.sliderReleased.connect(self.on_slider_released)
        self.dihedral_slider.valueChanged.connect(
            lambda v: self.on_slider_value_changed_click(v, self.dihedral_input, 1.0)
        )

        layout.addLayout(dihedral_layout)
        layout.addWidget(self.dihedral_slider)

        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Rotation Options:"))

        self.rotate_group_radio = QRadioButton(
            "Atoms 1,2,3: Fixed, Atom 4: Rotate connected group"
        )
        self.rotate_group_radio.setChecked(True)
        group_layout.addWidget(self.rotate_group_radio)

        self.rotate_atom_radio = QRadioButton(
            "Atoms 1,2,3: Fixed, Atom 4: Rotate atom only"
        )
        group_layout.addWidget(self.rotate_atom_radio)

        self.both_groups_radio = QRadioButton(
            "Bond (2-3) fixed: Both ends rotate equally"
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
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
        elif self.atom4_idx is None:
            self.atom4_idx = atom_idx
            # Take a fresh snapshot immediately upon completing the selection
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self.atom4_idx = None
            self._snapshot_positions = None

        self.update_display()

    def clear_selection(self) -> None:
        """Clear the current atom selection."""
        self.atom1_idx = None
        self.atom2_idx = None
        self.atom3_idx = None
        self.atom4_idx = None
        self._snapshot_positions = None
        self.clear_selection_labels()
        self.update_display()

    def show_atom_labels(self) -> None:
        """Display labels on the selected atoms."""
        selected_atoms = [
            self.atom1_idx,
            self.atom2_idx,
            self.atom3_idx,
            self.atom4_idx,
        ]
        labels = ["1st", "2nd", "3rd", "4th"]
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
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            self._snapshot_positions = None
            # Clear input
            try:
                self.dihedral_input.blockSignals(True)
                self.dihedral_input.clear()
                self.dihedral_input.blockSignals(False)
                self.dihedral_slider.blockSignals(True)
                self.dihedral_slider.setValue(180)
                self.dihedral_slider.setEnabled(False)
                self.dihedral_slider.blockSignals(False)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"Selected: {symbol1}({self.atom1_idx}) - ?")
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            self.add_selection_label(self.atom1_idx, "1")
        elif self.atom3_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(
                f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - ?"
            )
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")
        elif self.atom4_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            self.selection_label.setText(
                f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - {symbol3}({self.atom3_idx}) - ?"
            )
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")
            self.add_selection_label(self.atom3_idx, "3")
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            symbol4 = self.mol.GetAtomWithIdx(self.atom4_idx).GetSymbol()
            self.selection_label.setText(
                f"Dihedral: {symbol1}({self.atom1_idx})-{symbol2}({self.atom2_idx})-{symbol3}({self.atom3_idx})-{symbol4}({self.atom4_idx})"
            )

            # Calculate current dihedral
            current_dihedral = self.calculate_dihedral()
            self.dihedral_label.setText(f"Current dihedral: {current_dihedral:.2f}°")
            self.apply_button.setEnabled(True)
            # Update input box and slider
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
            except (AttributeError, RuntimeError, TypeError):
                pass

            # Add labels
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")
            self.add_selection_label(self.atom3_idx, "3")
            self.add_selection_label(self.atom4_idx, "4")

    def calculate_dihedral(self) -> float:
        """Calculate the current dihedral angle."""
        if not self._is_selection_complete():
            return 0.0
        assert self.atom1_idx is not None
        assert self.atom2_idx is not None
        assert self.atom3_idx is not None
        assert self.atom4_idx is not None
        return calculate_dihedral(
            self.mol.GetConformer().GetPositions(),
            self.atom1_idx,
            self.atom2_idx,
            self.atom3_idx,
            self.atom4_idx,
        )

    def _is_selection_complete(self) -> bool:
        """Check if all four atoms required for a dihedral are selected."""
        return (
            self.atom1_idx is not None
            and self.atom2_idx is not None
            and self.atom3_idx is not None
            and self.atom4_idx is not None
        )

    def on_dihedral_input_changed(self, text: str) -> None:
        """Line edit text changed, update slider."""
        if not self.dihedral_input.isEnabled() or not self.apply_button.isEnabled():
            return
        self._sync_input_to_slider(text, self.dihedral_slider, 1.0, wrap=True)

    def apply_changes(self) -> None:
        """Apply the dihedral changes to the molecule."""
        if not self._is_selection_complete():
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

        # Apply the update
        self.apply_geometry_update(new_dihedral)

        # Push Undo state AFTER modification
        self._push_undo()

    def apply_geometry_update(self, new_dihedral_deg: float) -> None:  # pylint: disable=arguments-renamed
        """Adjust the dihedral angle."""
        if not self._is_selection_complete():
            return

        conf = self.mol.GetConformer()

        # Use snapshot if available (slider dragging) to keep the rotation axis stable
        snapshot = self._snapshot_positions
        if snapshot is not None:
            positions = snapshot.copy()
        else:
            positions = conf.GetPositions()

        assert self.atom1_idx is not None
        assert self.atom2_idx is not None
        assert self.atom3_idx is not None
        assert self.atom4_idx is not None
        idx1: int = self.atom1_idx
        idx2: int = self.atom2_idx
        idx3: int = self.atom3_idx
        idx4: int = self.atom4_idx

        if self.both_groups_radio.isChecked():
            # Both ends rotate equally.
            # We use adjust_dihedral twice: once for the 4th-atom side, once for the 1st-atom side.
            current_dihedral = calculate_dihedral(positions, idx1, idx2, idx3, idx4)
            delta = new_dihedral_deg - current_dihedral

            # Shortest Path Wrapping for total delta
            if delta > 180:
                delta -= 360
            elif delta < -180:
                delta += 360

            # 1. Rotate group 4 by +half delta around 2->3
            adjust_dihedral(
                positions,
                idx1,
                idx2,
                idx3,
                idx4,
                current_dihedral + delta / 2.0,
                get_connected_group(self.mol, idx3, exclude=idx2),
            )
            # 2. Rotate group 1 by -half delta around 2->3 (which is +half around 3->2)
            # Recalculate Dihedral(4,3,2,1) which is now (current + delta/2)
            # Target (current + delta)
            adjust_dihedral(
                positions,
                idx4,
                idx3,
                idx2,
                idx1,
                current_dihedral + delta,
                get_connected_group(self.mol, idx2, exclude=idx3),
            )
        elif self.rotate_atom_radio.isChecked():
            # Move only atom 4
            adjust_dihedral(positions, idx1, idx2, idx3, idx4, new_dihedral_deg, {idx4})
        else:
            # Default: rotate atom 4 and its connected sub-structure
            atoms_to_move = get_connected_group(self.mol, idx3, exclude=idx2)
            adjust_dihedral(
                positions, idx1, idx2, idx3, idx4, new_dihedral_deg, atoms_to_move
            )

        # Write updated positions back using inherited helper
        self._update_molecule_geometry(positions)
