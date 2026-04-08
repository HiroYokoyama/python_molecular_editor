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
from PyQt6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

try:
    from .base_picking_dialog import BasePickingDialog
except ImportError:
    from moleditpy_linux.ui.base_picking_dialog import BasePickingDialog

_TAB_ABSOLUTE = 0
_TAB_DELTA = 1


class TranslationDialog(BasePickingDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        super().__init__(mol, main_window, parent)
        self.selected_atoms = set()

        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)

        self.init_ui()

        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def init_ui(self):
        self.setWindowTitle("Translate Atoms")
        self.setModal(False)
        layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_absolute_tab(), "Absolute")
        self.tabs.addTab(self._build_delta_tab(), "Delta")
        self.tabs.currentChanged.connect(self._on_tab_changed)
        layout.addWidget(self.tabs)

        # Shared Close button row
        close_row = QHBoxLayout()
        close_row.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        close_row.addWidget(close_btn)
        layout.addLayout(close_row)

        self.enable_picking()

    def _build_absolute_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        instr = QLabel(
            "Click one atom to select it, then specify its target absolute coordinates (Å)."
        )
        instr.setWordWrap(True)
        layout.addWidget(instr)

        self.abs_selection_label = QLabel("No atom selected")
        layout.addWidget(self.abs_selection_label)

        coord_row = QHBoxLayout()
        coord_row.addWidget(QLabel("X:"))
        self.abs_x_input = QLineEdit("0.000")
        coord_row.addWidget(self.abs_x_input)
        coord_row.addWidget(QLabel("Y:"))
        self.abs_y_input = QLineEdit("0.000")
        coord_row.addWidget(self.abs_y_input)
        coord_row.addWidget(QLabel("Z:"))
        self.abs_z_input = QLineEdit("0.000")
        coord_row.addWidget(self.abs_z_input)
        layout.addLayout(coord_row)

        self.move_mol_checkbox = QCheckBox("Move entire molecule")
        self.move_mol_checkbox.setChecked(True)
        self.move_mol_checkbox.stateChanged.connect(self._on_move_mol_toggled)
        layout.addWidget(self.move_mol_checkbox)

        origin_btn = QPushButton("Set to Origin")
        origin_btn.setToolTip("Set target coordinates to the origin")
        origin_btn.clicked.connect(self._set_origin)
        layout.addWidget(origin_btn)

        btn_row = QHBoxLayout()
        abs_clear_btn = QPushButton("Clear Selection")
        abs_clear_btn.clicked.connect(self._abs_clear_selection)
        btn_row.addWidget(abs_clear_btn)
        btn_row.addStretch()
        self.abs_apply_btn = QPushButton("Move Molecule")
        self.abs_apply_btn.clicked.connect(self.apply_absolute)
        self.abs_apply_btn.setEnabled(False)
        btn_row.addWidget(self.abs_apply_btn)
        layout.addLayout(btn_row)

        layout.addStretch()
        return widget

    def _build_delta_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        instr = QLabel(
            "Click atoms in the 3D view to select them, then specify the translation vector (Å)."
        )
        instr.setWordWrap(True)
        layout.addWidget(instr)

        self.delta_selection_label = QLabel("No atoms selected")
        layout.addWidget(self.delta_selection_label)

        vector_row = QHBoxLayout()
        vector_row.addWidget(QLabel("dX:"))
        self.dx_input = QLineEdit("0.0")
        vector_row.addWidget(self.dx_input)
        vector_row.addWidget(QLabel("dY:"))
        self.dy_input = QLineEdit("0.0")
        vector_row.addWidget(self.dy_input)
        vector_row.addWidget(QLabel("dZ:"))
        self.dz_input = QLineEdit("0.0")
        vector_row.addWidget(self.dz_input)
        layout.addLayout(vector_row)

        btn_row = QHBoxLayout()
        clear_btn = QPushButton("Clear Selection")
        clear_btn.clicked.connect(self.clear_selection)
        btn_row.addWidget(clear_btn)

        select_all_btn = QPushButton("Select All Atoms")
        select_all_btn.setToolTip("Select all atoms in the molecule")
        select_all_btn.clicked.connect(self.select_all_atoms)
        btn_row.addWidget(select_all_btn)

        btn_row.addStretch()
        self.apply_button = QPushButton("Apply Translation")
        self.apply_button.clicked.connect(self.apply_translation)
        self.apply_button.setEnabled(False)
        btn_row.addWidget(self.apply_button)
        layout.addLayout(btn_row)

        layout.addStretch()
        return widget

    # ------------------------------------------------------------------
    # Tab switching
    # ------------------------------------------------------------------

    def _on_tab_changed(self, index):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    # ------------------------------------------------------------------
    # Atom picking dispatch
    # ------------------------------------------------------------------

    def on_atom_picked(self, atom_idx):
        if self.tabs.currentIndex() == _TAB_ABSOLUTE:
            self._abs_on_atom_picked(atom_idx)
        else:
            self._delta_on_atom_picked(atom_idx)

    def _abs_on_atom_picked(self, atom_idx):
        # Enforce single selection: replace previous atom
        self.selected_atoms = {atom_idx}
        self._populate_abs_inputs_from_atom(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def _delta_on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    # ------------------------------------------------------------------
    # Absolute tab helpers
    # ------------------------------------------------------------------

    def _populate_abs_inputs_from_atom(self, atom_idx):
        pos = self.main_window.view_3d_manager.current_mol.GetConformer().GetPositions()[atom_idx]
        self.abs_x_input.setText(f"{pos[0]:.4f}")
        self.abs_y_input.setText(f"{pos[1]:.4f}")
        self.abs_z_input.setText(f"{pos[2]:.4f}")

    def _abs_clear_selection(self):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.abs_x_input.setText("0.000")
        self.abs_y_input.setText("0.000")
        self.abs_z_input.setText("0.000")
        self.update_display()

    def _set_origin(self):
        self.abs_x_input.setText("0.0000")
        self.abs_y_input.setText("0.0000")
        self.abs_z_input.setText("0.0000")

    def _on_move_mol_toggled(self, state):
        label = "Move Molecule" if self.move_mol_checkbox.isChecked() else "Move Atom"
        self.abs_apply_btn.setText(label)

    def apply_absolute(self):
        self.mol = self.main_window.view_3d_manager.current_mol
        if len(self.selected_atoms) != 1:
            QMessageBox.warning(self, "Warning", "Please select exactly one atom.")
            return

        try:
            tx = float(self.abs_x_input.text())
            ty = float(self.abs_y_input.text())
            tz = float(self.abs_z_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid numbers for X, Y, Z.")
            return

        atom_idx = next(iter(self.selected_atoms))
        positions = self.mol.GetConformer().GetPositions()
        current = positions[atom_idx]
        delta = np.array([tx, ty, tz]) - current

        if np.allclose(delta, 0):
            return

        if self.move_mol_checkbox.isChecked():
            positions += delta
        else:
            positions[atom_idx] += delta

        self._update_molecule_geometry(positions)
        self._push_undo()
        self.show_atom_labels()

    # ------------------------------------------------------------------
    # Delta tab methods (unchanged logic)
    # ------------------------------------------------------------------

    def clear_selection(self):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def select_all_atoms(self):
        try:
            if hasattr(self, "mol") and self.mol is not None:
                self.selected_atoms = set(range(self.mol.GetNumAtoms()))
            else:
                self.selected_atoms = (
                    set(self.main_window.state_manager.data.atoms.keys())
                    if hasattr(self.main_window.state_manager, "data")
                    else set()
                )
            self.show_atom_labels()
            self.update_display()
        except (AttributeError, RuntimeError, TypeError, KeyError) as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def apply_translation(self):
        self.mol = self.main_window.view_3d_manager.current_mol
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

        translation_vec = np.array([dx, dy, dz])
        positions = self.mol.GetConformer().GetPositions()
        for atom_idx in self.selected_atoms:
            positions[atom_idx] += translation_vec

        self._update_molecule_geometry(positions)
        self._push_undo()
        self.show_atom_labels()

    # ------------------------------------------------------------------
    # Shared display update
    # ------------------------------------------------------------------

    def update_display(self):
        tab = self.tabs.currentIndex()
        count = len(self.selected_atoms)

        if tab == _TAB_ABSOLUTE:
            if count == 0:
                self.abs_selection_label.setText("Click one atom to select it")
                self.abs_apply_btn.setEnabled(False)
            else:
                atom_idx = next(iter(self.selected_atoms))
                sym = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                self.abs_selection_label.setText(f"Selected: atom {atom_idx} ({sym})")
                self.abs_apply_btn.setEnabled(True)
        else:
            if count == 0:
                self.delta_selection_label.setText("Click atoms to select (minimum 1 required)")
                self.apply_button.setEnabled(False)
            else:
                self.delta_selection_label.setText(f"Selected {count} atom{'s' if count != 1 else ''}")
                self.apply_button.setEnabled(True)

    def show_atom_labels(self):
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            pairs = [(idx, str(i + 1)) for i, idx in enumerate(sorted_atoms)]
            self.show_atom_labels_for(pairs)
        else:
            self.clear_atom_labels()
