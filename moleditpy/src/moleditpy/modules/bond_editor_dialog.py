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
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QPushButton, QMessageBox
)
from PyQt6.QtCore import Qt

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

class BondEditorDialog(Dialog3DPickingMixin, QDialog):
    """
    A dialog for manually overriding the molecular topology.
    Useful for transition states or hypervalent coordination spheres.
    Click two atoms in the 3D view to select them, then choose a bond order.
    """
    def __init__(self, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.main_window = main_window
        self.mol = main_window.current_mol
        self.atom1_idx = None
        self.atom2_idx = None
        self.setWindowTitle("Edit Bond")

        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click two atoms in the 3D view to select a bond, then choose the bond order."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Bond Order
        order_layout = QHBoxLayout()
        order_layout.addWidget(QLabel("Bond Order:"))
        self.order_combo = QComboBox()
        self.order_combo.addItems(["Single", "Double", "Triple", "Delete (0)"])
        order_layout.addWidget(self.order_combo)
        layout.addLayout(order_layout)

        # Pre-populate from existing 3D selection
        self._populate_from_selection()

        # Buttons
        btn_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        btn_layout.addWidget(self.clear_button)

        btn_layout.addStretch()

        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        btn_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        btn_layout.addWidget(close_button)

        layout.addLayout(btn_layout)

        # Connect to main window's picker
        self.picker_connection = None
        if self.mol:
            self.enable_picking()

    def _populate_from_selection(self):
        """Pre-populate from existing 3D selection if available."""
        if hasattr(self.main_window, "selected_atoms_3d") and len(self.main_window.selected_atoms_3d) == 2:
            sel_list = list(self.main_window.selected_atoms_3d)
            if self.mol:
                self.atom1_idx = sel_list[0]
                self.atom2_idx = sel_list[1]
                self.update_display()
                self.show_atom_labels()

    def apply_changes(self):
        if self.atom1_idx is None or self.atom2_idx is None:
            QMessageBox.warning(self, "Invalid Input", "Please select two atoms first.")
            return

        if self.atom1_idx == self.atom2_idx:
            QMessageBox.warning(self, "Invalid Input", "Atom 1 and Atom 2 cannot be the same.")
            return

        # Get original 2D atom IDs for the data model
        try:
            atom1 = self.mol.GetAtomWithIdx(self.atom1_idx)
            atom2 = self.mol.GetAtomWithIdx(self.atom2_idx)
            id1 = atom1.GetIntProp("_original_atom_id")
            id2 = atom2.GetIntProp("_original_atom_id")
        except Exception:
            QMessageBox.warning(self, "Error", "Could not resolve original atom IDs. Hydrogen atoms added during 3D conversion cannot be used for bond editing.")
            return

        order_idx = self.order_combo.currentIndex()
        order_map = [1.0, 2.0, 3.0, 0.0]
        order = order_map[order_idx]

        # Dispatch to data logic safely
        try:
            self.main_window.apply_bond_editor_changes(id1, id2, order)
        except AttributeError:
            # Fallback: call force_bond directly on data
            try:
                self.main_window.data.force_bond(id1, id2, order)
                self.main_window.scene.reinitialize_items()
                self.main_window.scene.update()
                try:
                    self.main_window.trigger_conversion()
                except Exception as e2:
                    self.main_window.statusBar().showMessage(f"Failed to generate 3D: {e2}")
                self.main_window.push_undo_state()
                self.main_window.statusBar().showMessage(f"Applied manual bond edit between atoms.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to apply bond edit: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply bond edit: {e}")

    def on_atom_picked(self, atom_idx):
        if self.mol is None:
            return

        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        else:
            self.atom1_idx = atom_idx
            self.atom2_idx = None

        self.show_atom_labels()
        self.update_display()

    def update_display(self):
        labels = []
        if self.atom1_idx is not None and self.atom1_idx < self.mol.GetNumAtoms():
            try:
                a1 = self.mol.GetAtomWithIdx(self.atom1_idx)
                sym = a1.GetSymbol()
                try:
                    orig_id = a1.GetIntProp('_original_atom_id')
                    labels.append(f"Atom 1: {sym} (ID {orig_id + 1})")
                except Exception:
                    labels.append(f"Atom 1: {sym} (3D idx {self.atom1_idx})")
            except Exception:
                labels.append("Atom 1: ?")

        if self.atom2_idx is not None and self.atom2_idx < self.mol.GetNumAtoms():
            try:
                a2 = self.mol.GetAtomWithIdx(self.atom2_idx)
                sym = a2.GetSymbol()
                try:
                    orig_id = a2.GetIntProp('_original_atom_id')
                    labels.append(f"Atom 2: {sym} (ID {orig_id + 1})")
                except Exception:
                    labels.append(f"Atom 2: {sym} (3D idx {self.atom2_idx})")
            except Exception:
                labels.append("Atom 2: ?")

        if labels:
            self.selection_label.setText("Selected: " + ", ".join(labels))
        else:
            self.selection_label.setText("No atoms selected")

        # Enable/disable apply button
        self.apply_button.setEnabled(
            self.atom1_idx is not None and self.atom2_idx is not None
        )

    def show_atom_labels(self):
        if self.mol is None:
            return

        atoms_and_labels = []
        if self.atom1_idx is not None and self.atom1_idx < self.mol.GetNumAtoms():
            atoms_and_labels.append((self.atom1_idx, "1"))
        if self.atom2_idx is not None and self.atom2_idx < self.mol.GetNumAtoms():
            atoms_and_labels.append((self.atom2_idx, "2"))

        try:
            self.show_atom_labels_for(atoms_and_labels)
        except Exception:
            pass

    def clear_selection(self):
        self.atom1_idx = None
        self.atom2_idx = None
        self.update_display()
        try:
            self.clear_atom_labels()
        except Exception:
            pass

    def reject(self):
        self.disable_picking()
        try:
            self.clear_atom_labels()
        except Exception:
            pass
        super().reject()
        try:
            if self.main_window.current_mol:
                self.main_window.draw_molecule_3d(self.main_window.current_mol)
        except Exception:
            pass
