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

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
from rdkit.Chem import AllChem, rdMolTransforms

from .dialog_3d_picking_mixin import Dialog3DPickingMixin


class ConstrainedOptimizationDialog(Dialog3DPickingMixin, QDialog):
    """Dialog for constrained optimization."""

    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = []  # Using a list because order matters
        self.constraints = []  # (type, atoms_indices, value)
        self.constraint_labels = []  # 3D label actors
        self.init_ui()
        self.enable_picking()

        # Load existing constraints from MainWindow
        if self.main_window.constraints_3d:
            self.constraint_table.blockSignals(True)  # Block signals during loading
            try:
                # Load into self.constraints as (Type, (Idx...), Value, Force) tuples
                for const_data in self.main_window.constraints_3d:
                    # Support 3- or 4-element constraints for backward compatibility
                    if len(const_data) == 4:
                        const_type, atom_indices, value, force_const = const_data
                    else:
                        const_type, atom_indices, value = const_data
                        force_const = 1.0e5  # Default value

                    # Convert to tuple and add to internal list
                    self.constraints.append(
                        (const_type, tuple(atom_indices), value, force_const)
                    )

                    row_count = self.constraint_table.rowCount()
                    self.constraint_table.insertRow(row_count)

                    value_str = ""
                    if const_type == "Distance":
                        value_str = f"{value:.3f}"
                    else:
                        value_str = f"{value:.2f}"

                    # Column 0 (Type)
                    item_type = QTableWidgetItem(const_type)
                    item_type.setFlags(
                        Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                    )
                    item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 0, item_type)

                    # Column 1 (Atom Indices)
                    item_indices = QTableWidgetItem(str(atom_indices))
                    item_indices.setFlags(
                        Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                    )
                    item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 1, item_indices)

                    # Column 2 (Value)
                    item_value = QTableWidgetItem(value_str)
                    item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 2, item_value)

                    # Column 3 (Force)
                    item_force = QTableWidgetItem(f"{force_const:.2e}")
                    item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 3, item_force)
            finally:
                self.constraint_table.blockSignals(False)

            # <<< Load current optimization settings from MainWindow as defaults >>>
        try:
            # Add fallback for None case
            current_method_str = self.main_window.optimization_method or "MMFF_RDKIT"
            current_method = current_method_str.upper()

            # 1. UFF_RDKIT
            if current_method == "UFF_RDKIT":
                self.ff_combo.setCurrentText("UFF")

            # 2. MMFF94_RDKIT (MMFF94)
            elif current_method == "MMFF94_RDKIT":
                self.ff_combo.setCurrentText("MMFF94")

            # 3. MMFF_RDKIT (MMFF94s) - This is the default
            elif current_method == "MMFF_RDKIT":
                self.ff_combo.setCurrentText("MMFF94s")

            # 4. (Fallback from old config files, etc.)
            elif "UFF" in current_method:
                self.ff_combo.setCurrentText("UFF")
            elif "MMFF94S" in current_method:
                self.ff_combo.setCurrentText("MMFF94s")
            elif "MMFF94" in current_method:  # Already handled MMFF94_RDKIT above
                self.ff_combo.setCurrentText("MMFF94")

            # 5. Default
            else:
                self.ff_combo.setCurrentText("MMFF94s")

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            print(f"Could not set default force field: {e}")

    def init_ui(self):
        self.setWindowTitle("Constrained Optimization")
        self.setModal(False)
        self.resize(450, 500)
        layout = QVBoxLayout(self)

        # 1. Instructions
        instruction_label = QLabel(
            "Select 2-4 atoms to add a constraint. Select constraints in the table to remove them."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # 2. Optimization method and Force Constant
        form_layout = QFormLayout()
        self.ff_combo = QComboBox()
        self.ff_combo.addItems(["MMFF94s", "MMFF94", "UFF"])
        form_layout.addRow("Force Field:", self.ff_combo)

        # Force Constant setting
        self.force_const_input = QLineEdit("1.0e5")
        self.force_const_input.setToolTip(
            "Force constant for constraints (default: 1.0e5)"
        )
        form_layout.addRow("Force Constant:", self.force_const_input)

        layout.addLayout(form_layout)

        # 3. Selected atoms
        self.selection_label = QLabel("Selected atoms: None")
        layout.addWidget(self.selection_label)

        # 4. Constraints table
        self.constraint_table = QTableWidget()
        self.constraint_table.setColumnCount(4)
        self.constraint_table.setHorizontalHeaderLabels(
            ["Type", "Atom Indices", "Value (Å or °)", "Force"]
        )
        self.constraint_table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows
        )
        # Change edit triggers (e.g., double click)
        self.constraint_table.setEditTriggers(
            QTableWidget.EditTrigger.DoubleClicked
            | QTableWidget.EditTrigger.SelectedClicked
            | QTableWidget.EditTrigger.EditKeyPressed
        )
        self.constraint_table.itemSelectionChanged.connect(self.show_constraint_labels)
        self.constraint_table.cellChanged.connect(self.on_cell_changed)

        self.constraint_table.setStyleSheet("""
            QTableWidget QLineEdit {
                background-color: white;
                color: black;
                border: none;
            }
        """)

        layout.addWidget(self.constraint_table)

        # 5. Buttons (Add / Remove)
        button_layout = QHBoxLayout()
        self.add_button = QPushButton("Add Constraint")
        self.add_button.clicked.connect(self.add_constraint)
        self.add_button.setEnabled(False)
        button_layout.addWidget(self.add_button)

        self.remove_button = QPushButton("Remove Selected")
        self.remove_button.clicked.connect(self.remove_constraint)
        self.remove_button.setEnabled(False)
        button_layout.addWidget(self.remove_button)
        layout.addLayout(button_layout)

        self.remove_all_button = QPushButton("Remove All")
        self.remove_all_button.clicked.connect(self.remove_all_constraints)
        button_layout.addWidget(self.remove_all_button)

        # 6. Main buttons (Optimize / Close)
        main_buttons = QHBoxLayout()
        main_buttons.addStretch()
        self.optimize_button = QPushButton("Optimize")
        self.optimize_button.clicked.connect(self.apply_optimization)
        main_buttons.addWidget(self.optimize_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        main_buttons.addWidget(close_button)
        layout.addLayout(main_buttons)

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            if len(self.selected_atoms) >= 4:
                self.selected_atoms.pop(0)  # Up to 4 atoms
            self.selected_atoms.append(atom_idx)

        self.show_selection_labels()
        self.update_selection_display()

    def update_selection_display(self):
        self.show_selection_labels()
        n = len(self.selected_atoms)

        atom_str = ", ".join(map(str, self.selected_atoms))
        prefix = ""
        can_add = False

        if n == 0:
            prefix = "Selected atoms: None"
            atom_str = ""  # Set atom_str to empty
        elif n == 1:
            prefix = "Selected atoms: "
        elif n == 2:
            prefix = "Selected atoms: <b>Distance</b> "
            can_add = True
        elif n == 3:
            prefix = "Selected atoms: <b>Angle</b> "
            can_add = True
        elif n == 4:
            prefix = "Selected atoms: <b>Torsion</b> "
            can_add = True
        else:  # n > 4
            prefix = "Selected atoms (max 4): "

        # Set label text
        if n == 0:
            self.selection_label.setText(prefix)
        else:
            self.selection_label.setText(f"{prefix}[{atom_str}]")

        # Button text is fixed
        self.add_button.setText("Add Constraint")
        # Set button enabled state
        self.add_button.setEnabled(can_add)

    def add_constraint(self):
        n = len(self.selected_atoms)
        conf = self.mol.GetConformer()

        # Retrieve Force Constant
        try:
            force_const = float(self.force_const_input.text())
        except ValueError:
            QMessageBox.warning(
                self, "Warning", "Invalid Force Constant. Using default 1.0e5."
            )
            force_const = 1.0e5

        if n == 2:
            constraint_type = "Distance"
            value = conf.GetAtomPosition(self.selected_atoms[0]).Distance(
                conf.GetAtomPosition(self.selected_atoms[1])
            )
            value_str = f"{value:.3f}"
        elif n == 3:
            constraint_type = "Angle"
            value = rdMolTransforms.GetAngleDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        elif n == 4:
            constraint_type = "Torsion"
            value = rdMolTransforms.GetDihedralDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        else:
            return

        atom_indices = tuple(self.selected_atoms)

        # Check for duplicate constraints (same atom indices)
        for const in self.constraints:
            if const[0] == constraint_type and const[1] == atom_indices:
                QMessageBox.warning(
                    self, "Warning", "This exact constraint already exists."
                )
                return

        self.constraints.append((constraint_type, atom_indices, value, force_const))

        # Update table
        row_count = self.constraint_table.rowCount()
        self.constraint_table.insertRow(row_count)

        # --- Column 0 (Type) ---
        item_type = QTableWidgetItem(constraint_type)
        # Set non-editable flags
        item_type.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 0, item_type)

        # --- Column 1 (Atom Indices) ---
        item_indices = QTableWidgetItem(str(atom_indices))
        # Set non-editable flags
        item_indices.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 1, item_indices)

        # --- Column 2 (Value) ---
        item_value = QTableWidgetItem(value_str)
        item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # Editable flags (enabled by default)
        self.constraint_table.setItem(row_count, 2, item_value)

        # --- Column 3 (Force) ---
        item_force = QTableWidgetItem(f"{force_const:.2e}")
        item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # Editable
        self.constraint_table.setItem(row_count, 3, item_force)

        # Clear selection
        self.selected_atoms.clear()
        self.update_selection_display()

    def remove_constraint(self):
        selected_rows = sorted(
            list(set(index.row() for index in self.constraint_table.selectedIndexes())),
            reverse=True,
        )
        if not selected_rows:
            return

        self.constraint_table.blockSignals(True)

        for row in selected_rows:
            self.constraints.pop(row)
            self.constraint_table.removeRow(row)

        self.constraint_table.blockSignals(False)

        self.clear_constraint_labels()

    def remove_all_constraints(self):
        """Clear all constraints."""
        if not self.constraints:
            return

        # Clear internal list
        self.constraints.clear()

        # Remove all rows from the table
        self.constraint_table.blockSignals(True)
        self.constraint_table.setRowCount(0)
        self.constraint_table.blockSignals(False)

        # Clear 3D labels
        self.clear_constraint_labels()

        # Disable selection-based button
        self.remove_button.setEnabled(False)

    def show_constraint_labels(self):
        self.clear_constraint_labels()
        selected_items = self.constraint_table.selectedItems()
        if not selected_items:
            self.remove_button.setEnabled(False)
            return

        self.remove_button.setEnabled(True)

        # Get the constraint from the selected row (first selection only)
        try:
            row = selected_items[0].row()
            constraint_type, atom_indices, value, force_const = self.constraints[row]
        except (IndexError, TypeError, ValueError):
            # Unpack 3-element tuple for old format constraints
            try:
                constraint_type, atom_indices, value = self.constraints[row]
            except (IndexError, TypeError):
                return

        labels = []
        if constraint_type == "Distance":
            labels = ["A1", "A2"]
        elif constraint_type == "Angle":
            labels = ["A1", "A2 (V)", "A3"]
        elif constraint_type == "Torsion":
            labels = ["A1", "A2", "A3", "A4"]

        positions = []
        texts = []
        for i, atom_idx in enumerate(atom_indices):
            positions.append(self.main_window.atom_positions_3d[atom_idx])
            texts.append(labels[i])

        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions,
                texts,
                point_size=20,
                font_size=12,
                text_color="cyan",
                always_visible=True,
            )
            self.constraint_labels.append(label_actor)

    def clear_constraint_labels(self):
        for label_actor in self.constraint_labels:
            try:
                self.main_window.plotter.remove_actor(label_actor)
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress errors during constraint label removal

        self.constraint_labels = []

    def apply_optimization(self):
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No valid 3D molecule found.")
            return

        ff_name = self.ff_combo.currentText()
        conf = self.mol.GetConformer()

        try:
            ignore_interfrag = not self.main_window.settings.get(
                "optimize_intermolecular_interaction_rdkit", True
            )
            if ff_name.startswith("MMFF"):
                props = AllChem.MMFFGetMoleculeProperties(self.mol, mmffVariant=ff_name)
                ff = AllChem.MMFFGetMoleculeForceField(
                    self.mol,
                    props,
                    confId=0,
                    ignoreInterfragInteractions=ignore_interfrag,
                )
                add_dist_constraint = ff.MMFFAddDistanceConstraint
                add_angle_constraint = ff.MMFFAddAngleConstraint
                add_torsion_constraint = ff.MMFFAddTorsionConstraint
            else:  # UFF
                ff = AllChem.UFFGetMoleculeForceField(
                    self.mol, confId=0, ignoreInterfragInteractions=ignore_interfrag
                )
                add_dist_constraint = ff.UFFAddDistanceConstraint
                add_angle_constraint = ff.UFFAddAngleConstraint
                add_torsion_constraint = ff.UFFAddTorsionConstraint

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(
                self, "Error", f"Failed to initialize force field {ff_name}: {e}"
            )
            return

        # Add constraints
        try:
            for constraint in self.constraints:
                # Support 4- or 3-element constraints for backward compatibility
                if len(constraint) == 4:
                    const_type, atoms, value, force_const = constraint
                else:
                    const_type, atoms, value = constraint
                    force_const = 1.0e5  # Default value

                if const_type == "Distance":
                    # C++ signature: (self, idx1, idx2, bool relative, minLen, maxLen, forceConst)
                    add_dist_constraint(
                        int(atoms[0]),
                        int(atoms[1]),
                        False,
                        float(value),
                        float(value),
                        float(force_const),
                    )
                elif const_type == "Angle":
                    # C++ signature: (self, idx1, idx2, idx3, bool relative, minDeg, maxDeg, forceConst)
                    add_angle_constraint(
                        int(atoms[0]),
                        int(atoms[1]),
                        int(atoms[2]),
                        False,
                        float(value),
                        float(value),
                        float(force_const),
                    )
                elif const_type == "Torsion":
                    # C++ signature: (self, idx1, idx2, idx3, idx4, bool relative, minDeg, maxDeg, forceConst)
                    add_torsion_constraint(
                        int(atoms[0]),
                        int(atoms[1]),
                        int(atoms[2]),
                        int(atoms[3]),
                        False,
                        float(value),
                        float(value),
                        float(force_const),
                    )

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(self, "Error", f"Failed to add constraints: {e}")
            print(e)
            return

        # Execute optimization
        try:
            self.main_window.statusBar().showMessage(
                f"Running constrained {ff_name} optimization..."
            )
            ff.Minimize(maxIts=20000)

            # Apply optimized coordinates to the main window's numpy array
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                self.main_window.atom_positions_3d[i] = [pos.x, pos.y, pos.z]

            # Update 3D view
            self.main_window.draw_molecule_3d(self.mol)
            self.main_window.update_chiral_labels()
            self.main_window.push_undo_state()
            self.main_window.statusBar().showMessage(
                "Constrained optimization finished."
            )

            try:
                constrained_method_name = f"Constrained_{ff_name}"
                self.main_window.last_successful_optimization_method = (
                    constrained_method_name
                )
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Failed to set last_successful_optimization_method: {e}")

            # Save constraints list to MainWindow on success (same logic as reject)
            try:
                # Save as list for JSON compatibility
                json_safe_constraints = []
                for const in self.constraints:
                    # 4-element constraint
                    if len(const) == 4:
                        json_safe_constraints.append(
                            [const[0], list(const[1]), const[2], const[3]]
                        )
                    else:
                        # Add default force for old 3-element format
                        json_safe_constraints.append(
                            [const[0], list(const[1]), const[2], 1.0e5]
                        )

                # Update MainWindow only if changed
                if self.main_window.constraints_3d != json_safe_constraints:
                    self.main_window.constraints_3d = json_safe_constraints
                    self.main_window.has_unsaved_changes = (
                        True  # Mark as unsaved changes
                    )
                    self.main_window.update_window_title()

            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Failed to save constraints post-optimization: {e}")

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(self, "Error", f"Optimization failed: {e}")

    def closeEvent(self, event):
        self.reject()
        event.accept()

    def reject(self):
        self.clear_constraint_labels()
        self.clear_selection_labels()
        self.disable_picking()

        # Save constraints list to MainWindow when closing the dialog
        try:
            # Save as list for JSON compatibility
            json_safe_constraints = []
            for const in self.constraints:
                # Tuple to list conversion
                if len(const) == 4:
                    json_safe_constraints.append(
                        [const[0], list(const[1]), const[2], const[3]]
                    )
                else:
                    # Add default force for old 3-element format
                    json_safe_constraints.append(
                        [const[0], list(const[1]), const[2], 1.0e5]
                    )

            # Update MainWindow only if changed
            if self.main_window.constraints_3d != json_safe_constraints:
                self.main_window.constraints_3d = json_safe_constraints
                self.main_window.has_unsaved_changes = True  # Mark as unsaved changes
                self.main_window.update_window_title()

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            print(f"Failed to save constraints to main window: {e}")

        super().reject()

    def clear_selection(self):
        """Clear selection (called from Mixin when clicking outside an atom)"""
        self.selected_atoms.clear()
        self.clear_selection_labels()
        self.update_selection_display()

    def show_selection_labels(self):
        """Display labels on selected atoms."""
        self.clear_selection_labels()

        if not hasattr(self, "selection_labels"):
            self.selection_labels = []

        if (
            not hasattr(self.main_window, "atom_positions_3d")
            or self.main_window.atom_positions_3d is None
        ):
            return  # Do nothing if 3D coordinates are missing

        max_idx = len(self.main_window.atom_positions_3d) - 1
        positions = []
        texts = []

        for i, atom_idx in enumerate(self.selected_atoms):
            if atom_idx is not None and 0 <= atom_idx <= max_idx:
                positions.append(self.main_window.atom_positions_3d[atom_idx])
                texts.append(f"A{i + 1}")
            elif atom_idx is not None:
                # Log invalid index (for debugging)
                print(
                    f"Warning: Invalid atom index {atom_idx} in show_selection_labels"
                )

        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions,
                texts,
                point_size=20,
                font_size=12,
                text_color="yellow",
                always_visible=True,
            )
            # Consider case where add_point_labels returns a list
            if isinstance(label_actor, list):
                self.selection_labels.extend(label_actor)
            else:
                self.selection_labels.append(label_actor)

    def on_cell_changed(self, row, column):
        """Update internal data when a table cell is edited."""

        # Handle only Value (col 2) and Force (col 3) columns
        if column not in [2, 3]:
            return

        try:
            # Retrieve text from the modified item
            item = self.constraint_table.item(row, column)
            if not item:
                return

            new_value_str = item.text()
            new_value = float(new_value_str)

            # Update internal constraints list
            old_constraint = self.constraints[row]

            # Support 3- or 4-element constraints for backward compatibility
            if len(old_constraint) == 4:
                if column == 2:  # Value column
                    self.constraints[row] = (
                        old_constraint[0],
                        old_constraint[1],
                        new_value,
                        old_constraint[3],
                    )
                elif column == 3:  # Force column
                    self.constraints[row] = (
                        old_constraint[0],
                        old_constraint[1],
                        old_constraint[2],
                        new_value,
                    )
            else:
                # Old 3-element format
                if column == 2:  # Value column
                    self.constraints[row] = (
                        old_constraint[0],
                        old_constraint[1],
                        new_value,
                        1.0e5,
                    )
                elif column == 3:  # Force column (new addition)
                    self.constraints[row] = (
                        old_constraint[0],
                        old_constraint[1],
                        old_constraint[2],
                        new_value,
                    )

        except (ValueError, TypeError):
            # Case of invalid input (non-numeric)
            # Restore original value to table
            self.constraint_table.blockSignals(True)

            if column == 2:  # Value column
                old_value = self.constraints[row][2]
                if self.constraints[row][0] == "Distance":
                    item.setText(f"{old_value:.3f}")
                else:
                    item.setText(f"{old_value:.2f}")
            elif column == 3:  # Force column
                old_force = (
                    self.constraints[row][3]
                    if len(self.constraints[row]) == 4
                    else 1.0e5
                )
                item.setText(f"{old_force:.2e}")

            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.constraint_table.blockSignals(False)

            QMessageBox.warning(
                self, "Invalid Value", "Please enter a valid floating-point number."
            )
        except IndexError as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress sync errors between table and constraints list

    def keyPressEvent(self, event):
        """Handle keyboard events (Delete/Backspace to remove, Enter to optimize)."""
        key = event.key()

        # Check if Delete or Backspace was pressed
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            # Check if table has focus or items are selected
            if (
                self.constraint_table.hasFocus()
                or len(self.constraint_table.selectedIndexes()) > 0
            ):
                self.remove_constraint()
                event.accept()
                return

        # Check if Enter/Return was pressed (execute optimization)
        if key == Qt.Key.Key_Return or key == Qt.Key.Key_Enter:
            # Ensure table is not in editing state
            if self.constraint_table.state() != QAbstractItemView.State.EditingState:
                if self.optimize_button.isEnabled():
                    self.apply_optimization()
                event.accept()
                return

        # Default processing for other keys
        super().keyPressEvent(event)
