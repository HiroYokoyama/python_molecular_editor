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
import pyvista as pv
from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtWidgets import (
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

try:
    from .constants import VDW_RADII, pt
except ImportError:
    from modules.constants import VDW_RADII, pt

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin


class MoveGroupDialog(Dialog3DPickingMixin, QDialog):  # pragma: no cover
    """Dialog to select a connected molecular group and perform translation/rotation."""

    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()
        self.group_atoms = set()  # All atoms connected to selected atoms
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Move Group")
        self.setModal(False)
        self.resize(300, 400)
        layout = QVBoxLayout(self)

        # Drag state management
        self.is_dragging_group = False
        self.drag_start_pos = None
        self.mouse_moved_during_drag = False

        # Instructions
        instruction_label = QLabel(
            "Click an atom in the 3D view to select its connected molecule group.\n"
            "Left-drag: Move the group\n"
            "Right-drag: Rotate the group around its center"
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected group display
        self.selection_label = QLabel("No group selected")
        layout.addWidget(self.selection_label)

        # Translation controls
        trans_group = QLabel("Translation (Å):")
        trans_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(trans_group)

        trans_layout = QGridLayout()
        self.x_trans_input = QLineEdit("0.0")
        self.y_trans_input = QLineEdit("0.0")
        self.z_trans_input = QLineEdit("0.0")

        # Execute apply_translation on Enter key
        self.x_trans_input.returnPressed.connect(self.apply_translation)
        self.y_trans_input.returnPressed.connect(self.apply_translation)
        self.z_trans_input.returnPressed.connect(self.apply_translation)

        trans_layout.addWidget(QLabel("X:"), 0, 0)
        trans_layout.addWidget(self.x_trans_input, 0, 1)
        trans_layout.addWidget(QLabel("Y:"), 1, 0)
        trans_layout.addWidget(self.y_trans_input, 1, 1)
        trans_layout.addWidget(QLabel("Z:"), 2, 0)
        trans_layout.addWidget(self.z_trans_input, 2, 1)

        trans_button_layout = QHBoxLayout()
        reset_trans_button = QPushButton("Reset")
        reset_trans_button.clicked.connect(self.reset_translation_inputs)
        trans_button_layout.addWidget(reset_trans_button)

        apply_trans_button = QPushButton("Apply Translation")
        apply_trans_button.clicked.connect(self.apply_translation)
        trans_button_layout.addWidget(apply_trans_button)

        trans_layout.addLayout(trans_button_layout, 3, 0, 1, 2)

        layout.addLayout(trans_layout)

        layout.addSpacing(10)

        # Rotation controls
        rot_group = QLabel("Rotation (degrees):")
        rot_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(rot_group)

        rot_layout = QGridLayout()
        self.x_rot_input = QLineEdit("0.0")
        self.y_rot_input = QLineEdit("0.0")
        self.z_rot_input = QLineEdit("0.0")

        # Execute apply_rotation on Enter key
        self.x_rot_input.returnPressed.connect(self.apply_rotation)
        self.y_rot_input.returnPressed.connect(self.apply_rotation)
        self.z_rot_input.returnPressed.connect(self.apply_rotation)

        rot_layout.addWidget(QLabel("Around X:"), 0, 0)
        rot_layout.addWidget(self.x_rot_input, 0, 1)
        rot_layout.addWidget(QLabel("Around Y:"), 1, 0)
        rot_layout.addWidget(self.y_rot_input, 1, 1)
        rot_layout.addWidget(QLabel("Around Z:"), 2, 0)
        rot_layout.addWidget(self.z_rot_input, 2, 1)

        rot_button_layout = QHBoxLayout()
        reset_rot_button = QPushButton("Reset")
        reset_rot_button.clicked.connect(self.reset_rotation_inputs)
        rot_button_layout.addWidget(reset_rot_button)

        apply_rot_button = QPushButton("Apply Rotation")
        apply_rot_button.clicked.connect(self.apply_rotation)
        rot_button_layout.addWidget(apply_rot_button)

        rot_layout.addLayout(rot_button_layout, 3, 0, 1, 2)

        layout.addLayout(rot_layout)

        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)

        button_layout.addStretch()

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        # Enable picking to handle atom selection
        self.enable_picking()

    def eventFilter(self, obj, event):
        """Mouse event handling in 3D view - delegate to CustomInteractorStyle if a group is selected."""
        if obj == self.main_window.plotter.interactor:
            # Prevent state confusion from double/triple clicks
            if event.type() == QEvent.Type.MouseButtonDblClick:
                # Ignore double clicks and reset state
                self.is_dragging_group = False
                self.drag_start_pos = None
                self.mouse_moved_during_drag = False
                self.potential_drag = False
                if hasattr(self, "clicked_atom_for_toggle"):
                    delattr(self, "clicked_atom_for_toggle")
                return False

            if (
                event.type() == QEvent.Type.MouseButtonPress
                and event.button() == Qt.MouseButton.LeftButton
            ):
                # Clean up previous state (triple-click countermeasure)
                self.is_dragging_group = False
                self.potential_drag = False
                if hasattr(self, "clicked_atom_for_toggle"):
                    delattr(self, "clicked_atom_for_toggle")
                # Delegate to CustomInteractorStyle if a group is already selected
                if self.group_atoms:
                    return False

                # Mouse press handling
                try:
                    interactor = self.main_window.plotter.interactor
                    click_pos = interactor.GetEventPosition()

                    # Pick to identify which atom was clicked
                    picker = self.main_window.plotter.picker
                    picker.Pick(
                        click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer
                    )

                    clicked_atom_idx = None
                    if picker.GetActor() is self.main_window.atom_actor:
                        picked_position = np.array(picker.GetPickPosition())
                        distances = np.linalg.norm(
                            self.main_window.atom_positions_3d - picked_position, axis=1
                        )
                        closest_atom_idx = np.argmin(distances)

                        # Threshold check
                        if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                            atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                            if atom:
                                try:
                                    atomic_num = atom.GetAtomicNum()
                                    vdw_radius = pt.GetRvdw(atomic_num)
                                    if vdw_radius < 0.1:
                                        vdw_radius = 1.5
                                except (AttributeError, RuntimeError):
                                    vdw_radius = 1.5
                                click_threshold = vdw_radius * 1.5

                                if distances[closest_atom_idx] < click_threshold:
                                    clicked_atom_idx = int(closest_atom_idx)

                    # Handle clicked atom
                    if clicked_atom_idx is not None:
                        if self.group_atoms and clicked_atom_idx in self.group_atoms:
                            # Atom within existing group - prepare for drag
                            self.is_dragging_group = False
                            self.drag_start_pos = click_pos
                            self.drag_atom_idx = clicked_atom_idx
                            self.mouse_moved_during_drag = False
                            self.potential_drag = True  # Potential drag start
                            self.clicked_atom_for_toggle = (
                                clicked_atom_idx  # Save for toggling
                            )
                            # Allow camera operation (start drag if threshold exceeded)
                            return False
                        else:
                            # Atom outside group - select new group
                            # Manually call on_atom_picked from parent Mixin
                            self.on_atom_picked(clicked_atom_idx)
                            return True
                    else:
                        # Clicked outside atoms
                        # Allow normal camera operation even if a group exists
                        return False

                except (AttributeError, RuntimeError, ValueError) as e:
                    print(f"Error in mouse press: {e}")
                    return False

            elif event.type() == QEvent.Type.MouseMove:
                # Mouse move handling
                if (
                    getattr(self, "potential_drag", False)
                    and self.drag_start_pos
                    and not self.is_dragging_group
                ):
                    # potential_drag state: threshold check
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]

                        # Start drag if threshold is exceeded
                        drag_threshold = 5  # pixels
                        if abs(dx) > drag_threshold or abs(dy) > drag_threshold:
                            # Confirm drag start
                            self.is_dragging_group = True
                            self.potential_drag = False
                            try:
                                self.main_window.plotter.setCursor(
                                    Qt.CursorShape.ClosedHandCursor
                                )
                            except (AttributeError, RuntimeError):  # pragma: no cover
                                import traceback
                                traceback.print_exc()
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        import traceback
                        traceback.print_exc()

                    # Allow camera operation if below threshold
                    if not self.is_dragging_group:
                        return False

                if self.is_dragging_group and self.drag_start_pos:
                    # In drag mode - just record distance (no real-time update)
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()

                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]

                        if abs(dx) > 2 or abs(dy) > 2:
                            self.mouse_moved_during_drag = True
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        import traceback
                        traceback.print_exc()

                    # Consume event during drag to prevent camera rotation
                    return True

                # Hover handling (when not dragging)
                if self.group_atoms:
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        picker = self.main_window.plotter.picker
                        picker.Pick(
                            current_pos[0],
                            current_pos[1],
                            0,
                            self.main_window.plotter.renderer,
                        )

                        if picker.GetActor() is self.main_window.atom_actor:
                            picked_position = np.array(picker.GetPickPosition())
                            distances = np.linalg.norm(
                                self.main_window.atom_positions_3d - picked_position,
                                axis=1,
                            )
                            closest_atom_idx = np.argmin(distances)

                            if closest_atom_idx in self.group_atoms:
                                self.main_window.plotter.setCursor(
                                    Qt.CursorShape.OpenHandCursor
                                )
                            else:
                                self.main_window.plotter.setCursor(
                                    Qt.CursorShape.ArrowCursor
                                )
                        else:
                            self.main_window.plotter.setCursor(
                                Qt.CursorShape.ArrowCursor
                            )
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        import traceback
                        traceback.print_exc()

                # Allow camera rotation if not dragging
                return False

            elif (
                event.type() == QEvent.Type.MouseButtonRelease
                and event.button() == Qt.MouseButton.LeftButton
            ):
                # Mouse release handling
                if getattr(self, "potential_drag", False) or (
                    self.is_dragging_group and self.drag_start_pos
                ):
                    try:
                        if self.is_dragging_group and self.mouse_moved_during_drag:
                            # Drag executed - delegate to CustomInteractorStyle (do nothing)
                            pass
                        else:
                            # Mouse move below threshold = simple click
                            # Toggle selection if an atom within the group is clicked
                            if hasattr(self, "clicked_atom_for_toggle"):
                                clicked_atom = self.clicked_atom_for_toggle
                                delattr(self, "clicked_atom_for_toggle")
                                # Reset drag state before toggling
                                self.is_dragging_group = False
                                self.drag_start_pos = None
                                self.mouse_moved_during_drag = False
                                self.potential_drag = False
                                if hasattr(self, "last_drag_positions"):
                                    delattr(self, "last_drag_positions")
                                # Execute toggle process
                                self.on_atom_picked(clicked_atom)
                                try:
                                    self.main_window.plotter.setCursor(
                                        Qt.CursorShape.ArrowCursor
                                    )
                                except (AttributeError, RuntimeError):  # pragma: no cover
                                    import traceback
                                    traceback.print_exc()
                                return True

                    except (AttributeError, RuntimeError):  # pragma: no cover
                        import traceback
                        traceback.print_exc()
                    finally:
                        # Reset drag state
                        self.is_dragging_group = False
                        self.drag_start_pos = None
                        self.mouse_moved_during_drag = False
                        self.potential_drag = False
                        # Clear saved position data
                        if hasattr(self, "last_drag_positions"):
                            delattr(self, "last_drag_positions")
                        try:
                            self.main_window.plotter.setCursor(
                                Qt.CursorShape.ArrowCursor
                            )
                        except (AttributeError, RuntimeError):  # pragma: no cover
                            import traceback
                            traceback.print_exc()

                    return True  # Consume event

                # Normal release processing if not dragging
                return False

        # Pass other events to parent class
        return super().eventFilter(obj, event)

    def on_atom_picked(self, atom_idx):
        """Select the entire connected component the atom belongs to (supports multiple groups)."""
        # Do not change selection during drag (toggling on release is allowed)
        if getattr(self, "is_dragging_group", False):
            return

        # Search connected component via BFS/DFS
        visited = set()
        queue = [atom_idx]
        visited.add(atom_idx)

        while queue:
            current_idx = queue.pop(0)
            for bond_idx in range(self.mol.GetNumBonds()):
                bond = self.mol.GetBondWithIdx(bond_idx)
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()

                if begin_idx == current_idx and end_idx not in visited:
                    visited.add(end_idx)
                    queue.append(end_idx)
                elif end_idx == current_idx and begin_idx not in visited:
                    visited.add(begin_idx)
                    queue.append(begin_idx)

        # Add or remove as a new group
        if visited.issubset(self.group_atoms):
            # Already selected - remove
            self.group_atoms -= visited
        else:
            # Add new group
            self.group_atoms |= visited

        self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def update_display(self):
        if not self.group_atoms:
            self.selection_label.setText("No group selected")
        else:
            atom_info = []
            for atom_idx in sorted(self.group_atoms):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_info.append(f"{symbol}({atom_idx})")

            self.selection_label.setText(
                f"Selected group: {len(self.group_atoms)} atoms - {', '.join(atom_info[:5])}{' ...' if len(atom_info) > 5 else ''}"
            )

    def show_atom_labels(self):
        """Highlight atoms in the selected group (same style as Ctrl+Click)."""
        self.clear_atom_labels()

        if not self.group_atoms:
            return

        # Create list of selected atom indices
        selected_indices = list(self.group_atoms)

        # Get positions of selected atoms
        selected_positions = self.main_window.atom_positions_3d[selected_indices]

        # Highlight atoms with slightly larger radii
        selected_radii = np.array(
            [
                VDW_RADII.get(self.mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3
                for i in selected_indices
            ]
        )

        # Create dataset for highlighting
        highlight_source = pv.PolyData(selected_positions)
        highlight_source["radii"] = selected_radii

        # Highlight with semi-transparent yellow spheres
        highlight_glyphs = highlight_source.glyph(
            scale="radii",
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
            orient=False,
        )

        # Add and save highlight actor (set as non-pickable)
        self.highlight_actor = self.main_window.plotter.add_mesh(
            highlight_glyphs,
            color="yellow",
            opacity=0.3,
            name="move_group_highlight",
            pickable=False,  # Disable picking
        )

        self.main_window.plotter.render()

    def clear_atom_labels(self):
        """Clear atom highlights (including MoveGroup specific highlights)."""
        super().clear_atom_labels()
        try:
            self.main_window.plotter.remove_actor("move_group_highlight")
        except (AttributeError, RuntimeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        if hasattr(self, "highlight_actor"):
            try:
                self.main_window.plotter.remove_actor(self.highlight_actor)
            except (AttributeError, RuntimeError):  # pragma: no cover
                import traceback
                traceback.print_exc()

            self.highlight_actor = None
        try:
            self.main_window.plotter.render()
        except (AttributeError, RuntimeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

    def reset_translation_inputs(self):
        """Reset Translation input fields."""
        self.x_trans_input.setText("0.0")
        self.y_trans_input.setText("0.0")
        self.z_trans_input.setText("0.0")

    def apply_translation(self):
        """Translate the selected group."""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return

        try:
            dx = float(self.x_trans_input.text())
            dy = float(self.y_trans_input.text())
            dz = float(self.z_trans_input.text())
        except ValueError:
            QMessageBox.warning(
                self, "Warning", "Please enter valid translation values."
            )
            return

        translation_vector = np.array([dx, dy, dz])

        conf = self.mol.GetConformer()
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            new_pos = atom_pos + translation_vector
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos

        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # Redraw labels
        self.main_window.push_undo_state()

    def reset_rotation_inputs(self):
        """Reset Rotation input fields."""
        self.x_rot_input.setText("0.0")
        self.y_rot_input.setText("0.0")
        self.z_rot_input.setText("0.0")

    def apply_rotation(self):
        """Rotate the selected group."""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return

        try:
            rx = float(self.x_rot_input.text())
            ry = float(self.y_rot_input.text())
            rz = float(self.z_rot_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid rotation values.")
            return

        # Convert degrees to radians
        rx_rad = np.radians(rx)
        ry_rad = np.radians(ry)
        rz_rad = np.radians(rz)

        # Calculate group centroid
        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.group_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        centroid = np.mean(positions, axis=0)

        # Create rotation matrices
        # Around X-axis
        Rx = np.array(
            [
                [1, 0, 0],
                [0, np.cos(rx_rad), -np.sin(rx_rad)],
                [0, np.sin(rx_rad), np.cos(rx_rad)],
            ]
        )
        # Around Y-axis
        Ry = np.array(
            [
                [np.cos(ry_rad), 0, np.sin(ry_rad)],
                [0, 1, 0],
                [-np.sin(ry_rad), 0, np.cos(ry_rad)],
            ]
        )
        # Around Z-axis
        Rz = np.array(
            [
                [np.cos(rz_rad), -np.sin(rz_rad), 0],
                [np.sin(rz_rad), np.cos(rz_rad), 0],
                [0, 0, 1],
            ]
        )

        # Composite rotation matrix (Z * Y * X)
        R = Rz @ Ry @ Rx

        # Rotate each atom
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            # Move centroid to origin
            centered_pos = atom_pos - centroid
            # Rotate
            rotated_pos = R @ centered_pos
            # Restore centroid
            new_pos = rotated_pos + centroid
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos

        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # Redraw labels
        self.main_window.push_undo_state()

    def clear_selection(self):
        """Clear selection."""
        self.selected_atoms.clear()
        self.group_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        # Reset drag-related flags
        self.is_dragging_group = False
        self.drag_start_pos = None
        if hasattr(self, "last_drag_positions"):
            delattr(self, "last_drag_positions")

    def closeEvent(self, event):
        """Handle dialog close event."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """Handle cancel action."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
