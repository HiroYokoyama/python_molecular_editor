#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import TYPE_CHECKING, Optional, Union, Any
import logging
import numpy as np
import pyvista as pv
from PyQt6.QtCore import QEvent, Qt, QObject
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import (
    QApplication,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

try:
    from .base_picking_dialog import BasePickingDialog
    from ..utils.constants import VDW_RADII, pt
except ImportError:
    from moleditpy.ui.base_picking_dialog import BasePickingDialog
    from moleditpy.utils.constants import VDW_RADII, pt


class MoveGroupDialog(BasePickingDialog):
    """Dialog to select a connected molecular group and perform translation/rotation."""

    def __init__(
        self,
        mol: Any,
        main_window: Any,
        preselected_atoms: Any = None,
        parent: Any = None,
    ) -> None:
        super().__init__(mol, main_window, parent)
        self.selected_atoms: set[int] = set()
        self.group_atoms: set[int] = set()  # All atoms connected to selected atoms

        # Add preselected atoms
        if preselected_atoms:
            # For MoveGroup, we pick the first atom and select its connected group
            self.on_atom_picked(preselected_atoms[0])

        self.clicked_atom_for_toggle: Optional[int] = None
        # State for group movement (used by CustomInteractorStyle)
        self._initial_positions: dict = {}
        self._is_dragging_group_vtk = False
        self._is_rotating_group_vtk = False
        self._drag_atom_idx: Optional[int] = None
        self._drag_start_pos: Optional[Any] = None
        self._mouse_moved: bool = False
        self._rotation_start_pos: Optional[Any] = None
        self._rotation_mouse_moved: bool = False
        self._rotation_atom_idx: Optional[int] = None
        self._group_centroid: Optional[np.ndarray] = None

        # State for group movement (used by dialog's own event filter)
        self.drag_atom_idx: Optional[int] = None
        self.potential_drag: bool = False
        self.is_dragging_group: bool = False
        self.drag_start_pos: Optional[Any] = None
        self.mouse_moved_during_drag: bool = False
        self.highlight_actor: Optional[pv.Actor] = None

        self.init_ui()

    def init_ui(self) -> None:
        self.setWindowTitle("Move Group")
        self.setModal(False)
        self.resize(300, 400)
        layout = QVBoxLayout(self)

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

    def eventFilter(self, obj: Any, event: Any) -> bool:
        """Mouse event handling in 3D view - delegate to CustomInteractorStyle if a group is selected."""
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None or self.mol is None:
            return False

        if obj == plotter.interactor:
            # Prevent state confusion from double/triple clicks
            if event.type() == QEvent.Type.MouseButtonDblClick:
                # Ignore double clicks and reset state
                self.is_dragging_group = False
                self.drag_start_pos = None
                self.mouse_moved_during_drag = False
                self.potential_drag = False
                self.clicked_atom_for_toggle = None
                return False

            if (
                event.type() == QEvent.Type.MouseButtonPress
                and isinstance(event, QMouseEvent)
                and event.button() == Qt.MouseButton.LeftButton
            ):
                # Clean up previous state (triple-click countermeasure)
                self.is_dragging_group = False
                self.potential_drag = False
                self.clicked_atom_for_toggle = None
                # Delegate to CustomInteractorStyle if a group is already selected
                if self.group_atoms:
                    return False

                # Mouse press handling
                try:
                    interactor = plotter.interactor
                    if interactor is None:
                        return False
                    click_pos = interactor.GetEventPosition()

                    # Pick via plotter
                    picker = plotter.picker
                    if picker is None:
                        return False
                    picker.Pick(
                        click_pos[0],
                        click_pos[1],
                        0,
                        plotter.renderer,
                    )

                    clicked_atom_idx = None
                    if picker.GetActor() is self.main_window.view_3d_manager.atom_actor:
                        picked_position = np.array(picker.GetPickPosition())
                        if self.main_window.view_3d_manager.atom_positions_3d is None:
                            return False
                        distances = np.linalg.norm(
                            self.main_window.view_3d_manager.atom_positions_3d
                            - picked_position,
                            axis=1,
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
                                except (
                                    AttributeError,
                                    RuntimeError,
                                    ValueError,
                                    TypeError,
                                ):
                                    vdw_radius = 1.5
                                click_threshold = self._get_click_threshold(vdw_radius)

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
                            self.potential_drag = True
                            self.clicked_atom_for_toggle = clicked_atom_idx
                            return False
                        else:
                            # Atom outside group - select new group
                            self.on_atom_picked(clicked_atom_idx)
                            return True
                    else:
                        # Clicked outside atoms
                        return False

                except (AttributeError, RuntimeError, ValueError) as e:
                    logging.debug(f"Error in mouse press: {e}")
                    return False

            elif event.type() == QEvent.Type.MouseMove and isinstance(
                event, QMouseEvent
            ):
                # Mouse move handling
                if (
                    getattr(self, "potential_drag", False)
                    and self.drag_start_pos
                    and not self.is_dragging_group
                ):
                    try:
                        plotter_ref = self.main_window.view_3d_manager.plotter
                        if plotter_ref is None or plotter_ref.interactor is None:
                            return False
                        interactor = plotter_ref.interactor
                        current_pos = interactor.GetEventPosition()
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]

                        # Start drag if threshold is exceeded
                        drag_threshold = 5  # pixels
                        if abs(dx) > drag_threshold or abs(dy) > drag_threshold:
                            self.is_dragging_group = True
                            self.potential_drag = False
                            try:
                                plotter_ptr = self.main_window.view_3d_manager.plotter
                                if plotter_ptr is not None:
                                    plotter_ptr.setCursor(
                                        Qt.CursorShape.ClosedHandCursor
                                    )
                            except (
                                AttributeError,
                                RuntimeError,
                                ValueError,
                                TypeError,
                            ):
                                pass
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass

                    if not self.is_dragging_group:
                        return False

                if self.is_dragging_group and self.drag_start_pos:
                    try:
                        plotter_ref = self.main_window.view_3d_manager.plotter
                        if plotter_ref is None or plotter_ref.interactor is None:
                            return False
                        interactor = plotter_ref.interactor
                        current_pos = interactor.GetEventPosition()
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]
                        if abs(dx) > 2 or abs(dy) > 2:
                            self.mouse_moved_during_drag = True
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass
                    return True

                # Hover handling
                if self.group_atoms:
                    try:
                        plotter_ref = self.main_window.view_3d_manager.plotter
                        if plotter_ref is None or plotter_ref.interactor is None:
                            return False
                        interactor = plotter_ref.interactor
                        current_pos = interactor.GetEventPosition()
                        picker = plotter_ref.picker
                        if picker is None:
                            return False
                        picker.Pick(
                            current_pos[0],
                            current_pos[1],
                            0,
                            plotter_ref.renderer,
                        )

                        if (
                            picker.GetActor()
                            is self.main_window.view_3d_manager.atom_actor
                        ):
                            picked_position = np.array(picker.GetPickPosition())
                            if (
                                self.main_window.view_3d_manager.atom_positions_3d
                                is None
                            ):
                                return False
                            distances = np.linalg.norm(
                                self.main_window.view_3d_manager.atom_positions_3d
                                - picked_position,
                                axis=1,
                            )
                            closest_atom_idx = np.argmin(distances)

                            if closest_atom_idx in self.group_atoms:
                                plotter_ref.setCursor(Qt.CursorShape.OpenHandCursor)
                            else:
                                plotter_ref.setCursor(Qt.CursorShape.ArrowCursor)
                        else:
                            plotter_ref.setCursor(Qt.CursorShape.ArrowCursor)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass

                return False

            elif (
                event.type() == QEvent.Type.MouseButtonRelease
                and isinstance(event, QMouseEvent)
                and event.button() == Qt.MouseButton.LeftButton
            ):
                if getattr(self, "potential_drag", False) or (
                    self.is_dragging_group and self.drag_start_pos
                ):
                    try:
                        if not (
                            self.is_dragging_group and self.mouse_moved_during_drag
                        ):
                            # Mouse move below threshold = simple click (toggle)
                            if self.clicked_atom_for_toggle is not None:
                                clicked_atom = self.clicked_atom_for_toggle
                                self.clicked_atom_for_toggle = None
                                self.is_dragging_group = False
                                self.drag_start_pos = None
                                self.mouse_moved_during_drag = False
                                self.potential_drag = False
                                if clicked_atom is not None:
                                    self.on_atom_picked(clicked_atom)
                                try:
                                    plotter_ptr = (
                                        self.main_window.view_3d_manager.plotter
                                    )
                                    if plotter_ptr is not None:
                                        plotter_ptr.setCursor(
                                            Qt.CursorShape.ArrowCursor
                                        )
                                except (
                                    AttributeError,
                                    RuntimeError,
                                    ValueError,
                                    TypeError,
                                ):
                                    pass
                                return True
                            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                                logging.error(
                                    "REPORT ERROR: Missing attribute 'clicked_atom_for_toggle' on self"
                                )

                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass
                    finally:
                        self.is_dragging_group = False
                        self.drag_start_pos = None
                        self.mouse_moved_during_drag = False
                        self.potential_drag = False
                        try:
                            plotter_ptr = self.main_window.view_3d_manager.plotter
                            if plotter_ptr is not None:
                                plotter_ptr.setCursor(Qt.CursorShape.ArrowCursor)
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            pass

                    return True

                return False

        return super().eventFilter(obj, event)

    def on_atom_picked(self, atom_idx: int) -> None:
        """Select the entire connected component the atom belongs to."""
        if getattr(self, "is_dragging_group", False):
            return

        # BFS for connected atoms
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

        # Toggle group
        if visited.issubset(self.group_atoms):
            self.group_atoms -= visited
        else:
            self.group_atoms |= visited

        self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def update_display(self) -> None:
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

    def show_atom_labels(self) -> None:
        """Highlight atoms in the selected group."""
        self.clear_atom_labels()

        if not self.group_atoms:
            return

        selected_indices = list(self.group_atoms)
        plotter = self.main_window.view_3d_manager.plotter
        if self.main_window.view_3d_manager.atom_positions_3d is None:
            logging.error("atom_positions_3d is None in update_atom_labels")
            return
        selected_positions = self.main_window.view_3d_manager.atom_positions_3d[
            selected_indices
        ]
        selected_radii = np.array(
            [
                VDW_RADII.get(self.mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3
                for i in selected_indices
            ]
        )

        highlight_source = pv.PolyData(selected_positions)
        highlight_source["radii"] = selected_radii
        highlight_glyphs = highlight_source.glyph(
            scale="radii",
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
            orient=False,
        )

        if plotter is None:
            return
        self.highlight_actor = plotter.add_mesh(
            highlight_glyphs,
            color="yellow",
            opacity=0.3,
            name="move_group_highlight",
            pickable=False,
        )

        plotter.render()

    def clear_atom_labels(self) -> None:
        """Clear highlights."""
        # Call base which clears selection_labels (standard labels)
        super().clear_atom_labels()

        # Clear MoveGroup specific highlight
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is not None:
            try:
                plotter.remove_actor("move_group_highlight")
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass

        if self.highlight_actor:
            if plotter is not None:
                try:
                    plotter.remove_actor(self.highlight_actor)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    pass
            self.highlight_actor = None

        if plotter is not None:
            try:
                plotter.render()
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass

    def reset_translation_inputs(self) -> None:
        self.x_trans_input.setText("0.0")
        self.y_trans_input.setText("0.0")
        self.z_trans_input.setText("0.0")

    def apply_translation(self) -> None:
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
        positions = self.mol.GetConformer().GetPositions()
        for atom_idx in self.group_atoms:
            positions[atom_idx] += translation_vector

        # Write updated positions back using inherited helper
        self._update_molecule_geometry(positions)

        # Push Undo state AFTER modification
        self._push_undo()
        self.show_atom_labels()

    def reset_rotation_inputs(self) -> None:
        self.x_rot_input.setText("0.0")
        self.y_rot_input.setText("0.0")
        self.z_rot_input.setText("0.0")

    def apply_rotation(self) -> None:
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

        rx_rad, ry_rad, rz_rad = np.radians([rx, ry, rz])
        positions = self.mol.GetConformer().GetPositions()

        # Calculate centroid of the group
        group_indices = list(self.group_atoms)
        group_positions = positions[group_indices]
        centroid = np.mean(group_positions, axis=0)

        # Rotation matrices
        Rx = np.array(
            [
                [1, 0, 0],
                [0, np.cos(rx_rad), -np.sin(rx_rad)],
                [0, np.sin(rx_rad), np.cos(rx_rad)],
            ]
        )
        Ry = np.array(
            [
                [np.cos(ry_rad), 0, np.sin(ry_rad)],
                [0, 1, 0],
                [-np.sin(ry_rad), 0, np.cos(ry_rad)],
            ]
        )
        Rz = np.array(
            [
                [np.cos(rz_rad), -np.sin(rz_rad), 0],
                [np.sin(rz_rad), np.cos(rz_rad), 0],
                [0, 0, 1],
            ]
        )
        R = Rz @ Ry @ Rx

        for atom_idx in self.group_atoms:
            pos = positions[atom_idx]
            new_pos = R @ (pos - centroid) + centroid
            positions[atom_idx] = new_pos

        # Write updated positions back using inherited helper
        self._update_molecule_geometry(positions)

        # Push Undo state AFTER modification
        self._push_undo()
        self.show_atom_labels()

    def clear_selection(self) -> None:
        """Clear selection."""
        self.selected_atoms.clear()
        self.group_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        self.is_dragging_group = False
        self.drag_start_pos = None
