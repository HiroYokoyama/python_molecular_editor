#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Any, Dict, Optional
import logging
import numpy as np
import pyvista as pv
from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

from .atom_picking import pick_atom_index_from_screen
from .base_picking_dialog import BasePickingDialog
from .move_dialog_mixin import MoveDialogMixin


class MoveGroupDialog(MoveDialogMixin, BasePickingDialog):
    """Dialog to select a connected molecular group and perform translation/rotation."""

    def __init__(
        self,
        mol: Any,
        main_window: Any,
        preselected_atoms: Any = None,
        parent: Any = None,
    ) -> None:
        """Initialize move group transformation dialog."""
        super().__init__(mol, main_window, parent)
        self.clear_button: Any = None
        self.selection_label: Any = None
        self.x_rot_input: Any = None
        self.x_trans_input: Any = None
        self.y_rot_input: Any = None
        self.y_trans_input: Any = None
        self.z_rot_input: Any = None
        self.z_trans_input: Any = None
        self.selected_atoms: set[int] = set()
        self.group_atoms: set[int] = set()  # All atoms connected to selected atoms

        self.clicked_atom_for_toggle: Optional[int] = None
        # State for group movement (used by CustomInteractorStyle)
        self.initial_positions: Dict[int, np.ndarray] = {}
        self.is_dragging_group_vtk = False
        self.is_rotating_group_vtk = False
        self.drag_atom_idx_vtk: Optional[int] = None
        self.drag_start_pos_vtk: Optional[Any] = None
        self.mouse_moved_vtk: bool = False
        self.rotation_start_pos: Optional[Any] = None
        self.rotation_mouse_moved: bool = False
        self.rotation_atom_idx: Optional[int] = None
        self.group_centroid: Optional[np.ndarray] = None

        # State for group movement (used by dialog's own event filter)
        self.drag_atom_idx: Optional[int] = None
        self.potential_drag: bool = False
        self.is_dragging_group: bool = False
        self.drag_start_pos: Optional[Any] = None
        self.mouse_moved_during_drag: bool = False
        self._consume_next_left_release = False
        self.highlight_actor: Optional[pv.Actor] = None

        self.init_ui()

        # After init_ui: picking an atom highlights it and updates the labels,
        # which needs both the actor state and the widgets to exist.
        if preselected_atoms:
            self.on_atom_picked(preselected_atoms[0])

    def init_ui(self) -> None:
        """Build the move-group dialog with atom picker, translate/rotate inputs, and controls."""
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

                    clicked_atom_idx = pick_atom_index_from_screen(
                        self.main_window.view_3d_manager,
                        (int(click_pos[0]), int(click_pos[1])),
                        self.mol,
                    )

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
                            self._consume_next_left_release = True
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
                            ) as e:
                                logging.debug(f"Failed to set closed hand cursor: {e}")
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Error initiating drag on move: {e}")

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
                        if abs(dx) > 5 or abs(dy) > 5:
                            self.mouse_moved_during_drag = True
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Error tracking drag movement: {e}")
                    return True

                # Hover handling
                if self.group_atoms:
                    try:
                        plotter_ref = self.main_window.view_3d_manager.plotter
                        if plotter_ref is None or plotter_ref.interactor is None:
                            return False
                        interactor = plotter_ref.interactor
                        current_pos = interactor.GetEventPosition()
                        closest_atom_idx = pick_atom_index_from_screen(
                            self.main_window.view_3d_manager,
                            (int(current_pos[0]), int(current_pos[1])),
                            self.mol,
                        )

                        if closest_atom_idx in self.group_atoms:
                            plotter_ref.setCursor(Qt.CursorShape.OpenHandCursor)
                        else:
                            plotter_ref.setCursor(Qt.CursorShape.ArrowCursor)
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Error updating hover cursor: {e}")

                return False

            elif (
                event.type() == QEvent.Type.MouseButtonRelease
                and isinstance(event, QMouseEvent)
                and event.button() == Qt.MouseButton.LeftButton
            ):
                if self._consume_next_left_release:
                    self._consume_next_left_release = False
                    return True

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
                                ) as e:
                                    logging.debug(
                                        f"Failed to reset cursor to arrow: {e}"
                                    )
                                return True

                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Error in mouse release handling: {e}")
                    finally:
                        self.is_dragging_group = False
                        self.drag_start_pos = None
                        self.mouse_moved_during_drag = False
                        self.potential_drag = False
                        try:
                            plotter_ptr = self.main_window.view_3d_manager.plotter
                            if plotter_ptr is not None:
                                plotter_ptr.setCursor(Qt.CursorShape.ArrowCursor)
                        except (
                            AttributeError,
                            RuntimeError,
                            ValueError,
                            TypeError,
                        ) as e:
                            logging.debug(
                                f"Failed to reset cursor in release finally: {e}"
                            )

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
            if atom_idx in self.selected_atoms:
                self.selected_atoms.remove(atom_idx)
        else:
            self.group_atoms |= visited
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()
