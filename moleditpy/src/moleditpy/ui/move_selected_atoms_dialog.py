#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Optional, Any
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
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

try:
    from .atom_picking import pick_atom_index_from_screen
    from .base_picking_dialog import BasePickingDialog
    from ..utils.constants import VDW_RADII
except ImportError:
    from moleditpy.ui.atom_picking import pick_atom_index_from_screen
    from moleditpy.ui.base_picking_dialog import BasePickingDialog
    from moleditpy.utils.constants import VDW_RADII


class MoveSelectedAtomsDialog(BasePickingDialog):
    """Dialog to select specific atoms and perform translation/rotation on them."""

    def __init__(
        self,
        mol: Any,
        main_window: Any,
        preselected_atoms: Any = None,
        parent: Any = None,
    ) -> None:
        super().__init__(mol, main_window, parent)
        self.selected_atoms: set[int] = set()

        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)

        self.clicked_atom_for_toggle: Optional[int] = None
        # State for group movement (used by dialog's own event filter)
        self.drag_atom_idx: Optional[int] = None
        self.potential_drag: bool = False
        self.is_dragging_group: bool = False
        self.drag_start_pos: Optional[Any] = None
        self.mouse_moved_during_drag: bool = False
        self._consume_next_left_release: bool = False
        self.highlight_actor: Optional[pv.Actor] = None

        self.widgets: dict[str, Any] = {}

        self.init_ui()
        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    @property
    def group_atoms(self) -> set[int]:
        """Expose selected_atoms as group_atoms for CustomInteractorStyle compatibility."""
        return self.selected_atoms

    @group_atoms.setter
    def group_atoms(self, value: set[int]) -> None:
        self.selected_atoms = value

    @property
    def x_trans_input(self) -> QLineEdit:
        """Expose x_trans_input widget."""
        return self.widgets["x_trans_input"]

    @property
    def y_trans_input(self) -> QLineEdit:
        """Expose y_trans_input widget."""
        return self.widgets["y_trans_input"]

    @property
    def z_trans_input(self) -> QLineEdit:
        """Expose z_trans_input widget."""
        return self.widgets["z_trans_input"]

    @property
    def x_rot_input(self) -> QLineEdit:
        """Expose x_rot_input widget."""
        return self.widgets["x_rot_input"]

    @property
    def y_rot_input(self) -> QLineEdit:
        """Expose y_rot_input widget."""
        return self.widgets["y_rot_input"]

    @property
    def z_rot_input(self) -> QLineEdit:
        """Expose z_rot_input widget."""
        return self.widgets["z_rot_input"]

    def init_ui(self) -> None:
        """Initialize UI widgets and layout."""
        self.setWindowTitle("Move Selected Atoms")
        self.setModal(False)
        self.resize(300, 400)
        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click atoms in the 3D view to select/deselect them.\n"
            "Left-drag: Move the selected atoms\n"
            "Right-drag: Rotate the selected atoms around their center"
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.widgets["selection_label"] = QLabel("No atoms selected")
        layout.addWidget(self.widgets["selection_label"])

        self._init_translation_ui(layout)
        self._init_rotation_ui(layout)
        self._init_buttons_ui(layout)

        # Enable picking to handle atom selection
        self.enable_picking()

    def _init_translation_ui(self, layout: QVBoxLayout) -> None:
        """Initialize translation widgets."""
        trans_group = QLabel("Translation (Å):")
        trans_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(trans_group)

        trans_layout = QGridLayout()
        x_in = QLineEdit("0.0")
        y_in = QLineEdit("0.0")
        z_in = QLineEdit("0.0")

        self.widgets["x_trans_input"] = x_in
        self.widgets["y_trans_input"] = y_in
        self.widgets["z_trans_input"] = z_in

        # Execute apply_translation on Enter key
        x_in.returnPressed.connect(self.apply_translation)
        y_in.returnPressed.connect(self.apply_translation)
        z_in.returnPressed.connect(self.apply_translation)

        trans_layout.addWidget(QLabel("X:"), 0, 0)
        trans_layout.addWidget(x_in, 0, 1)
        trans_layout.addWidget(QLabel("Y:"), 1, 0)
        trans_layout.addWidget(y_in, 1, 1)
        trans_layout.addWidget(QLabel("Z:"), 2, 0)
        trans_layout.addWidget(z_in, 2, 1)

        trans_button_layout = QHBoxLayout()
        reset_trans_button = QPushButton("Reset")
        reset_trans_button.clicked.connect(self.reset_translation_inputs)
        trans_button_layout.addWidget(reset_trans_button)

        apply_trans_button = QPushButton("Apply Translation")
        apply_trans_button.clicked.connect(self.apply_translation)
        trans_button_layout.addWidget(apply_trans_button)

        trans_layout.addLayout(trans_button_layout, 3, 0, 1, 2)
        layout.addLayout(trans_layout)

    def _init_rotation_ui(self, layout: QVBoxLayout) -> None:
        """Initialize rotation widgets."""
        rot_group = QLabel("Rotation (degrees):")
        rot_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(rot_group)

        rot_layout = QGridLayout()
        x_rot = QLineEdit("0.0")
        y_rot = QLineEdit("0.0")
        z_rot = QLineEdit("0.0")

        self.widgets["x_rot_input"] = x_rot
        self.widgets["y_rot_input"] = y_rot
        self.widgets["z_rot_input"] = z_rot

        # Execute apply_rotation on Enter key
        x_rot.returnPressed.connect(self.apply_rotation)
        y_rot.returnPressed.connect(self.apply_rotation)
        z_rot.returnPressed.connect(self.apply_rotation)

        rot_layout.addWidget(QLabel("Around X:"), 0, 0)
        rot_layout.addWidget(x_rot, 0, 1)
        rot_layout.addWidget(QLabel("Around Y:"), 1, 0)
        rot_layout.addWidget(y_rot, 1, 1)
        rot_layout.addWidget(QLabel("Around Z:"), 2, 0)
        rot_layout.addWidget(z_rot, 2, 1)

        rot_button_layout = QHBoxLayout()
        reset_rot_button = QPushButton("Reset")
        reset_rot_button.clicked.connect(self.reset_rotation_inputs)
        rot_button_layout.addWidget(reset_rot_button)

        apply_rot_button = QPushButton("Apply Rotation")
        apply_rot_button.clicked.connect(self.apply_rotation)
        rot_button_layout.addWidget(apply_rot_button)

        rot_layout.addLayout(rot_button_layout, 3, 0, 1, 2)
        layout.addLayout(rot_layout)

    def _init_buttons_ui(self, layout: QVBoxLayout) -> None:
        """Initialize bottom buttons."""
        button_layout = QHBoxLayout()
        clear_btn = QPushButton("Clear Selection")
        clear_btn.clicked.connect(self.clear_selection)
        self.widgets["clear_button"] = clear_btn
        button_layout.addWidget(clear_btn)

        button_layout.addStretch()

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

    def eventFilter(self, obj: Any, event: Any) -> bool:
        """Mouse event handling in 3D view.

        Delegates to CustomInteractorStyle if atoms are selected.
        """
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None or self.mol is None:
            return False

        if obj != plotter.interactor:
            return super().eventFilter(obj, event)

        e_type = event.type()

        if e_type == QEvent.Type.MouseButtonDblClick:
            # Ignore double clicks and reset state
            self.is_dragging_group = False
            self.drag_start_pos = None
            self.mouse_moved_during_drag = False
            self.potential_drag = False
            self.clicked_atom_for_toggle = None
            return False

        if (
            e_type == QEvent.Type.MouseButtonPress
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            # Clean up previous state (triple-click countermeasure)
            self.is_dragging_group = False
            self.potential_drag = False
            self.clicked_atom_for_toggle = None
            # Delegate to CustomInteractorStyle if atoms are already selected
            if self.selected_atoms:
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
                    if self.selected_atoms and clicked_atom_idx in self.selected_atoms:
                        # Atom within existing group - prepare for drag
                        self.is_dragging_group = False
                        self.drag_start_pos = click_pos
                        self.drag_atom_idx = clicked_atom_idx
                        self.mouse_moved_during_drag = False
                        self.potential_drag = True
                        self.clicked_atom_for_toggle = clicked_atom_idx
                        return False
                    else:
                        # Atom outside group - select
                        self.on_atom_picked(clicked_atom_idx)
                        self._consume_next_left_release = True
                        return True
                else:
                    # Clicked outside atoms
                    return False

            except (AttributeError, RuntimeError, ValueError) as e:
                logging.debug(f"Error in mouse press: {e}")
                return False

        elif e_type == QEvent.Type.MouseMove and isinstance(event, QMouseEvent):
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
                                plotter_ptr.setCursor(Qt.CursorShape.ClosedHandCursor)
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
            if self.selected_atoms:
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

                    if closest_atom_idx in self.selected_atoms:
                        plotter_ref.setCursor(Qt.CursorShape.OpenHandCursor)
                    else:
                        plotter_ref.setCursor(Qt.CursorShape.ArrowCursor)
                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    logging.debug(f"Error updating hover cursor: {e}")

            return False

        elif (
            e_type == QEvent.Type.MouseButtonRelease
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
                    if not (self.is_dragging_group and self.mouse_moved_during_drag):
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
                                plotter_ptr = self.main_window.view_3d_manager.plotter
                                if plotter_ptr is not None:
                                    plotter_ptr.setCursor(Qt.CursorShape.ArrowCursor)
                            except (
                                AttributeError,
                                RuntimeError,
                                ValueError,
                                TypeError,
                            ) as e:
                                logging.debug(f"Failed to reset cursor to arrow: {e}")
                            return True
                        else:
                            logging.error(
                                "REPORT ERROR: Missing attribute 'clicked_atom_for_toggle' on self"
                            )

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
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Failed to reset cursor in release finally: {e}")

                return True

            return False

        return super().eventFilter(obj, event)

    def on_atom_picked(self, atom_idx: int) -> None:
        """Select or deselect the clicked atom."""
        if getattr(self, "is_dragging_group", False):
            return

        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)

        self.show_atom_labels()
        self.update_display()

    def update_display(self) -> None:
        """Update the selected atoms text label."""
        if not self.selected_atoms:
            self.widgets["selection_label"].setText("No atoms selected")
        else:
            atom_info = []
            for atom_idx in sorted(self.selected_atoms):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_info.append(f"{symbol}({atom_idx})")

            info_str = ", ".join(atom_info[:5])
            if len(atom_info) > 5:
                info_str += " ..."
            self.widgets["selection_label"].setText(
                f"Selected: {len(self.selected_atoms)} atoms - {info_str}"
            )

    def show_atom_labels(self) -> None:
        """Highlight selected atoms."""
        self.clear_atom_labels()

        if not self.selected_atoms:
            return

        selected_indices = list(self.selected_atoms)
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
            name="move_selected_atoms_highlight",
            pickable=False,
        )

        plotter.render()

    def clear_atom_labels(self) -> None:
        """Clear highlights."""
        super().clear_atom_labels()

        plotter = self.main_window.view_3d_manager.plotter
        if plotter is not None:
            try:
                plotter.remove_actor("move_selected_atoms_highlight")
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
        """Reset translation entry fields to 0.0."""
        self.widgets["x_trans_input"].setText("0.0")
        self.widgets["y_trans_input"].setText("0.0")
        self.widgets["z_trans_input"].setText("0.0")

    def apply_translation(self) -> None:
        """Translate the selected atoms."""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select atoms first.")
            return

        try:
            dx = float(self.widgets["x_trans_input"].text())
            dy = float(self.widgets["y_trans_input"].text())
            dz = float(self.widgets["z_trans_input"].text())
        except ValueError:
            QMessageBox.warning(
                self, "Warning", "Please enter valid translation values."
            )
            return

        translation_vector = np.array([dx, dy, dz])
        positions = self.mol.GetConformer().GetPositions()
        for atom_idx in self.selected_atoms:
            positions[atom_idx] += translation_vector

        self._update_molecule_geometry(positions)
        self._push_undo()
        self.show_atom_labels()

    def reset_rotation_inputs(self) -> None:
        """Reset rotation entry fields to 0.0."""
        self.widgets["x_rot_input"].setText("0.0")
        self.widgets["y_rot_input"].setText("0.0")
        self.widgets["z_rot_input"].setText("0.0")

    def apply_rotation(self) -> None:
        """Rotate the selected atoms around their centroid."""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select atoms first.")
            return

        try:
            rx_rad, ry_rad, rz_rad = np.radians(
                [
                    float(self.widgets["x_rot_input"].text()),
                    float(self.widgets["y_rot_input"].text()),
                    float(self.widgets["z_rot_input"].text()),
                ]
            )
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid rotation values.")
            return

        positions = self.mol.GetConformer().GetPositions()

        # Calculate centroid of the selected atoms
        selected_indices = list(self.selected_atoms)
        selected_positions = positions[selected_indices]
        centroid = np.mean(selected_positions, axis=0)

        # Rotation matrices
        r_x = np.array(
            [
                [1, 0, 0],
                [0, np.cos(rx_rad), -np.sin(rx_rad)],
                [0, np.sin(rx_rad), np.cos(rx_rad)],
            ]
        )
        r_y = np.array(
            [
                [np.cos(ry_rad), 0, np.sin(ry_rad)],
                [0, 1, 0],
                [-np.sin(ry_rad), 0, np.cos(ry_rad)],
            ]
        )
        r_z = np.array(
            [
                [np.cos(rz_rad), -np.sin(rz_rad), 0],
                [np.sin(rz_rad), np.cos(rz_rad), 0],
                [0, 0, 1],
            ]
        )
        rot_matrix = r_z @ r_y @ r_x

        for atom_idx in self.selected_atoms:
            pos = positions[atom_idx]
            new_pos = rot_matrix @ (pos - centroid) + centroid
            positions[atom_idx] = new_pos

        self._update_molecule_geometry(positions)
        self._push_undo()
        self.show_atom_labels()

    def clear_selection(self) -> None:
        """Clear selection."""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        self.is_dragging_group = False
        self.drag_start_pos = None
