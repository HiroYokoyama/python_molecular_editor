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
    from moleditpy_linux.ui.atom_picking import pick_atom_index_from_screen
    from moleditpy_linux.ui.base_picking_dialog import BasePickingDialog
    from moleditpy_linux.utils.constants import VDW_RADII


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
        self.highlight_actor: Optional[pv.Actor] = None

        # Grouped states to comply with Pylint instance attribute limit
        self.drag_state: dict[str, Any] = {
            "drag_atom_idx": None,
            "potential_drag": False,
            "is_dragging_group": False,
            "drag_start_pos": (0, 0),
            "mouse_moved_during_drag": False,
            "consume_next_left_release": False,
        }
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

    def _handle_double_click(self) -> bool:
        """Handle MouseButtonDblClick event."""
        self.drag_state["is_dragging_group"] = False
        self.drag_state["drag_start_pos"] = (0, 0)
        self.drag_state["mouse_moved_during_drag"] = False
        self.drag_state["potential_drag"] = False
        self.clicked_atom_for_toggle = None
        return False

    def _handle_mouse_press(self, plotter: Any) -> bool:
        """Handle MouseButtonPress event for LeftButton."""
        self.drag_state["is_dragging_group"] = False
        self.drag_state["potential_drag"] = False
        self.clicked_atom_for_toggle = None

        if self.selected_atoms:
            return False

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

            if clicked_atom_idx is not None:
                if self.selected_atoms and clicked_atom_idx in self.selected_atoms:
                    self.drag_state["is_dragging_group"] = False
                    self.drag_state["drag_start_pos"] = click_pos
                    self.drag_state["drag_atom_idx"] = clicked_atom_idx
                    self.drag_state["mouse_moved_during_drag"] = False
                    self.drag_state["potential_drag"] = True
                    self.clicked_atom_for_toggle = clicked_atom_idx
                    return False

                self.on_atom_picked(clicked_atom_idx)
                self.drag_state["consume_next_left_release"] = True
                return True

            return False

        except (AttributeError, RuntimeError, ValueError) as e:
            logging.debug("Error in mouse press: %s", e)
            return False

    def _handle_potential_drag(self, plotter: Any) -> bool:
        """Handle potential drag checking and threshold transition."""
        start_pos = self.drag_state["drag_start_pos"]
        try:
            interactor = plotter.interactor
            if interactor is None:
                return False
            current_pos = interactor.GetEventPosition()
            dx = current_pos[0] - start_pos[0]
            dy = current_pos[1] - start_pos[1]

            if abs(dx) > 5 or abs(dy) > 5:
                self.drag_state["is_dragging_group"] = True
                self.drag_state["potential_drag"] = False
                try:
                    plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    pass
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass
        return self.drag_state["is_dragging_group"]

    def _handle_actual_drag(self, plotter: Any) -> bool:
        """Update drag movement flags during active dragging."""
        start_pos = self.drag_state["drag_start_pos"]
        try:
            interactor = plotter.interactor
            if interactor is None:
                return False
            current_pos = interactor.GetEventPosition()
            dx = current_pos[0] - start_pos[0]
            dy = current_pos[1] - start_pos[1]
            if abs(dx) > 2 or abs(dy) > 2:
                self.drag_state["mouse_moved_during_drag"] = True
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass
        return True

    def _handle_hover_cursor(self, plotter: Any) -> bool:
        """Set appropriate hover cursor style over selected atoms."""
        try:
            interactor = plotter.interactor
            if interactor is None:
                return False
            current_pos = interactor.GetEventPosition()
            closest_atom_idx = pick_atom_index_from_screen(
                self.main_window.view_3d_manager,
                (int(current_pos[0]), int(current_pos[1])),
                self.mol,
            )

            if closest_atom_idx in self.selected_atoms:
                plotter.setCursor(Qt.CursorShape.OpenHandCursor)
            else:
                plotter.setCursor(Qt.CursorShape.ArrowCursor)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass
        return False

    def _handle_mouse_move(self, plotter: Any) -> bool:
        """Handle MouseMove event."""
        if (
            self.drag_state["potential_drag"]
            and not self.drag_state["is_dragging_group"]
        ):
            if not self._handle_potential_drag(plotter):
                return False

        if self.drag_state["is_dragging_group"]:
            return self._handle_actual_drag(plotter)

        if self.selected_atoms:
            return self._handle_hover_cursor(plotter)

        return False

    def _reset_drag_state(self) -> None:
        """Reset internal drag flags."""
        self.drag_state["is_dragging_group"] = False
        self.drag_state["drag_start_pos"] = (0, 0)
        self.drag_state["mouse_moved_during_drag"] = False
        self.drag_state["potential_drag"] = False

    def _restore_arrow_cursor(self) -> None:
        """Restore standard arrow cursor in the 3D viewer."""
        try:
            plotter_ptr = self.main_window.view_3d_manager.plotter
            if plotter_ptr is not None:
                plotter_ptr.setCursor(Qt.CursorShape.ArrowCursor)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            pass

    def _handle_mouse_release(self) -> bool:
        """Handle MouseButtonRelease event for LeftButton."""
        if self.drag_state["consume_next_left_release"]:
            self.drag_state["consume_next_left_release"] = False
            return True

        is_drag_active = (
            self.drag_state["potential_drag"] or self.drag_state["is_dragging_group"]
        )
        if not is_drag_active:
            return False

        has_moved = (
            self.drag_state["is_dragging_group"]
            and self.drag_state["mouse_moved_during_drag"]
        )
        clicked_atom = self.clicked_atom_for_toggle

        if not has_moved and clicked_atom is not None:
            self.clicked_atom_for_toggle = None
            self._reset_drag_state()
            try:
                self.on_atom_picked(clicked_atom)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                pass
            self._restore_arrow_cursor()
            return True

        self._reset_drag_state()
        self._restore_arrow_cursor()
        return True

    def eventFilter(self, obj: Any, event: Any) -> bool:
        """Mouse event handling in 3D view.

        Delegates to CustomInteractorStyle if atoms are selected.
        """
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None or self.mol is None:
            return False

        if obj != plotter.interactor:
            return super().eventFilter(obj, event)

        result = False
        e_type = event.type()

        if e_type == QEvent.Type.MouseButtonDblClick:
            result = self._handle_double_click()
        elif (
            e_type == QEvent.Type.MouseButtonPress
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            result = self._handle_mouse_press(plotter)
        elif e_type == QEvent.Type.MouseMove and isinstance(event, QMouseEvent):
            result = self._handle_mouse_move(plotter)
        elif (
            e_type == QEvent.Type.MouseButtonRelease
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            result = self._handle_mouse_release()

        return result

    def on_atom_picked(self, atom_idx: int) -> None:
        """Select or deselect the clicked atom."""
        if self.drag_state["is_dragging_group"]:
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
        self.drag_state["is_dragging_group"] = False
        self.drag_state["drag_start_pos"] = (0, 0)
