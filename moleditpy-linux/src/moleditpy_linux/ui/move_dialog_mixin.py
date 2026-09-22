#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Any, Callable, Optional, Set
import logging

import numpy as np
import pyvista as pv
from PyQt6.QtWidgets import QLabel, QLineEdit, QMessageBox

from ..utils.constants import VDW_DISPLAY_RADII


class MoveDialogMixin:
    """Translation, rotation and highlighting shared by the two move dialogs.

    Written against ``group_atoms``: MoveGroupDialog grows it from the clicked
    atom by BFS, MoveSelectedAtomsDialog aliases it onto ``selected_atoms``.
    """

    # Set by each dialog; they differ only in wording and in the actor name.
    HIGHLIGHT_ACTOR = "move_group_highlight"
    EMPTY_SELECTION_WARNING = "Please select a group first."
    NO_SELECTION_TEXT = "No group selected"
    SELECTION_PREFIX = "Selected group"

    # Provided by BasePickingDialog and the concrete dialogs.
    mol: Any
    main_window: Any
    group_atoms: Set[int]
    selected_atoms: Set[int]
    highlight_actor: Any
    selection_label: QLabel
    x_trans_input: QLineEdit
    y_trans_input: QLineEdit
    z_trans_input: QLineEdit
    x_rot_input: QLineEdit
    y_rot_input: QLineEdit
    z_rot_input: QLineEdit
    is_dragging_group: bool
    drag_start_pos: Optional[Any]
    # Annotations only -- assigning these would shadow BasePickingDialog's.
    _update_molecule_geometry: Callable[[Any], None]
    _push_undo: Callable[[], None]

    def _warn(self, message: str) -> None:
        """Show a modal warning parented to this dialog."""
        QMessageBox.warning(self, "Warning", message)  # type: ignore[arg-type]

    def update_display(self) -> None:
        """Refresh the selection label with the current atom count and symbols."""
        if not self.group_atoms:
            self.selection_label.setText(self.NO_SELECTION_TEXT)
            return

        atom_info = [
            f"{self.mol.GetAtomWithIdx(idx).GetSymbol()}({idx})"
            for idx in sorted(self.group_atoms)
        ]
        shown = ", ".join(atom_info[:5])
        if len(atom_info) > 5:
            shown += " ..."
        self.selection_label.setText(
            f"{self.SELECTION_PREFIX}: {len(self.group_atoms)} atoms - {shown}"
        )

    def show_atom_labels(self) -> None:
        """Highlight the selected atoms with translucent spheres."""
        plotter = self.main_window.view_3d_manager.plotter
        try:
            cam = plotter.camera_position if plotter else None
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        self.clear_atom_labels()

        if not self.group_atoms:
            return

        selected_indices = list(self.group_atoms)
        if self.main_window.view_3d_manager.atom_positions_3d is None:
            logging.warning("atom_positions_3d is None in update_atom_labels")
            return
        selected_positions = self.main_window.view_3d_manager.atom_positions_3d[
            selected_indices
        ]
        selected_radii = np.array(
            [
                VDW_DISPLAY_RADII.get(self.mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3
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
            name=self.HIGHLIGHT_ACTOR,
            pickable=False,
            reset_camera=False,
        )

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                logging.debug("Could not restore the camera", exc_info=True)

        plotter.render()

    def clear_atom_labels(self) -> None:
        """Clear the highlight spheres and the base class's own labels."""
        super().clear_atom_labels()  # type: ignore[misc]

        plotter = self.main_window.view_3d_manager.plotter
        if plotter is not None:
            try:
                plotter.remove_actor(self.HIGHLIGHT_ACTOR)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                logging.debug("Could not remove the highlight by name", exc_info=True)

        if self.highlight_actor:
            if plotter is not None:
                try:
                    plotter.remove_actor(self.highlight_actor)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    logging.debug("Could not remove the highlight actor", exc_info=True)
            self.highlight_actor = None

        if plotter is not None:
            try:
                plotter.render()
            except (AttributeError, RuntimeError, ValueError, TypeError):
                logging.debug("Could not render after clearing", exc_info=True)

    def reset_translation_inputs(self) -> None:
        """Reset all translation input fields to zero."""
        self.x_trans_input.setText("0.0")
        self.y_trans_input.setText("0.0")
        self.z_trans_input.setText("0.0")

    def reset_rotation_inputs(self) -> None:
        """Reset all rotation input fields to zero."""
        self.x_rot_input.setText("0.0")
        self.y_rot_input.setText("0.0")
        self.z_rot_input.setText("0.0")

    def apply_translation(self) -> None:
        """Translate the selected atoms."""
        if not self.group_atoms:
            self._warn(self.EMPTY_SELECTION_WARNING)
            return

        try:
            dx = float(self.x_trans_input.text())
            dy = float(self.y_trans_input.text())
            dz = float(self.z_trans_input.text())
        except ValueError:
            self._warn("Please enter valid translation values.")
            return

        translation_vector = np.array([dx, dy, dz])
        positions = self.mol.GetConformer().GetPositions()
        for atom_idx in self.group_atoms:
            positions[atom_idx] += translation_vector

        self._update_molecule_geometry(positions)
        self._push_undo()  # after the modification, not before
        self.show_atom_labels()

    def apply_rotation(self) -> None:
        """Rotate the selected atoms about their centroid."""
        if not self.group_atoms:
            self._warn(self.EMPTY_SELECTION_WARNING)
            return

        try:
            rx_rad, ry_rad, rz_rad = np.radians(
                [
                    float(self.x_rot_input.text()),
                    float(self.y_rot_input.text()),
                    float(self.z_rot_input.text()),
                ]
            )
        except ValueError:
            self._warn("Please enter valid rotation values.")
            return

        positions = self.mol.GetConformer().GetPositions()
        selected_indices = list(self.group_atoms)
        centroid = np.mean(positions[selected_indices], axis=0)

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

        for atom_idx in self.group_atoms:
            pos = positions[atom_idx]
            positions[atom_idx] = rot_matrix @ (pos - centroid) + centroid

        self._update_molecule_geometry(positions)
        self._push_undo()  # after the modification, not before
        self.show_atom_labels()

    def clear_selection(self) -> None:
        """Drop the selection and the highlight that draws it."""
        self.selected_atoms.clear()
        self.group_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        self.is_dragging_group = False
        self.drag_start_pos = None
