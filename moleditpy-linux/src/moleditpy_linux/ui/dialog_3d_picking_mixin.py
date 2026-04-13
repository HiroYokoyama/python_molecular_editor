#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging  # [REPORT ERROR MISSING ATTRIBUTE]

import numpy as np
from PyQt6.QtCore import QEvent, Qt

try:
    from ..utils.constants import pt
except ImportError:
    from moleditpy_linux.utils.constants import pt


class Dialog3DPickingMixin:
    """Mixin providing common functionality for 3D atom selection."""

    def __init__(self):
        """Initialize the Mixin."""
        self.picking_enabled = False
        self._mouse_press_pos = None
        self._mouse_moved = False
        self.selection_labels = []

    def eventFilter(self, obj, event):
        """Capture mouse clicks in the 3D view (reproducibly mimicking the original 3D edit logic)."""
        if (
            obj == self.main_window.view_3d_manager.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonPress
            and event.button() == Qt.MouseButton.LeftButton
        ):
            # Start tracking for smart selection (click vs drag)
            self._mouse_press_pos = event.pos()
            self._mouse_moved = False

            try:
                # Retrieve VTK event coordinates (matches original logic)
                interactor = self.main_window.view_3d_manager.plotter.interactor
                click_pos = interactor.GetEventPosition()
                picker = self.main_window.view_3d_manager.plotter.picker
                picker.Pick(
                    click_pos[0],
                    click_pos[1],
                    0,
                    self.main_window.view_3d_manager.plotter.renderer,
                )

                if picker.GetActor() is self.main_window.view_3d_manager.atom_actor:
                    picked_position = np.array(picker.GetPickPosition())
                    distances = np.linalg.norm(
                        self.main_window.view_3d_manager.atom_positions_3d
                        - picked_position,
                        axis=1,
                    )
                    closest_atom_idx = np.argmin(distances)

                    # Add range check
                    if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                        # Click threshold check (matches original logic)
                        atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                        if atom:
                            try:
                                atomic_num = atom.GetAtomicNum()
                                vdw_radius = pt.GetRvdw(atomic_num)
                                if vdw_radius < 0.1:
                                    vdw_radius = 1.5
                            except (AttributeError, RuntimeError, TypeError):
                                vdw_radius = 1.5
                            click_threshold = vdw_radius * 1.5

                            if distances[closest_atom_idx] < click_threshold:
                                if hasattr(self.main_window, "_picking_consumed"):
                                    self.main_window._picking_consumed = True
                                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                                    logging.error(
                                        "REPORT ERROR: Missing attribute '_picking_consumed' on self.main_window"
                                    )
                                self.on_atom_picked(int(closest_atom_idx))

                                # We picked an atom, so stop tracking for background click
                                self._mouse_press_pos = None
                                return True

                # Clicked something other than an atom
                # Permitting rotation (drag) without immediate selection clearing.
                # Actual clearing is handled in the MouseButtonRelease event.
                return False

            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Error in eventFilter: {e}")
                # On exception, don't swallow the event either — let the normal
                # event pipeline continue so the UI remains responsive.
                return False

        # Add movement tracking for smart selection
        elif (
            obj == self.main_window.view_3d_manager.plotter.interactor
            and event.type() == QEvent.Type.MouseMove
        ):
            if self._mouse_press_pos is not None:
                # Check if moved significantly
                diff = event.pos() - self._mouse_press_pos
                if diff.manhattanLength() > 3:
                    self._mouse_moved = True

        # Add release handling for smart selection
        elif (
            obj == self.main_window.view_3d_manager.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonRelease
            and event.button() == Qt.MouseButton.LeftButton
        ):
            if self._mouse_press_pos is not None:
                if not self._mouse_moved:
                    # Pure click (no drag) on background -> Clear selection
                    if hasattr(self, "clear_selection"):
                        self.clear_selection()
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error(
                            "REPORT ERROR: Missing attribute 'clear_selection' on self"
                        )

                # Reset state
                self._mouse_press_pos = None
                self._mouse_moved = False

        return super().eventFilter(obj, event)

    def enable_picking(self):
        """Enable atom selection in the 3D view."""
        self.main_window.view_3d_manager.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
        # Ensure the main window flag exists
        if hasattr(self.main_window, "_picking_consumed"):
            self.main_window._picking_consumed = False
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute '_picking_consumed' on self.main_window"
            )

    def disable_picking(self):
        """Disable atom selection in the 3D view."""
        if self.picking_enabled:
            self.main_window.view_3d_manager.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        if hasattr(self.main_window, "_picking_consumed"):
            self.main_window._picking_consumed = False
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute '_picking_consumed' on self.main_window"
            )

    def try_alternative_picking(self, x, y):
        """Alternative picking method (unused)."""

    # ------------------------------------------------------------------
    # Label management (shared across dialogs)
    # ------------------------------------------------------------------

    def clear_atom_labels(self):
        """Remove all label actors from the plotter."""
        for label_actor in self.selection_labels:
            try:
                if label_actor is not None:
                    self.main_window.view_3d_manager.plotter.remove_actor(label_actor)
            except (AttributeError, RuntimeError, TypeError):
                # Ignore actor removal failure on stale plotter
                pass

        self.selection_labels = []

    # Alias — some dialogs use this name instead.
    clear_selection_labels = clear_atom_labels

    def add_selection_label(self, atom_idx, label_text, color="yellow"):
        """Add a point label at the position of *atom_idx*.

        Parameters
        ----------
        atom_idx : int
            Index into ``self.main_window.view_3d_manager.atom_positions_3d``.
        label_text : str
            Text shown next to the atom.
        color : str, optional
            Label colour (default ``'yellow'``).
        """
        plotter = self.main_window.view_3d_manager.plotter
        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        pos = self.main_window.view_3d_manager.atom_positions_3d[atom_idx]
        label_actor = plotter.add_point_labels(
            [pos],
            [label_text],
            point_size=20,
            font_size=12,
            text_color=color,
            always_visible=True,
        )
        self.selection_labels.append(label_actor)

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                pass

    def show_atom_labels_for(self, atoms_and_labels, color="yellow"):
        """Clear existing labels and add new ones for each *(idx, text)* pair.

        Parameters
        ----------
        atoms_and_labels : list[tuple[int, str]]
            Each element is ``(atom_idx, label_text)``.
        color : str, optional
            Label colour (default ``'yellow'``).
        """
        plotter = self.main_window.view_3d_manager.plotter
        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        self.clear_atom_labels()

        for atom_idx, label_text in atoms_and_labels:
            pos = self.main_window.view_3d_manager.atom_positions_3d[atom_idx]
            label_actor = plotter.add_point_labels(
                [pos],
                [label_text],
                point_size=20,
                font_size=12,
                text_color=color,
                always_visible=True,
            )
            self.selection_labels.append(label_actor)

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                pass
