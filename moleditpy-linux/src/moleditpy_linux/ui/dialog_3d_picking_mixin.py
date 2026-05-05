#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

from PyQt6.QtCore import QEvent, Qt, QObject, QPoint
from PyQt6.QtGui import QMouseEvent
from typing import Any, Optional, TYPE_CHECKING

try:
    from .atom_picking import pick_atom_index_from_screen
except ImportError:
    from moleditpy_linux.ui.atom_picking import pick_atom_index_from_screen

if TYPE_CHECKING:
    from .main_window import MainWindow
    from rdkit import Chem


class Dialog3DPickingMixin:
    """Mixin providing common functionality for 3D atom selection."""

    main_window: MainWindow
    mol: Chem.Mol

    def __init__(self) -> None:
        """Initialize the Mixin."""
        self.picking_enabled = False
        self._mouse_press_pos: Optional[QPoint] = None
        self._mouse_moved = False
        self._consume_next_left_release = False
        self.selection_labels: list[Any] = []

    def eventFilter(self, obj: Optional[QObject], event: Optional[QEvent]) -> bool:
        """Capture mouse clicks in the 3D view (reproducibly mimicking the original 3D edit logic)."""
        if event is None:
            return False

        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None or self.mol is None:
            return False

        if (
            obj == plotter.interactor
            and event.type() == QEvent.Type.MouseButtonPress
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            # Start tracking for smart selection (click vs drag)
            self._mouse_press_pos = event.pos()
            self._mouse_moved = False

            try:
                # Retrieve VTK event coordinates (matches original logic)
                interactor = plotter.interactor
                if interactor is None:
                    return False
                click_pos = interactor.GetEventPosition()
                closest_atom_idx = pick_atom_index_from_screen(
                    self.main_window.view_3d_manager,
                    (int(click_pos[0]), int(click_pos[1])),
                    self.mol,
                )

                if (
                    closest_atom_idx is not None
                    and 0 <= closest_atom_idx < self.mol.GetNumAtoms()
                ):
                    atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        if hasattr(self.main_window, "_picking_consumed"):
                            self.main_window._picking_consumed = True

                        def _deferred_pick(idx=int(closest_atom_idx), target=self):
                            try:
                                target.on_atom_picked(idx)
                            except (AttributeError, RuntimeError):
                                pass

                        from PyQt6.QtCore import QTimer

                        QTimer.singleShot(0, _deferred_pick)

                        # We picked an atom, so stop tracking for background click
                        self._mouse_press_pos = None
                        self._consume_next_left_release = True
                        return True

                # Clicked something other than an atom
                return False

            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Error in eventFilter: {e}")
                return False

        # Add movement tracking for smart selection
        elif (
            obj == plotter.interactor
            and event.type() == QEvent.Type.MouseMove
            and isinstance(event, QMouseEvent)
        ):
            if self._mouse_press_pos is not None:
                # Check if moved significantly
                diff = event.pos() - self._mouse_press_pos
                if diff.manhattanLength() > 3:
                    self._mouse_moved = True

        # Add release handling for smart selection
        elif (
            obj == plotter.interactor
            and event.type() == QEvent.Type.MouseButtonRelease
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            if self._consume_next_left_release:
                self._consume_next_left_release = False
                self._mouse_press_pos = None
                self._mouse_moved = False
                return True

            if self._mouse_press_pos is not None:
                if not self._mouse_moved:
                    # Pure click (no drag) on background -> Clear selection
                    if hasattr(self, "clear_selection"):

                        def _deferred_clear(target=self):
                            try:
                                target.clear_selection()
                            except (AttributeError, RuntimeError):
                                pass

                        from PyQt6.QtCore import QTimer

                        QTimer.singleShot(0, _deferred_clear)

                # Reset state
                self._mouse_press_pos = None
                self._mouse_moved = False

        return False

    def on_atom_picked(self, atom_idx: int) -> None:
        """Handle atom picking. Must be implemented by the dialog."""
        raise NotImplementedError("on_atom_picked must be implemented by subclasses")

    def enable_picking(self) -> None:
        """Enable atom selection in the 3D view."""
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is not None and plotter.interactor is not None:
            plotter.interactor.installEventFilter(self)
            self.picking_enabled = True
        if hasattr(self.main_window, "_picking_consumed"):
            self.main_window._picking_consumed = False

    def disable_picking(self) -> None:
        """Disable atom selection in the 3D view."""
        if self.picking_enabled:
            plotter = self.main_window.view_3d_manager.plotter
            if plotter is not None and plotter.interactor is not None:
                plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        self._consume_next_left_release = False
        if hasattr(self.main_window, "_picking_consumed"):
            self.main_window._picking_consumed = False

    def clear_atom_labels(self) -> None:
        """Remove all label actors from the plotter."""
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is not None:
            for label_actor in self.selection_labels:
                try:
                    if label_actor is not None:
                        plotter.remove_actor(label_actor)
                except (AttributeError, RuntimeError, TypeError):
                    pass
        self.selection_labels = []

    clear_selection_labels = clear_atom_labels

    def add_selection_label(
        self, atom_idx: int, label_text: str, color: str = "yellow"
    ) -> None:
        """Add a point label at the position of *atom_idx*."""
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None:
            return

        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        if self.main_window.view_3d_manager.atom_positions_3d is not None:
            pos = self.main_window.view_3d_manager.atom_positions_3d[atom_idx]
            label_actor = plotter.add_point_labels(
                [pos],
                [label_text],
                point_size=0,
                font_size=12,
                text_color=color,
                always_visible=True,
                show_points=False,
                shape="rect",
                shape_color="gray",
                shape_opacity=0.5,
            )
            self.selection_labels.append(label_actor)

            if cam is not None:
                try:
                    plotter.camera_position = cam
                except (AttributeError, RuntimeError, TypeError):
                    pass

    def show_atom_labels_for(
        self, atoms_and_labels: list[tuple[int, str]], color: str = "yellow"
    ) -> None:
        """Clear existing labels and add new ones for each *(idx, text)* pair."""
        plotter = self.main_window.view_3d_manager.plotter
        if plotter is None:
            return

        try:
            cam = plotter.camera_position
        except (AttributeError, RuntimeError, TypeError):
            cam = None

        self.clear_atom_labels()

        if self.main_window.view_3d_manager.atom_positions_3d is not None:
            for atom_idx, label_text in atoms_and_labels:
                pos = self.main_window.view_3d_manager.atom_positions_3d[atom_idx]
                label_actor = plotter.add_point_labels(
                    [pos],
                    [label_text],
                    point_size=0,
                    font_size=12,
                    text_color=color,
                    always_visible=True,
                    show_points=False,
                    shape="rect",
                    shape_color="gray",
                    shape_opacity=0.5,
                )
                self.selection_labels.append(label_actor)

        if cam is not None:
            try:
                plotter.camera_position = cam
            except (AttributeError, RuntimeError, TypeError):
                pass
