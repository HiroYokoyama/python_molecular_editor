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

import logging
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtWidgets import QApplication
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera  # pylint: disable=no-name-in-module


from .atom_picking import pick_atom_index_from_screen
from ..utils.constants import MOVE_DIALOG_TYPES

from rdkit import Geometry

# VTK's trackball state id for camera rotation (VTKIS_ROTATE); the enum is not
# exposed to Python, so the literal is used with this name for clarity.
_VTKIS_ROTATE = 1

# Reference 3D-canvas size (px) that camera rotation speed is normalized to.
# VTK's stock Rotate divides motion by the live render size, so rotation slows
# as the canvas grows; normalizing to a fixed reference makes the speed
# window-size independent and matches the feel at the default window size.
_ROTATION_REFERENCE_SIZE = 640.0

# Minimum seconds between real-time drag redraws (~30 fps).
_DRAG_REDRAW_INTERVAL = 0.033

# Largest molecule that still gets real-time drag feedback (each frame is a
# full scene rebuild).
_REALTIME_DRAG_MAX_ATOMS = 300

# Radians of group rotation per pixel of mouse travel, before sensitivity.
_GROUP_ROTATION_RADIANS_PER_PIXEL = 0.008


def _rodrigues_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Return the rotation matrix for `angle` radians about `axis` (Rodrigues)."""
    norm = float(np.linalg.norm(axis))
    if norm < 1e-6 or abs(angle) < 1e-6:
        return np.eye(3)
    unit = axis / norm
    k_mat = np.array(
        [
            [0.0, -unit[2], unit[1]],
            [unit[2], 0.0, -unit[0]],
            [-unit[1], unit[0], 0.0],
        ]
    )
    return np.eye(3) + np.sin(angle) * k_mat + (1.0 - np.cos(angle)) * (k_mat @ k_mat)


class CustomInteractorStyle(vtkInteractorStyleTrackballCamera):
    """VTK interactor style extending trackball-camera with 3D atom drag and measurement."""

    def __init__(self, main_window: Any = None, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.main_window = main_window
        # Custom state flags
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos: Optional[tuple[int, int]] = None
        self._suppress_next_left_button_up = False
        # Throttle for real-time drag redraws
        self._last_drag_redraw_time: float = 0.0
        self._drag_redraw_pending = False
        self._active_drag_atoms: Optional[List[int]] = None

        self.AddObserver("LeftButtonPressEvent", self.on_left_button_down)  # type: ignore[arg-type]
        # self.AddObserver("LeftButtonDoubleClickEvent", self.on_left_button_down)
        self.AddObserver("RightButtonPressEvent", self.on_right_button_down)  # type: ignore[arg-type]
        self.AddObserver("MouseMoveEvent", self.on_mouse_move)  # type: ignore[arg-type]
        self.AddObserver("LeftButtonReleaseEvent", self.on_left_button_up)  # type: ignore[arg-type]
        self.AddObserver("RightButtonReleaseEvent", self.on_right_button_up)  # type: ignore[arg-type]

    def reset_interactor_state(self) -> None:
        """Force-reset all VTK and custom drag/selection state.

        Called when 3D-edit or measurement mode is toggled off, and
        defensively at the start of every new mouse-press, so that a
        previously corrupted VTK state machine cannot cause a permanent lock.
        """
        # Reset VTK's internal state machine (clears ROTATE / PAN / etc.)
        try:
            self.StopState()
        except (AttributeError, RuntimeError):
            # Safe defensive fallback catching AttributeError, RuntimeError
            logging.debug("Suppressed non-critical error", exc_info=True)

        self._end_drag_event()

        # Reset all custom flags
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None
        self._suppress_next_left_button_up = False

        # Clear dangling atom drag info on main window
        mw = self.main_window
        if mw is not None:
            try:
                mw.dragged_atom_info = None
            except (AttributeError, RuntimeError):
                # Safe defensive fallback catching AttributeError, RuntimeError
                logging.debug("Suppressed non-critical error", exc_info=True)
            try:
                mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except (AttributeError, RuntimeError):
                # Safe defensive fallback catching AttributeError, RuntimeError
                logging.debug("Suppressed non-critical error", exc_info=True)

    def _find_move_dialog(self) -> Any:
        """Return the visible Move Group / Move Selected Atoms dialog, if any."""
        try:
            for widget in QApplication.topLevelWidgets():
                if type(widget).__name__ in MOVE_DIALOG_TYPES and widget.isVisible():
                    return widget
        except (AttributeError, RuntimeError, TypeError):
            logging.debug("Suppressed non-critical error", exc_info=True)
        return None

    def _display_to_world(
        self,
        renderer: Any,
        screen_x: float,
        screen_y: float,
        depth_z: float,
    ) -> Optional[Tuple[float, float, float]]:
        """Convert screen (display) coordinates to 3D world coordinates.

        Uses the given depth_z (already in display-space) so the result
        lies on the same depth plane as the reference atom.
        """
        try:
            renderer.SetDisplayPoint(screen_x, screen_y, depth_z)
            renderer.DisplayToWorld()
            wp = renderer.GetWorldPoint()
            return (wp[0], wp[1], wp[2])
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
            return None

    def _world_to_display_depth(
        self, renderer: Any, x: float, y: float, z: float
    ) -> Optional[float]:
        """Return the display-space Z depth for a world-space point."""
        try:
            renderer.SetWorldPoint(x, y, z, 1.0)
            renderer.WorldToDisplay()
            return float(renderer.GetDisplayPoint()[2])
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
            return None

    def _invoke_drag_handlers(
        self,
        event_type: str,
        atom_indices: List[int],
        positions: Dict[int, Tuple[float, float, float]],
    ) -> None:
        """Forward a drag event to all registered plugin handlers."""
        mw = self.main_window
        if mw is None:
            return
        pm = getattr(mw, "plugin_manager", None)
        if pm is None:
            return
        invoke = getattr(pm, "invoke_atom_drag_handlers", None)
        if invoke is not None:
            try:
                invoke(event_type, atom_indices, positions)
            # Plugin faults are already isolated per handler by the manager.
            except (AttributeError, RuntimeError, TypeError, ValueError):
                logging.warning("invoke_atom_drag_handlers raised", exc_info=True)

    def _current_positions(
        self, atom_indices: List[int]
    ) -> Dict[int, Tuple[float, float, float]]:
        """Read the current conformer coordinates of the given atoms."""
        positions: Dict[int, Tuple[float, float, float]] = {}
        try:
            mol = self.main_window.view_3d_manager.current_mol
            if mol is None or mol.GetNumConformers() == 0:
                return positions
            conf = mol.GetConformer()
            for idx in atom_indices:
                pos = conf.GetAtomPosition(idx)
                positions[idx] = (float(pos.x), float(pos.y), float(pos.z))
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
        return positions

    def _begin_drag_event(self, atom_indices: Any) -> None:
        """Report the start of a drag gesture and arm the redraw throttle.

        A right-button press never passes through reset_interactor_state, so
        close any gesture whose end was missed before opening a new one.
        """
        self._end_drag_event()
        atoms = [int(idx) for idx in atom_indices]
        self._active_drag_atoms = atoms
        self._last_drag_redraw_time = 0.0
        self._drag_redraw_pending = False
        self._invoke_drag_handlers("start", list(atoms), {})

    def _end_drag_event(self) -> None:
        """Report the end of the active drag gesture exactly once.

        Every exit path (normal release, aborted gesture, forced reset) routes
        through here so a plugin can never be left believing a drag is still
        in progress.
        """
        atoms = self._active_drag_atoms
        self._active_drag_atoms = None
        self._drag_redraw_pending = False
        if atoms is None:
            return
        self._invoke_drag_handlers("end", list(atoms), self._current_positions(atoms))

    def _realtime_drag_active(self) -> bool:
        """Return True when real-time drag is both enabled and affordable here.

        Every frame rebuilds the whole 3D scene, so past a few hundred atoms
        the rebuild costs more than the throttle interval and dragging feels
        worse with it than without; the move is then applied on release only.
        """
        if not self._realtime_drag_enabled():
            return False
        return not self._molecule_too_large_for_realtime()

    def _molecule_too_large_for_realtime(self) -> bool:
        """Return True when the molecule exceeds the real-time drag size limit."""
        try:
            mol = self.main_window.view_3d_manager.current_mol
            if mol is None:
                return False
            return int(mol.GetNumAtoms()) > _REALTIME_DRAG_MAX_ATOMS
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
            return False

    def _should_drag_redraw(self) -> bool:
        """Return True when a real-time drag frame may be rendered now.

        Skips while a previous deferred redraw is still queued so a scene that
        rebuilds slower than the throttle interval cannot pile up frames.
        """
        if not self._realtime_drag_active():
            return False
        if self._drag_redraw_pending:
            return False
        return time.monotonic() - self._last_drag_redraw_time >= _DRAG_REDRAW_INTERVAL

    def _finish_drag_redraw(self) -> None:
        """Release the redraw gate once a deferred drag frame has rendered."""
        self._drag_redraw_pending = False
        self._last_drag_redraw_time = time.monotonic()

    def _abort_active_drag(self, moved: bool) -> None:
        """Close out a gesture whose release event was lost.

        Real-time dragging has already written the new geometry, so a moved
        gesture still needs an undo entry even though it never reached the
        release handler.
        """
        if self._active_drag_atoms is None:
            return
        if moved and self._realtime_drag_active():
            self._push_undo_after_aborted_drag()
        self._end_drag_event()

    def _push_undo_after_aborted_drag(self) -> None:
        """Record geometry left behind by a gesture whose release was lost.

        Real-time dragging mutates the conformer on every frame, so an aborted
        gesture would otherwise leave the molecule moved with nothing on the
        undo stack.
        """
        try:
            self.main_window.edit_actions_manager.push_undo_state()
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)

    def _realtime_drag_enabled(self) -> bool:
        """Return True if real-time drag redraw is enabled in settings."""
        mw = self.main_window
        if mw is None:
            return True
        try:
            return bool(mw.get_settings().get("realtime_3d_drag", True))
        except (AttributeError, RuntimeError, TypeError):
            return True

    def _rotation_sensitivity(self) -> float:
        """Return the user's rotation-speed multiplier, shared by camera and group rotation."""
        mw = self.main_window
        if mw is None:
            return 1.0
        try:
            return float(mw.get_settings().get("mouse_rotation_sensitivity", 1.0))
        except (AttributeError, TypeError, ValueError):
            return 1.0

    def _rotate_size_independent(self) -> None:
        """Rotate the camera at a speed that does not depend on canvas size.

        Mirrors vtkInteractorStyleTrackballCamera::Rotate but normalizes the
        per-pixel azimuth/elevation to a fixed reference size instead of the
        live render size, and scales by the user's rotation-sensitivity setting.
        """
        renderer = self.GetCurrentRenderer()
        interactor = self.GetInteractor()
        if renderer is None or interactor is None:
            return

        dx = interactor.GetEventPosition()[0] - interactor.GetLastEventPosition()[0]
        dy = interactor.GetEventPosition()[1] - interactor.GetLastEventPosition()[1]

        delta = (
            -20.0
            / _ROTATION_REFERENCE_SIZE
            * self.GetMotionFactor()
            * self._rotation_sensitivity()
        )
        camera = renderer.GetActiveCamera()
        camera.Azimuth(dx * delta)
        camera.Elevation(dy * delta)
        camera.OrthogonalizeViewUp()

        if self.GetAutoAdjustCameraClippingRange():
            renderer.ResetCameraClippingRange()
        if interactor.GetLightFollowCamera():
            renderer.UpdateLightsGeometryToFollowCamera()
        interactor.Render()

    def _stop_vtk_left_button_state(self) -> None:
        """Clear VTK's button/drag state after a custom-handled left click."""
        try:
            self.StopState()
        except (AttributeError, RuntimeError):
            # Safe defensive fallback catching AttributeError, RuntimeError
            logging.debug("Suppressed non-critical error", exc_info=True)

    def on_left_button_down(self, obj: Any, event: Any) -> None:
        """
        Dispatch click events.
        Use custom action if atom handles, else camera rotation.
        """
        # Defensively reset VTK state at the start of every new press so that
        # any previously stuck state (e.g. from a missed OnLeftButtonDown call)
        # is cleared before we process the new event.
        self.reset_interactor_state()

        mw = self.main_window

        # Clear previous drag state
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None

        move_group_dialog = self._find_move_dialog()

        if move_group_dialog and move_group_dialog.group_atoms:
            # Group drag if selected
            click_pos = self.GetInteractor().GetEventPosition()
            clicked_atom_idx = pick_atom_index_from_screen(
                mw.view_3d_manager,
                (int(click_pos[0]), int(click_pos[1])),
                mw.view_3d_manager.current_mol,
            )

            # If an atom in the group is clicked
            if clicked_atom_idx is not None:
                if clicked_atom_idx in move_group_dialog.group_atoms:
                    # Preparation for group drag
                    move_group_dialog.is_dragging_group_vtk = True
                    move_group_dialog.drag_atom_idx_vtk = clicked_atom_idx
                    move_group_dialog.drag_start_pos_vtk = click_pos
                    move_group_dialog.mouse_moved_vtk = False
                    # Save initial positions
                    move_group_dialog.initial_positions = {}
                    conf = mw.view_3d_manager.current_mol.GetConformer()
                    for atom_idx in move_group_dialog.group_atoms:
                        pos = conf.GetAtomPosition(atom_idx)
                        move_group_dialog.initial_positions[atom_idx] = np.array(
                            [pos.x, pos.y, pos.z]
                        )
                    mw.view_3d_manager.plotter.setCursor(
                        Qt.CursorShape.ClosedHandCursor
                    )
                    self._begin_drag_event(move_group_dialog.group_atoms)
                    self._suppress_next_left_button_up = True
                    return  # Disable camera rotation
                else:
                    if type(move_group_dialog).__name__ == "MoveSelectedAtomsDialog":
                        # For MoveSelectedAtomsDialog, we toggle ONLY the clicked atom, no BFS!
                        def _deferred_toggle() -> None:
                            try:
                                move_group_dialog.on_atom_picked(clicked_atom_idx)
                            except (AttributeError, RuntimeError):
                                logging.warning(
                                    "Move-dialog atom toggle failed", exc_info=True
                                )

                        QTimer.singleShot(0, _deferred_toggle)
                        self._suppress_next_left_button_up = True
                        return

                    # Clicked outside group - Search connected component
                    visited = set()
                    queue = [clicked_atom_idx]
                    visited.add(clicked_atom_idx)

                    while queue:
                        current_idx = queue.pop(0)
                        for bond_idx in range(
                            mw.view_3d_manager.current_mol.GetNumBonds()
                        ):
                            bond = mw.view_3d_manager.current_mol.GetBondWithIdx(
                                bond_idx
                            )
                            begin_idx = bond.GetBeginAtomIdx()
                            end_idx = bond.GetEndAtomIdx()

                            if begin_idx == current_idx and end_idx not in visited:
                                visited.add(end_idx)
                                queue.append(end_idx)
                            elif end_idx == current_idx and begin_idx not in visited:
                                visited.add(begin_idx)
                                queue.append(begin_idx)

                    # Multi-selection with Ctrl
                    is_ctrl_pressed = bool(
                        (
                            QApplication.keyboardModifiers()
                            & Qt.KeyboardModifier.ControlModifier
                        )
                        or (
                            self.GetInteractor()
                            and self.GetInteractor().GetControlKey()
                        )
                    )

                    if is_ctrl_pressed:
                        # Ctrl + Click: toggle selection
                        if visited.issubset(move_group_dialog.group_atoms):
                            # Already selected - deselect
                            move_group_dialog.group_atoms -= visited
                        else:
                            # Add new group
                            move_group_dialog.group_atoms |= visited
                    else:
                        # Replace existing selection
                        move_group_dialog.group_atoms = visited.copy()

                    move_group_dialog.selected_atoms.add(clicked_atom_idx)

                    def _deferred_move_group_update() -> None:
                        try:
                            move_group_dialog.show_atom_labels()
                            move_group_dialog.update_display()
                        except (AttributeError, RuntimeError):
                            # Safe defensive fallback catching AttributeError, RuntimeError
                            logging.debug(
                                "Suppressed non-critical error", exc_info=True
                            )

                    QTimer.singleShot(0, _deferred_move_group_update)
                    return
            else:
                # Track mouse event to distinguish rotation from click
                self._mouse_press_pos = self.GetInteractor().GetEventPosition()
                self._mouse_moved_during_drag = False

                # Allow camera rotation
                super(CustomInteractorStyle, self).OnLeftButtonDown()
                return

        interactor = self.GetInteractor()
        is_temp_mode = bool(
            (QApplication.keyboardModifiers() & Qt.KeyboardModifier.AltModifier)
            or (interactor and interactor.GetAltKey())
        )
        is_edit_active = mw.edit_3d_manager.is_3d_edit_mode or is_temp_mode

        # Handle measurement mode
        if mw.edit_3d_manager.measurement_mode and mw.view_3d_manager.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            self._mouse_moved_during_drag = False

            closest_atom_idx = pick_atom_index_from_screen(
                mw.view_3d_manager,
                (int(click_pos[0]), int(click_pos[1])),
                mw.view_3d_manager.current_mol,
            )

            # Special handling if atom clicked
            if closest_atom_idx is not None:
                # Add range check
                if 0 <= closest_atom_idx < mw.view_3d_manager.current_mol.GetNumAtoms():
                    # Check click threshold
                    atom = mw.view_3d_manager.current_mol.GetAtomWithIdx(
                        int(closest_atom_idx)
                    )
                    if atom:
                        if True:

                            def _deferred_measure() -> None:
                                try:
                                    mw.edit_3d_manager.handle_measurement_atom_selection(
                                        closest_atom_idx
                                    )
                                except (AttributeError, RuntimeError):
                                    logging.warning(
                                        "Measurement selection failed", exc_info=True
                                    )

                            QTimer.singleShot(0, _deferred_measure)
                            self._suppress_next_left_button_up = True
                            return  # Selection complete, disable camera rotation

            # Clear measurement if not dragging
            self._is_dragging_atom = False
            self._mouse_press_pos = click_pos
            super().OnLeftButtonDown()
            return

        # Handle selection if 3D mol exists
        if is_edit_active and mw.view_3d_manager.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            closest_atom_idx = pick_atom_index_from_screen(
                mw.view_3d_manager,
                (int(click_pos[0]), int(click_pos[1])),
                mw.view_3d_manager.current_mol,
            )

            if closest_atom_idx is not None:
                # Add range check
                if 0 <= closest_atom_idx < mw.view_3d_manager.current_mol.GetNumAtoms():
                    # Get atom safely from RDKit Mol
                    atom = mw.view_3d_manager.current_mol.GetAtomWithIdx(
                        int(closest_atom_idx)
                    )
                    if atom:
                        if True:
                            # Successfully grabbed atom
                            self._is_dragging_atom = True
                            self.is_dragging = False
                            mw.dragged_atom_info = {"id": int(closest_atom_idx)}
                            mw.view_3d_manager.plotter.setCursor(
                                Qt.CursorShape.ClosedHandCursor
                            )
                            self._begin_drag_event([int(closest_atom_idx)])
                            self._suppress_next_left_button_up = True
                            return  # Prevent camera rotation

        # Track mouse event to distinguish rotation from click
        self._mouse_press_pos = self.GetInteractor().GetEventPosition()
        self._mouse_moved_during_drag = False
        self._is_dragging_atom = False
        super().OnLeftButtonDown()

    def on_right_button_down(self, obj: Any, event: Any) -> None:
        """
        Right-click: Start group rotation if dialog open.
        """
        mw = self.main_window

        move_group_dialog = self._find_move_dialog()

        if move_group_dialog and move_group_dialog.group_atoms:
            # Start rotation drag if group selected
            click_pos = self.GetInteractor().GetEventPosition()
            clicked_atom_idx = pick_atom_index_from_screen(
                mw.view_3d_manager,
                (int(click_pos[0]), int(click_pos[1])),
                mw.view_3d_manager.current_mol,
            )

            follow_mouse = False
            try:
                follow_mouse = bool(
                    mw.get_settings().get("rotate_group_follow_mouse", False)
                )
            except (AttributeError, RuntimeError, TypeError):
                follow_mouse = False

            # If follow_mouse is True, require clicking directly on an atom in the group.
            # If follow_mouse is False (default), allow right-clicking anywhere on display.
            if not follow_mouse or (
                clicked_atom_idx is not None
                and clicked_atom_idx in move_group_dialog.group_atoms
            ):
                move_group_dialog.is_rotating_group_vtk = True
                move_group_dialog.rotation_start_pos = click_pos
                move_group_dialog.rotation_mouse_moved = False
                if clicked_atom_idx is not None:
                    move_group_dialog.rotation_atom_idx = (
                        clicked_atom_idx  # Record grabbed atom
                    )
                elif hasattr(move_group_dialog, "rotation_atom_idx"):
                    delattr(move_group_dialog, "rotation_atom_idx")

                # Save initial positions and centroid
                move_group_dialog.initial_positions = {}
                conf = mw.view_3d_manager.current_mol.GetConformer()
                centroid = np.zeros(3)
                for atom_idx in move_group_dialog.group_atoms:
                    pos = conf.GetAtomPosition(atom_idx)
                    pos_array = np.array([pos.x, pos.y, pos.z])
                    move_group_dialog.initial_positions[atom_idx] = pos_array
                    centroid += pos_array
                centroid /= len(move_group_dialog.group_atoms)
                move_group_dialog.group_centroid = centroid

                mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                self._begin_drag_event(move_group_dialog.group_atoms)
                return  # Disable camera rotation

        # Standard right-click
        super().OnRightButtonDown()

    def _heal_stuck_pointer_state(self, move_group_dialog: Any) -> None:
        """Self-heal drag/rotate state whose release event was lost.

        A release can be lost mid-gesture (e.g. pyvista temporarily swapping
        the interactor style, a dialog grabbing events), leaving this style or
        a Move Group dialog in a permanent drag/rotate state that blocks
        interaction. Verify against the real button state on every move and
        reset when the button is not actually held.
        """
        try:
            buttons = QApplication.mouseButtons()
            left_held = bool(buttons & Qt.MouseButton.LeftButton)
            right_held = bool(buttons & Qt.MouseButton.RightButton)
            any_held = buttons != Qt.MouseButton.NoButton
        except (AttributeError, RuntimeError, TypeError):
            return

        if self._is_dragging_atom and not left_held:
            self._abort_active_drag(self.is_dragging)
            self.reset_interactor_state()

        if not any_held:
            try:
                if self.GetState() != 0:  # stuck ROTATE/PAN/etc. without a button
                    self.reset_interactor_state()
            except (AttributeError, RuntimeError, TypeError):
                # Safe defensive fallback catching AttributeError, RuntimeError, TypeError
                logging.debug("Suppressed non-critical error", exc_info=True)

        if move_group_dialog is not None:
            try:
                if (
                    getattr(move_group_dialog, "is_dragging_group_vtk", False)
                    and not left_held
                ):
                    moved = getattr(move_group_dialog, "mouse_moved_vtk", False)
                    move_group_dialog.is_dragging_group_vtk = False
                    move_group_dialog.drag_start_pos_vtk = None
                    move_group_dialog.mouse_moved_vtk = False
                    self._abort_active_drag(moved)
                if (
                    getattr(move_group_dialog, "is_rotating_group_vtk", False)
                    and not right_held
                ):
                    moved = getattr(move_group_dialog, "rotation_mouse_moved", False)
                    move_group_dialog.is_rotating_group_vtk = False
                    move_group_dialog.rotation_start_pos = None
                    move_group_dialog.rotation_mouse_moved = False
                    self._abort_active_drag(moved)
            except (AttributeError, RuntimeError, TypeError):
                # Safe defensive fallback catching AttributeError, RuntimeError, TypeError
                logging.debug("Suppressed non-critical error", exc_info=True)

    def on_mouse_move(self, obj: Any, event: Any) -> None:
        """
        Handle mouse move (drag vs camera/hover).
        """
        mw = self.main_window

        move_group_dialog = self._find_move_dialog()

        self._heal_stuck_pointer_state(move_group_dialog)
        if move_group_dialog and getattr(
            move_group_dialog, "is_dragging_group_vtk", False
        ):
            # Dragging group - record offset and update positions in real-time
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()

            if move_group_dialog.drag_start_pos_vtk is None:
                return

            dx = current_pos[0] - move_group_dialog.drag_start_pos_vtk[0]
            dy = current_pos[1] - move_group_dialog.drag_start_pos_vtk[1]

            if abs(dx) > 5 or abs(dy) > 5:
                move_group_dialog.mouse_moved_vtk = True

            if (
                move_group_dialog.mouse_moved_vtk
                and hasattr(move_group_dialog, "initial_positions")
                and self._should_drag_redraw()
            ):
                self._do_realtime_group_translate(mw, move_group_dialog, current_pos)

            return  # Disable camera rotation

        # Group rotation handling
        if move_group_dialog and getattr(
            move_group_dialog, "is_rotating_group_vtk", False
        ):
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()

            if move_group_dialog.rotation_start_pos is None:
                return

            dx = current_pos[0] - move_group_dialog.rotation_start_pos[0]
            dy = current_pos[1] - move_group_dialog.rotation_start_pos[1]

            if abs(dx) > 5 or abs(dy) > 5:
                move_group_dialog.rotation_mouse_moved = True

            if (
                move_group_dialog.rotation_mouse_moved
                and hasattr(move_group_dialog, "initial_positions")
                and hasattr(move_group_dialog, "group_centroid")
                and self._should_drag_redraw()
            ):
                self._do_realtime_group_rotate(mw, move_group_dialog, current_pos)

            return  # Disable camera rotation

        interactor = self.GetInteractor()

        # Record mouse movement
        if self._mouse_press_pos is not None:
            current_pos = interactor.GetEventPosition()
            if (
                abs(current_pos[0] - self._mouse_press_pos[0]) > 3
                or abs(current_pos[1] - self._mouse_press_pos[1]) > 3
            ):
                self._mouse_moved_during_drag = True

        if self._is_dragging_atom and mw.dragged_atom_info is not None:
            # Custom atom drag
            self.is_dragging = True
            if self._should_drag_redraw():
                self._do_realtime_atom_drag(mw)
        else:
            # Camera interaction. VTK's built-in Rotate divides mouse motion by
            # the live render size, so rotation slows as the canvas grows. Handle
            # the rotate state ourselves for window-size-independent speed;
            # delegate pan/dolly/spin to the parent unchanged.
            if self.GetState() == _VTKIS_ROTATE:
                self._rotate_size_independent()
            else:
                super().OnMouseMove()

            # Update cursor display
            is_edit_active = (
                mw.edit_3d_manager.is_3d_edit_mode or interactor.GetAltKey()
            )
            if is_edit_active:
                # Hover check if edit active
                click_pos = interactor.GetEventPosition()
                atom_under_cursor = (
                    pick_atom_index_from_screen(
                        mw.view_3d_manager,
                        (int(click_pos[0]), int(click_pos[1])),
                        mw.view_3d_manager.current_mol,
                    )
                    is not None
                )

                if atom_under_cursor:
                    mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                else:
                    mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            else:
                mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.ArrowCursor)

    def _do_realtime_atom_drag(self, mw: Any) -> None:
        """Update the dragged atom's position and redraw during mouse-move."""
        try:
            atom_id = mw.dragged_atom_info.get("id") if mw.dragged_atom_info else None
            if atom_id is None:
                return
            if not (
                mw.view_3d_manager.current_mol
                and mw.view_3d_manager.current_mol.GetNumConformers() > 0
            ):
                return

            interactor = self.GetInteractor()
            renderer = mw.view_3d_manager.plotter.renderer
            current_pos = interactor.GetEventPosition()
            conf = mw.view_3d_manager.current_mol.GetConformer()
            pos_3d = conf.GetAtomPosition(atom_id)

            depth_z = self._world_to_display_depth(
                renderer, pos_3d.x, pos_3d.y, pos_3d.z
            )
            if depth_z is None:
                return
            new_world = self._display_to_world(
                renderer, current_pos[0], current_pos[1], depth_z
            )
            if new_world is None:
                return

            conf.SetAtomPosition(
                atom_id,
                Geometry.Point3D(
                    float(new_world[0]), float(new_world[1]), float(new_world[2])
                ),
            )
            if isinstance(
                mw.view_3d_manager.atom_positions_3d, (list, np.ndarray)
            ) and atom_id < len(mw.view_3d_manager.atom_positions_3d):
                mw.view_3d_manager.atom_positions_3d[atom_id] = list(new_world)

            _rt_mol = mw.view_3d_manager.current_mol

            def _deferred_rt_atom() -> None:
                try:
                    mw.view_3d_manager.draw_molecule_3d(_rt_mol)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    logging.debug("Suppressed non-critical error", exc_info=True)
                finally:
                    self._finish_drag_redraw()

            self._drag_redraw_pending = True
            QTimer.singleShot(0, _deferred_rt_atom)
            self._invoke_drag_handlers(
                "move",
                [atom_id],
                {atom_id: new_world},
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)

    def _do_realtime_group_translate(
        self, mw: Any, move_group_dialog: Any, current_pos: Any
    ) -> None:
        """Translate group atoms and redraw during mouse-move."""
        try:
            renderer = mw.view_3d_manager.plotter.renderer
            drag_atom_idx = move_group_dialog.drag_atom_idx_vtk
            drag_initial_pos = move_group_dialog.initial_positions[drag_atom_idx]

            depth_z = self._world_to_display_depth(
                renderer, drag_initial_pos[0], drag_initial_pos[1], drag_initial_pos[2]
            )
            if depth_z is None:
                return
            new_world = self._display_to_world(
                renderer, current_pos[0], current_pos[1], depth_z
            )
            if new_world is None:
                return

            conf = mw.view_3d_manager.current_mol.GetConformer()
            translation = np.array(new_world) - drag_initial_pos
            new_positions: Dict[int, Tuple[float, float, float]] = {}
            for atom_idx in move_group_dialog.group_atoms:
                initial_pos = move_group_dialog.initial_positions[atom_idx]
                new_pos = initial_pos + translation
                mw.view_3d_manager.atom_positions_3d[atom_idx] = new_pos
                conf.SetAtomPosition(
                    atom_idx,
                    Geometry.Point3D(
                        float(new_pos[0]), float(new_pos[1]), float(new_pos[2])
                    ),
                )
                new_positions[atom_idx] = (
                    float(new_pos[0]),
                    float(new_pos[1]),
                    float(new_pos[2]),
                )

            _rt_mol = mw.view_3d_manager.current_mol
            _rt_dlg = move_group_dialog

            def _deferred_rt_grp() -> None:
                try:
                    mw.view_3d_manager.draw_molecule_3d(_rt_mol)
                    _rt_dlg.show_atom_labels()
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    logging.debug("Suppressed non-critical error", exc_info=True)
                finally:
                    self._finish_drag_redraw()

            self._drag_redraw_pending = True
            QTimer.singleShot(0, _deferred_rt_grp)
            self._invoke_drag_handlers(
                "move", list(move_group_dialog.group_atoms), new_positions
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)

    def _compute_delta_rotation_matrix(
        self, renderer: Any, start_pos: tuple, current_pos: tuple
    ) -> Optional[np.ndarray]:
        """Compute 3D rotation matrix from screen mouse displacement (dx, dy)."""
        try:
            dx = current_pos[0] - start_pos[0]
            dy = current_pos[1] - start_pos[1]
            if abs(dx) < 1e-3 and abs(dy) < 1e-3:
                return None

            camera = renderer.GetActiveCamera()
            view_up = np.array(camera.GetViewUp())
            v_dir = np.array(camera.GetDirectionOfProjection())
            right = np.cross(v_dir, view_up)

            right_norm = np.linalg.norm(right)
            if right_norm > 1e-6:
                right = right / right_norm

            view_up_norm = np.linalg.norm(view_up)
            if view_up_norm > 1e-6:
                view_up = view_up / view_up_norm

            speed = _GROUP_ROTATION_RADIANS_PER_PIXEL * self._rotation_sensitivity()
            rot_y = _rodrigues_matrix(view_up, float(dx) * speed)
            rot_x = _rodrigues_matrix(right, -float(dy) * speed)
            return rot_x @ rot_y
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
            return None

    def _get_group_rotation_matrix(
        self, mw: Any, move_group_dialog: Any, current_pos: Any
    ) -> Optional[np.ndarray]:
        """Compute group rotation matrix using delta or follow-mouse method."""
        try:
            renderer = mw.view_3d_manager.plotter.renderer
            centroid = move_group_dialog.group_centroid
            start_pos = getattr(move_group_dialog, "rotation_start_pos", None)
            if start_pos is None:
                return None

            follow_mouse = bool(
                mw.get_settings().get("rotate_group_follow_mouse", False)
            )

            if follow_mouse:
                if not hasattr(move_group_dialog, "rotation_atom_idx"):
                    move_group_dialog.rotation_atom_idx = next(
                        iter(move_group_dialog.group_atoms)
                    )
                grabbed_atom_idx = move_group_dialog.rotation_atom_idx
                grabbed_initial_pos = move_group_dialog.initial_positions[
                    grabbed_atom_idx
                ]

                depth_z = self._world_to_display_depth(
                    renderer,
                    grabbed_initial_pos[0],
                    grabbed_initial_pos[1],
                    grabbed_initial_pos[2],
                )
                if depth_z is None:
                    return None
                target_world = self._display_to_world(
                    renderer, current_pos[0], current_pos[1], depth_z
                )
                if target_world is None:
                    return None
                target_pos = np.array(target_world)

                v1 = grabbed_initial_pos - centroid
                v2 = target_pos - centroid
                v1_norm = np.linalg.norm(v1)
                v2_norm = np.linalg.norm(v2)
                if v1_norm < 1e-6 or v2_norm < 1e-6:
                    return None

                v1_n = v1 / v1_norm
                v2_n = v2 / v2_norm
                rotation_axis = np.cross(v1_n, v2_n)
                axis_norm = np.linalg.norm(rotation_axis)
                if axis_norm < 1e-6:
                    return None
                rotation_axis /= axis_norm

                cos_angle = float(np.clip(np.dot(v1_n, v2_n), -1.0, 1.0))
                return _rodrigues_matrix(rotation_axis, float(np.arccos(cos_angle)))
            return self._compute_delta_rotation_matrix(renderer, start_pos, current_pos)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)
            return None

    def _do_realtime_group_rotate(
        self, mw: Any, move_group_dialog: Any, current_pos: Any
    ) -> None:
        """Rotate group atoms around their centroid and redraw during mouse-move."""
        try:
            rot_matrix = self._get_group_rotation_matrix(
                mw, move_group_dialog, current_pos
            )
            if rot_matrix is None:
                return

            centroid = move_group_dialog.group_centroid

            conf = mw.view_3d_manager.current_mol.GetConformer()
            new_positions: Dict[int, Tuple[float, float, float]] = {}
            for atom_idx in move_group_dialog.group_atoms:
                initial_pos = move_group_dialog.initial_positions[atom_idx]
                new_pos = rot_matrix @ (initial_pos - centroid) + centroid
                mw.view_3d_manager.atom_positions_3d[atom_idx] = new_pos
                conf.SetAtomPosition(
                    atom_idx,
                    Geometry.Point3D(
                        float(new_pos[0]), float(new_pos[1]), float(new_pos[2])
                    ),
                )
                new_positions[atom_idx] = (
                    float(new_pos[0]),
                    float(new_pos[1]),
                    float(new_pos[2]),
                )

            _rt_mol = mw.view_3d_manager.current_mol
            _rt_dlg = move_group_dialog

            def _deferred_rt_rot() -> None:
                try:
                    mw.view_3d_manager.draw_molecule_3d(_rt_mol)
                    _rt_dlg.show_atom_labels()
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    logging.debug("Suppressed non-critical error", exc_info=True)
                finally:
                    self._finish_drag_redraw()

            self._drag_redraw_pending = True
            QTimer.singleShot(0, _deferred_rt_rot)
            self._invoke_drag_handlers(
                "move", list(move_group_dialog.group_atoms), new_positions
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            logging.debug("Suppressed non-critical error", exc_info=True)

    def on_left_button_up(self, obj: Any, event: Any) -> None:
        """
        Handle click release and reset state.
        """
        mw = self.main_window

        move_group_dialog = self._find_move_dialog()
        # Prevent multi-click issues
        if move_group_dialog:
            if getattr(
                move_group_dialog, "is_dragging_group_vtk", False
            ) and not getattr(move_group_dialog, "mouse_moved_vtk", False):
                # No drag = click only -> toggle
                clicked_atom = getattr(move_group_dialog, "drag_atom_idx_vtk", None)
                if clicked_atom is not None:
                    try:
                        move_group_dialog.on_atom_picked(clicked_atom)
                    except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                        logging.warning(f"Error in toggle: {e}")
                # Reset if multi-clicked without drag
                move_group_dialog.is_dragging_group_vtk = False
                move_group_dialog.drag_start_pos_vtk = None
                move_group_dialog.mouse_moved_vtk = False
                if hasattr(move_group_dialog, "initial_positions"):  # Safe
                    delattr(move_group_dialog, "initial_positions")

        if move_group_dialog and getattr(
            move_group_dialog, "is_dragging_group_vtk", False
        ):
            if getattr(move_group_dialog, "mouse_moved_vtk", False):
                # Update coordinates on release if dragged
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.view_3d_manager.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.view_3d_manager.current_mol.GetConformer()

                    # Initial position of dragged atom
                    drag_atom_initial_pos = move_group_dialog.initial_positions[
                        move_group_dialog.drag_atom_idx_vtk
                    ]

                    # screen to world conversion
                    renderer.SetWorldPoint(
                        drag_atom_initial_pos[0],
                        drag_atom_initial_pos[1],
                        drag_atom_initial_pos[2],
                        1.0,
                    )
                    renderer.WorldToDisplay()
                    display_coords = renderer.GetDisplayPoint()

                    new_display_pos = (
                        current_pos[0],
                        current_pos[1],
                        display_coords[2],
                    )
                    renderer.SetDisplayPoint(
                        new_display_pos[0], new_display_pos[1], new_display_pos[2]
                    )
                    renderer.DisplayToWorld()
                    new_world_coords = renderer.GetWorldPoint()

                    # Translation vector
                    translation_vector = np.array(
                        [
                            new_world_coords[0] - drag_atom_initial_pos[0],
                            new_world_coords[1] - drag_atom_initial_pos[1],
                            new_world_coords[2] - drag_atom_initial_pos[2],
                        ]
                    )

                    # Move entire group
                    for atom_idx in move_group_dialog.group_atoms:
                        initial_pos = move_group_dialog.initial_positions[atom_idx]
                        new_pos = initial_pos + translation_vector
                        conf.SetAtomPosition(
                            atom_idx,
                            Geometry.Point3D(
                                float(new_pos[0]), float(new_pos[1]), float(new_pos[2])
                            ),
                        )
                        mw.view_3d_manager.atom_positions_3d[atom_idx] = new_pos

                    # Defer the full redraw out of the VTK observer callback to
                    # prevent re-entrant plotter.render() calls that can deadlock
                    # the VTK render window on some OS/driver combinations.
                    _grp_mol = mw.view_3d_manager.current_mol
                    _grp_dlg = move_group_dialog

                    def _deferred_group_redraw() -> None:
                        mw.view_3d_manager.draw_molecule_3d(_grp_mol)
                        mw.view_3d_manager.update_chiral_labels()
                        _grp_dlg.show_atom_labels()
                        mw.edit_actions_manager.push_undo_state()

                    QTimer.singleShot(0, _deferred_group_redraw)
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    logging.warning("Error finalizing group drag", exc_info=True)
            else:
                # No drag = click only -> toggle
                clicked_atom = getattr(move_group_dialog, "drag_atom_idx_vtk", None)
                if clicked_atom is not None:
                    try:
                        move_group_dialog.on_atom_picked(clicked_atom)
                    except (AttributeError, RuntimeError, TypeError, ValueError):
                        logging.warning("Error in toggle", exc_info=True)
            self._end_drag_event()

        # Background click: deselect Move Group
        if move_group_dialog and not getattr(
            move_group_dialog, "is_dragging_group_vtk", False
        ):
            if not self._mouse_moved_during_drag and self._mouse_press_pos is not None:
                # Background click: deselect
                def _deferred_clear_move_group() -> None:
                    try:
                        move_group_dialog.group_atoms.clear()
                        move_group_dialog.selected_atoms.clear()
                        move_group_dialog.clear_atom_labels()
                        move_group_dialog.update_display()
                    except (AttributeError, RuntimeError):
                        # Safe defensive fallback catching AttributeError, RuntimeError
                        logging.debug("Suppressed non-critical error", exc_info=True)

                QTimer.singleShot(0, _deferred_clear_move_group)

        # Measurement mode click handling
        if (
            mw.edit_3d_manager.measurement_mode
            and not self._mouse_moved_during_drag
            and self._mouse_press_pos is not None
        ):
            # Background click -> clear selection
            def _deferred_clear_measurement() -> None:
                try:
                    mw.edit_3d_manager.clear_measurement_selection()
                except (AttributeError, RuntimeError):
                    # Safe defensive fallback catching AttributeError, RuntimeError
                    logging.debug("Suppressed non-critical error", exc_info=True)

            QTimer.singleShot(0, _deferred_clear_measurement)

        if self._is_dragging_atom:
            # Finalize custom drag
            if self.is_dragging:
                if (
                    mw.view_3d_manager.current_mol
                    and mw.view_3d_manager.current_mol.GetNumConformers() > 0
                ):
                    try:
                        atom_id = None
                        try:
                            atom_id = (
                                mw.dragged_atom_info.get("id")
                                if mw.dragged_atom_info
                                else None
                            )
                        except (AttributeError, KeyError, TypeError, ValueError):
                            atom_id = None

                        if atom_id is not None:
                            try:
                                interactor = self.GetInteractor()
                                renderer = mw.view_3d_manager.plotter.renderer
                                current_display_pos = interactor.GetEventPosition()
                                conf = mw.view_3d_manager.current_mol.GetConformer()
                                pos_3d = conf.GetAtomPosition(atom_id)
                                renderer.SetWorldPoint(
                                    pos_3d.x, pos_3d.y, pos_3d.z, 1.0
                                )
                                renderer.WorldToDisplay()
                                display_coords = renderer.GetDisplayPoint()
                                new_display_pos = (
                                    current_display_pos[0],
                                    current_display_pos[1],
                                    display_coords[2],
                                )
                                renderer.SetDisplayPoint(
                                    new_display_pos[0],
                                    new_display_pos[1],
                                    new_display_pos[2],
                                )
                                renderer.DisplayToWorld()
                                new_world_coords_tuple = renderer.GetWorldPoint()
                                new_world_coords = list(new_world_coords_tuple)[:3]
                                # Ensure container supports assignment
                                if isinstance(
                                    mw.view_3d_manager.atom_positions_3d,
                                    (list, np.ndarray),
                                ) and atom_id < len(
                                    mw.view_3d_manager.atom_positions_3d
                                ):
                                    try:
                                        mw.view_3d_manager.atom_positions_3d[
                                            atom_id
                                        ] = new_world_coords
                                    except (
                                        AttributeError,
                                        RuntimeError,
                                        ValueError,
                                        TypeError,
                                    ):  # [VTK SYNC] atom_positions_3d update may race with VTK teardown; skip safely.
                                        logging.debug(
                                            "Suppressed non-critical error",
                                            exc_info=True,
                                        )
                            except (
                                AttributeError,
                                RuntimeError,
                                TypeError,
                                ValueError,
                            ):  # [VTK SYNC] Outer drag-loop coordinate sync may race with VTK teardown; skip safely.
                                logging.debug(
                                    "Suppressed non-critical error", exc_info=True
                                )
                        conf = mw.view_3d_manager.current_mol.GetConformer()
                        pos_count = (
                            len(mw.view_3d_manager.atom_positions_3d)
                            if isinstance(
                                mw.view_3d_manager.atom_positions_3d, (list, np.ndarray)
                            )
                            else 0
                        )
                        for i in range(
                            min(mw.view_3d_manager.current_mol.GetNumAtoms(), pos_count)
                        ):
                            pos = mw.view_3d_manager.atom_positions_3d[i]
                            conf.SetAtomPosition(
                                i,
                                Geometry.Point3D(
                                    float(pos[0]), float(pos[1]), float(pos[2])
                                ),
                            )
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        logging.warning(
                            "Caught exception in " + __file__, exc_info=True
                        )

                    # Defer the redraw + undo push out of the VTK observer
                    # callback to prevent re-entrant plotter.render() deadlock.
                    _atom_mol = mw.view_3d_manager.current_mol

                    def _deferred_atom_redraw() -> None:
                        try:
                            mw.view_3d_manager.draw_molecule_3d(_atom_mol)
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            logging.warning(
                                "Caught exception in " + __file__, exc_info=True
                            )
                        mw.edit_actions_manager.push_undo_state()

                    QTimer.singleShot(0, _deferred_atom_redraw)

            self._end_drag_event()
            mw.dragged_atom_info = None
            self._stop_vtk_left_button_state()

            if self.is_dragging:
                # Defer all display-update calls that internally call plotter.render()
                # out of the VTK LeftButtonReleaseEvent callback to avoid re-entrant
                # rendering, which deadlocks the VTK render window on some platforms.
                _update_calls = [
                    mw.edit_3d_manager.update_3d_selection_display,
                    mw.edit_3d_manager.update_measurement_labels_display,
                    mw.edit_3d_manager.update_2d_measurement_labels,
                    mw.view_3d_manager.show_all_atom_info,
                ]

                def _deferred_updates() -> None:
                    for fn in _update_calls:
                        try:
                            fn()
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            logging.warning(
                                "Caught exception in " + __file__, exc_info=True
                            )

                QTimer.singleShot(0, _deferred_updates)
        else:
            if self._suppress_next_left_button_up:
                self._stop_vtk_left_button_state()
            else:
                # Delegate cleanup to parent
                super().OnLeftButtonUp()

        # Handle click release and reset state.
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_press_pos = None
        self._mouse_moved_during_drag = False
        self._suppress_next_left_button_up = False

        # Clear Move Group state
        try:
            if move_group_dialog:
                move_group_dialog.is_dragging_group_vtk = False
                move_group_dialog.drag_start_pos_vtk = None
                move_group_dialog.mouse_moved_vtk = False
                if hasattr(move_group_dialog, "initial_positions"):  # Safe
                    delattr(move_group_dialog, "initial_positions")
                if hasattr(move_group_dialog, "drag_atom_idx_vtk"):  # Safe
                    delattr(move_group_dialog, "drag_atom_idx_vtk")
        except (AttributeError, RuntimeError, ValueError, TypeError):
            logging.warning("Caught exception in " + __file__, exc_info=True)

        # Update cursor after release
        try:
            mw.view_3d_manager.plotter.setCursor(Qt.CursorShape.ArrowCursor)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            logging.warning("Caught exception in " + __file__, exc_info=True)

        # Restore focus to 2D view
        if mw and mw.init_manager.view_2d:
            mw.init_manager.view_2d.setFocus()

    def on_right_button_up(self, obj: Any, event: Any) -> None:
        """
        Finalize group rotation on right-click release.
        """
        mw = self.main_window

        move_group_dialog = self._find_move_dialog()
        if move_group_dialog and getattr(
            move_group_dialog, "is_rotating_group_vtk", False
        ):
            # Maintain selection on rotate release
            if getattr(move_group_dialog, "rotation_mouse_moved", False):
                # Apply rotation on release if moved
                try:
                    interactor = self.GetInteractor()
                    current_pos = interactor.GetEventPosition()
                    conf = mw.view_3d_manager.current_mol.GetConformer()
                    centroid = move_group_dialog.group_centroid

                    rot_matrix = self._get_group_rotation_matrix(
                        mw, move_group_dialog, current_pos
                    )
                    if rot_matrix is not None:
                        for atom_idx in move_group_dialog.group_atoms:
                            initial_pos = move_group_dialog.initial_positions[atom_idx]
                            relative_pos = initial_pos - centroid
                            rotated_pos = rot_matrix @ relative_pos
                            new_pos = rotated_pos + centroid

                            conf.SetAtomPosition(
                                atom_idx,
                                Geometry.Point3D(
                                    float(new_pos[0]),
                                    float(new_pos[1]),
                                    float(new_pos[2]),
                                ),
                            )
                            mw.view_3d_manager.atom_positions_3d[atom_idx] = new_pos

                        mw.view_3d_manager.draw_molecule_3d(
                            mw.view_3d_manager.current_mol
                        )
                        mw.view_3d_manager.update_chiral_labels()
                        move_group_dialog.show_atom_labels()
                        mw.edit_actions_manager.push_undo_state()
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    logging.warning("Error finalizing group rotation", exc_info=True)
            self._end_drag_event()

            # Reset state
            move_group_dialog.is_rotating_group_vtk = False
            move_group_dialog.rotation_start_pos = None
            move_group_dialog.rotation_mouse_moved = False
            if hasattr(move_group_dialog, "initial_positions"):  # Safe
                delattr(move_group_dialog, "initial_positions")
            if hasattr(move_group_dialog, "group_centroid"):  # Safe
                delattr(move_group_dialog, "group_centroid")
            if hasattr(move_group_dialog, "rotation_atom_idx"):  # Safe
                delattr(move_group_dialog, "rotation_atom_idx")

            return

        # Standard right-click release
        super().OnRightButtonUp()
