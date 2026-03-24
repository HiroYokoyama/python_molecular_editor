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
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera  # pylint: disable=no-name-in-module

try:
    from .constants import pt
except ImportError:
    from modules.constants import pt
try:
    from .move_group_dialog import MoveGroupDialog
except ImportError:
    from modules.move_group_dialog import MoveGroupDialog


class CustomInteractorStyle(vtkInteractorStyleTrackballCamera):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        # Custom state flags
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None

        self.AddObserver("LeftButtonPressEvent", self.on_left_button_down)
        # self.AddObserver("LeftButtonDoubleClickEvent", self.on_left_button_down)
        self.AddObserver("RightButtonPressEvent", self.on_right_button_down)
        self.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.AddObserver("LeftButtonReleaseEvent", self.on_left_button_up)
        self.AddObserver("RightButtonReleaseEvent", self.on_right_button_up)

    def on_left_button_down(self, obj, event):
        """
        Dispatch click events.
        Use custom action if atom handles, else camera rotation.
        """
        mw = self.main_window

        # Clear previous drag state
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None

        # Check Move Group dialog
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        if move_group_dialog and move_group_dialog.group_atoms:
            # Group drag if selected
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(
                    mw.atom_positions_3d - picked_position, axis=1
                )
                closest_atom_idx = np.argmin(distances)

                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        try:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            if vdw_radius < 0.1:
                                vdw_radius = 1.5
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            vdw_radius = 1.5
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)

            # If an atom in the group is clicked
            if clicked_atom_idx is not None:
                if clicked_atom_idx in move_group_dialog.group_atoms:
                    # Preparation for group drag
                    move_group_dialog._is_dragging_group_vtk = True
                    move_group_dialog._drag_atom_idx = clicked_atom_idx
                    move_group_dialog._drag_start_pos = click_pos
                    move_group_dialog._mouse_moved = False
                    # Save initial positions
                    move_group_dialog._initial_positions = {}
                    conf = mw.current_mol.GetConformer()
                    for atom_idx in move_group_dialog.group_atoms:
                        pos = conf.GetAtomPosition(atom_idx)
                        move_group_dialog._initial_positions[atom_idx] = np.array(
                            [pos.x, pos.y, pos.z]
                        )
                    mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                    return  # Disable camera rotation
                else:
                    # Clicked outside group - Search connected component
                    visited = set()
                    queue = [clicked_atom_idx]
                    visited.add(clicked_atom_idx)

                    while queue:
                        current_idx = queue.pop(0)
                        for bond_idx in range(mw.current_mol.GetNumBonds()):
                            bond = mw.current_mol.GetBondWithIdx(bond_idx)
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
                        QApplication.keyboardModifiers()
                        & Qt.KeyboardModifier.ControlModifier
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
                    move_group_dialog.show_atom_labels()
                    move_group_dialog.update_display()
                    return
            else:
                # Track mouse event to distinguish rotation from click
                self._mouse_press_pos = self.GetInteractor().GetEventPosition()
                self._mouse_moved_during_drag = False

                # Allow camera rotation
                super(CustomInteractorStyle, self).OnLeftButtonDown()
                return

        is_temp_mode = bool(
            QApplication.keyboardModifiers() & Qt.KeyboardModifier.AltModifier
        )
        is_edit_active = mw.is_3d_edit_mode or is_temp_mode

        # Ctrl+Click for atom selection (3D edit)
        is_ctrl_click = bool(
            QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier
        )

        # Handle measurement mode
        if mw.measurement_mode and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            self._mouse_moved_during_drag = False

            picker = mw.plotter.picker

            # Run pick process
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            # Special handling if atom clicked
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(
                    mw.atom_positions_3d - picked_position, axis=1
                )
                closest_atom_idx = np.argmin(distances)

                # Add range check
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # Check click threshold
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        try:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            if vdw_radius < 0.1:
                                vdw_radius = 1.5
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            vdw_radius = 1.5
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            mw.handle_measurement_atom_selection(int(closest_atom_idx))
                            return  # Selection complete, disable camera rotation

            # Clear measurement if not dragging
            self._is_dragging_atom = False
            self._mouse_press_pos = click_pos
            super().OnLeftButtonDown()
            return

        # Handle selection if 3D mol exists
        if is_edit_active and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(
                    mw.atom_positions_3d - picked_position, axis=1
                )
                closest_atom_idx = np.argmin(distances)

                # Add range check
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # Get atom safely from RDKit Mol
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        try:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            if vdw_radius < 0.1:
                                vdw_radius = 1.5
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            vdw_radius = 1.5
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            # Successfully grabbed atom
                            self._is_dragging_atom = True
                        self.is_dragging = False
                        mw.dragged_atom_info = {"id": int(closest_atom_idx)}
                        mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                        return  # Prevent camera rotation

        self._is_dragging_atom = False
        super().OnLeftButtonDown()

    def on_right_button_down(self, obj, event):
        """
        Right-click: Start group rotation if dialog open.
        """
        mw = self.main_window

        # Check if Move Group dialog is open
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        if move_group_dialog and move_group_dialog.group_atoms:
            # Start rotation drag if group selected
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(
                    mw.atom_positions_3d - picked_position, axis=1
                )
                closest_atom_idx = np.argmin(distances)

                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        try:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            if vdw_radius < 0.1:
                                vdw_radius = 1.5
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            vdw_radius = 1.5
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)

            # Start rotation drag if atom inside group clicked
            if (
                clicked_atom_idx is not None
                and clicked_atom_idx in move_group_dialog.group_atoms
            ):
                move_group_dialog._is_rotating_group_vtk = True
                move_group_dialog._rotation_start_pos = click_pos
                move_group_dialog._rotation_mouse_moved = False
                move_group_dialog._rotation_atom_idx = (
                    clicked_atom_idx  # Record grabbed atom
                )

                # Save initial positions and centroid
                move_group_dialog._initial_positions = {}
                conf = mw.current_mol.GetConformer()
                centroid = np.zeros(3)
                for atom_idx in move_group_dialog.group_atoms:
                    pos = conf.GetAtomPosition(atom_idx)
                    pos_array = np.array([pos.x, pos.y, pos.z])
                    move_group_dialog._initial_positions[atom_idx] = pos_array
                    centroid += pos_array
                centroid /= len(move_group_dialog.group_atoms)
                move_group_dialog._group_centroid = centroid

                mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                return  # Disable camera rotation

        # Standard right-click
        super().OnRightButtonDown()

    def on_mouse_move(self, obj, event):
        """
        Handle mouse move (drag vs camera/hover).
        """
        mw = self.main_window

        # Move Group drag handling
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        if move_group_dialog and getattr(
            move_group_dialog, "_is_dragging_group_vtk", False
        ):
            # Dragging group - record offset
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()

            dx = current_pos[0] - move_group_dialog._drag_start_pos[0]
            dy = current_pos[1] - move_group_dialog._drag_start_pos[1]

            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._mouse_moved = True

            return  # Disable camera rotation

        # Group rotation handling
        if move_group_dialog and getattr(
            move_group_dialog, "_is_rotating_group_vtk", False
        ):
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()

            dx = current_pos[0] - move_group_dialog._rotation_start_pos[0]
            dy = current_pos[1] - move_group_dialog._rotation_start_pos[1]

            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._rotation_mouse_moved = True

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
            atom_id = mw.dragged_atom_info["id"]
        else:
            # Delegate camera rotation to parent
            super().OnMouseMove()

            # Update cursor display
            is_edit_active = mw.is_3d_edit_mode or interactor.GetAltKey()
            if is_edit_active:
                # Hover check if edit active
                atom_under_cursor = False
                click_pos = interactor.GetEventPosition()
                picker = mw.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
                if picker.GetActor() is mw.atom_actor:
                    atom_under_cursor = True

                if atom_under_cursor:
                    mw.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                else:
                    mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            else:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)

    def on_left_button_up(self, obj, event):
        """
        Handle click release and reset state.
        """
        mw = self.main_window

        # Finalize Move Group drag
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        # Prevent multi-click issues
        if move_group_dialog:
            if getattr(
                move_group_dialog, "_is_dragging_group_vtk", False
            ) and not getattr(move_group_dialog, "_mouse_moved", False):
                # Reset if multi-clicked without drag
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, "_initial_positions"):
                    delattr(move_group_dialog, "_initial_positions")

        if move_group_dialog and getattr(
            move_group_dialog, "_is_dragging_group_vtk", False
        ):
            if getattr(move_group_dialog, "_mouse_moved", False):
                # Update coordinates on release if dragged
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()

                    # Initial position of dragged atom
                    drag_atom_initial_pos = move_group_dialog._initial_positions[
                        move_group_dialog._drag_atom_idx
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
                        initial_pos = move_group_dialog._initial_positions[atom_idx]
                        new_pos = initial_pos + translation_vector
                        conf.SetAtomPosition(atom_idx, new_pos.tolist())
                        mw.atom_positions_3d[atom_idx] = new_pos

                    # Update 3D display
                    mw.draw_molecule_3d(mw.current_mol)
                    mw.update_chiral_labels()
                    move_group_dialog.show_atom_labels()
                    mw.push_undo_state()
                except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                    print(f"Error finalizing group drag: {e}")
            else:
                # No drag = click only -> toggle
                if hasattr(move_group_dialog, "_drag_atom_idx"):
                    clicked_atom = move_group_dialog._drag_atom_idx
                    try:
                        move_group_dialog.on_atom_picked(clicked_atom)
                    except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                        print(f"Error in toggle: {e}")

        # Background click: deselect
        if move_group_dialog and not getattr(
            move_group_dialog, "_is_dragging_group_vtk", False
        ):
            if not self._mouse_moved_during_drag and self._mouse_press_pos is not None:
                # Background click: deselect
                move_group_dialog.group_atoms.clear()
                move_group_dialog.selected_atoms.clear()
                move_group_dialog.clear_atom_labels()
                move_group_dialog.update_display()

        # Measurement mode click handling
        if (
            mw.measurement_mode
            and not self._mouse_moved_during_drag
            and self._mouse_press_pos is not None
        ):
            # Background click -> clear selection
            mw.clear_measurement_selection()

        if self._is_dragging_atom:
            # Finalize custom drag
            if self.is_dragging:
                if mw.current_mol and mw.current_mol.GetNumConformers() > 0:
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
                                renderer = mw.plotter.renderer
                                current_display_pos = interactor.GetEventPosition()
                                conf = mw.current_mol.GetConformer()
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
                                try:
                                    mw.atom_positions_3d[atom_id] = new_world_coords
                                except (AttributeError, KeyError, TypeError, ValueError, IndexError):
                                    try:
                                        ap = list(mw.atom_positions_3d)
                                        ap[atom_id] = new_world_coords
                                        mw.atom_positions_3d = ap
                                    except (AttributeError, RuntimeError):  # pragma: no cover
                                        import traceback
                                        traceback.print_exc()
                            except (AttributeError, RuntimeError, TypeError, ValueError):  # pragma: no cover
                                import traceback
                                traceback.print_exc()
                        conf = mw.current_mol.GetConformer()
                        for i in range(mw.current_mol.GetNumAtoms()):
                            try:
                                pos = mw.atom_positions_3d[i]
                                conf.SetAtomPosition(i, pos.tolist())
                            except (AttributeError, KeyError, TypeError, ValueError):
                                # Skip individual failures but continue applying
                                # other atom positions.
                                pass
                    except (AttributeError, RuntimeError):
                        # If applying positions fails, continue to redraw from
                        # whatever authoritative state is available.
                        pass

                    # Redraw and push undo state
                    try:
                        mw.draw_molecule_3d(mw.current_mol)
                    except (AttributeError, RuntimeError):  # pragma: no cover
                        import traceback
                        traceback.print_exc()

                    mw.push_undo_state()
            mw.dragged_atom_info = None
            # Update all relevant UI displays and labels
            for update_call in [
                mw.update_3d_selection_display,
                mw.update_measurement_labels_display,
                mw.update_2d_measurement_labels,
                mw.show_all_atom_info,
            ]:
                try:
                    update_call()
                except (AttributeError, RuntimeError):  # pragma: no cover
                    pass
        else:
            # Delegate cleanup to parent
            super().OnLeftButtonUp()

        # Handle click release and reset state.
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_press_pos = None
        self._mouse_moved_during_drag = False

        # Clear Move Group state
        try:
            if move_group_dialog:
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, "_initial_positions"):
                    delattr(move_group_dialog, "_initial_positions")
                if hasattr(move_group_dialog, "_drag_atom_idx"):
                    delattr(move_group_dialog, "_drag_atom_idx")
        except (AttributeError, RuntimeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        # Update cursor after release
        try:
            mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
        except (AttributeError, RuntimeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        # Restore focus to 2D view
        if mw and mw.view_2d:
            mw.view_2d.setFocus()

    def on_right_button_up(self, obj, event):
        """
        Finalize group rotation on right-click release.
        """
        mw = self.main_window

        # Finalize Move Group rotation
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except (AttributeError, RuntimeError, TypeError):  # pragma: no cover
            import traceback
            traceback.print_exc()

        if move_group_dialog and getattr(
            move_group_dialog, "_is_rotating_group_vtk", False
        ):
            # Maintain selection on rotate release
            if getattr(move_group_dialog, "_rotation_mouse_moved", False):
                # Apply rotation on release if moved
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()
                    centroid = move_group_dialog._group_centroid

                    # Save initial grabbed atom index
                    if not hasattr(move_group_dialog, "_rotation_atom_idx"):
                        move_group_dialog._rotation_atom_idx = next(
                            iter(move_group_dialog.group_atoms)
                        )

                    grabbed_atom_idx = move_group_dialog._rotation_atom_idx
                    grabbed_initial_pos = move_group_dialog._initial_positions[
                        grabbed_atom_idx
                    ]

                    # Get start screen coordinates
                    renderer.SetWorldPoint(
                        grabbed_initial_pos[0],
                        grabbed_initial_pos[1],
                        grabbed_initial_pos[2],
                        1.0,
                    )
                    renderer.WorldToDisplay()
                    start_display = renderer.GetDisplayPoint()

                    # Convert current mouse pos to world (same depth)
                    renderer.SetDisplayPoint(
                        current_pos[0], current_pos[1], start_display[2]
                    )
                    renderer.DisplayToWorld()
                    target_world = renderer.GetWorldPoint()
                    target_pos = np.array(
                        [target_world[0], target_world[1], target_world[2]]
                    )

                    # Vectors relative to centroid
                    v1 = grabbed_initial_pos - centroid
                    v2 = target_pos - centroid

                    # Normalize vectors
                    v1_norm = np.linalg.norm(v1)
                    v2_norm = np.linalg.norm(v2)

                    if v1_norm > 1e-6 and v2_norm > 1e-6:
                        v1_normalized = v1 / v1_norm
                        v2_normalized = v2 / v2_norm

                        # Rotation axis
                        rotation_axis = np.cross(v1_normalized, v2_normalized)
                        axis_norm = np.linalg.norm(rotation_axis)

                        if axis_norm > 1e-6:
                            rotation_axis = rotation_axis / axis_norm

                            # Rotation angle
                            cos_angle = np.clip(
                                np.dot(v1_normalized, v2_normalized), -1.0, 1.0
                            )
                            angle = np.arccos(cos_angle)

                            # Create rotation matrix
                            K = np.array(
                                [
                                    [0, -rotation_axis[2], rotation_axis[1]],
                                    [rotation_axis[2], 0, -rotation_axis[0]],
                                    [-rotation_axis[1], rotation_axis[0], 0],
                                ]
                            )

                            rot_matrix = (
                                np.eye(3)
                                + np.sin(angle) * K
                                + (1 - np.cos(angle)) * (K @ K)
                            )

                            # Rotate group around centroid
                            for atom_idx in move_group_dialog.group_atoms:
                                initial_pos = move_group_dialog._initial_positions[
                                    atom_idx
                                ]
                                # Relative position from centroid
                                relative_pos = initial_pos - centroid
                                # Apply rotation
                                rotated_pos = rot_matrix @ relative_pos
                                # Restore absolute position
                                new_pos = rotated_pos + centroid

                                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                                mw.atom_positions_3d[atom_idx] = new_pos

                            # Update 3D display
                            mw.draw_molecule_3d(mw.current_mol)
                            mw.update_chiral_labels()
                            move_group_dialog.show_atom_labels()
                            mw.push_undo_state()
                except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                    print(f"Error finalizing group rotation: {e}")

            # Reset state
            move_group_dialog._is_rotating_group_vtk = False
            move_group_dialog._rotation_start_pos = None
            move_group_dialog._rotation_mouse_moved = False
            if hasattr(move_group_dialog, "_initial_positions"):
                delattr(move_group_dialog, "_initial_positions")
            if hasattr(move_group_dialog, "_group_centroid"):
                delattr(move_group_dialog, "_group_centroid")
            if hasattr(move_group_dialog, "_rotation_atom_idx"):
                delattr(move_group_dialog, "_rotation_atom_idx")

            try:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except (AttributeError, RuntimeError):  # pragma: no cover
                import traceback
                traceback.print_exc()

            return

        # Standard right-click release
        super().OnRightButtonUp()
