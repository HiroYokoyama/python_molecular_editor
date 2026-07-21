"""Unit tests for CustomInteractorStyle 3D event handling."""

import numpy as np
import pytest

from unittest.mock import MagicMock, patch
from moleditpy.ui.custom_interactor_style import CustomInteractorStyle
from PyQt6.QtCore import Qt

_VTK_BASE = "vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera"


def _style_at(host, pos=(100, 100)):
    """Style with a mocked interactor at *pos* and neutral VTK state."""
    style = CustomInteractorStyle(host)
    style.GetState = MagicMock(return_value=0)
    style.StopState = MagicMock()
    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = pos
    mock_interactor.GetAltKey.return_value = False
    style.GetInteractor = MagicMock(return_value=mock_interactor)
    return style


def _move_dialog(**flags):
    dlg = MagicMock()
    dlg.isVisible.return_value = True
    type(dlg).__name__ = "MoveGroupDialog"
    for name, value in flags.items():
        setattr(dlg, name, value)
    return dlg


def test_custom_interactor_style_left_click_atom_selection(app, mock_parser_host):
    """Verify that left-clicking an atom in 3D edit mode successfully selects it."""
    # Ensure no MoveGroupDialog interferes
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = True

    interactor_style = CustomInteractorStyle(mock_parser_host)

    # Mocking VTK methods
    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    mock_parser_host.edit_3d_manager.measurement_mode = False

    # Setup mock RDKit Mol
    mock_mol = MagicMock()
    mock_mol.GetNumAtoms.return_value = 1
    mock_atom = MagicMock()
    mock_atom.GetAtomicNum.return_value = 6  # Carbon
    mock_mol.GetAtomWithIdx.return_value = mock_atom
    mock_parser_host.view_3d_manager.current_mol = mock_mol

    # Mock VTK Picker
    mock_picker = MagicMock()
    mock_parser_host.plotter.picker = mock_picker
    mock_parser_host.view_3d_manager.atom_actor = MagicMock()
    mock_picker.GetActor.return_value = mock_parser_host.view_3d_manager.atom_actor
    mock_picker.GetPickPosition.return_value = (0, 0, 0)

    import numpy as np

    mock_parser_host.view_3d_manager.atom_positions_3d = np.array([[0.0, 0.0, 0.0]])

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=0,
        ),
        patch(
            "vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera.OnLeftButtonDown"
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier

        interactor_style.on_left_button_down(None, None)

        # Verify the atom was selected.
        # Note: reset_interactor_state() resets cursor to ArrowCursor at the start
        # of every press, so setCursor is called twice. We verify the final call
        # is ClosedHandCursor (the meaningful state after atom grab).
        assert interactor_style._is_dragging_atom is True
        mock_parser_host.plotter.setCursor.assert_called_with(
            Qt.CursorShape.ClosedHandCursor
        )
        assert mock_parser_host.dragged_atom_info["id"] == 0


def test_custom_interactor_style_background_click(app, mock_parser_host):
    """Verify that clicking the background (no atom selected) allows VTK trackball rotation."""
    mock_parser_host.is_3d_edit_mode = True
    interactor_style = CustomInteractorStyle(mock_parser_host)

    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    # Mock VTK Picker returning None (background click)
    mock_picker = MagicMock()
    mock_parser_host.plotter.picker = mock_picker
    mock_picker.GetActor.return_value = None

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera.OnLeftButtonDown"
        ) as mock_super_down,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier

        interactor_style.on_left_button_down(None, None)

        # Verify drag is false and superclass rotation method was called
        assert interactor_style._is_dragging_atom is False
        mock_super_down.assert_called_once()


def test_measurement_atom_click_suppresses_vtk_release(app, mock_parser_host):
    """Atom measurement clicks consume press and must not leak release to VTK."""
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = False
    mock_parser_host.edit_3d_manager.measurement_mode = True

    interactor_style = CustomInteractorStyle(mock_parser_host)

    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    mock_mol = MagicMock()
    mock_mol.GetNumAtoms.return_value = 1
    mock_mol.GetAtomWithIdx.return_value = MagicMock()
    mock_parser_host.view_3d_manager.current_mol = mock_mol

    mock_picker = MagicMock()
    mock_parser_host.plotter.picker = mock_picker
    mock_parser_host.view_3d_manager.atom_actor = MagicMock()
    mock_picker.GetActor.return_value = mock_parser_host.view_3d_manager.atom_actor
    mock_picker.GetPickPosition.return_value = (0, 0, 0)

    import numpy as np

    mock_parser_host.view_3d_manager.atom_positions_3d = np.array([[0.0, 0.0, 0.0]])

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=0,
        ),
        patch(
            "vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera.OnLeftButtonUp"
        ) as mock_super_up,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier

        interactor_style.on_left_button_down(None, None)
        interactor_style.on_left_button_up(None, None)

        mock_super_up.assert_not_called()


def test_atom_click_without_drag_resets_without_render_updates(app, mock_parser_host):
    """Simple atom clicks should not run drag redraw/update work."""
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = True
    mock_parser_host.edit_3d_manager.measurement_mode = False
    mock_parser_host.edit_3d_manager.update_3d_selection_display = MagicMock()
    mock_parser_host.edit_3d_manager.update_measurement_labels_display = MagicMock()
    mock_parser_host.edit_3d_manager.update_2d_measurement_labels = MagicMock()
    mock_parser_host.view_3d_manager.show_all_atom_info = MagicMock()

    interactor_style = CustomInteractorStyle(mock_parser_host)
    interactor_style.StopState = MagicMock()

    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    mock_mol = MagicMock()
    mock_mol.GetNumAtoms.return_value = 1
    mock_mol.GetAtomWithIdx.return_value = MagicMock()
    mock_parser_host.view_3d_manager.current_mol = mock_mol

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=0,
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier

        interactor_style.on_left_button_down(None, None)
        interactor_style.on_left_button_up(None, None)

    interactor_style.StopState.assert_called()
    mock_parser_host.edit_3d_manager.update_3d_selection_display.assert_not_called()
    mock_parser_host.edit_3d_manager.update_measurement_labels_display.assert_not_called()
    mock_parser_host.edit_3d_manager.update_2d_measurement_labels.assert_not_called()
    mock_parser_host.view_3d_manager.show_all_atom_info.assert_not_called()


def test_custom_interactor_style_move_selected_atoms_no_bfs(app, mock_parser_host):
    """Verify that clicking outside selection in MoveSelectedAtomsDialog toggles only that atom (no BFS)."""
    interactor_style = CustomInteractorStyle(mock_parser_host)

    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    mock_dialog = MagicMock()
    mock_dialog.isVisible.return_value = True
    type(mock_dialog).__name__ = "MoveSelectedAtomsDialog"
    mock_dialog.group_atoms = {0}

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch("moleditpy.ui.custom_interactor_style.QTimer.singleShot") as mock_timer,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=1,
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = [mock_dialog]
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier

        interactor_style.on_left_button_down(None, None)

        mock_timer.assert_called_once()
        callback = mock_timer.call_args[0][1]
        callback()

        mock_dialog.on_atom_picked.assert_called_once_with(1)


def test_heal_resets_stuck_atom_drag_when_left_button_not_held(app, mock_parser_host):
    """A lost release must not leave _is_dragging_atom stuck when no button is held."""
    interactor_style = CustomInteractorStyle(mock_parser_host)
    interactor_style.GetState = MagicMock(return_value=0)
    interactor_style.StopState = MagicMock()
    interactor_style._is_dragging_atom = True

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        interactor_style._heal_stuck_pointer_state(None)

    assert interactor_style._is_dragging_atom is False


def test_heal_resets_stuck_camera_state_when_no_button_held(app, mock_parser_host):
    """A stuck VTK ROTATE/PAN state with no button held must be stopped."""
    interactor_style = CustomInteractorStyle(mock_parser_host)
    interactor_style.GetState = MagicMock(return_value=1)  # VTKIS_ROTATE
    interactor_style.StopState = MagicMock()

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        interactor_style._heal_stuck_pointer_state(None)

    interactor_style.StopState.assert_called()


def test_heal_keeps_active_drag_while_left_button_held(app, mock_parser_host):
    """A genuine in-progress atom drag must not be reset."""
    interactor_style = CustomInteractorStyle(mock_parser_host)
    interactor_style.GetState = MagicMock(return_value=0)
    interactor_style.StopState = MagicMock()
    interactor_style._is_dragging_atom = True

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        interactor_style._heal_stuck_pointer_state(None)

    assert interactor_style._is_dragging_atom is True
    interactor_style.StopState.assert_not_called()


def test_heal_clears_stuck_move_group_drag_flags(app, mock_parser_host):
    """Stuck Move Group drag/rotate latches are cleared when buttons are not held."""
    interactor_style = CustomInteractorStyle(mock_parser_host)
    interactor_style.GetState = MagicMock(return_value=0)
    interactor_style.StopState = MagicMock()

    mock_dialog = MagicMock()
    mock_dialog.is_dragging_group_vtk = True
    mock_dialog.is_rotating_group_vtk = True

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        interactor_style._heal_stuck_pointer_state(mock_dialog)

    assert mock_dialog.is_dragging_group_vtk is False
    assert mock_dialog.is_rotating_group_vtk is False


def test_custom_interactor_style_right_click_rotation(app, mock_parser_host):
    """Verify right-click triggers group rotation and release finalizes it."""
    interactor_style = CustomInteractorStyle(mock_parser_host)

    mock_interactor = MagicMock()
    mock_interactor.GetEventPosition.return_value = (100, 100)
    interactor_style.GetInteractor = MagicMock(return_value=mock_interactor)

    mock_dialog = MagicMock()
    mock_dialog.isVisible.return_value = True
    type(mock_dialog).__name__ = "MoveGroupDialog"
    mock_dialog.group_atoms = {0, 1}
    mock_dialog.is_rotating_group_vtk = False
    mock_dialog.rotation_mouse_moved = False
    mock_dialog.rotation_start_pos = None

    # Mock RDKit Mol and Conformer
    mock_mol = MagicMock()
    mock_conformer = MagicMock()
    mock_pos0 = MagicMock(x=0.0, y=0.0, z=0.0)
    mock_pos1 = MagicMock(x=2.0, y=0.0, z=0.0)
    mock_conformer.GetAtomPosition.side_effect = (
        lambda idx: mock_pos0 if idx == 0 else mock_pos1
    )
    mock_mol.GetConformer.return_value = mock_conformer
    mock_parser_host.view_3d_manager.current_mol = mock_mol

    import numpy as np

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=0,
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = [mock_dialog]

        # 1. Test right-click press starts rotation
        interactor_style.on_right_button_down(None, None)
        assert mock_dialog.is_rotating_group_vtk is True
        assert mock_dialog.rotation_start_pos == (100, 100)
        assert mock_dialog.rotation_atom_idx == 0
        assert mock_dialog.group_centroid is not None

        # 2. Test right-click release with no movement resets state
        mock_dialog.rotation_mouse_moved = False
        interactor_style.on_right_button_up(None, None)
        assert mock_dialog.is_rotating_group_vtk is False
        assert mock_dialog.rotation_start_pos is None

        # 3. Test right-click release with movement applies rotation
        mock_dialog.is_rotating_group_vtk = True
        mock_dialog.rotation_mouse_moved = True
        mock_dialog.rotation_start_pos = (100, 100)
        mock_dialog.rotation_atom_idx = 0
        mock_dialog.initial_positions = {
            0: np.array([0.0, 0.0, 0.0]),
            1: np.array([2.0, 0.0, 0.0]),
        }
        mock_dialog.group_centroid = np.array([1.0, 0.0, 0.0])

        # Mock renderer methods to return coordinate tuples
        mock_renderer = mock_parser_host.view_3d_manager.plotter.renderer
        mock_renderer.GetDisplayPoint.return_value = (100.0, 100.0, 0.5)
        mock_renderer.GetWorldPoint.return_value = (1.0, 1.0, 0.0, 1.0)

        interactor_style.on_right_button_up(None, None)
        assert mock_dialog.is_rotating_group_vtk is False
        assert mock_dialog.rotation_start_pos is None
        mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()


# ---------------------------------------------------------------------------
# on_mouse_move
# ---------------------------------------------------------------------------


def test_mouse_move_group_drag_past_threshold_marks_moved(app, mock_parser_host):
    style = _style_at(mock_parser_host, pos=(110, 110))
    dlg = _move_dialog(
        is_dragging_group_vtk=True,
        drag_start_pos_vtk=(100, 100),
        mouse_moved_vtk=False,
    )

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove") as mock_super_move,
    ):
        mock_qapp.topLevelWidgets.return_value = [dlg]
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert dlg.mouse_moved_vtk is True
    mock_super_move.assert_not_called()  # camera rotation stays disabled


def test_mouse_move_group_drag_below_threshold_not_moved(app, mock_parser_host):
    style = _style_at(mock_parser_host, pos=(103, 103))
    dlg = _move_dialog(
        is_dragging_group_vtk=True,
        drag_start_pos_vtk=(100, 100),
        mouse_moved_vtk=False,
    )

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.topLevelWidgets.return_value = [dlg]
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert dlg.mouse_moved_vtk is False


def test_mouse_move_group_drag_without_start_pos_returns_early(app, mock_parser_host):
    style = _style_at(mock_parser_host)
    dlg = _move_dialog(
        is_dragging_group_vtk=True, drag_start_pos_vtk=None, mouse_moved_vtk=False
    )

    with patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp:
        mock_qapp.topLevelWidgets.return_value = [dlg]
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert dlg.mouse_moved_vtk is False


def test_mouse_move_group_rotation_past_threshold_marks_moved(app, mock_parser_host):
    style = _style_at(mock_parser_host, pos=(90, 120))
    dlg = _move_dialog(
        is_dragging_group_vtk=False,
        is_rotating_group_vtk=True,
        rotation_start_pos=(100, 100),
        rotation_mouse_moved=False,
    )

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove") as mock_super_move,
    ):
        mock_qapp.topLevelWidgets.return_value = [dlg]
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.RightButton
        style.on_mouse_move(None, None)

    assert dlg.rotation_mouse_moved is True
    mock_super_move.assert_not_called()


def test_mouse_move_past_press_threshold_sets_moved_during_drag(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = False
    style = _style_at(mock_parser_host, pos=(104, 104))
    style._mouse_press_pos = (100, 100)

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove"),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert style._mouse_moved_during_drag is True


def test_mouse_move_within_press_threshold_keeps_flag_clear(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = False
    style = _style_at(mock_parser_host, pos=(102, 102))
    style._mouse_press_pos = (100, 100)

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove"),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert style._mouse_moved_during_drag is False


def test_mouse_move_during_atom_drag_sets_is_dragging_and_skips_camera(
    app, mock_parser_host
):
    style = _style_at(mock_parser_host)
    style._is_dragging_atom = True
    mock_parser_host.dragged_atom_info = {"id": 0}

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove") as mock_super_move,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        style.on_mouse_move(None, None)

    assert style.is_dragging is True
    mock_super_move.assert_not_called()


def test_mouse_move_hover_over_atom_shows_open_hand(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = True
    style = _style_at(mock_parser_host)

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove") as mock_super_move,
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=0,
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        style.on_mouse_move(None, None)

    mock_super_move.assert_called_once()
    mock_parser_host.view_3d_manager.plotter.setCursor.assert_called_with(
        Qt.CursorShape.OpenHandCursor
    )


def test_mouse_move_hover_over_background_shows_arrow(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = True
    style = _style_at(mock_parser_host)

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove"),
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen",
            return_value=None,
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        style.on_mouse_move(None, None)

    mock_parser_host.view_3d_manager.plotter.setCursor.assert_called_with(
        Qt.CursorShape.ArrowCursor
    )


def test_mouse_move_outside_edit_mode_shows_arrow_without_picking(
    app, mock_parser_host
):
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = False
    style = _style_at(mock_parser_host)

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnMouseMove"),
        patch(
            "moleditpy.ui.custom_interactor_style.pick_atom_index_from_screen"
        ) as mock_pick,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        style.on_mouse_move(None, None)

    mock_pick.assert_not_called()
    mock_parser_host.view_3d_manager.plotter.setCursor.assert_called_with(
        Qt.CursorShape.ArrowCursor
    )


# ---------------------------------------------------------------------------
# on_left_button_up
# ---------------------------------------------------------------------------


def test_release_click_only_toggles_group_atom_and_resets_latch(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host)
    dlg = _move_dialog(
        is_dragging_group_vtk=True,
        mouse_moved_vtk=False,
        drag_atom_idx_vtk=3,
        initial_positions={3: np.array([0.0, 0.0, 0.0])},
    )

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnLeftButtonUp"),
    ):
        mock_qapp.topLevelWidgets.return_value = [dlg]
        style.on_left_button_up(None, None)

    dlg.on_atom_picked.assert_called_once_with(3)
    assert dlg.is_dragging_group_vtk is False
    assert dlg.drag_start_pos_vtk is None
    assert dlg.mouse_moved_vtk is False
    assert not hasattr(dlg, "initial_positions")


def test_release_after_group_drag_translates_whole_group(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host, pos=(150, 150))
    dlg = _move_dialog(
        is_dragging_group_vtk=True,
        mouse_moved_vtk=True,
        drag_atom_idx_vtk=0,
        group_atoms={0, 1},
        initial_positions={
            0: np.array([0.0, 0.0, 0.0]),
            1: np.array([2.0, 0.0, 0.0]),
        },
    )

    renderer = mock_parser_host.view_3d_manager.plotter.renderer
    renderer.GetDisplayPoint.return_value = (100.0, 100.0, 0.5)
    renderer.GetWorldPoint.return_value = (1.0, 1.0, 0.0, 1.0)

    mock_conf = MagicMock()
    mock_mol = MagicMock()
    mock_mol.GetConformer.return_value = mock_conf
    mock_parser_host.view_3d_manager.current_mol = mock_mol
    mock_parser_host.view_3d_manager.atom_positions_3d = np.zeros((2, 3))

    deferred = []
    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.QTimer.singleShot",
            side_effect=lambda _ms, fn: deferred.append(fn),
        ),
        patch(f"{_VTK_BASE}.OnLeftButtonUp"),
    ):
        mock_qapp.topLevelWidgets.return_value = [dlg]
        style.on_left_button_up(None, None)

    # Translation vector is (1, 1, 0): both atoms moved by it
    assert mock_conf.SetAtomPosition.call_count == 2
    assert np.allclose(
        mock_parser_host.view_3d_manager.atom_positions_3d,
        [[1.0, 1.0, 0.0], [3.0, 1.0, 0.0]],
    )

    # Redraw + undo push are deferred out of the VTK callback
    assert len(deferred) == 1
    deferred[0]()
    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()
    mock_parser_host.view_3d_manager.update_chiral_labels.assert_called_once()
    dlg.show_atom_labels.assert_called_once()
    mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()


def test_release_background_click_deselects_move_group_deferred(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host)
    style._mouse_press_pos = (100, 100)
    style._mouse_moved_during_drag = False
    dlg = _move_dialog(is_dragging_group_vtk=False)

    deferred = []
    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.QTimer.singleShot",
            side_effect=lambda _ms, fn: deferred.append(fn),
        ),
        patch(f"{_VTK_BASE}.OnLeftButtonUp"),
    ):
        mock_qapp.topLevelWidgets.return_value = [dlg]
        style.on_left_button_up(None, None)

    assert len(deferred) == 1
    deferred[0]()
    dlg.group_atoms.clear.assert_called_once()
    dlg.selected_atoms.clear.assert_called_once()
    dlg.clear_atom_labels.assert_called_once()
    dlg.update_display.assert_called_once()


def test_release_background_click_in_measurement_mode_clears_selection(
    app, mock_parser_host
):
    mock_parser_host.edit_3d_manager.measurement_mode = True
    style = _style_at(mock_parser_host)
    style._mouse_press_pos = (100, 100)
    style._mouse_moved_during_drag = False

    deferred = []
    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.QTimer.singleShot",
            side_effect=lambda _ms, fn: deferred.append(fn),
        ),
        patch(f"{_VTK_BASE}.OnLeftButtonUp"),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        style.on_left_button_up(None, None)

    assert len(deferred) == 1
    deferred[0]()
    mock_parser_host.edit_3d_manager.clear_measurement_selection.assert_called_once()


def test_release_after_camera_drag_does_not_clear_measurement(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = True
    style = _style_at(mock_parser_host)
    style._mouse_press_pos = (100, 100)
    style._mouse_moved_during_drag = True  # camera was rotated, not a click

    deferred = []
    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.QTimer.singleShot",
            side_effect=lambda _ms, fn: deferred.append(fn),
        ),
        patch(f"{_VTK_BASE}.OnLeftButtonUp"),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        style.on_left_button_up(None, None)

    assert deferred == []


def test_release_after_atom_drag_writes_position_and_defers_redraw(
    app, mock_parser_host
):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host, pos=(150, 150))
    style._is_dragging_atom = True
    style.is_dragging = True
    mock_parser_host.dragged_atom_info = {"id": 0}

    mol = MagicMock()
    mock_parser_host.view_3d_manager.current_mol = mol
    mol.GetNumConformers.return_value = 1
    mol.GetNumAtoms.return_value = 1
    mock_conf = MagicMock()
    mock_conf.GetAtomPosition.return_value = MagicMock(x=0.0, y=0.0, z=0.0)
    mol.GetConformer.return_value = mock_conf

    renderer = mock_parser_host.view_3d_manager.plotter.renderer
    renderer.GetDisplayPoint.return_value = (100.0, 100.0, 0.5)
    renderer.GetWorldPoint.return_value = (1.0, 1.0, 0.0, 1.0)
    mock_parser_host.view_3d_manager.atom_positions_3d = np.zeros((1, 3))

    deferred = []
    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(
            "moleditpy.ui.custom_interactor_style.QTimer.singleShot",
            side_effect=lambda _ms, fn: deferred.append(fn),
        ),
    ):
        mock_qapp.topLevelWidgets.return_value = []
        style.on_left_button_up(None, None)

    assert np.allclose(
        mock_parser_host.view_3d_manager.atom_positions_3d, [[1.0, 1.0, 0.0]]
    )
    mock_conf.SetAtomPosition.assert_called_once()
    assert mock_parser_host.dragged_atom_info is None
    assert style._is_dragging_atom is False
    assert style.is_dragging is False

    # Two deferred jobs: redraw+undo, then display updates
    assert len(deferred) == 2
    for fn in deferred:
        fn()
    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()
    mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()
    mock_parser_host.edit_3d_manager.update_3d_selection_display.assert_called_once()
    mock_parser_host.view_3d_manager.show_all_atom_info.assert_called_once()


def test_release_suppressed_after_group_grab_skips_vtk_cleanup(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host)
    style._suppress_next_left_button_up = True

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnLeftButtonUp") as mock_super_up,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        style.on_left_button_up(None, None)

    mock_super_up.assert_not_called()
    style.StopState.assert_called()
    assert style._suppress_next_left_button_up is False


def test_release_delegates_to_vtk_and_restores_focus(app, mock_parser_host):
    mock_parser_host.edit_3d_manager.measurement_mode = False
    style = _style_at(mock_parser_host)
    style._mouse_press_pos = None

    with (
        patch("moleditpy.ui.custom_interactor_style.QApplication") as mock_qapp,
        patch(f"{_VTK_BASE}.OnLeftButtonUp") as mock_super_up,
    ):
        mock_qapp.topLevelWidgets.return_value = []
        style.on_left_button_up(None, None)

    mock_super_up.assert_called_once()
    mock_parser_host.view_3d_manager.plotter.setCursor.assert_called_with(
        Qt.CursorShape.ArrowCursor
    )
    mock_parser_host.init_manager.view_2d.setFocus.assert_called_once()


# ---------------------------------------------------------------------------
# reset_interactor_state
# ---------------------------------------------------------------------------


def test_reset_interactor_state_clears_all_flags_and_cursor(app, mock_parser_host):
    style = _style_at(mock_parser_host)
    style._is_dragging_atom = True
    style.is_dragging = True
    style._mouse_moved_during_drag = True
    style._mouse_press_pos = (10, 10)
    style._suppress_next_left_button_up = True

    style.reset_interactor_state()

    assert style._is_dragging_atom is False
    assert style.is_dragging is False
    assert style._mouse_moved_during_drag is False
    assert style._mouse_press_pos is None
    assert style._suppress_next_left_button_up is False
    assert mock_parser_host.dragged_atom_info is None
    mock_parser_host.view_3d_manager.plotter.setCursor.assert_called_with(
        Qt.CursorShape.ArrowCursor
    )


def test_reset_interactor_state_survives_broken_plotter(app, mock_parser_host):
    style = _style_at(mock_parser_host)
    mock_parser_host.view_3d_manager.plotter.setCursor.side_effect = RuntimeError(
        "render window gone"
    )

    style.reset_interactor_state()  # must not raise

    assert style._is_dragging_atom is False


def test_reset_interactor_state_without_main_window(app):
    style = CustomInteractorStyle(None)
    style.StopState = MagicMock()
    style._is_dragging_atom = True

    style.reset_interactor_state()  # must not raise

    assert style._is_dragging_atom is False


def _rotation_style(sensitivity=1.0, motion_factor=10.0):
    """A style wired for _rotate_size_independent with a fake camera."""
    host = MagicMock()
    host.get_settings.return_value = {"mouse_rotation_sensitivity": sensitivity}
    style = CustomInteractorStyle(host)

    camera = MagicMock()
    renderer = MagicMock()
    renderer.GetActiveCamera.return_value = camera
    interactor = MagicMock()
    interactor.GetEventPosition.return_value = (130, 120)
    interactor.GetLastEventPosition.return_value = (100, 100)  # dx=30, dy=20
    interactor.GetLightFollowCamera.return_value = False

    style.GetCurrentRenderer = MagicMock(return_value=renderer)
    style.GetInteractor = MagicMock(return_value=interactor)
    style.GetMotionFactor = MagicMock(return_value=motion_factor)
    style.GetAutoAdjustCameraClippingRange = MagicMock(return_value=True)
    return style, camera


def test_rotation_speed_is_window_size_independent(app):
    """The applied azimuth/elevation depend only on mouse delta, motion factor,
    sensitivity and a fixed reference size — never on the live render size."""
    from moleditpy.ui.custom_interactor_style import _ROTATION_REFERENCE_SIZE

    style, camera = _rotation_style(sensitivity=1.0, motion_factor=10.0)
    style._rotate_size_independent()

    delta = -20.0 / _ROTATION_REFERENCE_SIZE * 10.0 * 1.0
    camera.Azimuth.assert_called_once_with(30 * delta)
    camera.Elevation.assert_called_once_with(20 * delta)
    camera.OrthogonalizeViewUp.assert_called_once()


def test_rotation_speed_scales_with_sensitivity(app):
    """Doubling the sensitivity setting doubles the rotation applied."""
    s1, cam1 = _rotation_style(sensitivity=1.0)
    s1._rotate_size_independent()
    s2, cam2 = _rotation_style(sensitivity=2.0)
    s2._rotate_size_independent()

    az1 = cam1.Azimuth.call_args[0][0]
    az2 = cam2.Azimuth.call_args[0][0]
    assert az2 == pytest.approx(az1 * 2.0)


def test_rotation_handles_missing_renderer(app):
    """No renderer yet (early startup) must not raise."""
    style, _ = _rotation_style()
    style.GetCurrentRenderer = MagicMock(return_value=None)
    style._rotate_size_independent()  # must not raise
