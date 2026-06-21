"""Unit tests for CustomInteractorStyle 3D event handling."""
from unittest.mock import MagicMock, patch
from moleditpy.ui.custom_interactor_style import CustomInteractorStyle
from PyQt6.QtCore import Qt


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
