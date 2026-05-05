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
