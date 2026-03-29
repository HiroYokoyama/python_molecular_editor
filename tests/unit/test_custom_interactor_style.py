import pytest
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
    mock_parser_host.current_mol = mock_mol
    
    # Mock VTK Picker
    mock_picker = MagicMock()
    mock_parser_host.plotter.picker = mock_picker
    mock_parser_host.atom_actor = MagicMock()
    mock_picker.GetActor.return_value = mock_parser_host.atom_actor
    mock_picker.GetPickPosition.return_value = (0, 0, 0)
    
    import numpy as np
    mock_parser_host.view_3d_manager.atom_positions_3d = np.array([[0.0, 0.0, 0.0]])
    
    with patch('moleditpy.ui.custom_interactor_style.QApplication') as mock_qapp, \
         patch('vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera.OnLeftButtonDown'):
         
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier
        
        interactor_style.on_left_button_down(None, None)
        
        # Verify the atom was selected
        assert interactor_style._is_dragging_atom is True
        mock_parser_host.plotter.setCursor.assert_called_once_with(Qt.CursorShape.ClosedHandCursor)
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
    
    with patch('moleditpy.ui.custom_interactor_style.QApplication') as mock_qapp, \
         patch('vtkmodules.vtkInteractionStyle.vtkInteractorStyleTrackballCamera.OnLeftButtonDown') as mock_super_down:
         
        mock_qapp.topLevelWidgets.return_value = []
        mock_qapp.keyboardModifiers.return_value = Qt.KeyboardModifier.NoModifier
        
        interactor_style.on_left_button_down(None, None)
        
        # Verify drag is false and superclass rotation method was called
        assert interactor_style._is_dragging_atom is False
        mock_super_down.assert_called_once()
