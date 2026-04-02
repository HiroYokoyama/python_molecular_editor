import pytest
from PyQt6.QtCore import Qt
from unittest.mock import patch, MagicMock
from moleditpy.ui.about_dialog import AboutDialog

def test_about_dialog_initialization(app, mock_parser_host):
    """Verify AboutDialog initializes correctly with the injected host."""
    dialog = AboutDialog(main_window=mock_parser_host)
    assert dialog.windowTitle() == "About MoleditPy"
    assert dialog.image_label is not None

def test_about_dialog_easter_egg(app, mock_parser_host):
    """Verify the easter egg triggers correctly on right click."""
    dialog = AboutDialog(main_window=mock_parser_host)
    
    # Mock a right-click mouse event
    event = MagicMock()
    event.button.return_value = Qt.MouseButton.RightButton
    
    dialog.image_mouse_press_event(event)
    
    # Assert that it cleared the scene and loaded the specific SMILES
    mock_parser_host.edit_actions_manager.clear_all.assert_called_once()
    mock_parser_host.string_importer_manager.load_from_smiles.assert_called_once_with("C1=CN=C(N=C1)C2=NC=CC=N2")

def test_about_dialog_ignore_left_click(app, mock_parser_host):
    """Verify left clicks do not trigger the easter egg."""
    dialog = AboutDialog(main_window=mock_parser_host)
    
    event = MagicMock()
    event.button.return_value = Qt.MouseButton.LeftButton
    
    dialog.image_mouse_press_event(event)
    
    mock_parser_host.clear_all.assert_not_called()
    event.ignore.assert_called_once()
