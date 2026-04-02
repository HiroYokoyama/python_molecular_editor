import pytest
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QColorDialog
from unittest.mock import patch, MagicMock
from moleditpy.ui.color_settings_dialog import ColorSettingsDialog

def test_color_settings_dialog_initialization(app, mock_parser_host):
    """Verify ColorSettingsDialog initializes and parses settings correctly."""
    settings = {"cpk_colors": {"C": "#123456"}, "ball_stick_bond_color": "#aabbcc"}
    dialog = ColorSettingsDialog(current_settings=settings, parent=None)
    mock_parser_host.init_manager.settings.clear()
    mock_parser_host.default_settings = {}
    dialog.parent_window = mock_parser_host
    
    assert dialog.windowTitle() == "CPK Colors"
    assert "C" in dialog.element_buttons
    # Initial color parsing check usually relies on btn styles, we just ensure it didn't crash
    assert dialog.changed_cpk == {}

@patch('moleditpy.ui.color_settings_dialog.QColorDialog.getColor')
def test_color_settings_dialog_pick_color(mock_get_color, app, mock_parser_host):
    """Verify that clicking an element button updates the changed_cpk dictionary."""
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_parser_host.init_manager.settings.clear()
    mock_parser_host.default_settings = {}
    dialog.parent_window = mock_parser_host
    
    # Mock valid color return from QColorDialog
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#ff0000"
    mock_color.red.return_value = 255
    mock_color.green.return_value = 0
    mock_color.blue.return_value = 0
    mock_get_color.return_value = mock_color
    
    btn = dialog.element_buttons["C"]
    
    with patch.object(dialog, 'sender', return_value=btn):
        dialog.on_element_clicked()
        
    assert dialog.changed_cpk["C"] == "#ff0000"
    assert "background-color: #ff0000" in btn.styleSheet()

def test_color_settings_dialog_reset_all(app, mock_parser_host):
    """Verify reset_all clears changes and sets the reset flag."""
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_parser_host.init_manager.settings.clear()
    mock_parser_host.default_settings = {}
    dialog.parent_window = mock_parser_host
    dialog.changed_cpk["C"] = "#ff0000"
    
    dialog.reset_all()
    
    assert dialog._reset_all_flag is True
    assert dialog.changed_cpk == {}

def test_color_settings_dialog_apply_changes(app, mock_parser_host):
    """Verify that apply_changes pushes updates to the parent window settings."""
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_parser_host.init_manager.settings.clear()
    mock_parser_host.default_settings = {}
    dialog.parent_window = mock_parser_host
    dialog.changed_cpk["O"] = "#00ff00"
    dialog.changed_bs_color = "#112233"
    
    # Ensure parent settings is a mock dictionary
    mock_parser_host.init_manager.settings.clear()
    mock_parser_host.default_settings = {}
    
    dialog.apply_changes()
    
    assert mock_parser_host.settings_dirty is True
    assert mock_parser_host.settings["cpk_colors"]["O"] == "#00ff00"
    assert mock_parser_host.settings["ball_stick_bond_color"] == "#112233"
