# -*- coding: utf-8 -*-
import os
from unittest.mock import MagicMock
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QMessageBox
from moleditpy.utils.constants import DEFAULT_CPK_COLORS

def test_load_command_line_file_with_plugin(window, monkeypatch):
    """Test that load_command_line_file uses plugin openers when available."""
    mock_callback = MagicMock()
    window.plugin_manager.file_openers = {".testext": [{"callback": mock_callback, "plugin": "TestPlugin"}]}
    
    test_file = "test.testext"
    # Mock os.path.exists to return True for our dummy file
    monkeypatch.setattr(os.path, "exists", lambda x: True if x == test_file else os.path.exists(x))
    
    window.init_manager.load_command_line_file(test_file)
    
    mock_callback.assert_called_once_with(test_file)
    assert window.init_manager.current_file_path == test_file

def test_load_command_line_file_default_extensions(window, monkeypatch):
    """Test that load_command_line_file handles standard extensions."""
    # Mock the specific loading methods called by load_command_line_file
    window.io_manager.load_mol_file_for_3d_viewing = MagicMock()
    window.io_manager.load_xyz_for_3d_viewing = MagicMock()
    window.io_manager.open_project_file = MagicMock()
    
    # Mock os.path.exists
    monkeypatch.setattr(os.path, "exists", lambda x: True)
    
    window.init_manager.load_command_line_file("test.mol")
    window.io_manager.load_mol_file_for_3d_viewing.assert_called_once_with("test.mol")
    
    window.init_manager.load_command_line_file("test.xyz")
    window.io_manager.load_xyz_for_3d_viewing.assert_called_once_with("test.xyz")
    
    window.init_manager.load_command_line_file("test.pmeprj")
    window.io_manager.open_project_file.assert_called_once_with(file_path="test.pmeprj")

def test_update_cpk_colors_from_settings(window):
    """Test that CPK colors are updated correctly from settings overrides."""
    from moleditpy.utils import constants
    
    # Define an override
    overrides = {"C": "#FF0000"}  # Carbon as Red
    window.init_manager.settings["cpk_colors"] = overrides
    
    window.init_manager.update_cpk_colors_from_settings()
    
    # Verify constants.CPK_COLORS["C"] is now Red
    assert constants.CPK_COLORS["C"] == QColor("#FF0000")
    # Verify CPK_COLORS_PV for "C" is [1.0, 0.0, 0.0]
    assert constants.CPK_COLORS_PV["C"] == [1.0, 0.0, 0.0]
    
    # Reset to defaults for other tests (important since constants are global)
    window.init_manager.settings["cpk_colors"] = {}
    window.init_manager.update_cpk_colors_from_settings()
    assert constants.CPK_COLORS["C"] == DEFAULT_CPK_COLORS.get("C", constants.CPK_COLORS["C"])

def test_apply_initial_settings(window, monkeypatch):
    """Test that apply_initial_settings updates scene background and style."""
    window.init_manager.settings["background_color_2d"] = "#FFEECC"
    window.init_manager.settings["background_color"] = "#112233"
    
    window.view_3d_manager.plotter.set_background = MagicMock()
    window.view_3d_manager.apply_3d_settings = MagicMock()

    window.init_manager.apply_initial_settings()

    # Check 2D scene background
    expected_color = QColor("#FFEECC")
    assert window.init_manager.scene.backgroundBrush().color() == expected_color

    # Check 3D background
    window.view_3d_manager.plotter.set_background.assert_called_with("#112233")

def test_reset_all_settings_flow(window, monkeypatch):
    """Test the complete settings reset flow."""
    # Mock confirmation to return Yes
    monkeypatch.setattr(QMessageBox, "question", lambda *args: QMessageBox.StandardButton.Yes)
    monkeypatch.setattr(QMessageBox, "information", lambda *args: None)
    
    # Mock internal reset methods
    window.init_manager._perform_settings_reset = MagicMock()
    window.init_manager._refresh_ui_after_reset = MagicMock()
    
    window.init_manager.reset_all_settings_menu()
    
    window.init_manager._perform_settings_reset.assert_called_once()
    window.init_manager._refresh_ui_after_reset.assert_called_once()

def test_perform_settings_reset_logic(window, monkeypatch, tmp_path):
    """Test the low-level settings reset (file deletion and reload)."""
    # Create a dummy settings file
    settings_dir = tmp_path / ".moleditpy"
    settings_dir.mkdir()
    settings_file = settings_dir / "settings.json"
    settings_file.write_text("{}")
    
    window.init_manager.settings_dir = str(settings_dir)
    window.init_manager.settings_file = str(settings_file)
    window.init_manager.load_settings = MagicMock()
    
    window.init_manager._perform_settings_reset()
    
    assert not settings_file.exists()
    window.init_manager.load_settings.assert_called_once()
    assert window.init_manager.settings_dirty is True
