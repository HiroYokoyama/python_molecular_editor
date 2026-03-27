# -*- coding: utf-8 -*-
import pytest
import os
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QColor, QBrush
from PyQt6.QtWidgets import QMessageBox, QApplication
from moleditpy.modules.constants import DEFAULT_CPK_COLORS

def test_load_command_line_file_with_plugin(window, monkeypatch):
    """Test that load_command_line_file uses plugin openers when available."""
    mock_callback = MagicMock()
    window.plugin_manager.file_openers = {".testext": [{"callback": mock_callback, "plugin": "TestPlugin"}]}
    
    test_file = "test.testext"
    # Mock os.path.exists to return True for our dummy file
    monkeypatch.setattr(os.path, "exists", lambda x: True if x == test_file else os.path.exists(x))
    
    window.load_command_line_file(test_file)
    
    mock_callback.assert_called_once_with(test_file)
    assert window.current_file_path == test_file

def test_load_command_line_file_default_extensions(window, monkeypatch):
    """Test that load_command_line_file handles standard extensions."""
    # Mock the specific loading methods called by load_command_line_file
    window.load_mol_file_for_3d_viewing = MagicMock()
    window.load_xyz_for_3d_viewing = MagicMock()
    window.open_project_file = MagicMock()
    
    # Mock os.path.exists
    monkeypatch.setattr(os.path, "exists", lambda x: True)
    
    window.load_command_line_file("test.mol")
    window.load_mol_file_for_3d_viewing.assert_called_once_with("test.mol")
    
    window.load_command_line_file("test.xyz")
    window.load_xyz_for_3d_viewing.assert_called_once_with("test.xyz")
    
    window.load_command_line_file("test.pmeprj")
    window.open_project_file.assert_called_once_with(file_path="test.pmeprj")

def test_update_cpk_colors_from_settings(window):
    """Test that CPK colors are updated correctly from settings overrides."""
    from moleditpy.modules import constants
    
    # Define an override
    overrides = {"C": "#FF0000"}  # Carbon as Red
    window.settings["cpk_colors"] = overrides
    
    window.update_cpk_colors_from_settings()
    
    # Verify constants.CPK_COLORS["C"] is now Red
    assert constants.CPK_COLORS["C"] == QColor("#FF0000")
    # Verify CPK_COLORS_PV for "C" is [1.0, 0.0, 0.0]
    assert constants.CPK_COLORS_PV["C"] == [1.0, 0.0, 0.0]
    
    # Reset to defaults for other tests (important since constants are global)
    window.settings["cpk_colors"] = {}
    window.update_cpk_colors_from_settings()

def test_apply_initial_settings(window, monkeypatch):
    """Test that apply_initial_settings updates scene background and style."""
    window.settings["background_color_2d"] = "#FFEECC"
    window.settings["background_color"] = "#112233"
    
    # Mock plotter if it exists
    if window.plotter:
        window.plotter.set_background = MagicMock()
        window.apply_3d_settings = MagicMock()
    
    window.apply_initial_settings()
    
    # Check 2D scene background
    expected_color = QColor("#FFEECC")
    assert window.scene.backgroundBrush().color() == expected_color
    
    # Check 3D background if plotter was mocked
    if window.plotter:
        window.plotter.set_background.assert_called_with("#112233")

def test_reset_all_settings_flow(window, monkeypatch):
    """Test the complete settings reset flow."""
    # Mock confirmation to return Yes
    monkeypatch.setattr(QMessageBox, "question", lambda *args: QMessageBox.StandardButton.Yes)
    monkeypatch.setattr(QMessageBox, "information", lambda *args: None)
    
    # Mock internal reset methods
    window._perform_settings_reset = MagicMock()
    window._refresh_ui_after_reset = MagicMock()
    
    window.reset_all_settings_menu()
    
    window._perform_settings_reset.assert_called_once()
    window._refresh_ui_after_reset.assert_called_once()

def test_perform_settings_reset_logic(window, monkeypatch, tmp_path):
    """Test the low-level settings reset (file deletion and reload)."""
    # Create a dummy settings file
    settings_dir = tmp_path / ".moleditpy"
    settings_dir.mkdir()
    settings_file = settings_dir / "settings.json"
    settings_file.write_text("{}")
    
    window.settings_dir = str(settings_dir)
    window.settings_file = str(settings_file)
    window.load_settings = MagicMock()
    
    window._perform_settings_reset()
    
    assert not settings_file.exists()
    window.load_settings.assert_called_once()
    assert window.settings_dirty is True
