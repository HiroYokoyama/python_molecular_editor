import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtWidgets import QMainWindow, QApplication
from moleditpy.ui.ui_manager import MainWindowUiManager

class MockUiManager(MainWindowUiManager, QMainWindow):
    def __init__(self):
        super().__init__()
        self.scene = MagicMock()
        self.view_2d = MagicMock()
        self.statusBar_mock = MagicMock()
        self.splitter = MagicMock()
        self.tool_group = MagicMock()
        self.mode_actions = {}
        self.settings = {}
        self.initial_settings = {}
        
    def statusBar(self):
        return self.statusBar_mock

def test_close_event_robustness(app):
    """Test that closeEvent handles missing attributes gracefully."""
    win = MockUiManager()
    event = MagicMock()
    
    # These should not raise even if attributes are missing or malformed
    win.closeEvent(event)
    assert event.accept.called

def test_enable_3d_features_robustness():
    """Test that _enable_3d_features handles missing widgets."""
    win = MockUiManager()
    # Delete some common attributes to test robustness
    if hasattr(win, "optimize_3d_button"): del win.optimize_3d_button
    
    # Should not raise
    win._enable_3d_features(True)

def test_set_mode_robustness():
    """Test that set_mode handles template preview visibility."""
    win = MockUiManager()
    win.scene.mode = "select"
    win.scene.template_preview = MagicMock()
    
    win.set_mode("atom_C")
    assert win.scene.mode == "atom_C"
    assert win.scene.template_preview.hide.called
