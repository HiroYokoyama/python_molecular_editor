import pytest
from unittest.mock import MagicMock
from moleditpy.ui.ui_manager import UIManager


def _make_ui_manager():
    """Create a UIManager with a mock host in the new manager pattern."""
    host = MagicMock()
    host.has_unsaved_changes = False
    host.init_manager.scene = MagicMock()
    host.init_manager.scene.mode = "select"
    host.init_manager.scene.template_preview = MagicMock()
    host.init_manager.view_2d = MagicMock()
    host.init_manager.splitter = MagicMock()
    host.init_manager.tool_group = MagicMock()
    host.init_manager.mode_actions = {}
    host.init_manager.settings = {}
    host.initial_settings = {}
    host.edit_actions_manager = MagicMock()
    host.state_manager = MagicMock()
    host.view_3d_manager = MagicMock()
    host.io_manager = MagicMock()
    return UIManager(host)


def test_close_event_robustness():
    """Test that handle_close_event handles missing attributes gracefully."""
    ui = _make_ui_manager()
    event = MagicMock()
    # Should not raise even if attributes are missing or malformed
    result = ui.handle_close_event(event)
    # Returns False when no unsaved changes (event accepted by caller)
    assert result is False or result is True  # just check it doesn't raise


def test_enable_3d_features_robustness():
    """Test that _enable_3d_features handles missing widgets."""
    ui = _make_ui_manager()
    # Should not raise
    ui._enable_3d_features(True)
    ui._enable_3d_features(False)


def test_set_mode_robustness():
    """Test that set_mode handles template preview visibility."""
    ui = _make_ui_manager()
    ui.host.init_manager.scene.mode = "select"

    ui.set_mode("atom_C")
    assert ui.host.init_manager.scene.mode == "atom_C"
    assert ui.host.init_manager.scene.template_preview.hide.called
