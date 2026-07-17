"""Unit tests for UIManager robustness and edge cases."""

from unittest.mock import MagicMock
from moleditpy.ui.ui_manager import UIManager


def _make_ui_manager():
    """Create a UIManager with a mock host in the new manager pattern."""
    host = MagicMock()
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
    host.state_manager.has_unsaved_changes = False
    host.view_3d_manager = MagicMock()
    host.io_manager = MagicMock()

    def set_scene_mode(mode):
        host.init_manager.scene.mode = mode

    host.set_scene_mode.side_effect = set_scene_mode

    def set_scene_atom_symbol(symbol):
        host.init_manager.scene.current_atom_symbol = symbol

    host.set_scene_atom_symbol.side_effect = set_scene_atom_symbol

    def set_scene_bond_properties(order, stereo=0):
        host.init_manager.scene.bond_order = order
        host.init_manager.scene.bond_stereo = stereo

    host.set_scene_bond_properties.side_effect = set_scene_bond_properties

    def update_status_message(msg):
        pass

    host.update_status_message.side_effect = update_status_message

    def set_settings_dirty(val):
        host.init_manager.settings_dirty = val

    host.set_settings_dirty.side_effect = set_settings_dirty

    return UIManager(host)


def test_close_event_robustness():
    """Test that handle_close_event handles missing attributes gracefully."""
    ui = _make_ui_manager()
    event = MagicMock()
    # Should not raise even if attributes are missing or malformed
    result = ui.handle_close_event(event)
    # No unsaved changes in mock setup — function should allow close
    assert result is True


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


# ---------------------------------------------------------------------------
# Interactor style installation / watchdog (fast-click freeze root cause)
# ---------------------------------------------------------------------------


class _NoIrenPlotter:
    """Plotter stand-in without pyvista's iren attribute."""

    def __init__(self):
        self.interactor = MagicMock()

    @property
    def iren(self):
        raise AttributeError("no iren")


def test_install_interactor_style_registers_via_pyvista():
    """Style must go through iren.style so pyvista's update_style() keeps it."""
    ui = _make_ui_manager()
    style = object()

    ui._install_interactor_style(style)

    plotter = ui.host.view_3d_manager.plotter
    assert plotter.iren.style is style
    plotter.interactor.SetInteractorStyle.assert_not_called()


def test_install_interactor_style_falls_back_to_raw_vtk():
    ui = _make_ui_manager()
    plotter = _NoIrenPlotter()
    ui.host.view_3d_manager.plotter = plotter
    style = object()

    ui._install_interactor_style(style)

    plotter.interactor.SetInteractorStyle.assert_called_once_with(style)


def test_style_watchdog_reinstalls_replaced_style():
    """A style evicted by pyvista (e.g. its double-click chart handler) is reinstalled."""
    ui = _make_ui_manager()
    expected = object()
    ui._expected_style = expected
    plotter = ui.host.view_3d_manager.plotter
    plotter.interactor.GetInteractorStyle.return_value = MagicMock()

    ui._check_interactor_style()

    assert plotter.iren.style is expected


def test_style_watchdog_keeps_installed_style():
    ui = _make_ui_manager()
    expected = object()
    ui._expected_style = expected
    plotter = ui.host.view_3d_manager.plotter
    plotter.interactor.GetInteractorStyle.return_value = expected

    ui._check_interactor_style()

    assert plotter.iren.style is not expected  # no reinstall happened


def test_style_watchdog_ignores_rubberband_style():
    """pyvista's temporary box-selection style must not be fought by the watchdog."""
    ui = _make_ui_manager()
    ui._expected_style = object()
    rubber = type("InteractorStyleRubberBandPick", (), {})()
    plotter = ui.host.view_3d_manager.plotter
    plotter.interactor.GetInteractorStyle.return_value = rubber

    ui._check_interactor_style()

    assert plotter.iren.style is not ui._expected_style
