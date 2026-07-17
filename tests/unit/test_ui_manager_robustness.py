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

from unittest.mock import patch


# ---------------------------------------------------------------------------
# set_mode branches
# ---------------------------------------------------------------------------

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QGraphicsView


def test_set_mode_tuple_becomes_bond_mode():
    ui = _make_ui_manager()
    ui.set_mode((2, 1))
    assert ui.host.init_manager.scene.mode == "bond_2_1"
    assert ui.host.init_manager.scene.bond_order == 2
    assert ui.host.init_manager.scene.bond_stereo == 1


def test_set_mode_atom_sets_symbol_and_cross_cursor():
    ui = _make_ui_manager()
    ui.set_mode("atom_N")
    assert ui.host.init_manager.scene.current_atom_symbol == "N"
    ui.host.init_manager.view_2d.setCursor.assert_called_with(
        Qt.CursorShape.CrossCursor
    )
    ui.host.init_manager.view_2d.setDragMode.assert_called_with(
        QGraphicsView.DragMode.NoDrag
    )
    ui.host.update_status_message.assert_any_call("Mode: Draw Atom (N)")


def test_set_mode_bond_parses_order_and_stereo():
    ui = _make_ui_manager()
    ui.set_mode("bond_1_2")
    assert ui.host.init_manager.scene.bond_order == 1
    assert ui.host.init_manager.scene.bond_stereo == 2
    ui.host.update_status_message.assert_any_call(
        "Mode: Draw Bond (Order: 1 (Dash))"
    )


def test_set_mode_builtin_template_status():
    ui = _make_ui_manager()
    ui.set_mode("template_benzene")
    ui.host.statusBar().showMessage.assert_any_call("Mode: Benzene Template")


def test_set_mode_user_template_status():
    ui = _make_ui_manager()
    ui.set_mode("template_user_MyFrag")
    ui.host.statusBar().showMessage.assert_any_call(
        "Mode: User Template (MyFrag)"
    )


def test_set_mode_leaving_template_clears_preview():
    ui = _make_ui_manager()
    ui.host.init_manager.scene.mode = "template_benzene"
    ui.set_mode("select")
    ui.host.init_manager.scene.clear_template_preview.assert_called_once()


def test_set_mode_charge_modes_status():
    ui = _make_ui_manager()
    ui.set_mode("charge_plus")
    ui.host.statusBar().showMessage.assert_any_call(
        "Mode: Increase Charge (Click on Atom)"
    )
    ui.set_mode("charge_minus")
    ui.host.statusBar().showMessage.assert_any_call(
        "Mode: Decrease Charge (Click on Atom)"
    )


def test_set_mode_radical_status():
    ui = _make_ui_manager()
    ui.set_mode("radical")
    ui.host.update_status_message.assert_any_call(
        "Mode: Toggle Radical (Click on Atom)"
    )


def test_set_mode_select_uses_rubber_band():
    ui = _make_ui_manager()
    ui.set_mode("select")
    ui.host.init_manager.view_2d.setDragMode.assert_called_with(
        QGraphicsView.DragMode.RubberBandDrag
    )
    ui.host.update_status_message.assert_any_call("Mode: Select")


# ---------------------------------------------------------------------------
# Toolbar sync + shortcuts
# ---------------------------------------------------------------------------


def test_set_mode_and_update_toolbar_checks_matching_action():
    ui = _make_ui_manager()
    select_action = MagicMock()
    atom_action = MagicMock()
    ui.host.init_manager.mode_actions = {
        "select": select_action,
        "atom_C": atom_action,
    }
    btn = MagicMock()
    ui.host.init_manager.toolbar.widgetForAction.return_value = btn

    ui.set_mode_and_update_toolbar("atom_C")

    atom_action.setChecked.assert_called_with(True)
    select_action.setChecked.assert_called_with(False)
    btn.setStyleSheet.assert_called_with("")


def test_set_mode_and_update_toolbar_user_template_highlight():
    ui = _make_ui_manager()
    template_action = MagicMock()
    ui.host.init_manager.mode_actions = {"template_user": template_action}
    btn = MagicMock()
    ui.host.init_manager.toolbar.widgetForAction.return_value = btn

    ui.set_mode_and_update_toolbar("template_user_MyFrag")

    template_action.setChecked.assert_called_with(True)
    assert "background-color" in btn.setStyleSheet.call_args.args[0]


def test_activate_select_mode_checks_select_action():
    ui = _make_ui_manager()
    select_action = MagicMock()
    ui.host.init_manager.mode_actions = {"select": select_action}
    ui.activate_select_mode()
    assert ui.host.init_manager.scene.mode == "select"
    select_action.setChecked.assert_called_with(True)


def test_set_atom_from_periodic_table():
    ui = _make_ui_manager()
    ui.set_atom_from_periodic_table("Fe")
    assert ui.host.init_manager.scene.mode == "atom_Fe"
    assert ui.host.init_manager.scene.current_atom_symbol == "Fe"


# ---------------------------------------------------------------------------
# Drag & drop
# ---------------------------------------------------------------------------


def _drag_event_with_files(*paths, local=True):
    event = MagicMock()
    urls = []
    for p in paths:
        url = MagicMock()
        url.isLocalFile.return_value = local
        url.toLocalFile.return_value = p
        urls.append(url)
    event.mimeData.return_value.hasUrls.return_value = bool(urls)
    event.mimeData.return_value.urls.return_value = urls
    return event


def test_drag_enter_accepts_known_extension():
    ui = _make_ui_manager()
    event = _drag_event_with_files("C:/mols/foo.mol")
    ui.handle_drag_enter_event(event)
    event.acceptProposedAction.assert_called_once()


def test_drag_enter_ignores_unknown_extension_without_plugins():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    event = _drag_event_with_files("C:/mols/foo.docx")
    ui.handle_drag_enter_event(event)
    event.ignore.assert_called_once()
    event.acceptProposedAction.assert_not_called()


def test_drag_enter_accepts_unknown_extension_with_plugin_handler():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = [{"callback": MagicMock()}]
    event = _drag_event_with_files("C:/mols/foo.custom")
    ui.handle_drag_enter_event(event)
    event.acceptProposedAction.assert_called_once()


def test_drag_enter_ignores_null_mime():
    ui = _make_ui_manager()
    event = MagicMock()
    event.mimeData.return_value = None
    ui.handle_drag_enter_event(event)
    event.ignore.assert_called_once()


def test_drop_project_file_opens_project():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    event = _drag_event_with_files("C:/mols/session.pmeprj")
    with patch("moleditpy.ui.ui_manager.QTimer.singleShot"):
        ui.handle_drop_event(event)
    ui.host.io_manager.open_project_file.assert_called_once_with(
        file_path="C:/mols/session.pmeprj"
    )
    event.acceptProposedAction.assert_called_once()


def test_drop_mol_file_on_2d_loads_into_editor():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    ui.host.view_3d_manager.plotter = None  # no 3D target check possible
    event = _drag_event_with_files("C:/mols/mol.sdf")
    with patch("moleditpy.ui.ui_manager.QTimer.singleShot"):
        ui.handle_drop_event(event)
    ui.host.io_manager.load_mol_file.assert_called_once_with(
        file_path="C:/mols/mol.sdf"
    )


def test_drop_mol_file_on_3d_loads_for_viewing():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    plotter_widget = MagicMock()
    ui.host.init_manager.splitter.widget.return_value = plotter_widget
    ui.host.childAt.return_value = plotter_widget  # dropped onto the plotter
    event = _drag_event_with_files("C:/mols/mol.mol")
    with patch("moleditpy.ui.ui_manager.QTimer.singleShot"):
        ui.handle_drop_event(event)
    ui.host.io_manager.load_mol_file_for_3d_viewing.assert_called_once_with(
        file_path="C:/mols/mol.mol"
    )
    ui.host.io_manager.load_mol_file.assert_not_called()


def test_drop_xyz_loads_for_3d_viewing():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    event = _drag_event_with_files("C:/mols/geom.xyz")
    with patch("moleditpy.ui.ui_manager.QTimer.singleShot"):
        ui.handle_drop_event(event)
    ui.host.io_manager.load_xyz_for_3d_viewing.assert_called_once_with(
        file_path="C:/mols/geom.xyz"
    )


def test_drop_plugin_handler_wins_over_builtin():
    ui = _make_ui_manager()
    handler = MagicMock(return_value=True)
    ui.host.plugin_manager.drop_handlers = [{"callback": handler}]
    event = _drag_event_with_files("C:/mols/mol.mol")
    ui.handle_drop_event(event)
    handler.assert_called_once_with("C:/mols/mol.mol")
    ui.host.io_manager.load_mol_file.assert_not_called()
    event.acceptProposedAction.assert_called_once()


def test_drop_broken_plugin_handler_falls_through_to_builtin():
    ui = _make_ui_manager()
    handler = MagicMock(side_effect=RuntimeError("plugin bug"))
    ui.host.plugin_manager.drop_handlers = [{"callback": handler}]
    event = _drag_event_with_files("C:/mols/geom.xyz")
    with patch("moleditpy.ui.ui_manager.QTimer.singleShot"):
        ui.handle_drop_event(event)
    ui.host.io_manager.load_xyz_for_3d_viewing.assert_called_once()


def test_drop_unsupported_extension_reports_status():
    ui = _make_ui_manager()
    ui.host.plugin_manager.drop_handlers = []
    event = _drag_event_with_files("C:/mols/notes.txt")
    ui.handle_drop_event(event)
    event.ignore.assert_called_once()
    assert any(
        "Unsupported file type" in str(c.args[0])
        for c in ui.host.statusBar().showMessage.call_args_list
    )


def test_drop_without_local_file_ignored():
    ui = _make_ui_manager()
    event = _drag_event_with_files("http://x/foo.mol", local=False)
    ui.handle_drop_event(event)
    event.ignore.assert_called_once()


# ---------------------------------------------------------------------------
# Panel layout
# ---------------------------------------------------------------------------


def test_minimize_2d_panel_collapses_left():
    ui = _make_ui_manager()
    ui.host.init_manager.splitter.sizes.return_value = [400, 800]
    ui.minimize_2d_panel()
    ui.host.init_manager.splitter.setSizes.assert_called_once_with([0, 1200])


def test_minimize_2d_panel_noop_when_already_minimized():
    ui = _make_ui_manager()
    ui.host.init_manager.splitter.sizes.return_value = [0, 1200]
    ui.minimize_2d_panel()
    ui.host.init_manager.splitter.setSizes.assert_not_called()


def test_restore_2d_panel_only_when_hidden():
    ui = _make_ui_manager()
    ui.host.init_manager.splitter.sizes.return_value = [0, 1200]
    ui.restore_2d_panel()
    ui.host.init_manager.splitter.setSizes.assert_called_once_with([600, 600])

    ui.host.init_manager.splitter.setSizes.reset_mock()
    ui.host.init_manager.splitter.sizes.return_value = [500, 700]
    ui.restore_2d_panel()
    ui.host.init_manager.splitter.setSizes.assert_not_called()


def test_set_panel_layout_splits_width_by_percent():
    ui = _make_ui_manager()
    ui.host.init_manager.splitter.width.return_value = 1000
    ui.set_panel_layout(30, 70)
    ui.host.init_manager.splitter.setSizes.assert_called_once_with([300, 700])
    assert any(
        "Panel layout set to 30% : 70%" in str(c.args[0])
        for c in ui.host.statusBar().showMessage.call_args_list
    )


def test_set_panel_layout_rejects_bad_percentages():
    ui = _make_ui_manager()
    ui.set_panel_layout(30, 60)  # 90 != 100
    ui.host.init_manager.splitter.setSizes.assert_not_called()


def test_set_panel_layout_uses_default_width_when_zero():
    ui = _make_ui_manager()
    ui.host.init_manager.splitter.width.return_value = 0
    ui.set_panel_layout(50, 50)
    ui.host.init_manager.splitter.setSizes.assert_called_once_with([600, 600])


def test_toggle_2d_panel_restores_when_hidden():
    ui = _make_ui_manager()
    ui.host.ui_manager = MagicMock()
    ui.host.init_manager.splitter.sizes.return_value = [0, 1200]
    ui.toggle_2d_panel()
    ui.host.ui_manager.restore_2d_panel.assert_called_once()


def test_toggle_2d_panel_minimizes_when_visible():
    ui = _make_ui_manager()
    ui.host.ui_manager = MagicMock()
    ui.host.init_manager.splitter.sizes.return_value = [600, 600]
    ui.toggle_2d_panel()
    ui.host.ui_manager.minimize_2d_panel.assert_called_once()


def test_toggle_2d_panel_empty_sizes_noop():
    ui = _make_ui_manager()
    ui.host.ui_manager = MagicMock()
    ui.host.init_manager.splitter.sizes.return_value = []
    ui.toggle_2d_panel()
    ui.host.ui_manager.restore_2d_panel.assert_not_called()
    ui.host.ui_manager.minimize_2d_panel.assert_not_called()


# ---------------------------------------------------------------------------
# 3D viewer UI mode transitions
# ---------------------------------------------------------------------------


def test_enter_3d_viewer_ui_mode_disables_2d_tools():
    ui = _make_ui_manager()
    action = MagicMock()
    ui.host.init_manager.tool_group.actions.return_value = [action]
    ui.host.init_manager.splitter.sizes.return_value = [600, 600]
    ui.enable_3d_features = MagicMock()

    ui.enter_3d_viewer_ui_mode()

    assert ui.is_2d_editable is False
    ui.host.init_manager.cleanup_button.setEnabled.assert_called_with(False)
    ui.host.init_manager.convert_button.setEnabled.assert_called_with(False)
    action.setEnabled.assert_called_with(False)
    ui.enable_3d_features.assert_called_once_with(True)


def test_restore_ui_for_editing_reenables_2d_tools():
    ui = _make_ui_manager()
    action = MagicMock()
    ui.host.init_manager.tool_group.actions.return_value = [action]
    ui.host.init_manager.splitter.sizes.return_value = [0, 1200]
    ui._enable_3d_edit_actions = MagicMock()

    ui.restore_ui_for_editing()

    assert ui.is_2d_editable is True
    ui.host.init_manager.cleanup_button.setEnabled.assert_called_with(True)
    action.setEnabled.assert_called_with(True)
    ui._enable_3d_edit_actions.assert_called_once_with(False)


def test_enable_3d_edit_actions_toggles_all_actions_and_menus():
    ui = _make_ui_manager()
    ui._enable_3d_edit_actions(True)
    ui.host.translation_action.setEnabled.assert_called_with(True)
    ui.host.align_menu.setEnabled.assert_called_with(True)
    ui._enable_3d_edit_actions(False)
    ui.host.dihedral_action.setEnabled.assert_called_with(False)
    ui.host.align_menu.setEnabled.assert_called_with(False)
