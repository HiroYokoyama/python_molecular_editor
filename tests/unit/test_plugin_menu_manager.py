"""Unit tests for PluginMenuManager — plugin UI lifecycle management."""
import pytest
from unittest.mock import MagicMock, patch, call
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QMenu

from moleditpy.ui.plugin_menu_manager import PluginMenuManager


# ---------------------------------------------------------------------------
# Module-level patch: QAction(text, MagicMock) fails because PyQt6 strictly
# validates the parent type. Strip the mock parent so tests stay lightweight.
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _patch_qaction(monkeypatch):
    """Allow QAction to be constructed with a MagicMock host in unit tests."""
    from PyQt6.QtGui import QAction as _RealQAction

    def _safe_qaction(text: str, parent=None) -> QAction:
        return _RealQAction(text)

    monkeypatch.setattr("moleditpy.ui.plugin_menu_manager.QAction", _safe_qaction)


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def make_init_manager(
    *,
    has_toolbar: bool = True,
    has_export_button: bool = True,
    has_style_button: bool = True,
    has_import_menu: bool = True,
    plugin_manager: MagicMock | None = None,
) -> MagicMock:
    """Return a minimal MainInitManager mock suitable for PluginMenuManager."""
    im = MagicMock()
    im.host.menuBar.return_value = MagicMock()
    im.host.menuBar.return_value.actions.return_value = []
    im.plugin_menubar_separator_added = False

    if plugin_manager is None:
        plugin_manager = MagicMock()
        plugin_manager.menu_actions = []
        plugin_manager.toolbar_actions = []
        plugin_manager.export_actions = []
        plugin_manager.file_openers = {}
        plugin_manager.analysis_tools = []
        plugin_manager.custom_3d_styles = []
    im.host.plugin_manager = plugin_manager

    if has_toolbar:
        im.plugin_toolbar = MagicMock()
    else:
        del im.plugin_toolbar

    if has_export_button:
        im.export_button = MagicMock()
        im.export_button.menu.return_value = MagicMock()
    else:
        im.export_button = None

    if has_style_button:
        im.style_button = MagicMock()
        im.style_button.menu.return_value = MagicMock()
    else:
        im.style_button = None

    if has_import_menu:
        im.import_menu = MagicMock()
    else:
        im.import_menu = None

    return im


@pytest.fixture
def im():
    return make_init_manager()


@pytest.fixture
def pmm(im) -> PluginMenuManager:
    return PluginMenuManager(im)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestConstruction:
    def test_holds_init_manager_reference(self, im):
        mgr = PluginMenuManager(im)
        assert mgr._im is im


# ---------------------------------------------------------------------------
# update_plugin_menu
# ---------------------------------------------------------------------------

class TestUpdatePluginMenu:
    def test_does_nothing_when_no_plugin_manager(self, im, pmm):
        im.host.plugin_manager = None
        menu = MagicMock(spec=QMenu)
        pmm.update_plugin_menu(menu)
        menu.clear.assert_not_called()

    def test_clears_and_adds_manage_action(self, im, pmm):
        menu = MagicMock(spec=QMenu)
        pmm.update_plugin_menu(menu)
        menu.clear.assert_called_once()
        # addAction called at least once (the "Plugin Manager..." entry)
        assert menu.addAction.called or menu.addSeparator.called

    def test_discover_plugins_called(self, im, pmm):
        menu = MagicMock(spec=QMenu)
        im.host.plugin_manager.discover_plugins.return_value = []
        pmm.update_plugin_menu(menu)
        im.host.plugin_manager.discover_plugins.assert_called_once_with(im.host)

    def test_all_integration_hooks_called(self, im):
        """update_plugin_menu must call all 6 integration methods."""
        pmm = PluginMenuManager(im)
        im.host.plugin_manager.discover_plugins.return_value = []

        methods = [
            "update_style_menu_with_plugins",
            "add_registered_plugin_actions",
            "add_plugin_toolbar_actions",
            "_add_legacy_plugin_actions",
            "integrate_plugin_export_actions",
            "integrate_plugin_file_openers",
            "integrate_plugin_analysis_tools",
        ]
        mocks = {}
        for m in methods:
            mocks[m] = MagicMock()
            setattr(pmm, m, mocks[m])

        pmm.update_plugin_menu(MagicMock(spec=QMenu))

        for m in methods:
            mocks[m].assert_called_once()


# ---------------------------------------------------------------------------
# rebuild_plugin_menus
# ---------------------------------------------------------------------------

class TestRebuildPluginMenus:
    def test_calls_all_six_rebuild_steps(self, im):
        pmm = PluginMenuManager(im)
        im.host.menuBar.return_value.actions.return_value = []

        steps = [
            "add_registered_plugin_actions",
            "add_plugin_toolbar_actions",
            "integrate_plugin_export_actions",
            "integrate_plugin_file_openers",
            "integrate_plugin_analysis_tools",
            "update_style_menu_with_plugins",
        ]
        mocks = {s: MagicMock() for s in steps}
        for s, m in mocks.items():
            setattr(pmm, s, m)

        pmm.rebuild_plugin_menus()

        for s, m in mocks.items():
            m.assert_called_once(), f"{s} was not called"

    def test_resets_separator_flag(self, im, pmm):
        im.plugin_menubar_separator_added = True
        im.host.menuBar.return_value.actions.return_value = []
        pmm.rebuild_plugin_menus()
        assert im.plugin_menubar_separator_added is False

    def test_one_step_exception_does_not_stop_others(self, im):
        """A failing step must not prevent subsequent steps from running."""
        pmm = PluginMenuManager(im)
        im.host.menuBar.return_value.actions.return_value = []

        call_order = []

        def boom():
            call_order.append("boom")
            raise RuntimeError("step failed")

        def ok_step(name):
            def _step():
                call_order.append(name)
            return _step

        pmm.add_registered_plugin_actions = boom
        pmm.add_plugin_toolbar_actions = ok_step("toolbar")
        pmm.integrate_plugin_export_actions = ok_step("export")
        pmm.integrate_plugin_file_openers = ok_step("file_openers")
        pmm.integrate_plugin_analysis_tools = ok_step("analysis")
        pmm.update_style_menu_with_plugins = ok_step("style")

        pmm.rebuild_plugin_menus()

        assert "boom" in call_order
        assert "toolbar" in call_order
        assert "export" in call_order
        assert "style" in call_order

    def test_menu_cleanup_removes_tagged_actions(self, im, pmm):
        """Tagged plugin actions are stripped before rebuild."""
        tagged_action = MagicMock()
        tagged_action.data.return_value = "plugin_managed"
        tagged_action.menu.return_value = None

        top_menu = MagicMock()
        top_menu.actions.return_value = [tagged_action]

        top_action = MagicMock()
        top_action.menu.return_value = top_menu

        im.host.menuBar.return_value.actions.return_value = [top_action]

        pmm.rebuild_plugin_menus()

        top_menu.removeAction.assert_called_with(tagged_action)


# ---------------------------------------------------------------------------
# add_registered_plugin_actions
# ---------------------------------------------------------------------------

class TestAddRegisteredPluginActions:
    def test_no_menu_actions_does_nothing(self, im, pmm):
        im.host.plugin_manager.menu_actions = []
        pmm.add_registered_plugin_actions()
        im.host.menuBar.return_value.addMenu.assert_not_called()

    def test_creates_new_top_level_menu_with_separator(self, im, pmm):
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {"path": "MyPlugin/Action", "callback": callback, "text": "Run It", "shortcut": None}
        ]
        im.host.menuBar.return_value.actions.return_value = []  # no existing menus
        new_menu = MagicMock()
        new_menu.actions.return_value = []
        im.host.menuBar.return_value.addMenu.return_value = new_menu

        pmm.add_registered_plugin_actions()

        # A separator should have been added to menubar before the new menu
        im.host.menuBar.return_value.addSeparator.assert_called_once()
        im.host.menuBar.return_value.addMenu.assert_called_once_with("MyPlugin")
        new_menu.addAction.assert_called_once()

    def test_reuses_existing_top_level_menu(self, im, pmm):
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {"path": "File/ExportXYZ", "callback": callback, "text": "Export XYZ", "shortcut": None}
        ]
        existing_menu = MagicMock()
        existing_menu.actions.return_value = []
        file_action = MagicMock()
        file_action.menu.return_value = existing_menu
        file_action.text.return_value = "File"

        im.host.menuBar.return_value.actions.return_value = [file_action]

        pmm.add_registered_plugin_actions()

        # Should NOT create a new top-level menu
        im.host.menuBar.return_value.addMenu.assert_not_called()
        existing_menu.addAction.assert_called_once()

    def test_separator_added_only_once_for_multiple_plugins(self, im, pmm):
        """The menubar separator is added exactly once even with multiple new menus."""
        im.plugin_menubar_separator_added = False
        im.host.plugin_manager.menu_actions = [
            {"path": "PlugA/Action1", "callback": MagicMock(), "text": "A1", "shortcut": None},
            {"path": "PlugB/Action2", "callback": MagicMock(), "text": "B1", "shortcut": None},
        ]
        im.host.menuBar.return_value.actions.return_value = []
        im.host.menuBar.return_value.addMenu.return_value = MagicMock(
            **{"actions.return_value": []}
        )

        pmm.add_registered_plugin_actions()

        assert im.host.menuBar.return_value.addSeparator.call_count == 1

    def test_shortcut_applied_when_present(self, im, pmm):
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {"path": "Plug/Act", "callback": callback, "text": "Act", "shortcut": "Ctrl+P"}
        ]
        im.host.menuBar.return_value.actions.return_value = []
        new_menu = MagicMock(**{"actions.return_value": []})
        im.host.menuBar.return_value.addMenu.return_value = new_menu

        pmm.add_registered_plugin_actions()

        added_action = new_menu.addAction.call_args[0][0]
        assert isinstance(added_action, QAction)
        # shortcut is set — just verify no exception was raised


# ---------------------------------------------------------------------------
# add_plugin_toolbar_actions
# ---------------------------------------------------------------------------

class TestAddPluginToolbarActions:
    def test_no_toolbar_attribute_does_nothing(self):
        im = make_init_manager(has_toolbar=False)
        pmm = PluginMenuManager(im)
        pmm.add_plugin_toolbar_actions()  # must not raise

    def test_hides_toolbar_when_no_actions(self, im, pmm):
        im.host.plugin_manager.toolbar_actions = []
        pmm.add_plugin_toolbar_actions()
        im.plugin_toolbar.hide.assert_called_once()
        im.plugin_toolbar.clear.assert_called_once()

    def test_shows_toolbar_and_adds_actions(self, im, pmm):
        im.host.plugin_manager.toolbar_actions = [
            {"text": "Run", "callback": MagicMock(), "icon": "", "tooltip": "Run it"},
        ]
        pmm.add_plugin_toolbar_actions()
        im.plugin_toolbar.show.assert_called_once()
        im.plugin_toolbar.addAction.assert_called_once()

    def test_icon_set_when_file_exists(self, im, pmm, tmp_path):
        icon_file = tmp_path / "icon.png"
        icon_file.write_bytes(b"")  # create the file
        im.host.plugin_manager.toolbar_actions = [
            {"text": "Icon Action", "callback": MagicMock(), "icon": str(icon_file), "tooltip": ""}
        ]
        pmm.add_plugin_toolbar_actions()
        added = im.plugin_toolbar.addAction.call_args[0][0]
        assert isinstance(added, QAction)


# ---------------------------------------------------------------------------
# integrate_plugin_export_actions
# ---------------------------------------------------------------------------

class TestIntegratePluginExportActions:
    def test_no_export_actions_does_nothing(self, im, pmm):
        im.host.plugin_manager.export_actions = []
        pmm.integrate_plugin_export_actions()
        im.export_button.menu.return_value.addSeparator.assert_not_called()

    def test_adds_actions_to_export_button_menu(self, im, pmm):
        im.host.plugin_manager.export_actions = [
            {"label": "Export CSV", "callback": MagicMock()}
        ]
        im.host.menuBar.return_value.actions.return_value = []  # no File menu

        pmm.integrate_plugin_export_actions()

        export_menu = im.export_button.menu.return_value
        export_menu.addSeparator.assert_called_once()
        assert export_menu.addAction.called

    def test_adds_actions_to_both_export_button_and_file_menu(self, im, pmm):
        im.host.plugin_manager.export_actions = [
            {"label": "Export PDF", "callback": MagicMock()}
        ]

        export_submenu = MagicMock()
        export_submenu.actions.return_value = []
        export_sub_action = MagicMock()
        export_sub_action.menu.return_value = export_submenu
        export_sub_action.text.return_value = "Export"

        file_menu = MagicMock()
        file_menu.actions.return_value = [export_sub_action]
        file_top_action = MagicMock()
        file_top_action.menu.return_value = file_menu
        file_top_action.text.return_value = "File"

        im.host.menuBar.return_value.actions.return_value = [file_top_action]

        pmm.integrate_plugin_export_actions()

        # Both menus received a separator + action
        assert im.export_button.menu.return_value.addAction.called
        assert export_submenu.addAction.called


# ---------------------------------------------------------------------------
# integrate_plugin_analysis_tools
# ---------------------------------------------------------------------------

class TestIntegratePluginAnalysisTools:
    def test_no_analysis_menu_does_nothing(self, im, pmm):
        im.host.menuBar.return_value.actions.return_value = []
        im.host.plugin_manager.analysis_tools = [
            {"label": "NMR", "callback": MagicMock(), "plugin": "NMRPlugin"}
        ]
        pmm.integrate_plugin_analysis_tools()  # must not raise

    def test_no_tools_skips_separator(self, im, pmm):
        analysis_menu = MagicMock()
        analysis_action = MagicMock()
        analysis_action.text.return_value = "Analysis"
        analysis_action.menu.return_value = analysis_menu
        im.host.menuBar.return_value.actions.return_value = [analysis_action]
        im.host.plugin_manager.analysis_tools = []

        pmm.integrate_plugin_analysis_tools()

        analysis_menu.addSeparator.assert_not_called()

    def test_adds_tools_to_analysis_menu(self, im, pmm):
        analysis_menu = MagicMock()
        analysis_action = MagicMock()
        analysis_action.text.return_value = "Analysis"
        analysis_action.menu.return_value = analysis_menu
        im.host.menuBar.return_value.actions.return_value = [analysis_action]
        im.host.plugin_manager.analysis_tools = [
            {"label": "NMR Predict", "callback": MagicMock(), "plugin": "NMRPlugin"}
        ]

        pmm.integrate_plugin_analysis_tools()

        analysis_menu.addSeparator.assert_called_once()
        analysis_menu.addAction.assert_called_once()


# ---------------------------------------------------------------------------
# update_style_menu_with_plugins
# ---------------------------------------------------------------------------

class TestUpdateStyleMenuWithPlugins:
    def test_no_style_button_does_nothing(self):
        im = make_init_manager(has_style_button=False)
        pmm = PluginMenuManager(im)
        pmm.update_style_menu_with_plugins()  # must not raise

    def test_no_custom_styles_does_nothing(self, im, pmm):
        im.host.plugin_manager.custom_3d_styles = []
        pmm.update_style_menu_with_plugins()
        im.style_button.menu.return_value.addAction.assert_not_called()

    def test_adds_custom_style_actions(self, im, pmm):
        im.host.plugin_manager.custom_3d_styles = ["Ball", "Licorice"]

        style_menu = MagicMock()
        existing_action = MagicMock()
        existing_action.actionGroup.return_value = MagicMock()  # truthy group
        existing_action.isSeparator.return_value = False
        existing_action.text.return_value = "Existing"
        style_menu.actions.return_value = [existing_action]
        style_menu.addAction = MagicMock()
        im.style_button.menu.return_value = style_menu

        pmm.update_style_menu_with_plugins()

        assert style_menu.addAction.call_count == 2  # two new styles
        assert style_menu.addSeparator.called  # separator before new styles

    def test_does_not_duplicate_existing_style(self, im, pmm):
        im.host.plugin_manager.custom_3d_styles = ["Ball"]

        style_menu = MagicMock()
        group_action = MagicMock()
        group_action.actionGroup.return_value = MagicMock()
        group_action.isSeparator.return_value = False

        existing_style_action = MagicMock()
        existing_style_action.actionGroup.return_value = None
        existing_style_action.text.return_value = "Ball"
        existing_style_action.isSeparator.return_value = False

        style_menu.actions.return_value = [group_action, existing_style_action]
        im.style_button.menu.return_value = style_menu

        pmm.update_style_menu_with_plugins()

        style_menu.addAction.assert_not_called()


# ---------------------------------------------------------------------------
# integrate_plugin_file_openers
# ---------------------------------------------------------------------------

class TestIntegratePluginFileOpeners:
    def test_no_file_openers_does_nothing(self, im, pmm):
        im.host.plugin_manager.file_openers = {}
        pmm.integrate_plugin_file_openers()
        im.import_menu.addSeparator.assert_not_called()

    def test_no_import_menu_does_nothing(self):
        im = make_init_manager(has_import_menu=False)
        pm = MagicMock()
        pm.file_openers = {".xyz": [{"callback": MagicMock(), "plugin": "XYZPlugin"}]}
        im.host.plugin_manager = pm
        pmm = PluginMenuManager(im)
        pmm.integrate_plugin_file_openers()  # must not raise

    def test_adds_opener_actions(self, im, pmm):
        cb = MagicMock()
        im.host.plugin_manager.file_openers = {
            ".xyz": [{"callback": cb, "plugin": "XYZPlugin"}]
        }

        pmm.integrate_plugin_file_openers()

        im.import_menu.addSeparator.assert_called_once()
        im.import_menu.addAction.assert_called_once()


# ---------------------------------------------------------------------------
# _clear_all_plugin_actions
# ---------------------------------------------------------------------------

class TestClearAllPluginActions:
    def test_clears_plugin_menu(self, im, pmm):
        plugin_menu = MagicMock(spec=QMenu)
        plugin_menu.actions.return_value = []
        im.host.menuBar.return_value.actions.return_value = []
        pmm._clear_all_plugin_actions(plugin_menu)
        plugin_menu.clear.assert_called_once()

    def test_removes_tagged_actions_from_all_menus(self, im, pmm):
        tagged = MagicMock()
        tagged.data.return_value = "plugin_managed"
        tagged.menu.return_value = None

        sub_menu = MagicMock()
        sub_menu.actions.return_value = [tagged]
        sub_action = MagicMock()
        sub_action.data.return_value = None
        sub_action.menu.return_value = sub_menu

        top_menu = MagicMock()
        top_menu.actions.return_value = [sub_action]
        top_action = MagicMock()
        top_action.menu.return_value = top_menu

        im.host.menuBar.return_value.actions.return_value = [top_action]
        plugin_menu = MagicMock(spec=QMenu)
        plugin_menu.actions.return_value = []

        pmm._clear_all_plugin_actions(plugin_menu)

        sub_menu.removeAction.assert_called_with(tagged)
