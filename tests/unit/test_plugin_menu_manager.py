"""Unit tests for PluginMenuManager — plugin UI lifecycle management."""

import pytest
from typing import Optional
from unittest.mock import MagicMock
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QMenu

from moleditpy.ui.plugin_menu_manager import PluginMenuManager


# ---------------------------------------------------------------------------
# Module-level patch: QAction(text, MagicMock) fails because PyQt6 strictly
# validates the parent type. Strip the mock parent so tests stay lightweight.
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _patch_qaction(monkeypatch, app):
    """Allow QAction to be constructed with a MagicMock host in unit tests.

    Parent the actions to a real QObject holder (not None): unparented QActions
    get GC'd after QApplication teardown, causing an access-violation crash on
    Windows when this file runs in isolation. The holder gives them an owner
    destroyed in order while the app is alive.
    """
    from PyQt6.QtCore import QObject
    from PyQt6.QtGui import QAction as _RealQAction

    holder = QObject()

    def _safe_qaction(text: str, parent=None) -> QAction:
        return _RealQAction(text, holder)

    monkeypatch.setattr("moleditpy.ui.plugin_menu_manager.QAction", _safe_qaction)
    yield
    holder.deleteLater()


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def make_init_manager(
    *,
    has_toolbar: bool = True,
    has_export_button: bool = True,
    has_style_button: bool = True,
    has_import_menu: bool = True,
    plugin_manager: Optional[MagicMock] = None,
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
        """PluginMenuManager stores a reference to the init manager it received."""
        mgr = PluginMenuManager(im)
        assert mgr._im is im


# ---------------------------------------------------------------------------
# update_plugin_menu
# ---------------------------------------------------------------------------


class TestUpdatePluginMenu:
    def test_does_nothing_when_no_plugin_manager(self, im, pmm):
        """update_plugin_menu does nothing when plugin_manager is None."""
        im.host.plugin_manager = None
        menu = MagicMock(spec=QMenu)
        pmm.update_plugin_menu(menu)
        menu.clear.assert_not_called()

    def test_clears_and_adds_manage_action(self, im, pmm):
        """update_plugin_menu clears the menu and adds at least the manage action."""
        menu = MagicMock(spec=QMenu)
        pmm.update_plugin_menu(menu)
        menu.clear.assert_called_once()
        # addAction called at least once (the "Plugin Manager..." entry)
        assert menu.addAction.called or menu.addSeparator.called

    def test_discover_plugins_called(self, im, pmm):
        """update_plugin_menu calls discover_plugins on the plugin_manager."""
        menu = MagicMock(spec=QMenu)
        im.host.plugin_manager.discover_plugins.return_value = []
        pmm.update_plugin_menu(menu)
        im.host.plugin_manager.discover_plugins.assert_called_once_with(im.host)

    def test_menu_action_appears_after_update(self, im):
        """A registered menu action is present in the menu after update_plugin_menu."""
        pmm = PluginMenuManager(im)
        im.host.plugin_manager.discover_plugins.return_value = []
        im.host.plugin_manager.menu_actions = [
            {
                "path": "MyPlugin/DoThing",
                "callback": MagicMock(),
                "text": "Do Thing",
                "shortcut": None,
            }
        ]
        im.host.menuBar.return_value.actions.return_value = []
        added_menu = MagicMock(**{"actions.return_value": []})
        im.host.menuBar.return_value.addMenu.return_value = added_menu

        pmm.update_plugin_menu(MagicMock(spec=QMenu))

        im.host.menuBar.return_value.addMenu.assert_called_once_with("MyPlugin")
        added_menu.addAction.assert_called_once()


# ---------------------------------------------------------------------------
# rebuild_plugin_menus
# ---------------------------------------------------------------------------


class TestRebuildPluginMenus:
    def test_toolbar_and_export_populated_after_rebuild(self, im):
        """rebuild_plugin_menus wires toolbar and export actions into UI."""
        pmm = PluginMenuManager(im)
        im.host.menuBar.return_value.actions.return_value = []
        im.host.plugin_manager.toolbar_actions = [
            {"text": "Run", "callback": MagicMock(), "icon": "", "tooltip": ""}
        ]
        im.host.plugin_manager.export_actions = [
            {"label": "Export CSV", "callback": MagicMock()}
        ]

        pmm.rebuild_plugin_menus()

        im.plugin_toolbar.addAction.assert_called_once()
        im.export_button.menu.return_value.addAction.assert_called_once()

    def test_resets_separator_flag(self, im, pmm):
        """rebuild_plugin_menus resets the plugin_menubar_separator_added flag."""
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
        """add_registered_plugin_actions does nothing when no menu actions are registered."""
        im.host.plugin_manager.menu_actions = []
        pmm.add_registered_plugin_actions()
        im.host.menuBar.return_value.addMenu.assert_not_called()

    def test_creates_new_top_level_menu_with_separator(self, im, pmm):
        """add_registered_plugin_actions creates a new top-level menu and adds a separator."""
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {
                "path": "MyPlugin/Action",
                "callback": callback,
                "text": "Run It",
                "shortcut": None,
            }
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
        """add_registered_plugin_actions reuses an existing top-level menu rather than creating a new one."""
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {
                "path": "File/ExportXYZ",
                "callback": callback,
                "text": "Export XYZ",
                "shortcut": None,
            }
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
            {
                "path": "PlugA/Action1",
                "callback": MagicMock(),
                "text": "A1",
                "shortcut": None,
            },
            {
                "path": "PlugB/Action2",
                "callback": MagicMock(),
                "text": "B1",
                "shortcut": None,
            },
        ]
        im.host.menuBar.return_value.actions.return_value = []
        im.host.menuBar.return_value.addMenu.return_value = MagicMock(
            **{"actions.return_value": []}
        )

        pmm.add_registered_plugin_actions()

        assert im.host.menuBar.return_value.addSeparator.call_count == 1

    def test_shortcut_applied_when_present(self, im, pmm):
        """add_registered_plugin_actions applies a keyboard shortcut when one is specified."""
        callback = MagicMock()
        im.host.plugin_manager.menu_actions = [
            {
                "path": "Plug/Act",
                "callback": callback,
                "text": "Act",
                "shortcut": "Ctrl+P",
            }
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
        """add_plugin_toolbar_actions does not raise when plugin_toolbar is missing."""
        im = make_init_manager(has_toolbar=False)
        pmm = PluginMenuManager(im)
        pmm.add_plugin_toolbar_actions()  # must not raise

    def test_hides_toolbar_when_no_actions(self, im, pmm):
        """add_plugin_toolbar_actions hides and clears the toolbar when no actions are registered."""
        im.host.plugin_manager.toolbar_actions = []
        pmm.add_plugin_toolbar_actions()
        im.plugin_toolbar.hide.assert_called_once()
        im.plugin_toolbar.clear.assert_called_once()

    def test_shows_toolbar_and_adds_actions(self, im, pmm):
        """add_plugin_toolbar_actions shows the toolbar and adds an action for each registered entry."""
        im.host.plugin_manager.toolbar_actions = [
            {"text": "Run", "callback": MagicMock(), "icon": "", "tooltip": "Run it"},
        ]
        pmm.add_plugin_toolbar_actions()
        im.plugin_toolbar.show.assert_called_once()
        im.plugin_toolbar.addAction.assert_called_once()

    def test_icon_set_when_file_exists(self, im, pmm, tmp_path):
        """add_plugin_toolbar_actions creates a QAction with an icon when the icon file exists."""
        icon_file = tmp_path / "icon.png"
        icon_file.write_bytes(b"")  # create the file
        im.host.plugin_manager.toolbar_actions = [
            {
                "text": "Icon Action",
                "callback": MagicMock(),
                "icon": str(icon_file),
                "tooltip": "",
            }
        ]
        pmm.add_plugin_toolbar_actions()
        added = im.plugin_toolbar.addAction.call_args[0][0]
        assert isinstance(added, QAction)


# ---------------------------------------------------------------------------
# integrate_plugin_export_actions
# ---------------------------------------------------------------------------


class TestIntegratePluginExportActions:
    def test_no_export_actions_does_nothing(self, im, pmm):
        """integrate_plugin_export_actions does nothing when no export actions are registered."""
        im.host.plugin_manager.export_actions = []
        pmm.integrate_plugin_export_actions()
        im.export_button.menu.return_value.addSeparator.assert_not_called()

    def test_adds_actions_to_export_button_menu(self, im, pmm):
        """integrate_plugin_export_actions adds a separator and action to the export button menu."""
        im.host.plugin_manager.export_actions = [
            {"label": "Export CSV", "callback": MagicMock()}
        ]
        im.host.menuBar.return_value.actions.return_value = []  # no File menu

        pmm.integrate_plugin_export_actions()

        export_menu = im.export_button.menu.return_value
        export_menu.addSeparator.assert_called_once()
        assert export_menu.addAction.called

    def test_adds_actions_to_both_export_button_and_file_menu(self, im, pmm):
        """integrate_plugin_export_actions adds actions to both the export button and the File/Export submenu."""
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
        """integrate_plugin_analysis_tools does not raise when no Analysis menu exists."""
        im.host.menuBar.return_value.actions.return_value = []
        im.host.plugin_manager.analysis_tools = [
            {"label": "NMR", "callback": MagicMock(), "plugin": "NMRPlugin"}
        ]
        pmm.integrate_plugin_analysis_tools()  # must not raise

    def test_no_tools_skips_separator(self, im, pmm):
        """integrate_plugin_analysis_tools skips the separator when no tools are registered."""
        analysis_menu = MagicMock()
        analysis_action = MagicMock()
        analysis_action.text.return_value = "Analysis"
        analysis_action.menu.return_value = analysis_menu
        im.host.menuBar.return_value.actions.return_value = [analysis_action]
        im.host.plugin_manager.analysis_tools = []

        pmm.integrate_plugin_analysis_tools()

        analysis_menu.addSeparator.assert_not_called()

    def test_adds_tools_to_analysis_menu(self, im, pmm):
        """integrate_plugin_analysis_tools adds a separator and action to the Analysis menu."""
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
        """update_style_menu_with_plugins does not raise when style_button is None."""
        im = make_init_manager(has_style_button=False)
        pmm = PluginMenuManager(im)
        pmm.update_style_menu_with_plugins()  # must not raise

    def test_no_custom_styles_does_nothing(self, im, pmm):
        """update_style_menu_with_plugins does nothing when no custom styles are registered."""
        im.host.plugin_manager.custom_3d_styles = []
        pmm.update_style_menu_with_plugins()
        im.style_button.menu.return_value.addAction.assert_not_called()

    def test_adds_custom_style_actions(self, im, pmm):
        """update_style_menu_with_plugins adds a separator and one action per custom style."""
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
        """update_style_menu_with_plugins skips a style already present in the menu."""
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
        """integrate_plugin_file_openers does nothing when no file openers are registered."""
        im.host.plugin_manager.file_openers = {}
        pmm.integrate_plugin_file_openers()
        im.import_menu.addSeparator.assert_not_called()

    def test_no_import_menu_does_nothing(self):
        """integrate_plugin_file_openers does not raise when import_menu is None."""
        im = make_init_manager(has_import_menu=False)
        pm = MagicMock()
        pm.file_openers = {".xyz": [{"callback": MagicMock(), "plugin": "XYZPlugin"}]}
        im.host.plugin_manager = pm
        pmm = PluginMenuManager(im)
        pmm.integrate_plugin_file_openers()  # must not raise

    def test_adds_opener_actions(self, im, pmm):
        """integrate_plugin_file_openers adds a separator and action for each file opener."""
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
    def test_reset_plugin_menu_clears_and_re_adds_the_header(self, im, pmm):
        """_reset_plugin_menu empties the Plugin menu and restores its header."""
        plugin_menu = MagicMock(spec=QMenu)
        plugin_menu.actions.return_value = []

        pmm._reset_plugin_menu(plugin_menu)

        plugin_menu.clear.assert_called_once()
        plugin_menu.addAction.assert_called_once()
        plugin_menu.addSeparator.assert_called_once()

    def test_removes_tagged_actions_from_all_menus(self, im, pmm):
        """_clear_all_plugin_actions removes all plugin_managed tagged actions from all menus."""
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

        pmm._clear_all_plugin_actions()

        sub_menu.removeAction.assert_called_with(tagged)


# ---------------------------------------------------------------------------
# _add_legacy_plugin_actions — categorized run()-plugin menu building
# ---------------------------------------------------------------------------

import os


def _legacy_plugin(name, category=None):
    mod = MagicMock()  # MagicMock has a `run` attribute -> treated as legacy
    p = {"name": name, "module": mod}
    if category is not None:
        p["category"] = category
    return p


class TestAddLegacyPluginActions:
    def test_empty_list_adds_disabled_placeholder(self, im, pmm):
        menu = QMenu()
        pmm._add_legacy_plugin_actions(menu, [])
        actions = menu.actions()
        assert len(actions) == 1
        assert actions[0].text() == "(No plugins found)"
        assert not actions[0].isEnabled()

    def test_root_plugins_added_sorted_by_name(self, im, pmm):
        menu = QMenu()
        pmm._add_legacy_plugin_actions(
            menu, [_legacy_plugin("Zeta"), _legacy_plugin("Alpha")]
        )
        texts = [a.text() for a in menu.actions()]
        assert texts == ["Alpha", "Zeta"]

    def test_categorized_plugin_goes_into_submenu(self, im, pmm):
        menu = QMenu()
        pmm._add_legacy_plugin_actions(menu, [_legacy_plugin("Doer", "Tools")])
        submenus = [a.menu() for a in menu.actions() if a.menu()]
        assert len(submenus) == 1
        assert submenus[0].title() == "Tools"
        assert [a.text() for a in submenus[0].actions()] == ["Doer"]

    def test_nested_category_builds_nested_submenus(self, im, pmm):
        menu = QMenu()
        pmm._add_legacy_plugin_actions(
            menu, [_legacy_plugin("Deep", f"Outer{os.sep}Inner")]
        )
        outer = next(a.menu() for a in menu.actions() if a.menu())
        assert outer.title() == "Outer"
        inner = next(a.menu() for a in outer.actions() if a.menu())
        assert inner.title() == "Inner"
        assert [a.text() for a in inner.actions()] == ["Deep"]

    def test_same_category_reuses_single_submenu(self, im, pmm):
        menu = QMenu()
        pmm._add_legacy_plugin_actions(
            menu,
            [_legacy_plugin("Beta", "Tools"), _legacy_plugin("Alpha", "Tools")],
        )
        submenus = [a.menu() for a in menu.actions() if a.menu()]
        assert len(submenus) == 1  # one shared submenu, not two
        assert [a.text() for a in submenus[0].actions()] == ["Alpha", "Beta"]

    def test_triggering_action_runs_plugin(self, im, pmm):
        menu = QMenu()
        plugin = _legacy_plugin("Runner")
        pmm._add_legacy_plugin_actions(menu, [plugin])
        menu.actions()[0].trigger()
        im.host.plugin_manager.run_plugin.assert_called_once_with(
            plugin["module"], im.host
        )


# ---------------------------------------------------------------------------
# add_registered_plugin_actions — divider placement in tree-structured menus
# ---------------------------------------------------------------------------

TAG = PluginMenuManager._PLUGIN_ACTION_TAG


def _menubar_with_tools(im, app):
    """Give *im* a real menu bar holding a native 'Tools' menu with two entries."""
    from PyQt6.QtWidgets import QMenuBar

    bar = QMenuBar()
    tools = bar.addMenu("Tools")
    tools.addAction("Native A")
    tools.addAction("Native B")
    im.host.menuBar = MagicMock(return_value=bar)
    return bar, tools


def _entry(path, text):
    return {
        "plugin": "P",
        "path": path,
        "callback": lambda: None,
        "text": text,
        "icon": None,
        "shortcut": None,
    }


class TestPluginSeparatorPlacement:
    def test_divider_precedes_a_nested_plugin_submenu(self, im, pmm, app):
        """A nested path puts the divider before the submenu, not after it.

        The divider used to be checked against the freshly created (empty) leaf
        submenu, so no divider was added at all and the plugin submenu sat
        flush against the native entries.
        """
        _bar, tools = _menubar_with_tools(im, app)
        im.host.plugin_manager.menu_actions = [_entry("Tools/Sub/Item", "Item")]

        pmm.add_registered_plugin_actions()

        kinds = [
            "sep" if a.isSeparator() else ("menu" if a.menu() else a.text())
            for a in tools.actions()
        ]
        assert kinds == ["Native A", "Native B", "sep", "menu"]

    def test_no_extra_divider_between_plugin_siblings(self, im, pmm, app):
        """An entry following a plugin submenu does not get its own divider.

        The submenu action is tagged, so the dedupe check recognises it as
        plugin-owned instead of mistaking it for a native entry.
        """
        _bar, tools = _menubar_with_tools(im, app)
        im.host.plugin_manager.menu_actions = [
            _entry("Tools/Sub/Item", "Item"),
            _entry("Tools/Direct", "Direct"),
        ]

        pmm.add_registered_plugin_actions()

        kinds = [
            "sep" if a.isSeparator() else ("menu" if a.menu() else a.text())
            for a in tools.actions()
        ]
        assert kinds == ["Native A", "Native B", "sep", "menu", "Direct"]

    def test_plugin_created_submenu_is_tagged(self, im, pmm, app):
        """Submenus built for a plugin path carry the plugin tag for cleanup."""
        _bar, tools = _menubar_with_tools(im, app)
        im.host.plugin_manager.menu_actions = [_entry("Tools/Sub/Item", "Item")]

        pmm.add_registered_plugin_actions()

        sub_action = next(a for a in tools.actions() if a.menu())
        assert sub_action.data() == TAG

    def test_shared_submenu_gets_one_divider_and_no_duplicate(self, im, pmm, app):
        """Two plugins under the same submenu reuse it without extra dividers."""
        _bar, tools = _menubar_with_tools(im, app)
        im.host.plugin_manager.menu_actions = [
            _entry("Tools/Sub/One", "One"),
            _entry("Tools/Sub/Two", "Two"),
        ]

        pmm.add_registered_plugin_actions()

        submenus = [a.menu() for a in tools.actions() if a.menu()]
        assert len(submenus) == 1
        assert [a.text() for a in submenus[0].actions()] == ["One", "Two"]
        assert sum(1 for a in tools.actions() if a.isSeparator()) == 1


# ---------------------------------------------------------------------------
# Menu-bar cleanup — shared by both clear paths
# ---------------------------------------------------------------------------


def _menubar_with_plugin_top_level(im):
    """Build a bar holding a native menu plus a tagged divider and plugin menu."""
    from PyQt6.QtWidgets import QMenuBar

    bar = QMenuBar()
    bar.addMenu("Tools").addAction("Native A")

    sep = bar.addSeparator()
    sep.setData(TAG)

    plugin_top = bar.addMenu("MyPlugin")
    plugin_top.menuAction().setData(TAG)
    entry = plugin_top.addAction("Gone")
    entry.setData(TAG)

    im.host.menuBar = MagicMock(return_value=bar)
    im.plugin_menubar_separator_added = True
    return bar


class TestClearPluginMenubarEntries:
    def test_clear_all_drops_stale_plugin_top_level_menu(self, im, pmm, app):
        """_clear_all_plugin_actions reclaims the menu bar, not just the menus.

        _get_plugin_target_menus walks into the bar's menus but never the bar
        itself, so an uninstalled plugin's own top-level menu used to survive
        as an empty stub with an orphaned divider in front of it.
        """
        bar = _menubar_with_plugin_top_level(im)

        pmm._clear_all_plugin_actions()

        titles = [a.text() for a in bar.actions() if a.menu()]
        assert titles == ["Tools"]
        assert not any(a.isSeparator() for a in bar.actions())

    def test_clear_all_resets_separator_flag(self, im, pmm, app):
        """The flag is reset so the next populate pass re-adds the divider."""
        _menubar_with_plugin_top_level(im)

        pmm._clear_all_plugin_actions()

        assert im.plugin_menubar_separator_added is False

    def test_native_top_level_menu_is_never_removed(self, im, pmm, app):
        """An untagged top-level menu survives even when it is empty."""
        from PyQt6.QtWidgets import QMenuBar

        bar = QMenuBar()
        bar.addMenu("Empty Native")
        im.host.menuBar = MagicMock(return_value=bar)

        pmm._clear_all_plugin_actions()

        assert [a.text() for a in bar.actions() if a.menu()] == ["Empty Native"]

    def test_rebuild_uses_the_same_menubar_cleanup(self, im, pmm, app):
        """rebuild_plugin_menus drops the stale top-level menu too."""
        bar = _menubar_with_plugin_top_level(im)

        pmm.rebuild_plugin_menus()

        assert [a.text() for a in bar.actions() if a.menu()] == ["Tools"]
        assert im.plugin_menubar_separator_added is False


# ---------------------------------------------------------------------------
# Plugin menu composition — context-injected vs folder-derived groups
# ---------------------------------------------------------------------------


def _menubar_with_plugin_menu(im):
    """Give *im* a real menu bar whose only menu is the Plugin menu."""
    from PyQt6.QtWidgets import QMenuBar

    bar = QMenuBar()
    plugin_menu = bar.addMenu("&Plugin")
    im.host.menuBar = MagicMock(return_value=bar)
    im.plugin_menu = plugin_menu
    return bar, plugin_menu


def _shape(menu):
    """Describe a menu as a list of 'sep' / submenu titles / action labels."""
    return [
        "sep" if a.isSeparator() else (a.menu().title() + "/" if a.menu() else a.text())
        for a in menu.actions()
    ]


class TestPluginMenuComposition:
    def test_update_puts_context_entries_above_folder_ones(self, im, pmm, app):
        """A fresh build lists context-injected entries, a divider, then folders."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_entry("Plugin/Injected", "Injected")]
        im.host.plugin_manager.discover_plugins.return_value = [
            _legacy_plugin("Filed", "Folder")
        ]

        pmm.update_plugin_menu(plugin_menu)

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "sep",
            "Injected",
            "sep",
            "Folder/",
        ]

    def test_rebuild_preserves_that_order(self, im, pmm, app):
        """A reload must not flip the two groups.

        rebuild_plugin_menus used to strip only the tagged (context-injected)
        entries and append the fresh ones after the surviving folder entries,
        so every reload swapped the order.
        """
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_entry("Plugin/Injected", "Injected")]
        im.host.plugin_manager.plugins = [_legacy_plugin("Filed", "Folder")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "sep",
            "Injected",
            "sep",
            "Folder/",
        ]

    def test_repeated_rebuilds_are_stable(self, im, pmm, app):
        """Reloading twice yields the same menu, with no accumulation."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_entry("Plugin/Injected", "Injected")]
        im.host.plugin_manager.plugins = [_legacy_plugin("Filed", "Folder")]

        pmm.rebuild_plugin_menus()
        first = _shape(plugin_menu)
        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == first

    def test_no_trailing_divider_without_folder_plugins(self, im, pmm, app):
        """With only context-injected entries the menu ends on a real entry."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_entry("Plugin/Injected", "Injected")]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == ["Plugin Manager...", "sep", "Injected"]

    def test_folder_plugins_alone_get_no_leading_divider(self, im, pmm, app):
        """Without context-injected entries the header divider is enough."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = []
        im.host.plugin_manager.plugins = [_legacy_plugin("Filed", "Folder")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == ["Plugin Manager...", "sep", "Folder/"]


# ---------------------------------------------------------------------------
# Plugin Installer — pinned beside the Plugin Manager
# ---------------------------------------------------------------------------


def _v4_plugin(name):
    """A discovered initialize()-style plugin: present, but with no run()."""
    p = _legacy_plugin(name)
    del p["module"].run
    return p


def _installer_entry(plugin="Plugin Installer", label="Plugin Installer..."):
    entry = _entry("Plugin/" + label, label)
    entry["plugin"] = plugin
    return entry


class TestPluginInstallerPinning:
    def test_installer_sits_directly_below_the_manager(self, im, pmm, app):
        """The installer joins the header pair with no divider between them."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_installer_entry()]
        im.host.plugin_manager.plugins = [_legacy_plugin("Filed", "Folder")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "Plugin Installer...",
            "sep",
            "Folder/",
        ]

    def test_installer_stays_above_other_injected_plugins(self, im, pmm, app):
        """Ordinary context-injected entries stay below the divider."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [
            _entry("Plugin/Injected", "Injected"),
            _installer_entry(),
        ]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "Plugin Installer...",
            "sep",
            "Injected",
        ]

    def test_registration_order_does_not_matter(self, im, pmm, app):
        """Registering the installer first gives the same layout."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [
            _installer_entry(),
            _entry("Plugin/Injected", "Injected"),
        ]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "Plugin Installer...",
            "sep",
            "Injected",
        ]

    def test_installer_alone_leaves_no_dangling_divider(self, im, pmm, app):
        """With nothing under it, the header divider is dropped."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [_installer_entry()]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == ["Plugin Manager...", "Plugin Installer..."]

    def test_matched_case_insensitively_by_label(self, im, pmm, app):
        """A differently-cased PLUGIN_NAME or label still pins."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [
            _installer_entry(plugin="Some Author Bundle", label="Plugin installer")
        ]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == ["Plugin Manager...", "Plugin installer"]

    def test_unrelated_plugin_is_not_pinned(self, im, pmm, app):
        """A plugin merely mentioning install is left below the divider."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        im.host.plugin_manager.menu_actions = [
            _installer_entry(plugin="Installer Helper", label="Install Helper")
        ]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == [
            "Plugin Manager...",
            "sep",
            "Install Helper",
        ]

    def test_nested_installer_path_is_not_pinned(self, im, pmm, app):
        """Only a direct Plugin/<entry> registration is pinned."""
        _bar, plugin_menu = _menubar_with_plugin_menu(im)
        entry = _entry("Plugin/Tools/Plugin Installer...", "Plugin Installer...")
        entry["plugin"] = "Plugin Installer"
        im.host.plugin_manager.menu_actions = [entry]
        im.host.plugin_manager.plugins = [_v4_plugin("Modern")]

        pmm.rebuild_plugin_menus()

        assert _shape(plugin_menu) == ["Plugin Manager...", "sep", "Tools/"]
