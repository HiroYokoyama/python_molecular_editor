"""Integration tests for plugin menu reload completeness.

Verifies that PluginMenuManager.rebuild_plugin_menus fully rebuilds all
six integration points: menus, toolbar, export, file-openers, analysis, and
3D styles — and that MainInitManager correctly delegates to it.
"""

import pytest
from typing import Optional
from unittest.mock import MagicMock

from moleditpy.ui.plugin_menu_manager import PluginMenuManager


# QAction(text, MagicMock) fails because PyQt6 strictly validates parent type.
# Patch it to parent actions to a real QObject instead of the mock. Parenting
# (rather than dropping the parent) gives the QActions a deterministic owner, so
# Qt destroys them in order with the holder while QApplication is still alive —
# leaving them unparented causes an access-violation crash at GC on Windows when
# this file runs in isolation.
@pytest.fixture(autouse=True)
def _patch_qaction(monkeypatch, app):
    from PyQt6.QtCore import QObject
    from PyQt6.QtGui import QAction as _RealQAction

    holder = QObject()

    def _safe_qaction(text: str, parent=None) -> _RealQAction:
        return _RealQAction(text, holder)

    monkeypatch.setattr("moleditpy.ui.plugin_menu_manager.QAction", _safe_qaction)
    yield
    holder.deleteLater()


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_plugin_manager(
    *,
    menu_actions=None,
    toolbar_actions=None,
    export_actions=None,
    file_openers=None,
    analysis_tools=None,
    custom_3d_styles=None,
) -> MagicMock:
    pm = MagicMock()
    pm.menu_actions = menu_actions or []
    pm.toolbar_actions = toolbar_actions or []
    pm.export_actions = export_actions or []
    pm.file_openers = file_openers or {}
    pm.analysis_tools = analysis_tools or []
    pm.custom_3d_styles = custom_3d_styles or []
    pm.optimization_methods = {}
    return pm


def _make_im(plugin_manager: Optional[MagicMock] = None) -> MagicMock:
    im = MagicMock()
    im.plugin_menubar_separator_added = False
    im.host.menuBar.return_value.actions.return_value = []
    im.host.plugin_manager = plugin_manager or _make_plugin_manager()
    im.plugin_toolbar = MagicMock()
    im.export_button = MagicMock()
    im.export_button.menu.return_value = MagicMock()
    im.import_menu = MagicMock()
    im.style_button = MagicMock()
    im.style_button.menu.return_value = MagicMock()
    return im


# ---------------------------------------------------------------------------
# rebuild_plugin_menus — completeness
# ---------------------------------------------------------------------------


class TestRebuildCompleteness:
    """rebuild_plugin_menus must always call all six integration points."""

    def test_all_six_steps_run_when_all_registries_empty(self):
        """Even with empty registries, all six rebuild methods execute."""
        im = _make_im()
        pmm = PluginMenuManager(im)

        called = []

        def track(name):
            def _fn():
                called.append(name)

            return _fn

        pmm.add_registered_plugin_actions = track("menu_actions")
        pmm.add_plugin_toolbar_actions = track("toolbar_actions")
        pmm.integrate_plugin_export_actions = track("export_actions")
        pmm.integrate_plugin_file_openers = track("file_openers")
        pmm.integrate_plugin_analysis_tools = track("analysis_tools")
        pmm.update_style_menu_with_plugins = track("style_menu")
        pmm.integrate_plugin_optimization_methods = track("optimization_methods")

        pmm.rebuild_plugin_menus()

        assert called == [
            "menu_actions",
            "toolbar_actions",
            "export_actions",
            "file_openers",
            "analysis_tools",
            "style_menu",
            "optimization_methods",
        ], f"Expected all six steps in order, got: {called}"

    def test_all_six_steps_run_when_all_registries_populated(self):
        """With all registries populated, all six steps still run."""
        pm = _make_plugin_manager(
            menu_actions=[
                {"path": "P/A", "callback": MagicMock(), "text": "A", "shortcut": None}
            ],
            toolbar_actions=[
                {"text": "TB", "callback": MagicMock(), "icon": "", "tooltip": ""}
            ],
            export_actions=[{"label": "EXP", "callback": MagicMock()}],
            file_openers={".xyz": [{"callback": MagicMock(), "plugin": "Xyz"}]},
            analysis_tools=[{"label": "NMR", "callback": MagicMock(), "plugin": "NMR"}],
            custom_3d_styles=["NewStyle"],
        )
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        called = []

        def track(name):
            def _fn():
                called.append(name)

            return _fn

        pmm.add_registered_plugin_actions = track("menu_actions")
        pmm.add_plugin_toolbar_actions = track("toolbar_actions")
        pmm.integrate_plugin_export_actions = track("export_actions")
        pmm.integrate_plugin_file_openers = track("file_openers")
        pmm.integrate_plugin_analysis_tools = track("analysis_tools")
        pmm.update_style_menu_with_plugins = track("style_menu")

        pmm.rebuild_plugin_menus()

        assert set(called) == {
            "menu_actions",
            "toolbar_actions",
            "export_actions",
            "file_openers",
            "analysis_tools",
            "style_menu",
        }

    def test_separator_flag_reset_before_rebuild(self):
        """rebuild_plugin_menus resets the separator flag to False before running."""
        im = _make_im()
        im.plugin_menubar_separator_added = True
        pmm = PluginMenuManager(im)
        pmm.rebuild_plugin_menus()
        assert im.plugin_menubar_separator_added is False

    def test_step_exception_does_not_abort_remaining_steps(self):
        """If one step raises, all subsequent steps still execute."""
        im = _make_im()
        pmm = PluginMenuManager(im)

        survived = []

        def step_ok(name):
            def _fn():
                survived.append(name)

            return _fn

        pmm.add_registered_plugin_actions = lambda: (_ for _ in ()).throw(
            RuntimeError("boom")
        )
        pmm.add_plugin_toolbar_actions = step_ok("toolbar")
        pmm.integrate_plugin_export_actions = step_ok("export")
        pmm.integrate_plugin_file_openers = step_ok("file_openers")
        pmm.integrate_plugin_analysis_tools = step_ok("analysis")
        pmm.update_style_menu_with_plugins = step_ok("style")

        pmm.rebuild_plugin_menus()

        assert "toolbar" in survived
        assert "export" in survived
        assert "style" in survived


class TestRebuildCleanup:
    """Verify that rebuild_plugin_menus cleans up existing tagged plugin actions."""

    def test_rebuild_cleans_tagged_actions_from_target_menus(self):
        """rebuild_plugin_menus must remove actions tagged with 'plugin_managed' from target menus."""
        im = _make_im()

        # Setup a mock menu for export button
        export_menu = MagicMock()
        im.export_button.menu.return_value = export_menu

        # Create a mock action that is tagged
        tagged_action = MagicMock()
        tagged_action.data.return_value = "plugin_managed"
        tagged_action.menu.return_value = None

        # Create a mock action that is NOT tagged
        untagged_action = MagicMock()
        untagged_action.data.return_value = None
        untagged_action.menu.return_value = None

        # Export menu has both actions
        export_menu.actions.return_value = [tagged_action, untagged_action]

        pmm = PluginMenuManager(im)
        pmm.rebuild_plugin_menus()

        # The tagged action should be removed
        export_menu.removeAction.assert_called_once_with(tagged_action)

    def test_rebuild_detaches_tagged_action_from_its_action_group(self):
        """A tagged action in a QActionGroup (e.g. a custom 3D style) must be
        removed from the group as well as the menu, otherwise it leaks into the
        exclusive group across rebuilds and can drop the active check state."""
        im = _make_im()

        style_menu = MagicMock()
        im.style_button.menu.return_value = style_menu

        group = MagicMock()
        tagged_action = MagicMock()
        tagged_action.data.return_value = "plugin_managed"
        tagged_action.menu.return_value = None
        tagged_action.actionGroup.return_value = group
        style_menu.actions.return_value = [tagged_action]

        pmm = PluginMenuManager(im)
        pmm.rebuild_plugin_menus()

        group.removeAction.assert_called_once_with(tagged_action)
        style_menu.removeAction.assert_called_once_with(tagged_action)

    def test_untagged_action_in_group_is_left_alone(self):
        """A non-plugin action sharing the style group must not be touched."""
        im = _make_im()

        style_menu = MagicMock()
        im.style_button.menu.return_value = style_menu

        group = MagicMock()
        builtin_action = MagicMock()
        builtin_action.data.return_value = None
        builtin_action.menu.return_value = None
        builtin_action.actionGroup.return_value = group
        style_menu.actions.return_value = [builtin_action]

        pmm = PluginMenuManager(im)
        pmm.rebuild_plugin_menus()

        group.removeAction.assert_not_called()
        style_menu.removeAction.assert_not_called()


# ---------------------------------------------------------------------------
# MainWindow.plugin_menu_manager proxy and PluginManager call path
# ---------------------------------------------------------------------------


class TestPluginMenuManagerRouting:
    """Verify the direct-routing architecture: callers go through
    MainWindow.plugin_menu_manager, not through MainInitManager wrappers."""

    def test_main_window_proxy_returns_init_manager_pmm(self):
        """MainWindow.plugin_menu_manager property reads from init_manager."""
        from moleditpy.ui.main_window import MainWindow

        # Read the property descriptor directly without constructing MainWindow
        # (QMainWindow requires a QApplication and full Qt setup to instantiate).
        prop = MainWindow.__dict__.get("plugin_menu_manager")
        assert prop is not None and isinstance(prop, property), (
            "MainWindow.plugin_menu_manager must be a @property"
        )
        # Verify the getter delegates to init_manager.plugin_menu_manager
        pmm_mock = MagicMock()
        fake_mw = MagicMock()
        fake_mw.init_manager.plugin_menu_manager = pmm_mock
        result = prop.fget(fake_mw)
        assert result is pmm_mock

    def test_plugin_manager_rebuild_calls_main_window_pmm(self):
        """PluginManager.rebuild_plugin_menus() calls main_window.plugin_menu_manager."""
        from moleditpy.plugins.plugin_manager import PluginManager

        pm = PluginManager.__new__(PluginManager)
        pmm_mock = MagicMock()
        pm.main_window = MagicMock()
        pm.main_window.plugin_menu_manager = pmm_mock

        pm.rebuild_plugin_menus()

        pmm_mock.rebuild_plugin_menus.assert_called_once()

    def test_plugin_manager_rebuild_no_main_window_does_nothing(self):
        """PluginManager.rebuild_plugin_menus() with no main_window is a no-op."""
        from moleditpy.plugins.plugin_manager import PluginManager

        pm = PluginManager.__new__(PluginManager)
        pm.main_window = None

        pm.rebuild_plugin_menus()  # must not raise

    def test_main_init_manager_has_no_wrapper_methods(self):
        """MainInitManager no longer exposes the old wrapper methods directly."""
        from moleditpy.ui.main_window_init import MainInitManager

        for method in (
            "rebuild_plugin_menus",
            "update_plugin_menu",
            "add_registered_plugin_actions",
            "add_plugin_toolbar_actions",
        ):
            assert not hasattr(MainInitManager, method), (
                f"MainInitManager should not have {method!r} — "
                "callers must use host.plugin_menu_manager directly"
            )


# ---------------------------------------------------------------------------
# update_plugin_menu — discovers and integrates
# ---------------------------------------------------------------------------


class TestUpdatePluginMenuIntegration:
    def test_discover_plugins_called_with_host(self):
        """update_plugin_menu calls discover_plugins with the init_manager host."""
        pm = _make_plugin_manager()
        pm.discover_plugins.return_value = []
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        plugin_menu = MagicMock()
        pmm.update_plugin_menu(plugin_menu)

        pm.discover_plugins.assert_called_once_with(im.host)

    def test_calls_seven_integration_methods(self):
        """update_plugin_menu invokes all seven integration methods exactly once."""
        pm = _make_plugin_manager()
        pm.discover_plugins.return_value = []
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        integration_methods = [
            "update_style_menu_with_plugins",
            "add_registered_plugin_actions",
            "add_plugin_toolbar_actions",
            "_add_legacy_plugin_actions",
            "integrate_plugin_export_actions",
            "integrate_plugin_file_openers",
            "integrate_plugin_analysis_tools",
            "integrate_plugin_optimization_methods",
        ]
        mocks = {m: MagicMock() for m in integration_methods}
        for m, mock in mocks.items():
            setattr(pmm, m, mock)

        pmm.update_plugin_menu(MagicMock())

        for m, mock in mocks.items():
            mock.assert_called_once(), f"Expected {m} to be called once"

    def test_does_nothing_when_plugin_manager_is_none(self):
        """update_plugin_menu is a no-op when plugin_manager is None."""
        im = _make_im()
        im.host.plugin_manager = None
        pmm = PluginMenuManager(im)
        plugin_menu = MagicMock()

        pmm.update_plugin_menu(plugin_menu)

        plugin_menu.clear.assert_not_called()


# ---------------------------------------------------------------------------
# Export action propagation end-to-end
# ---------------------------------------------------------------------------


class TestExportActionEndToEnd:
    """Export actions from the plugin manager appear in the export button menu."""

    def test_export_action_appears_in_button_menu(self):
        """A registered export action is added to the export button menu."""
        callback = MagicMock()
        pm = _make_plugin_manager(
            export_actions=[{"label": "Export XYZ", "callback": callback}]
        )
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        pmm.integrate_plugin_export_actions()

        export_menu = im.export_button.menu.return_value
        # A separator and at least one action were added
        export_menu.addSeparator.assert_called_once()
        export_menu.addAction.assert_called_once()
        added_action = export_menu.addAction.call_args[0][0]
        # The added action has the right label
        assert added_action.text() == "Export XYZ"

    def test_multiple_export_actions_all_added(self):
        """All registered export actions are each added to the export menu."""
        pm = _make_plugin_manager(
            export_actions=[
                {"label": "Export XYZ", "callback": MagicMock()},
                {"label": "Export CIF", "callback": MagicMock()},
            ]
        )
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        pmm.integrate_plugin_export_actions()

        assert im.export_button.menu.return_value.addAction.call_count == 2


# ---------------------------------------------------------------------------
# Analysis tool propagation end-to-end
# ---------------------------------------------------------------------------


class TestAnalysisToolEndToEnd:
    def test_analysis_tools_appear_in_analysis_menu(self):
        """Registered analysis tools appear as actions in the Analysis menu after rebuild."""
        callback = MagicMock()
        pm = _make_plugin_manager(
            analysis_tools=[
                {"label": "NMR Predict", "callback": callback, "plugin": "NMRPro"}
            ]
        )
        im = _make_im(pm)

        # Wire up an Analysis menu in the menubar
        analysis_menu = MagicMock()
        analysis_action = MagicMock()
        analysis_action.text.return_value = "Analysis"
        analysis_action.menu.return_value = analysis_menu
        im.host.menuBar.return_value.actions.return_value = [analysis_action]

        pmm = PluginMenuManager(im)
        pmm.integrate_plugin_analysis_tools()

        analysis_menu.addSeparator.assert_called_once()
        assert analysis_menu.addAction.call_count == 1
        added = analysis_menu.addAction.call_args[0][0]
        assert "NMR Predict" in added.text()
        assert "NMRPro" in added.text()


# ---------------------------------------------------------------------------
# rebuild_plugin_menus — plugin-created top-level menu removal (dev-4.3.1)
# ---------------------------------------------------------------------------


class TestTopLevelMenuCleanup:
    def _make_im_with_real_menubar(self, pm):
        from PyQt6.QtWidgets import QMenuBar

        bar = QMenuBar()
        im = _make_im(pm)
        im.host.menuBar = lambda: bar
        return im, bar

    def _titles(self, bar):
        return [a.text() for a in bar.actions() if a.menu() is not None]

    def test_emptied_plugin_toplevel_menu_removed_on_rebuild(self, app):
        """Regression: a plugin-created top-level menu stayed in the menu bar
        forever (dead and empty) after the plugin was uninstalled/reloaded."""
        pm = _make_plugin_manager(
            menu_actions=[
                {
                    "plugin": "TestPlugin",
                    "path": "MyPluginMenu/Do Thing",
                    "callback": lambda: None,
                    "text": "Do Thing",
                    "icon": None,
                    "shortcut": None,
                }
            ]
        )
        im, bar = self._make_im_with_real_menubar(pm)
        mgr = PluginMenuManager(im)

        mgr.add_registered_plugin_actions()
        assert "MyPluginMenu" in self._titles(bar)

        pm.menu_actions = []  # plugin uninstalled, registries re-discovered
        mgr.rebuild_plugin_menus()

        assert "MyPluginMenu" not in self._titles(bar)
        # No orphaned menu-bar separator either
        assert bar.actions() == []

    def test_still_registered_toplevel_menu_survives_rebuild(self, app):
        pm = _make_plugin_manager(
            menu_actions=[
                {
                    "plugin": "TestPlugin",
                    "path": "MyPluginMenu/Do Thing",
                    "callback": lambda: None,
                    "text": "Do Thing",
                    "icon": None,
                    "shortcut": None,
                }
            ]
        )
        im, bar = self._make_im_with_real_menubar(pm)
        mgr = PluginMenuManager(im)

        mgr.add_registered_plugin_actions()
        mgr.rebuild_plugin_menus()

        assert self._titles(bar) == ["MyPluginMenu"]
        menu = bar.actions()[-1].menu()
        labels = [a.text() for a in menu.actions() if not a.isSeparator()]
        assert labels == ["Do Thing"]

    def test_builtin_menus_never_removed(self, app):
        """Untagged (application-owned) menus must survive even when empty."""
        pm = _make_plugin_manager()
        im, bar = self._make_im_with_real_menubar(pm)
        bar.addMenu("File")  # built-in, untagged, currently empty
        mgr = PluginMenuManager(im)

        mgr.rebuild_plugin_menus()

        assert "File" in self._titles(bar)
