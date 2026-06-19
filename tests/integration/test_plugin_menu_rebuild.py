"""Integration tests for plugin menu reload completeness.

Verifies that PluginMenuManager.rebuild_plugin_menus fully rebuilds all
six integration points: menus, toolbar, export, file-openers, analysis, and
3D styles — and that MainInitManager correctly delegates to it.
"""
import pytest
from unittest.mock import MagicMock, patch, call

from moleditpy.ui.plugin_menu_manager import PluginMenuManager


# QAction(text, MagicMock) fails because PyQt6 strictly validates parent type.
# Patch it to drop the mock parent so tests run headlessly (including CI).
@pytest.fixture(autouse=True)
def _patch_qaction(monkeypatch):
    from PyQt6.QtGui import QAction as _RealQAction

    def _safe_qaction(text: str, parent=None) -> _RealQAction:
        return _RealQAction(text)

    monkeypatch.setattr("moleditpy.ui.plugin_menu_manager.QAction", _safe_qaction)


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
    return pm


def _make_im(plugin_manager: MagicMock | None = None) -> MagicMock:
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
        pmm._integrate_plugin_export_actions = track("export_actions")
        pmm._integrate_plugin_file_openers = track("file_openers")
        pmm._integrate_plugin_analysis_tools = track("analysis_tools")
        pmm._update_style_menu_with_plugins = track("style_menu")

        pmm.rebuild_plugin_menus()

        assert called == [
            "menu_actions",
            "toolbar_actions",
            "export_actions",
            "file_openers",
            "analysis_tools",
            "style_menu",
        ], f"Expected all six steps in order, got: {called}"

    def test_all_six_steps_run_when_all_registries_populated(self):
        """With all registries populated, all six steps still run."""
        pm = _make_plugin_manager(
            menu_actions=[{"path": "P/A", "callback": MagicMock(), "text": "A", "shortcut": None}],
            toolbar_actions=[{"text": "TB", "callback": MagicMock(), "icon": "", "tooltip": ""}],
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
        pmm._integrate_plugin_export_actions = track("export_actions")
        pmm._integrate_plugin_file_openers = track("file_openers")
        pmm._integrate_plugin_analysis_tools = track("analysis_tools")
        pmm._update_style_menu_with_plugins = track("style_menu")

        pmm.rebuild_plugin_menus()

        assert set(called) == {
            "menu_actions", "toolbar_actions", "export_actions",
            "file_openers", "analysis_tools", "style_menu",
        }

    def test_separator_flag_reset_before_rebuild(self):
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

        pmm.add_registered_plugin_actions = lambda: (_ for _ in ()).throw(RuntimeError("boom"))
        pmm.add_plugin_toolbar_actions = step_ok("toolbar")
        pmm._integrate_plugin_export_actions = step_ok("export")
        pmm._integrate_plugin_file_openers = step_ok("file_openers")
        pmm._integrate_plugin_analysis_tools = step_ok("analysis")
        pmm._update_style_menu_with_plugins = step_ok("style")

        pmm.rebuild_plugin_menus()

        assert "toolbar" in survived
        assert "export" in survived
        assert "style" in survived


# ---------------------------------------------------------------------------
# MainInitManager delegates to PluginMenuManager
# ---------------------------------------------------------------------------

class TestMainInitManagerDelegation:
    """Verify MainInitManager.rebuild_plugin_menus and friends delegate to PluginMenuManager."""

    def _make_real_init_manager_stub(self):
        """Return a MainInitManager with plugin_menu_manager injected."""
        from moleditpy.ui.main_window_init import MainInitManager
        # Construct without running __init__ to keep it lightweight
        im = object.__new__(MainInitManager)
        im.plugin_menu_manager = MagicMock()
        return im

    def test_rebuild_plugin_menus_delegates(self):
        from moleditpy.ui.main_window_init import MainInitManager
        im = object.__new__(MainInitManager)
        pmm_mock = MagicMock()
        im.plugin_menu_manager = pmm_mock

        im.rebuild_plugin_menus()

        pmm_mock.rebuild_plugin_menus.assert_called_once()

    def test_update_plugin_menu_delegates(self):
        from moleditpy.ui.main_window_init import MainInitManager
        im = object.__new__(MainInitManager)
        pmm_mock = MagicMock()
        im.plugin_menu_manager = pmm_mock

        sentinel_menu = MagicMock()
        im.update_plugin_menu(sentinel_menu)

        pmm_mock.update_plugin_menu.assert_called_once_with(sentinel_menu)

    def test_add_registered_plugin_actions_delegates(self):
        from moleditpy.ui.main_window_init import MainInitManager
        im = object.__new__(MainInitManager)
        pmm_mock = MagicMock()
        im.plugin_menu_manager = pmm_mock

        im.add_registered_plugin_actions()

        pmm_mock.add_registered_plugin_actions.assert_called_once()

    def test_add_plugin_toolbar_actions_delegates(self):
        from moleditpy.ui.main_window_init import MainInitManager
        im = object.__new__(MainInitManager)
        pmm_mock = MagicMock()
        im.plugin_menu_manager = pmm_mock

        im.add_plugin_toolbar_actions()

        pmm_mock.add_plugin_toolbar_actions.assert_called_once()


# ---------------------------------------------------------------------------
# update_plugin_menu — discovers and integrates
# ---------------------------------------------------------------------------

class TestUpdatePluginMenuIntegration:
    def test_discover_plugins_called_with_host(self):
        pm = _make_plugin_manager()
        pm.discover_plugins.return_value = []
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        plugin_menu = MagicMock()
        pmm.update_plugin_menu(plugin_menu)

        pm.discover_plugins.assert_called_once_with(im.host)

    def test_calls_seven_integration_methods(self):
        pm = _make_plugin_manager()
        pm.discover_plugins.return_value = []
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        integration_methods = [
            "_update_style_menu_with_plugins",
            "add_registered_plugin_actions",
            "add_plugin_toolbar_actions",
            "_add_legacy_plugin_actions",
            "_integrate_plugin_export_actions",
            "_integrate_plugin_file_openers",
            "_integrate_plugin_analysis_tools",
        ]
        mocks = {m: MagicMock() for m in integration_methods}
        for m, mock in mocks.items():
            setattr(pmm, m, mock)

        pmm.update_plugin_menu(MagicMock())

        for m, mock in mocks.items():
            mock.assert_called_once(), f"Expected {m} to be called once"

    def test_does_nothing_when_plugin_manager_is_none(self):
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
        callback = MagicMock()
        pm = _make_plugin_manager(export_actions=[{"label": "Export XYZ", "callback": callback}])
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        pmm._integrate_plugin_export_actions()

        export_menu = im.export_button.menu.return_value
        # A separator and at least one action were added
        export_menu.addSeparator.assert_called_once()
        export_menu.addAction.assert_called_once()
        added_action = export_menu.addAction.call_args[0][0]
        # The added action has the right label
        assert added_action.text() == "Export XYZ"

    def test_multiple_export_actions_all_added(self):
        pm = _make_plugin_manager(export_actions=[
            {"label": "Export XYZ", "callback": MagicMock()},
            {"label": "Export CIF", "callback": MagicMock()},
        ])
        im = _make_im(pm)
        pmm = PluginMenuManager(im)

        pmm._integrate_plugin_export_actions()

        assert im.export_button.menu.return_value.addAction.call_count == 2


# ---------------------------------------------------------------------------
# Analysis tool propagation end-to-end
# ---------------------------------------------------------------------------

class TestAnalysisToolEndToEnd:
    def test_analysis_tools_appear_in_analysis_menu(self):
        callback = MagicMock()
        pm = _make_plugin_manager(
            analysis_tools=[{"label": "NMR Predict", "callback": callback, "plugin": "NMRPro"}]
        )
        im = _make_im(pm)

        # Wire up an Analysis menu in the menubar
        analysis_menu = MagicMock()
        analysis_action = MagicMock()
        analysis_action.text.return_value = "Analysis"
        analysis_action.menu.return_value = analysis_menu
        im.host.menuBar.return_value.actions.return_value = [analysis_action]

        pmm = PluginMenuManager(im)
        pmm._integrate_plugin_analysis_tools()

        analysis_menu.addSeparator.assert_called_once()
        assert analysis_menu.addAction.call_count == 1
        added = analysis_menu.addAction.call_args[0][0]
        assert "NMR Predict" in added.text()
        assert "NMRPro" in added.text()
