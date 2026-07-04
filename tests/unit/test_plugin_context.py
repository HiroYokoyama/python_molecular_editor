"""
tests/unit/test_plugin_context.py

Unit tests for PluginContext.mark_project_modified() — the V4 API
that notifies the application of unsaved changes via a stable context method
instead of direct mw.state_manager access.
"""

import sys
import os
import unittest
import pytest
from unittest.mock import MagicMock

# Ensure moleditpy is importable
_SRC = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from moleditpy.plugins.plugin_interface import PluginContext


def _make_context(mw=None):
    """Return a PluginContext whose manager.get_main_window() returns mw."""
    manager = MagicMock()
    manager.get_main_window.return_value = mw
    return PluginContext(manager, "TestPlugin")


class TestMarkProjectModified(unittest.TestCase):
    def test_sets_has_unsaved_changes(self):
        """mark_project_modified sets has_unsaved_changes to True on the state manager."""
        mw = MagicMock()
        mw.state_manager.has_unsaved_changes = False
        ctx = _make_context(mw)
        ctx.mark_project_modified()
        self.assertTrue(mw.state_manager.has_unsaved_changes)

    def test_calls_update_window_title(self):
        """mark_project_modified calls update_window_title on the state manager."""
        mw = MagicMock()
        ctx = _make_context(mw)
        ctx.mark_project_modified()
        mw.state_manager.update_window_title.assert_called_once()

    def test_no_crash_when_state_manager_missing(self):
        """mark_project_modified does not raise when the main window has no state_manager."""
        mw = MagicMock(spec=[])  # no attributes
        ctx = _make_context(mw)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(f"mark_project_modified() raised without state_manager: {exc}")

    def test_no_crash_when_update_window_title_missing(self):
        """mark_project_modified does not raise when update_window_title is absent."""
        mw = MagicMock()
        del mw.state_manager.update_window_title
        ctx = _make_context(mw)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(
                f"mark_project_modified() raised without update_window_title: {exc}"
            )

    def test_no_crash_when_state_manager_raises(self):
        """mark_project_modified does not propagate exceptions from state_manager."""
        mw = MagicMock()
        mw.state_manager.update_window_title.side_effect = RuntimeError("boom")
        ctx = _make_context(mw)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(f"mark_project_modified() propagated internal exception: {exc}")

    def test_method_exists_on_plugincontext(self):
        """PluginContext exposes a callable mark_project_modified method."""
        self.assertTrue(
            hasattr(PluginContext, "mark_project_modified"),
            "PluginContext must expose mark_project_modified()",
        )
        self.assertTrue(callable(getattr(PluginContext, "mark_project_modified")))


class TestRefreshUi(unittest.TestCase):
    def test_calls_all_required_managers(self):
        """refresh_ui calls update_realtime_info, update_undo_redo_actions, and update_window_title."""
        mw = MagicMock()
        _make_context(mw).refresh_ui()
        mw.state_manager.update_realtime_info.assert_called_once()
        mw.edit_actions_manager.update_undo_redo_actions.assert_called_once()
        mw.state_manager.update_window_title.assert_called_once()

    def test_no_crash_when_managers_missing(self):
        """refresh_ui does not raise when the main window has no manager attributes."""
        _make_context(MagicMock(spec=[])).refresh_ui()


class TestFit2dView(unittest.TestCase):
    def test_calls_fit_to_view(self):
        """fit_2d_view delegates to view_3d_manager.fit_to_view."""
        mw = MagicMock()
        _make_context(mw).fit_2d_view()
        mw.view_3d_manager.fit_to_view.assert_called_once()

    def test_no_crash_when_fit_to_view_missing(self):
        """fit_2d_view does not raise when fit_to_view is absent on the view manager."""
        mw = MagicMock(spec=["view_3d_manager"])
        mw.view_3d_manager = MagicMock(spec=[])
        _make_context(mw).fit_2d_view()


@pytest.mark.parametrize("push_to_undo", [True, False])
def test_clear_canvas_delegates(push_to_undo):
    """clear_canvas delegates to edit_actions_manager.clear_2d_editor with push_to_undo."""
    mw = MagicMock()
    _make_context(mw).clear_canvas(push_to_undo=push_to_undo)
    mw.edit_actions_manager.clear_2d_editor.assert_called_once_with(
        push_to_undo=push_to_undo
    )


def test_clear_canvas_no_crash_when_manager_missing():
    """clear_canvas does not raise when edit_actions_manager is absent."""
    _make_context(MagicMock(spec=[])).clear_canvas()


@pytest.mark.parametrize("enabled", [True, False])
def test_set_3d_features_enabled_delegates(enabled):
    """set_3d_features_enabled delegates to ui_manager.enable_3d_features."""
    mw = MagicMock()
    _make_context(mw).set_3d_features_enabled(enabled)
    mw.ui_manager.enable_3d_features.assert_called_once_with(enabled)


def test_set_3d_features_enabled_no_crash_when_ui_manager_missing():
    """set_3d_features_enabled does not raise when ui_manager is absent."""
    _make_context(MagicMock(spec=[])).set_3d_features_enabled(False)


@pytest.mark.parametrize("enabled", [True, False])
def test_set_analysis_enabled_delegates(enabled):
    """set_analysis_enabled delegates to init_manager.analysis_action.setEnabled."""
    mw = MagicMock()
    _make_context(mw).set_analysis_enabled(enabled)
    mw.init_manager.analysis_action.setEnabled.assert_called_once_with(enabled)


def test_set_analysis_enabled_no_crash_when_action_missing():
    """set_analysis_enabled does not raise when analysis_action is absent."""
    mw = MagicMock(spec=["init_manager"])
    mw.init_manager = MagicMock(spec=[])
    _make_context(mw).set_analysis_enabled(True)


class TestCheckChemistryProblems(unittest.TestCase):
    def test_calls_fallback(self):
        """check_chemistry_problems delegates to compute_manager.check_chemistry_problems_fallback."""
        mw = MagicMock()
        _make_context(mw).check_chemistry_problems()
        mw.compute_manager.check_chemistry_problems_fallback.assert_called_once()

    def test_no_crash_when_compute_manager_missing(self):
        """check_chemistry_problems does not raise when compute_manager is absent."""
        _make_context(MagicMock(spec=[])).check_chemistry_problems()


class TestRefresh2dScene(unittest.TestCase):
    def test_calls_update_all_items(self):
        """refresh_2d_scene delegates to init_manager.scene.update_all_items."""
        mw = MagicMock()
        _make_context(mw).refresh_2d_scene()
        mw.init_manager.scene.update_all_items.assert_called_once()

    def test_no_crash_when_init_manager_missing(self):
        """refresh_2d_scene does not raise when init_manager is absent."""
        _make_context(MagicMock(spec=[])).refresh_2d_scene()

    def test_no_crash_when_scene_missing(self):
        """refresh_2d_scene does not raise when scene is absent on init_manager."""
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=[])
        _make_context(mw).refresh_2d_scene()

    def test_no_crash_when_update_all_items_missing(self):
        """refresh_2d_scene does not raise when update_all_items is absent on scene."""
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=["scene"])
        mw.init_manager.scene = MagicMock(spec=[])
        _make_context(mw).refresh_2d_scene()


@pytest.mark.parametrize(
    "call",
    [
        lambda ctx: ctx.mark_project_modified(),
        lambda ctx: ctx.refresh_ui(),
        lambda ctx: ctx.fit_2d_view(),
        lambda ctx: ctx.clear_canvas(),
        lambda ctx: ctx.set_3d_features_enabled(True),
        lambda ctx: ctx.set_analysis_enabled(True),
        lambda ctx: ctx.check_chemistry_problems(),
        lambda ctx: ctx.refresh_2d_scene(),
    ],
    ids=[
        "mark_project_modified",
        "refresh_ui",
        "fit_2d_view",
        "clear_canvas",
        "set_3d_features_enabled",
        "set_analysis_enabled",
        "check_chemistry_problems",
        "refresh_2d_scene",
    ],
)
def test_no_crash_when_mw_is_none(call):
    """All PluginContext API methods are no-ops when the main window is None."""
    call(_make_context(None))


if __name__ == "__main__":
    unittest.main()
