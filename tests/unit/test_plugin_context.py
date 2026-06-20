"""
tests/unit/test_plugin_context.py

Unit tests for PluginContext.mark_project_modified() — the V4 API
that notifies the application of unsaved changes via a stable context method
instead of direct mw.state_manager access.
"""

import sys
import os
import unittest
from unittest.mock import MagicMock, patch, PropertyMock

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
        mw = MagicMock()
        mw.state_manager.has_unsaved_changes = False
        ctx = _make_context(mw)
        ctx.mark_project_modified()
        self.assertTrue(mw.state_manager.has_unsaved_changes)

    def test_calls_update_window_title(self):
        mw = MagicMock()
        ctx = _make_context(mw)
        ctx.mark_project_modified()
        mw.state_manager.update_window_title.assert_called_once()

    def test_no_crash_when_mw_is_none(self):
        ctx = _make_context(mw=None)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(f"mark_project_modified() raised with mw=None: {exc}")

    def test_no_crash_when_state_manager_missing(self):
        mw = MagicMock(spec=[])  # no attributes
        ctx = _make_context(mw)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(f"mark_project_modified() raised without state_manager: {exc}")

    def test_no_crash_when_update_window_title_missing(self):
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
        mw = MagicMock()
        mw.state_manager.update_window_title.side_effect = RuntimeError("boom")
        ctx = _make_context(mw)
        try:
            ctx.mark_project_modified()
        except Exception as exc:
            self.fail(f"mark_project_modified() propagated internal exception: {exc}")

    def test_method_exists_on_plugincontext(self):
        self.assertTrue(
            hasattr(PluginContext, "mark_project_modified"),
            "PluginContext must expose mark_project_modified()",
        )
        self.assertTrue(callable(getattr(PluginContext, "mark_project_modified")))


class TestRefreshUi(unittest.TestCase):
    def test_calls_update_realtime_info(self):
        mw = MagicMock()
        _make_context(mw).refresh_ui()
        mw.state_manager.update_realtime_info.assert_called_once()

    def test_calls_update_undo_redo_actions(self):
        mw = MagicMock()
        _make_context(mw).refresh_ui()
        mw.edit_actions_manager.update_undo_redo_actions.assert_called_once()

    def test_calls_update_window_title(self):
        mw = MagicMock()
        _make_context(mw).refresh_ui()
        mw.state_manager.update_window_title.assert_called_once()

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).refresh_ui()

    def test_no_crash_when_managers_missing(self):
        _make_context(MagicMock(spec=[])).refresh_ui()


class TestFit3dView(unittest.TestCase):
    def test_calls_fit_to_view(self):
        mw = MagicMock()
        _make_context(mw).fit_3d_view()
        mw.view_3d_manager.fit_to_view.assert_called_once()

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).fit_3d_view()

    def test_no_crash_when_fit_to_view_missing(self):
        mw = MagicMock(spec=["view_3d_manager"])
        mw.view_3d_manager = MagicMock(spec=[])
        _make_context(mw).fit_3d_view()


class TestClearCanvas(unittest.TestCase):
    def test_calls_clear_2d_editor_default(self):
        mw = MagicMock()
        _make_context(mw).clear_canvas()
        mw.edit_actions_manager.clear_2d_editor.assert_called_once_with(
            push_to_undo=True
        )

    def test_calls_clear_2d_editor_no_undo(self):
        mw = MagicMock()
        _make_context(mw).clear_canvas(push_to_undo=False)
        mw.edit_actions_manager.clear_2d_editor.assert_called_once_with(
            push_to_undo=False
        )

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).clear_canvas()

    def test_no_crash_when_manager_missing(self):
        _make_context(MagicMock(spec=[])).clear_canvas()


class TestSet3dFeaturesEnabled(unittest.TestCase):
    def test_calls_enable_3d_features_true(self):
        mw = MagicMock()
        _make_context(mw).set_3d_features_enabled(True)
        mw.ui_manager.enable_3d_features.assert_called_once_with(True)

    def test_calls_enable_3d_features_false(self):
        mw = MagicMock()
        _make_context(mw).set_3d_features_enabled(False)
        mw.ui_manager.enable_3d_features.assert_called_once_with(False)

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).set_3d_features_enabled(True)

    def test_no_crash_when_ui_manager_missing(self):
        _make_context(MagicMock(spec=[])).set_3d_features_enabled(False)


class TestSetAnalysisEnabled(unittest.TestCase):
    def test_calls_set_enabled_true(self):
        mw = MagicMock()
        _make_context(mw).set_analysis_enabled(True)
        mw.init_manager.analysis_action.setEnabled.assert_called_once_with(True)

    def test_calls_set_enabled_false(self):
        mw = MagicMock()
        _make_context(mw).set_analysis_enabled(False)
        mw.init_manager.analysis_action.setEnabled.assert_called_once_with(False)

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).set_analysis_enabled(True)

    def test_no_crash_when_analysis_action_missing(self):
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=[])
        _make_context(mw).set_analysis_enabled(True)


class TestCheckChemistryProblems(unittest.TestCase):
    def test_calls_fallback(self):
        mw = MagicMock()
        _make_context(mw).check_chemistry_problems()
        mw.compute_manager.check_chemistry_problems_fallback.assert_called_once()

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).check_chemistry_problems()

    def test_no_crash_when_compute_manager_missing(self):
        _make_context(MagicMock(spec=[])).check_chemistry_problems()


class TestRefresh2dScene(unittest.TestCase):
    def test_calls_update_all_items(self):
        mw = MagicMock()
        _make_context(mw).refresh_2d_scene()
        mw.init_manager.scene.update_all_items.assert_called_once()

    def test_no_crash_when_mw_is_none(self):
        _make_context(None).refresh_2d_scene()

    def test_no_crash_when_init_manager_missing(self):
        _make_context(MagicMock(spec=[])).refresh_2d_scene()

    def test_no_crash_when_scene_missing(self):
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=[])
        _make_context(mw).refresh_2d_scene()

    def test_no_crash_when_update_all_items_missing(self):
        mw = MagicMock(spec=["init_manager"])
        mw.init_manager = MagicMock(spec=["scene"])
        mw.init_manager.scene = MagicMock(spec=[])
        _make_context(mw).refresh_2d_scene()


if __name__ == "__main__":
    unittest.main()
