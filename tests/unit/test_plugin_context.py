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


if __name__ == "__main__":
    unittest.main()
