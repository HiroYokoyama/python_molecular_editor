"""Prevent an installed legacy plugin from loading a second native Qt binding."""

import sys

import pytest

from moleditpy.plugins.plugin_manager import PluginManager


@pytest.mark.parametrize(
    "statement",
    [
        "from PyQt6.QtWidgets import QWidget",
        "import PyQt6.QtCore",
        "__import__('PyQt6.QtGui')",
        "import importlib; importlib.import_module('PyQt6')",
    ],
)
def test_legacy_plugin_is_rejected_before_any_code_runs(tmp_path, statement):
    marker = tmp_path / "executed"
    plugin = tmp_path / "legacy.py"
    plugin.write_text(
        f"from pathlib import Path\nPath({str(marker)!r}).touch()\n{statement}\n",
        encoding="utf-8",
    )
    manager = PluginManager()
    manager._load_single_plugin(str(plugin), "legacy_qt_plugin", "")
    assert not marker.exists()
    assert "PyQt6" not in sys.modules
    assert manager.plugins[0]["module"] is None
    assert "migration to PySide6" in manager.plugins[0]["status"]


def test_package_submodules_are_checked_before_entry_point(tmp_path):
    plugin = tmp_path / "plugin"
    plugin.mkdir()
    entry = plugin / "__init__.py"
    entry.write_text("raise RuntimeError('must not execute')", encoding="utf-8")
    (plugin / "dialog.py").write_text("from PyQt6 import QtWidgets", encoding="utf-8")
    manager = PluginManager()
    manager._load_single_plugin(str(entry), "legacy_qt_package", "")
    assert "migration to PySide6" in manager.plugins[0]["status"]


def test_pyside_plugin_is_allowed_and_pyqt_documentation_is_not_an_import(tmp_path):
    plugin = tmp_path / "modern.py"
    plugin.write_text(
        '"""Migrated from PyQt6."""\nfrom PySide6.QtCore import QObject\nPLUGIN_NAME = "Modern"\n',
        encoding="utf-8",
    )
    manager = PluginManager()
    manager._load_single_plugin(str(plugin), "modern_qt_plugin", "")
    assert manager.plugins[0]["module"] is not None
