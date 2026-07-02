# -*- coding: utf-8 -*-
"""E2E plugin workflow tests using the REAL PluginManager.

Covers install -> discover -> initialize for .py and .zip plugins against a
temporary plugin directory, and plugin state persistence through the project
save/load cycle on a live MainWindow.
"""

import importlib
import importlib.util
import textwrap
import zipfile

import pytest

_PKG = (
    "moleditpy_linux"
    if importlib.util.find_spec("moleditpy_linux") is not None
    else "moleditpy"
)
PluginManager = importlib.import_module(f"{_PKG}.plugins.plugin_manager").PluginManager

PLUGIN_SRC = textwrap.dedent(
    """\
    PLUGIN_NAME = "E2E Test Plugin"
    PLUGIN_VERSION = "1.0"

    def initialize(context):
        context.add_menu_action("E2E/Test", lambda: None, text="E2E Action")
    """
)


def _make_pm(tmp_path):
    pm = PluginManager(main_window=None)
    pm.plugin_dir = str(tmp_path / "plugins")
    return pm


def test_install_and_discover_py_plugin(tmp_path):
    """A .py plugin installs, loads, and registers its menu action."""
    src = tmp_path / "e2e_plugin.py"
    src.write_text(PLUGIN_SRC, encoding="utf-8")

    pm = _make_pm(tmp_path)
    ok, msg = pm.install_plugin(str(src))
    assert ok, msg

    plugins = pm.discover_plugins()
    match = [p for p in plugins if p["name"] == "E2E Test Plugin"]
    assert len(match) == 1
    assert match[0]["status"] == "Loaded"
    assert match[0]["version"] == "1.0"
    assert any(a["path"] == "E2E/Test" for a in pm.menu_actions)


def test_zip_single_file_install_gets_wrapper_folder(tmp_path):
    """A ZIP containing only one .py file extracts into a wrapper folder."""
    zip_path = tmp_path / "soloplugin.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("solo.py", PLUGIN_SRC)

    pm = _make_pm(tmp_path)
    ok, msg = pm.install_plugin(str(zip_path))
    assert ok, msg

    extracted = tmp_path / "plugins" / "soloplugin" / "solo.py"
    assert extracted.exists()
    plugins = pm.discover_plugins()
    assert any(p["name"] == "E2E Test Plugin" for p in plugins)


@pytest.mark.gui
def test_plugin_state_persists_through_project_save_load(window):
    """Plugin save/load handlers round-trip data through create/load_json_data."""
    pm = PluginManager(main_window=window)
    window.plugin_manager = pm

    received = {}
    pm.register_save_handler("StatefulPlugin", lambda: {"counter": 42})
    pm.register_load_handler("StatefulPlugin", received.update)

    json_data = window.state_manager.create_json_data()
    assert json_data["plugins"]["StatefulPlugin"] == {"counter": 42}

    window.state_manager.load_from_json_data(json_data)
    assert received == {"counter": 42}


@pytest.mark.gui
def test_missing_plugin_state_preserved_on_resave(window):
    """State from an absent plugin survives a load + save cycle unchanged."""
    pm = PluginManager(main_window=window)
    window.plugin_manager = pm

    json_data = window.state_manager.create_json_data()
    json_data["plugins"] = {"AbsentPlugin": {"keep": "me"}}

    window.state_manager.load_from_json_data(json_data)  # no handler registered
    resaved = window.state_manager.create_json_data()
    assert resaved["plugins"]["AbsentPlugin"] == {"keep": "me"}
