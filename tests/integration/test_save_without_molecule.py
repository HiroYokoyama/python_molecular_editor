"""A project with no molecule still saves, so plugin-only state is kept."""

import json
from unittest.mock import MagicMock

import pytest


@pytest.fixture
def empty_window_with_plugin(window):
    window.edit_actions_manager.clear_all(skip_check=True)
    assert not window.state_manager.data.atoms
    assert window.view_3d_manager.current_mol is None
    loaded = []
    window.plugin_manager.save_handlers["demo"] = lambda: {"notes": "plugin only"}
    window.plugin_manager.load_handlers["demo"] = loaded.append
    return window, loaded


@pytest.mark.parametrize("action", ["save_project", "save_project_as", "save_as_json"])
def test_an_empty_document_saves_and_reloads_its_plugin_data(
    empty_window_with_plugin, monkeypatch, tmp_path, action
):
    window, loaded = empty_window_with_plugin
    path = tmp_path / "plugin_only.pmeprj"
    monkeypatch.setattr(
        "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
        MagicMock(return_value=(str(path), "")),
    )
    window.statusBar().clearMessage()

    getattr(window.io_manager, action)()

    assert window.statusBar().currentMessage() != "Error: Nothing to save."
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["format"] == "PME Project"
    assert data["plugins"] == {"demo": {"notes": "plugin only"}}
    assert "2d_structure" not in data
    assert data.get("3d_structure") is None

    window.io_manager.load_json_data(str(path))
    assert loaded == [{"notes": "plugin only"}]
    assert not window.state_manager.data.atoms


def test_saving_over_the_open_project_needs_no_molecule(
    empty_window_with_plugin, tmp_path
):
    window, _ = empty_window_with_plugin
    path = tmp_path / "open.pmeprj"
    window.set_current_file_path(str(path))

    window.io_manager.save_project()

    assert json.loads(path.read_text(encoding="utf-8"))["plugins"] == {
        "demo": {"notes": "plugin only"}
    }
