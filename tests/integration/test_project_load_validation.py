"""Malformed project contents must not destroy the open document."""

import copy
import json
import pickle
from unittest.mock import MagicMock

import pytest
from PySide6.QtCore import QPointF


@pytest.mark.parametrize(
    "structure",
    [
        None,
        [],
        {"atoms": None},
        {"atoms": [None]},
        {"atoms": [{"id": 0, "symbol": "C", "x": "bad", "y": 0}]},
        {"atoms": [{"id": 0, "symbol": "C", "x": float("nan"), "y": 0}]},
        {"atoms": [{"id": 0, "symbol": "C", "x": 10**400, "y": 0}]},
        {"atoms": [{"id": [], "symbol": "C", "x": 0, "y": 0}]},
        {"atoms": [{"id": 0, "symbol": [], "x": 0, "y": 0}]},
        {"atoms": [{"id": 0, "symbol": "C", "x": 0, "y": 0, "charge": []}]},
        {"atoms": [{"id": 0, "symbol": "C", "x": 0, "y": 0, "radical": -1}]},
        {"atoms": [{"id": 0, "symbol": "C", "x": 0, "y": 0, "rs": "bad"}]},
        {"atoms": [], "bonds": None},
        {"atoms": [], "bonds": [{"atom1": 0, "atom2": 1, "order": 1}]},
        {"atoms": [], "next_atom_id": "bad"},
    ],
)
def test_invalid_2d_project_keeps_open_document(
    window, monkeypatch, tmp_path, structure
):
    _assert_failed_load_keeps_document(
        window,
        monkeypatch,
        tmp_path,
        {"format": "PME Project", "2d_structure": structure},
    )


@pytest.mark.parametrize(
    "structure", [[], "bad", {"mol_binary_base64": "not base64"}, {"atoms": None}]
)
def test_invalid_3d_project_keeps_open_document(
    window, monkeypatch, tmp_path, structure
):
    _assert_failed_load_keeps_document(
        window,
        monkeypatch,
        tmp_path,
        {"format": "PME Project", "3d_structure": structure},
    )


@pytest.mark.parametrize("state", [[], None, {"atoms": None}, {"atoms": {0: None}}])
def test_invalid_raw_project_keeps_open_document(window, monkeypatch, tmp_path, state):
    _assert_failed_load_keeps_document(window, monkeypatch, tmp_path, state, raw=True)


def _assert_failed_load_keeps_document(
    window, monkeypatch, tmp_path, payload, raw=False
):
    scene = window.init_manager.scene
    scene.create_atom("N", QPointF(10, 20), charge=1)
    window.edit_actions_manager.push_undo_state()
    original_path = str(tmp_path / "original.pmeprj")
    window.set_current_file_path(original_path)
    window.set_has_unsaved_changes(True)
    previous = copy.deepcopy(window.state_manager.get_current_state())
    history = copy.deepcopy(window.edit_actions_manager.undo_stack)
    monkeypatch.setattr(window.state_manager, "check_unsaved_changes", lambda: True)
    reset = MagicMock()
    monkeypatch.setattr(window.plugin_manager, "invoke_document_reset_handlers", reset)
    report = MagicMock()
    monkeypatch.setattr(window.io_manager, "_report_load_error", report)
    path = tmp_path / ("broken.pmeraw" if raw else "broken.pmeprj")
    if raw:
        path.write_bytes(pickle.dumps(payload))
        monkeypatch.setattr(window.io_manager, "_confirm_pickle_load", lambda _: True)
        window.io_manager.load_raw_data(str(path))
    else:
        path.write_text(json.dumps(payload), encoding="utf-8")
        window.io_manager.load_json_data(str(path))
    assert window.state_manager.get_current_state() == previous
    assert window.edit_actions_manager.undo_stack == history
    assert window.init_manager.current_file_path == original_path
    assert window.state_manager.has_unsaved_changes
    assert len(scene.atom_items) == 1
    assert next(iter(scene.atom_items.values())).symbol == "N"
    reset.assert_not_called()
    report.assert_called_once()
