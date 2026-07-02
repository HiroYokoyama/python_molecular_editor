# -*- coding: utf-8 -*-
"""E2E project save/load round-trip tests: real MainWindow, real IOManager.

Covers the .pmeprj (JSON) and .pmeraw (pickle) project formats end-to-end:
draw -> save -> clear -> load -> verify nothing is lost, plus the
unsaved-changes prompt and Save-As extension handling.
"""

import importlib.util
import json

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QMessageBox

_PKG = (
    "moleditpy_linux"
    if importlib.util.find_spec("moleditpy_linux") is not None
    else "moleditpy"
)


def _draw_ethanol_like(window):
    """Draw a small 3-atom test molecule (C-C-O with a double bond)."""
    scene = window.init_manager.scene
    c1 = scene.create_atom("C", QPointF(0, 0))
    c2 = scene.create_atom("C", QPointF(50, 0))
    o1 = scene.create_atom("O", QPointF(100, 0))
    scene.create_bond(scene.atom_items[c1], scene.atom_items[c2], bond_order=1)
    scene.create_bond(scene.atom_items[c2], scene.atom_items[o1], bond_order=2)
    return c1, c2, o1


def _snapshot(window):
    atoms = {
        aid: (d["symbol"], d.get("charge", 0), tuple(d["pos"]))
        for aid, d in window.data.atoms.items()
    }
    bonds = {
        tuple(sorted(k)): (v["order"], v.get("stereo", 0))
        for k, v in window.data.bonds.items()
    }
    return atoms, bonds, window.data.next_atom_id


def _patch_unsaved_prompt(monkeypatch, button):
    monkeypatch.setattr(
        f"{_PKG}.ui.app_state.QMessageBox.question",
        lambda *a, **k: button,
        raising=True,
    )


@pytest.mark.gui
def test_pmeprj_roundtrip(window, qtbot, tmp_path, monkeypatch):
    """Draw -> serialize -> clear -> load .pmeprj -> identical document."""
    _draw_ethanol_like(window)
    atoms_before, bonds_before, next_id_before = _snapshot(window)

    json_data = window.state_manager.create_json_data()
    path = tmp_path / "roundtrip.pmeprj"
    path.write_text(json.dumps(json_data), encoding="utf-8")

    window.edit_actions_manager.clear_all(skip_check=True)
    assert window.data.atoms == {}

    _patch_unsaved_prompt(monkeypatch, QMessageBox.StandardButton.No)
    window.io_manager.open_project_file(str(path))
    qtbot.wait(200)  # let deferred QTimer plotter calls fire (None-guarded)

    atoms_after, bonds_after, next_id_after = _snapshot(window)
    assert set(atoms_after) == set(atoms_before)
    for aid, (sym, chg, pos) in atoms_before.items():
        sym2, chg2, pos2 = atoms_after[aid]
        assert (sym2, chg2) == (sym, chg)
        assert pos2 == pytest.approx(pos)
    assert bonds_after == bonds_before
    assert next_id_after == next_id_before
    assert window.state_manager.has_unsaved_changes is False


@pytest.mark.gui
def test_pmeprj_save_via_dialog_appends_extension(window, tmp_path, monkeypatch):
    """Save As with a bare filename creates <name>.pmeprj and clears dirty flag."""
    _draw_ethanol_like(window)
    target = tmp_path / "noext"
    monkeypatch.setattr(
        f"{_PKG}.ui.io_logic.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(target), None),
        raising=True,
    )

    window.io_manager.save_project_as()

    saved = tmp_path / "noext.pmeprj"
    assert saved.exists()
    data = json.loads(saved.read_text(encoding="utf-8"))
    assert data["format"] == "PME Project"
    assert len(data["2d_structure"]["atoms"]) == 3
    assert window.state_manager.has_unsaved_changes is False


@pytest.mark.gui
def test_pmeraw_roundtrip(window, qtbot, tmp_path, monkeypatch):
    """Draw -> Save as .pmeraw (pickle) via dialog -> clear -> load -> identical."""
    _draw_ethanol_like(window)
    atoms_before, bonds_before, next_id_before = _snapshot(window)

    target = tmp_path / "raw_project"
    monkeypatch.setattr(
        f"{_PKG}.ui.io_logic.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(target), None),
        raising=True,
    )
    window.io_manager.save_raw_data()
    saved = tmp_path / "raw_project.pmeraw"
    assert saved.exists()

    window.edit_actions_manager.clear_all(skip_check=True)
    assert window.data.atoms == {}

    _patch_unsaved_prompt(monkeypatch, QMessageBox.StandardButton.No)
    window.io_manager.load_raw_data(str(saved))
    qtbot.wait(200)

    atoms_after, bonds_after, next_id_after = _snapshot(window)
    assert set(atoms_after) == set(atoms_before)
    assert bonds_after == bonds_before
    assert next_id_after == next_id_before


@pytest.mark.gui
def test_load_cancelled_on_unsaved_changes_keeps_document(
    window, qtbot, tmp_path, monkeypatch
):
    """Cancelling the unsaved-changes prompt aborts the load, keeping the doc."""
    # Prepare a valid project file first (1-atom document)
    scene = window.init_manager.scene
    scene.create_atom("N", QPointF(0, 0))
    path = tmp_path / "other.pmeprj"
    path.write_text(
        json.dumps(window.state_manager.create_json_data()), encoding="utf-8"
    )

    # Now build a different unsaved document
    window.edit_actions_manager.clear_all(skip_check=True)
    _draw_ethanol_like(window)
    window.state_manager.has_unsaved_changes = True
    atoms_before, bonds_before, _ = _snapshot(window)

    _patch_unsaved_prompt(monkeypatch, QMessageBox.StandardButton.Cancel)
    window.io_manager.open_project_file(str(path))
    qtbot.wait(100)

    atoms_after, bonds_after, _ = _snapshot(window)
    assert atoms_after == atoms_before
    assert bonds_after == bonds_before
