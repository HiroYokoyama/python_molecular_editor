# -*- coding: utf-8 -*-
"""Save the benzene session to a .pmeprj project and load it back.

Both halves go through the application's own menu handlers -- `save_project_as`
and `open_project_file` -- with only the file dialogs stubbed, so the test
covers the real serialisation, the real "clear then rebuild" load path, and the
real scene reconstruction.
"""

import json

import numpy as np
import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QFileDialog

from benzene import assert_is_2d_benzene, draw_benzene_by_hand, place_benzene_template
from test_benzene_conversion import (
    CONVERT_TIMEOUT_MS,
    _convert_via_button,
    assert_is_3d_benzene,
)

pytestmark = pytest.mark.full_gui


def _save_project_as(window, monkeypatch, path):
    """Drive File > Save Project As with the save dialog answering `path`."""
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: (str(path), "PME Project Files (*.pmeprj)")),
    )
    window.io_manager.save_project_as()
    assert path.exists(), f"save_project_as wrote nothing to {path}"
    return path


def _positions(mol):
    return np.array(mol.GetConformer().GetPositions())


# ---------------------------------------------------------------------------
# The file itself
# ---------------------------------------------------------------------------


def test_saved_project_is_a_pme_project_json(full_window, qtbot, monkeypatch, tmp_path):
    place_benzene_template(full_window, qtbot)
    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")

    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload.get("format") == "PME Project", (
        f"unexpected format tag: {payload.get('format')!r}"
    )
    assert payload.get("2d_structure", {}).get("atoms"), "project stores no 2D atoms"
    assert payload.get("application") == "MoleditPy"


def test_save_clears_the_unsaved_marker(full_window, qtbot, monkeypatch, tmp_path):
    place_benzene_template(full_window, qtbot)
    assert full_window.state_manager.has_unsaved_changes

    _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")
    assert not full_window.state_manager.has_unsaved_changes


# ---------------------------------------------------------------------------
# Round trips
# ---------------------------------------------------------------------------


def test_2d_benzene_survives_a_round_trip(full_window, qtbot, monkeypatch, tmp_path):
    """Save a drawn ring, wipe the canvas, load it back unchanged."""
    place_benzene_template(full_window, qtbot)
    before = {
        aid: (a["symbol"], tuple(round(c, 6) for c in a["pos"]))
        for aid, a in full_window.state_manager.data.atoms.items()
    }
    before_bonds = {
        k: v["order"] for k, v in full_window.state_manager.data.bonds.items()
    }

    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False
    assert len(full_window.state_manager.data.atoms) == 0

    full_window.io_manager.open_project_file(str(target))
    qtbot.waitUntil(
        lambda: len(full_window.state_manager.data.atoms) == 6, timeout=30000
    )

    data = full_window.state_manager.data
    assert_is_2d_benzene(data)

    after = {
        aid: (a["symbol"], tuple(round(c, 6) for c in a["pos"]))
        for aid, a in data.atoms.items()
    }
    assert after == before, "2D atom identities/positions changed across the round trip"
    assert {k: v["order"] for k, v in data.bonds.items()} == before_bonds, (
        "bond orders changed across the round trip"
    )


def test_round_trip_rebuilds_the_scene_items(full_window, qtbot, monkeypatch, tmp_path):
    """Loading restores drawable items, not just the data model."""
    place_benzene_template(full_window, qtbot)
    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False

    full_window.io_manager.open_project_file(str(target))
    scene = full_window.init_manager.scene
    qtbot.waitUntil(lambda: len(scene.atom_items) == 6, timeout=30000)
    assert len(scene.bond_items) == 6
    assert all(item.scene() is scene for item in scene.atom_items.values())


def test_3d_geometry_survives_a_round_trip(full_window, qtbot, monkeypatch, tmp_path):
    """The converted 3D structure comes back with identical coordinates."""
    draw_benzene_by_hand(full_window)
    mol = _convert_via_button(full_window, qtbot)
    before = _positions(mol)

    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene3d.pmeprj")
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload.get("3d_structure"), "project kept no 3D structure"

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False
    assert full_window.view_3d_manager.current_mol is None

    full_window.io_manager.open_project_file(str(target))
    qtbot.waitUntil(
        lambda: full_window.view_3d_manager.current_mol is not None,
        timeout=CONVERT_TIMEOUT_MS,
    )

    restored = full_window.view_3d_manager.current_mol
    assert_is_3d_benzene(restored, "pmeprj-round-trip")

    after = _positions(restored)
    assert after.shape == before.shape, (
        f"atom count changed: {before.shape} -> {after.shape}"
    )
    # MOL blocks carry four decimals, so this is exact within the format.
    assert np.allclose(after, before, atol=1e-3), (
        f"3D coordinates drifted; max delta {np.abs(after - before).max():.5f} A"
    )


def test_charges_and_radicals_survive_a_round_trip(
    full_window, qtbot, monkeypatch, tmp_path
):
    """Per-atom charge and radical state is part of the project format."""
    ids = draw_benzene_by_hand(full_window)
    scene = full_window.init_manager.scene

    charged, radical = scene.atom_items[ids[0]], scene.atom_items[ids[3]]
    charged.charge = 1
    full_window.state_manager.data.atoms[ids[0]]["charge"] = 1
    radical.radical = 1
    full_window.state_manager.data.atoms[ids[3]]["radical"] = 1
    scene.update_all_items()

    target = _save_project_as(full_window, monkeypatch, tmp_path / "ions.pmeprj")

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False

    full_window.io_manager.open_project_file(str(target))
    qtbot.waitUntil(
        lambda: len(full_window.state_manager.data.atoms) == 6, timeout=30000
    )

    atoms = full_window.state_manager.data.atoms
    assert sum(a.get("charge", 0) for a in atoms.values()) == 1, (
        "formal charge was lost"
    )
    assert sum(a.get("radical", 0) for a in atoms.values()) == 1, (
        "radical state was lost"
    )


def test_round_trip_leaves_the_document_clean(
    full_window, qtbot, monkeypatch, tmp_path
):
    """A freshly loaded project must not look like it has pending edits."""
    place_benzene_template(full_window, qtbot)
    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False
    full_window.io_manager.open_project_file(str(target))
    qtbot.waitUntil(
        lambda: len(full_window.state_manager.data.atoms) == 6, timeout=30000
    )

    assert not full_window.state_manager.has_unsaved_changes
    assert full_window.init_manager.current_file_path == str(target)


def test_reloaded_project_can_be_converted_again(
    full_window, qtbot, monkeypatch, tmp_path
):
    """A round-tripped 2D drawing still drives a correct 2D->3D conversion."""
    place_benzene_template(full_window, qtbot)
    target = _save_project_as(full_window, monkeypatch, tmp_path / "benzene.pmeprj")

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False
    full_window.io_manager.open_project_file(str(target))
    qtbot.waitUntil(
        lambda: len(full_window.state_manager.data.atoms) == 6, timeout=30000
    )

    qtbot.mouseClick(full_window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: full_window.view_3d_manager.current_mol is not None,
        timeout=CONVERT_TIMEOUT_MS,
    )
    assert_is_3d_benzene(full_window.view_3d_manager.current_mol, "reload+convert")
