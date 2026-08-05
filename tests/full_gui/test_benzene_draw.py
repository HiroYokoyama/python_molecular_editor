# -*- coding: utf-8 -*-
"""Draw benzene in the running application and check the 2D model it produced."""

import math

import pytest
from PyQt6.QtCore import QPoint, QPointF, Qt

from benzene import (
    RING_RADIUS,
    assert_is_2d_benzene,
    draw_benzene_by_hand,
    place_benzene_template,
)

pytestmark = pytest.mark.full_gui


def _centroid(points):
    n = len(points)
    return QPointF(
        sum(p.x() for p in points) / n,
        sum(p.y() for p in points) / n,
    )


def test_benzene_template_places_a_ring(full_window, qtbot):
    """Toolbar benzene tool + hover + click puts a Kekule ring in the model."""
    place_benzene_template(full_window, qtbot)
    assert_is_2d_benzene(full_window.state_manager.data)


def test_benzene_template_creates_scene_items(full_window, qtbot):
    """The ring is not just data -- six AtomItems and six BondItems exist."""
    place_benzene_template(full_window, qtbot)
    scene = full_window.init_manager.scene
    assert len(scene.atom_items) == 6
    assert len(scene.bond_items) == 6


def test_bond_drawn_by_real_mouse_drag(full_window, qtbot):
    """Press-drag-release on the canvas creates two atoms and a bond.

    The one drawing path driven entirely by real Qt mouse events.
    """
    view = full_window.init_manager.view_2d
    viewport = view.viewport()
    full_window.ui_manager.set_mode("atom_C")

    centre = viewport.rect().center()
    start = QPoint(centre.x() - 60, centre.y())
    end = QPoint(centre.x() + 60, centre.y())

    qtbot.mousePress(viewport, Qt.MouseButton.LeftButton, pos=start)
    qtbot.mouseMove(viewport, end)
    qtbot.mouseRelease(viewport, Qt.MouseButton.LeftButton, pos=end)

    data = full_window.state_manager.data
    qtbot.waitUntil(lambda: len(data.bonds) == 1, timeout=10000)
    assert len(data.atoms) == 2
    assert [a["symbol"] for a in data.atoms.values()] == ["C", "C"]


def test_benzene_by_hand_matches_template(full_window):
    """Hand-built benzene yields the same model the template does."""
    draw_benzene_by_hand(full_window)
    assert_is_2d_benzene(full_window.state_manager.data)


def test_benzene_is_aromatic_to_rdkit(full_window, qtbot):
    """RDKit perceives the drawn ring as an aromatic six-ring."""
    from rdkit import Chem

    place_benzene_template(full_window, qtbot)
    mol = full_window.state_manager.data.to_rdkit_mol()
    assert mol is not None

    Chem.SanitizeMol(mol)
    assert mol.GetNumAtoms() == 6
    assert all(a.GetIsAromatic() for a in mol.GetAtoms()), "ring is not aromatic"
    assert Chem.MolToSmiles(mol) == "c1ccccc1"


def test_benzene_implicit_hydrogens_are_one_per_carbon(full_window, qtbot):
    """Each aromatic CH shows exactly one implicit hydrogen in the scene."""
    place_benzene_template(full_window, qtbot)
    full_window.edit_actions_manager.update_implicit_hydrogens()
    qtbot.wait(200)

    counts = sorted(
        a.implicit_h_count for a in full_window.init_manager.scene.atom_items.values()
    )
    assert counts == [1] * 6, f"expected one H per carbon, got {counts}"


def test_no_chemistry_problems_flagged(full_window, qtbot):
    """A valid benzene must not raise the red 'problem' marker on any atom."""
    place_benzene_template(full_window, qtbot)
    full_window.edit_actions_manager.update_implicit_hydrogens()
    qtbot.wait(200)

    flagged = [
        a.atom_id
        for a in full_window.init_manager.scene.atom_items.values()
        if a.has_problem
    ]
    assert not flagged, f"atoms wrongly flagged as problematic: {flagged}"


def test_benzene_undo_restores_empty_canvas(full_window, qtbot):
    """Ctrl+Z after placing the template empties the model again."""
    place_benzene_template(full_window, qtbot)
    assert len(full_window.state_manager.data.atoms) == 6

    full_window.edit_actions_manager.undo()
    qtbot.wait(100)
    assert len(full_window.state_manager.data.atoms) == 0

    full_window.edit_actions_manager.redo()
    qtbot.wait(100)
    assert_is_2d_benzene(full_window.state_manager.data)


def test_second_ring_fuses_into_naphthalene(full_window, qtbot):
    """Clicking the template on an existing bond fuses rather than overlaps."""
    place_benzene_template(full_window, qtbot)

    scene = full_window.init_manager.scene
    view = full_window.init_manager.view_2d
    viewport = view.viewport()

    bond_item = next(iter(scene.bond_items.values()))
    mid = (bond_item.atom1.pos() + bond_item.atom2.pos()) / 2.0

    # Exactly on the midpoint the new ring's side is undefined and it lands on
    # top of the old one, fusing every vertex; nudge outward from the centre.
    centre = _centroid([a.pos() for a in scene.atom_items.values()])
    outward = mid - centre
    scale = 3.0 / math.hypot(outward.x(), outward.y())
    probe = mid + outward * scale
    pos = viewport.mapFrom(view, view.mapFromScene(probe))

    full_window.init_manager.mode_actions["template_benzene"].trigger()

    # Arm the preview the way the hover handler does; the fusion is what is
    # under test here, not Qt's hit-testing of a one-pixel-wide bond.
    scene.update_template_preview(probe)
    ctx = scene.template_context
    # "items" holds the two atoms of the bond the new ring will share.
    assert len(ctx.get("items", [])) == 2, (
        f"template did not latch onto the bond; context items={ctx.get('items')}"
    )
    assert ctx.get("points"), "template preview produced no ring points"

    qtbot.mouseClick(viewport, Qt.MouseButton.LeftButton, pos=pos)

    data = full_window.state_manager.data
    qtbot.waitUntil(lambda: len(data.atoms) > 6, timeout=10000)
    assert len(data.atoms) == 10, (
        f"fused bicyclic should have 10 carbons, got {len(data.atoms)}"
    )
    assert len(data.bonds) == 11, f"expected 11 bonds, got {len(data.bonds)}"


def test_ring_geometry_is_a_regular_hexagon(full_window):
    """Hand-drawn ring keeps equal bond lengths, so 2D layout is sane."""
    ids = draw_benzene_by_hand(full_window)
    scene = full_window.init_manager.scene
    lengths = []
    for i in range(6):
        p = scene.atom_items[ids[i]].pos()
        q = scene.atom_items[ids[(i + 1) % 6]].pos()
        lengths.append(math.hypot(q.x() - p.x(), q.y() - p.y()))

    # A regular hexagon's side equals its circumradius.
    assert max(lengths) - min(lengths) < 1e-6, f"uneven ring: {lengths}"
    assert abs(lengths[0] - RING_RADIUS) < 1e-6
