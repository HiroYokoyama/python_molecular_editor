# -*- coding: utf-8 -*-
"""Shared benzene helpers for the full-GUI tier."""

import math

from PyQt6.QtCore import QPoint, QPointF, Qt

RING_RADIUS = 45.0
CENTER = QPointF(0.0, 0.0)


def ring_positions(center: QPointF = CENTER, radius: float = RING_RADIUS):
    """Six vertices of a regular hexagon, first vertex pointing up."""
    return [
        QPointF(
            center.x() + radius * math.cos(-math.pi / 2 + i * math.pi / 3),
            center.y() + radius * math.sin(-math.pi / 2 + i * math.pi / 3),
        )
        for i in range(6)
    ]


def draw_benzene_by_hand(window):
    """Build benzene atom-by-atom through the scene's own creation API.

    Deterministic (no hit-testing), so tests that are really about conversion
    or rendering do not fail for drawing reasons.
    """
    scene = window.init_manager.scene
    window.ui_manager.set_mode("atom_C")

    ids = [scene.create_atom("C", p) for p in ring_positions()]
    for i in range(6):
        a, b = scene.atom_items[ids[i]], scene.atom_items[ids[(i + 1) % 6]]
        scene.create_bond(a, b, bond_order=2 if i % 2 == 0 else 1)

    scene.update_all_items()
    return ids


def place_benzene_template(window, qtbot, viewport_pos: QPoint = None):
    """Place the benzene ring the way a user does: pick the tool, hover, click.

    Arming the preview is a real step, not cosmetic: the release handler only
    commits a template when `template_context` already holds points, and it is
    `MoleculeScene.mouseMoveEvent` that fills it. Synthesized hover moves do
    not reach that handler dependably without a physical cursor, so the hover's
    own callee -- `update_template_preview` -- is invoked directly and the
    commit is then driven by a genuine mouse click.
    """
    view = window.init_manager.view_2d
    viewport = view.viewport()
    if viewport_pos is None:
        viewport_pos = viewport.rect().center()

    scene = window.init_manager.scene
    window.init_manager.mode_actions["template_benzene"].trigger()
    assert scene.mode == "template_benzene"

    scene_pos = view.mapToScene(viewport.mapTo(view, viewport_pos))
    scene.update_template_preview(scene_pos)
    assert scene.template_context.get("points"), "template preview was not armed"

    qtbot.mouseClick(viewport, Qt.MouseButton.LeftButton, pos=viewport_pos)
    qtbot.waitUntil(lambda: len(window.state_manager.data.atoms) > 0, timeout=10000)
    return scene_pos


def assert_is_2d_benzene(data):
    """The 2D model holds a six-carbon ring with 3 double + 3 single bonds."""
    symbols = [a["symbol"] for a in data.atoms.values()]
    assert len(data.atoms) == 6, f"expected 6 atoms, got {len(data.atoms)}: {symbols}"
    assert symbols == ["C"] * 6, f"expected all carbon, got {symbols}"
    assert len(data.bonds) == 6, f"expected 6 bonds, got {len(data.bonds)}"

    orders = sorted(b["order"] for b in data.bonds.values())
    assert orders == [1, 1, 1, 2, 2, 2], f"expected Kekule benzene, got {orders}"

    # Every carbon must carry exactly two ring bonds.
    degree = {aid: 0 for aid in data.atoms}
    for a1, a2 in data.bonds:
        degree[a1] += 1
        degree[a2] += 1
    assert set(degree.values()) == {2}, f"not a simple ring: {degree}"
