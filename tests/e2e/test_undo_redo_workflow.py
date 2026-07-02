# -*- coding: utf-8 -*-
"""E2E undo/redo workflow tests through the real StateManager/EditActionsManager.

Covers full user cycles: draw -> undo -> redo, bond-data restoration, and the
undo-history depth cap.
"""

import importlib
import importlib.util

import pytest
from PyQt6.QtCore import QPointF

_PKG = (
    "moleditpy_linux"
    if importlib.util.find_spec("moleditpy_linux") is not None
    else "moleditpy"
)
UNDO_STACK_MAX_DEPTH = importlib.import_module(
    f"{_PKG}.utils.constants"
).UNDO_STACK_MAX_DEPTH


@pytest.mark.gui
def test_draw_undo_redo_cycle(window):
    """Two draw steps can be fully undone and redone with exact atom counts."""
    mgr = window.edit_actions_manager
    mgr.reset_history()  # empty baseline entry

    scene = window.init_manager.scene
    scene.create_atom("C", QPointF(0, 0))
    mgr.push_undo_state()
    scene.create_atom("O", QPointF(50, 0))
    mgr.push_undo_state()
    assert len(window.data.atoms) == 2

    mgr.undo()
    assert len(window.data.atoms) == 1
    mgr.undo()
    assert len(window.data.atoms) == 0
    mgr.undo()  # at stack floor: no-op
    assert len(window.data.atoms) == 0

    mgr.redo()
    assert len(window.data.atoms) == 1
    mgr.redo()
    assert len(window.data.atoms) == 2
    symbols = sorted(d["symbol"] for d in window.data.atoms.values())
    assert symbols == ["C", "O"]


@pytest.mark.gui
def test_undo_restores_bond_order(window):
    """Undoing a bond-order change restores the original order."""
    mgr = window.edit_actions_manager
    scene = window.init_manager.scene
    c1 = scene.create_atom("C", QPointF(0, 0))
    c2 = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[c1], scene.atom_items[c2], bond_order=1)
    mgr.reset_history()

    key = next(iter(window.data.bonds))
    assert window.data.bonds[key]["order"] == 1

    window.data.add_bond(c1, c2, order=2)
    mgr.push_undo_state()
    assert window.data.bonds[key]["order"] == 2

    mgr.undo()
    key_after = next(iter(window.data.bonds))
    assert window.data.bonds[key_after]["order"] == 1
    assert len(window.data.bonds) == 1


@pytest.mark.gui
def test_undo_stack_capped_e2e(window):
    """The undo history never exceeds UNDO_STACK_MAX_DEPTH entries."""
    mgr = window.edit_actions_manager
    scene = window.init_manager.scene
    aid = scene.create_atom("C", QPointF(0, 0))
    mgr.reset_history()

    for i in range(1, UNDO_STACK_MAX_DEPTH + 6):
        window.data.set_atom_pos(aid, (float(i), 0.0))
        mgr.push_undo_state()

    assert len(mgr.undo_stack) == UNDO_STACK_MAX_DEPTH
    # Newest state is the last position pushed
    newest = mgr.undo_stack[-1]["atoms"][aid]["pos"]
    assert tuple(newest) == pytest.approx((float(UNDO_STACK_MAX_DEPTH + 5), 0.0))
