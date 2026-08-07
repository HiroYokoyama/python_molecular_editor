# -*- coding: utf-8 -*-
"""GUI tests for MoleculeScene mouse and keyboard events."""

import math

import pytest
from PyQt6.QtCore import Qt, QPoint, QPointF
from PyQt6.QtGui import QKeyEvent
from PyQt6.QtTest import QTest

from moleditpy.ui.chain_mixin import CHAIN_HALF_ANGLE_DEG
from moleditpy.utils.constants import DEFAULT_BOND_LENGTH

CHAIN_AXIS_STEP = DEFAULT_BOND_LENGTH * math.cos(math.radians(CHAIN_HALF_ANGLE_DEG))


class FakeSceneMouseEvent:
    """QGraphicsSceneMouseEvent stand-in; PyQt6 forbids instantiating the real one."""

    def __init__(self, pos, button=Qt.MouseButton.LeftButton):
        self._pos = pos
        self._button = button
        self.accepted = False

    def scenePos(self):
        return self._pos

    def pos(self):
        return self._pos

    def button(self):
        return self._button

    def buttons(self):
        return self._button

    def accept(self):
        self.accepted = True

    def ignore(self):
        self.accepted = False


def _drag_chain(scene, start, end):
    scene.mousePressEvent(FakeSceneMouseEvent(start))
    scene.mouseMoveEvent(FakeSceneMouseEvent(end))
    scene.mouseReleaseEvent(FakeSceneMouseEvent(end))


def test_bond_stereo_toggle_keys(window, qtbot):
    """Test Z and E keys toggle double bond stereochemistry."""
    scene = window.init_manager.scene
    data = window.state_manager.data

    # Create two atoms and a double bond
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a2_id = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[a1_id], scene.atom_items[a2_id], bond_order=2)

    bond_key = (min(a1_id, a2_id), max(a1_id, a2_id))
    bond_item = scene.bond_items[bond_key]

    # Simulate hover
    scene.hovered_item = bond_item

    # Press 'Z'
    event_z = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_Z, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_z)
    assert data.bonds[bond_key]["stereo"] == 3
    assert bond_item.stereo == 3

    # Press 'E'
    event_e = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_E, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_e)
    assert data.bonds[bond_key]["stereo"] == 4
    assert bond_item.stereo == 4


def test_atom_addition_keys(window, qtbot):
    """Test 1, 2, 3 keys add atoms/bonds from selected atom."""
    scene = window.init_manager.scene
    data = window.state_manager.data

    # Place initial atom and select it
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a1_item = scene.atom_items[a1_id]
    a1_item.setSelected(True)

    # Press '1' to add a single bond C
    event_1 = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_1, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_1)

    assert len(data.atoms) == 2
    assert len(data.bonds) == 1
    # Check if the new bond is single
    bond_key = list(data.bonds.keys())[0]
    assert data.bonds[bond_key]["order"] == 1

    # Select the new atom
    other_atom_ids = [aid for aid in data.atoms if aid != a1_id]
    assert len(other_atom_ids) == 1, "Expected exactly one new atom after pressing '1'"
    new_atom_id = other_atom_ids[0]
    new_atom_item = scene.atom_items[new_atom_id]
    scene.clearSelection()
    new_atom_item.setSelected(True)

    # Press '2' to add a double bond C
    event_2 = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_2, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_2)

    assert len(data.atoms) == 3
    assert len(data.bonds) == 2

    # Find the newest bond
    for b_key, b_data in data.bonds.items():
        if b_key != bond_key:
            assert b_data["order"] == 2


def test_delete_items_keys(window, qtbot):
    """Test Delete and Backspace keys remove selected items."""
    scene = window.init_manager.scene
    data = window.state_manager.data

    # Create some items
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a2_id = scene.create_atom("C", QPointF(50, 0))
    a1_item = scene.atom_items[a1_id]
    a2_item = scene.atom_items[a2_id]

    a1_item.setSelected(True)

    # Press 'Delete'
    event_del = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_del)

    assert a1_id not in data.atoms
    assert a2_id in data.atoms

    # Select a2 and press 'Backspace'
    a2_item.setSelected(True)
    event_back = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_Backspace, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_back)

    assert a2_id not in data.atoms
    assert len(data.atoms) == 0


def test_temp_line_cancellation(window, qtbot):
    """Test Delete key cancels an active temp_line (bond drawing)."""
    scene = window.init_manager.scene

    a1_id = scene.create_atom("C", QPointF(0, 0))
    a1_item = scene.atom_items[a1_id]

    # Simulate start of bond drawing
    scene.start_atom = a1_item
    scene.start_pos = a1_item.pos()
    from PyQt6.QtWidgets import QGraphicsLineItem

    scene.temp_line = QGraphicsLineItem()
    scene.addItem(scene.temp_line)

    # Press 'Delete'
    event_del = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_del)

    assert scene.temp_line is None


def test_chain_drag_creates_zigzag_chain(window, qtbot):
    """Dragging in chain mode creates one atom per bond at fixed bond length."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    window.ui_manager.set_mode("chain")

    _drag_chain(scene, QPointF(0, 0), QPointF(4 * CHAIN_AXIS_STEP, 0))

    assert len(data.atoms) == 5
    assert len(data.bonds) == 4
    for (id1, id2), bond in data.bonds.items():
        assert bond["order"] == 1
        p1 = scene.atom_items[id1].pos()
        p2 = scene.atom_items[id2].pos()
        assert math.hypot(p2.x() - p1.x(), p2.y() - p1.y()) == pytest.approx(
            DEFAULT_BOND_LENGTH, abs=1e-6
        )


def test_chain_drag_grows_from_an_existing_atom(window, qtbot):
    """A chain started on an existing atom extends it instead of duplicating it."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    anchor_id = scene.create_atom("C", QPointF(0, 0))
    window.ui_manager.set_mode("chain")

    _drag_chain(scene, QPointF(0, 0), QPointF(3 * CHAIN_AXIS_STEP, 0))

    assert anchor_id in data.atoms
    assert len(data.atoms) == 4
    assert len(data.bonds) == 3
    assert any(anchor_id in key for key in data.bonds)


def _drag_and_capture_label(scene, start, end):
    """Drag a chain, returning the label shown mid-drag and the atoms added."""
    data = scene.data
    before = len(data.atoms)
    scene.mousePressEvent(FakeSceneMouseEvent(start))
    scene.mouseMoveEvent(FakeSceneMouseEvent(end))
    label = scene.template_preview.chain_label
    scene.mouseReleaseEvent(FakeSceneMouseEvent(end))
    return label, len(data.atoms) - before


def test_previewed_count_matches_the_atoms_drawn_on_empty_canvas(window, qtbot):
    """The n badge promises a number of atoms; releasing must deliver exactly that."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("chain")

    label, added = _drag_and_capture_label(
        scene, QPointF(0, 0), QPointF(4 * CHAIN_AXIS_STEP, 0)
    )

    assert added == 5
    assert label == "n = 5"


def test_previewed_count_matches_the_atoms_drawn_from_an_existing_atom(window, qtbot):
    """The reused start atom is not newly drawn, so it must not inflate the badge."""
    scene = window.init_manager.scene
    scene.create_atom("C", QPointF(0, 0))
    window.ui_manager.set_mode("chain")

    label, added = _drag_and_capture_label(
        scene, QPointF(0, 0), QPointF(4 * CHAIN_AXIS_STEP, 0)
    )

    assert added == 4
    assert label == "n = 4"


def _bonded_to_symbol(data, symbol):
    return any(
        symbol in (data.atoms[i]["symbol"], data.atoms[j]["symbol"])
        for i, j in data.bonds
    )


def test_chain_end_joins_the_atom_under_the_cursor(window, qtbot):
    """Releasing on an atom joins it, stretching only the tail bond."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    window.ui_manager.set_mode("chain")

    # Sits past where the geometry alone would end, so only a cursor snap reaches it.
    target = QPointF(4 * CHAIN_AXIS_STEP + 8, 0)
    target_id = scene.create_atom("O", target)

    _drag_chain(scene, QPointF(0, 0), target)

    assert target_id in data.atoms
    assert len(data.atoms) == 5, "the hovered atom must be reused, not duplicated"
    assert _bonded_to_symbol(data, "O")

    # Identify the tail by the joined atom, not by sorting: the snap can shorten
    # it as easily as lengthen it, so sorted() would silently pick a different bond.
    def length(i, j):
        return math.hypot(
            scene.atom_items[j].pos().x() - scene.atom_items[i].pos().x(),
            scene.atom_items[j].pos().y() - scene.atom_items[i].pos().y(),
        )

    tail = [length(i, j) for i, j in data.bonds if target_id in (i, j)]
    rest = [length(i, j) for i, j in data.bonds if target_id not in (i, j)]
    assert len(tail) == 1
    assert rest == pytest.approx([DEFAULT_BOND_LENGTH] * 3, abs=1e-6)
    assert tail[0] != pytest.approx(DEFAULT_BOND_LENGTH, abs=1e-6)


def test_chain_end_ignores_an_atom_the_cursor_is_not_on(window, qtbot):
    """An atom well away from the cursor must not capture the chain end."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    window.ui_manager.set_mode("chain")

    end = QPointF(4 * CHAIN_AXIS_STEP, 0)
    scene.create_atom("O", QPointF(end.x() + 60, end.y()))

    _drag_chain(scene, QPointF(0, 0), end)

    assert len(data.atoms) == 6
    assert not _bonded_to_symbol(data, "O")


def test_both_ends_of_a_chain_use_the_bond_snapping_setting(window, qtbot):
    """Start and end share the setting bond drawing uses, not a private radius."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("chain")
    tol = scene.get_setting("bond_snapping_distance_2d", 14.0)

    start_id = scene.create_atom("N", QPointF(0, 0))
    end_id = scene.create_atom("O", QPointF(4 * CHAIN_AXIS_STEP, 0))

    # Both ends are offset by less than the shared snapping distance.
    _drag_chain(
        scene,
        QPointF(tol - 2, 0),
        QPointF(4 * CHAIN_AXIS_STEP + tol - 2, 0),
    )

    data = window.state_manager.data
    assert len(data.atoms) == 5, "both existing atoms must be reused, not duplicated"
    assert any(start_id in key for key in data.bonds)
    assert any(end_id in key for key in data.bonds)


def test_previewed_count_matches_the_atoms_drawn_when_the_end_fuses(window, qtbot):
    """A joined end atom already exists, so the badge must not promise it."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("chain")
    target = QPointF(4 * CHAIN_AXIS_STEP + 8, 0)
    scene.create_atom("O", target)

    label, added = _drag_and_capture_label(scene, QPointF(0, 0), target)

    assert added == 4
    assert label == "n = 4"


def test_chain_preview_follows_a_real_mouse_drag(window, qtbot):
    """Real Qt events through the view must reach the live preview, not just fakes."""
    view = window.init_manager.view_2d
    scene = window.init_manager.scene
    window.ui_manager.set_mode("chain")
    view.show()
    qtbot.waitExposed(view)

    start = view.viewport().rect().center()
    far = QPoint(start.x() + 250, start.y())

    QTest.mousePress(view.viewport(), Qt.MouseButton.LeftButton, pos=start)
    assert scene.chain_active is True

    QTest.mouseMove(view.viewport(), far)
    qtbot.wait(20)
    previewed_bonds = len(scene.chain_points) - 1
    assert previewed_bonds > 0
    assert scene.template_preview.isVisible()

    QTest.mouseRelease(view.viewport(), Qt.MouseButton.LeftButton, pos=far)
    assert scene.template_preview.isVisible() is False
    assert len(window.state_manager.data.bonds) == previewed_bonds


def test_chain_tool_button_sits_between_the_9_ring_and_user_templates(window, qtbot):
    """The chain tool shares the exclusive tool group and closes the template row."""
    init = window.init_manager
    chain_action = init.mode_actions["chain"]
    toolbar_actions = init.toolbar_bottom.actions()

    ring9_index = toolbar_actions.index(init.mode_actions["template_9"])
    assert toolbar_actions[ring9_index + 1] is chain_action
    assert toolbar_actions[ring9_index + 2] is init.mode_actions["template_user"]
    assert chain_action.isCheckable()
    assert chain_action in init.tool_group.actions()
    assert not chain_action.icon().isNull()


def test_chain_click_without_drag_adds_nothing(window, qtbot):
    """A bare click in chain mode is not long enough to place a bond."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    window.ui_manager.set_mode("chain")

    _drag_chain(scene, QPointF(0, 0), QPointF(2, 0))

    assert len(data.atoms) == 0
    assert len(data.bonds) == 0
    assert scene.chain_active is False


def test_bonding_to_existing_atom(window, qtbot):
    """Test that pressing 1, 2, 3 bonds to an existing atom if it's nearby."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    from moleditpy.utils.constants import DEFAULT_BOND_LENGTH

    # Create two atoms separated by DEFAULT_BOND_LENGTH
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a2_id = scene.create_atom("C", QPointF(0, -DEFAULT_BOND_LENGTH))

    a1_item = scene.atom_items[a1_id]
    a1_item.setSelected(True)

    # Press '1'. Since (0, -L) is exactly where a2 is, it should bond to a2
    # instead of creating a new atom.
    event_1 = QKeyEvent(
        QKeyEvent.Type.KeyPress, Qt.Key.Key_1, Qt.KeyboardModifier.NoModifier
    )
    scene.keyPressEvent(event_1)

    assert len(data.atoms) == 2  # No new atom created
    assert len(data.bonds) == 1  # Bond created between a1 and a2

    bond_key = (min(a1_id, a2_id), max(a1_id, a2_id))
    assert bond_key in data.bonds
