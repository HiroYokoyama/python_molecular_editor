import pytest
import math
from PyQt6.QtCore import Qt, QPointF, QLineF
from PyQt6.QtGui import QKeyEvent, QTransform, QMouseEvent
from PyQt6.QtWidgets import QApplication, QGraphicsLineItem
from moleditpy.modules.molecule_scene import MoleculeScene
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from unittest.mock import MagicMock, patch


@pytest.fixture(autouse=True)
def mock_graphics_scene_events():
    """Patch QGraphicsScene mouse events to accept MagicMock."""
    with (
        patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent"),
        patch("PyQt6.QtWidgets.QGraphicsScene.mouseMoveEvent"),
        patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent"),
        patch("PyQt6.QtWidgets.QGraphicsScene.keyPressEvent"),
    ):
        yield


def setup_scene_with_view(mock_parser_host):
    scene = MoleculeScene(mock_parser_host.data, mock_parser_host)
    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()
    mock_view.mapToScene.return_value = QPointF(0, 0)
    mock_view.mapFromGlobal.return_value = QPointF(0, 0)
    scene.views = MagicMock(return_value=[mock_view])
    scene.template_preview = MagicMock()
    scene.window.is_2d_editable = True
    scene.mode = "select"
    return scene


def create_mock_mouse_event(pos, button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.scenePos.return_value = pos
    event.button.return_value = button
    return event


def test_scene_keypress_modes(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)

    def mock_set_mode(mode):
        scene.mode = mode

    scene.window.set_mode_and_update_toolbar.side_effect = mock_set_mode

    # Test space -> activate_select_mode
    scene.mode = "atom_C"
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch("PyQt6.QtGui.QCursor.pos", return_value=QPointF(0, 0)),
    ):
        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_Space, Qt.KeyboardModifier.NoModifier
        )
        scene.keyPressEvent(event)
        assert scene.window.activate_select_mode.called

    key_map = {
        Qt.Key.Key_C: "atom_C",
        Qt.Key.Key_H: "atom_H",
        Qt.Key.Key_O: "atom_O",
        Qt.Key.Key_N: "atom_N",
        Qt.Key.Key_S: "atom_S",
        Qt.Key.Key_P: "atom_P",
        Qt.Key.Key_F: "atom_F",
        Qt.Key.Key_B: "atom_B",
        Qt.Key.Key_I: "atom_I",
    }

    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch("PyQt6.QtGui.QCursor.pos", return_value=QPointF(0, 0)),
    ):
        for key, expected_mode in key_map.items():
            event = QKeyEvent(
                QKeyEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier
            )
            scene.keyPressEvent(event)
            assert scene.mode == expected_mode


def test_scene_keypress_special_symbols(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)

    def mock_set_mode(mode):
        scene.mode = mode

    scene.window.set_mode_and_update_toolbar.side_effect = mock_set_mode

    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch("PyQt6.QtGui.QCursor.pos", return_value=QPointF(0, 0)),
    ):
        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_C, Qt.KeyboardModifier.ShiftModifier
        )
        scene.keyPressEvent(event)
        assert scene.mode == "atom_Cl"

        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_B, Qt.KeyboardModifier.ShiftModifier
        )
        scene.keyPressEvent(event)
        assert scene.mode == "atom_Br"

        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_S, Qt.KeyboardModifier.ShiftModifier
        )
        scene.keyPressEvent(event)
        assert scene.mode == "atom_Si"


def test_scene_keypress_delete(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]["item"]
    atom_item.setSelected(True)

    with (
        patch.object(MoleculeScene, "delete_items") as mock_delete,
        patch.object(MoleculeScene, "selectedItems", return_value=[atom_item]),
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch("PyQt6.QtGui.QCursor.pos", return_value=QPointF(0, 0)),
    ):
        event = QKeyEvent(
            QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier
        )
        scene.keyPressEvent(event)
        mock_delete.assert_called()


def test_scene_maintenance_methods(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    aid = scene.create_atom("C", QPointF(0, 0))
    # mock_parser_host's create_atom should populate data.atoms
    atom_item = mock_parser_host.data.atoms[aid]["item"]

    atom_item.has_problem = True
    scene.clear_all_problem_flags()
    assert atom_item.has_problem is False
    scene.purge_deleted_items()
    scene.reinitialize_items()
    scene.update_all_items()


def test_scene_queries(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    a1 = MagicMock(spec=AtomItem)
    a2 = MagicMock(spec=AtomItem)
    a1.pos.return_value = QPointF(0, 0)
    a1.scenePos.return_value = QPointF(0, 0)
    a1.__class__ = AtomItem
    a2.pos.return_value = QPointF(50, 0)
    a2.scenePos.return_value = QPointF(50, 0)
    a2.__class__ = AtomItem

    with patch.object(scene, "items", return_value=[a1, a2]):
        assert scene.find_atom_near(QPointF(50, 0)) == a2


def test_scene_update_connected_bonds(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    a = MagicMock(spec=AtomItem)
    a.bonds = [MagicMock()]
    scene.update_connected_bonds([a])
    assert a.bonds[0].update_position.called


def test_scene_leave_event(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.leaveEvent(None)
    assert scene.template_preview.hide.called


def test_scene_update_bond_stereo(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    bond = MagicMock()
    bond.order = 2
    bond.stereo = 0
    bond.atom1.atom_id = 10
    bond.atom2.atom_id = 20
    scene.data.bonds[(10, 20)] = {"stereo": 0, "item": bond}

    scene.update_bond_stereo(bond, 1)
    assert scene.data.bonds[(10, 20)]["stereo"] == 1
    assert bond.set_stereo.called


def test_scene_mouse_drag_create_bond_existing_atoms(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "atom_C"
    scene.current_atom_symbol = "C"

    aid1 = scene.create_atom("C", QPointF(0, 0))
    aid2 = scene.create_atom("C", QPointF(50, 0))
    a1 = mock_parser_host.data.atoms[aid1]["item"]
    a2 = mock_parser_host.data.atoms[aid2]["item"]

    # 1. Press on A1
    press_event = create_mock_mouse_event(QPointF(0, 0))
    with patch.object(scene, "itemAt", return_value=a1):
        scene.mousePressEvent(press_event)

    # 2. Release on A2
    release_event = create_mock_mouse_event(QPointF(50, 0))
    with (
        patch("PyQt6.QtWidgets.QApplication.startDragDistance", return_value=10),
        patch.object(scene, "itemAt", return_value=a2),
    ):
        scene.mouseReleaseEvent(release_event)

    # Verify bond creation via mock_parser_host's data.add_bond
    assert (aid1, aid2) in scene.data.bonds or (aid2, aid1) in scene.data.bonds


def test_scene_mouse_click_create_single_atom(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "atom_O"
    scene.current_atom_symbol = "O"

    # Click in blank space
    click_event = create_mock_mouse_event(QPointF(100, 100))
    with (
        patch.object(scene, "itemAt", return_value=None),
        patch("PyQt6.QtWidgets.QApplication.startDragDistance", return_value=10),
    ):
        # mousePressEvent sets self.start_pos and self.press_pos
        scene.mousePressEvent(click_event)
        # mouseReleaseEvent calculates is_click and triggers atom creation if start_pos matches current pos
        scene.mouseReleaseEvent(click_event)

    # Verify an 'O' atom was created and undo state pushed
    assert any(a["symbol"] == "O" for a in scene.data.atoms.values())
    assert mock_parser_host.push_undo_state.called


def test_scene_right_click_bond_delete(qtbot, mock_parser_host):
    """Test right-click deletion on a bond."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "select"
    aid1 = scene.create_atom("C", QPointF(0, 0))
    aid2 = scene.create_atom("C", QPointF(50, 0))
    a1 = mock_parser_host.data.atoms[aid1]["item"]
    a2 = mock_parser_host.data.atoms[aid2]["item"]
    scene.create_bond(a1, a2)
    bond_item = list(mock_parser_host.data.bonds.values())[0]["item"]

    event = create_mock_mouse_event(QPointF(25, 0), button=Qt.MouseButton.RightButton)
    with (
        patch.object(scene, "itemAt", return_value=bond_item),
        patch.object(
            MoleculeScene, "delete_items", wraps=scene.delete_items
        ) as mock_del,
    ):
        scene.mousePressEvent(event)
        mock_del.assert_called()
        assert len(mock_parser_host.data.bonds) == 0
        assert mock_parser_host.push_undo_state.called

    # Safety for headless environments: ensure any popup menus are closed
    # before the test finishes to avoid segmentation faults during teardown.
    for widget in QApplication.topLevelWidgets():
        try:
            if widget.metaObject().className() == "QMenu":
                widget.close()
        except (RuntimeError, AttributeError):
            continue

    # Wait for the event loop to clear any pending events (e.g. deletion, menu fade out)
    qtbot.wait(50)


def test_scene_drag_and_drop_atom(mock_parser_host):
    """Test moving an atom via drag-and-drop."""
    scene = setup_scene_with_view(mock_parser_host)
    start_pos = QPointF(0, 0)
    aid = scene.create_atom("C", start_pos)
    atom_item = mock_parser_host.data.atoms[aid]["item"]

    scene.mode = "select"
    new_pos = QPointF(50, 50)

    # 1. Press
    press_event = create_mock_mouse_event(start_pos)
    with (
        patch.object(scene, "items", return_value=[atom_item]),
        patch.object(scene, "itemAt", return_value=atom_item),
    ):
        scene.mousePressEvent(press_event)

        # 2. Move
        atom_item.setPos(new_pos)
        move_event = create_mock_mouse_event(new_pos)
        with patch("PyQt6.QtWidgets.QApplication.startDragDistance", return_value=5):
            scene.mouseMoveEvent(move_event)

            # 3. Release
            scene.mouseReleaseEvent(move_event)

    assert atom_item.pos() == new_pos
    assert mock_parser_host.data.atoms[aid]["pos"] == new_pos
    assert mock_parser_host.push_undo_state.called


def test_scene_delete_mixed_selection(mock_parser_host):
    """Test deleting a selection containing both atoms and bonds."""
    scene = setup_scene_with_view(mock_parser_host)
    aid1 = scene.create_atom("C", QPointF(0, 0))
    aid2 = scene.create_atom("C", QPointF(50, 0))
    a1 = mock_parser_host.data.atoms[aid1]["item"]
    a2 = mock_parser_host.data.atoms[aid2]["item"]
    scene.create_bond(a1, a2)
    bond_item = list(mock_parser_host.data.bonds.values())[0]["item"]

    a1.setSelected(True)
    bond_item.setSelected(True)

    scene.delete_items({a1, bond_item})

    assert aid1 not in mock_parser_host.data.atoms
    assert len(mock_parser_host.data.bonds) == 0
    assert aid2 in mock_parser_host.data.atoms
