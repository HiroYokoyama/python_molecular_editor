"""Unit tests for MoleculeScene mouse event handling."""

import math
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QTransform
from moleditpy.ui.molecule_scene import MoleculeScene
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem
from unittest.mock import MagicMock, patch


def setup_scene_with_view(mock_parser_host):
    scene = MoleculeScene(mock_parser_host.data, mock_parser_host)
    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()
    mock_view.mapToScene.return_value = QPointF(0, 0)
    mock_view.mapFromGlobal.return_value = QPointF(0, 0)
    scene.views = MagicMock(return_value=[mock_view])
    return scene


def create_mock_event(pos=None, button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.scenePos.return_value = QPointF(100, 100) if pos is None else pos
    event.button.return_value = button
    return event


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_scene_toggle_radical(mock_drag, mock_release, mock_press, mock_parser_host):
    """Radical mode cycles atom radical count through 0→1→2→0 on successive clicks."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "radical"
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = scene.atom_items[aid]
    atom_item.radical = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()

    with (
        patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
        patch.object(MoleculeScene, "itemAt", return_value=atom_item),
    ):
        event = create_mock_event(QPointF(100, 100))
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.radical == 1
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.radical == 2
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.radical == 0


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_scene_toggle_charge(mock_drag, mock_release, mock_press, mock_parser_host):
    """Charge plus/minus modes increment and decrement atom charge correctly."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "charge_plus"
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = scene.atom_items[aid]
    atom_item.charge = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()

    with (
        patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
        patch.object(MoleculeScene, "itemAt", return_value=atom_item),
    ):
        event = create_mock_event(QPointF(100, 100))
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.charge == 1
        scene.mode = "charge_minus"
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.charge == 0


def test_add_benzene_fragment(mock_parser_host):
    """add_molecule_fragment creates 6 atoms and 6 bonds with correct alternating order."""
    scene = setup_scene_with_view(mock_parser_host)
    points = [
        QPointF(
            20 * math.cos(math.radians(i * 60)), 20 * math.sin(math.radians(i * 60))
        )
        for i in range(6)
    ]
    bonds_info = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
    scene.add_molecule_fragment(points, bonds_info)
    assert len(scene.data.atoms) == 6
    assert len(scene.data.bonds) == 6
    assert len([b for b in scene.data.bonds.values() if b["order"] == 2]) == 3


def test_benzene_fusion_rotation(mock_parser_host):
    """Benzene fusion reuses existing edge atoms and preserves bond order."""
    scene = setup_scene_with_view(mock_parser_host)
    id1 = scene.create_atom("C", QPointF(-10, 0))
    id2 = scene.create_atom("C", QPointF(10, 0))
    scene.create_bond(scene.atom_items[id1], scene.atom_items[id2], bond_order=2)
    points = [QPointF(10, 0), QPointF(-10, 0)] + [
        QPointF(
            20 * math.cos(math.radians(i * 60)), 20 * math.sin(math.radians(i * 60))
        )
        for i in range(2, 6)
    ]
    bonds_info = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
    scene.add_molecule_fragment(
        points,
        bonds_info,
        existing_items=[scene.atom_items[id2], scene.atom_items[id1]],
    )
    assert len(scene.data.atoms) == 6
    eb = scene.find_bond_between(scene.atom_items[id1], scene.atom_items[id2])
    assert eb.order == 2


def test_delete_selected_items(mock_parser_host):
    """Right-click on a selected atom triggers delete_items."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "select"
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = scene.atom_items[aid]
    atom_item.setSelected(True)
    with patch.object(MoleculeScene, "delete_items", return_value=True) as mock_delete:
        with patch.object(MoleculeScene, "selectedItems", return_value=[atom_item]):
            with patch.object(MoleculeScene, "itemAt", return_value=atom_item):
                event = MagicMock()
                event.button.return_value = Qt.MouseButton.RightButton
                event.scenePos.return_value = QPointF(0, 0)
                scene.mousePressEvent(event)
                mock_delete.assert_called()


@patch("PyQt6.QtWidgets.QGraphicsScene.mouseDoubleClickEvent")
def test_double_click_select_component(mock_dbl, mock_parser_host):
    """Double-click selects all atoms in the connected component, not isolated ones."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "select"
    id1 = scene.create_atom("C", QPointF(0, 0))
    id2 = scene.create_atom("C", QPointF(20, 0))
    id3 = scene.create_atom("C", QPointF(40, 0))
    scene.create_bond(scene.atom_items[id1], scene.atom_items[id2])
    scene.create_bond(scene.atom_items[id2], scene.atom_items[id3])
    id_iso = scene.create_atom("O", QPointF(100, 100))
    atom_iso = scene.atom_items[id_iso]
    for it in scene.atom_items.values():
        it.__class__ = AtomItem
        it.setSelected = MagicMock()
    with patch.object(MoleculeScene, "itemAt", return_value=scene.atom_items[id1]):
        _orig_isinstance = isinstance

        def mock_isinstance(obj, types):
            tpl = types if _orig_isinstance(types, tuple) else (types,)
            if hasattr(obj, "atom_id") and AtomItem in tpl:
                return True
            if hasattr(obj, "atom1") and BondItem in tpl:
                return True
            return _orig_isinstance(obj, types)

        with patch(
            "moleditpy.ui.molecule_scene.isinstance", side_effect=mock_isinstance
        ):
            event = create_mock_event(QPointF(0, 0))
            scene.mouseDoubleClickEvent(event)
            assert scene.atom_items[id1].setSelected.called
            assert scene.atom_items[id2].setSelected.called
            assert scene.atom_items[id3].setSelected.called
            assert not atom_iso.setSelected.called


@patch("PyQt6.QtWidgets.QGraphicsScene.keyPressEvent")
def test_scene_key_event_dispatch(mock_kp, mock_parser_host):
    """Test key events like '4' (template), '.' (radical), +/- (charge)."""
    scene = setup_scene_with_view(mock_parser_host)

    # Case 1: Key '4' -> Template mode when hovering nothing
    # Mock itemAt returning None
    with patch.object(MoleculeScene, "itemAt", return_value=None):
        event = MagicMock()
        event.key.return_value = Qt.Key.Key_4
        event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

        scene.keyPressEvent(event)

        # Verify mode switched to 'template_benzene'
        mock_parser_host.ui_manager.set_mode_and_update_toolbar.assert_called_with(
            "template_benzene"
        )

    # Case 2: Key '+' -> Charge increase on hovered atom
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = scene.atom_items[aid]
    atom_item.charge = 0

    with patch.object(MoleculeScene, "itemAt", return_value=atom_item):
        event = MagicMock()
        event.key.return_value = Qt.Key.Key_Plus
        event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

        scene.keyPressEvent(event)

        assert atom_item.charge == 1
        assert scene.data.atoms[aid]["charge"] == 1


def test_scene_update_template_preview_logic(mock_parser_host):
    """Test update_template_preview with different targets."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "template_benzene"  # 6-ring

    scene.template_preview = MagicMock()

    # Case 1: Hovering nothing
    with patch.object(MoleculeScene, "items", return_value=[]):
        scene.update_template_preview(QPointF(100, 100))
        # Should show preview at pos
        scene.template_preview.set_geometry.assert_called()
        args, _ = scene.template_preview.set_geometry.call_args
        points = args[0]
        assert len(points) == 6

    # Case 2: Hovering Atom
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = scene.atom_items[aid]

    with patch.object(MoleculeScene, "items", return_value=[atom_item]):
        scene.update_template_preview(QPointF(20, 0))  # Mouse strictly to right

        # Should fuse to atom
        scene.template_preview.set_geometry.assert_called()
        # Verify context was set
        assert "items" in scene.template_context
        assert scene.template_context["items"][0] == atom_item


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseMoveEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
def test_scene_drag_create_bond_sequence(
    mock_release, mock_move, mock_press, mock_parser_host
):
    """Test the full mouse press -> move -> release sequence for creating a bond."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "bond_1"  # Single bond mode

    # Atom A at 0,0
    id1 = scene.create_atom("C", QPointF(0, 0))
    a1 = scene.atom_items[id1]

    # Atom B at 100,100
    id2 = scene.create_atom("C", QPointF(100, 100))
    a2 = scene.atom_items[id2]

    # 1. Press on A1
    with patch.object(MoleculeScene, "itemAt", return_value=a1):
        event_press = create_mock_event(QPointF(0, 0))
        scene.mousePressEvent(event_press)
        # Verify start_atom is set
        assert getattr(scene, "start_atom", None) == a1

    # 2. Move to A2
    # Logic might update temp line
    event_move = create_mock_event(QPointF(50, 50))
    scene.mouseMoveEvent(event_move)

    # 3. Release on A2
    with patch.object(MoleculeScene, "itemAt", return_value=a2):
        event_release = create_mock_event(QPointF(100, 100))
        scene.mouseReleaseEvent(event_release)

    # Verify bond created
    bond = scene.find_bond_between(a1, a2)
    assert bond.order == 1


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseMoveEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
def test_scene_bond_snapping_distance(
    mock_release, mock_move, mock_press, mock_parser_host
):
    """Verify that mouse interactions respect the configured bond_snapping_distance_2d setting."""
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "bond_1"

    # Configure settings
    mock_parser_host.init_manager.settings = {"bond_snapping_distance_2d": 25.0}

    id1 = scene.create_atom("C", QPointF(0, 0))
    scene.atom_items[id1]

    # Press near A1 but not exactly on it (within 25.0 tolerance)
    with patch.object(scene, "find_atom_near") as mock_find_near:
        event_press = create_mock_event(QPointF(10, 10))
        scene.mousePressEvent(event_press)
        mock_find_near.assert_called_with(QPointF(10, 10), tol=25.0)


def _keys_in_sync(scene):
    """bond_items must be keyed exactly like data.bonds, one item per key."""
    return set(scene.bond_items) == set(scene.data.bonds)


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_clicking_a_wedge_again_keeps_bond_items_keyed_like_the_data(
    mock_drag, mock_release, mock_press, mock_parser_host
):
    """Flipping a wedge re-keys data.bonds; bond_items used to keep the old key."""
    scene = setup_scene_with_view(mock_parser_host)
    a = scene.atom_items[scene.create_atom("C", QPointF(0, 0))]
    b = scene.atom_items[scene.create_atom("C", QPointF(50, 0))]
    scene.create_bond(a, b, bond_order=1, bond_stereo=1)
    bond = scene.find_bond_between(a, b)
    scene.mode, scene.bond_order, scene.bond_stereo = "bond_1_1", 1, 1

    with patch.object(MoleculeScene, "itemAt", return_value=bond):
        event = create_mock_event(QPointF(25, 0))
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)

    assert (b.atom_id, a.atom_id) in scene.data.bonds
    assert _keys_in_sync(scene)


def test_keyboard_wedge_flip_keeps_bond_items_keyed_like_the_data(mock_parser_host):
    """W on a wedge flips it; the scene's bond index follows the data key."""
    scene = setup_scene_with_view(mock_parser_host)
    a = scene.atom_items[scene.create_atom("C", QPointF(0, 0))]
    b = scene.atom_items[scene.create_atom("C", QPointF(50, 0))]
    scene.create_bond(a, b, bond_order=1, bond_stereo=1)
    bond = scene.find_bond_between(a, b)

    event = MagicMock()
    event.key.return_value = Qt.Key.Key_W
    event.modifiers.return_value = Qt.KeyboardModifier.NoModifier
    with (
        patch.object(MoleculeScene, "find_atom_near", return_value=None),
        patch.object(MoleculeScene, "itemAt", return_value=bond),
    ):
        scene.keyPressEvent(event)
    assert (b.atom_id, a.atom_id) in scene.data.bonds
    assert _keys_in_sync(scene)

    # Plain bond again ("1"): data sorts the key; bond_items must follow.
    event.key.return_value = Qt.Key.Key_1
    with (
        patch.object(MoleculeScene, "find_atom_near", return_value=None),
        patch.object(MoleculeScene, "itemAt", return_value=bond),
    ):
        scene.keyPressEvent(event)
    assert _keys_in_sync(scene)


@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_right_click_clears_ez_label_whatever_the_key_direction(
    mock_drag, mock_parser_host
):
    """Clearing an E/Z label finds the data key without trusting bond_items keys."""
    scene = setup_scene_with_view(mock_parser_host)
    a = scene.atom_items[scene.create_atom("C", QPointF(0, 0))]
    b = scene.atom_items[scene.create_atom("C", QPointF(50, 0))]
    scene.create_bond(a, b, bond_order=2, bond_stereo=0)
    bond = scene.find_bond_between(a, b)
    key = next(iter(scene.data.bonds))
    scene.data.bonds[key]["stereo"] = 3
    bond.stereo = 3
    # A stale reversed key, as older flips left behind, listed first.
    scene.bond_items = {key[::-1]: bond, **scene.bond_items}
    scene.mode = "bond_2_5"

    with patch.object(MoleculeScene, "itemAt", return_value=bond):
        scene.mousePressEvent(
            create_mock_event(QPointF(25, 0), Qt.MouseButton.RightButton)
        )

    assert scene.data.bonds[key]["stereo"] == 0


def test_purge_on_quit_is_connected_once(mock_parser_host):
    """reinitialize_items runs on every clear but must not stack connections."""
    from PyQt6.QtWidgets import QApplication

    scene = setup_scene_with_view(mock_parser_host)
    app = QApplication.instance()
    with patch.object(type(app), "aboutToQuit") as signal:
        scene._purge_connected = False
        scene.reinitialize_items()
        scene.reinitialize_items()
        scene.reinitialize_items()
    assert signal.connect.call_count == 1
