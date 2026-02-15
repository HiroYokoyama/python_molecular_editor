import pytest
import math
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QMouseEvent, QTransform
from moleditpy.modules.molecule_scene import MoleculeScene
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from unittest.mock import MagicMock, patch


def setup_scene_with_view(mock_parser_host):
    scene = MoleculeScene(mock_parser_host.data, mock_parser_host)
    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()
    scene.views = MagicMock(return_value=[mock_view])
    return scene


def create_mock_event(pos=QPointF(100, 100), button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.scenePos.return_value = pos
    event.button.return_value = button
    return event


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch(
    "moleditpy.modules.molecule_scene.QApplication.startDragDistance", return_value=1000
)
def test_scene_toggle_radical(mock_drag, mock_release, mock_press, mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "radical"
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = mock_parser_host.data.atoms[aid]["item"]
    atom_item.radical = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()

    with (
        patch("moleditpy.modules.molecule_scene.isinstance", return_value=True),
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
@patch(
    "moleditpy.modules.molecule_scene.QApplication.startDragDistance", return_value=1000
)
def test_scene_toggle_charge(mock_drag, mock_release, mock_press, mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "charge_plus"
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = mock_parser_host.data.atoms[aid]["item"]
    atom_item.charge = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()

    with (
        patch("moleditpy.modules.molecule_scene.isinstance", return_value=True),
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
    scene = setup_scene_with_view(mock_parser_host)
    id1 = scene.create_atom("C", QPointF(-10, 0))
    id2 = scene.create_atom("C", QPointF(10, 0))
    scene.create_bond(
        scene.data.atoms[id1]["item"], scene.data.atoms[id2]["item"], bond_order=2
    )
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
        existing_items=[scene.data.atoms[id2]["item"], scene.data.atoms[id1]["item"]],
    )
    assert len(scene.data.atoms) == 6
    eb = scene.find_bond_between(
        scene.data.atoms[id1]["item"], scene.data.atoms[id2]["item"]
    )
    assert eb.order == 2


def test_delete_selected_items(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "select"
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]["item"]
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
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = "select"
    id1 = scene.create_atom("C", QPointF(0, 0))
    id2 = scene.create_atom("C", QPointF(20, 0))
    id3 = scene.create_atom("C", QPointF(40, 0))
    scene.create_bond(scene.data.atoms[id1]["item"], scene.data.atoms[id2]["item"])
    scene.create_bond(scene.data.atoms[id2]["item"], scene.data.atoms[id3]["item"])
    id_iso = scene.create_atom("O", QPointF(100, 100))
    atom_iso = scene.data.atoms[id_iso]["item"]
    for atom_data in scene.data.atoms.values():
        it = atom_data["item"]
        it.__class__ = AtomItem
        it.setSelected = MagicMock()
    with patch.object(
        MoleculeScene, "itemAt", return_value=scene.data.atoms[id1]["item"]
    ):
        _orig_isinstance = isinstance

        def mock_isinstance(obj, types):
            tpl = types if _orig_isinstance(types, tuple) else (types,)
            if hasattr(obj, "atom_id") and AtomItem in tpl:
                return True
            if hasattr(obj, "atom1") and BondItem in tpl:
                return True
            return _orig_isinstance(obj, types)

        with patch(
            "moleditpy.modules.molecule_scene.isinstance", side_effect=mock_isinstance
        ):
            event = create_mock_event(QPointF(0, 0))
            scene.mouseDoubleClickEvent(event)
            assert scene.data.atoms[id1]["item"].setSelected.called
            assert scene.data.atoms[id2]["item"].setSelected.called
            assert scene.data.atoms[id3]["item"].setSelected.called
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
        mock_parser_host.set_mode_and_update_toolbar.assert_called_with(
            "template_benzene"
        )

    # Case 2: Key '+' -> Charge increase on hovered atom
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = scene.data.atoms[aid]["item"]
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
    atom_item = scene.data.atoms[aid]["item"]

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
    a1 = scene.data.atoms[id1]["item"]

    # Atom B at 100,100
    id2 = scene.create_atom("C", QPointF(100, 100))
    a2 = scene.data.atoms[id2]["item"]

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
    assert bond is not None
    assert bond.order == 1
