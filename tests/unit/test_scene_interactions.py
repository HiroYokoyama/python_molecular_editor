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

@patch('PyQt6.QtWidgets.QGraphicsScene.mousePressEvent')
@patch('PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent')
@patch('moleditpy.modules.molecule_scene.QApplication.startDragDistance', return_value=1000)
def test_scene_toggle_radical(mock_drag, mock_release, mock_press, mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'radical'
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    atom_item.radical = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()
    
    with patch('moleditpy.modules.molecule_scene.isinstance', return_value=True), \
         patch.object(MoleculeScene, 'itemAt', return_value=atom_item):
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

@patch('PyQt6.QtWidgets.QGraphicsScene.mousePressEvent')
@patch('PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent')
@patch('moleditpy.modules.molecule_scene.QApplication.startDragDistance', return_value=1000)
def test_scene_toggle_charge(mock_drag, mock_release, mock_press, mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'charge_plus'
    aid = scene.create_atom("C", QPointF(100, 100))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    atom_item.charge = 0
    atom_item.prepareGeometryChange = MagicMock()
    atom_item.update_style = MagicMock()
    
    with patch('moleditpy.modules.molecule_scene.isinstance', return_value=True), \
         patch.object(MoleculeScene, 'itemAt', return_value=atom_item):
        event = create_mock_event(QPointF(100, 100))
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.charge == 1
        scene.mode = 'charge_minus'
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
        assert atom_item.charge == 0

def test_add_benzene_fragment(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    points = [QPointF(20 * math.cos(math.radians(i * 60)), 20 * math.sin(math.radians(i * 60))) for i in range(6)]
    bonds_info = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
    scene.add_molecule_fragment(points, bonds_info)
    assert len(scene.data.atoms) == 6
    assert len(scene.data.bonds) == 6
    assert len([b for b in scene.data.bonds.values() if b['order'] == 2]) == 3

def test_benzene_fusion_rotation(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    id1 = scene.create_atom("C", QPointF(-10, 0))
    id2 = scene.create_atom("C", QPointF(10, 0))
    scene.create_bond(scene.data.atoms[id1]['item'], scene.data.atoms[id2]['item'], bond_order=2)
    points = [QPointF(10, 0), QPointF(-10, 0)] + [QPointF(20 * math.cos(math.radians(i * 60)), 20 * math.sin(math.radians(i * 60))) for i in range(2, 6)]
    bonds_info = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
    scene.add_molecule_fragment(points, bonds_info, existing_items=[scene.data.atoms[id2]['item'], scene.data.atoms[id1]['item']])
    assert len(scene.data.atoms) == 6
    eb = scene.find_bond_between(scene.data.atoms[id1]['item'], scene.data.atoms[id2]['item'])
    assert eb.order == 2

def test_delete_selected_items(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    atom_item.setSelected(True)
    with patch.object(MoleculeScene, 'delete_items', return_value=True) as mock_delete:
        with patch.object(MoleculeScene, 'selectedItems', return_value=[atom_item]):
            with patch.object(MoleculeScene, 'itemAt', return_value=atom_item):
                event = MagicMock()
                event.button.return_value = Qt.MouseButton.RightButton
                event.scenePos.return_value = QPointF(0, 0)
                scene.mousePressEvent(event)
                mock_delete.assert_called()

@patch('PyQt6.QtWidgets.QGraphicsScene.mouseDoubleClickEvent')
def test_double_click_select_component(mock_dbl, mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    id1 = scene.create_atom("C", QPointF(0, 0))
    id2 = scene.create_atom("C", QPointF(20, 0))
    id3 = scene.create_atom("C", QPointF(40, 0))
    scene.create_bond(scene.data.atoms[id1]['item'], scene.data.atoms[id2]['item'])
    scene.create_bond(scene.data.atoms[id2]['item'], scene.data.atoms[id3]['item'])
    id_iso = scene.create_atom("O", QPointF(100, 100))
    atom_iso = scene.data.atoms[id_iso]['item']
    for atom_data in scene.data.atoms.values():
        it = atom_data['item']
        it.__class__ = AtomItem
        it.setSelected = MagicMock()
    with patch.object(MoleculeScene, 'itemAt', return_value=scene.data.atoms[id1]['item']):
        _orig_isinstance = isinstance
        def mock_isinstance(obj, types):
            tpl = types if _orig_isinstance(types, tuple) else (types,)
            if hasattr(obj, 'atom_id') and AtomItem in tpl: return True
            if hasattr(obj, 'atom1') and BondItem in tpl: return True
            return _orig_isinstance(obj, types)
        with patch('moleditpy.modules.molecule_scene.isinstance', side_effect=mock_isinstance):
            event = create_mock_event(QPointF(0, 0))
            scene.mouseDoubleClickEvent(event)
            assert scene.data.atoms[id1]['item'].setSelected.called
            assert scene.data.atoms[id2]['item'].setSelected.called
            assert scene.data.atoms[id3]['item'].setSelected.called
            assert not atom_iso.setSelected.called
