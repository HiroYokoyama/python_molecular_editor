import pytest
import math
from PyQt6.QtCore import Qt, QPointF, QLineF
from PyQt6.QtGui import QKeyEvent, QTransform, QMouseEvent
from PyQt6.QtWidgets import QApplication, QGraphicsLineItem
from moleditpy.modules.molecule_scene import MoleculeScene
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from unittest.mock import MagicMock, patch

def setup_scene_with_view(mock_parser_host):
    scene = MoleculeScene(mock_parser_host.data, mock_parser_host)
    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()
    mock_view.mapToScene.return_value = QPointF(0, 0)
    mock_view.mapFromGlobal.return_value = QPointF(0, 0)
    scene.views = MagicMock(return_value=[mock_view])
    scene.template_preview = MagicMock()
    scene.window.is_2d_editable = True
    scene.mode = 'select'
    return scene

def create_mock_mouse_event(pos, button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.scenePos.return_value = pos
    event.button.return_value = button
    return event

def test_scene_keypress_modes(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    def mock_set_mode(mode): scene.mode = mode
    scene.window.set_mode_and_update_toolbar.side_effect = mock_set_mode

    # Test space -> activate_select_mode
    scene.mode = 'atom_C'
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Space, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert scene.window.activate_select_mode.called

    key_map = {
        Qt.Key.Key_C: 'atom_C',
        Qt.Key.Key_H: 'atom_H',
        Qt.Key.Key_O: 'atom_O',
        Qt.Key.Key_N: 'atom_N',
        Qt.Key.Key_S: 'atom_S',
        Qt.Key.Key_P: 'atom_P',
        Qt.Key.Key_F: 'atom_F',
        Qt.Key.Key_B: 'atom_B',
        Qt.Key.Key_I: 'atom_I'
    }
    
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        for key, expected_mode in key_map.items():
            event = QKeyEvent(QKeyEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier)
            scene.keyPressEvent(event)
            assert scene.mode == expected_mode

def test_scene_keypress_special_symbols(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    def mock_set_mode(mode): scene.mode = mode
    scene.window.set_mode_and_update_toolbar.side_effect = mock_set_mode
    
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_C, Qt.KeyboardModifier.ShiftModifier)
        scene.keyPressEvent(event)
        assert scene.mode == 'atom_Cl'
        
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_B, Qt.KeyboardModifier.ShiftModifier)
        scene.keyPressEvent(event)
        assert scene.mode == 'atom_Br'
        
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_S, Qt.KeyboardModifier.ShiftModifier)
        scene.keyPressEvent(event)
        assert scene.mode == 'atom_Si'

def test_scene_keypress_delete(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    atom_item.setSelected(True)
    
    with patch.object(MoleculeScene, 'delete_items') as mock_delete, \
         patch.object(MoleculeScene, 'selectedItems', return_value=[atom_item]), \
         patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
            event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier)
            scene.keyPressEvent(event)
            mock_delete.assert_called()

def test_scene_maintenance_methods(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    aid = scene.create_atom("C", QPointF(0,0))
    # mock_parser_host's create_atom should populate data.atoms
    atom_item = mock_parser_host.data.atoms[aid]['item']
    
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
    a1.pos.return_value = QPointF(0,0)
    a1.scenePos.return_value = QPointF(0,0)
    a1.__class__ = AtomItem
    a2.pos.return_value = QPointF(50,0)
    a2.scenePos.return_value = QPointF(50,0)
    a2.__class__ = AtomItem
    
    with patch.object(scene, 'items', return_value=[a1, a2]):
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
    scene.data.bonds[(10, 20)] = {'stereo': 0, 'item': bond}
    
    scene.update_bond_stereo(bond, 1)
    assert scene.data.bonds[(10, 20)]['stereo'] == 1
    assert bond.set_stereo.called

def test_scene_mouse_drag_create_bond_existing_atoms(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'atom_C'
    scene.current_atom_symbol = 'C'
    
    aid1 = scene.create_atom("C", QPointF(0, 0))
    aid2 = scene.create_atom("C", QPointF(50, 0))
    a1 = mock_parser_host.data.atoms[aid1]['item']
    a2 = mock_parser_host.data.atoms[aid2]['item']
    
    # 1. Press on A1
    press_event = create_mock_mouse_event(QPointF(0, 0))
    with patch.object(scene, 'itemAt', return_value=a1):
        scene.mousePressEvent(press_event)
    
    # 2. Release on A2
    release_event = create_mock_mouse_event(QPointF(50, 0))
    with patch('PyQt6.QtWidgets.QApplication.startDragDistance', return_value=10), \
         patch.object(scene, 'itemAt', return_value=a2):
        scene.mouseReleaseEvent(release_event)
    
    # Verify bond creation via mock_parser_host's data.add_bond
    assert (aid1, aid2) in scene.data.bonds or (aid2, aid1) in scene.data.bonds

def test_scene_mouse_click_create_single_atom(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'atom_O'
    scene.current_atom_symbol = 'O'
    
    # Click in blank space
    click_event = create_mock_mouse_event(QPointF(100, 100))
    with patch.object(scene, 'itemAt', return_value=None), \
         patch('PyQt6.QtWidgets.QApplication.startDragDistance', return_value=10):
        # mousePressEvent sets self.start_pos and self.press_pos
        scene.mousePressEvent(click_event)
        # mouseReleaseEvent calculates is_click and triggers atom creation if start_pos matches current pos
        scene.mouseReleaseEvent(click_event)
    
    # Verify an 'O' atom was created
    assert any(a['symbol'] == 'O' for a in scene.data.atoms.values())

def test_scene_right_click_delete(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    aid = scene.create_atom("C", QPointF(0, 0))
    a1 = mock_parser_host.data.atoms[aid]['item']
    
    event = create_mock_mouse_event(QPointF(0, 0), button=Qt.MouseButton.RightButton)
    with patch.object(scene, 'itemAt', return_value=a1), \
         patch.object(MoleculeScene, 'delete_items') as mock_del:
        scene.mousePressEvent(event)
        mock_del.assert_called_once_with({a1})

def test_scene_bond_complex_interactions(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    
    # Create bond
    aid1 = scene.create_atom("C", QPointF(0, 0))
    aid2 = scene.create_atom("C", QPointF(50, 0))
    a1 = mock_parser_host.data.atoms[aid1]['item']
    a2 = mock_parser_host.data.atoms[aid2]['item']
    scene.create_bond(a1, a2, bond_order=2)
    bond_key = (aid1, aid2) if (aid1, aid2) in scene.data.bonds else (aid2, aid1)
    bond_item = scene.data.bonds[bond_key]['item']
    
    # Test Key_E / Key_Z on hovered bond
    scene.hovered_item = bond_item
    
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(25,0)):
        
        # Test Key_E (Stereo 4)
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_E, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert scene.data.bonds[bond_key]['stereo'] == 4
        
        # Test Key_Z (Stereo 3)
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Z, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert scene.data.bonds[bond_key]['stereo'] == 3

    # Test Key_W / Key_D on selected bond
    scene.hovered_item = None
    bond_item.setSelected(True)
    
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        # Key_W (Stereo 1 - Up)
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_W, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        # Bond order should become 1, stereo 1
        current_bond_data = scene.data.bonds.get(bond_key) or scene.data.bonds.get((aid2, aid1))
        assert current_bond_data['order'] == 1
        assert current_bond_data['stereo'] == 1

def test_scene_atom_keyboard_properties(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    atom_item.setSelected(True)
    
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        
        # Charge +
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Plus, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert atom_item.charge == 1
        
        # Charge -
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Minus, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert atom_item.charge == 0
        
        # Radical .
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Period, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        assert atom_item.radical == 1

def test_scene_template_interaction(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    scene.mode = 'select'
    
    # Test Key_4 on blank space -> triggers mode change
    with patch.object(MoleculeScene, 'itemAt', return_value=None), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(0,0)):
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_4, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        scene.window.set_mode_and_update_toolbar.assert_called_with('template_benzene')

    # Test Key_4 on existing atom -> triggers add_molecule_fragment
    aid = scene.create_atom("C", QPointF(0, 0))
    atom_item = mock_parser_host.data.atoms[aid]['item']
    
    with patch.object(MoleculeScene, 'itemAt', return_value=atom_item), \
         patch('PyQt6.QtGui.QCursor.pos', return_value=QPointF(20,0)), \
         patch.object(scene, 'add_molecule_fragment') as mock_add_frag:
        event = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_4, Qt.KeyboardModifier.NoModifier)
        scene.keyPressEvent(event)
        mock_add_frag.assert_called()

def test_scene_user_template_logic(mock_parser_host):
    scene = setup_scene_with_view(mock_parser_host)
    
    # Mock template data
    scene.user_template_data = {
        'atoms': [{'id': 'a1', 'symbol': 'C', 'x': 0, 'y': 0}],
        'bonds': []
    }
    
    # Test update_user_template_preview
    scene.update_user_template_preview(QPointF(10, 10))
    scene.template_preview.set_user_template_geometry.assert_called()
    scene.template_preview.show.assert_called()
    
    # Test add_user_template_fragment
    context = {
        'points': [QPointF(0, 0)],
        'atoms_data': [{'id': 'a1', 'symbol': 'N', 'charge': 1}],
        'bonds_info': []
    }
    scene.add_user_template_fragment(context)
    # Verify N+ atom created
    assert any(a['symbol'] == 'N' and a['charge'] == 1 for a in scene.data.atoms.values())

