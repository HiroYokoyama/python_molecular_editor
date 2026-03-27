import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QKeyEvent
from moleditpy.modules.molecular_scene_handler import KeyboardMixin
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem

class MockScene(KeyboardMixin):
    def __init__(self):
        self.window = MagicMock()
        self.window.is_2d_editable = True
        self.data = MagicMock()
        self._selected_items = []
        self._item_at_cursor = None
        self.view_mock = MagicMock()
        self.view_mock.mapToScene.return_value = QPointF(0,0)
        self.view_mock.mapFromGlobal.return_value = QPointF(0,0)
        self.view_mock.transform.return_value = MagicMock()
        self.views = MagicMock(return_value=[self.view_mock])
        self.update_all_items = MagicMock()
        # Mocking calculation method that is used when Key_4 is pressed
        self._calculate_polygon_from_edge = MagicMock(return_value=[QPointF(0,0)]*6)
        self.add_molecule_fragment = MagicMock()

    def selectedItems(self):
        return self._selected_items

    def itemAt(self, pos, transform):
        return self._item_at_cursor

@pytest.fixture
def scene():
    return MockScene()

def test_benzene_shortcut_on_atom(scene):
    """Test one-shot benzene placement when cursor is over an atom."""
    atom = MagicMock(spec=AtomItem)
    atom.pos.return_value = QPointF(100, 100)
    scene._item_at_cursor = atom
    
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_4
    
    # Needs to match the view's mapToScene etc. logic
    # In MockScene, itemAt returns self._item_at_cursor regardless of pos argument
    with patch("moleditpy.modules.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(150, 100) # Simulated global cursor pos
        scene.keyPressEvent(event)
    
    scene._calculate_polygon_from_edge.assert_called_once()
    scene.add_molecule_fragment.assert_called_once()
    scene.window.push_undo_state.assert_called_once()
    event.accept.assert_called_once()

def test_benzene_shortcut_on_bond(scene):
    """Test one-shot benzene placement when cursor is over a bond."""
    bond = MagicMock(spec=BondItem)
    bond.atom1 = MagicMock(spec=AtomItem)
    bond.atom2 = MagicMock(spec=AtomItem)
    bond.atom1.pos.return_value = QPointF(100, 100)
    bond.atom2.pos.return_value = QPointF(100, 120)
    scene._item_at_cursor = bond
    
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_4
    
    with patch("moleditpy.modules.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(120, 110)
        print(f"DEBUG: item_at_cursor type: {type(scene._item_at_cursor)}")
        scene.keyPressEvent(event)
    
    # Check that use_existing_length=True was passed for bond placement
    # We use ANY for variable geometric args
    called_args, called_kwargs = scene._calculate_polygon_from_edge.call_args
    assert called_kwargs.get("use_existing_length") is True
    scene.add_molecule_fragment.assert_called_once()

def test_benzene_shortcut_empty_space(scene):
    """Test mode switch when cursor is over empty space."""
    scene._item_at_cursor = None
    
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_4
    
    with patch("moleditpy.modules.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(200, 200)
        scene.keyPressEvent(event)
    
    scene.window.set_mode_and_update_toolbar.assert_called_with("template_benzene")
    event.accept.assert_called_once()
    scene.add_molecule_fragment.assert_not_called()
