import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QKeyEvent
from moleditpy.ui.molecular_scene_handler import KeyboardMixin
from moleditpy.ui.atom_item import AtomItem

def MockAtom(atom_id, radical=0):
    atom = MagicMock(spec=AtomItem)
    atom.atom_id = atom_id
    atom.radical = radical
    atom.prepareGeometryChange = MagicMock()
    atom.update_style = MagicMock()
    atom.pos = MagicMock(return_value=QPointF(0,0))
    return atom

class MockScene(KeyboardMixin):
    def __init__(self):
        self.window = MagicMock()
        self.window.ui_manager.is_2d_editable = True
        self.data = MagicMock()
        self.data.atoms = {}
        self._selected_items = []
        self._item_at_cursor = None
        self.views = MagicMock(return_value=[MagicMock()])
        self.update_all_items = MagicMock()

    def selectedItems(self):
        return self._selected_items

    def itemAt(self, pos, transformCase):
        return self._item_at_cursor

@pytest.fixture
def scene():
    return MockScene()

def test_radical_toggle_selected(scene):
    """Test toggling radical on multiple selected atoms."""
    a1 = MockAtom(1, 0)
    a2 = MockAtom(2, 1)
    scene._selected_items = [a1, a2]
    scene.data.atoms = {1: {"radical": 0}, 2: {"radical": 1}}
    
    # Create Period key event
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_Period
    
    # We need to patch QCursor.pos() inside KeyboardMixin.keyPressEvent
    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0,0)
        scene.keyPressEvent(event)
    
    assert a1.radical == 1
    assert a2.radical == 2
    assert scene.data.atoms[1]["radical"] == 1
    assert scene.data.atoms[2]["radical"] == 2
    a1.update_style.assert_called_once()
    a2.update_style.assert_called_once()
    scene.window.edit_actions_manager.push_undo_state.assert_called_once()
    event.accept.assert_called_once()

def test_radical_toggle_at_cursor(scene):
    """Test toggling radical on an atom at the cursor when nothing is selected."""
    a1 = MockAtom(1, 2)
    scene._selected_items = []
    scene._item_at_cursor = a1
    scene.data.atoms = {1: {"radical": 2}}
    
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_Period
    
    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0,0)
        scene.keyPressEvent(event)
    
    assert a1.radical == 0
    assert scene.data.atoms[1]["radical"] == 0
    scene.window.edit_actions_manager.push_undo_state.assert_called_once()

def test_radical_toggle_no_target(scene):
    """Test that nothing happens if no atoms are selected or at the cursor."""
    scene._selected_items = []
    scene._item_at_cursor = None
    
    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_Period
    
    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0,0)
        scene.keyPressEvent(event)
    
    scene.window.edit_actions_manager.push_undo_state.assert_not_called()
    event.accept.assert_not_called()
