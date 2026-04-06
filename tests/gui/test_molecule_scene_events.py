# -*- coding: utf-8 -*-
import pytest
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QKeyEvent
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem

def test_bond_stereo_toggle_keys(window, qtbot):
    """Test Z and E keys toggle double bond stereochemistry."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    
    # Create two atoms and a double bond
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a2_id = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(data.atoms[a1_id]["item"], data.atoms[a2_id]["item"], bond_order=2)
    
    bond_key = (min(a1_id, a2_id), max(a1_id, a2_id))
    bond_item = data.bonds[bond_key]["item"]
    
    # Simulate hover
    scene.hovered_item = bond_item
    
    # Press 'Z'
    event_z = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Z, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_z)
    assert data.bonds[bond_key]["stereo"] == 3
    assert bond_item.stereo == 3
    
    # Press 'E'
    event_e = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_E, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_e)
    assert data.bonds[bond_key]["stereo"] == 4
    assert bond_item.stereo == 4

def test_atom_addition_keys(window, qtbot):
    """Test 1, 2, 3 keys add atoms/bonds from selected atom."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    
    # Place initial atom and select it
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a1_item = data.atoms[a1_id]["item"]
    a1_item.setSelected(True)
    
    # Press '1' to add a single bond C
    event_1 = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_1, Qt.KeyboardModifier.NoModifier)
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
    new_atom_item = data.atoms[new_atom_id]["item"]
    scene.clearSelection()
    new_atom_item.setSelected(True)
    
    # Press '2' to add a double bond C
    event_2 = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_2, Qt.KeyboardModifier.NoModifier)
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
    a1_item = data.atoms[a1_id]["item"]
    a2_item = data.atoms[a2_id]["item"]
    
    a1_item.setSelected(True)
    
    # Press 'Delete'
    event_del = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_del)
    
    assert a1_id not in data.atoms
    assert a2_id in data.atoms
    
    # Select a2 and press 'Backspace'
    a2_item.setSelected(True)
    event_back = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Backspace, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_back)
    
    assert a2_id not in data.atoms
    assert len(data.atoms) == 0

def test_temp_line_cancellation(window, qtbot):
    """Test Delete key cancels an active temp_line (bond drawing)."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a1_item = data.atoms[a1_id]["item"]
    
    # Simulate start of bond drawing
    scene.start_atom = a1_item
    scene.start_pos = a1_item.pos()
    from PyQt6.QtWidgets import QGraphicsLineItem
    scene.temp_line = QGraphicsLineItem()
    scene.addItem(scene.temp_line)
    
    # Press 'Delete'
    event_del = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_Delete, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_del)
    
    assert scene.temp_line is None

def test_bonding_to_existing_atom(window, qtbot):
    """Test that pressing 1, 2, 3 bonds to an existing atom if it's nearby."""
    scene = window.init_manager.scene
    data = window.state_manager.data
    from moleditpy.utils.constants import DEFAULT_BOND_LENGTH
    
    # Create two atoms separated by DEFAULT_BOND_LENGTH
    a1_id = scene.create_atom("C", QPointF(0, 0))
    a2_id = scene.create_atom("C", QPointF(0, -DEFAULT_BOND_LENGTH))
    
    a1_item = data.atoms[a1_id]["item"]
    a1_item.setSelected(True)
    
    # Press '1'. Since (0, -L) is exactly where a2 is, it should bond to a2
    # instead of creating a new atom.
    event_1 = QKeyEvent(QKeyEvent.Type.KeyPress, Qt.Key.Key_1, Qt.KeyboardModifier.NoModifier)
    scene.keyPressEvent(event_1)
    
    assert len(data.atoms) == 2 # No new atom created
    assert len(data.bonds) == 1 # Bond created between a1 and a2
    
    bond_key = (min(a1_id, a2_id), max(a1_id, a2_id))
    assert bond_key in data.bonds
