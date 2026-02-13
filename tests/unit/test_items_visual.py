import pytest
from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtGui import QPainter, QTransform, QPainterPath, QBrush
from PyQt6.QtWidgets import QGraphicsScene, QStyleOptionGraphicsItem, QGraphicsItem
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from unittest.mock import MagicMock, patch
from rdkit import Chem

def test_atom_item_visual_states(mock_parser_host):
    """Test AtomItem paint and boundingRect with different states."""
    atom = AtomItem(0, "C", QPointF(0, 0))
    atom.charge = 1
    atom.radical = 1
    atom.implicit_h_count = 1
    
    # Mock painter and option
    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()
    
    # Test boundingRect handles charges and radicals
    rect = atom.boundingRect()
    assert isinstance(rect, QRectF)
    assert not rect.isEmpty()
    
    # Test paint
    atom.paint(painter, option, None)
    assert painter.drawText.called
    assert painter.drawEllipse.called # Radical dot

def test_atom_item_h_label_flip(mock_parser_host):
    """Test AtomItem H-label flipping logic based on bond direction."""
    atom_center = AtomItem(0, "N", QPointF(100, 100))
    atom_center.implicit_h_count = 1
    
    # Partner atom to the right
    atom_right = AtomItem(1, "C", QPointF(200, 100))
    
    bond = MagicMock()
    bond.atom1 = atom_center
    bond.atom2 = atom_right
    atom_center.bonds = [bond]
    
    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()
    
    # With partner to the right, H should be on the left (flipped)
    # The logic in paint() checks total_dx = partner.x - my.x
    # 200 - 100 = 100 > 0 -> flip_text = True
    atom_center.paint(painter, option, None)
    
    # Verify alignment or text order if possible, but mainly ensure it runs
    assert painter.drawText.called

def test_bond_item_render_complex_types(mock_parser_host):
    """Test BondItem rendering for triple bonds, wedge/dash, and E/Z labels."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    
    # Triple bond
    bond = BondItem(a1, a2, order=3)
    painter = MagicMock(spec=QPainter)
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    # Triple bond draws 3 lines
    assert painter.drawLine.call_count == 3
    
    # Wedge
    painter.reset_mock()
    bond.order = 1
    bond.stereo = 1 # Wedge
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    assert painter.drawPolygon.called
    
    # Dash
    painter.reset_mock()
    bond.stereo = 2 # Dash
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    # Dash draws multiple lines
    assert painter.drawLine.call_count > 1
    
    # E/Z Label
    painter.reset_mock()
    bond.order = 2
    bond.stereo = 3 # Z
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    assert painter.drawPath.called

def test_bond_item_ring_logic(mock_parser_host):
    """Test BondItem ring rendering logic by mocking RDKit mol integration."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    bond = BondItem(a1, a2, order=2)
    
    # Set up mock scene/window environment
    scene = MagicMock(spec=QGraphicsScene)
    view = MagicMock()
    window = MagicMock()
    scene.views.return_value = [view]
    view.window.return_value = window
    bond.scene = MagicMock(return_value=scene)
    
    # Create a mock RDKit mol with a ring
    mol = Chem.MolFromSmiles("C1=CC=CC=C1")
    # Need to set _original_atom_id for the mapping logic
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i)
    
    # Mock to_rdkit_mol
    window.data.to_rdkit_mol.return_value = mol
    
    painter = MagicMock(spec=QPainter)
    # The paint logic tries to find atoms with atom_id 0 and 1
    a1.atom_id = 0
    a2.atom_id = 1
    
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    
    # In ring, it draws a center line and a shortened inner line (2 lines total)
    assert painter.drawLine.call_count == 2

def test_atom_item_item_change_updates_bonds(mock_parser_host):
    """Test that moving an atom updates connected bonds."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "O", QPointF(10, 0))
    
    # Mock bond
    bond = MagicMock(spec=BondItem)
    bond.scene.return_value = True # ensure scene() check passes
    a1.bonds = [bond]
    
    # Trigger itemChange by moving
    # itemChange is protected, usually called by C++ or internal logic.
    # We can call it directly or trigger via setPos if we were in a full scene (but headless sometimes tricky).
    # Calling itemChange directly is robust for unit testing the logic.
    a1.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable, True)
    a1.itemChange(QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged, QPointF(5, 5))
    
    bond.update_position.assert_called()

def test_atom_item_paint_transparent_bg(mock_parser_host):
    """Test AtomItem.paint with transparent background logic."""
    atom = AtomItem(0, "C", QPointF(0, 0))
    atom.is_visible = True
    
    scene = MagicMock(spec=QGraphicsScene)
    # Mock background brush to be NoBrush (transparent)
    scene.backgroundBrush.return_value = QBrush(Qt.BrushStyle.NoBrush)
    
    # Mock helper to allow .scene() call
    atom.scene = MagicMock(return_value=scene)
    
    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()
    
    atom.paint(painter, option, None)
    
    # Needs to verify it set composition mode to Clear
    painter.setCompositionMode.assert_called_with(QPainter.CompositionMode.CompositionMode_Clear)

def test_atom_item_paint_resilience_to_deleted_bond(mock_parser_host):
    """Test paint doesn't crash if a bond refers to a deleted atom."""
    atom = AtomItem(0, "C", QPointF(0, 0))
    atom.is_visible = True
    
    # Create a bond where the PARTNER is "deleted"
    bond = MagicMock()
    bond.atom1 = atom
    # atom2 is the partner
    bond.atom2 = MagicMock()
    
    atom.bonds = [bond]
    
    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()
    
    with patch('moleditpy.modules.atom_item.sip_isdeleted_safe', return_value=True):
        # Should not crash
        atom.paint(painter, option, None)
        # Should still draw text
        assert painter.drawText.called

def test_atom_item_shape_collision(mock_parser_host):
    """Test atom shape returns a path for collision detection."""
    atom = AtomItem(0, "C", QPointF(0, 0))
    
    # Case 1: No scene/view (fallback)
    path = atom.shape()
    assert isinstance(path, QPainterPath)
    assert not path.isEmpty()
    
    # Case 2: With scene/view (scaled)
    scene = MagicMock(spec=QGraphicsScene)
    view = MagicMock()
    view.transform.return_value.m11.return_value = 2.0 # Scale
    scene.views.return_value = [view]
    atom.scene = MagicMock(return_value=scene)
    
    path_scaled = atom.shape()
    assert not path_scaled.isEmpty()

def test_bond_item_shape_stroked(mock_parser_host):
    """Test bond shape is a stroked path (wider than line)."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(10, 0))
    bond = BondItem(a1, a2)
    
    # Set mocked scene for settings (though shape might not use settings directly, it uses width)
    # BondItem.shape usually creates a path and uses QPainterPathStroker
    
    path = bond.shape()
    assert isinstance(path, QPainterPath)
    # A simple line has 0 area, but a stroked shape should have area
    assert not path.isEmpty()
    rect = path.boundingRect()
    assert rect.width() > 0
    assert rect.height() > 0 # Should have thickness
