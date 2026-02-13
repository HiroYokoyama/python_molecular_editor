import pytest
from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtGui import QPainter, QTransform
from PyQt6.QtWidgets import QGraphicsScene, QStyleOptionGraphicsItem
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
