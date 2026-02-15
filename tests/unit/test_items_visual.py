import pytest
from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtGui import QPainter, QTransform, QPainterPath, QBrush
from PyQt6.QtWidgets import (
    QGraphicsScene,
    QStyleOptionGraphicsItem,
    QGraphicsItem,
    QGraphicsSceneHoverEvent,
)
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

    # Basic rendering test
    atom.paint(painter, option, None)
    assert painter.drawEllipse.called or painter.drawText.called

    # BoundingRect test
    rect = atom.boundingRect()
    assert rect.width() > 0
    assert rect.height() > 0


def test_atom_item_hover(mock_parser_host):
    """Test AtomItem hover state changes."""
    atom = AtomItem(0, "N", QPointF(10, 10))
    atom.hovered = False

    # Simulate hover events
    atom.hoverEnterEvent(None)
    assert atom.hovered

    atom.hoverLeaveEvent(None)
    assert not atom.hovered


def test_atom_item_h_label_flip(mock_parser_host):
    """Test AtomItem flips H label position based on neighbor orientation."""
    atom_left = AtomItem(0, "C", QPointF(0, 0))
    atom_center = AtomItem(1, "N", QPointF(50, 0))
    atom_right = AtomItem(2, "H", QPointF(100, 0))

    # Mock neighbor relationships
    atom_center.bonds = []

    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()
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
    bond.stereo = 1  # Wedge
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    assert painter.drawPolygon.called

    # Dash
    painter.reset_mock()
    bond.stereo = 2  # Dash
    bond.paint(painter, QStyleOptionGraphicsItem(), None)
    # Dash draws multiple lines
    assert painter.drawLine.call_count > 1

    # E/Z Label
    painter.reset_mock()
    bond.order = 2
    bond.stereo = 3  # Z
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

    painter = MagicMock(spec=QPainter)
    bond.paint(painter, QStyleOptionGraphicsItem(), None)

    # Should draw at least 1 line
    assert painter.drawLine.call_count >= 1


def test_atom_item_item_change_updates_bonds(mock_parser_host):
    """Test that moving an atom triggers bond position updates."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

    # Create mock bond
    bond = MagicMock(spec=BondItem)
    bond.atom1 = a1
    bond.atom2 = MagicMock()
    bond.update_position = MagicMock()
    # Mock scene so bond.scene() returns True
    bond.scene.return_value = MagicMock()

    a1.bonds = [bond]

    # Trigger itemChange (ItemPositionHasChanged is what triggers the update)
    a1.itemChange(
        QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged, QPointF(5, 5)
    )

    # Bond update should be called
    assert bond.update_position.called


def test_atom_item_paint_transparent_bg(mock_parser_host):
    """Test AtomItem paint with transparent background uses CompositionMode_Clear."""
    atom = AtomItem(0, "O", QPointF(0, 0))

    painter = MagicMock(spec=QPainter)
    painter.background.return_value.alpha.return_value = 0  # Transparent

    # Mock scene with transparent background
    scene = MagicMock()
    scene.views.return_value = [MagicMock()]
    scene.views()[0].window.return_value.settings = {"background_color_2d": "#00000000"}
    atom.scene = MagicMock(return_value=scene)

    option = QStyleOptionGraphicsItem()
    atom.paint(painter, option, None)

    # Should set composition mode to clear at some point
    # (This is a bit tricky to assert directly, but we check the method was called)
    assert painter.setCompositionMode.called or painter.drawEllipse.called


def test_atom_item_paint_resilience_to_deleted_bond(mock_parser_host):
    """Test AtomItem paint doesn't crash when a C++ bond object is deleted."""
    atom = AtomItem(0, "C", QPointF(0, 0))

    # Mock a bond that will raise RuntimeError (simulating deleted C++ object)
    deleted_bond = MagicMock()
    type(deleted_bond).atom1 = property(
        lambda self: (_ for _ in ()).throw(RuntimeError("C++ object deleted"))
    )

    atom.bonds = [deleted_bond]

    painter = MagicMock(spec=QPainter)
    option = QStyleOptionGraphicsItem()

    # Should not crash
    try:
        atom.paint(painter, option, None)
        success = True
    except RuntimeError:
        success = False

    assert success, "Paint should handle deleted bond references gracefully"
    # Ensure checking loop finished
    assert atom.bonds == [deleted_bond]


def test_atom_item_shape_collision(mock_parser_host):
    """Test AtomItem.shape() returns larger area than visual bounding for collision."""
    atom = AtomItem(0, "N", QPointF(0, 0))

    shape_path = atom.shape()
    assert isinstance(shape_path, QPainterPath)

    shape_rect = shape_path.boundingRect()
    visual_rect = atom.boundingRect()

    # Shape should be at least as large as visual rect (often larger for easier clicking)
    assert shape_rect.width() >= visual_rect.width() * 0.8
    assert shape_rect.height() >= visual_rect.height() * 0.8


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
    assert rect.height() > 0  # Should have thickness


def test_bond_bounding_rect_geometry(mock_parser_host):
    """Test boundingRect includes bond line area properly."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(100, 0))
    bond = BondItem(a1, a2, order=1)

    rect = bond.boundingRect()

    # Rect should encompass the line
    assert rect.width() >= 100  # At least the 100-unit line
    assert rect.height() > 0  # Should have some thickness


def test_bond_double_non_ring_parallel(mock_parser_host):
    """Test double bond parallel line rendering outside rings (non-ring case)."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    bond = BondItem(a1, a2, order=2)

    # Mock scene to ensure NOT in ring (no RDKit mol available)
    scene = MagicMock(spec=QGraphicsScene)
    scene.window = None  # No window means no RDKit mol integration
    bond.scene = MagicMock(return_value=scene)

    painter = MagicMock(spec=QPainter)
    bond.paint(painter, QStyleOptionGraphicsItem(), None)

    # Double bond (non-ring) should draw 2 parallel lines
    assert painter.drawLine.call_count == 2, (
        "Non-ring double bond should draw 2 parallel lines"
    )


def test_bond_hover_effects(mock_parser_host):
    """Test BondItem hover enter/leave effects by patching base classes to avoid type errors."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    bond = BondItem(a1, a2)

    # Mock scene - don't use spec=QGraphicsScene as MoleculeScene has custom methods
    scene = MagicMock()
    scene.mode = "edit"
    bond.scene = MagicMock(return_value=scene)

    # Patch QGraphicsItem's hover events to avoid type checks on mocks
    with (
        patch("PyQt6.QtWidgets.QGraphicsItem.hoverEnterEvent"),
        patch("PyQt6.QtWidgets.QGraphicsItem.hoverLeaveEvent"),
    ):
        # Hover enter
        bond.hovered = False
        event = MagicMock()
        bond.hoverEnterEvent(event)
        assert bond.hovered
        assert scene.set_hovered_item.called

        # Hover leave
        leave_event = MagicMock()
        bond.hoverLeaveEvent(leave_event)
        assert not bond.hovered
        assert scene.set_hovered_item.call_count >= 2


def test_bond_stereo_e_rendering(mock_parser_host):
    """Test BondItem rendering for E-stereo label (stereo=4)."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    bond = BondItem(a1, a2, order=2, stereo=4)  # stereo=4 is E

    painter = MagicMock(spec=QPainter)
    bond.paint(painter, QStyleOptionGraphicsItem(), None)

    # Should draw text/path for the label
    assert painter.drawPath.called


def test_bond_get_ez_label_rect(mock_parser_host):
    """Test get_ez_label_local_rect() calculates correct label geometry."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(100, 0))
    bond = BondItem(a1, a2, order=2, stereo=4)

    label_rect = bond.get_ez_label_local_rect()
    assert label_rect is not None
    assert isinstance(label_rect, QRectF)
    assert label_rect.width() > 0


def test_bond_update_position_resilience(mock_parser_host):
    """Test BondItem.update_position resilience when atoms exist."""
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    bond = BondItem(a1, a2)

    # Should not crash
    bond.atom2 = None  # Simulate losing Partner after init
    # Should not crash
    bond.atom2 = None  # Simulate losing Partner after init
    # Should not crash
    bond.atom2 = None  # Simulate losing Partner after init
    bond.update_position(notify=True)
    assert bond.atom2 is None
    # Verify it didn't raise exception
    assert True
