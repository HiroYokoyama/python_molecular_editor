import pytest
from unittest.mock import MagicMock
from PyQt6.QtCore import QPointF
from moleditpy.ui.molecular_scene_handler import KeyboardMixin
from moleditpy.ui.atom_item import AtomItem


def MockAtom(pos=QPointF(0, 0), symbol="C"):
    atom = MagicMock(spec=AtomItem)
    atom.pos = MagicMock(return_value=pos)
    atom.symbol = symbol
    atom.bonds = []
    return atom


class MockBond:
    def __init__(self, a1, a2):
        self.atom1 = a1
        self.atom2 = a2


class MockScene(KeyboardMixin):
    def __init__(self):
        pass


@pytest.fixture
def scene():
    return MockScene()


L = 20.0  # Default bond length reference


def test_placement_0_neighbors(scene):
    """Zero bonds: should place default (up)."""
    start_atom = MockAtom(QPointF(100, 100))
    offset = scene._calculate_new_atom_position(start_atom, L)
    assert offset.x() == pytest.approx(0)
    assert offset.y() == pytest.approx(-L)


def test_placement_1_neighbor(scene):
    """One bond: should rotate 60 degrees clockwise from existing bond vector."""
    start_atom = MockAtom(QPointF(100, 100))
    other_atom = MockAtom(
        QPointF(80, 100)
    )  # Neighbor at Left (existing vector is Right 20, 0)
    bond = MockBond(start_atom, other_atom)
    start_atom.bonds = [bond]

    offset = scene._calculate_new_atom_position(start_atom, L)
    # Existing vector (start - other): (20, 0)
    # 60 deg CW rotation: (20*cos60 - 0*sin60, 20*sin60 + 0*cos60) = (10, 17.32)
    # Then scaled to L=20: (10, 17.32)
    assert offset.x() == pytest.approx(10.0)
    assert offset.y() == pytest.approx(17.32, abs=0.01)


def test_placement_3_neighbors_balanced(scene):
    """Three balanced neighbors: sum is near zero, should use fallback (45 deg offset)."""
    start_atom = MockAtom(QPointF(0, 0))
    # Neighbors at ~trigonal planar positions
    a1 = MockAtom(QPointF(1, 0))
    a2 = MockAtom(QPointF(-0.5, 0.866))
    a3 = MockAtom(QPointF(-0.5, -0.866))

    start_atom.bonds = [
        MockBond(start_atom, a1),
        MockBond(start_atom, a2),
        MockBond(start_atom, a3),
    ]

    offset = scene._calculate_new_atom_position(start_atom, L)
    # Sum is ~0. Fallback is (L*0.7071, -L*0.7071)
    assert offset.x() == pytest.approx(L * 0.7071)
    assert offset.y() == pytest.approx(-L * 0.7071)


def test_placement_3_neighbors_unbalanced(scene):
    """Three unbalanced neighbors: should place in direction opposite to the sum."""
    start_atom = MockAtom(QPointF(0, 0))
    # Neighbors at Right, Bottom, Left
    a1 = MockAtom(QPointF(1, 0))
    a2 = MockAtom(QPointF(0, 1))
    a3 = MockAtom(QPointF(-1, 0))

    start_atom.bonds = [
        MockBond(start_atom, a1),
        MockBond(start_atom, a2),
        MockBond(start_atom, a3),
    ]

    offset = scene._calculate_new_atom_position(start_atom, L)
    # Vector sum: (1,0) + (0,1) + (-1,0) = (0, 1)
    # Opposite: (0, -1)
    # Scaled to L: (0, -L)
    assert offset.x() == pytest.approx(0)
    assert offset.y() == pytest.approx(-L)


def test_placement_2_neighbors_skeleton(scene):
    """Two neighbors: should continue skeleton (opposite to average bond vector)."""
    start_atom = MockAtom(QPointF(0, 0))
    # "V" shape at bottom
    a1 = MockAtom(QPointF(-0.866, 0.5))
    a2 = MockAtom(QPointF(0.866, 0.5))

    start_atom.bonds = [MockBond(start_atom, a1), MockBond(start_atom, a2)]

    offset = scene._calculate_new_atom_position(start_atom, L)
    # Unit vectors to neighbors: (-0.866, 0.5) and (0.866, 0.5)
    # Sum: (0, 1.0). Opposite: (0, -1.0). Scaled: (0, -L)
    assert offset.x() == pytest.approx(0)
    assert offset.y() == pytest.approx(-L)
