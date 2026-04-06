import pytest
from unittest.mock import MagicMock
from moleditpy.ui.molecular_scene_handler import TemplateMixin

class MockScene(TemplateMixin):
    def __init__(self):
        self.data = MagicMock()
        self.data.atoms = {}
        self.data.bonds = {}

    def find_bond_between(self, atom1, atom2):
        """Helper to find a bond based on atom.bonds lists."""
        for b in atom1.bonds:
            if (b.atom1 is atom1 and b.atom2 is atom2) or (
                b.atom1 is atom2 and b.atom2 is atom1
            ):
                return b
        return None

@pytest.fixture
def scene():
    return MockScene()

def test_calculate_6ring_rotation_empty(scene):
    """Test with no existing bonds."""
    num_points = 6
    bonds_info = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    atom_items = [MagicMock() for _ in range(6)]
    for a in atom_items:
        a.bonds = []

    rot = scene._calculate_6ring_rotation(num_points, bonds_info, atom_items)
    assert rot == 0

def test_calculate_6ring_rotation_single_edge_single(scene):
    """Test fusing on a single bond (order 1). Should prefer alternating (template double)."""
    num_points = 6
    # Template: 0:D, 1:S, 2:D, 3:S, 4:D, 5:S
    bonds_info = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    atom_items = [MagicMock() for _ in range(6)]
    for a in atom_items:
        a.bonds = []

    # Fuse on edge k=0 (between atoms 0 and 1) with existing single bond
    bond = MagicMock()
    bond.order = 1
    bond.atom1, bond.atom2 = atom_items[0], atom_items[1]
    atom_items[0].bonds.append(bond)
    atom_items[1].bonds.append(bond)

    rot = scene._calculate_6ring_rotation(num_points, bonds_info, atom_items)
    # Case B logic: 1-edge fuse.
    # exist_order=1. template_ord at (k_fuse+rot)%6.
    # Score 100 if (exist_order=1 and template_ord=2).
    # rot=0 -> template_ord = bonds_info[0].order = 2. Score 100.
    # rot=1 -> template_ord = bonds_info[1].order = 1. Score 0.
    assert rot % 2 == 0  # Should be 0, 2, or 4

def test_calculate_6ring_rotation_single_edge_double(scene):
    """Test fusing on a double bond. Should prefer 1 or 2 (alternating or matching)."""
    num_points = 6
    bonds_info = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    atom_items = [MagicMock() for _ in range(6)]
    for a in atom_items:
        a.bonds = []

    # Existing double bond on edge k=0
    bond = MagicMock()
    bond.order = 2
    bond.atom1, bond.atom2 = atom_items[0], atom_items[1]
    atom_items[0].bonds.append(bond)
    atom_items[1].bonds.append(bond)

    rot = scene._calculate_6ring_rotation(num_points, bonds_info, atom_items)
    # Case B logic: if exist_order == 2:
    # 100 if template_ord == 1 or 2
    # 50 if adj template orders match alternating preference
    # If rot=0: template_ord=2. score 100. adj are 1 and 1. (50+50 bonus for adjacent single bonds if exist_order=2)
    # Total rot=0: 100 + 50 + 50 = 200.
    assert rot % 2 == 0

def test_calculate_6ring_rotation_multi_edge_fused(scene):
    """Test fusing on two adjacent edges (naphthalene-like)."""
    num_points = 6
    # Template: 0:D, 1:S, 2:D, 3:S, 4:D, 5:S
    bonds_info = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    atom_items = [MagicMock() for _ in range(6)]
    for a in atom_items:
        a.bonds = []

    # Edge k=0 (atoms 0,1) and k=1 (atoms 1,2)
    b1 = MagicMock()
    b1.order = 2
    b1.atom1, b1.atom2 = atom_items[0], atom_items[1]
    atom_items[0].bonds.append(b1)
    atom_items[1].bonds.append(b1)

    b2 = MagicMock()
    b2.order = 1
    b2.atom1, b2.atom2 = atom_items[1], atom_items[2]
    atom_items[1].bonds.append(b2)
    atom_items[2].bonds.append(b2)

    rot = scene._calculate_6ring_rotation(num_points, bonds_info, atom_items)

    # Case A logic: >= 2 edges.
    # Base scoring: 1000 per matched double bond, 100 bonus for matches, -50 for mismatch.
    # rot=0:
    # k=0 (exist 2): template_ord = bonds_info[0]=2. match_double_count=1, match_bonus=100.
    # k=1 (exist 1): template_ord = bonds_info[1]=1. match_bonus=100.
    # Total = 1000 + 200 = 1200.

    # rot=1:
    # k=0 (exist 2): template_ord = bonds_info[1]=1. mismatch_penalty=-50.
    # k=1 (exist 1): template_ord = bonds_info[2]=2. mismatch_penalty=-50.

    # Connection safety factor (5000) also applies.
    assert rot == 0

def test_calculate_6ring_rotation_connection_safety(scene):
    """Verify that 'safe connection' scoring prioritizes rotations where template single bonds connect to fusion points."""
    num_points = 6
    # Template: 0:D, 1:S, 2:D, 3:S, 4:D, 5:S
    bonds_info = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    atom_items = [MagicMock() for _ in range(6)]
    for a in atom_items:
        a.bonds = []

    # Multi-fuse on 0-1 and 2-3 (non-adjacent)
    # k=0 (0-1) and k=2 (2-3)
    b1 = MagicMock(); b1.order = 1; b1.atom1 = atom_items[0]; b1.atom2 = atom_items[1]
    atom_items[0].bonds.append(b1); atom_items[1].bonds.append(b1)

    b2 = MagicMock(); b2.order = 1; b2.atom1 = atom_items[2]; b2.atom2 = atom_items[3]
    atom_items[2].bonds.append(b2); atom_items[3].bonds.append(b2)

    rot = scene._calculate_6ring_rotation(num_points, bonds_info, atom_items)
    # Safe connections (adj template edges are single bonds) have massive 5000 pts bonus.
    # In benzene (2-1-2-1-2-1), double bonds (0, 2, 4) have single-bond neighbors.
    # Single bonds (1, 3, 5) have double-bond neighbors.
    # Thus, the score heavily favors aligning fusion edges k=0,2 with template double bonds at 0,2,4.
    assert rot % 2 == 0
