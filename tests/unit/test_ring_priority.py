import pytest
from unittest.mock import MagicMock
from moleditpy.core.molecular_data import MolecularData

def test_ring_priority_smaller_wins():
    """Verify that smaller rings are prioritized for double bond shift logic."""
    data = MolecularData()
    
    # Create a 5/6 fused system (Indole-like skeleton)
    # 6-ring: 0-1-2-3-4-5
    # 5-ring: 6-4-3-8-7
    # Shared bond: 3-4
    
    # 6-ring atoms
    data.add_atom("C", (0, 0))   # 0
    data.add_atom("C", (10, 0))  # 1
    data.add_atom("C", (20, 10)) # 2
    data.add_atom("C", (20, 20)) # 3 (Shared)
    data.add_atom("C", (10, 30)) # 4 (Shared)
    data.add_atom("C", (0, 30))  # 5
    
    # 5-ring atoms
    data.add_atom("C", (30, 25)) # 6
    data.add_atom("C", (30, 15)) # 7
    data.add_atom("C", (25, 20)) # 8
    
    # Bonds for 6-ring
    data.add_bond(0, 1)
    data.add_bond(1, 2)
    data.add_bond(2, 3)
    data.add_bond(3, 4, order=2) # Double bond in shared position
    data.add_bond(4, 5)
    data.add_bond(5, 0)
    
    # Bonds for 5-ring
    data.add_bond(3, 8)
    data.add_bond(8, 7)
    data.add_bond(7, 6)
    data.add_bond(6, 4)
    
    # Mock items for bonds and atoms
    for atom_id in data.atoms:
        item = MagicMock()
        item.atom_id = atom_id # Set atom_id for the mapping logic
        # Set up pos() to return a mock that responds to x() and y() with the actual coordinates
        pos_mock = MagicMock()
        pos_mock.x.return_value = float(data.atoms[atom_id]["pos"][0])
        pos_mock.y.return_value = float(data.atoms[atom_id]["pos"][1])
        item.pos.return_value = pos_mock
        data.atoms[atom_id]["item"] = item
        
    for bond_key in data.bonds:
        bond_item = MagicMock()
        data.bonds[bond_key]["item"] = bond_item
        
    # Run ring info update
    data.update_ring_info_2d()
    
    # Check shared bond (3, 4)
    shared_bond_key = (3, 4)
    if shared_bond_key not in data.bonds:
        shared_bond_key = (4, 3)
    
    shared_bond_item = data.bonds[shared_bond_key]["item"]
    assert shared_bond_item.is_in_ring is True
    
    # Calculate ring centers manually
    # Ring 0 (atoms 0,1,2,3,4,5) - 6-ring
    six_ring_atoms = [0, 1, 2, 3, 4, 5]
    six_center_x = sum(data.atoms[i]["pos"][0] for i in six_ring_atoms) / 6
    six_center_y = sum(data.atoms[i]["pos"][1] for i in six_ring_atoms) / 6
    six_center = (six_center_x, six_center_y)
    
    # Ring 1 (atoms 3,4,6,7,8) - 5-ring
    five_ring_atoms = [3, 4, 6, 7, 8]
    five_center_x = sum(data.atoms[i]["pos"][0] for i in five_ring_atoms) / 5
    five_center_y = sum(data.atoms[i]["pos"][1] for i in five_ring_atoms) / 5
    five_center = (five_center_x, five_center_y)
    
    # We expect the 5-ring center to be used because it's smaller
    assert shared_bond_item.ring_center == five_center
