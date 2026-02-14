import pytest
from rdkit import Chem
from moleditpy.modules.calculation_worker import CalculationWorker
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF


def test_ez_preservation_logic(qtbot):
    """Verify worker preserves explicit E/Z labels even if RDKit might lose them."""
    # Create a MOL block with explicit Z stereo (stereo=3)
    data = MolecularData()
    # Cis-2-butene structure
    data.add_atom("C", QPointF(-10, 10))
    data.add_atom("C", QPointF(0, 0))
    data.add_atom("C", QPointF(10, 0))
    data.add_atom("C", QPointF(20, 10))
    data.add_bond(1, 2, order=2, stereo=3) # Z
    data.add_bond(0, 1, order=1)
    data.add_bond(2, 3, order=1)
    
    # We use our own SERIALIZATION logic to ensure M CFG or V2000 stereo code is present
    mol_block = data.to_mol_block() 
    
    worker = CalculationWorker()
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, {'conversion_mode': 'rdkit'})
        
    res = blocker.args[0]
    mol = res[1] if isinstance(res, tuple) else res
    
    # The double bond should persist as Z
    bond = mol.GetBondBetweenAtoms(1, 2)
    assert bond.GetStereo() == Chem.BondStereo.STEREOZ

def test_molecular_data_fallback_serialization():
    """Verify MolecularData uses manual string construction if RDKit fails (fallback logic)."""
    data = MolecularData()
    # Intentionally broken data that might confuse RDKit but we want to serialize anyway
    data.add_atom("X", QPointF(0, 0)) # Non-standard element
    data.add_atom("C", QPointF(10, 0))
    data.add_bond(0, 1, order=1)
    
    # RDKit might fail to sanitize "X", so to_rdkit_mol returns None
    from unittest.mock import MagicMock
    item1 = MagicMock()
    item1.pos.return_value = QPointF(0,0)
    data.atoms[0]['item'] = item1
    
    item2 = MagicMock()
    item2.pos.return_value = QPointF(10,0)
    data.atoms[1]['item'] = item2

    # to_mol_block should then use the fallback path
    mol_block = data.to_mol_block()
    
    assert mol_block is not None
    assert "MoleditPy" in mol_block
    assert "X  " in mol_block
    assert "C  " in mol_block

def test_coordinate_mapping_primary():
    """Verify primary RDKit conversion logic uses 'pos' attribute correctly."""
    from moleditpy.modules.constants import ANGSTROM_PER_PIXEL
    
    data = MolecularData()
    # RDKit path uses 'pos' passed in add_atom
    data.add_atom("C", QPointF(0, 0))
    data.add_atom("C", QPointF(100, 0)) 
    data.add_bond(0, 1, order=1)
    
    mol_block = data.to_mol_block()
    
    expected_x = 100.0 * ANGSTROM_PER_PIXEL
    
    lines = mol_block.split('\n')
    atom_lines = [l for l in lines if " C " in l]
    assert len(atom_lines) == 2
    
    parts = atom_lines[1].split()
    x_coord = float(parts[0])
    
    print(f"DEBUG: Primary x_coord={x_coord}, expected={expected_x}")
    assert abs(x_coord - expected_x) < 0.0001

def test_coordinate_mapping_fallback():
    """Verify fallback serialization (reverted logic) uses atom['item'].pos()."""
    from moleditpy.modules.constants import ANGSTROM_PER_PIXEL
    from unittest.mock import MagicMock
    
    data = MolecularData()
    # Use invalid atom symbol to force RDKit serialization to fail/return None
    # to_rdkit_mol sanitization will likely fail or return None for "X" depending on config
    # We'll use "X" and ensure to_rdkit_mol returns None or raises to trigger fallback
    data.add_atom("X", QPointF(0, 0))
    data.add_atom("X", QPointF(0, 0)) 
    data.add_bond(0, 1, order=1)
    
    # Mock items for fallback logic
    item1 = MagicMock()
    item1.pos.return_value = QPointF(0.0, 0.0)
    data.atoms[0]['item'] = item1
    
    item2 = MagicMock()
    item2.pos.return_value = QPointF(100.0, 0.0)
    data.atoms[1]['item'] = item2
    
    mol_block = data.to_mol_block()
    assert mol_block is not None
    assert "MoleditPy" in mol_block # Check header of fallback block
    
    expected_x = 100.0 * ANGSTROM_PER_PIXEL
    
    lines = mol_block.split('\n')
    atom_lines = [l for l in lines if " X " in l]
    assert len(atom_lines) == 2
    
    parts = atom_lines[1].split()
    x_coord = float(parts[0])
    
    print(f"DEBUG: Fallback x_coord={x_coord}, expected={expected_x}")
    assert abs(x_coord - expected_x) < 0.0001
    
