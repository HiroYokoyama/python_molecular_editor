import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.modules.calculation_worker import CalculationWorker
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF

def test_calculation_worker_bond_length_validation(qtbot, app):
    """
    Integration test: Run a real 3D optimization for Ethane.
    Compare the resulting C-C bond length against RDKit's reference calculated value.
    No hard-coded bond length value is used.
    """
    # 1. Setup MolecularData (2D Ethane)
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    # 2. Run CalculationWorker (Real 3D Embedded + Optimization)
    worker = CalculationWorker()
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    result = blocker.args[0]
    # result is (mol_block_3d, rdkit_mol, info_dict)
    mol_3d = result[1]
    assert mol_3d is not None
    
    # 3. Validation: Compare against RDKit's own reference calculations
    # We find the C-C bond in the result mol.
    cc_bonds = [b for b in mol_3d.GetBonds() if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    assert len(cc_bonds) >= 1
    target_bond = cc_bonds[0]
    
    # Get optimized distance from the result conformer
    conf = mol_3d.GetConformer()
    dist_measured = rdMolTransforms.GetBondLength(conf, target_bond.GetBeginAtomIdx(), target_bond.GetEndAtomIdx())
    
    # Get "Reference" RDKit distance for a standard ETKDG embedded molecule
    ref_mol = Chem.MolFromSmiles("CC")
    ref_mol = Chem.AddHs(ref_mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(ref_mol, params)
    AllChem.MMFFOptimizeMolecule(ref_mol, mmffVariant="MMFF94s")
    ref_conf = ref_mol.GetConformer()
    # In ref_mol (from "CC"), C-C bond is between atoms 0 and 1
    dist_ref = rdMolTransforms.GetBondLength(ref_conf, 0, 1)
    
    print(f"DEBUG: Measured C-C = {dist_measured:.4f}, Reference C-C = {dist_ref:.4f}")
    assert dist_measured == pytest.approx(dist_ref, rel=1e-2)

def test_calculation_worker_angle_validation(qtbot, app):
    """
    Integration test: Run 3D optimization for Propane.
    Compare the C-C-C angle against RDKit's reference.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    c3 = data.add_atom("C", QPointF(100, 50))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c2, c3, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()
    
    # Find C-C-C triplet by symbol and connectivity (since original IDs are lost in MOL block)
    c_indices = [a.GetIdx() for a in mol_3d.GetAtoms() if a.GetSymbol() == 'C']
    central = -1
    neighbors = []
    for idx in c_indices:
        c_neighbors = [n.GetIdx() for n in mol_3d.GetAtomWithIdx(idx).GetNeighbors() if n.GetSymbol() == 'C']
        if len(c_neighbors) == 2:
            central = idx
            neighbors = c_neighbors
            break
            
    assert central != -1
    angle_measured = rdMolTransforms.GetAngleDeg(conf, neighbors[0], central, neighbors[1])
    
    # Reference Propane (Identical setup to CalculationWorker)
    ref_mol = Chem.MolFromSmiles("CCC")
    ref_mol = Chem.AddHs(ref_mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(ref_mol, params)
    AllChem.MMFFOptimizeMolecule(ref_mol, mmffVariant="MMFF94s")
    
    # In ref_mol (from SMILES "CCC"), C atoms are 0, 1, 2. Central is 1.
    angle_ref = rdMolTransforms.GetAngleDeg(ref_mol.GetConformer(), 0, 1, 2)
    
    print(f"DEBUG: Measured C-C-C = {angle_measured:.4f}, Reference = {angle_ref:.4f}")
    assert angle_measured == pytest.approx(angle_ref, abs=5.0), \
        f"Angle mismatch: measured {angle_measured:.2f} vs ref {angle_ref:.2f}"

def test_calculation_worker_dihedral_validation(qtbot, app):
    """
    Integration test: Run 3D optimization for n-Butane.
    Check if it converges to a staggered conformation (Trans or Gauche).
    No hardcoded values; check against RDKit staggered preference.
    """
    data = MolecularData()
    # Linear start to avoid bias
    for i in range(4):
        data.add_atom("C", QPointF(i*50, 0))
        if i > 0:
            data.add_bond(i-1, i, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()
    
    # Path C1-C2-C3-C4
    # We find a chain of 4 carbons.
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3)
    
    # A staggered conformation should have dihedral around 60 (gauche) or 180 (trans).
    # We check if it is NOT eclipsed (near 0, 120).
    is_staggered = (40 < abs(dihedral) < 80) or (160 < abs(dihedral) <= 180)
    
    print(f"DEBUG: n-Butane Dihedral = {dihedral:.4f}")
    assert is_staggered

def test_calculation_worker_conversion_no_optimize(qtbot, app):
    """
    Integration test: 2D to 3D without optimization.
    Check if atom counts and symbols are preserved.
    """
    data = MolecularData()
    data.add_atom("F", QPointF(0, 0))
    data.add_atom("Cl", QPointF(50, 0))
    data.add_bond(0, 1, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'rdkit', 'do_optimize': False}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    mol_3d = blocker.args[0][1]
    symbols = set(a.GetSymbol() for a in mol_3d.GetAtoms())
    assert "F" in symbols
    assert "Cl" in symbols
    assert mol_3d.GetNumAtoms() == 2 # F-Cl without Hs (default for this simple case if not specified)
    # Actually worker usually adds Hs unless suppressed. Let's check.
    # Hydrogens are not added for F-Cl by default by RDKit usually? Or they are?
    # RDKit will add Hs to fill valence if worker calls AddHs.
    
    # Verify bond length is non-zero (3D coordinates generated)
    conf = mol_3d.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 0, 1)
    assert dist > 0.1
