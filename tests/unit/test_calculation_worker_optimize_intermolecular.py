import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from moleditpy.ui.calculation_worker import _iterative_optimize

def test_intermolecular_interaction_toggle():
    """Test that ignoreInterfragInteractions is correctly toggled via options."""
    # Create two methane molecules
    mol1 = Chem.MolFromSmiles("C")
    mol1 = Chem.AddHs(mol1)
    AllChem.EmbedMolecule(mol1, randomSeed=42)
    
    mol2 = Chem.MolFromSmiles("C")
    mol2 = Chem.AddHs(mol2)
    AllChem.EmbedMolecule(mol2, randomSeed=43)
    
    # Combine them into one molecule with two fragments and SANITIZE
    combined = Chem.CombineMols(mol1, mol2)
    Chem.SanitizeMol(combined)
    
    # Move the second fragment 6.0A away along X axis
    conf = combined.GetConformer()
    for i in range(mol1.GetNumAtoms(), combined.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 6.0, pos.y, pos.z))
    
    def check_halted(): return False
    def safe_status(msg): pass

    # --- Case 1: Intermolecular Interaction ON (Default) ---
    mol_on = Chem.Mol(combined)
    options_on = {"optimize_intermolecular_interaction_rdkit": True}
    _iterative_optimize(mol_on, "MMFF94s", check_halted, safe_status, options=options_on, max_iters=500)
    
    dist_on = np.linalg.norm(
        np.array(mol_on.GetConformer().GetAtomPosition(0)) - 
        np.array(mol_on.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )
    
    # --- Case 2: Intermolecular Interaction OFF ---
    mol_off = Chem.Mol(combined)
    options_off = {"optimize_intermolecular_interaction_rdkit": False}
    _iterative_optimize(mol_off, "MMFF94s", check_halted, safe_status, options=options_off, max_iters=500)
    
    dist_off = np.linalg.norm(
        np.array(mol_off.GetConformer().GetAtomPosition(0)) - 
        np.array(mol_off.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )
    
    print(f"Dist ON: {dist_on:.6f}")
    print(f"Dist OFF: {dist_off:.6f}")
    
    # With interaction OFF, the distance should be exactly the initial 6.0
    assert abs(dist_off - 6.0) < 1e-5
    
    # With interaction ON, it should have changed (pulled together)
    assert dist_on < 5.0

def test_intermolecular_interaction_uff():
    """Test UFF path as well."""
    mol1 = Chem.MolFromSmiles("C")
    mol1 = Chem.AddHs(mol1)
    AllChem.EmbedMolecule(mol1, randomSeed=42)
    
    mol2 = Chem.MolFromSmiles("C")
    mol2 = Chem.AddHs(mol2)
    AllChem.EmbedMolecule(mol2, randomSeed=43)
    
    combined = Chem.CombineMols(mol1, mol2)
    Chem.SanitizeMol(combined)
    
    conf = combined.GetConformer()
    for i in range(mol1.GetNumAtoms(), combined.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 6.0, pos.y, pos.z))
    
    def check_halted(): return False
    def safe_status(msg): pass

    # OFF
    mol_off = Chem.Mol(combined)
    options_off = {"optimize_intermolecular_interaction_rdkit": False}
    _iterative_optimize(mol_off, "UFF", check_halted, safe_status, options=options_off, max_iters=500)
    
    dist_off = np.linalg.norm(
        np.array(mol_off.GetConformer().GetAtomPosition(0)) - 
        np.array(mol_off.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )
    
    assert abs(dist_off - 6.0) < 1e-5
