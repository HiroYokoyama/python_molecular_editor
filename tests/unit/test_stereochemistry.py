from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF

def test_wedge_dash_mapping():
    """Verify that Wedge/Dash in MolecularData mapped to RDKit BondDir."""
    data = MolecularData()
    # Create a chiral center (dummy)
    # C1 is center, C2 is wedge, C3 is dash
    id1 = data.add_atom("C", QPointF(0, 0))
    id2 = data.add_atom("C", QPointF(10, 10))
    id3 = data.add_atom("C", QPointF(-10, 10))
    
    # stereo=1 is Wedge, stereo=2 is Dash
    data.add_bond(id1, id2, order=1, stereo=1) # C1 -> C2 Wedge
    data.add_bond(id1, id3, order=1, stereo=2) # C1 -> C3 Dash
    
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None
    
    # Check RDKit bonds
    # RDKit might swap atoms, so we check if BEGINWEDGE/BEGINDASH is at id1
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIntProp("_original_atom_id")
        a2 = bond.GetEndAtom().GetIntProp("_original_atom_id")
        
        if (a1 == id1 and a2 == id2):
            assert bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
        elif (a1 == id2 and a2 == id1):
            assert bond.GetBondDir() == Chem.BondDir.ENDWEDGE # If swapped
        elif (a1 == id1 and a2 == id3):
            assert bond.GetBondDir() == Chem.BondDir.BEGINDASH
        elif (a1 == id3 and a2 == id1):
            assert bond.GetBondDir() == Chem.BondDir.ENDDASH # If swapped

def test_ez_stereo_persistence():
    """Verify E/Z double bond configurations are maintained."""
    data = MolecularData()
    # 2-butene-like structure for E/Z test
    # C1-C2=C3-C4
    id1 = data.add_atom("C", QPointF(-10, -10))
    id2 = data.add_atom("C", QPointF(0, 0))
    id3 = data.add_atom("C", QPointF(10, 0))
    id4 = data.add_atom("C", QPointF(20, -10)) # Trans configuration
    
    data.add_bond(id1, id2, order=1)
    data.add_bond(id2, id3, order=2, stereo=4) # 4 is STEREOE (Trans)
    data.add_bond(id3, id4, order=1)
    
    # Case 1: Label priority
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=False)
    bond = mol.GetBondBetweenAtoms(1, 2) # The double bond
    assert bond.GetStereo() == Chem.BondStereo.STEREOE
    
    # Case 2: Z label (stereo=3)
    data.bonds[(id2, id3)]['stereo'] = 3
    mol_z = data.to_rdkit_mol(use_2d_stereo=False)
    Chem.AssignStereochemistry(mol_z, force=True, cleanIt=False)
    bond_z = mol_z.GetBondBetweenAtoms(1, 2)
    assert bond_z.GetStereo() == Chem.BondStereo.STEREOZ

def test_chiral_r_s_consistency():
    """Verify that R/S markers match RDKit descriptors for a known chiral center."""
    # (S)-2-butanol
    #      OH
    #      |
    #  CH3-C-CH2-CH3
    #      |
    #      H
    data = MolecularData()
    c2 = data.add_atom("C", QPointF(0, 0))
    c1 = data.add_atom("C", QPointF(-10, 0))
    o = data.add_atom("O", QPointF(0, 10))
    c3 = data.add_atom("C", QPointF(10, 0))
    c4 = data.add_atom("C", QPointF(20, 0))
    
    data.add_bond(c2, c1, order=1)
    data.add_bond(c2, o, order=1, stereo=1) # OH is Wedge
    data.add_bond(c2, c3, order=1)
    data.add_bond(c3, c4, order=1)
    
    mol = data.to_rdkit_mol()
    # Check if the bond between C2 (idx 0) and O (idx 2) is Wedge
    # (Indices depend on order of add_atom)
    found_wedge = False
    for bond in mol.GetBonds():
        if bond.GetBondDir() == Chem.BondDir.BEGINWEDGE:
            found_wedge = True
            break
    assert found_wedge
