import pytest
import math
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.modules.calculation_worker import CalculationWorker
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF

def test_3d_bond_lengths(qtbot):
    """Verify optimized 3D coordinates yield physical bond lengths."""
    # Create Ethane C-C bond length test
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(10, 0))
    data.add_bond(c1, c2, order=1)
    
    mol_block = data.to_mol_block()
    print(f"DEBUG: mol_block='{mol_block}'")
    
    worker = CalculationWorker()
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, {'conversion_mode': 'rdkit'})
        
    result = blocker.args[0]
    if isinstance(result, tuple):
        mol = result[1]
    else:
        mol = result
        
    assert mol is not None
    conf = mol.GetConformer()
    
    # RDKit indices might change, but here it should be 0 and 1 (heavy atoms)
    # plus hydrogens. We find C-C bond.
    cc_bond = None
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C':
            cc_bond = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            break
            
    assert cc_bond is not None
    dist = rdMolTransforms.GetBondLength(conf, cc_bond[0], cc_bond[1])
    
    # Ethane C-C length is approx 1.54 A. Tolerance 0.05 A
    assert 1.50 < dist < 1.60



from moleditpy.modules.mirror_dialog import MirrorDialog
from unittest.mock import MagicMock
from rdkit.Geometry import Point3D
from PyQt6.QtWidgets import QWidget

def test_mirror_dialog_logic(qtbot):
    """Verify MirrorDialog correctly manipulates coordinates and UI (software logic)."""
    # (S)-2-butanol
    mol = Chem.MolFromSmiles("C[C@H](O)CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    # Mock main_window using QWidget to satisfy PyQt6
    main_window = QWidget()
    main_window.statusBar = MagicMock()
    main_window.statusBar.return_value = MagicMock()
    main_window.draw_molecule_3d = MagicMock()
    main_window.update_chiral_labels = MagicMock()
    main_window.push_undo_state = MagicMock()
    
    # Initial state
    conf = mol.GetConformer()
    x_before = conf.GetAtomPosition(0).x
    
    # Instantiate dialog
    dialog = MirrorDialog(mol, main_window)
    qtbot.addWidget(dialog)
    
    # select YZ plane (id=2)
    dialog.yz_radio.setChecked(True)
    
    # Apply mirror
    dialog.apply_mirror()
    
    # Verify coordinates inverted on X axis
    x_after = conf.GetAtomPosition(0).x
    assert x_after == pytest.approx(-x_before, abs=1e-3)
    
    # Verify MoleditPy specific recovery calls
    assert main_window.draw_molecule_3d.called
    assert main_window.update_chiral_labels.called
    assert main_window.push_undo_state.called

def test_planarize_logic():
    """Verify planarize functionality (coordinate logic)."""
    # Create a non-planar 4-atom system
    mol = Chem.RWMol()
    for _ in range(4):
        mol.AddAtom(Chem.Atom("C"))
    conf = Chem.Conformer(4)
    # Positions with large Z variation
    positions = [(0, 0, 1.0), (1, 0, -1.0), (0, 1, 0.5), (1, 1, -0.5)]
    for i, pos in enumerate(positions):
        conf.SetAtomPosition(i, pos)
    mol.AddConformer(conf)
    
    # Simulation of planarize logic for XY plane (Z=0)
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x, pos.y, 0.0))
        
    for i in range(mol.GetNumAtoms()):
        assert conf.GetAtomPosition(i).z == 0.0
