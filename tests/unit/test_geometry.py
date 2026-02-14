import pytest
import math
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.modules.calculation_worker import CalculationWorker
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest import mock as _mock

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

from moleditpy.modules.planarize_dialog import PlanarizeDialog
import numpy as np

def test_planarize_logic(qtbot):
    """Verify planarize functionality using the actual PlanarizeDialog logic."""
    # (R)-2-butanol
    mol = Chem.MolFromSmiles("C[C@H](O)CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()
    
    # Mock main_window
    main_window = QWidget()
    main_window.atom_positions_3d = mol.GetConformer().GetPositions()
    main_window.plotter = MagicMock()
    main_window.draw_molecule_3d = MagicMock()
    main_window.update_chiral_labels = MagicMock()
    main_window.push_undo_state = MagicMock()
    
    # Instantiate dialog with all atoms selected
    dialog = PlanarizeDialog(mol, main_window, preselected_atoms=range(mol.GetNumAtoms()))
    qtbot.addWidget(dialog)
    
    # Apply planarize
    # We mock QMessageBox to prevent blocking
    with _mock.patch('PyQt6.QtWidgets.QMessageBox.information'):
        dialog.apply_planarize()
    
    new_positions = conf.GetPositions()
    centroid = np.mean(new_positions, axis=0)
    centered = new_positions - centroid
    
    u, s, vh = np.linalg.svd(centered)
    # The smallest singular value s[-1] should be near 0 if they are planar
    assert s[-1] < 1e-10
    
    # Also verify that main_window methods were called
    assert main_window.draw_molecule_3d.called
    assert main_window.update_chiral_labels.called
    assert main_window.push_undo_state.called
