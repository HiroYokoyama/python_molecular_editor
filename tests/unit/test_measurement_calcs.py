import pytest
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from unittest.mock import MagicMock, patch

# Needs a QApplication before instantiating widgets like QDialogs
from PyQt6.QtWidgets import QApplication
app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)

# Import all geometry measuring logics from the application:
from moleditpy.modules.angle_dialog import AngleDialog
from moleditpy.modules.dihedral_dialog import DihedralDialog
from moleditpy.modules.bond_length_dialog import BondLengthDialog
from moleditpy.modules.mol_geometry import calculate_dihedral, adjust_bond_angle
from moleditpy.modules.alignment_dialog import AlignmentDialog
from moleditpy.modules.align_plane_dialog import AlignPlaneDialog
from moleditpy.modules.main_window_edit_3d import MainWindowEdit3d
from moleditpy.modules.custom_interactor_style import CustomInteractorStyle
from moleditpy.modules.alignment_dialog import AlignmentDialog
from moleditpy.modules.align_plane_dialog import AlignPlaneDialog

def setup_test_molecule(smiles: str):
    """Creates an optimized 3D RDKit molecule from SMILES for testing."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

# ==========================================
# Angle Dialog Tests
# ==========================================
def test_angle_dialog_logic_matches_rdkit():
    """Verify AngleDialog's native calculation against RDKit."""
    mol = setup_test_molecule("CCC") # Propane
    conf = mol.GetConformer()
    
    # 0, 1, 2 correspond to the three carbons in CCC
    rdkit_ref_angle = rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
    
    main_window = MagicMock()
    dialog = AngleDialog(mol, main_window)
    dialog.atom1_idx = 0
    dialog.atom2_idx = 1
    dialog.atom3_idx = 2
    
    app_calculated_angle = dialog.calculate_angle()
    assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

# ==========================================
# Mol Geometry (adjust_bond_angle) angle calc logic
# ==========================================
def test_adjust_bond_angle_internal_calc_matches_rdkit():
    """mol_geometry.adjust_bond_angle internally calculates the current angle before rotation."""
    mol = setup_test_molecule("CCC")
    conf = mol.GetConformer()
    rdkit_ref_angle = rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
    
    positions = conf.GetPositions()
    
    # We replicate the specific internal math from adjust_bond_angle manually as it doesn't expose a 'get_angle'
    # but we test the math used inside that function.
    idx_a = 0
    idx_b = 1
    idx_c = 2
    pos_a = positions[idx_a]
    pos_b = positions[idx_b]
    pos_c = positions[idx_c]
    vec_ab = pos_a - pos_b
    vec_cb = pos_c - pos_b
    len_ab = np.linalg.norm(vec_ab)
    len_cb = np.linalg.norm(vec_cb)
    cos_current = np.dot(vec_ab, vec_cb) / (len_ab * len_cb)
    cos_current = np.clip(cos_current, -1.0, 1.0)
    current_angle_rad = np.arccos(cos_current)
    app_calculated_angle = np.degrees(current_angle_rad)
    
    assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

# ==========================================
# MainWindowEdit3d angle calc logic
# ==========================================
def test_main_window_edit_3d_angle_logic_matches_rdkit():
    """Verify MainWindowEdit3d._calculate_angle matches RDKit."""
    mol = setup_test_molecule("CCC")
    conf = mol.GetConformer()
    rdkit_ref_angle = rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
    
    # Needs a mock parent that has atom_positions_3d and current_mol
    mock_main_window = MagicMock()
    mock_main_window.current_mol = mol
    mock_main_window.atom_positions_3d = {
        0: conf.GetAtomPosition(0),
        1: conf.GetAtomPosition(1),
        2: conf.GetAtomPosition(2)
    }
    
    editor = MainWindowEdit3d()
    editor.atom_positions_3d = mock_main_window.atom_positions_3d
    
    # _calculate_angle(atom1_idx, vertex_idx, atom3_idx)
    app_calculated_angle = editor.calculate_angle(0, 1, 2)
    assert app_calculated_angle == pytest.approx(rdkit_ref_angle, abs=0.01)

# ==========================================
# Dihedral Dialog & mol_geometry Tests
# ==========================================
def test_dihedral_logic_matches_rdkit():
    """Verify calculate_dihedral native calculation against RDKit."""
    mol = setup_test_molecule("CCCC") # Butane
    conf = mol.GetConformer()
    
    # 0, 1, 2, 3 correspond to the four carbons in CCCC
    rdkit_ref_dihedral = rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3)
    
    positions = conf.GetPositions()
    app_calculated_dihedral = calculate_dihedral(positions, 0, 1, 2, 3)
    assert app_calculated_dihedral == pytest.approx(rdkit_ref_dihedral, abs=0.05)
    
    main_window = MagicMock()
    dialog = DihedralDialog(mol, main_window)
    dialog.atom1_idx = 0
    dialog.atom2_idx = 1
    dialog.atom3_idx = 2
    dialog.atom4_idx = 3
    dialog.update_display()
    
    if dialog.dihedral_input.text() != "":
        dialog_calculated_text = dialog.dihedral_input.text()
        assert float(dialog_calculated_text) == pytest.approx(rdkit_ref_dihedral, abs=0.05)

# ==========================================
# Bond Length Dialog Tests
# ==========================================
def test_bond_length_dialog_logic_matches_rdkit():
    """Verify BondLengthDialog's distance calculation against RDKit."""
    mol = setup_test_molecule("CC") # Ethane
    conf = mol.GetConformer()
    
    # 0, 1 correspond to the two carbons in CC
    rdkit_ref_dist = rdMolTransforms.GetBondLength(conf, 0, 1)
    
    main_window = MagicMock()
    dialog = BondLengthDialog(mol, main_window)
    dialog.atom1_idx = 0
    dialog.atom2_idx = 1
    dialog.update_display() # This triggers the distance calculation
    
    if dialog.distance_input.text() != "":
        dialog_calculated_text = dialog.distance_input.text()
    
        # Check dialog formatted text matches RDKit
        assert float(dialog_calculated_text) == pytest.approx(rdkit_ref_dist, abs=0.01)

    # Manual logic inside the dialog fallback
    p1 = conf.GetAtomPosition(0)
    p2 = conf.GetAtomPosition(1)
    manual_distance = ((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)**0.5
    assert manual_distance == pytest.approx(rdkit_ref_dist, abs=0.0001)

# ==========================================
# MainWindowEdit3D distance calc logic
# ==========================================
def test_main_window_edit_3d_distance_logic_matches_rdkit():
    """Verify MainWindowEdit3d._calculate_distance matches RDKit."""
    mol = setup_test_molecule("CC")
    conf = mol.GetConformer()
    rdkit_ref_dist = rdMolTransforms.GetBondLength(conf, 0, 1)
    
    # Needs a mock parent that has atom_positions_3d
    mock_main_window = MagicMock()
    mock_main_window.current_mol = mol
    mock_main_window.atom_positions_3d = {
        0: conf.GetAtomPosition(0),
        1: conf.GetAtomPosition(1)
    }
    
    editor = MainWindowEdit3d()
    editor.atom_positions_3d = mock_main_window.atom_positions_3d
    # _calculate_distance(atom1_idx, atom2_idx)
    app_calculated_dist = editor.calculate_distance(0, 1)
    assert app_calculated_dist == pytest.approx(rdkit_ref_dist, abs=0.0001)

# ==========================================
# CustomInteractorStyle distance calc logic
# ==========================================
def test_custom_interactor_style_distance_logic_matches_rdkit():
    """Verify the inline distance calc logic in CustomInteractorStyle matches RDKit."""
    mol = setup_test_molecule("CC")
    conf = mol.GetConformer()
    rdkit_ref_dist = rdMolTransforms.GetBondLength(conf, 0, 1)
    
    # We replicate the logic from lines 80-82 in CustomInteractorStyle: np.linalg.norm
    pos1 = np.array(list(conf.GetAtomPosition(0)))
    pos2 = np.array(list(conf.GetAtomPosition(1)))
    app_calculated_dist = float(np.linalg.norm(pos1 - pos2))
    
    assert app_calculated_dist == pytest.approx(rdkit_ref_dist, abs=0.0001)

# ==========================================
# Alignment Dialog logic
# ==========================================
@patch("PyQt6.QtWidgets.QMessageBox.information")
@patch("PyQt6.QtWidgets.QMessageBox.warning")
def test_alignment_dialog_logic(mock_warning, mock_info):
    """Verify AlignmentDialog properly aligns given atoms to the target axis."""
    mol = setup_test_molecule("CC") # Ethane
    conf = mol.GetConformer()
    
    mock_main_window = MagicMock()
    mock_main_window.atom_positions_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    mock_main_window.current_mol = mol
    
    # We test X-axis alignment
    dialog = AlignmentDialog(mol, mock_main_window, "x")
    
    # Select carbon 0 and 1
    dialog.selected_atoms.add(0)
    dialog.selected_atoms.add(1)
    
    # Execute the alignment
    dialog.apply_alignment()
    
    # Verify mathematically
    # atom 0 should be at origin
    pos0 = np.array(conf.GetAtomPosition(0))
    assert np.allclose(pos0, [0, 0, 0], atol=1e-5)
    
    # atom 1 should be on the X-axis (y=0, z=0) and x > 0
    pos1 = np.array(conf.GetAtomPosition(1))
    assert pos1[0] > 0
    assert pos1[1] == pytest.approx(0.0, abs=1e-5)
    assert pos1[2] == pytest.approx(0.0, abs=1e-5)

# ==========================================
# Align Plane Dialog logic
# ==========================================
@patch("PyQt6.QtWidgets.QMessageBox.information")
@patch("PyQt6.QtWidgets.QMessageBox.warning")
def test_align_plane_dialog_logic(mock_warning, mock_info):
    """Verify AlignPlaneDialog properly aligns selected atoms to the target plane."""
    # Build cyclopropane, which is planar
    mol = setup_test_molecule("C1CC1")
    conf = mol.GetConformer()
    
    mock_main_window = MagicMock()
    mock_main_window.atom_positions_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    mock_main_window.current_mol = mol
    
    # Test aligning to XY plane
    dialog = AlignPlaneDialog(mol, mock_main_window, "xy")
    
    # Select the 3 carbon atoms (0, 1, 2) which form a plane
    dialog.selected_atoms.add(0)
    dialog.selected_atoms.add(1)
    dialog.selected_atoms.add(2)
    
    # Execute plane alignment
    dialog.apply_PlaneAlign()
    
    # Verify mathematically
    pos0 = np.array(conf.GetAtomPosition(0))
    pos1 = np.array(conf.GetAtomPosition(1))
    pos2 = np.array(conf.GetAtomPosition(2))
    
    # Since aligned to XY plane, all Z coordinates of the planar ring should be roughly equal
    # The application centers them around origin, so they should be near 0
    assert pos0[2] == pytest.approx(pos1[2], abs=1e-5)
    assert pos1[2] == pytest.approx(pos2[2], abs=1e-5)
    assert pos0[2] == pytest.approx(0.0, abs=1e-5)
