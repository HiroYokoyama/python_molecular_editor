import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from rdkit import Chem
from rdkit.Chem import AllChem

from moleditpy.ui.view_3d_logic import View3DManager

def _make_view3d(mock_host):
    """Create a View3DManager instance with the given mock host."""
    view3d = View3DManager(mock_host)
    # Ensure standard defaults
    view3d._drawing_3d = False
    view3d.current_3d_style = "ball_and_stick"
    return view3d

@pytest.fixture
def mock_pv():
    with patch('moleditpy.ui.view_3d_logic.pv') as mock:
        # Mock PolyData and its glyph/tube methods
        mock_poly = MagicMock()
        mock.PolyData.return_value = mock_poly
        mock_poly.glyph.return_value = MagicMock()
        mock_poly.tube.return_value = MagicMock()

        # Mock Spline
        mock_spline = MagicMock()
        mock.Spline.return_value = mock_spline
        mock_spline.tube.return_value = MagicMock()

        # Mock Light
        mock.Light.return_value = MagicMock()

        yield mock

def test_add_3d_atom_glyphs_styles(app, mock_parser_host, mock_pv):
    """Verify radii and resolutions for different styles in _add_3d_atom_glyphs."""
    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    sym = ["C"]
    col = np.array([[0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    # Check Ball and Stick
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "ball_and_stick", True, mesh_props)
    assert mock_pv.PolyData.call_count >= 1
    # Check that radii were set (VDW_RADII for C is 0.4 * 1.0 = 0.4)
    mock_pv.PolyData.call_args[0][0]
    # Wait, the radii are stored in the mesh, not passed as first arg to PolyData
    # Actually, lines 330-332:
    # self.glyph_source = pv.PolyData(self.atom_positions_3d)
    # self.glyph_source["colors"] = col
    # self.glyph_source["radii"] = rad

    # Let's verify the rad array calculation directly or by checking what's assigned to the mock
    mock_poly = mock_pv.PolyData.return_value
    assert "radii" in mock_poly.__setitem__.call_args_list[1][0]
    rad_array = mock_poly.__setitem__.call_args_list[1][0][1]
    assert np.isclose(rad_array[0], 0.51) # Default VDW for C (1.7 * 0.3)

    # Check CPK
    mock_pv.PolyData.reset_mock()
    mock_parser_host.init_manager.settings["cpk_atom_scale"] = 1.0
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "cpk", True, mesh_props)
    rad_array = mock_pv.PolyData.return_value.__setitem__.call_args_list[1][0][1]
    # pt.GetRvdw(6) for C is ~1.7
    assert rad_array[0] > 1.5

    # Check Stick
    mock_pv.PolyData.reset_mock()
    mock_parser_host.init_manager.settings["stick_bond_radius"] = 0.15
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "stick", True, mesh_props)
    rad_array = mock_pv.PolyData.return_value.__setitem__.call_args_list[1][0][1]
    assert np.isclose(rad_array[0], 0.15)

def test_add_3d_atom_glyphs_stick_split(app, mock_parser_host, mock_pv):
    """Verify that terminal multiple bonds lead to atom splitting in stick mode."""
    view3d = _make_view3d(mock_parser_host)

    # Ethylene (C=C) - two terminal carbons with a double bond
    mol = Chem.MolFromSmiles("C=C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    sym = ["C", "C"]
    col = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    # Mock positions for better control
    view3d.atom_positions_3d = np.array([[0.0, 0.0, 0.0], [1.33, 0.0, 0.0]])

    mock_pv.PolyData.reset_mock()
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "stick", True, mesh_props)

    # In stick mode, a double bond terminal atom should be split into 2
    # So for 2 atoms with a double bond, they are both terminal, so we expect 2*2 = 4 spheres
    # Wait, the code says:
    # if atom.GetDegree() == 1: # Terminal atom
    # for C=C, both carbons have degree 1 (only 1 bond each).
    # So both should be split.
    # Total positions should be 4.

    # Check PolyData calls
    # First call is to create glyph_source = pv.PolyData(self.atom_positions_3d) at line 330
    # Second call is at line 472: glyph_source = pv.PolyData(np.array(new_positions))
    assert mock_pv.PolyData.call_count >= 2
    new_positions = mock_pv.PolyData.call_args_list[-1][0][0]
    assert len(new_positions) == 4

def test_add_3d_bond_cylinders_basic(app, mock_parser_host, mock_pv):
    """Verify single, double, and triple bond generation."""
    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("C#CC=CC")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    ["C"] * 5
    col = np.array([[0.5, 0.5, 0.5]] * 5)
    mesh_props = {"smooth_shading": True}

    # Mock plotter.add_mesh
    mock_parser_host.view_3d_manager.plotter.add_mesh = MagicMock()

    view3d._add_3d_bond_cylinders(mol, conf, col, "ball_and_stick", mesh_props)

    # C#CC=CC has 4 bonds: 1 Triple, 1 Single, 1 Double, 1 Single
    # Triple -> 3 segments
    # Single -> 1 segment (Uniform B&S)
    # Double -> 2 segments
    # Single -> 1 segment
    # Total segments = 3 + 1 + 2 + 1 = 7 segments

    # Check that segments were added to the lists and then PolyData created
    # Only one PolyData call expected since we didn't call _add_3d_atom_glyphs
    assert mock_pv.PolyData.call_count >= 1
    points = mock_pv.PolyData.call_args_list[-1][0][0]
    # Each segment has 2 points
    assert len(points) == 14

def test_add_3d_bond_cylinders_styles(app, mock_parser_host, mock_pv):
    """Check style-dependent factors (radius/offset factors) in bond drawing."""
    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("C=C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    col = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    # 1. Ball and Stick
    mock_parser_host.init_manager.settings.update({
        "ball_stick_bond_radius": 0.1,
        "ball_stick_double_bond_radius_factor": 0.8,
        "ball_stick_double_bond_offset_factor": 2.0
    })
    view3d._add_3d_bond_cylinders(mol, conf, col, "ball_and_stick", mesh_props)
    # Double bond radius should be 0.1 * 0.8 = 0.08
    # radii are set on PolyData.point_data["radii"]
    mock_pd = mock_pv.PolyData.return_value
    radii = mock_pd.point_data.__setitem__.call_args_list[0][0][1]
    assert np.allclose(radii, 0.08)

    # 2. Stick
    mock_pv.PolyData.reset_mock()
    mock_pd.point_data.reset_mock()
    mock_parser_host.init_manager.settings.update({
        "stick_bond_radius": 0.15,
        "stick_double_bond_radius_factor": 0.6,
        "stick_double_bond_offset_factor": 1.5
    })
    view3d._add_3d_bond_cylinders(mol, conf, col, "stick", mesh_props)
    # Double bond radius should be 0.15 * 0.6 = 0.09
    radii = mock_pd.point_data.__setitem__.call_args_list[0][0][1]
    assert np.allclose(radii, 0.09)

def test_add_3d_bond_cylinders_overrides(app, mock_parser_host, mock_pv):
    """Verify plugin bond color overrides."""
    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("CC")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    col = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]) # White atoms
    mesh_props = {"smooth_shading": True}

    # Set override for bond 0
    view3d._plugin_bond_color_overrides = {0: "#FF0000"} # Red

    view3d._add_3d_bond_cylinders(mol, conf, col, "stick", mesh_props)

    # colors are set on PolyData.cell_data["colors"]
    mock_pd = mock_pv.PolyData.return_value
    colors = mock_pd.cell_data.__setitem__.call_args_list[0][0][1]
    # Red is [255, 0, 0]
    assert np.array_equal(colors[0], [255, 0, 0])

def test_add_3d_aromatic_rings(app, mock_parser_host, mock_pv):
    """Verify aromatic torus generation."""
    view3d = _make_view3d(mock_parser_host)

    # Benzene
    mol = Chem.MolFromSmiles("c1ccccc1")
    AllChem.EmbedMolecule(mol)
    view3d.atom_positions_3d = np.array([list(mol.GetConformer().GetAtomPosition(i)) for i in range(6)])

    mock_parser_host.init_manager.settings["display_aromatic_circles_3d"] = True
    mock_parser_host.view_3d_manager.plotter.add_mesh = MagicMock()

    mesh_props = {"smooth_shading": True}
    view3d._add_3d_aromatic_rings(mol, "ball_and_stick", mesh_props)

    # Should call pv.Spline and then tube()
    assert mock_pv.Spline.call_count == 1
    # Check that add_mesh was called for the torus
    assert mock_parser_host.view_3d_manager.plotter.add_mesh.call_count == 1

def test_calculate_double_bond_offset(app, mock_parser_host):
    """Verify neighbor-based plane calculation for double bond offset."""
    view3d = _make_view3d(mock_parser_host)

    # Ethene
    mol = Chem.MolFromSmiles("C=C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    bond = mol.GetBondWithIdx(0)

    # Set custom positions to define a plane
    # C0 at (0,0,0), C1 at (1,0,0)
    # H on C0 at (0,1,0), H on C1 at (1,1,0)
    # This should put the double bond in XY plane, offset should be along Z or Y?

    # Creating a molecule with neighbors to test neighbor logic
    mol = Chem.MolFromSmiles("C=C(C)C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    bond = mol.GetBondWithIdx(0) # The C=C bond

    offset = view3d._calculate_double_bond_offset(mol, bond, conf)
    assert len(offset) == 3
    assert np.isclose(np.linalg.norm(offset), 1.0)

def test_show_ez_labels_3d(app, mock_parser_host):
    """Verify EZ label detection and discrepancy marking."""
    view3d = _make_view3d(mock_parser_host)

    # trans-2-butene (E)
    mol = Chem.MolFromSmiles("C/C=C/C")
    AllChem.EmbedMolecule(mol)

    mock_parser_host.view_3d_manager.plotter.add_point_labels = MagicMock()

    view3d.show_ez_labels_3d(mol)

    # Check if add_point_labels was called with "E"
    args, kwargs = mock_parser_host.view_3d_manager.plotter.add_point_labels.call_args
    assert "E" in args[1]

    # Test discrepancy
    # Mock 2D metadata to say it's Z (3)
    mock_parser_host.state_manager.data.bonds = { (1, 2): {"stereo": 3} }
    # We need to make sure the atoms have _original_atom_id
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i)

    mock_parser_host.view_3d_manager.plotter.add_point_labels.reset_mock()
    view3d.show_ez_labels_3d(mol)
    args, kwargs = mock_parser_host.view_3d_manager.plotter.add_point_labels.call_args
    assert "?" in args[1]

def test_chiral_labels_logic(app, mock_parser_host, mock_pv):
    """Verify chiral label toggling and update."""
    view3d = _make_view3d(mock_parser_host)
    # Ensure drawing logic is safe
    mock_parser_host.statusBar.return_value = MagicMock()
    mock_parser_host.view_3d_manager.plotter = MagicMock()

    # Chiral carbon
    mol = Chem.MolFromSmiles("C[C@H](F)Cl")
    AllChem.EmbedMolecule(mol)
    # Set original atom IDs for mapping
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i + 1)

    # Ensure view3d has current_mol and host has it too
    # Mocking the plotter on the manager itself because it acts as host.view_3d_manager
    view3d.plotter = MagicMock()
    view3d.host.view_3d_manager = view3d
    view3d.current_mol = mol

    view3d.plotter.add_point_labels = MagicMock()

    view3d.toggle_chiral_labels_display(True)
    assert view3d.show_chiral_labels is True

    # Directly test update_chiral_labels
    atom_item = MagicMock()
    # RDKit idx 1 is the chiral center in C[C@H](F)Cl (atoms: C0, C1, F2, Cl3)
    # Its _original_atom_id is 1+1 = 2
    mock_parser_host.state_manager.data.atoms = { 2: {"item": atom_item, "symbol": "C"} }

    with patch('moleditpy.ui.view_3d_logic.Chem.FindMolChiralCenters', return_value=[(1, 'S')]):
        view3d.update_chiral_labels()
        assert atom_item.chiral_label == 'S'

def test_color_overrides(app, mock_parser_host, mock_pv):
    """Verify color override API functions."""
    view3d = _make_view3d(mock_parser_host)
    view3d.draw_molecule_3d = MagicMock()
    mock_parser_host.view_3d_manager.current_mol = Chem.MolFromSmiles("C")

    # Bond override
    view3d.update_bond_color_override(0, "#FF0000")
    assert view3d._plugin_bond_color_overrides[0] == "#FF0000"
    view3d.draw_molecule_3d.assert_called()

    # Atom override
    view3d.update_atom_color_override(0, "#00FF00")
    assert view3d._plugin_color_overrides[0] == "#00FF00"
    assert view3d.draw_molecule_3d.call_count == 2
