"""Unit tests for View3DManager rendering and camera operations."""

import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.view_3d_logic import View3DManager


def _make_view3d(mock_host):
    """Create a View3DManager instance with the given mock host."""
    from moleditpy.ui.view_3d_logic import View3DManager

    view3d = View3DManager.__new__(View3DManager)
    view3d.host = mock_host
    view3d._drawing_3d = False
    view3d.current_3d_style = "Ball and Stick"
    view3d.atom_info_display_mode = None
    view3d._3d_color_map = {}
    view3d._plugin_color_overrides = {}
    view3d._plugin_bond_color_overrides = {}
    view3d.atom_positions_3d = None
    view3d.atom_label_legend_names = []
    view3d.current_atom_info_labels = None
    view3d.axes_widget = None
    view3d.current_mol = None
    view3d.plotter = MagicMock()
    return view3d


def test_view_3d_draw_standard_3d_style(app, mock_parser_host):
    """Verify that draw_standard_3d_style clears the plotter and constructs the correct VTK meshes."""
    mock_parser_host.init_manager.settings.update(
        {
            "projection_mode": "Perspective",
            "background_color": "#ffffff",
            "display_kekule_3d": False,
        }
    )
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)
    view3d.plotter = MagicMock()

    mol = Chem.MolFromSmiles("CCO")
    AllChem.EmbedMolecule(mol)

    with patch("moleditpy.ui.view_3d_logic.pv") as mock_pv:
        mock_pv.PolyData.return_value.glyph.return_value = MagicMock()
        mock_pv.PolyData.return_value.tube.return_value = MagicMock()
        mock_pv.Light.return_value = MagicMock()

        view3d.draw_standard_3d_style(mol)

        view3d.plotter.clear.assert_called()
        view3d.plotter.set_background.assert_called_with("#ffffff")
        assert view3d.plotter.add_mesh.call_count >= 1
        view3d.plotter.render.assert_called()

        import numpy as np

        assert hasattr(view3d, "atom_positions_3d")
        assert isinstance(view3d.atom_positions_3d, np.ndarray)
        assert len(view3d.atom_positions_3d) == 3


def test_view_3d_draw_none(app, mock_parser_host):
    """Verify that calling draw with None safely clears the renderer."""
    mock_parser_host.init_manager.settings.update({"background_color": "#000000"})
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)
    view3d.plotter = MagicMock()

    view3d.draw_standard_3d_style(None)

    view3d.plotter.clear.assert_called()
    view3d.plotter.render.assert_called()
    assert view3d.current_mol is None


def test_original_id_mode_hides_rdkit_id_labels(app, mock_parser_host):
    """Original-ID labels should not fall back to RDKit index labels."""
    view3d = _make_view3d(mock_parser_host)
    mock_parser_host.view_3d_manager = view3d
    view3d.plotter = MagicMock()
    view3d.plotter.add_point_labels.return_value = MagicMock()
    view3d.plotter.add_text.return_value = MagicMock()
    view3d.plotter.remove_actor = MagicMock()
    view3d.current_atom_info_labels = None
    view3d.atom_label_legend_names = []
    view3d.atom_info_display_mode = "original_id"
    view3d.atom_index_base = 0
    view3d.atom_positions_3d = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    mol = Chem.RWMol()
    first_idx = mol.AddAtom(Chem.Atom("C"))
    second_idx = mol.AddAtom(Chem.Atom("O"))
    mol.GetAtomWithIdx(first_idx).SetIntProp("_original_atom_id", 42)
    view3d.current_mol = mol

    view3d.show_all_atom_info()

    view3d.plotter.add_point_labels.assert_called_once()
    _, label_texts = view3d.plotter.add_point_labels.call_args.args[:2]
    assert label_texts == ["42"]
    view3d.plotter.add_text.assert_called_once()
    assert view3d.plotter.add_text.call_args.args[0] == "ID"


def test_fit_to_view_empty(app, mock_parser_host):
    """Verify that fit_to_view resets zoom when the scene has no items."""
    view3d = _make_view3d(mock_parser_host)
    mock_parser_host.init_manager.scene.items.return_value = []

    # Mock reset_zoom to verify it gets called
    view3d.reset_zoom = MagicMock()

    view3d.fit_to_view()

    view3d.reset_zoom.assert_called_once()


def test_fit_to_view_with_items(app, mock_parser_host):
    """Verify that fit_to_view calculates bounding box and fits view when scene has items."""
    view3d = _make_view3d(mock_parser_host)

    # Mock items in scene
    mock_item = MagicMock()
    mock_item.isVisible.return_value = True
    from PyQt6.QtCore import QRectF

    mock_item.sceneBoundingRect.return_value = QRectF(10, 10, 50, 50)

    mock_parser_host.init_manager.scene.items.return_value = [mock_item]

    # Mock anchors and fitInView
    view_2d = mock_parser_host.init_manager.view_2d
    view_2d.transformationAnchor.return_value = 1
    view_2d.resizeAnchor.return_value = 2

    view3d.fit_to_view()

    view_2d.setTransformationAnchor.assert_called()
    view_2d.setResizeAnchor.assert_called()
    view_2d.fitInView.assert_called_once()


# ---------------------------------------------------------------------------
# Extended tests (use _make_view3d_ext — normal constructor with mock_host)
# ---------------------------------------------------------------------------


def _make_view3d_ext(mock_host):
    """Create a View3DManager using the normal constructor with a mock host."""
    view3d = View3DManager(mock_host)
    view3d._drawing_3d = False
    view3d.current_3d_style = "ball_and_stick"
    view3d.plotter = MagicMock()
    return view3d


@pytest.fixture
def mock_pv():
    with patch("moleditpy.ui.view_3d_logic.pv") as mock:
        mock_poly = MagicMock()
        mock.PolyData.return_value = mock_poly
        mock_poly.glyph.return_value = MagicMock()
        mock_poly.tube.return_value = MagicMock()

        mock_spline = MagicMock()
        mock.Spline.return_value = mock_spline
        mock_spline.tube.return_value = MagicMock()

        mock.Light.return_value = MagicMock()

        yield mock


def test_add_3d_atom_glyphs_styles(app, mock_parser_host, mock_pv):
    """Verify radii and resolutions for different styles in _add_3d_atom_glyphs."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    sym = ["C"]
    col = np.array([[0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "ball_and_stick", True, mesh_props)
    assert mock_pv.PolyData.call_count >= 1

    mock_poly = mock_pv.PolyData.return_value
    assert "radii" in mock_poly.__setitem__.call_args_list[1][0]
    rad_array = mock_poly.__setitem__.call_args_list[1][0][1]
    assert np.isclose(rad_array[0], 0.51)  # Default VDW for C (1.7 * 0.3)

    mock_pv.PolyData.reset_mock()
    mock_parser_host.init_manager.settings["cpk_atom_scale"] = 1.0
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "cpk", True, mesh_props)
    rad_array = mock_pv.PolyData.return_value.__setitem__.call_args_list[1][0][1]
    assert rad_array[0] > 1.5

    mock_pv.PolyData.reset_mock()
    mock_parser_host.init_manager.settings["stick_bond_radius"] = 0.15
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "stick", True, mesh_props)
    rad_array = mock_pv.PolyData.return_value.__setitem__.call_args_list[1][0][1]
    assert np.isclose(rad_array[0], 0.15)


def test_add_3d_atom_glyphs_stick_split(app, mock_parser_host, mock_pv):
    """Verify that terminal multiple bonds lead to atom splitting in stick mode."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C=C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    sym = ["C", "C"]
    col = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    view3d.atom_positions_3d = np.array([[0.0, 0.0, 0.0], [1.33, 0.0, 0.0]])

    mock_pv.PolyData.reset_mock()
    view3d._add_3d_atom_glyphs(mol, conf, sym, col, "stick", True, mesh_props)

    assert mock_pv.PolyData.call_count >= 2
    new_positions = mock_pv.PolyData.call_args_list[-1][0][0]
    assert len(new_positions) == 4


def test_add_3d_bond_cylinders_basic(app, mock_parser_host, mock_pv):
    """Verify single, double, and triple bond generation."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C#CC=CC")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    col = np.array([[0.5, 0.5, 0.5]] * 5)
    mesh_props = {"smooth_shading": True}

    view3d.plotter.add_mesh = MagicMock()

    view3d._add_3d_bond_cylinders(mol, conf, col, "ball_and_stick", mesh_props)

    assert mock_pv.PolyData.call_count >= 1
    points = mock_pv.PolyData.call_args_list[-1][0][0]
    assert len(points) == 14


def test_add_3d_bond_cylinders_styles(app, mock_parser_host, mock_pv):
    """Check style-dependent factors (radius/offset factors) in bond drawing."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C=C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    col = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])
    mesh_props = {"smooth_shading": True}

    mock_parser_host.init_manager.settings.update(
        {
            "ball_stick_bond_radius": 0.1,
            "ball_stick_double_bond_radius_factor": 0.8,
            "ball_stick_double_bond_offset_factor": 2.0,
        }
    )
    view3d._add_3d_bond_cylinders(mol, conf, col, "ball_and_stick", mesh_props)
    mock_pd = mock_pv.PolyData.return_value
    radii = mock_pd.point_data.__setitem__.call_args_list[0][0][1]
    assert np.allclose(radii, 0.08)

    mock_pv.PolyData.reset_mock()
    mock_pd.point_data.reset_mock()
    mock_parser_host.init_manager.settings.update(
        {
            "stick_bond_radius": 0.15,
            "stick_double_bond_radius_factor": 0.6,
            "stick_double_bond_offset_factor": 1.5,
        }
    )
    view3d._add_3d_bond_cylinders(mol, conf, col, "stick", mesh_props)
    radii = mock_pd.point_data.__setitem__.call_args_list[0][0][1]
    assert np.allclose(radii, 0.09)


def test_add_3d_bond_cylinders_overrides(app, mock_parser_host, mock_pv):
    """Verify plugin bond color overrides."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("CC")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    col = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    mesh_props = {"smooth_shading": True}

    view3d._plugin_bond_color_overrides = {0: "#FF0000"}

    view3d._add_3d_bond_cylinders(mol, conf, col, "stick", mesh_props)

    mock_pd = mock_pv.PolyData.return_value
    colors = mock_pd.cell_data.__setitem__.call_args_list[0][0][1]
    assert np.array_equal(colors[0], [255, 0, 0])


def test_add_3d_aromatic_rings(app, mock_parser_host, mock_pv):
    """Verify aromatic torus generation."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("c1ccccc1")
    AllChem.EmbedMolecule(mol)
    view3d.atom_positions_3d = np.array(
        [list(mol.GetConformer().GetAtomPosition(i)) for i in range(6)]
    )

    mock_parser_host.init_manager.settings["display_aromatic_circles_3d"] = True
    view3d.plotter.add_mesh = MagicMock()

    mesh_props = {"smooth_shading": True}
    view3d._add_3d_aromatic_rings(mol, "ball_and_stick", mesh_props)

    assert mock_pv.Spline.call_count == 1
    assert view3d.plotter.add_mesh.call_count == 1


def test_calculate_double_bond_offset(app, mock_parser_host):
    """Verify neighbor-based plane calculation for double bond offset."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C=C(C)C")
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    bond = mol.GetBondWithIdx(0)

    offset = view3d._calculate_double_bond_offset(mol, bond, conf)
    assert len(offset) == 3
    assert np.isclose(np.linalg.norm(offset), 1.0)


def test_show_ez_labels_3d(app, mock_parser_host):
    """Verify EZ label detection and discrepancy marking."""
    view3d = _make_view3d_ext(mock_parser_host)

    mol = Chem.MolFromSmiles("C/C=C/C")
    AllChem.EmbedMolecule(mol)

    view3d.plotter.add_point_labels = MagicMock()

    view3d.show_ez_labels_3d(mol)

    args, kwargs = view3d.plotter.add_point_labels.call_args
    assert "E" in args[1]

    mock_parser_host.state_manager.data.bonds = {(1, 2): {"stereo": 3}}
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i)

    view3d.plotter.add_point_labels.reset_mock()
    view3d.show_ez_labels_3d(mol)
    args, kwargs = view3d.plotter.add_point_labels.call_args
    assert "?" in args[1]


def test_chiral_labels_logic(app, mock_parser_host, mock_pv):
    """Verify chiral label toggling and update."""
    view3d = _make_view3d_ext(mock_parser_host)
    mock_parser_host.statusBar.return_value = MagicMock()
    view3d.plotter = MagicMock()

    mol = Chem.MolFromSmiles("C[C@H](F)Cl")
    AllChem.EmbedMolecule(mol)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i + 1)

    view3d.plotter = MagicMock()
    view3d.host.view_3d_manager = view3d
    view3d.current_mol = mol

    view3d.plotter.add_point_labels = MagicMock()

    view3d.toggle_chiral_labels_display(True)
    assert view3d.show_chiral_labels is True

    atom_item = MagicMock()
    mock_parser_host.state_manager.data.atoms = {2: {"symbol": "C"}}
    mock_parser_host.init_manager.scene.atom_items[2] = atom_item

    with patch(
        "moleditpy.ui.view_3d_logic.Chem.FindMolChiralCenters", return_value=[(1, "S")]
    ):
        view3d.update_chiral_labels()
        assert atom_item.chiral_label == "S"


def test_color_overrides(app, mock_parser_host, mock_pv):
    """Verify color override API functions."""
    view3d = _make_view3d_ext(mock_parser_host)
    view3d.draw_molecule_3d = MagicMock()
    view3d.current_mol = Chem.MolFromSmiles("C")

    view3d.update_bond_color_override(0, "#FF0000")
    assert view3d._plugin_bond_color_overrides[0] == "#FF0000"
    view3d.draw_molecule_3d.assert_called()

    view3d.update_atom_color_override(0, "#00FF00")
    assert view3d._plugin_color_overrides[0] == "#00FF00"
    assert view3d.draw_molecule_3d.call_count == 2
