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


# ---------------------------------------------------------------------------
# plotter property / cleanup / style switching
# ---------------------------------------------------------------------------


def test_plotter_property_roundtrip(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    sentinel = MagicMock()
    view3d.plotter = sentinel
    assert view3d.plotter is sentinel
    view3d.plotter = None
    assert view3d.plotter is None


def test_cleanup_clears_and_closes_plotter(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    plotter = MagicMock()
    view3d.plotter = plotter
    view3d.current_mol = MagicMock()
    view3d.cleanup()
    plotter.clear.assert_called_once()
    plotter.close.assert_called_once()
    assert view3d.current_mol is None


def test_cleanup_without_plotter(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.plotter = None
    view3d.current_mol = MagicMock()
    view3d.cleanup()  # must not raise
    assert view3d.current_mol is None


def test_set_3d_style_same_style_is_noop(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "stick"
    view3d.set_3d_style("stick")
    mock_parser_host.edit_3d_manager.clear_3d_selection.assert_not_called()


def test_set_3d_style_disables_edit_modes_and_redraws(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "ball_and_stick"
    view3d.current_mol = MagicMock()
    view3d.draw_molecule_3d = MagicMock()
    mock_parser_host.edit_3d_manager.measurement_mode = True
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = True

    view3d.set_3d_style("stick")

    mock_parser_host.init_manager.measurement_action.setChecked.assert_called_with(
        False
    )
    mock_parser_host.edit_3d_manager.toggle_measurement_mode.assert_called_once_with(
        False
    )
    mock_parser_host.ui_manager.toggle_3d_edit_mode.assert_called_once_with(False)
    mock_parser_host.edit_3d_manager.clear_3d_selection.assert_called_once()
    assert view3d.current_3d_style == "stick"
    view3d.draw_molecule_3d.assert_called_once_with(view3d.current_mol)


def test_set_3d_style_without_molecule_skips_redraw(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "ball_and_stick"
    view3d.current_mol = None
    view3d.draw_molecule_3d = MagicMock()
    mock_parser_host.edit_3d_manager.measurement_mode = False
    mock_parser_host.edit_3d_manager.is_3d_edit_mode = False

    view3d.set_3d_style("stick")

    view3d.draw_molecule_3d.assert_not_called()


def test_draw_molecule_3d_uses_plugin_custom_style(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "plugin_style"
    handler = MagicMock()
    mock_parser_host.plugin_manager.custom_3d_styles = {
        "plugin_style": {"callback": handler}
    }
    view3d.draw_standard_3d_style = MagicMock()
    mol = MagicMock()

    view3d.draw_molecule_3d(mol)

    handler.assert_called_once_with(mock_parser_host, mol)
    view3d.draw_standard_3d_style.assert_not_called()


def test_draw_molecule_3d_broken_plugin_falls_back_to_standard(
    mock_parser_host,
):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "plugin_style"
    mock_parser_host.plugin_manager.custom_3d_styles = {
        "plugin_style": {"callback": MagicMock(side_effect=RuntimeError("bug"))}
    }
    view3d.draw_standard_3d_style = MagicMock()
    mol = MagicMock()

    view3d.draw_molecule_3d(mol)

    view3d.draw_standard_3d_style.assert_called_once_with(mol)


def test_draw_molecule_3d_standard_when_style_not_custom(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_3d_style = "ball_and_stick"
    mock_parser_host.plugin_manager.custom_3d_styles = {}
    view3d.draw_standard_3d_style = MagicMock()
    mol = MagicMock()

    view3d.draw_molecule_3d(mol)

    view3d.draw_standard_3d_style.assert_called_once_with(mol)


# ---------------------------------------------------------------------------
# Chiral labels
# ---------------------------------------------------------------------------


def test_toggle_chiral_labels_on_redraws_and_reports(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_mol = MagicMock()
    view3d.draw_molecule_3d = MagicMock()

    view3d.toggle_chiral_labels_display(True)

    assert view3d.show_chiral_labels is True
    view3d.draw_molecule_3d.assert_called_once()
    assert any(
        "Chiral labels" in str(c.args[0])
        for c in mock_parser_host.statusBar().showMessage.call_args_list
    )


def test_toggle_chiral_labels_off_reports_disabled(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_mol = None
    view3d.toggle_chiral_labels_display(False)
    assert view3d.show_chiral_labels is False
    assert any(
        "Chiral labels disabled." in str(c.args[0])
        for c in mock_parser_host.statusBar().showMessage.call_args_list
    )


def test_update_chiral_labels_disabled_clears_items(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.show_chiral_labels = False
    item = MagicMock()
    item.chiral_label = "R"
    mock_parser_host.init_manager.scene.atom_items = {1: item}

    view3d.update_chiral_labels()

    assert item.chiral_label is None
    mock_parser_host.init_manager.scene.update.assert_called()


def test_update_chiral_labels_assigns_r_or_s(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.show_chiral_labels = True

    mol = Chem.MolFromSmiles("[C@H](F)(Cl)Br")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    for atom in mol.GetAtoms():
        atom.SetIntProp("_original_atom_id", atom.GetIdx() + 100)
    view3d.current_mol = mol

    items = {idx + 100: MagicMock(chiral_label=None) for idx in range(mol.GetNumAtoms())}
    mock_parser_host.init_manager.scene.atom_items = items

    view3d.update_chiral_labels()

    labels = {i.chiral_label for i in items.values()}
    assert labels & {"R", "S"}, f"expected a chiral label, got {labels}"


# ---------------------------------------------------------------------------
# Atom info display
# ---------------------------------------------------------------------------


def test_toggle_atom_info_same_mode_turns_off(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.atom_info_display_mode = "rdkit_index"
    view3d.clear_all_atom_info_labels = MagicMock()

    view3d.toggle_atom_info_display("rdkit_index")

    assert view3d.atom_info_display_mode is None
    mock_parser_host.init_manager.show_index_action.setChecked.assert_called_with(
        False
    )
    mock_parser_host.init_manager.atom_index_base_menu.setEnabled.assert_called_with(
        False
    )
    assert any(
        "Atom info display disabled." in str(c.args[0])
        for c in mock_parser_host.statusBar().showMessage.call_args_list
    )


def test_toggle_atom_info_new_index_mode_enables_base_menu(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.atom_info_display_mode = None
    view3d.clear_all_atom_info_labels = MagicMock()
    view3d.show_all_atom_info = MagicMock()

    view3d.toggle_atom_info_display("rdkit_index")

    assert view3d.atom_info_display_mode == "rdkit_index"
    mock_parser_host.init_manager.show_index_action.setChecked.assert_called_with(
        True
    )
    mock_parser_host.init_manager.show_atom_symbol_action.setChecked.assert_called_with(
        False
    )
    mock_parser_host.init_manager.atom_index_base_menu.setEnabled.assert_called_with(
        True
    )
    view3d.show_all_atom_info.assert_called_once()
    view3d.plotter.render.assert_called_once()


def test_toggle_atom_info_symbol_mode_disables_base_menu(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.atom_info_display_mode = None
    view3d.clear_all_atom_info_labels = MagicMock()
    view3d.show_all_atom_info = MagicMock()

    view3d.toggle_atom_info_display("symbol")

    mock_parser_host.init_manager.atom_index_base_menu.setEnabled.assert_called_with(
        False
    )
    assert any(
        "Displaying: Element Symbol" in str(c.args[0])
        for c in mock_parser_host.statusBar().showMessage.call_args_list
    )


def test_set_atom_index_base_syncs_checkmarks(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.atom_info_display_mode = None
    view3d.set_atom_index_base(1)
    assert view3d.atom_index_base == 1
    mock_parser_host.init_manager.atom_index_base_0_action.setChecked.assert_called_with(
        False
    )
    mock_parser_host.init_manager.atom_index_base_1_action.setChecked.assert_called_with(
        True
    )


def test_set_atom_index_base_refreshes_active_index_labels(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.atom_info_display_mode = "xyz_index"
    view3d.clear_all_atom_info_labels = MagicMock()
    view3d.show_all_atom_info = MagicMock()

    view3d.set_atom_index_base(0)

    view3d.clear_all_atom_info_labels.assert_called_once()
    view3d.show_all_atom_info.assert_called_once()
    view3d.plotter.render.assert_called_once()


# ---------------------------------------------------------------------------
# Molecule kind detection + menu state
# ---------------------------------------------------------------------------


def test_is_xyz_derived_molecule_detects_property(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol.GetAtomWithIdx(0).SetProp("xyz_unique_id", "7")
    view3d.current_mol = mol
    assert view3d.is_xyz_derived_molecule() is True


def test_is_xyz_derived_molecule_false_cases(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_mol = None
    assert view3d.is_xyz_derived_molecule() is False
    view3d.current_mol = Chem.MolFromSmiles("C")
    assert view3d.is_xyz_derived_molecule() is False


def test_has_original_atom_ids(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    view3d.current_mol = mol
    assert view3d.has_original_atom_ids() is False
    mol.GetAtomWithIdx(1).SetIntProp("_original_atom_id", 5)
    assert view3d.has_original_atom_ids() is True


def test_update_atom_id_menu_state_enables_matching_actions(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol.GetAtomWithIdx(0).SetProp("xyz_unique_id", "1")
    view3d.current_mol = mol

    view3d.update_atom_id_menu_state()

    mock_parser_host.init_manager.show_original_id_action.setEnabled.assert_called_with(
        False
    )
    mock_parser_host.init_manager.show_xyz_index_action.setEnabled.assert_called_with(
        True
    )


# ---------------------------------------------------------------------------
# Zoom / fit / hover
# ---------------------------------------------------------------------------


def test_zoom_in_and_out_scale_view(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.zoom_in()
    mock_parser_host.init_manager.view_2d.scale.assert_called_with(1.2, 1.2)
    view3d.zoom_out()
    mock_parser_host.init_manager.view_2d.scale.assert_called_with(
        1 / 1.2, 1 / 1.2
    )


def test_reset_zoom_sets_075_transform(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.reset_zoom()
    transform = mock_parser_host.init_manager.view_2d.setTransform.call_args.args[0]
    assert transform.m11() == 0.75
    assert transform.m22() == 0.75


def test_fit_to_view_empty_scene_resets_zoom(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    mock_parser_host.init_manager.scene.items.return_value = []
    view3d.reset_zoom = MagicMock()
    view3d.fit_to_view()
    view3d.reset_zoom.assert_called_once()


def test_fit_to_view_fits_visible_items(mock_parser_host):
    from PyQt6.QtCore import QRectF

    view3d = _make_view3d(mock_parser_host)
    item = MagicMock()
    item.isVisible.return_value = True
    item.sceneBoundingRect.return_value = QRectF(0, 0, 100, 50)
    mock_parser_host.init_manager.scene.items.return_value = [item]
    view3d.reset_zoom = MagicMock()

    view3d.fit_to_view()

    view3d.reset_zoom.assert_not_called()
    mock_parser_host.init_manager.view_2d.fitInView.assert_called_once()
    padded = mock_parser_host.init_manager.view_2d.fitInView.call_args.args[0]
    assert padded.width() == pytest.approx(110.0)  # 10% padding
    assert padded.height() == pytest.approx(55.0)


def test_setup_3d_hover_only_refreshes_when_mode_active(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.show_all_atom_info = MagicMock()
    view3d.atom_info_display_mode = None
    view3d.setup_3d_hover()
    view3d.show_all_atom_info.assert_not_called()
    view3d.atom_info_display_mode = "symbol"
    view3d.setup_3d_hover()
    view3d.show_all_atom_info.assert_called_once()


# ---------------------------------------------------------------------------
# Plugin color overrides
# ---------------------------------------------------------------------------


def test_update_atom_color_override_set_and_clear(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_mol = MagicMock()
    view3d.draw_molecule_3d = MagicMock()

    view3d.update_atom_color_override(2, "#ff0000")
    assert view3d._plugin_color_overrides == {2: "#ff0000"}
    view3d.update_atom_color_override(2, None)
    assert view3d._plugin_color_overrides == {}
    assert view3d.draw_molecule_3d.call_count == 2


def test_update_bond_color_override_set_and_clear(mock_parser_host):
    view3d = _make_view3d(mock_parser_host)
    view3d.current_mol = None
    view3d.draw_molecule_3d = MagicMock()

    view3d.update_bond_color_override(1, "#00ff00")
    assert view3d._plugin_bond_color_overrides == {1: "#00ff00"}
    view3d.update_bond_color_override(1, None)
    assert view3d._plugin_bond_color_overrides == {}
    view3d.draw_molecule_3d.assert_not_called()  # no molecule loaded
