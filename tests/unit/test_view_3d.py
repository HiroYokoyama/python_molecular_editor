from unittest.mock import MagicMock, patch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def _make_view3d(mock_host):
    """Create a View3DManager instance with the given mock host."""
    from moleditpy.ui.view_3d_logic import View3DManager

    view3d = View3DManager.__new__(View3DManager)
    view3d.host = mock_host
    view3d._drawing_3d = False
    view3d.current_3d_style = "Ball and Stick"
    view3d.atom_info_display_mode = None
    return view3d


def test_view_3d_draw_standard_3d_style(app, mock_parser_host):
    """Verify that draw_standard_3d_style clears the plotter and constructs the correct VTK meshes."""
    mock_parser_host.view_3d_manager.plotter = MagicMock()
    mock_parser_host.init_manager.settings.update(
        {
            "projection_mode": "Perspective",
            "background_color": "#ffffff",
            "display_kekule_3d": False,
        }
    )
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("CCO")
    AllChem.EmbedMolecule(mol)

    with patch("moleditpy.ui.view_3d_logic.pv") as mock_pv:
        mock_pv.PolyData.return_value.glyph.return_value = MagicMock()
        mock_pv.PolyData.return_value.tube.return_value = MagicMock()
        mock_pv.Light.return_value = MagicMock()

        view3d.draw_standard_3d_style(mol)

        mock_parser_host.view_3d_manager.plotter.clear.assert_called()
        mock_parser_host.view_3d_manager.plotter.set_background.assert_called_with(
            "#ffffff"
        )
        assert mock_parser_host.view_3d_manager.plotter.add_mesh.call_count >= 1
        mock_parser_host.view_3d_manager.plotter.render.assert_called()

        import numpy as np

        assert hasattr(view3d, "atom_positions_3d")
        assert isinstance(view3d.atom_positions_3d, np.ndarray)
        assert len(view3d.atom_positions_3d) == 3


def test_view_3d_draw_none(app, mock_parser_host):
    """Verify that calling draw with None safely clears the renderer."""
    mock_parser_host.view_3d_manager.plotter = MagicMock()
    mock_parser_host.init_manager.settings.update({"background_color": "#000000"})
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)

    view3d.draw_standard_3d_style(None)

    mock_parser_host.view_3d_manager.plotter.clear.assert_called()
    mock_parser_host.view_3d_manager.plotter.render.assert_called()
    assert mock_parser_host.view_3d_manager.current_mol is None


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
