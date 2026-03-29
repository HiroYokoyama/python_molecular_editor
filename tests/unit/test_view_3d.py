import pytest
from unittest.mock import MagicMock, patch
from rdkit import Chem
from rdkit.Chem import AllChem


def _make_view3d(mock_host):
    """Create a View3DManager instance with the given mock host."""
    from moleditpy.ui.view_3d_logic import View3DManager
    view3d = View3DManager.__new__(View3DManager)
    view3d.host = mock_host
    view3d._drawing_3d = False
    view3d.current_3d_style = "Ball and Stick"
    return view3d


def test_view_3d_draw_standard_3d_style(app, mock_parser_host):
    """Verify that draw_standard_3d_style clears the plotter and constructs the correct VTK meshes."""
    mock_parser_host.plotter = MagicMock()
    mock_parser_host.settings = {
        "projection_mode": "Perspective",
        "background_color": "#ffffff",
        "display_kekule_3d": False
    }
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)

    mol = Chem.MolFromSmiles("CCO")
    AllChem.EmbedMolecule(mol)

    with patch('moleditpy.ui.view_3d_logic.pv') as mock_pv:
        mock_pv.PolyData.return_value.glyph.return_value = MagicMock()
        mock_pv.PolyData.return_value.tube.return_value = MagicMock()
        mock_pv.Light.return_value = MagicMock()

        view3d.draw_standard_3d_style(mol)

        mock_parser_host.plotter.clear.assert_called()
        mock_parser_host.plotter.set_background.assert_called_with("#ffffff")
        assert mock_parser_host.plotter.add_mesh.call_count >= 1
        mock_parser_host.plotter.render.assert_called()

        import numpy as np
        assert hasattr(view3d, "atom_positions_3d")
        assert isinstance(view3d.atom_positions_3d, np.ndarray)
        assert len(view3d.atom_positions_3d) == 3


def test_view_3d_draw_none(app, mock_parser_host):
    """Verify that calling draw with None safely clears the renderer."""
    mock_parser_host.plotter = MagicMock()
    mock_parser_host.settings = {"background_color": "#000000"}
    mock_parser_host.edit_3d_manager = MagicMock()

    view3d = _make_view3d(mock_parser_host)

    view3d.draw_standard_3d_style(None)

    mock_parser_host.plotter.clear.assert_called()
    mock_parser_host.plotter.render.assert_called()
    assert mock_parser_host.current_mol is None
