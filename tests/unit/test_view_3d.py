import pytest
import types
from unittest.mock import MagicMock, patch
from rdkit import Chem
from rdkit.Chem import AllChem

def test_view_3d_draw_standard_3d_style(app, mock_parser_host):
    """Verify that draw_standard_3d_style clears the plotter and constructs the correct VTK meshes."""
    from moleditpy.ui.view_3d import MainWindowView3d
    
    # Mixin the methods to our mock host
    for attr in dir(MainWindowView3d):
        if not attr.startswith('__') and callable(getattr(MainWindowView3d, attr)):
            setattr(mock_parser_host, attr, types.MethodType(getattr(MainWindowView3d, attr), mock_parser_host))
    
    # Provide a simple 3D molecule
    mol = Chem.MolFromSmiles("CCO")
    AllChem.EmbedMolecule(mol)
    
    mock_parser_host._drawing_3d = False
    mock_parser_host.statusBar = MagicMock()
    mock_parser_host.measurement_mode = False
    mock_parser_host.plotter = MagicMock()
    mock_parser_host.settings = {
        "projection_mode": "Perspective",
        "background_color": "#ffffff",
        "display_kekule_3d": False
    }
    
    # We must patch pyvista PolyData so it doesn't try to open X11/OpenGL contexts
    with patch('moleditpy.ui.view_3d.pv') as mock_pv:
        mock_pv.PolyData.return_value.glyph.return_value = MagicMock()
        mock_pv.PolyData.return_value.tube.return_value = MagicMock()
        mock_pv.Light.return_value = MagicMock()
        
        mock_parser_host.draw_standard_3d_style(mol)
        
        # Check high level rendering calls
        mock_parser_host.plotter.clear.assert_called()
        mock_parser_host.plotter.set_background.assert_called_once_with("#ffffff")
        assert mock_parser_host.plotter.add_mesh.call_count >= 1
        mock_parser_host.plotter.render.assert_called()
        
        # Ensure our properties were created correctly based on the RDKit Mol
        import numpy as np
        assert hasattr(mock_parser_host, "atom_positions_3d")
        assert isinstance(mock_parser_host.atom_positions_3d, np.ndarray)
        assert len(mock_parser_host.atom_positions_3d) == 3

def test_view_3d_draw_none(app, mock_parser_host):
    mock_parser_host._drawing_3d = False

    """Verify that calling draw with None safely cleats the renderer."""
    from moleditpy.ui.view_3d import MainWindowView3d
    for attr in dir(MainWindowView3d):
        if not attr.startswith('__') and callable(getattr(MainWindowView3d, attr)):
            setattr(mock_parser_host, attr, types.MethodType(getattr(MainWindowView3d, attr), mock_parser_host))

    mock_parser_host.plotter = MagicMock()
    mock_parser_host.settings = {"background_color": "#000000"}

    mock_parser_host.draw_standard_3d_style(None)
    
    mock_parser_host.plotter.clear.assert_called()
    mock_parser_host.plotter.render.assert_called()
    assert mock_parser_host.current_mol is None
