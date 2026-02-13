import pytest
import os
from rdkit import Chem
from moleditpy.modules.main_window_compute import MainWindowCompute
from PyQt6.QtCore import Qt
from unittest.mock import MagicMock, patch

class DummyCompute(MainWindowCompute):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.view_2d = host.view_2d
        self.view_3d = host.view_3d
        self.settings = host.settings
        self.current_mol = host.current_mol
        self.active_worker_ids = host.active_worker_ids
        self.halt_ids = host.halt_ids
        self.convert_button = MagicMock()
        self.cleanup_button = MagicMock()
    
    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self): return self._host.statusBar()
    def update_from_rdkit_mol(self, mol): pass
    def create_atom_id_mapping(self): return {}
    def draw_molecule_3d(self, mol): pass
    def update_chiral_labels(self): pass
    def setup_3d_hover(self): pass
    def update_atom_id_menu_text(self): pass
    def update_atom_id_menu_state(self): pass

def test_compute_set_optimization_method(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.set_optimization_method("GAFF_OBABEL")
    assert compute.settings['optimization_method'] == 'GAFF_OBABEL'
    assert compute.statusBar().showMessage.called

def test_compute_halt_logic(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids.add("test_id")
    compute.halt_conversion()
    assert "test_id" in compute.halt_ids
    assert len(compute.active_worker_ids) == 0
    assert compute.statusBar().showMessage.called

def test_on_calculation_finished_basic(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    worker_id = "test_worker"
    compute.active_worker_ids.add(worker_id)
    mol = Chem.MolFromSmiles("C")
    result = (worker_id, mol)
    
    with patch.object(compute, 'draw_molecule_3d') as mock_draw:
        compute.on_calculation_finished(result)
        assert compute.current_mol == mol
        assert worker_id not in compute.active_worker_ids
        assert mock_draw.called

def test_on_calculation_error_basic(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.on_calculation_error("Error Occurred")
    assert compute.statusBar().showMessage.called

def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    """Test manual chemistry problem detection in fallback mode."""
    compute = DummyCompute(mock_parser_host)
    c_id = mock_parser_host.scene.create_atom("C", [0,0,0])
    c_item = mock_parser_host.data.atoms[c_id]['item']
    
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", [i+1, 0, 0])
        h_item = mock_parser_host.data.atoms[h_id]['item']
        mock_parser_host.scene.create_bond(c_item, h_item, bond_order=1)
    
    compute.check_chemistry_problems_fallback()
    assert c_item.has_problem is True
    assert compute.statusBar().showMessage.called
