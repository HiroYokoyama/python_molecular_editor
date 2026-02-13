import pytest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_compute import MainWindowCompute
from PyQt6.QtCore import Qt, QPointF
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
        self.active_worker_ids = set() 
        self.halt_ids = set()          
        self.convert_button = MagicMock()
        self.cleanup_button = MagicMock()
        self.optimize_3d_button = MagicMock()
        self.export_button = MagicMock()
        self.molecule_3d_action = MagicMock()
        self.analysis_action = MagicMock()
        self.edit_3d_action = MagicMock()
        self.is_xyz_derived = False
        self.plotter = MagicMock()
        self.waiting_worker_id = None
        self.optimization_method = 'MMFF_RDKIT'
        self.opt3d_method_labels = {'MMFF_RDKIT': 'MMFF94s (RDKit)', 'UFF_RDKIT': 'UFF (RDKit)'}
    
    def __getattr__(self, name):
        if name in ('_host', 'data', 'scene', 'view_2d', 'view_3d', 'settings', 'current_mol', 'active_worker_ids', 'halt_ids', 'convert_button', 'cleanup_button', 'optimize_3d_button', 'export_button', 'molecule_3d_action', 'analysis_action', 'edit_3d_action', 'is_xyz_derived', 'plotter', 'waiting_worker_id', 'optimization_method', 'opt3d_method_labels'):
            return object.__getattribute__(self, name)
        return getattr(self._host, name)

    def statusBar(self): return self._host.statusBar()
    def update_from_rdkit_mol(self, mol): pass
    def create_atom_id_mapping(self): return {}
    def draw_molecule_3d(self, mol): pass
    def update_chiral_labels(self): pass
    def setup_3d_hover(self): pass
    def update_atom_id_menu_text(self): pass
    def update_atom_id_menu_state(self): pass
    def update_window_title(self): pass
    def setEnabled(self, val): pass
    def _enable_3d_features(self, val): pass
    def _enable_3d_edit_actions(self, val): pass
    def push_undo_state(self): pass

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

def test_on_calculation_error_basic(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    worker_id = "worker_1"
    compute.active_worker_ids.add(worker_id)
    compute.on_calculation_error((worker_id, "Error Occurred"))
    assert compute.statusBar().showMessage.called
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "Error: Error Occurred" in msg

def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    """Test manual chemistry problem detection in fallback mode."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.data.bonds = {}
    
    c_id = mock_parser_host.scene.create_atom("C", QPointF(0,0))
    c_entry = mock_parser_host.data.atoms[c_id]
    c_item = c_entry['item']
    
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", QPointF(i+1, 0))
        h_entry = mock_parser_host.data.atoms[h_id]
        h_item = h_entry['item']
        mock_parser_host.scene.create_bond(c_item, h_item, bond_order=1)
    
    compute.check_chemistry_problems_fallback()
    assert c_item.has_problem is True
    assert compute.statusBar().showMessage.called

def test_trigger_conversion_empty(mock_parser_host):
    """Test trigger_conversion when there are no atoms."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.trigger_conversion()
    assert compute.statusBar().showMessage.called

def test_trigger_conversion_with_atoms(mock_parser_host):
    """Test trigger_conversion with atoms."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    compute.settings['conversion_target'] = 'all'
    
    with patch('moleditpy.modules.main_window_compute.CalculationWorker') as mock_worker, \
         patch('PyQt6.QtCore.QThread') as mock_thread:
        compute.trigger_conversion()
        assert mock_worker.called

def test_optimize_3d_structure_logic(mock_parser_host):
    """Test optimize_3d_structure basic flow."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol) 
    compute.current_mol = mol
    
    # Just verify it shows a status message and doesn't crash
    compute.optimize_3d_structure()
    assert compute.statusBar().showMessage.called
    # It should show "Optimizing..." initially
    
def test_on_calculation_finished_worker_id_mismatch(mock_parser_host):
    """Test on_calculation_finished with a worker_id not in active_worker_ids."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"valid_id"}
    mol = Chem.MolFromSmiles("C")
    result = ("invalid_id", mol)
    
    compute.on_calculation_finished(result)
    assert compute.current_mol is None

def test_on_calculation_error_stale(mock_parser_host):
    """Test on_calculation_error when the worker is stale (not in active set)."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"current_id"}
    compute.on_calculation_error(("stale_id", "Some Error"))
    assert not compute.statusBar().showMessage.called
