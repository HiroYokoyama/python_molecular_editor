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
        self._active_calc_threads = []
    
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
    def update_window_title(self): pass
    def setEnabled(self, val): pass
    def _enable_3d_features(self, val): pass
    def _enable_3d_edit_actions(self, val): pass
    def push_undo_state(self): pass

def test_on_calculation_error_stale(mock_parser_host):
    """Test on_calculation_error when the worker is stale (not in active set)."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"new_worker_id"}
    compute.on_calculation_error(("stale_id", "Ignore this error"))
    assert not compute.statusBar().showMessage.called

def test_on_calculation_error_basic(mock_parser_host):
    """Test on_calculation_error for an ACTIVE worker."""
    compute = DummyCompute(mock_parser_host)
    worker_id = "active_id"
    compute.active_worker_ids = {worker_id}
    compute.on_calculation_error((worker_id, "Real Error"))
    assert compute.statusBar().showMessage.called
    assert "Real Error" in compute.statusBar().showMessage.call_args[0][0]

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

def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.data.bonds = {}
    c_id = mock_parser_host.scene.create_atom("C", QPointF(0,0))
    c_item = mock_parser_host.data.atoms[c_id]['item']
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", QPointF(i+1, 0))
        h_item = mock_parser_host.data.atoms[h_id]['item']
        mock_parser_host.scene.create_bond(c_item, h_item, bond_order=1)
    compute.check_chemistry_problems_fallback()
    assert c_item.has_problem is True
    assert compute.statusBar().showMessage.called

def test_trigger_conversion_empty(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.trigger_conversion()
    assert compute.statusBar().showMessage.called

def test_trigger_conversion_with_atoms(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    compute.settings['conversion_target'] = 'all'
    with patch('moleditpy.modules.main_window_compute.CalculationWorker'), \
         patch('moleditpy.modules.main_window_compute.QThread'):
        compute.trigger_conversion()
        assert compute.statusBar().showMessage.called

def test_optimize_3d_structure_logic(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol) 
    compute.current_mol = mol
    compute.optimize_3d_structure()
    assert compute.statusBar().showMessage.called

def test_on_calculation_finished_worker_id_mismatch(mock_parser_host):
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"valid_id"}
    mol = Chem.MolFromSmiles("C")
    result = ("invalid_id", mol)
    compute.on_calculation_finished(result)
    assert compute.current_mol is None

def test_trigger_conversion_chemistry_problems(mock_parser_host):
    """Test trigger_conversion when Chem.DetectChemistryProblems finds issues."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    atom = mol.GetAtomWithIdx(0)
    atom.SetIntProp("_original_atom_id", 1)
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    
    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0
    with patch('rdkit.Chem.DetectChemistryProblems', return_value=[problem]), \
         patch.object(compute.data, 'to_rdkit_mol', return_value=mol):
        compute.trigger_conversion()
        all_messages = [str(call[0][0]) for call in compute.statusBar().showMessage.call_args_list if call[0]]
        assert any("chemistry problem(s) found" in msg for msg in all_messages)

def test_trigger_conversion_sanitize_error(mock_parser_host):
    """Test trigger_conversion when Chem.SanitizeMol fails."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    # MUST populate atoms to avoid empty trigger return
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    
    with patch('rdkit.Chem.DetectChemistryProblems', return_value=[]), \
         patch('rdkit.Chem.SanitizeMol', side_effect=ValueError("Sanitize failed")), \
         patch.object(compute.data, 'to_rdkit_mol', return_value=mol):
        compute.trigger_conversion()
        all_messages = [str(call[0][0]) for call in compute.statusBar().showMessage.call_args_list if call[0]]
        assert any("Error: Invalid chemical structure." in msg for msg in all_messages)

def test_trigger_conversion_multiple_frags(mock_parser_host):
    """Test trigger_conversion with multiple fragments."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C.C")
    # MUST populate atoms to avoid empty trigger return
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}, 2: {'symbol': 'C', 'item': MagicMock()}}
    
    with patch('rdkit.Chem.DetectChemistryProblems', return_value=[]), \
         patch('rdkit.Chem.SanitizeMol'), \
         patch.object(compute.data, 'to_rdkit_mol', return_value=mol), \
         patch('moleditpy.modules.main_window_compute.CalculationWorker'), \
         patch('moleditpy.modules.main_window_compute.QThread'):
        compute.trigger_conversion()
        all_messages = [str(call[0][0]) for call in compute.statusBar().showMessage.call_args_list if call[0]]
        assert any("collision detection" in msg for msg in all_messages)

def test_on_calculation_finished_single_mol_legacy(mock_parser_host):
    """Test on_calculation_finished with a single mol (legacy result format)."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    compute.on_calculation_finished(mol)
    assert compute.current_mol == mol

def test_on_calculation_error_legacy_payload(mock_parser_host):
    """Test on_calculation_error with a string (legacy error format)."""
    compute = DummyCompute(mock_parser_host)
    compute.on_calculation_error("Fatal Error")
    assert any("Fatal Error" in str(call[0][0]) for call in compute.statusBar().showMessage.call_args_list if call[0])
