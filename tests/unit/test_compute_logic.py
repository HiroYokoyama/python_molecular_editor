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
        self.active_worker_ids = set() # Ensure it's a set
        self.halt_ids = set()          # Ensure it's a set
        self.convert_button = MagicMock()
        self.cleanup_button = MagicMock()
        self.optimize_3d_button = MagicMock()
        self.molecule_3d_action = MagicMock()
        self.is_xyz_derived = False
    
    def __getattr__(self, name):
        # Avoid recursion for standard members
        if name in ('_host', 'data', 'scene', 'view_2d', 'view_3d', 'settings', 'current_mol'):
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
    # Payload should be (worker_id, error_msg)
    compute.on_calculation_error(("worker_1", "Error Occurred"))
    assert compute.statusBar().showMessage.called

def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    """Test manual chemistry problem detection in fallback mode."""
    compute = DummyCompute(mock_parser_host)
    c_id = mock_parser_host.scene.create_atom("C", QPointF(0,0))
    c_item = mock_parser_host.data.atoms[c_id]['item']
    
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", QPointF(i+1, 0))
        h_item = mock_parser_host.data.atoms[h_id]['item']
        mock_parser_host.scene.create_bond(c_item, h_item, bond_order=1)
    
    compute.check_chemistry_problems_fallback()
    assert c_item.has_problem is True
    assert compute.statusBar().showMessage.called

def test_trigger_conversion_no_atoms(mock_parser_host):
    """Test trigger_conversion when there are no atoms."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.trigger_conversion()
    assert compute.statusBar().showMessage.called
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "No atoms" in msg

def test_trigger_conversion_all_mode(mock_parser_host):
    """Test trigger_conversion with 'all' mode."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    compute.settings['conversion_target'] = 'all'
    
    with patch('moleditpy.modules.main_window_compute.CalculationWorker') as mock_worker, \
         patch('PyQt6.QtCore.QThread') as mock_thread:
        compute.trigger_conversion()
        # It should at least attempt to create a worker
        assert mock_worker.called

def test_optimize_3d_structure_no_mol(mock_parser_host):
    """Test optimize_3d_structure with no current molecule."""
    compute = DummyCompute(mock_parser_host)
    compute.current_mol = None
    compute.optimize_3d_structure()
    assert compute.statusBar().showMessage.called
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "No 3D structure" in msg

def test_optimize_3d_structure_with_selection(mock_parser_host):
    """Test optimize_3d_structure when some atoms are selected."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    compute.current_mol = mol
    
    # Mock selection
    mock_item = MagicMock()
    mock_item.isSelected.return_value = True
    compute.data.atoms = {1: {'item': mock_item, 'rdkit_idx': 0}}
    
    with patch('moleditpy.modules.main_window_compute.CalculationWorker') as mock_worker, \
         patch('PyQt6.QtCore.QThread') as mock_thread:
        compute.optimize_3d_structure()
        assert mock_worker.called

def test_on_calculation_finished_worker_id_mismatch(mock_parser_host):
    """Test on_calculation_finished with a worker_id not in active_worker_ids."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"valid_id"}
    mol = Chem.MolFromSmiles("C")
    result = ("invalid_id", mol)
    
    compute.on_calculation_finished(result)
    assert compute.current_mol is None

def test_on_calculation_error_halted(mock_parser_host):
    """Test on_calculation_error when the worker was halted."""
    compute = DummyCompute(mock_parser_host)
    compute.halt_ids = {"id123"}
    compute.on_calculation_error(("id123", "Halted"))
    assert compute.statusBar().showMessage.called
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "halted" in msg.lower()

def test_trigger_conversion_no_atoms(mock_parser_host):
    """Test trigger_conversion when there are no atoms."""
    compute = DummyCompute(mock_parser_host)
    mock_parser_host.data.atoms = {}
    compute.trigger_conversion()
    assert compute.statusBar().showMessage.called
    assert "No atoms to convert" in compute.statusBar().showMessage.call_args[0][0]

def test_trigger_conversion_all_mode(mock_parser_host):
    """Test trigger_conversion with 'all' mode."""
    compute = DummyCompute(mock_parser_host)
    mock_parser_host.data.atoms = {1: {'symbol': 'C', 'item': MagicMock()}}
    compute.settings['conversion_target'] = 'all'
    
    with patch('moleditpy.modules.main_window_compute.CalculationWorker') as mock_worker, \
         patch('PyQt6.QtCore.QThread') as mock_thread:
        compute.trigger_conversion()
        assert mock_worker.called
        assert mock_thread.called

def test_optimize_3d_structure_no_mol(mock_parser_host):
    """Test optimize_3d_structure with no current molecule."""
    compute = DummyCompute(mock_parser_host)
    compute.current_mol = None
    compute.optimize_3d_structure()
    assert compute.statusBar().showMessage.called
    assert "No 3D structure" in compute.statusBar().showMessage.call_args[0][0]

def test_optimize_3d_structure_with_selection(mock_parser_host):
    """Test optimize_3d_structure when some atoms are selected."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    compute.current_mol = mol
    
    # Mock selection
    mock_item = MagicMock()
    mock_item.isSelected.return_value = True
    compute._host.data.atoms = {1: {'item': mock_item, 'rdkit_idx': 0}}
    
    with patch('moleditpy.modules.main_window_compute.CalculationWorker') as mock_worker:
        compute.optimize_3d_structure()
        assert mock_worker.called
        # Check if constraints were passed (simplified check)
        args, kwargs = mock_worker.call_args
        assert 'options' in kwargs
        # The logic adds constraints for UNSELECTED atoms if some are selected

def test_on_calculation_finished_worker_id_mismatch(mock_parser_host):
    """Test on_calculation_finished with a worker_id not in active_worker_ids."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"valid_id"}
    mol = Chem.MolFromSmiles("C")
    result = ("invalid_id", mol)
    
    compute.on_calculation_finished(result)
    # Should not update anything significantly if ID is unknown/halted
    assert compute.current_mol is None

def test_on_calculation_error_halted(mock_parser_host):
    """Test on_calculation_error when the worker was halted."""
    compute = DummyCompute(mock_parser_host)
    compute.halt_ids = {123}
    compute.on_calculation_error((123, "Halted message"))
    # Should show "Conversion halted"
    assert "halted" in compute.statusBar().showMessage.call_args[0][0].lower()
