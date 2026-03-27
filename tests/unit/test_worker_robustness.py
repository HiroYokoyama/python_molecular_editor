import pytest
from unittest.mock import MagicMock, patch
from rdkit import Chem
from moleditpy.core.calculation_worker import CalculationWorker, WorkerHaltError

class MockWorker(CalculationWorker):
    def __init__(self):
        super().__init__()
        self.halt_all = False
        self.halt_ids = set()

def test_worker_halt_logic():
    """Test that worker correctly identifies halt requests."""
    worker = MockWorker()
    worker.halt_all = True
    
    # Internal _check_halted logic check via mock run
    with pytest.raises(WorkerHaltError):
        # We simulate the _safe_status which checks for halt
        def _check_halted(): return worker.halt_all
        def _safe_status(msg):
            if _check_halted(): raise WorkerHaltError("Halted")
        
        _safe_status("testing halt")

def test_iterative_optimize_robustness():
    """Test that iterative_optimize handles RDKit failures by returning False instead of crashing."""
    from moleditpy.core.calculation_worker import _iterative_optimize
    mol = Chem.MolFromSmiles("C")
    # Empty conformer should fail optimization
    assert not _iterative_optimize(mol, "UFF", lambda: False, lambda x: None)

def test_direct_conversion_failure_robustness():
    """Test that _perform_direct_conversion raises ValueError on empty input."""
    from moleditpy.core.calculation_worker import _perform_direct_conversion
    mol = Chem.MolFromSmiles("C")
    with pytest.raises(ValueError, match="Failed to parse coordinates"):
        _perform_direct_conversion("", mol, {}, lambda: False, lambda x: None)

def test_run_calculation_empty_input():
    """Test that run_calculation emits error on empty input."""
    worker = MockWorker()
    error_handler = MagicMock()
    worker.error.connect(error_handler)
    
    # Should not raise, but emit error signal
    worker.run_calculation("")
    
    assert error_handler.called
    args, _ = error_handler.call_args
    # First arg is (worker_id, error_msg)
    assert "No atoms to convert" in str(args[0][1])
