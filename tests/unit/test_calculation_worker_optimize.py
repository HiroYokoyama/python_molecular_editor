import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.calculation_worker import CalculationWorker
from unittest.mock import patch

# Helper to capture signal emissions
class SignalCaptor:
    def __init__(self):
        self.emitted_values = []

    def capture(self, *args):
        if len(args) == 1:
            self.emitted_values.append(args[0])
        else:
            self.emitted_values.append(args)

@pytest.fixture
def worker():
    return CalculationWorker()

def test_optimize_only_mmff94s(worker):
    """Test optimize_only mode with MMFF94s."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "MMFF94s_RDKIT",
        "worker_id": 1
    }

    worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result

    assert res_mol.HasProp("_pme_optimization_method")
    assert res_mol.GetProp("_pme_optimization_method").upper() == "MMFF94S_RDKIT"

def test_optimize_only_uff(worker):
    """Test optimize_only mode with UFF."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "UFF_RDKIT",
        "worker_id": 2
    }

    worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result

    assert res_mol.GetProp("_pme_optimization_method") == "UFF_RDKIT"

def test_collision_avoidance_trigger(worker):
    """Test that collision avoidance is called in direct mode."""
    # Use two separate molecules (fragments)
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    # Place them very close
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (0, 0, 0))
    mol.AddConformer(conf)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    # direct mode triggers collision avoidance ONLY IF do_optimize is True
    options = {"conversion_mode": "direct", "do_optimize": True}

    # Mock _iterative_optimize to avoid actual optimization and just check collision avoidance
    with patch("moleditpy.ui.calculation_worker._iterative_optimize", return_value=True):
        worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    res_mol = finish_captor.emitted_values[0]
    if isinstance(res_mol, tuple): res_mol = res_mol[1]

    # Atoms should have been moved apart
    conf = res_mol.GetConformer()
    p1 = conf.GetAtomPosition(0)
    p2 = conf.GetAtomPosition(1)
    dist = (p1.x-p2.x)**2 + (p1.y-p2.y)**2 + (p1.z-p2.z)**2
    assert dist > 0.01

def test_iterative_optimize_halt(worker):
    """Test that iterative optimization respects halt signals."""
    mol = Chem.MolFromSmiles("CCCC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    worker.halt_ids = {1}

    from moleditpy.ui.calculation_worker import _iterative_optimize, WorkerHaltError

    def check_halted():
        return 1 in worker.halt_ids

    def safe_status(msg):
        pass

    with pytest.raises(WorkerHaltError):
        _iterative_optimize(mol, "MMFF94s", check_halted, safe_status)

def test_obabel_optimization_flow(worker):
    """Test the flow of OpenBabel optimization (mocked)."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "UFF_OBABEL",
        "worker_id": 3
    }

    # Mock both availability and the iterative function
    with patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", True), \
         patch("moleditpy.ui.calculation_worker._iterative_optimize_obabel", return_value=True) as mock_opt:
        worker.run_calculation(mol_block, options)

    assert mock_opt.called
    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result
    assert res_mol.GetProp("_pme_optimization_method") == "UFF_OBABEL"
