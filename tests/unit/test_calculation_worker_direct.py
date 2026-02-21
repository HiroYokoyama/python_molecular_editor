import pytest
from PyQt6.QtCore import QObject, pyqtSignal
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.calculation_worker import CalculationWorker
from unittest.mock import MagicMock, patch


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
    w = CalculationWorker()
    return w


def test_calculation_worker_init(worker):
    """Verify the initial state of the CalculationWorker."""
    from PyQt6.QtCore import QObject

    # Verify it is a QObject subclass
    assert isinstance(worker, QObject)

    # Verify the start_work signal exists and is connectable
    assert hasattr(worker, "start_work")

    # By default, halt_ids is not set until shared by the host.
    assert getattr(worker, "halt_ids", None) is None

    # halt_all should default to falsy to avoid premature halting
    assert not getattr(worker, "halt_all", False)


def test_calculation_worker_halt_logic(worker):
    """Test the internal _check_halted logic via run_calculation."""
    mol = Chem.MolFromSmiles("C")
    mol_block = Chem.MolToMolBlock(mol)

    # Mock error signal to catch the "Halted" message
    error_captor = SignalCaptor()
    worker.error.connect(error_captor.capture)

    # Set halt state
    worker.halt_ids = {123}
    options = {"worker_id": 123}

    worker.run_calculation(mol_block, options)

    # Should emit error with "Halted"
    assert any("Halted" in str(val) for val in error_captor.emitted_values)
    assert error_captor.emitted_values[0] == (123, "Halted")


def test_calculation_worker_direct_mode(worker):
    """Test 'direct' conversion mode which avoids RDKit 3D embedding."""
    # Simple ethanol-like structure in 2D
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {"conversion_mode": "direct", "worker_id": 1}
    worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    # Payload should be (worker_id, mol) or just mol
    result = finish_captor.emitted_values[0]
    if isinstance(result, tuple):
        res_mol = result[1]
    else:
        res_mol = result

    assert res_mol.GetNumAtoms() > 2  # Should have added Hydrogens
    # Check that Z is not all 0 (due to H-offset logic)
    conf = res_mol.GetConformer()
    z_coords = [conf.GetAtomPosition(i).z for i in range(res_mol.GetNumAtoms())]
    assert any(z > 0 for z in z_coords)


def test_calculation_worker_explicit_stereo_m_cfg(worker):
    """Test parsing of M CFG labels in MOL block."""
    # Trans-2-butene with explicit CFG 2 (E)
    mol_block = """
  2-Butene
  MoleditPy

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  1  4  1  0
M  CFG   1   1   2   2
M  END
"""
    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    # Use direct mode to verify processing
    worker.run_calculation(mol_block.strip() + "\n", {"conversion_mode": "direct"})

    assert len(finish_captor.emitted_values) > 0
    res_mol = finish_captor.emitted_values[0]
    if isinstance(res_mol, tuple):
        res_mol = res_mol[1]

    # Find double bond and check stereo
    for bond in res_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            assert bond.GetStereo() == Chem.BondStereo.STEREOE


def test_calculation_worker_error_empty_input(worker):
    error_captor = SignalCaptor()
    worker.error.connect(error_captor.capture)

    worker.run_calculation("", None)
    assert any("No atoms to convert" in str(val) for val in error_captor.emitted_values)


def test_calculation_worker_safe_helpers_halted(worker):
    """Test that safe helpers don't emit finished if halted."""
    # We need to set halt_all on the worker object
    setattr(worker, "halt_all", True)

    status_captor = SignalCaptor()
    worker.status_update.connect(status_captor.capture)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    mol = Chem.MolFromSmiles("C")
    # run_calculation uses _check_halted which checks self.halt_all if worker_id is None
    worker.run_calculation(Chem.MolToMolBlock(mol), None)

    # Finished should NOT be emitted
    assert len(finish_captor.emitted_values) == 0
    # Status should NOT be emitted (except maybe the first "Creating 3D structure..." if it happens before halt check?)


def test_calculation_worker_rdkit_embedding_fail_fallback(worker):
    """Test that embedding failure triggers fallback status or error messages."""
    status_captor = SignalCaptor()
    worker.status_update.connect(status_captor.capture)
    error_captor = SignalCaptor()
    worker.error.connect(error_captor.capture)

    mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles("C"))

    # If EmbedMolecule returns -1, it doesn't raise Exception, so it proceeds to constraint embedding
    with (
        patch("moleditpy.modules.calculation_worker.AllChem.EmbedMolecule", return_value=-1),
        patch(
            "moleditpy.modules.calculation_worker.AllChem.GetMoleculeBoundsMatrix",
            side_effect=Exception("Bounds fail"),
        ),
    ):
        worker.run_calculation(mol_block, {"conversion_mode": "rdkit"})

    all_msgs = [str(m).lower() for m in status_captor.emitted_values] + [
        str(m).lower() for m in error_captor.emitted_values
    ]
    all_msgs_str = " | ".join(all_msgs)

    # Check for specific failure indication from the embedding/conversion path
    assert any(
        keyword in all_msgs_str
        for keyword in ["embedding failed", "conversion failed", "bounds fail"]
    ), f"Expected specific embedding failure message, got: {all_msgs_str}"
