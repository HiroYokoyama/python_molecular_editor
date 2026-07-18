"""Unit tests for CalculationWorker optimization and intermolecular modes."""

import pytest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.calculation_worker import CalculationWorker
from unittest.mock import patch


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


# ---------------------------------------------------------------------------
# Direct / core worker tests
# ---------------------------------------------------------------------------


def test_calculation_worker_init(worker):
    """CalculationWorker initialises as a QObject with halt attributes unset."""
    from PyQt6.QtCore import QObject

    assert isinstance(worker, QObject)
    assert getattr(worker, "halt_ids", None) is None
    assert not getattr(worker, "halt_all", False)


def test_calculation_worker_explicit_stereo_m_cfg(worker):
    """M CFG stereo block in V2000 molfile is preserved through direct conversion."""
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

    worker.run_calculation(mol_block.strip() + "\n", {"conversion_mode": "direct"})

    assert len(finish_captor.emitted_values) > 0
    res_mol = finish_captor.emitted_values[0]
    if isinstance(res_mol, tuple):
        res_mol = res_mol[1]

    for bond in res_mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            assert bond.GetStereo() == Chem.BondStereo.STEREOE


def test_calculation_worker_error_empty_input(worker):
    """Empty mol block emits an error signal indicating no atoms to convert."""
    error_captor = SignalCaptor()
    worker.error.connect(error_captor.capture)

    worker.run_calculation("", None)
    assert any("No atoms to convert" in str(val) for val in error_captor.emitted_values)


def test_calculation_worker_none_conversion_mode_defaults_to_fallback(worker):
    """conversion_mode=None falls through to the RDKit workflow path."""
    mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles("C"))

    with patch.object(worker, "_run_rdkit_workflow", return_value=True) as mock_rdkit:
        worker.run_calculation(mol_block, {"conversion_mode": None})

    assert mock_rdkit.called


def test_calculation_worker_safe_helpers_halted(worker):
    """halt_all=True prevents the finished signal from being emitted."""
    setattr(worker, "halt_all", True)

    status_captor = SignalCaptor()
    worker.status_update.connect(status_captor.capture)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    mol = Chem.MolFromSmiles("C")
    worker.run_calculation(Chem.MolToMolBlock(mol), None)

    assert len(finish_captor.emitted_values) == 0


def test_calculation_worker_rdkit_embedding_fail_fallback(worker):
    """Embedding failure is reported via status/error signal and does not raise."""
    status_captor = SignalCaptor()
    worker.status_update.connect(status_captor.capture)
    error_captor = SignalCaptor()
    worker.error.connect(error_captor.capture)

    mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles("C"))

    with (
        patch("moleditpy.ui.calculation_worker.AllChem.EmbedMolecule", return_value=-1),
        patch(
            "moleditpy.ui.calculation_worker.AllChem.GetMoleculeBoundsMatrix",
            side_effect=RuntimeError("Bounds fail"),
        ),
    ):
        worker.run_calculation(mol_block, {"conversion_mode": "rdkit"})

    all_msgs = [str(m).lower() for m in status_captor.emitted_values] + [
        str(m).lower() for m in error_captor.emitted_values
    ]
    all_msgs_str = " | ".join(all_msgs)

    assert any(
        keyword in all_msgs_str
        for keyword in ["embedding failed", "conversion failed", "bounds fail"]
    ), f"Expected specific embedding failure message, got: {all_msgs_str}"


# ---------------------------------------------------------------------------
# Optimization tests
# ---------------------------------------------------------------------------


def test_optimize_only_mmff94s(worker):
    """optimize_only with MMFF94s sets the _pme_optimization_method property."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "MMFF94s_RDKIT",
        "worker_id": 1,
    }

    worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result

    assert res_mol.HasProp("_pme_optimization_method")
    assert res_mol.GetProp("_pme_optimization_method").upper() == "MMFF94S_RDKIT"


def test_optimize_only_uff(worker):
    """optimize_only with UFF_RDKIT sets the _pme_optimization_method property."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "UFF_RDKIT",
        "worker_id": 2,
    }

    worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result

    assert res_mol.GetProp("_pme_optimization_method") == "UFF_RDKIT"


def test_collision_avoidance_trigger(worker):
    """Overlapping atoms in direct mode trigger collision avoidance separation."""
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (0, 0, 0))
    mol.AddConformer(conf)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {"conversion_mode": "direct", "do_optimize": True}

    with patch(
        "moleditpy.ui.calculation_worker._iterative_optimize", return_value=True
    ):
        worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) > 0
    res_mol = finish_captor.emitted_values[0]
    if isinstance(res_mol, tuple):
        res_mol = res_mol[1]

    conf = res_mol.GetConformer()
    p1 = conf.GetAtomPosition(0)
    p2 = conf.GetAtomPosition(1)
    dist = (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2
    assert dist > 0.01


def test_iterative_optimize_halt(worker):
    """_iterative_optimize raises WorkerHaltError when halt_ids contains the worker id."""
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
    """optimize_only with UFF_OBABEL calls the obabel optimization helper."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "UFF_OBABEL",
        "worker_id": 3,
    }

    with (
        patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", True),
        patch(
            "moleditpy.ui.calculation_worker._iterative_optimize_obabel",
            return_value=True,
        ) as mock_opt,
    ):
        worker.run_calculation(mol_block, options)

    assert mock_opt.called
    assert len(finish_captor.emitted_values) > 0
    result = finish_captor.emitted_values[0]
    res_mol = result[1] if isinstance(result, tuple) else result
    assert res_mol.GetProp("_pme_optimization_method") == "UFF_OBABEL"


# ---------------------------------------------------------------------------
# Intermolecular interaction tests
# ---------------------------------------------------------------------------


def test_intermolecular_interaction_toggle():
    """Enabling intermolecular interaction optimization draws separate molecules together."""
    mol1 = Chem.MolFromSmiles("C")
    mol1 = Chem.AddHs(mol1)
    AllChem.EmbedMolecule(mol1, randomSeed=42)

    mol2 = Chem.MolFromSmiles("C")
    mol2 = Chem.AddHs(mol2)
    AllChem.EmbedMolecule(mol2, randomSeed=43)

    combined = Chem.CombineMols(mol1, mol2)
    Chem.SanitizeMol(combined)

    conf = combined.GetConformer()
    for i in range(mol1.GetNumAtoms(), combined.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 6.0, pos.y, pos.z))

    from moleditpy.ui.calculation_worker import _iterative_optimize

    def check_halted():
        return False

    def safe_status(msg):
        pass

    mol_on = Chem.Mol(combined)
    _iterative_optimize(
        mol_on,
        "MMFF94s",
        check_halted,
        safe_status,
        options={"optimize_intermolecular_interaction_rdkit": True},
        max_iters=500,
    )
    dist_on = np.linalg.norm(
        np.array(mol_on.GetConformer().GetAtomPosition(0))
        - np.array(mol_on.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )

    mol_off = Chem.Mol(combined)
    _iterative_optimize(
        mol_off,
        "MMFF94s",
        check_halted,
        safe_status,
        options={"optimize_intermolecular_interaction_rdkit": False},
        max_iters=500,
    )
    dist_off = np.linalg.norm(
        np.array(mol_off.GetConformer().GetAtomPosition(0))
        - np.array(mol_off.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )

    assert abs(dist_off - 6.0) < 1e-5
    assert dist_on < 5.0


def test_intermolecular_interaction_uff():
    """UFF with intermolecular interaction disabled keeps molecules at their initial separation."""
    mol1 = Chem.MolFromSmiles("C")
    mol1 = Chem.AddHs(mol1)
    AllChem.EmbedMolecule(mol1, randomSeed=42)

    mol2 = Chem.MolFromSmiles("C")
    mol2 = Chem.AddHs(mol2)
    AllChem.EmbedMolecule(mol2, randomSeed=43)

    combined = Chem.CombineMols(mol1, mol2)
    Chem.SanitizeMol(combined)

    conf = combined.GetConformer()
    for i in range(mol1.GetNumAtoms(), combined.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 6.0, pos.y, pos.z))

    from moleditpy.ui.calculation_worker import _iterative_optimize

    def check_halted():
        return False

    def safe_status(msg):
        pass

    mol_off = Chem.Mol(combined)
    _iterative_optimize(
        mol_off,
        "UFF",
        check_halted,
        safe_status,
        options={"optimize_intermolecular_interaction_rdkit": False},
        max_iters=500,
    )
    dist_off = np.linalg.norm(
        np.array(mol_off.GetConformer().GetAtomPosition(0))
        - np.array(mol_off.GetConformer().GetAtomPosition(mol1.GetNumAtoms()))
    )

    assert abs(dist_off - 6.0) < 1e-5


# ---------------------------------------------------------------------------
# UFF fallback removal tests
# ---------------------------------------------------------------------------


def test_iterative_optimize_mmff_unsupported_returns_false():
    """_iterative_optimize returns False (no UFF fallback) when MMFF props are None."""
    from moleditpy.ui.calculation_worker import _iterative_optimize

    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    def check_halted():
        return False

    status_msgs = []

    def safe_status(msg):
        status_msgs.append(msg)

    with patch(
        "moleditpy.ui.calculation_worker.AllChem.MMFFGetMoleculeProperties",
        return_value=None,
    ):
        result = _iterative_optimize(mol, "MMFF94s", check_halted, safe_status)

    assert result is False
    assert any("Failed to setup" in m for m in status_msgs)


def test_optimize_only_mmff_unsupported_emits_error():
    """optimize_only with unsupported MMFF (props=None) emits an error signal, not success."""
    worker = CalculationWorker()

    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol_block = Chem.MolToMolBlock(mol)

    finish_captor = SignalCaptor()
    error_captor = SignalCaptor()
    worker.finished.connect(finish_captor.capture)
    worker.error.connect(error_captor.capture)

    options = {
        "conversion_mode": "optimize_only",
        "optimization_method": "MMFF_RDKIT",
        "worker_id": 99,
    }

    with patch(
        "moleditpy.ui.calculation_worker.AllChem.MMFFGetMoleculeProperties",
        return_value=None,
    ):
        worker.run_calculation(mol_block, options)

    assert len(finish_captor.emitted_values) == 0, "No success when MMFF unsupported"
    assert len(error_captor.emitted_values) > 0, "Error signal expected"
    err_msg = str(error_captor.emitted_values[0])
    assert "MMFF" in err_msg.upper() or "failed" in err_msg.lower()


# ---------------------------------------------------------------------------
# Pure module-level helpers
# ---------------------------------------------------------------------------

from moleditpy.ui.calculation_worker import (
    _resolve_method_key,
    _parse_explicit_stereo,
    _apply_explicit_stereo,
    _iterative_optimize_obabel,
)


@pytest.mark.parametrize(
    "method, expected",
    [
        ("UFF", "UFF"),
        ("uff optimize", "UFF"),
        ("GAFF", "GAFF"),
        ("Ghemical", "GHEMICAL"),
        ("MMFF94", "MMFF94"),
        ("MMFF94S", "MMFF94s"),
        ("mmff94s", "MMFF94s"),
        ("something else", "MMFF94s"),
    ],
)
def test_resolve_method_key(method, expected):
    assert _resolve_method_key(method) == expected


def test_parse_explicit_stereo_z_and_e():
    mol_block = "\n".join(
        [
            "header",
            "M  CFG  1   2   1",  # bond idx 1 -> Z
            "M  CFG  1   3   2",  # bond idx 2 -> E
            "M  END",
        ]
    )
    result = _parse_explicit_stereo(mol_block)
    assert result == {
        1: Chem.BondStereo.STEREOZ,
        2: Chem.BondStereo.STEREOE,
    }


def test_parse_explicit_stereo_ignores_short_line():
    result = _parse_explicit_stereo("M  CFG  1\nM  END")
    assert result == {}


def test_parse_explicit_stereo_skips_non_integer_fields():
    # int(parts[3]) raises ValueError -> the line is skipped
    result = _parse_explicit_stereo("M  CFG  1   x   y")
    assert result == {}


def test_apply_explicit_stereo_sets_double_bond_stereo():
    mol = Chem.MolFromSmiles("CC=CC")  # double bond at index 1
    _apply_explicit_stereo(mol, {1: Chem.BondStereo.STEREOZ})
    bond = mol.GetBondWithIdx(1)
    assert bond.GetStereo() == Chem.BondStereo.STEREOZ
    assert bond.GetStereoAtoms()[0] == 0  # methyl carbons chosen as stereo atoms
    assert bond.GetStereoAtoms()[1] == 3


def test_apply_explicit_stereo_ignores_out_of_range_bond():
    mol = Chem.MolFromSmiles("CC=CC")
    _apply_explicit_stereo(mol, {99: Chem.BondStereo.STEREOE})  # no such bond
    assert mol.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREONONE


def test_apply_explicit_stereo_skips_single_bond():
    mol = Chem.MolFromSmiles("CCC")  # only single bonds
    _apply_explicit_stereo(mol, {0: Chem.BondStereo.STEREOZ})
    assert mol.GetBondWithIdx(0).GetStereo() == Chem.BondStereo.STEREONONE


def test_iterative_optimize_obabel_unavailable_reports_and_returns_false():
    mol = Chem.MolFromSmiles("CC")
    AllChem.EmbedMolecule(Chem.AddHs(mol), randomSeed=1)
    messages = []
    with patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", False):
        result = _iterative_optimize_obabel(
            mol,
            "UFF",
            check_halted_cb=lambda: False,
            safe_status_cb=messages.append,
        )
    assert result is False
    assert any("OpenBabel" in m for m in messages)
