"""Exercise the installed Open Babel backend, including its bundled data."""

import importlib.util
import math
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from moleditpy.ui import calculation_worker as worker


def _run_isolated(script):
    env = os.environ.copy()
    source = str(Path(__file__).resolve().parents[2] / "moleditpy" / "src")
    env["PYTHONPATH"] = os.pathsep.join([source, env.get("PYTHONPATH", "")])
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(script)],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


def test_worker_does_not_load_openbabel_in_gui_process():
    result = _run_isolated("""
        import sys
        import moleditpy
        moleditpy.OBABEL_AVAILABLE = True
        sys.modules["openbabel"] = None
        from moleditpy.ui import calculation_worker
        assert calculation_worker.OBABEL_AVAILABLE
        assert sys.modules["openbabel"] is None
        assert "openbabel._openbabel" not in sys.modules
    """)
    assert result.returncode == 0, result.stderr


@pytest.fixture
def openbabel_backend():
    if importlib.util.find_spec("openbabel") is None:
        pytest.skip("Open Babel is not installed")
    assert worker.OBABEL_AVAILABLE
    return worker


def test_worker_configures_openbabel_runtime_paths(openbabel_backend):
    result = _run_isolated("""
        import os
        from pathlib import Path
        from rdkit import Chem
        from openbabel import pybel
        libdir = os.environ.get("BABEL_LIBDIR")
        from moleditpy.utils.openbabel_worker import _load_pybel
        pybel = _load_pybel()
        assert libdir == os.environ.get("BABEL_LIBDIR")
        assert (Path(os.environ["BABEL_DATADIR"]) / "UFF.prm").is_file()
        mol = pybel.readstring("smi", "CCO")
        assert mol.OBMol.NumAtoms() == 3
        assert mol.write("mol").strip().endswith("M  END")
    """)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("method", ["MMFF94s", "MMFF94", "UFF", "GAFF", "Ghemical"])
def test_openbabel_forcefield_parameters_load(openbabel_backend, method):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert rdDistGeom.EmbedMolecule(mol, randomSeed=42) >= 0
    messages = []
    assert openbabel_backend._iterative_optimize_obabel(
        mol, method, lambda: False, messages.append, max_iters=20, chunk_size=10
    ), messages
    positions = mol.GetConformer().GetPositions()
    assert all(math.isfinite(float(value)) for row in positions for value in row)
    distance = (
        mol.GetConformer()
        .GetAtomPosition(0)
        .Distance(mol.GetConformer().GetAtomPosition(1))
    )
    assert 1.2 < distance < 1.8


def test_openbabel_subprocess_conversion(openbabel_backend):
    mol = Chem.MolFromSmiles("CCO")
    results = []
    messages = []
    assert openbabel_backend._perform_obabel_conversion(
        Chem.MolToMolBlock(mol),
        "obabel",
        "UFF_OBABEL",
        "runtime-test",
        None,
        lambda: False,
        messages.append,
        results.append,
    ), messages
    worker_id, converted = results[0]
    assert worker_id == "runtime-test"
    assert converted.GetNumAtoms() == 9
    assert converted.GetProp("_pme_conversion_backend") == "Open Babel"
    assert converted.GetProp("_pme_optimization_method") == "UFF_OBABEL"
    assert converted.GetConformer().Is3D()
    positions = converted.GetConformer().GetPositions()
    assert all(math.isfinite(float(value)) for row in positions for value in row)


def test_openbabel_cancellation_reaps_child(openbabel_backend, monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert rdDistGeom.EmbedMolecule(mol, randomSeed=42) >= 0
    positions = mol.GetConformer().GetPositions().copy()
    processes = []
    real_popen = subprocess.Popen

    def track_process(*args, **kwargs):
        process = real_popen(*args, **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr(worker.subprocess, "Popen", track_process)
    checks = iter([False, True])
    with pytest.raises(worker.WorkerHaltError):
        openbabel_backend._iterative_optimize_obabel(
            mol, "UFF", lambda: next(checks), lambda message: None
        )
    assert len(processes) == 1
    assert processes[0].poll() is not None
    assert (mol.GetConformer().GetPositions() == positions).all()


def test_openbabel_child_failure_preserves_coordinates(openbabel_backend):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert rdDistGeom.EmbedMolecule(mol, randomSeed=42) >= 0
    positions = mol.GetConformer().GetPositions().copy()
    messages = []
    assert not openbabel_backend._iterative_optimize_obabel(
        mol, "unknown-forcefield", lambda: False, messages.append
    )
    assert messages
    assert (mol.GetConformer().GetPositions() == positions).all()
