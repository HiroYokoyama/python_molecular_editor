# -*- coding: utf-8 -*-
"""Integration tests: ethane 2D structure + real 2D→3D conversion, all platforms.

On Linux the moleditpy-linux package is used (mirrors CI / installer behaviour).
On Windows/macOS the standard moleditpy package is used.
All tests are headless — no display or Qt window is required.
"""

import os
import sys
import platform
import importlib
import importlib.util
import pytest
from PyQt6.QtCore import QPointF

# ---------------------------------------------------------------------------
# Package selection: moleditpy_linux on Linux, moleditpy everywhere else.
# Set MOLEDITPY_USE_INSTALLED=1 to use whichever package is pip-installed
# instead of the local source tree (used by the pip-install CI workflow).
# ---------------------------------------------------------------------------
_IS_LINUX = sys.platform.startswith("linux")

_LINUX_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy-linux", "src")
)
_MAIN_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)

if os.environ.get("MOLEDITPY_USE_INSTALLED") == "1":
    # CI pip-install mode: use whichever package is installed; skip local source.
    if importlib.util.find_spec("moleditpy_linux") is not None:
        _PKG = "moleditpy_linux"
    else:
        _PKG = "moleditpy"
elif _IS_LINUX and os.path.isdir(_LINUX_SRC):
    if _LINUX_SRC not in sys.path:
        sys.path.insert(0, _LINUX_SRC)
    _PKG = "moleditpy_linux"
else:
    if _MAIN_SRC not in sys.path:
        sys.path.insert(0, _MAIN_SRC)
    _PKG = "moleditpy"

MolecularData = importlib.import_module(f"{_PKG}.core.molecular_data").MolecularData
CalculationWorker = importlib.import_module(
    f"{_PKG}.ui.calculation_worker"
).CalculationWorker
_pkg_mod = importlib.import_module(_PKG)
OBABEL_AVAILABLE = getattr(_pkg_mod, "OBABEL_AVAILABLE", False)

_OS_LABEL = platform.system()  # "Windows" / "Linux" / "Darwin"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_ethane() -> str:
    """Return a MOL-block string for ethane (2 C, 1 bond) via MolecularData."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    return data.to_mol_block()


def _run_worker_sync(mol_block: str, mode: str) -> object:
    """Run CalculationWorker synchronously and return the RDKit mol (or None)."""
    result_holder: list = []
    worker = CalculationWorker()
    worker.finished.connect(lambda result: result_holder.append(result))
    worker.run_calculation(mol_block, {"conversion_mode": mode, "do_optimize": True})
    return result_holder[0][1] if result_holder and result_holder[0] else None


def _assert_ethane_3d(mol, mode_label: str) -> None:
    """Common geometry checks for a 3D ethane mol."""
    import numpy as np

    assert mol is not None, f"[{_OS_LABEL}/{mode_label}] conversion returned None"
    assert mol.GetNumConformers() > 0, f"[{_OS_LABEL}/{mode_label}] no 3D conformer"

    carbons = [a for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    assert len(carbons) == 2, (
        f"[{_OS_LABEL}/{mode_label}] expected 2 carbons, got {len(carbons)}"
    )

    conf = mol.GetConformer()
    p0 = np.array(conf.GetAtomPosition(carbons[0].GetIdx()))
    p1 = np.array(conf.GetAtomPosition(carbons[1].GetIdx()))
    cc_dist = float(np.linalg.norm(p1 - p0))

    assert 1.3 < cc_dist < 1.7, (
        f"[{_OS_LABEL}/{mode_label}] C-C distance {cc_dist:.3f} Å outside [1.3, 1.7]"
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_ethane_2d_structure():
    """Ethane can be built as a 2D MolecularData on every platform."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)

    assert len(data.atoms) == 2
    assert len(data.bonds) == 1
    symbols = [v["symbol"] for v in data.atoms.values()]
    assert symbols.count("C") == 2

    mol_block = data.to_mol_block()
    assert "C" in mol_block
    assert mol_block.strip() != ""


def test_ethane_molblock_roundtrip():
    """MOL-block from MolecularData is parseable by RDKit on every platform."""
    from rdkit import Chem

    mol = Chem.MolFromMolBlock(_build_ethane(), removeHs=False)
    assert mol is not None
    carbons = [a for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    assert len(carbons) == 2


def test_ethane_conversion_rdkit(app, qtbot):
    """RDKit mode: ethane 2D→3D produces a valid 3D mol on this platform."""
    mol = _run_worker_sync(_build_ethane(), "rdkit")
    _assert_ethane_3d(mol, "rdkit")


def test_ethane_conversion_direct(app, qtbot):
    """Direct mode: ethane 2D→3D produces a valid 3D mol on this platform."""
    mol = _run_worker_sync(_build_ethane(), "direct")
    _assert_ethane_3d(mol, "direct")


def test_ethane_conversion_fallback(app, qtbot):
    """Fallback mode: RDKit → obabel chain must produce a valid 3D mol."""
    mol = _run_worker_sync(_build_ethane(), "fallback")
    _assert_ethane_3d(mol, "fallback")


@pytest.mark.skipif(
    not OBABEL_AVAILABLE,
    reason="Open Babel not available on this platform",
)
def test_ethane_conversion_obabel(app, qtbot):
    """Obabel mode: only runs when Open Babel is present (typically Linux)."""
    mol = _run_worker_sync(_build_ethane(), "obabel")
    _assert_ethane_3d(mol, "obabel")


@pytest.mark.skipif(not _IS_LINUX, reason="Linux-only: verify Linux package is used")
def test_linux_package_is_active():
    """On Linux, moleditpy_linux must be the active package."""
    assert _PKG == "moleditpy_linux"


@pytest.mark.skipif(
    _IS_LINUX or _PKG != "moleditpy",
    reason="Windows/macOS-only: verify main package is used (skipped when moleditpy_linux is active)",
)
def test_main_package_is_active():
    """On Windows/macOS, moleditpy must be the active package."""
    assert _PKG == "moleditpy"
