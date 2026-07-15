# -*- coding: utf-8 -*-
"""E2E import workflow tests: SMILES, MOL file, and XYZ through the real app.

Each test uses the real importer managers on a live MainWindow and asserts
the resulting document state (2D scene data or 3D viewer mode).
"""

import importlib.util

import pytest
from PyQt6.QtWidgets import QMessageBox
from rdkit import Chem
from rdkit.Chem import AllChem

_PKG = (
    "moleditpy_linux"
    if importlib.util.find_spec("moleditpy_linux") is not None
    else "moleditpy"
)

WATER_XYZ = """3
water test
O    0.000000    0.000000    0.117300
H    0.000000    0.757200   -0.469200
H    0.000000   -0.757200   -0.469200
"""


@pytest.mark.gui
def test_smiles_import_creates_scene_atoms(window, qtbot):
    """Importing ethanol SMILES adds 3 heavy atoms and 2 bonds to the 2D data."""
    window.string_importer_manager.load_from_smiles("CCO")
    qtbot.wait(150)  # drain deferred fit_to_view timer while window is alive

    assert len(window.data.atoms) == 3
    assert len(window.data.bonds) == 2
    symbols = sorted(d["symbol"] for d in window.data.atoms.values())
    assert symbols == ["C", "C", "O"]
    assert window.state_manager.has_unsaved_changes is True


@pytest.mark.gui
def test_invalid_smiles_does_not_crash_or_modify(window):
    """An invalid SMILES leaves the document untouched and raises nothing."""
    window.string_importer_manager.load_from_smiles("not_a_smiles((")

    assert len(window.data.atoms) == 0
    assert len(window.data.bonds) == 0


@pytest.mark.gui
def test_mol_file_import(window, qtbot, tmp_path):
    """Importing a benzene MOL file yields 6 kekulized ring atoms and bonds."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    AllChem.Compute2DCoords(mol)
    path = tmp_path / "benzene.mol"
    path.write_text(Chem.MolToMolBlock(mol), encoding="utf-8")

    window.io_manager.load_mol_file(str(path))
    qtbot.wait(150)  # drain deferred fit_to_view timer while window is alive

    assert len(window.data.atoms) == 6
    assert len(window.data.bonds) == 6
    orders = sorted(b["order"] for b in window.data.bonds.values())
    assert orders == [1, 1, 1, 2, 2, 2]  # kekulized benzene


@pytest.mark.gui
def test_xyz_load_enters_viewer_mode(window, qtbot, tmp_path, monkeypatch):
    """Loading an XYZ file produces a 3-atom 3D mol and switches to viewer mode."""
    path = tmp_path / "water.xyz"
    path.write_text(WATER_XYZ, encoding="utf-8")

    # Avoid both the unsaved-changes prompt and the charge dialog
    monkeypatch.setattr(
        f"{_PKG}.ui.app_state.QMessageBox.question",
        lambda *a, **k: QMessageBox.StandardButton.No,
        raising=True,
    )
    window.init_manager.settings["skip_chemistry_checks"] = True

    window.io_manager.load_xyz_for_3d_viewing(str(path))
    qtbot.wait(200)

    mol = window.view_3d_manager.current_mol
    assert mol is not None
    assert mol.GetNumAtoms() == 3
    assert sorted(a.GetSymbol() for a in mol.GetAtoms()) == ["H", "H", "O"]
    assert window.ui_manager.is_2d_editable is False


@pytest.mark.gui
def test_3d_imports_invoke_plugin_document_reset(window, qtbot, tmp_path, monkeypatch):
    """XYZ and MOL 3D-viewer imports replace the document, so the plugin
    document-reset dispatch must fire (plugins rely on it to drop state tied
    to the previous structure). The dispatch loop itself is unit-tested in
    test_plugin_manager.py; here we pin that the import paths reach it."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    monkeypatch.setattr(
        f"{_PKG}.ui.app_state.QMessageBox.question",
        lambda *a, **k: QMessageBox.StandardButton.No,
        raising=True,
    )
    window.init_manager.settings["skip_chemistry_checks"] = True

    resets = []
    monkeypatch.setattr(
        window.plugin_manager,
        "invoke_document_reset_handlers",
        lambda: resets.append(True),
        raising=True,
    )

    xyz_path = tmp_path / "water.xyz"
    xyz_path.write_text(WATER_XYZ, encoding="utf-8")
    window.io_manager.load_xyz_for_3d_viewing(str(xyz_path))
    qtbot.wait(200)
    assert len(resets) == 1, "XYZ 3D import must invoke document reset"

    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol_path = tmp_path / "methane.mol"
    mol_path.write_text(Chem.MolToMolBlock(mol), encoding="utf-8")
    window.io_manager.load_mol_file_for_3d_viewing(str(mol_path))
    qtbot.wait(200)
    assert len(resets) == 2, "MOL 3D import must invoke document reset"
