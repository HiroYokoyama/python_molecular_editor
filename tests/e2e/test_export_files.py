# -*- coding: utf-8 -*-
"""E2E file export tests: 2D MOL, 3D MOL, and XYZ through the real IOManager.

Each test drives the real save function with the file dialog patched to a
tmp_path target (with NO extension, to also cover extension appending), then
re-reads the produced file with RDKit / plain parsing to prove it is valid.
"""

import importlib.util

import pytest
from PyQt6.QtCore import QPointF, Qt
from rdkit import Chem

_PKG = (
    "moleditpy_linux"
    if importlib.util.find_spec("moleditpy_linux") is not None
    else "moleditpy"
)


def _draw_ethane(window):
    scene = window.init_manager.scene
    c1 = scene.create_atom("C", QPointF(0, 0))
    c2 = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[c1], scene.atom_items[c2])
    return c1, c2


def _convert_to_3d(window, qtbot):
    qtbot.mouseClick(window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None, timeout=5000
    )


def _patch_save_dialog(monkeypatch, target):
    monkeypatch.setattr(
        f"{_PKG}.ui.io_logic.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(target), None),
        raising=True,
    )


@pytest.mark.gui
def test_save_as_mol_2d(window, tmp_path, monkeypatch):
    """2D MOL export appends .mol, is RDKit-parseable, and carries the header."""
    _draw_ethane(window)
    _patch_save_dialog(monkeypatch, tmp_path / "out")

    window.io_manager.save_as_mol()

    saved = tmp_path / "out.mol"
    assert saved.exists()
    text = saved.read_text(encoding="utf-8")
    assert "MoleditPy" in text
    mol = Chem.MolFromMolFile(str(saved))
    assert mol is not None
    assert sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C") == 2


@pytest.mark.gui
def test_save_3d_as_mol_appends_extension(window, qtbot, tmp_path, monkeypatch):
    """3D MOL export appends .mol and produces a parseable 3D structure."""
    _draw_ethane(window)
    _convert_to_3d(window, qtbot)
    _patch_save_dialog(monkeypatch, tmp_path / "bare_name")

    window.io_manager.save_3d_as_mol()

    saved = tmp_path / "bare_name.mol"
    assert saved.exists()
    mol = Chem.MolFromMolFile(str(saved), removeHs=False)
    assert mol is not None
    assert mol.GetNumConformers() == 1
    conf = mol.GetConformer()
    zs = [abs(conf.GetAtomPosition(i).z) for i in range(mol.GetNumAtoms())]
    assert any(z > 1e-4 for z in zs), "expected non-planar 3D coordinates"


@pytest.mark.gui
def test_save_as_xyz(window, qtbot, tmp_path, monkeypatch):
    """XYZ export appends .xyz with count, chrg/mult comment, and 4-field rows."""
    _draw_ethane(window)
    _convert_to_3d(window, qtbot)
    _patch_save_dialog(monkeypatch, tmp_path / "structure")

    window.io_manager.save_as_xyz()

    saved = tmp_path / "structure.xyz"
    assert saved.exists()
    lines = saved.read_text(encoding="utf-8").strip().splitlines()
    num_atoms = int(lines[0])
    assert num_atoms == window.view_3d_manager.current_mol.GetNumAtoms()
    assert "chrg =" in lines[1] and "mult =" in lines[1]
    assert len(lines) == 2 + num_atoms
    for row in lines[2:]:
        parts = row.split()
        assert len(parts) == 4
        float(parts[1]), float(parts[2]), float(parts[3])  # coordinates parse
