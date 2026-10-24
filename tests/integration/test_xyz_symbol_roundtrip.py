"""XYZ pseudo-atom labels survive both export entry points."""

from unittest.mock import MagicMock

import pytest
from PySide6.QtWidgets import QFileDialog

from moleditpy.plugins.plugin_interface import PluginContext


@pytest.mark.parametrize("symbol", ["XX", "X:", "-", "Bq"])
def test_xyz_exports_preserve_pseudo_atom_symbols(
    window, monkeypatch, tmp_path, symbol
):
    mol = window.io_manager._mol_from_xyz_lines(
        ["2", "Ghost site", f"{symbol} 1 2 3", "C 4 5 6"]
    )
    window.set_current_molecule(mol)
    path = tmp_path / "ghost.xyz"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
    window.io_manager.save_as_xyz()
    exported = path.read_text(encoding="utf-8").splitlines()
    assert exported[2].split() == [symbol, "1.00000000", "2.00000000", "3.00000000"]
    assert exported[3].split()[0] == "C"
    manager = MagicMock()
    manager.get_main_window.return_value = window
    block = PluginContext(manager, "Test plugin").to_xyz_block()
    assert block.splitlines()[0].split() == exported[2].split()
    round_trip = window.io_manager._mol_from_xyz_lines(exported)
    assert round_trip.GetAtomWithIdx(0).GetProp("xyz_original_symbol") == symbol
    assert tuple(round_trip.GetConformer().GetAtomPosition(0)) == (1, 2, 3)
