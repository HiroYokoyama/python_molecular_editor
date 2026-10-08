"""Real Qt/RDKit regressions for stale dialogs and private optimization results."""

from unittest.mock import MagicMock

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from moleditpy.ui.bond_length_dialog import BondLengthDialog
from moleditpy.ui.constrained_optimization_dialog import (
    ConstrainedOptimizationDialog,
    ConstrainedOptimizationThread,
)


def install_molecule(window, monkeypatch):
    """Install real embedded butane while isolating GPU rendering."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    monkeypatch.setattr(window.view_3d_manager, "draw_molecule_3d", MagicMock())
    window.set_current_molecule(mol)
    window.set_3d_atom_positions(mol.GetConformer().GetPositions().copy())
    window.edit_actions_manager.reset_history()
    window.set_has_unsaved_changes(False)
    return mol


@pytest.mark.parametrize(
    "replacement", ["undo", "redo", "clear", "file", "calculation", "setter"]
)
def test_replacing_molecule_invalidates_geometry_dialog(
    window, monkeypatch, replacement
):
    """Every replacement closes the old tool and prevents later edits/history pushes."""
    mol = install_molecule(window, monkeypatch)
    mol.GetConformer().SetAtomPosition(0, (3, 4, 5))
    window.edit_actions_manager.push_undo_state()
    if replacement == "redo":
        window.edit_actions_manager.undo()
        mol = window.current_mol
    dialog = BondLengthDialog(mol, window, [0, 1], parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    dialog._molecule_modified = True
    old_positions = mol.GetConformer().GetPositions().copy()
    if replacement == "undo":
        window.edit_actions_manager.undo()
    elif replacement == "redo":
        window.edit_actions_manager.redo()
    elif replacement == "clear":
        window.edit_actions_manager.clear_all(skip_check=True)
    elif replacement == "file":
        monkeypatch.setattr(window.state_manager, "check_unsaved_changes", lambda: True)
        monkeypatch.setattr(
            window.io_manager, "load_xyz_block", lambda _: Chem.Mol(mol)
        )
        window.io_manager.show_xyz_data("data")
    elif replacement == "calculation":
        window.compute_manager.on_calculation_finished(Chem.Mol(mol))
    else:
        window.current_mol = Chem.Mol(mol)
    assert dialog._invalidated
    assert not dialog.isVisible()
    assert not window.edit_3d_manager.active_3d_dialogs
    history = list(window.edit_actions_manager.undo_stack)
    current = window.current_mol
    positions = current.GetConformer().GetPositions().copy() if current else None
    dialog._update_molecule_geometry(old_positions + 10)
    dialog._push_undo()
    dialog.on_slider_released()
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), old_positions)
    if current:
        np.testing.assert_array_equal(current.GetConformer().GetPositions(), positions)
    assert window.edit_actions_manager.undo_stack == history


@pytest.mark.parametrize("outcome", ["close", "replace", "geometry_changed", "commit"])
def test_constrained_optimization_only_commits_valid_private_result(
    window, monkeypatch, outcome
):
    """Real UFF never modifies the live conformer until a valid completion."""
    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    dialog.ff_combo.setCurrentText("UFF")
    threads = []
    monkeypatch.setattr(
        ConstrainedOptimizationThread, "start", lambda worker: threads.append(worker)
    )
    dialog.apply_optimization()
    worker = threads[0]
    assert worker.mol is not mol
    redraw = window.view_3d_manager.draw_molecule_3d
    redraw.reset_mock()
    if outcome == "close":
        dialog.reject()
    elif outcome == "replace":
        window.set_current_molecule(Chem.Mol(mol))
    elif outcome == "geometry_changed":
        mol.GetConformer().SetAtomPosition(0, (10, 20, 30))
    before_run = mol.GetConformer().GetPositions().copy()
    history_size = len(window.edit_actions_manager.undo_stack)
    worker.run()
    optimized = worker.mol.GetConformer().GetPositions()
    assert not np.allclose(optimized, original)
    if outcome == "commit":
        np.testing.assert_allclose(mol.GetConformer().GetPositions(), optimized)
        np.testing.assert_allclose(window.view_3d_manager.atom_positions_3d, optimized)
        assert window.state_manager.has_unsaved_changes
        assert len(window.edit_actions_manager.undo_stack) == history_size + 1
        redraw.assert_called_once_with(mol)
        window.edit_actions_manager.undo()
        np.testing.assert_allclose(
            window.current_mol.GetConformer().GetPositions(), original
        )
    else:
        np.testing.assert_array_equal(mol.GetConformer().GetPositions(), before_run)
        assert len(window.edit_actions_manager.undo_stack) == history_size
        assert not window.state_manager.has_unsaved_changes
        redraw.assert_not_called()
    dialog.reject()


def test_closing_running_optimization_joins_worker(window, monkeypatch):
    """Closing a real running QThread cannot destroy a worker still minimizing."""
    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    dialog.ff_combo.setCurrentText("UFF")
    dialog.apply_optimization()
    worker = dialog._opt_thread
    dialog.reject()
    assert not worker.isRunning()
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), original)
    assert not window.state_manager.has_unsaved_changes
