"""Real Qt/RDKit regressions for stale dialogs and private optimization results."""

from unittest.mock import MagicMock

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.plugins.plugin_interface import PluginContext

from moleditpy.ui.bond_length_dialog import BondLengthDialog
from moleditpy.ui.alignment_dialog import AlignmentDialog
from moleditpy.ui.constrained_optimization_dialog import (
    ConstrainedOptimizationDialog,
    ConstrainedOptimizationThread,
)


def test_close_defers_until_conversion_thread_finishes(window, monkeypatch):
    from PySide6.QtCore import QEvent

    thread = MagicMock()
    thread.isRunning.return_value = True
    compute = window.compute_manager
    compute._active_calc_threads.append(thread)
    compute.active_worker_ids.add(41)
    window.set_has_unsaved_changes(False)
    try:
        assert not window.ui_manager.handle_close_event(QEvent(QEvent.Type.Close))
        assert 41 in compute.halt_ids
        assert not compute.active_worker_ids
        thread.wait.assert_not_called()
        thread.quit.assert_not_called()
    finally:
        compute._active_calc_threads.clear()


def test_calculation_shutdown_resumes_close_after_native_thread_cleanup(
    window, monkeypatch, qtbot
):
    import threading
    from PySide6.QtCore import QCoreApplication, QEvent
    from moleditpy.ui.calculation_worker import CalculationWorker

    entered = threading.Event()
    release = threading.Event()
    completed = threading.Event()

    def controlled_calculation(worker, _block, options):
        entered.set()
        try:
            if not release.wait(10):
                worker.error.emit((options["worker_id"], "Test gate timed out"))
                return
            worker.finished.emit((options["worker_id"], None))
        finally:
            completed.set()

    monkeypatch.setattr(CalculationWorker, "run_calculation", controlled_calculation)
    close = MagicMock()
    monkeypatch.setattr(window, "close", close)
    compute = window.compute_manager
    compute.active_worker_ids.add(41)
    compute._start_calculation_worker("unused", {"worker_id": 41}, 41)
    try:
        qtbot.waitUntil(entered.is_set, timeout=5000)
        assert not window.ui_manager.handle_close_event(QEvent(QEvent.Type.Close))
        assert 41 in compute.halt_ids
        assert not compute.active_worker_ids
        close.assert_not_called()
    finally:
        release.set()
        qtbot.waitUntil(lambda: not compute._active_calc_threads, timeout=5000)
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    qtbot.waitUntil(lambda: close.call_count == 1, timeout=5000)
    assert completed.is_set()
    assert not compute._calc_objects
    assert not compute._close_after_calculation


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


def test_alignment_dialog_cannot_edit_replaced_molecule(window, monkeypatch):
    mol = install_molecule(window, monkeypatch)
    dialog = AlignmentDialog(mol, window, "x", [0, 1], parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    original = mol.GetConformer().GetPositions().copy()
    replacement = Chem.Mol(mol)

    window.set_current_molecule(replacement)

    assert dialog._invalidated
    assert not dialog.isEnabled()
    assert not dialog.isVisible()
    assert not window.edit_3d_manager.active_3d_dialogs
    history = list(window.edit_actions_manager.undo_stack)
    window.view_3d_manager.draw_molecule_3d.reset_mock()
    dialog.apply_alignment()
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), original)
    np.testing.assert_array_equal(replacement.GetConformer().GetPositions(), original)
    assert window.edit_actions_manager.undo_stack == history
    window.view_3d_manager.draw_molecule_3d.assert_not_called()


@pytest.mark.parametrize(
    "replacement", ["undo", "redo", "clear", "file", "calculation", "setter", "plugin"]
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
    elif replacement == "plugin":
        manager = MagicMock()
        manager.get_main_window.return_value = window
        PluginContext(manager, "Test plugin").current_mol = Chem.Mol(mol)
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


@pytest.mark.parametrize(
    "outcome",
    ["close", "replace", "unregistered_replace", "geometry_changed", "commit"],
)
def test_constrained_optimization_only_commits_valid_private_result(
    window, monkeypatch, outcome
):
    """Real UFF never modifies the live conformer until a valid completion."""
    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    if outcome != "unregistered_replace":
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
    elif outcome in ("replace", "unregistered_replace"):
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
    worker.finished.emit()
    dialog.reject()


def test_closing_running_optimization_cancels_worker(window, monkeypatch, qtbot):
    """Closing a real running QThread cannot destroy a worker still minimizing."""
    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    dialog.ff_combo.setCurrentText("UFF")
    dialog.apply_optimization()
    dialog.reject()
    qtbot.waitUntil(lambda: dialog._opt_thread is None, timeout=5000)
    assert not window.edit_3d_manager._optimization_threads
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), original)
    assert not window.state_manager.has_unsaved_changes


def test_running_optimization_commits_on_gui_thread(window, monkeypatch, qtbot):
    """An actual QThread result is applied and undoable through the Qt event loop."""
    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    dialog.ff_combo.setCurrentText("UFF")
    dialog.apply_optimization()
    qtbot.waitUntil(
        lambda: len(window.edit_actions_manager.undo_stack) == 2, timeout=5000
    )
    assert not np.allclose(mol.GetConformer().GetPositions(), original)
    assert window.state_manager.has_unsaved_changes
    dialog.reject()
    qtbot.waitUntil(lambda: dialog._opt_thread is None, timeout=5000)
    window.edit_actions_manager.undo()
    np.testing.assert_allclose(
        window.current_mol.GetConformer().GetPositions(), original
    )


def test_mirror_dialog_cannot_edit_replaced_molecule(window, monkeypatch):
    """Replacement disables a real mirror dialog and preserves both conformers."""
    mol = install_molecule(window, monkeypatch)
    window.dialog_manager.open_mirror_dialog()
    dialog = window.edit_3d_manager.active_3d_dialogs[0]
    replacement = Chem.Mol(mol)
    old_positions = mol.GetConformer().GetPositions().copy()
    window.set_current_molecule(replacement)
    assert dialog._invalidated
    assert not dialog.isEnabled()
    assert not dialog.isVisible()
    history = list(window.edit_actions_manager.undo_stack)
    window.view_3d_manager.draw_molecule_3d.reset_mock()
    dialog.apply_mirror()
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), old_positions)
    np.testing.assert_array_equal(
        replacement.GetConformer().GetPositions(), old_positions
    )
    assert window.edit_actions_manager.undo_stack == history
    window.view_3d_manager.draw_molecule_3d.assert_not_called()


def test_completed_optimization_can_restart_then_close(window, monkeypatch, qtbot):
    """A new run retires the finished QThread, and a closed dialog cannot restart."""
    mol = install_molecule(window, monkeypatch)
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    dialog.ff_combo.setCurrentText("UFF")
    dialog.apply_optimization()
    first_worker = dialog._opt_thread
    assert first_worker.wait(5000)
    held = []
    monkeypatch.setattr(
        ConstrainedOptimizationThread, "start", lambda worker: held.append(worker)
    )
    dialog.apply_optimization()
    second_worker = dialog._opt_thread
    assert second_worker is not first_worker
    qtbot.waitUntil(
        lambda: first_worker not in window.edit_3d_manager._optimization_threads,
        timeout=5000,
    )
    assert len(window.edit_actions_manager.undo_stack) == 1
    assert dialog._opt_thread is second_worker
    dialog.reject()
    positions = mol.GetConformer().GetPositions().copy()
    history = list(window.edit_actions_manager.undo_stack)
    dialog.apply_optimization()
    assert dialog._opt_thread is second_worker
    np.testing.assert_array_equal(mol.GetConformer().GetPositions(), positions)
    assert window.edit_actions_manager.undo_stack == history
    assert held == [second_worker]
    second_worker.run()
    second_worker.finished.emit()
    qtbot.waitUntil(lambda: dialog._opt_thread is None, timeout=5000)


@pytest.mark.parametrize(
    "action", ["close", "delete_on_close", "undo", "application_close"]
)
def test_slow_optimization_close_keeps_gui_responsive_and_worker_alive(
    window, monkeypatch, qtbot, action
):
    """A blocked minimization cannot freeze close/undo or die with its dialog."""
    from threading import Event
    from shiboken6 import isValid
    from PySide6.QtCore import Qt, QTimer
    from PySide6.QtGui import QCloseEvent
    from rdkit.Chem import rdForceFieldHelpers

    mol = install_molecule(window, monkeypatch)
    original = mol.GetConformer().GetPositions().copy()
    if action == "undo":
        mol.GetConformer().SetAtomPosition(0, (10, 20, 30))
        window.edit_actions_manager.push_undo_state()
    entered = Event()
    release = Event()
    pulse = Event()
    ff = MagicMock()

    def slow_minimization(maxIts):
        entered.set()
        if not release.wait(timeout=5):
            raise RuntimeError("test did not release minimization")
        return 0

    ff.Minimize.side_effect = slow_minimization
    monkeypatch.setattr(
        rdForceFieldHelpers, "UFFGetMoleculeForceField", lambda *a, **k: ff
    )
    dialog = ConstrainedOptimizationDialog(mol, window, parent=window)
    window.dialog_manager._open_3d_edit_dialog(dialog)
    dialog.ff_combo.setCurrentText("UFF")
    dialog.apply_optimization()
    worker = dialog._opt_thread
    retry_close = MagicMock()
    try:
        qtbot.waitUntil(entered.is_set, timeout=5000)
        assert worker.parent() is window
        if action == "delete_on_close":
            dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
            dialog.close()
            qtbot.waitUntil(lambda: not isValid(dialog), timeout=1000)
        elif action == "undo":
            window.edit_actions_manager.undo()
        elif action == "application_close":
            monkeypatch.setattr(window, "close", retry_close)
            assert not window.ui_manager.handle_close_event(QCloseEvent())
            retry_close.assert_not_called()
        else:
            dialog.close()
        assert worker.isRunning()
        assert worker in window.edit_3d_manager._optimization_threads
        QTimer.singleShot(0, pulse.set)
        qtbot.waitUntil(pulse.is_set, timeout=500)
    finally:
        release.set()
        qtbot.waitUntil(
            lambda: not window.edit_3d_manager._optimization_threads, timeout=5000
        )
    if action == "undo":
        np.testing.assert_allclose(
            window.current_mol.GetConformer().GetPositions(), original
        )
    else:
        np.testing.assert_array_equal(mol.GetConformer().GetPositions(), original)
        assert not window.state_manager.has_unsaved_changes
    if action == "application_close":
        qtbot.waitUntil(lambda: retry_close.call_count == 1, timeout=1000)
