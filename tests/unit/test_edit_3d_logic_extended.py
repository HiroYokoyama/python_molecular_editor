"""
Extended tests for Edit3DManager (ui/edit_3d_logic.py).

Covers uncovered paths:
  - toggle_measurement_mode when is_3d_edit_mode=True (disables edit, closes dialogs)
  - close_all_3d_edit_dialogs: closes all dialogs, clears list, handles errors
  - update_measurement_labels_display: calls plotter.add_point_labels with positions
  - update_measurement_labels_display: no labels / no mol → early return
  - clear_measurement_selection: clears labels, removes actor, renders
  - update_2d_measurement_labels: no mol → early return
  - update_2d_measurement_labels: maps atoms and adds labels
  - update_3d_selection_display: empty selection → just render
  - remove_dialog_from_list
  - calculate_and_display_measurements: 3-atom (angle) and 4-atom (dihedral) paths
"""

import numpy as np
import pytest
from unittest.mock import MagicMock, patch, call
from PyQt6.QtCore import QPointF, QRectF

from moleditpy.ui.edit_3d_logic import Edit3DManager


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_host():
    host = MagicMock()
    host.view_3d_manager.plotter = MagicMock()
    host.view_3d_manager.current_mol = None
    host.view_3d_manager.atom_positions_3d = None
    host.init_manager.settings = {"background_color": "#919191"}
    host.statusBar.return_value = MagicMock()
    return host


def _make_mgr(host=None):
    if host is None:
        host = _make_host()
    return Edit3DManager(host), host


def _positions_4():
    return np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [4.5, 0.0, 0.0],
    ])


# ---------------------------------------------------------------------------
# toggle_measurement_mode when is_3d_edit_mode=True
# ---------------------------------------------------------------------------

def test_toggle_on_when_edit_mode_active_disables_edit_mode(app):
    mgr, host = _make_mgr()
    mgr.is_3d_edit_mode = True
    host.init_manager.edit_3d_action = MagicMock()

    mgr.toggle_measurement_mode(True)

    host.init_manager.edit_3d_action.setChecked.assert_called_once_with(False)
    host.ui_manager.toggle_3d_edit_mode.assert_called_once_with(False)


def test_toggle_on_closes_active_dialogs(app):
    mgr, host = _make_mgr()
    mgr.is_3d_edit_mode = True
    host.init_manager.edit_3d_action = MagicMock()

    dlg = MagicMock()
    mgr.active_3d_dialogs.append(dlg)

    mgr.toggle_measurement_mode(True)

    dlg.close.assert_called_once()
    assert mgr.active_3d_dialogs == []


# ---------------------------------------------------------------------------
# close_all_3d_edit_dialogs
# ---------------------------------------------------------------------------

def test_close_all_closes_each_dialog(app):
    mgr, _ = _make_mgr()
    d1, d2 = MagicMock(), MagicMock()
    mgr.active_3d_dialogs = [d1, d2]

    mgr.close_all_3d_edit_dialogs()

    d1.close.assert_called_once()
    d2.close.assert_called_once()
    assert mgr.active_3d_dialogs == []


def test_close_all_handles_close_error(app):
    mgr, _ = _make_mgr()
    dlg = MagicMock()
    dlg.close.side_effect = RuntimeError("already closed")
    mgr.active_3d_dialogs = [dlg]

    mgr.close_all_3d_edit_dialogs()   # should not raise
    assert mgr.active_3d_dialogs == []


def test_close_all_empty_list_is_noop(app):
    mgr, _ = _make_mgr()
    mgr.close_all_3d_edit_dialogs()   # should not raise


# ---------------------------------------------------------------------------
# update_measurement_labels_display
# ---------------------------------------------------------------------------

def test_update_labels_display_adds_point_labels(app):
    mgr, host = _make_mgr()
    mol = MagicMock()
    mol.GetNumAtoms.return_value = 2
    host.view_3d_manager.current_mol = mol
    host.view_3d_manager.atom_positions_3d = _positions_4()

    mgr.measurement_labels = [(0, "1"), (1, "2")]
    mgr.update_measurement_labels_display()

    host.view_3d_manager.plotter.add_point_labels.assert_called_once()


def test_update_labels_display_no_labels_returns_early(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.current_mol = MagicMock()
    mgr.measurement_labels = []
    mgr.update_measurement_labels_display()
    host.view_3d_manager.plotter.add_point_labels.assert_not_called()


def test_update_labels_display_no_mol_returns_early(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.current_mol = None
    mgr.measurement_labels = [(0, "1")]
    mgr.update_measurement_labels_display()
    host.view_3d_manager.plotter.add_point_labels.assert_not_called()


# ---------------------------------------------------------------------------
# clear_measurement_selection
# ---------------------------------------------------------------------------

def test_clear_measurement_selection_clears_state(app):
    mgr, host = _make_mgr()
    mgr.selected_atoms_for_measurement = [0, 1]
    mgr.measurement_labels = [(0, "1")]
    mgr.measurement_text_actor = None

    mgr.clear_measurement_selection()

    assert mgr.selected_atoms_for_measurement == []
    assert mgr.measurement_labels == []
    host.view_3d_manager.plotter.render.assert_called()


def test_clear_measurement_selection_removes_text_actor(app):
    mgr, host = _make_mgr()
    actor = MagicMock()
    mgr.measurement_text_actor = actor

    mgr.clear_measurement_selection()

    host.view_3d_manager.plotter.remove_actor.assert_any_call(actor)
    assert mgr.measurement_text_actor is None


# ---------------------------------------------------------------------------
# update_2d_measurement_labels
# ---------------------------------------------------------------------------

def test_update_2d_labels_no_mol_returns_early(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.current_mol = None
    mgr.measurement_labels = [(0, "1")]

    with patch.object(mgr, "add_2d_measurement_label") as mock_add:
        mgr.update_2d_measurement_labels()
    mock_add.assert_not_called()


def test_update_2d_labels_no_atoms_data_returns_early(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.current_mol = MagicMock()
    host.state_manager.data.atoms = {}   # empty
    mgr.measurement_labels = [(0, "1")]

    with patch.object(mgr, "add_2d_measurement_label") as mock_add:
        mgr.update_2d_measurement_labels()
    mock_add.assert_not_called()


def test_update_2d_labels_maps_atom_and_adds_label(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.current_mol = MagicMock()
    host.state_manager.data.atoms = {1: {}}  # non-empty

    atom_item = MagicMock()
    atom_item.atom_id = 1
    atom_item.symbol = "C"

    scene = MagicMock()
    scene.items.return_value = [atom_item]
    host.init_manager.scene = scene

    mgr.measurement_labels = [(0, "1")]

    with patch.object(mgr, "find_rdkit_atom_index", return_value=0), \
         patch.object(mgr, "add_2d_measurement_label") as mock_add, \
         patch.object(mgr, "clear_2d_measurement_labels"):
        mgr.update_2d_measurement_labels()

    mock_add.assert_called_once_with(atom_item, "1")


# ---------------------------------------------------------------------------
# update_3d_selection_display
# ---------------------------------------------------------------------------

def test_update_3d_selection_empty_renders(app):
    mgr, host = _make_mgr()
    mgr.selected_atoms_3d = set()
    mgr.update_3d_selection_display()
    host.view_3d_manager.plotter.render.assert_called()


def test_update_3d_selection_no_mol_renders(app):
    mgr, host = _make_mgr()
    mgr.selected_atoms_3d = {0}
    host.view_3d_manager.current_mol = None
    mgr.update_3d_selection_display()
    host.view_3d_manager.plotter.render.assert_called()


# ---------------------------------------------------------------------------
# remove_dialog_from_list
# ---------------------------------------------------------------------------

def test_remove_dialog_present(app):
    mgr, _ = _make_mgr()
    dlg = MagicMock()
    mgr.active_3d_dialogs = [dlg]
    mgr.remove_dialog_from_list(dlg)
    assert dlg not in mgr.active_3d_dialogs


def test_remove_dialog_absent_is_noop(app):
    mgr, _ = _make_mgr()
    mgr.active_3d_dialogs = []
    mgr.remove_dialog_from_list(MagicMock())   # should not raise


# ---------------------------------------------------------------------------
# calculate_and_display_measurements — 3-atom and 4-atom paths
# ---------------------------------------------------------------------------

def test_calculate_and_display_3_atoms_includes_angle(app):
    mgr, host = _make_mgr()
    positions = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    host.view_3d_manager.atom_positions_3d = positions
    mgr.selected_atoms_for_measurement = [0, 1, 2]

    with patch.object(mgr, "display_measurement_text") as mock_disp:
        mgr.calculate_and_display_measurements()

    lines = mock_disp.call_args[0][0]
    assert any("Angle" in l for l in lines)
    assert any("Distance" in l for l in lines)


def test_calculate_and_display_4_atoms_includes_dihedral(app):
    mgr, host = _make_mgr()
    positions = np.array([
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 1.0, 1.0],
    ])
    host.view_3d_manager.atom_positions_3d = positions
    mgr.selected_atoms_for_measurement = [0, 1, 2, 3]

    with patch.object(mgr, "display_measurement_text") as mock_disp:
        mgr.calculate_and_display_measurements()

    lines = mock_disp.call_args[0][0]
    assert any("Dihedral" in l for l in lines)
    assert any("Angle" in l for l in lines)
    assert any("Distance" in l for l in lines)


def test_calculate_and_display_1_atom_does_nothing(app):
    mgr, host = _make_mgr()
    host.view_3d_manager.atom_positions_3d = _positions_4()
    mgr.selected_atoms_for_measurement = [0]

    with patch.object(mgr, "display_measurement_text") as mock_disp:
        mgr.calculate_and_display_measurements()

    mock_disp.assert_not_called()
