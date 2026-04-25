"""
Tests for BasePickingDialog (ui/base_picking_dialog.py).

Covers:
  - keyPressEvent: Enter/Return clicks apply_button if enabled
  - keyPressEvent: other keys delegate to super
  - closeEvent: calls clear_atom_labels and disable_picking
  - reject: calls clear_atom_labels and disable_picking
  - accept: calls clear_atom_labels and disable_picking
  - _update_molecule_geometry: updates conformer, cache, redraws
  - _update_molecule_geometry: dict-form positions
  - _push_undo: calls push_undo_state, clears _molecule_modified
  - done: pushes undo if _molecule_modified is True
  - done: skips undo if _molecule_modified is False
"""

import numpy as np
import pytest
from unittest.mock import MagicMock, patch, call
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QKeyEvent
from PyQt6.QtWidgets import QApplication
from rdkit import Chem
from rdkit.Chem import AllChem

from moleditpy.ui.base_picking_dialog import BasePickingDialog


# ---------------------------------------------------------------------------
# Minimal concrete subclass
# ---------------------------------------------------------------------------

class _ConcreteDialog(BasePickingDialog):
    def on_atom_picked(self, atom_idx: int) -> None:
        pass


def _ethane():
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_mw(mol):
    mw = MagicMock()
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    mw.view_3d_manager.plotter = MagicMock()
    return mw


def _make_dlg(app, mol=None):
    if mol is None:
        mol = _ethane()
    mw = _make_mw(mol)
    with patch.object(_ConcreteDialog, "enable_picking"):
        # Most subclasses call enable_picking in init_ui; avoid side effects
        dlg = _ConcreteDialog(mol, mw)
    return dlg, mol, mw


# ---------------------------------------------------------------------------
# keyPressEvent
# ---------------------------------------------------------------------------

def _key_event(key):
    ev = MagicMock(spec=QKeyEvent)
    ev.key.return_value = key
    return ev


def test_key_enter_clicks_apply_button_if_enabled(app):
    dlg, _, _ = _make_dlg(app)
    apply_btn = MagicMock()
    apply_btn.isEnabled.return_value = True
    dlg.apply_button = apply_btn

    ev = _key_event(Qt.Key.Key_Return)
    dlg.keyPressEvent(ev)
    apply_btn.click.assert_called_once()
    ev.accept.assert_called_once()


def test_key_enter_does_not_click_disabled_apply_button(app):
    dlg, _, _ = _make_dlg(app)
    apply_btn = MagicMock()
    apply_btn.isEnabled.return_value = False
    dlg.apply_button = apply_btn

    ev = _key_event(Qt.Key.Key_Return)
    dlg.keyPressEvent(ev)
    apply_btn.click.assert_not_called()


def test_key_enter_no_apply_button_does_not_raise(app):
    dlg, _, _ = _make_dlg(app)
    ev = _key_event(Qt.Key.Key_Return)
    dlg.keyPressEvent(ev)   # should not raise


def test_key_none_event_returns_early(app):
    dlg, _, _ = _make_dlg(app)
    dlg.keyPressEvent(None)   # should not raise


def test_key_other_key_passes_to_super(app):
    from PyQt6.QtGui import QKeyEvent
    from PyQt6.QtCore import QEvent
    dlg, _, _ = _make_dlg(app)
    ev = QKeyEvent(QEvent.Type.KeyPress, Qt.Key.Key_Escape, Qt.KeyboardModifier.NoModifier)
    # Should not raise; delegates to QDialog.keyPressEvent
    dlg.keyPressEvent(ev)


# ---------------------------------------------------------------------------
# closeEvent / reject / accept
# ---------------------------------------------------------------------------

def test_close_event_clears_labels_and_disables_picking(app):
    from PyQt6.QtGui import QCloseEvent
    dlg, _, _ = _make_dlg(app)
    with patch.object(dlg, "clear_atom_labels") as mock_clear, \
         patch.object(dlg, "disable_picking") as mock_disable:
        dlg.closeEvent(QCloseEvent())
    mock_clear.assert_called_once()
    mock_disable.assert_called_once()


def test_reject_clears_labels_and_disables_picking(app):
    dlg, _, _ = _make_dlg(app)
    with patch.object(dlg, "clear_atom_labels") as mock_clear, \
         patch.object(dlg, "disable_picking") as mock_disable, \
         patch("moleditpy.ui.base_picking_dialog.QDialog.reject"):
        dlg.reject()
    mock_clear.assert_called_once()
    mock_disable.assert_called_once()


def test_accept_clears_labels_and_disables_picking(app):
    dlg, _, _ = _make_dlg(app)
    with patch.object(dlg, "clear_atom_labels") as mock_clear, \
         patch.object(dlg, "disable_picking") as mock_disable, \
         patch("moleditpy.ui.base_picking_dialog.QDialog.accept"):
        dlg.accept()
    mock_clear.assert_called_once()
    mock_disable.assert_called_once()


# ---------------------------------------------------------------------------
# _update_molecule_geometry
# ---------------------------------------------------------------------------

def test_update_molecule_geometry_array_updates_conformer(app):
    dlg, mol, mw = _make_dlg(app)
    positions = mol.GetConformer().GetPositions().copy()
    positions[0] = [10.0, 20.0, 30.0]

    dlg._update_molecule_geometry(positions)

    updated = mol.GetConformer().GetPositions()
    assert updated[0] == pytest.approx([10.0, 20.0, 30.0], abs=1e-4)
    assert dlg._molecule_modified is True


def test_update_molecule_geometry_dict_form(app):
    dlg, mol, mw = _make_dlg(app)
    dlg._update_molecule_geometry({0: np.array([5.0, 6.0, 7.0])})

    updated = mol.GetConformer().GetPositions()
    assert updated[0] == pytest.approx([5.0, 6.0, 7.0], abs=1e-4)


def test_update_molecule_geometry_calls_draw_molecule_3d(app):
    dlg, mol, mw = _make_dlg(app)
    positions = mol.GetConformer().GetPositions().copy()
    dlg._update_molecule_geometry(positions)
    mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(mol)


def test_update_molecule_geometry_updates_cache(app):
    dlg, mol, mw = _make_dlg(app)
    positions = mol.GetConformer().GetPositions().copy()
    positions[1] = [99.0, 0.0, 0.0]
    dlg._update_molecule_geometry(positions)
    assert mw.view_3d_manager.atom_positions_3d[1] == pytest.approx([99.0, 0.0, 0.0], abs=1e-4)


# ---------------------------------------------------------------------------
# _push_undo
# ---------------------------------------------------------------------------

def test_push_undo_calls_push_undo_state(app):
    dlg, _, mw = _make_dlg(app)
    dlg._molecule_modified = True
    dlg._push_undo()
    mw.edit_actions_manager.push_undo_state.assert_called_once()
    assert dlg._molecule_modified is False


def test_push_undo_no_state_manager_no_crash(app):
    dlg, _, mw = _make_dlg(app)
    del mw.state_manager  # simulate missing attribute
    dlg._push_undo()   # should not raise


# ---------------------------------------------------------------------------
# done
# ---------------------------------------------------------------------------

def test_done_pushes_undo_if_molecule_modified(app):
    dlg, _, mw = _make_dlg(app)
    dlg._molecule_modified = True
    with patch.object(dlg, "_push_undo") as mock_undo, \
         patch("moleditpy.ui.base_picking_dialog.QDialog.done"):
        dlg.done(0)
    mock_undo.assert_called_once()


def test_done_skips_undo_if_not_modified(app):
    dlg, _, mw = _make_dlg(app)
    dlg._molecule_modified = False
    with patch.object(dlg, "_push_undo") as mock_undo, \
         patch("moleditpy.ui.base_picking_dialog.QDialog.done"):
        dlg.done(0)
    mock_undo.assert_not_called()
