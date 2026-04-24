"""
Unit tests for ui/translation_dialog.py (TranslationDialog).

Covers:
  - _on_tab_changed: clears selection (not during init)
  - Absolute tab (_TAB_ABSOLUTE=0):
      - _abs_on_atom_picked: single-atom enforce, populates coordinate inputs
      - _abs_clear_selection: resets inputs and label
      - _set_origin: sets inputs to 0.0000
      - _on_move_mol_toggled: button label changes
      - apply_absolute: guard (no atom), guard (bad input), moves whole molecule, moves atom-only
      - apply_absolute pushes undo
  - Delta tab (_TAB_DELTA=1):
      - _delta_on_atom_picked: toggle, accumulate
      - clear_selection: empties set
      - select_all_atoms: selects all N atoms
      - apply_translation: guard (no atoms), guard (bad input), zero-vector no-op, correct shift
      - apply_translation pushes undo
  - update_display: correct label and button state per tab and count
  - preselected_atoms: 1 atom → absolute tab populated, many atoms → delta tab
"""

import os
import sys
import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _ethane():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_mw(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    pos = np.array(mol.GetConformer().GetPositions(), dtype=float)
    mw.view_3d_manager.atom_positions_3d = pos
    mw.view_3d_manager.current_mol = mol
    return mw


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None, preselected_atoms=None):
        from moleditpy.ui.translation_dialog import TranslationDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_mw(_mol)
        with (
            patch.object(TranslationDialog, "show_atom_labels"),
            patch.object(TranslationDialog, "clear_atom_labels"),
            patch.object(TranslationDialog, "enable_picking"),
            patch.object(TranslationDialog, "disable_picking"),
        ):
            dlg = TranslationDialog(_mol, mw, preselected_atoms=preselected_atoms)
        created.append(dlg)
        return dlg, _mol, mw

    yield _factory

    for dlg in created:
        try:
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Tab switching
# ---------------------------------------------------------------------------


class TestTabSwitching:
    def test_tab_change_clears_selection(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg._on_tab_changed(1)
        assert len(dlg.selected_atoms) == 0

    def test_tab_change_during_init_is_noop(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0}
        dlg._is_initializing = True
        dlg._on_tab_changed(0)
        assert 0 in dlg.selected_atoms  # not cleared


# ---------------------------------------------------------------------------
# Absolute tab — picking
# ---------------------------------------------------------------------------


class TestAbsolutePicking:
    def test_abs_pick_enforces_single_atom(self, make_dialog):
        dlg, mol, mw = make_dialog()
        dlg.tabs.setCurrentIndex(0)
        mw.view_3d_manager.current_mol = mol
        with patch.object(type(dlg), "show_atom_labels"):
            dlg._abs_on_atom_picked(0)
            dlg._abs_on_atom_picked(1)  # replaces, not adds
        assert dlg.selected_atoms == {1}

    def test_abs_pick_populates_inputs(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        pos = mol.GetConformer().GetPositions()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg._abs_on_atom_picked(0)
        assert dlg.abs_x_input.text() == f"{pos[0][0]:.4f}"
        assert dlg.abs_y_input.text() == f"{pos[0][1]:.4f}"
        assert dlg.abs_z_input.text() == f"{pos[0][2]:.4f}"


# ---------------------------------------------------------------------------
# Absolute tab — helpers
# ---------------------------------------------------------------------------


class TestAbsoluteHelpers:
    def test_abs_clear_resets_inputs(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.abs_x_input.setText("1.5")
        dlg.abs_y_input.setText("-2.0")
        dlg.abs_z_input.setText("3.0")
        dlg.selected_atoms = {0}
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg._abs_clear_selection()
        assert dlg.abs_x_input.text() == "0.000"
        assert dlg.abs_y_input.text() == "0.000"
        assert dlg.abs_z_input.text() == "0.000"
        assert len(dlg.selected_atoms) == 0

    def test_set_origin_fills_zeros(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.abs_x_input.setText("5.0")
        dlg._set_origin()
        assert dlg.abs_x_input.text() == "0.0000"
        assert dlg.abs_y_input.text() == "0.0000"
        assert dlg.abs_z_input.text() == "0.0000"

    def test_move_mol_toggled_changes_button_label(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.move_mol_checkbox.setChecked(False)
        dlg._on_move_mol_toggled(0)
        assert dlg.abs_apply_btn.text() == "Move Atom"
        dlg.move_mol_checkbox.setChecked(True)
        dlg._on_move_mol_toggled(2)
        assert dlg.abs_apply_btn.text() == "Move Molecule"


# ---------------------------------------------------------------------------
# apply_absolute
# ---------------------------------------------------------------------------


class TestApplyAbsolute:
    def test_no_atom_shows_warning(self, make_dialog):
        dlg, _, mw = make_dialog()
        dlg.selected_atoms.clear()
        mw.view_3d_manager.current_mol = dlg.mol
        with patch("moleditpy.ui.translation_dialog.QMessageBox") as mb:
            dlg.apply_absolute()
        mb.warning.assert_called_once()

    def test_bad_input_shows_warning(self, make_dialog):
        dlg, mol, mw = make_dialog()
        dlg.selected_atoms = {0}
        mw.view_3d_manager.current_mol = mol
        dlg.abs_x_input.setText("bad")
        with patch("moleditpy.ui.translation_dialog.QMessageBox") as mb:
            dlg.apply_absolute()
        mb.warning.assert_called_once()

    def test_move_molecule_shifts_all_atoms(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms = {0}
        before = mol.GetConformer().GetPositions().copy()
        # Move atom 0 to absolute (0,0,0) with move_mol=True
        dlg.abs_x_input.setText("0.0000")
        dlg.abs_y_input.setText("0.0000")
        dlg.abs_z_input.setText("0.0000")
        dlg.move_mol_checkbox.setChecked(True)
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_absolute()
        after = mol.GetConformer().GetPositions()
        delta = np.array([0.0, 0.0, 0.0]) - before[0]
        for i in range(mol.GetNumAtoms()):
            assert after[i] == pytest.approx(before[i] + delta, abs=1e-4)

    def test_move_atom_only_shifts_one_atom(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms = {0}
        before = mol.GetConformer().GetPositions().copy()
        dlg.abs_x_input.setText("0.0000")
        dlg.abs_y_input.setText("0.0000")
        dlg.abs_z_input.setText("0.0000")
        dlg.move_mol_checkbox.setChecked(False)
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_absolute()
        after = mol.GetConformer().GetPositions()
        assert after[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-4)
        # Other atoms unchanged
        for i in range(1, mol.GetNumAtoms()):
            assert after[i] == pytest.approx(before[i], abs=1e-4)

    def test_apply_absolute_pushes_undo(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms = {0}
        dlg.abs_x_input.setText("0.0000")
        dlg.abs_y_input.setText("0.0000")
        dlg.abs_z_input.setText("0.0000")
        dlg.move_mol_checkbox.setChecked(True)
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_absolute()
        mw.edit_actions_manager.push_undo_state.assert_called()


# ---------------------------------------------------------------------------
# Delta tab — picking
# ---------------------------------------------------------------------------


class TestDeltaPicking:
    def test_delta_pick_adds_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(1)
        with patch.object(type(dlg), "show_atom_labels"):
            dlg._delta_on_atom_picked(0)
        assert 0 in dlg.selected_atoms

    def test_delta_repick_removes_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(1)
        with patch.object(type(dlg), "show_atom_labels"):
            dlg._delta_on_atom_picked(0)
            dlg._delta_on_atom_picked(0)
        assert 0 not in dlg.selected_atoms

    def test_delta_picks_accumulate(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(1)
        with patch.object(type(dlg), "show_atom_labels"):
            for i in range(3):
                dlg._delta_on_atom_picked(i)
        assert dlg.selected_atoms == {0, 1, 2}


# ---------------------------------------------------------------------------
# select_all_atoms (delta tab)
# ---------------------------------------------------------------------------


class TestSelectAllAtoms:
    def test_select_all_selects_every_atom(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.select_all_atoms()
        assert dlg.selected_atoms == set(range(mol.GetNumAtoms()))


# ---------------------------------------------------------------------------
# apply_translation (delta tab)
# ---------------------------------------------------------------------------


class TestApplyTranslation:
    def test_no_atoms_shows_warning(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.translation_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_bad_input_shows_warning(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms = {0}
        dlg.dx_input.setText("bad")
        with patch("moleditpy.ui.translation_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_zero_vector_is_noop(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        before = mol.GetConformer().GetPositions().copy()
        dlg.selected_atoms = {0}
        dlg.dx_input.setText("0.0")
        dlg.dy_input.setText("0.0")
        dlg.dz_input.setText("0.0")
        dlg.apply_translation()
        after = mol.GetConformer().GetPositions()
        assert after == pytest.approx(before, abs=1e-6)

    def test_delta_shifts_selected_atoms(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        before = mol.GetConformer().GetPositions().copy()
        dlg.selected_atoms = {0, 1}
        dlg.dx_input.setText("1.0")
        dlg.dy_input.setText("2.0")
        dlg.dz_input.setText("3.0")
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_translation()
        after = mol.GetConformer().GetPositions()
        for idx in [0, 1]:
            assert after[idx] == pytest.approx(before[idx] + [1.0, 2.0, 3.0], abs=1e-4)
        # Unselected atoms unchanged
        for idx in range(2, mol.GetNumAtoms()):
            assert after[idx] == pytest.approx(before[idx], abs=1e-4)

    def test_apply_translation_pushes_undo(self, make_dialog):
        dlg, mol, mw = make_dialog()
        mw.view_3d_manager.current_mol = mol
        dlg.selected_atoms = {0}
        dlg.dx_input.setText("1.0")
        dlg.dy_input.setText("0.0")
        dlg.dz_input.setText("0.0")
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_translation()
        mw.edit_actions_manager.push_undo_state.assert_called()


# ---------------------------------------------------------------------------
# update_display
# ---------------------------------------------------------------------------


class TestUpdateDisplay:
    def test_abs_tab_no_atom_disables_button(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(0)
        dlg.selected_atoms.clear()
        dlg.update_display()
        assert not dlg.abs_apply_btn.isEnabled()

    def test_abs_tab_one_atom_enables_button(self, make_dialog):
        dlg, mol, _ = make_dialog()
        dlg.tabs.setCurrentIndex(0)
        dlg.selected_atoms = {0}
        dlg.update_display()
        assert dlg.abs_apply_btn.isEnabled()

    def test_delta_tab_no_atom_disables_button(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(1)
        dlg.selected_atoms.clear()
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_delta_tab_one_atom_enables_button(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.tabs.setCurrentIndex(1)
        dlg.selected_atoms = {0}
        dlg.update_display()
        assert dlg.apply_button.isEnabled()


# ---------------------------------------------------------------------------
# preselected_atoms
# ---------------------------------------------------------------------------


class TestPreselectedAtoms:
    def test_single_preselected_goes_to_abs_tab(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0])
        assert dlg.tabs.currentIndex() == 0
        assert dlg.selected_atoms == {0}

    def test_multiple_preselected_goes_to_delta_tab(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0, 1, 2])
        assert dlg.tabs.currentIndex() == 1
        assert {0, 1, 2}.issubset(dlg.selected_atoms)
