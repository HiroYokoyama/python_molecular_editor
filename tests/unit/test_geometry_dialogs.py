"""
Unit tests for BondLengthDialog, AngleDialog, DihedralDialog.

All three share a common factory pattern.  Tests focus on:
  - Atom-picking sequences and state transitions
  - clear_selection resets all state
  - _is_selection_complete gate
  - update_display label text at each stage
  - apply_changes: wrapping, invalid-input guard, calls apply_geometry_update
  - apply_geometry_update: correct coordinate math (atom-only and group modes)
  - Preselected atoms loaded on init
  - on_*_input_changed syncs slider
"""

import os
import sys
import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from moleditpy.ui.bond_length_dialog import BondLengthDialog
from moleditpy.ui.angle_dialog import AngleDialog
from moleditpy.ui.dihedral_dialog import DihedralDialog


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _ethane():
    """Ethane with 3D coords: C0-C1, H2-H7."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_mw(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


def _patch_labels(dlg):
    """Suppress all 3D label calls on a dialog instance."""
    dlg.add_selection_label = MagicMock()
    dlg.clear_selection_labels = MagicMock()
    dlg.show_atom_labels_for = MagicMock()
    dlg.show_atom_labels = MagicMock()


# ---------------------------------------------------------------------------
# BondLengthDialog
# ---------------------------------------------------------------------------


@pytest.fixture
def bond_dlg(qapp):
    mol = _ethane()
    mw = _make_mw(mol)
    dlg = BondLengthDialog(mol, mw)
    _patch_labels(dlg)
    yield dlg, mol, mw
    try:
        dlg.picking_enabled = False
        dlg.close()
    except Exception:
        pass


class TestBondLengthPicking:
    def test_first_pick_sets_atom1(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.on_atom_picked(0)
        assert dlg.atom1_idx == 0
        assert dlg.atom2_idx is None

    def test_second_pick_sets_atom2(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        assert dlg.atom1_idx == 0
        assert dlg.atom2_idx == 1

    def test_third_pick_resets_to_new_atom1(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.on_atom_picked(2)
        assert dlg.atom1_idx == 2
        assert dlg.atom2_idx is None

    def test_clear_selection(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.clear_selection()
        assert dlg.atom1_idx is None
        assert dlg.atom2_idx is None

    def test_preselected_atoms_loaded(self, qapp):
        mol = _ethane()
        mw = _make_mw(mol)
        dlg = BondLengthDialog(mol, mw, preselected_atoms=[0, 1])
        _patch_labels(dlg)
        assert dlg.atom1_idx == 0
        assert dlg.atom2_idx == 1
        dlg.picking_enabled = False
        dlg.close()


class TestBondLengthIsComplete:
    def test_incomplete_with_one_atom(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        assert not dlg._is_selection_complete()

    def test_complete_with_two_atoms(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        assert dlg._is_selection_complete()


class TestBondLengthUpdateDisplay:
    def test_no_atoms_disables_apply(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()
        assert "No atoms" in dlg.selection_label.text()

    def test_one_atom_shows_symbol(self, bond_dlg):
        dlg, mol, _ = bond_dlg
        dlg.atom1_idx = 0
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()
        assert mol.GetAtomWithIdx(0).GetSymbol() in dlg.selection_label.text()

    def test_two_atoms_enables_apply_and_shows_distance(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.update_display()
        assert dlg.apply_button.isEnabled()
        assert "Å" in dlg.distance_label.text()


class TestBondLengthApplyChanges:
    def test_incomplete_selection_is_noop(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.apply_changes()  # must not raise

    def test_invalid_input_shows_warning(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.distance_input.setText("abc")
        with patch("moleditpy.ui.bond_length_dialog.QMessageBox") as mb:
            dlg.apply_changes()
        mb.warning.assert_called_once()

    def test_negative_distance_shows_warning(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.distance_input.setText("-1.0")
        with patch("moleditpy.ui.bond_length_dialog.QMessageBox") as mb:
            dlg.apply_changes()
        mb.warning.assert_called_once()

    def test_valid_apply_pushes_undo(self, bond_dlg):
        dlg, _, mw = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.distance_input.setText("1.54")
        dlg.apply_changes()
        mw.edit_actions_manager.push_undo_state.assert_called()


class TestBondLengthGeometry:
    def test_atom_only_mode_moves_atom2(self, bond_dlg):
        dlg, mol, _ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.atom1_fix_radio.setChecked(True)

        before_pos1 = np.array(mol.GetConformer().GetAtomPosition(0))
        dlg.apply_geometry_update(2.0)
        after_pos2 = np.array(mol.GetConformer().GetAtomPosition(1))
        after_pos1 = np.array(mol.GetConformer().GetAtomPosition(0))

        # atom1 must not move
        assert after_pos1 == pytest.approx(before_pos1, abs=1e-4)
        # distance atom1->atom2 must be ~2.0 Å
        dist = np.linalg.norm(after_pos2 - after_pos1)
        assert dist == pytest.approx(2.0, abs=1e-4)

    def test_default_mode_moves_group(self, bond_dlg):
        dlg, mol, _ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.atom1_fix_group_radio.setChecked(True)

        before_pos1 = np.array(mol.GetConformer().GetAtomPosition(0))
        dlg.apply_geometry_update(2.0)
        after_pos1 = np.array(mol.GetConformer().GetAtomPosition(0))
        after_pos2 = np.array(mol.GetConformer().GetAtomPosition(1))

        assert after_pos1 == pytest.approx(before_pos1, abs=1e-4)
        dist = np.linalg.norm(after_pos2 - after_pos1)
        assert dist == pytest.approx(2.0, abs=1e-4)

    def test_both_groups_mode(self, bond_dlg):
        dlg, mol, _ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.both_groups_radio.setChecked(True)

        dlg.apply_geometry_update(2.0)
        pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
        pos1 = np.array(mol.GetConformer().GetAtomPosition(1))
        dist = np.linalg.norm(pos1 - pos0)
        assert dist == pytest.approx(2.0, abs=1e-3)

    def test_on_distance_input_changed_syncs_slider(self, bond_dlg):
        dlg, *_ = bond_dlg
        dlg.atom1_idx = 0
        dlg.atom2_idx = 1
        dlg.apply_button.setEnabled(True)
        dlg.distance_slider.blockSignals(True)
        dlg.distance_slider.setValue(0)
        dlg.distance_slider.blockSignals(False)
        dlg.on_distance_input_changed("2.00")
        assert dlg.distance_slider.value() == 200


# ---------------------------------------------------------------------------
# AngleDialog
# ---------------------------------------------------------------------------


@pytest.fixture
def angle_dlg(qapp):
    mol = _ethane()
    mw = _make_mw(mol)
    dlg = AngleDialog(mol, mw)
    _patch_labels(dlg)
    yield dlg, mol, mw
    try:
        dlg.picking_enabled = False
        dlg.close()
    except Exception:
        pass


class TestAngleDialogPicking:
    def test_sequential_picking(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.on_atom_picked(2)
        assert dlg.atom1_idx == 2 and dlg.atom2_idx is None
        dlg.on_atom_picked(0)
        assert dlg.atom2_idx == 0 and dlg.atom3_idx is None
        dlg.on_atom_picked(1)
        assert dlg.atom3_idx == 1

    def test_fourth_pick_resets(self, angle_dlg):
        dlg, *_ = angle_dlg
        for i in (2, 0, 1):
            dlg.on_atom_picked(i)
        dlg.on_atom_picked(5)
        assert dlg.atom1_idx == 5
        assert dlg.atom2_idx is None
        assert dlg.atom3_idx is None

    def test_clear_selection(self, angle_dlg):
        dlg, *_ = angle_dlg
        for i in (2, 0, 1):
            dlg.on_atom_picked(i)
        dlg.clear_selection()
        assert dlg.atom1_idx is None
        assert dlg._snapshot_positions is None

    def test_preselected_atoms(self, qapp):
        mol = _ethane()
        mw = _make_mw(mol)
        dlg = AngleDialog(mol, mw, preselected_atoms=[2, 0, 1])
        _patch_labels(dlg)
        assert dlg.atom1_idx == 2
        assert dlg.atom2_idx == 0
        assert dlg.atom3_idx == 1
        dlg.picking_enabled = False
        dlg.close()


class TestAngleDialogIsComplete:
    def test_incomplete_with_two_atoms(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        assert not dlg._is_selection_complete()

    def test_complete_with_three_atoms(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        assert dlg._is_selection_complete()


class TestAngleDialogUpdateDisplay:
    def test_no_atoms_disables_apply(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_three_atoms_enables_apply(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.update_display()
        assert dlg.apply_button.isEnabled()
        assert "°" in dlg.angle_label.text()


class TestAngleDialogApplyChanges:
    def test_invalid_input_shows_warning(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.angle_input.setText("bad")
        with patch("moleditpy.ui.angle_dialog.QMessageBox") as mb:
            dlg.apply_changes()
        mb.warning.assert_called_once()

    def test_angle_wraps_over_180(self, angle_dlg):
        """270° input must wrap to -90°."""
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg._snapshot_positions = dlg.mol.GetConformer().GetPositions().copy()
        dlg.angle_input.setText("270")
        dlg.apply_changes()
        assert float(dlg.angle_input.text()) == pytest.approx(-90.0, abs=0.01)

    def test_apply_changes_pushes_undo(self, angle_dlg):
        dlg, _, mw = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg._snapshot_positions = dlg.mol.GetConformer().GetPositions().copy()
        dlg.angle_input.setText("109.5")
        dlg.apply_changes()
        mw.edit_actions_manager.push_undo_state.assert_called()

    def test_on_angle_input_changed_syncs_slider(self, angle_dlg):
        dlg, *_ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.apply_button.setEnabled(True)
        dlg.on_angle_input_changed("90")
        assert dlg.angle_slider.value() == 90


class TestAngleDialogGeometry:
    def test_rotate_atom_only_mode(self, angle_dlg):
        dlg, mol, _ = angle_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.rotate_atom_radio.setChecked(True)
        dlg._snapshot_positions = mol.GetConformer().GetPositions().copy()

        dlg.apply_geometry_update(120.0)

        from moleditpy.core.mol_geometry import calc_angle_deg

        conf = mol.GetConformer()
        new_angle = calc_angle_deg(
            conf.GetAtomPosition(2),
            conf.GetAtomPosition(0),
            conf.GetAtomPosition(1),
        )
        assert new_angle == pytest.approx(120.0, abs=0.5)


# ---------------------------------------------------------------------------
# DihedralDialog
# ---------------------------------------------------------------------------


@pytest.fixture
def dihedral_dlg(qapp):
    mol = _ethane()
    mw = _make_mw(mol)
    dlg = DihedralDialog(mol, mw)
    _patch_labels(dlg)
    yield dlg, mol, mw
    try:
        dlg.picking_enabled = False
        dlg.close()
    except Exception:
        pass


class TestDihedralDialogPicking:
    def test_sequential_picking(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        for i, (attr, expected) in enumerate(
            zip(["atom1_idx", "atom2_idx", "atom3_idx", "atom4_idx"], [2, 0, 1, 5])
        ):
            dlg.on_atom_picked(expected)
            assert getattr(dlg, attr) == expected

    def test_fifth_pick_resets(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        for idx in (2, 0, 1, 5):
            dlg.on_atom_picked(idx)
        dlg.on_atom_picked(3)
        assert dlg.atom1_idx == 3
        assert dlg.atom2_idx is None
        assert dlg.atom4_idx is None

    def test_clear_selection(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        for idx in (2, 0, 1, 5):
            dlg.on_atom_picked(idx)
        dlg.clear_selection()
        assert all(
            getattr(dlg, a) is None
            for a in ["atom1_idx", "atom2_idx", "atom3_idx", "atom4_idx"]
        )
        assert dlg._snapshot_positions is None

    def test_preselected_atoms(self, qapp):
        mol = _ethane()
        mw = _make_mw(mol)
        dlg = DihedralDialog(mol, mw, preselected_atoms=[2, 0, 1, 5])
        _patch_labels(dlg)
        assert dlg.atom1_idx == 2
        assert dlg.atom4_idx == 5
        dlg.picking_enabled = False
        dlg.close()


class TestDihedralDialogIsComplete:
    def test_incomplete_with_three_atoms(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        assert not dlg._is_selection_complete()

    def test_complete_with_four_atoms(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        assert dlg._is_selection_complete()


class TestDihedralDialogCalculate:
    def test_calculate_dihedral_incomplete_returns_zero(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        assert dlg.calculate_dihedral() == pytest.approx(0.0)

    def test_calculate_dihedral_complete_returns_value(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        val = dlg.calculate_dihedral()
        assert -180.0 <= val <= 180.0


class TestDihedralDialogUpdateDisplay:
    def test_no_atoms_disables_apply(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_four_atoms_enables_apply(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg.update_display()
        assert dlg.apply_button.isEnabled()
        assert "°" in dlg.dihedral_label.text()


class TestDihedralDialogApplyChanges:
    def test_invalid_input_shows_warning(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg.dihedral_input.setText("xyz")
        with patch("moleditpy.ui.dihedral_dialog.QMessageBox") as mb:
            dlg.apply_changes()
        mb.warning.assert_called_once()

    def test_dihedral_wraps_over_180(self, dihedral_dlg):
        """270° input must wrap to -90°."""
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg._snapshot_positions = dlg.mol.GetConformer().GetPositions().copy()
        dlg.dihedral_input.setText("270")
        dlg.apply_changes()
        assert float(dlg.dihedral_input.text()) == pytest.approx(-90.0, abs=0.01)

    def test_apply_changes_pushes_undo(self, dihedral_dlg):
        dlg, _, mw = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg._snapshot_positions = dlg.mol.GetConformer().GetPositions().copy()
        dlg.dihedral_input.setText("60.0")
        dlg.apply_changes()
        mw.edit_actions_manager.push_undo_state.assert_called()

    def test_on_dihedral_input_changed_syncs_slider(self, dihedral_dlg):
        dlg, *_ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg.apply_button.setEnabled(True)
        dlg.on_dihedral_input_changed("60")
        assert dlg.dihedral_slider.value() == 60


class TestDihedralDialogGeometry:
    def test_rotate_atom_only_sets_dihedral(self, dihedral_dlg):
        dlg, mol, _ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg.rotate_atom_radio.setChecked(True)
        dlg._snapshot_positions = mol.GetConformer().GetPositions().copy()

        dlg.apply_geometry_update(60.0)

        from moleditpy.core.mol_geometry import calculate_dihedral

        positions = mol.GetConformer().GetPositions()
        result = calculate_dihedral(positions, 2, 0, 1, 5)
        assert result == pytest.approx(60.0, abs=1.0)

    def test_default_group_mode_sets_dihedral(self, dihedral_dlg):
        dlg, mol, _ = dihedral_dlg
        dlg.atom1_idx = 2
        dlg.atom2_idx = 0
        dlg.atom3_idx = 1
        dlg.atom4_idx = 5
        dlg.rotate_group_radio.setChecked(True)
        dlg._snapshot_positions = mol.GetConformer().GetPositions().copy()

        dlg.apply_geometry_update(60.0)

        from moleditpy.core.mol_geometry import calculate_dihedral

        positions = mol.GetConformer().GetPositions()
        result = calculate_dihedral(positions, 2, 0, 1, 5)
        assert result == pytest.approx(60.0, abs=1.0)
