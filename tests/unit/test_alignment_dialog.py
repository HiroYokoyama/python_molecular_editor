"""
Unit tests for ui/alignment_dialog.py (AlignmentDialog).

Covers:
  - on_atom_picked: toggle (select/deselect), max-2 cap, update_display side-effects
  - preselected_atoms: loaded on init
  - clear_selection: empties set and disables apply
  - update_display: label text at 0/1/2 atoms, apply_button enabled state
  - apply_alignment: guard (< 2 atoms), Rodrigues rotation math (atom1→origin, atom2→axis)
  - apply_alignment pushes undo
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
    mw.view_3d_manager.current_mol = mol
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


def _patch_labels(dlg):
    """Context manager that suppresses pyvista label calls on the dialog."""
    from unittest.mock import patch as _patch

    return _patch.object(type(dlg), "add_selection_label"), _patch.object(
        type(dlg), "clear_selection_labels"
    )


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None, axis="x", preselected_atoms=None):
        from moleditpy.ui.alignment_dialog import AlignmentDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_mw(_mol)
        with (
            patch.object(AlignmentDialog, "add_selection_label"),
            patch.object(AlignmentDialog, "clear_selection_labels"),
            patch.object(AlignmentDialog, "enable_picking"),
            patch.object(AlignmentDialog, "disable_picking"),
        ):
            dlg = AlignmentDialog(_mol, mw, axis=axis, preselected_atoms=preselected_atoms)
        created.append(dlg)
        return dlg, _mol, mw

    yield _factory

    for dlg in created:
        try:
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# on_atom_picked
# ---------------------------------------------------------------------------


class TestOnAtomPicked:
    def test_first_pick_adds_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
        assert 0 in dlg.selected_atoms

    def test_second_pick_adds_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(1)
        assert dlg.selected_atoms == {0, 1}

    def test_third_pick_capped_at_two(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(1)
            dlg.on_atom_picked(2)  # must be ignored
        assert len(dlg.selected_atoms) == 2
        assert 2 not in dlg.selected_atoms

    def test_repick_deselects(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(0)
        assert 0 not in dlg.selected_atoms

    def test_enables_apply_when_two_atoms(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(1)
        assert dlg.apply_button.isEnabled()


# ---------------------------------------------------------------------------
# preselected_atoms
# ---------------------------------------------------------------------------


class TestPreselectedAtoms:
    def test_preselected_loaded(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0, 1])
        assert dlg.selected_atoms == {0, 1}

    def test_preselected_capped_at_two(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0, 1, 2])
        assert len(dlg.selected_atoms) == 2


# ---------------------------------------------------------------------------
# clear_selection
# ---------------------------------------------------------------------------


class TestClearSelection:
    def test_clears_atoms_and_disables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "add_selection_label"),
            patch.object(type(dlg), "clear_selection_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(1)
        with patch.object(type(dlg), "clear_selection_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0
        assert not dlg.apply_button.isEnabled()


# ---------------------------------------------------------------------------
# update_display
# ---------------------------------------------------------------------------


class TestUpdateDisplay:
    def test_zero_atoms_label(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_one_atom_label_contains_symbol(self, make_dialog):
        dlg, mol, _ = make_dialog()
        dlg.selected_atoms = {0}
        dlg.update_display()
        sym = mol.GetAtomWithIdx(0).GetSymbol()
        assert sym in dlg.selection_label.text()
        assert not dlg.apply_button.isEnabled()

    def test_two_atoms_enables_apply(self, make_dialog):
        dlg, mol, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        dlg.update_display()
        assert dlg.apply_button.isEnabled()
        assert "2" in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# apply_alignment — guard
# ---------------------------------------------------------------------------


class TestApplyAlignmentGuard:
    def test_fewer_than_two_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0}
        with patch("moleditpy.ui.alignment_dialog.QMessageBox") as mb:
            dlg.apply_alignment()
        mb.warning.assert_called_once()

    def test_zero_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.alignment_dialog.QMessageBox") as mb:
            dlg.apply_alignment()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# apply_alignment — geometry math
# ---------------------------------------------------------------------------


class TestApplyAlignmentMath:
    def _two_atom_mol(self, p1, p2):
        """Build a 2-atom H2 molecule with explicit positions."""
        from rdkit import Chem, Geometry
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles("[H][H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(*p1))
        conf.SetAtomPosition(1, Geometry.Point3D(*p2))
        return mol

    def _make_dlg(self, make_dialog, mol, axis):
        dlg, mol, mw = make_dialog(mol=mol, axis=axis)
        mw.view_3d_manager.current_mol = mol
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )
        dlg.selected_atoms = {0, 1}
        return dlg, mol, mw

    def test_x_axis_alignment_atom1_at_origin(self, make_dialog):
        """After X-alignment, atom1 must be at origin."""
        mol = self._two_atom_mol([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
        dlg, mol, _ = self._make_dlg(make_dialog, mol, "x")
        with patch("moleditpy.ui.alignment_dialog.QMessageBox"):
            dlg.apply_alignment()
        pos = mol.GetConformer().GetPositions()
        # atom1 (lower index after sort = idx 0) moved to origin
        assert pos[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-5)

    def test_x_axis_alignment_atom2_on_x_axis(self, make_dialog):
        """After X-alignment, atom2 must lie on positive X-axis (y=0, z=0)."""
        mol = self._two_atom_mol([0.0, 0.0, 0.0], [0.0, 3.0, 0.0])
        dlg, mol, _ = self._make_dlg(make_dialog, mol, "x")
        with patch("moleditpy.ui.alignment_dialog.QMessageBox"):
            dlg.apply_alignment()
        pos = mol.GetConformer().GetPositions()
        assert pos[1][1] == pytest.approx(0.0, abs=1e-5)
        assert pos[1][2] == pytest.approx(0.0, abs=1e-5)
        assert pos[1][0] > 0  # on positive X

    def test_z_axis_alignment_atom2_on_z_axis(self, make_dialog):
        """After Z-alignment, atom2 must lie on Z-axis (x=0, y=0)."""
        mol = self._two_atom_mol([0.0, 0.0, 0.0], [3.0, 0.0, 0.0])
        dlg, mol, _ = self._make_dlg(make_dialog, mol, "z")
        with patch("moleditpy.ui.alignment_dialog.QMessageBox"):
            dlg.apply_alignment()
        pos = mol.GetConformer().GetPositions()
        assert pos[1][0] == pytest.approx(0.0, abs=1e-5)
        assert pos[1][1] == pytest.approx(0.0, abs=1e-5)
        assert pos[1][2] > 0

    def test_apply_pushes_undo(self, make_dialog):
        mol = self._two_atom_mol([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        dlg, _, mw = self._make_dlg(make_dialog, mol, "x")
        with patch("moleditpy.ui.alignment_dialog.QMessageBox"):
            dlg.apply_alignment()
        mw.edit_actions_manager.push_undo_state.assert_called()
