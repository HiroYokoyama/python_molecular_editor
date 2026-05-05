"""
Unit tests for ui/planarize_dialog.py (PlanarizeDialog).

Covers:
  - on_atom_picked: toggle select/deselect
  - preselected_atoms: loaded on init
  - clear_selection: empties set, disables apply
  - select_all_atoms: selects all N atoms
  - update_display: label text and apply_button state at 0 / 1-2 / >=3 atoms
  - apply_planarize: guard (<3 atoms), geometry math (atoms in a plane stay planar),
    undo is pushed
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
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None, preselected_atoms=None):
        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_mw(_mol)
        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            dlg = PlanarizeDialog(_mol, mw, preselected_atoms=preselected_atoms)
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
    def test_pick_adds_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.on_atom_picked(0)
        assert 0 in dlg.selected_atoms

    def test_repick_removes_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(0)
        assert 0 not in dlg.selected_atoms

    def test_multiple_picks_accumulate(self, make_dialog):
        dlg, _, _ = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        assert dlg.selected_atoms == {0, 1, 2, 3}


# ---------------------------------------------------------------------------
# preselected_atoms
# ---------------------------------------------------------------------------


class TestPreselectedAtoms:
    def test_preselected_atoms_loaded(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0, 1, 2])
        assert dlg.selected_atoms == {0, 1, 2}


# ---------------------------------------------------------------------------
# clear_selection
# ---------------------------------------------------------------------------


class TestClearSelection:
    def test_clear_empties_and_disables(self, make_dialog):
        dlg, _, _ = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0
        assert not dlg.apply_button.isEnabled()


# ---------------------------------------------------------------------------
# select_all_atoms
# ---------------------------------------------------------------------------


class TestSelectAllAtoms:
    def test_select_all_selects_every_atom(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.select_all_atoms()
        assert dlg.selected_atoms == set(range(mol.GetNumAtoms()))

    def test_select_all_enables_apply(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.select_all_atoms()
        assert dlg.apply_button.isEnabled()


# ---------------------------------------------------------------------------
# update_display
# ---------------------------------------------------------------------------


class TestUpdateDisplay:
    def test_zero_atoms_disables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_two_atoms_disables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_three_atoms_enables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1, 2}
        dlg.update_display()
        assert dlg.apply_button.isEnabled()

    def test_count_shown_in_label(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1, 2, 3}
        dlg.update_display()
        assert "4" in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# apply_planarize — guard
# ---------------------------------------------------------------------------


class TestApplyPlanarizeGuard:
    def test_fewer_than_three_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox") as mb:
            dlg.apply_planarize()
        mb.warning.assert_called_once()

    def test_zero_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.planarize_dialog.QMessageBox") as mb:
            dlg.apply_planarize()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# apply_planarize — geometry math
# ---------------------------------------------------------------------------


class TestApplyPlanarizeGeometry:
    def _planar_mol(self):
        """4 atoms in XY plane with small Z noise — planarize should reduce Z variance."""
        from rdkit import Chem, Geometry
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles("C1CCC1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=1)
        conf = mol.GetConformer()
        pts = [
            (1.0, 0.0, 0.1),
            (0.0, 1.0, -0.1),
            (-1.0, 0.0, 0.1),
            (0.0, -1.0, -0.1),
        ]
        for i, (x, y, z) in enumerate(pts):
            conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
        return mol

    def test_planarize_reduces_z_spread(self, make_dialog):
        mol = self._planar_mol()
        mw = _make_mw(mol)
        # atom_positions_3d must reflect the mol positions
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )

        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            dlg = PlanarizeDialog(mol, mw)

        dlg.selected_atoms = {0, 1, 2, 3}
        before_z_var = np.var(
            [mol.GetConformer().GetAtomPosition(i).z for i in range(4)]
        )

        with patch("moleditpy.ui.planarize_dialog.QMessageBox"):
            dlg.apply_planarize()

        after_z_var = np.var(
            [mol.GetConformer().GetAtomPosition(i).z for i in range(4)]
        )
        assert after_z_var <= before_z_var + 1e-6

    def test_planarize_pushes_undo(self, make_dialog):
        mol = self._planar_mol()
        mw = _make_mw(mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )

        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            dlg = PlanarizeDialog(mol, mw)

        dlg.selected_atoms = {0, 1, 2, 3}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox"):
            dlg.apply_planarize()
        mw.edit_actions_manager.push_undo_state.assert_called()

    def test_apply_calls_draw_molecule_3d(self, make_dialog):
        mol = self._planar_mol()
        mw = _make_mw(mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )

        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            dlg = PlanarizeDialog(mol, mw)

        dlg.selected_atoms = {0, 1, 2, 3}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox"):
            dlg.apply_planarize()
        mw.view_3d_manager.draw_molecule_3d.assert_called()
