"""
Unit tests for ui/align_plane_dialog.py (AlignPlaneDialog).

Covers:
  - on_atom_picked: toggle select/deselect, update_display side-effects
  - preselected_atoms: loaded on init
  - clear_selection: empties set and disables apply
  - select_all_atoms: selects all N atoms from mol
  - update_display: label text and apply_button state at 0 / 1-2 / >=3 atoms
  - apply_PlaneAlign: guard (<3 atoms), plane-normal math (atoms in XY plane → XY align leaves z small)
  - apply_PlaneAlign pushes undo via _update_molecule_geometry path
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

    def _factory(mol=None, plane="xy", preselected_atoms=None):
        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_mw(_mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(
                _mol, mw, plane=plane, preselected_atoms=preselected_atoms
            )
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
    def test_clear_empties_selection(self, make_dialog):
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
# apply_PlaneAlign — guard
# ---------------------------------------------------------------------------


class TestApplyPlaneAlignGuard:
    def test_fewer_than_three_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        with patch("moleditpy.ui.align_plane_dialog.QMessageBox") as mb:
            dlg.apply_PlaneAlign()
        mb.warning.assert_called_once()

    def test_zero_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.align_plane_dialog.QMessageBox") as mb:
            dlg.apply_PlaneAlign()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# apply_PlaneAlign — geometry math
# ---------------------------------------------------------------------------


class TestApplyPlaneAlignMath:
    def _flat_mol_in_xy(self):
        """4-atom molecule lying exactly in the XY plane with noise."""
        from rdkit import Chem, Geometry
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles("C1CCC1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=1)
        # Override first 4 atoms to lie in XY plane (z=0)
        conf = mol.GetConformer()
        xy_pts = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
        ]
        for i, (x, y, z) in enumerate(xy_pts):
            conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
        return mol

    def test_xy_align_reduces_z_variance(self, make_dialog):
        """Atoms already in XY plane aligned to XY: z-coords should stay ~0."""
        mol = self._flat_mol_in_xy()
        mw = MagicMock()
        mw._picking_consumed = False
        positions = mol.GetConformer().GetPositions()
        mw.view_3d_manager.atom_positions_3d = positions.copy()

        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(mol, mw, plane="xy")

        dlg.selected_atoms = {0, 1, 2, 3}
        dlg.apply_PlaneAlign()

        after = mol.GetConformer().GetPositions()
        # All 4 atoms should have near-zero z after alignment
        for i in range(4):
            assert abs(after[i][2]) < 0.5

    def test_apply_calls_draw_molecule_3d(self, make_dialog):
        dlg, mol, mw = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"), patch.object(
            type(dlg), "select_all_atoms"
        ):
            dlg.selected_atoms = set(range(mol.GetNumAtoms()))
        dlg.apply_PlaneAlign()
        mw.view_3d_manager.draw_molecule_3d.assert_called()
