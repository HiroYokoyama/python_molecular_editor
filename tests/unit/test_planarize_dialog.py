"""
Unit tests for ui/planarize_dialog.py (PlanarizeDialog).

Shared BasePickingDialog contract tests (pick/select/display) live in
PickingDialogContractTests and run via TestPlanarizeDialogContract.
"""

import os
import sys
import numpy as np
import pytest
from unittest.mock import patch

from PyQt6.QtWidgets import QApplication
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from _picking_dialog_contract import (
    PickingDialogContractTests,
    make_ethane,
    make_mock_mw,
)


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _planar_mol():
    mol = Chem.MolFromSmiles("C1CCC1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate(
        [(1, 0, 0.1), (0, 1, -0.1), (-1, 0, 0.1), (0, -1, -0.1)]
    ):
        conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
    return mol


# ---------------------------------------------------------------------------
# Shared BasePickingDialog contract
# ---------------------------------------------------------------------------


class TestPlanarizeDialogContract(PickingDialogContractTests):
    @pytest.fixture
    def make_dialog(self, qapp):
        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        created = []

        def _factory(mol=None, preselected_atoms=None):
            _mol = mol if mol is not None else make_ethane()
            mw = make_mock_mw(_mol)
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
# apply_planarize — guard
# ---------------------------------------------------------------------------


class TestApplyPlanarizeGuard:
    @pytest.fixture
    def dlg(self, qapp):
        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        mol = make_ethane()
        mw = make_mock_mw(mol)
        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            d = PlanarizeDialog(mol, mw)
        yield d
        try:
            d.close()
        except Exception:
            pass

    def test_fewer_than_three_atoms_shows_warning(self, dlg):
        dlg.selected_atoms = {0, 1}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox") as mb:
            dlg.apply_planarize()
        mb.warning.assert_called_once()

    def test_zero_atoms_shows_warning(self, dlg):
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.planarize_dialog.QMessageBox") as mb:
            dlg.apply_planarize()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# apply_planarize — geometry math
# ---------------------------------------------------------------------------


class TestApplyPlanarizeGeometry:
    def _make_dlg(self, mol, qapp):
        from moleditpy.ui.planarize_dialog import PlanarizeDialog

        mw = make_mock_mw(mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )
        with (
            patch.object(PlanarizeDialog, "show_atom_labels"),
            patch.object(PlanarizeDialog, "clear_atom_labels"),
            patch.object(PlanarizeDialog, "enable_picking"),
            patch.object(PlanarizeDialog, "disable_picking"),
        ):
            return PlanarizeDialog(mol, mw), mw

    def test_planarize_reduces_z_spread(self, qapp):
        mol = _planar_mol()
        dlg, _ = self._make_dlg(mol, qapp)
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

    def test_planarize_pushes_undo(self, qapp):
        mol = _planar_mol()
        dlg, mw = self._make_dlg(mol, qapp)
        dlg.selected_atoms = {0, 1, 2, 3}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox"):
            dlg.apply_planarize()
        mw.edit_actions_manager.push_undo_state.assert_called()

    def test_apply_calls_draw_molecule_3d(self, qapp):
        mol = _planar_mol()
        dlg, mw = self._make_dlg(mol, qapp)
        dlg.selected_atoms = {0, 1, 2, 3}
        with patch("moleditpy.ui.planarize_dialog.QMessageBox"):
            dlg.apply_planarize()
        mw.view_3d_manager.draw_molecule_3d.assert_called()
