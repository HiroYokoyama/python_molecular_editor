"""
Unit tests for ui/align_plane_dialog.py (AlignPlaneDialog).

Shared BasePickingDialog contract tests (pick/select/display) live in
PickingDialogContractTests and run via TestAlignPlaneDialogContract.
"""

import os
import sys
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


def _flat_mol_in_xy():
    mol = Chem.MolFromSmiles("C1CCC1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate([(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]):
        conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
    return mol


# ---------------------------------------------------------------------------
# Shared BasePickingDialog contract
# ---------------------------------------------------------------------------


class TestAlignPlaneDialogContract(PickingDialogContractTests):
    @pytest.fixture
    def make_dialog(self, qapp):
        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        created = []

        def _factory(mol=None, plane="xy", preselected_atoms=None):
            _mol = mol if mol is not None else make_ethane()
            mw = make_mock_mw(_mol)
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
# apply_PlaneAlign — guard
# ---------------------------------------------------------------------------


class TestApplyPlaneAlignGuard:
    @pytest.fixture
    def dlg(self, qapp):
        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        mol = make_ethane()
        mw = make_mock_mw(mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            d = AlignPlaneDialog(mol, mw, plane="xy")
        yield d
        try:
            d.close()
        except Exception:
            pass

    def test_fewer_than_three_atoms_shows_warning(self, dlg):
        """Warning shown when fewer than 3 atoms are selected for plane alignment."""
        dlg.selected_atoms = {0, 1}
        with patch("moleditpy.ui.align_plane_dialog.QMessageBox") as mb:
            dlg.apply_PlaneAlign()
        mb.warning.assert_called_once()

    def test_zero_atoms_shows_warning(self, dlg):
        """Warning shown when no atoms are selected for plane alignment."""
        dlg.selected_atoms.clear()
        with patch("moleditpy.ui.align_plane_dialog.QMessageBox") as mb:
            dlg.apply_PlaneAlign()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# apply_PlaneAlign — geometry math
# ---------------------------------------------------------------------------


class TestApplyPlaneAlignMath:
    def test_xy_align_reduces_z_variance(self, qapp):
        """XY alignment brings all selected atom Z-coordinates close to zero."""
        mol = _flat_mol_in_xy()
        mw = make_mock_mw(mol)
        mw.view_3d_manager.atom_positions_3d = mol.GetConformer().GetPositions().copy()

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
        for i in range(4):
            assert abs(after[i][2]) < 0.5

    def test_apply_calls_draw_molecule_3d(self, qapp):
        """Plane alignment triggers a 3D redraw after modifying geometry."""
        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        mol = make_ethane()
        mw = make_mock_mw(mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(mol, mw, plane="xy")
        dlg.selected_atoms = set(range(mol.GetNumAtoms()))
        dlg.apply_PlaneAlign()
        mw.view_3d_manager.draw_molecule_3d.assert_called()

    def test_align_plane_with_move_to_zero_plane_true(self, qapp):
        """move_to_zero_plane=True translates selected atoms to Z=0 after alignment."""
        mol = _flat_mol_in_xy()
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x, p.y, p.z + 4.0))

        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        mw = make_mock_mw(mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(mol, mw, plane="xy")
        dlg.selected_atoms = [1, 2, 3]
        dlg.move_to_zero_plane_checkbox.setChecked(True)
        dlg.apply_PlaneAlign()
        pos = mol.GetConformer().GetPositions()
        for i in [1, 2, 3]:
            assert pos[i][2] == pytest.approx(0.0, abs=1e-5)

    def test_align_plane_with_move_to_zero_plane_false(self, qapp):
        """move_to_zero_plane=False keeps the centroid Z offset after alignment."""
        mol = _flat_mol_in_xy()
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x, p.y, p.z + 4.0))

        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        mw = make_mock_mw(mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(mol, mw, plane="xy")
        dlg.selected_atoms = [1, 2, 3]
        dlg.move_to_zero_plane_checkbox.setChecked(False)
        dlg.apply_PlaneAlign()
        pos = mol.GetConformer().GetPositions()
        for i in [1, 2, 3]:
            assert abs(pos[i][2]) > 1.0

    def test_align_plane_already_aligned_with_move_to_zero_plane_true(self, qapp):
        """Already-aligned plane with move_to_zero still zeroes the Z offset."""
        mol = _flat_mol_in_xy()
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x, p.y, 3.0))

        from moleditpy.ui.align_plane_dialog import AlignPlaneDialog

        mw = make_mock_mw(mol)
        with (
            patch.object(AlignPlaneDialog, "show_atom_labels"),
            patch.object(AlignPlaneDialog, "clear_atom_labels"),
            patch.object(AlignPlaneDialog, "enable_picking"),
            patch.object(AlignPlaneDialog, "disable_picking"),
        ):
            dlg = AlignPlaneDialog(mol, mw, plane="xy")
        dlg.selected_atoms = [1, 2, 3]
        dlg.move_to_zero_plane_checkbox.setChecked(True)
        dlg.apply_PlaneAlign()
        pos = mol.GetConformer().GetPositions()
        for i in [1, 2, 3]:
            assert pos[i][2] == pytest.approx(0.0, abs=1e-5)
