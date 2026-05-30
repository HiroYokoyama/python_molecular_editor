"""
Unit tests for MoveSelectedAtomsDialog (ui/move_selected_atoms_dialog.py).

Covers:
  - on_atom_picked: toggles clicked atom individually (no BFS)
  - apply_translation: no-selection guard, invalid input guard, correct coordinate update for selected atoms only
  - apply_rotation: no-selection guard, invalid input guard, rotation math (90° around centroid) for selected atoms only
  - reset_translation_inputs / reset_rotation_inputs
  - clear_selection: clears all state
"""

import os
import sys

import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

# ---------------------------------------------------------------------------
# Molecules
# ---------------------------------------------------------------------------


def _ethane():
    """Ethane with 3D coords: 2 C atoms (idx 0,1) + 6 H atoms (idx 2-7)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _make_main_window(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None):
        from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_main_window(_mol)
        # Patch show_atom_labels and clear_atom_labels to avoid pyvista calls
        with (
            patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
            patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
        ):
            dlg = MoveSelectedAtomsDialog(_mol, mw)
        entry = (dlg, _mol, mw)
        created.append(entry)
        return entry

    yield _factory

    for dlg, *_ in created:
        try:
            dlg.picking_enabled = False
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# on_atom_picked
# ---------------------------------------------------------------------------


class TestOnAtomPicked:
    def test_picks_atom_individually_no_bfs(self, make_dialog):
        """Picking an atom should select ONLY that atom (no BFS connected components)."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
        assert dlg.selected_atoms == {0}

    def test_picking_again_toggles_deselect(self, make_dialog):
        """Re-picking the same atom deselects it."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
            assert 0 in dlg.selected_atoms
            dlg.on_atom_picked(0)
            assert 0 not in dlg.selected_atoms


# ---------------------------------------------------------------------------
# apply_translation
# ---------------------------------------------------------------------------


class TestApplyTranslation:
    def test_no_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_trans_input.setText("bad")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_translation_updates_only_selected_atoms(self, make_dialog):
        dlg, mol, _ = make_dialog()
        dlg.selected_atoms.add(0)  # select carbon 0 only

        before = mol.GetConformer().GetPositions().copy()
        dlg.x_trans_input.setText("1.0")
        dlg.y_trans_input.setText("2.0")
        dlg.z_trans_input.setText("3.0")

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_translation()

        after = mol.GetConformer().GetPositions()
        # Atom 0 should have shifted
        assert after[0] == pytest.approx(before[0] + [1.0, 2.0, 3.0], abs=1e-4)
        # Atom 1 (unselected carbon) should NOT have shifted
        assert after[1] == pytest.approx(before[1], abs=1e-4)


# ---------------------------------------------------------------------------
# apply_rotation
# ---------------------------------------------------------------------------


class TestApplyRotation:
    def test_no_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_rot_input.setText("not_a_number")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_rotation_updates_only_selected_atoms_around_centroid(self, make_dialog):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import Geometry

        # Build a 2-atom molecule at known positions: [1,0,0] and [-1,0,0], and another at [5,5,5]
        mol = Chem.MolFromSmiles("[H][H].[H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(1.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Geometry.Point3D(-1.0, 0.0, 0.0))
        conf.SetAtomPosition(2, Geometry.Point3D(5.0, 5.0, 5.0))

        dlg, mol, mw = make_dialog(mol=mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [5.0, 5.0, 5.0]]
        )

        dlg.selected_atoms.update([0, 1])  # select atoms 0 and 1
        dlg.z_rot_input.setText("90.0")

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()

        after = mol.GetConformer().GetPositions()
        # Centroid of 0 and 1 is [0,0,0]. 90° rotation around Z moves them to [0,1,0] and [0,-1,0]
        assert after[0] == pytest.approx([0.0, 1.0, 0.0], abs=1e-5)
        assert after[1] == pytest.approx([0.0, -1.0, 0.0], abs=1e-5)
        # Atom 2 (unselected) should NOT rotate
        assert after[2] == pytest.approx([5.0, 5.0, 5.0], abs=1e-5)


# ---------------------------------------------------------------------------
# Reset inputs
# ---------------------------------------------------------------------------


class TestResetInputs:
    def test_reset_translation_inputs(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.x_trans_input.setText("5.5")
        dlg.y_trans_input.setText("-3.0")
        dlg.z_trans_input.setText("1.2")
        dlg.reset_translation_inputs()
        assert dlg.x_trans_input.text() == "0.0"
        assert dlg.y_trans_input.text() == "0.0"
        assert dlg.z_trans_input.text() == "0.0"

    def test_reset_rotation_inputs(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.x_rot_input.setText("45.0")
        dlg.y_rot_input.setText("90.0")
        dlg.z_rot_input.setText("180.0")
        dlg.reset_rotation_inputs()
        assert dlg.x_rot_input.text() == "0.0"
        assert dlg.y_rot_input.text() == "0.0"
        assert dlg.z_rot_input.text() == "0.0"


# ---------------------------------------------------------------------------
# clear_selection
# ---------------------------------------------------------------------------


class TestClearSelection:
    def test_clear_removes_selected_atoms(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.update([0, 1])
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0
