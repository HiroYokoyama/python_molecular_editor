"""
Unit tests for MoveGroupDialog (ui/move_group_dialog.py).

Covers:
  - on_atom_picked: BFS connected-component traversal
  - Toggle behaviour: re-picking the same group deselects it
  - Multi-fragment: BFS stays within the picked fragment
  - update_display: label text (no group vs group with atoms)
  - apply_translation: no-group guard, invalid input guard, correct coordinate update
  - apply_rotation: no-group guard, invalid input guard, rotation math (Z-axis 90°)
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
    """Ethane with 3D coords: 2 C atoms (idx 0,1) + 6 H atoms (idx 2-7), all connected."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _two_methanes():
    """Two disconnected CH4 molecules: fragment A = atoms 0-4, fragment B = atoms 5-9."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("[CH4].[CH4]")
    mol = Chem.AddHs(mol)
    # useRandomCoords=True is required for reliable embedding of disconnected molecules
    AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
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
        from moleditpy.ui.move_group_dialog import MoveGroupDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_main_window(_mol)
        # Patch show_atom_labels and clear_atom_labels to avoid pyvista calls
        with (
            patch.object(MoveGroupDialog, "show_atom_labels"),
            patch.object(MoveGroupDialog, "clear_atom_labels"),
        ):
            dlg = MoveGroupDialog(_mol, mw)
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
# BFS / on_atom_picked
# ---------------------------------------------------------------------------


class TestOnAtomPicked:
    def test_picks_entire_connected_component(self, make_dialog):
        """Picking any atom in ethane should BFS to all 8 atoms."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
        assert len(dlg.group_atoms) == mol.GetNumAtoms()

    def test_picking_again_toggles_deselect(self, make_dialog):
        """Re-picking any atom in the already-selected group deselects everything."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
            dlg.on_atom_picked(1)  # still same group → deselect
        assert len(dlg.group_atoms) == 0

    def test_bfs_stays_within_fragment(self, make_dialog):
        """In a two-fragment system, BFS must not cross to the other fragment."""
        from rdkit.Chem import GetMolFrags

        mol = _two_methanes()
        dlg, mol, _ = make_dialog(mol=mol)

        # Determine the two fragments by RDKit — do not assume fixed index ranges
        frags = GetMolFrags(mol)  # tuple of tuples of atom indices
        assert len(frags) == 2, "Expected 2 disconnected fragments"
        frag_a = set(frags[0])
        frag_b = set(frags[1])

        # Pick any atom from fragment A
        seed = next(iter(frag_a))
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(seed)

        assert len(dlg.group_atoms) == len(frag_a)
        assert dlg.group_atoms == frag_a
        # No atom from fragment B must be selected
        assert dlg.group_atoms.isdisjoint(frag_b)

    def test_selected_atoms_records_clicked_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(3)
        assert 3 in dlg.selected_atoms

    def test_skip_pick_when_dragging(self, make_dialog):
        """on_atom_picked must be a no-op while dragging."""
        dlg, _, _ = make_dialog()
        dlg.is_dragging_group = True
        with patch.object(type(dlg), "show_atom_labels") as mock_labels:
            dlg.on_atom_picked(0)
        mock_labels.assert_not_called()
        assert len(dlg.group_atoms) == 0


# ---------------------------------------------------------------------------
# update_display
# ---------------------------------------------------------------------------


class TestUpdateDisplay:
    def test_no_group_shows_placeholder(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.group_atoms.clear()
        dlg.update_display()
        assert "No group" in dlg.selection_label.text()

    def test_group_shows_count_and_symbols(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)  # selects all 8 ethane atoms

        dlg.update_display()
        text = dlg.selection_label.text()
        assert "8" in text or "atoms" in text.lower()

    def test_more_than_5_atoms_appended_with_ellipsis(self, make_dialog):
        """If >5 atoms selected, display must show '...' at the end."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)  # 8-atom ethane

        dlg.update_display()
        assert "..." in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# apply_translation
# ---------------------------------------------------------------------------


class TestApplyTranslation:
    def _pick_all(self, dlg):
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

    def test_no_group_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_group_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        self._pick_all(dlg)
        dlg.x_trans_input.setText("bad")
        with patch("moleditpy.ui.move_group_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_translation_updates_conformer_positions(self, make_dialog):
        dlg, mol, _ = make_dialog()
        self._pick_all(dlg)

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
        # Every atom in the group should have shifted by [1, 2, 3]
        for idx in dlg.group_atoms:
            assert after[idx] == pytest.approx(before[idx] + [1.0, 2.0, 3.0], abs=1e-4)

    def test_translation_pushes_undo(self, make_dialog):
        dlg, _, mw = make_dialog()
        self._pick_all(dlg)
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_translation()
        mw.edit_actions_manager.push_undo_state.assert_called()


# ---------------------------------------------------------------------------
# apply_rotation
# ---------------------------------------------------------------------------


class TestApplyRotation:
    def _pick_all(self, dlg):
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

    def test_no_group_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_group_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        self._pick_all(dlg)
        dlg.x_rot_input.setText("not_a_number")
        with patch("moleditpy.ui.move_group_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_zero_rotation_leaves_positions_unchanged(self, make_dialog):
        dlg, mol, _ = make_dialog()
        self._pick_all(dlg)
        before = mol.GetConformer().GetPositions().copy()
        # All rotation inputs default to "0.0"
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()
        after = mol.GetConformer().GetPositions()
        assert after == pytest.approx(before, abs=1e-5)

    def test_90deg_z_rotation_around_centroid(self, make_dialog):
        """90° rotation around Z maps (centroid+[1,0,0]) → (centroid+[0,1,0])."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import Geometry

        # Build a 2-atom molecule at known positions: [1,0,0] and [-1,0,0]
        mol = Chem.MolFromSmiles("[H][H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(1.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Geometry.Point3D(-1.0, 0.0, 0.0))

        dlg, mol, mw = make_dialog(mol=mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
        )

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

        dlg.z_rot_input.setText("90.0")
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()

        after = mol.GetConformer().GetPositions()
        # Centroid = [0,0,0]; 90° Rz: [1,0,0]→[0,1,0], [-1,0,0]→[0,-1,0]
        assert after[0] == pytest.approx([0.0, 1.0, 0.0], abs=1e-5)
        assert after[1] == pytest.approx([0.0, -1.0, 0.0], abs=1e-5)

    def test_rotation_pushes_undo(self, make_dialog):
        dlg, _, mw = make_dialog()
        self._pick_all(dlg)
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()
        mw.edit_actions_manager.push_undo_state.assert_called()


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
    def test_clear_removes_group_and_selected_atoms(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()

        assert len(dlg.group_atoms) == 0
        assert len(dlg.selected_atoms) == 0

    def test_clear_resets_drag_state(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.is_dragging_group = True
        dlg.drag_start_pos = (10, 20)
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert not dlg.is_dragging_group
        assert dlg.drag_start_pos is None

    def test_clear_updates_display(self, make_dialog):
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert "No group" in dlg.selection_label.text()
