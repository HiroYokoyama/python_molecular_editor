"""
Unit tests for ConstrainedOptimizationDialog.

Covers:
  - Atom selection logic: on_atom_picked (toggle, max-4 FIFO eviction)
  - update_selection_display: label text and add-button state for 0-4 atoms
  - add_constraint: Distance / Angle / Torsion type detection, duplicate guard
  - remove_constraint / remove_all_constraints: internal list + table sync
  - on_cell_changed: valid value edit, invalid input reverts to original
  - __init__ constraint loading: backward-compat 3-element (default force) and 4-element
  - reject: tuple→list JSON-safe serialisation saved to main_window
  - Force-field combo mapping from optimization_method string
"""

import os
import sys

import pytest
from unittest.mock import MagicMock, call, patch

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QTableWidgetItem

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ethane_mol():
    """Ethane with explicit H and a 3D conformer (8 atoms: C0, C1, H2-H7)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_main_window(constraints=None, opt_method="MMFF_RDKIT"):
    """Minimal mock of MainWindow for ConstrainedOptimizationDialog."""
    mw = MagicMock()
    mw.edit_3d_manager.constraints_3d = list(constraints or [])
    mw.init_manager.optimization_method = opt_method
    mw._picking_consumed = False
    return mw


# ---------------------------------------------------------------------------
# Session-scoped QApplication
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


# ---------------------------------------------------------------------------
# Dialog factory fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def make_dialog(qapp):
    """Return a factory that creates a fresh dialog with a given main_window."""
    created = []

    def _factory(constraints=None, opt_method="MMFF_RDKIT", mol=None):
        from moleditpy.ui.constrained_optimization_dialog import (
            ConstrainedOptimizationDialog,
        )

        _mol = mol if mol is not None else _make_ethane_mol()
        mw = _make_main_window(constraints=constraints, opt_method=opt_method)
        dlg = ConstrainedOptimizationDialog(_mol, mw)
        created.append(dlg)
        return dlg

    yield _factory

    for dlg in created:
        try:
            dlg.picking_enabled = False
            dlg.constraint_labels = []
            dlg.selection_labels = []
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Atom selection  (on_atom_picked / update_selection_display)
# ---------------------------------------------------------------------------


class TestAtomSelection:
    def test_pick_adds_atom(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(3)
        assert dlg.selected_atoms == [3]

    def test_pick_same_atom_toggles_off(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(3)
        dlg.on_atom_picked(3)
        assert dlg.selected_atoms == []

    def test_pick_four_atoms(self, make_dialog):
        dlg = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        assert dlg.selected_atoms == [0, 1, 2, 3]

    def test_fifth_atom_evicts_first_fifo(self, make_dialog):
        """Adding a 5th atom should drop the oldest (index 0)."""
        dlg = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        dlg.on_atom_picked(4)
        assert dlg.selected_atoms == [1, 2, 3, 4]

    def test_display_zero_atoms(self, make_dialog):
        dlg = make_dialog()
        assert not dlg.add_button.isEnabled()
        assert "None" in dlg.selection_label.text()

    def test_display_one_atom(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        assert not dlg.add_button.isEnabled()

    def test_display_two_atoms_enables_add(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        assert dlg.add_button.isEnabled()
        assert "Distance" in dlg.selection_label.text()

    def test_display_three_atoms(self, make_dialog):
        dlg = make_dialog()
        for i in range(3):
            dlg.on_atom_picked(i)
        assert dlg.add_button.isEnabled()
        assert "Angle" in dlg.selection_label.text()

    def test_display_four_atoms(self, make_dialog):
        dlg = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        assert dlg.add_button.isEnabled()
        assert "Torsion" in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# add_constraint
# ---------------------------------------------------------------------------


class TestAddConstraint:
    def test_distance_constraint_added(self, make_dialog):
        dlg = make_dialog()
        # atoms 0 and 1 are the bonded carbons in ethane
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, cforce = dlg.constraints[0]
        assert ctype == "Distance"
        assert cidx == (0, 1)
        assert cval > 0.0  # C–C bond length ~1.54 Å
        assert cforce == pytest.approx(1.0e5)

    def test_distance_constraint_clears_selection(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        assert dlg.selected_atoms == []

    def test_angle_constraint_added(self, make_dialog):
        dlg = make_dialog()
        # H2–C0–C1: a valid valence angle in ethane
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, _ = dlg.constraints[0]
        assert ctype == "Angle"
        assert cidx == (2, 0, 1)
        assert 90.0 <= cval <= 130.0  # tetrahedral angle ~109.5°

    def test_torsion_constraint_added(self, make_dialog):
        dlg = make_dialog()
        # H2–C0–C1–H5: dihedral in ethane
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.on_atom_picked(5)
        dlg.add_constraint()

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, _ = dlg.constraints[0]
        assert ctype == "Torsion"
        assert cidx == (2, 0, 1, 5)

    def test_duplicate_constraint_rejected(self, make_dialog):
        """Adding the same constraint twice should not append a second entry."""
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()

        # Same atoms again
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)

        with patch.object(dlg, "_QDialog__init__", create=True):
            with patch(
                "moleditpy.ui.constrained_optimization_dialog.QMessageBox"
            ) as mb:
                dlg.add_constraint()

        assert len(dlg.constraints) == 1

    def test_table_row_added_on_constraint(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        assert dlg.constraint_table.rowCount() == 1

    def test_custom_force_constant_used(self, make_dialog):
        dlg = make_dialog()
        dlg.force_const_input.setText("500.0")
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        _, _, _, cforce = dlg.constraints[0]
        assert cforce == pytest.approx(500.0)

    def test_invalid_force_constant_defaults(self, make_dialog):
        """Non-numeric force constant falls back to 1.0e5 (shows warning)."""
        dlg = make_dialog()
        dlg.force_const_input.setText("not_a_number")
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        with patch(
            "moleditpy.ui.constrained_optimization_dialog.QMessageBox"
        ):
            dlg.add_constraint()
        _, _, _, cforce = dlg.constraints[0]
        assert cforce == pytest.approx(1.0e5)


# ---------------------------------------------------------------------------
# remove_constraint / remove_all_constraints
# ---------------------------------------------------------------------------


class TestRemoveConstraint:
    def _add_distance(self, dlg, a, b):
        dlg.on_atom_picked(a)
        dlg.on_atom_picked(b)
        dlg.add_constraint()

    def test_remove_selected_row(self, make_dialog):
        dlg = make_dialog()
        self._add_distance(dlg, 0, 1)
        self._add_distance(dlg, 0, 2)

        dlg.constraint_table.selectRow(0)
        dlg.remove_constraint()

        assert len(dlg.constraints) == 1
        assert dlg.constraint_table.rowCount() == 1
        # Only the second constraint remains (0, 2)
        assert dlg.constraints[0][1] == (0, 2)

    def test_remove_all(self, make_dialog):
        dlg = make_dialog()
        self._add_distance(dlg, 0, 1)
        self._add_distance(dlg, 0, 2)
        dlg.remove_all_constraints()

        assert dlg.constraints == []
        assert dlg.constraint_table.rowCount() == 0

    def test_remove_all_on_empty_is_noop(self, make_dialog):
        dlg = make_dialog()
        dlg.remove_all_constraints()  # should not raise
        assert dlg.constraints == []

    def test_remove_button_disabled_after_remove_all(self, make_dialog):
        dlg = make_dialog()
        self._add_distance(dlg, 0, 1)
        dlg.remove_all_constraints()
        assert not dlg.remove_button.isEnabled()


# ---------------------------------------------------------------------------
# on_cell_changed — inline editing
# ---------------------------------------------------------------------------


class TestCellChanged:
    def _setup_one_constraint(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        return dlg

    def _set_cell_text(self, dlg, row, col, text):
        """Set cell text with signals blocked to avoid premature signal firing."""
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.item(row, col).setText(text)
        dlg.constraint_table.blockSignals(False)

    def test_edit_value_column_updates_constraint(self, make_dialog):
        dlg = self._setup_one_constraint(make_dialog)
        original_type = dlg.constraints[0][0]

        self._set_cell_text(dlg, 0, 2, "2.000")
        dlg.on_cell_changed(0, 2)

        assert dlg.constraints[0][2] == pytest.approx(2.0)
        assert dlg.constraints[0][0] == original_type  # type unchanged

    def test_edit_force_column_updates_constraint(self, make_dialog):
        dlg = self._setup_one_constraint(make_dialog)
        self._set_cell_text(dlg, 0, 3, "2.5e4")
        dlg.on_cell_changed(0, 3)

        assert dlg.constraints[0][3] == pytest.approx(2.5e4)

    def test_invalid_value_reverts_to_original(self, make_dialog):
        dlg = self._setup_one_constraint(make_dialog)
        original_val = dlg.constraints[0][2]

        # Block signals so setText doesn't fire QMessageBox before patch is active
        self._set_cell_text(dlg, 0, 2, "xyz")
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox"):
            dlg.on_cell_changed(0, 2)

        # Internal value must be unchanged
        assert dlg.constraints[0][2] == pytest.approx(original_val)

    def test_non_value_column_ignored(self, make_dialog):
        dlg = self._setup_one_constraint(make_dialog)
        original = dlg.constraints[0]
        dlg.on_cell_changed(0, 0)  # Type column — must be ignored
        assert dlg.constraints[0] == original

    def test_3element_constraint_upgraded_on_edit(self, make_dialog):
        """on_cell_changed must upgrade a legacy 3-element tuple to 4-element."""
        dlg = make_dialog()
        # Inject a 3-element constraint directly (old format)
        dlg.constraints = [("Distance", (0, 1), 1.5)]
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.insertRow(0)
        dlg.constraint_table.setItem(0, 2, QTableWidgetItem("1.5"))
        dlg.constraint_table.setItem(0, 3, QTableWidgetItem("1.00e+05"))
        dlg.constraint_table.blockSignals(False)

        self._set_cell_text(dlg, 0, 2, "1.8")
        dlg.on_cell_changed(0, 2)

        assert len(dlg.constraints[0]) == 4
        assert dlg.constraints[0][2] == pytest.approx(1.8)
        assert dlg.constraints[0][3] == pytest.approx(1.0e5)  # default force


# ---------------------------------------------------------------------------
# __init__  constraint loading
# ---------------------------------------------------------------------------


class TestConstraintLoading:
    def test_loads_4element_constraints(self, make_dialog):
        existing = [["Distance", [0, 1], 1.54, 2.5e4]]
        dlg = make_dialog(constraints=existing)

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, cforce = dlg.constraints[0]
        assert ctype == "Distance"
        assert cidx == (0, 1)
        assert cval == pytest.approx(1.54)
        assert cforce == pytest.approx(2.5e4)

    def test_loads_3element_constraints_with_default_force(self, make_dialog):
        """Legacy 3-element constraints should get default force 1.0e5."""
        existing = [["Angle", [2, 0, 1], 109.5]]
        dlg = make_dialog(constraints=existing)

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, cforce = dlg.constraints[0]
        assert ctype == "Angle"
        assert cval == pytest.approx(109.5)
        assert cforce == pytest.approx(1.0e5)

    def test_loads_multiple_constraint_types(self, make_dialog):
        existing = [
            ["Distance", [0, 1], 1.54, 1.0e5],
            ["Angle", [2, 0, 1], 109.5, 5.0e4],
            ["Torsion", [2, 0, 1, 5], 60.0, 1.0e5],
        ]
        dlg = make_dialog(constraints=existing)
        assert len(dlg.constraints) == 3
        assert dlg.constraint_table.rowCount() == 3

    def test_table_populated_on_load(self, make_dialog):
        existing = [["Distance", [0, 1], 1.54, 1.0e5]]
        dlg = make_dialog(constraints=existing)
        assert dlg.constraint_table.rowCount() == 1
        assert dlg.constraint_table.item(0, 0).text() == "Distance"


# ---------------------------------------------------------------------------
# Force-field combo mapping
# ---------------------------------------------------------------------------


class TestForceFieldMapping:
    @pytest.mark.parametrize(
        "opt_method,expected_text",
        [
            ("UFF_RDKIT", "UFF"),
            ("MMFF94_RDKIT", "MMFF94"),
            ("MMFF_RDKIT", "MMFF94s"),
        ],
    )
    def test_ff_combo_set_from_opt_method(self, make_dialog, opt_method, expected_text):
        dlg = make_dialog(opt_method=opt_method)
        assert dlg.ff_combo.currentText() == expected_text

    def test_unknown_method_defaults_to_mmff94s(self, make_dialog):
        dlg = make_dialog(opt_method="UNKNOWN_FORCE_FIELD")
        assert dlg.ff_combo.currentText() == "MMFF94s"


# ---------------------------------------------------------------------------
# reject — JSON-safe serialisation
# ---------------------------------------------------------------------------


class TestReject:
    def test_reject_converts_tuples_to_lists(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()

        # Disable picking to avoid side effects during reject
        dlg.picking_enabled = False
        dlg.constraint_labels = []
        dlg.selection_labels = []

        # Reset mock so we can inspect the assignment
        dlg.main_window.edit_3d_manager.constraints_3d = []

        with patch.object(
            type(dlg.main_window.state_manager), "has_unsaved_changes", create=True
        ):
            dlg.reject()

        saved = dlg.main_window.edit_3d_manager.constraints_3d
        assert isinstance(saved, list)
        assert isinstance(saved[0], list), "Constraint must be a list for JSON compat"
        assert isinstance(saved[0][1], list), "Atom indices must be a list"

    def test_reject_3element_constraint_gets_default_force(self, make_dialog):
        """Legacy 3-element tuples in constraints list must be serialised with default force."""
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54)]  # old 3-element format

        dlg.picking_enabled = False
        dlg.constraint_labels = []
        dlg.selection_labels = []
        dlg.main_window.edit_3d_manager.constraints_3d = []

        dlg.reject()

        saved = dlg.main_window.edit_3d_manager.constraints_3d
        assert saved[0][3] == pytest.approx(1.0e5)

    def test_reject_no_change_skips_state_update(self, make_dialog):
        """If saved constraints haven't changed, has_unsaved_changes must NOT be set."""
        dlg = make_dialog()
        # Pre-populate main_window with the same data that reject would produce
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = [
            ["Distance", [0, 1], 1.54, 1.0e5]
        ]

        dlg.picking_enabled = False
        dlg.constraint_labels = []
        dlg.selection_labels = []

        dlg.reject()

        # state_manager.has_unsaved_changes must NOT have been written
        dlg.main_window.state_manager.update_window_title.assert_not_called()
