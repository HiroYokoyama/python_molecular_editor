"""
Unit tests for ConstrainedOptimizationDialog.

Covers:
  - Atom selection logic: on_atom_picked (toggle, max-4 FIFO eviction)
  - update_selection_display: label text and add-button state for 0-4 atoms
  - add_constraint: Distance / Angle / Dihedral type detection, duplicate guard
  - remove_constraint / remove_all_constraints: internal list + table sync
  - on_cell_changed: valid value edit, invalid input reverts to original
  - __init__ constraint loading: backward-compat 3-element (default force) and 4-element
  - reject: tuple→list JSON-safe serialisation saved to main_window
  - Force-field combo mapping from optimization_method string
"""

import os
import sys

import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtGui import QKeyEvent
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
        """Clicking an atom appends its index to selected_atoms."""
        dlg = make_dialog()
        dlg.on_atom_picked(3)
        assert dlg.selected_atoms == [3]

    def test_pick_same_atom_toggles_off(self, make_dialog):
        """Clicking the same atom twice removes it from the selection."""
        dlg = make_dialog()
        dlg.on_atom_picked(3)
        dlg.on_atom_picked(3)
        assert dlg.selected_atoms == []

    def test_pick_four_atoms(self, make_dialog):
        """Picking four distinct atoms fills selected_atoms to capacity."""
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
        """With no atoms selected the Add button is disabled and label shows 'None'."""
        dlg = make_dialog()
        assert not dlg.add_button.isEnabled()
        assert "None" in dlg.selection_label.text()

    def test_display_one_atom(self, make_dialog):
        """With one atom selected the Add button remains disabled."""
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        assert not dlg.add_button.isEnabled()

    def test_display_two_atoms_enables_add(self, make_dialog):
        """With two atoms selected the Add button is enabled and label shows 'Distance'."""
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        assert dlg.add_button.isEnabled()
        assert "Distance" in dlg.selection_label.text()

    def test_display_three_atoms(self, make_dialog):
        """With three atoms selected the label shows 'Angle'."""
        dlg = make_dialog()
        for i in range(3):
            dlg.on_atom_picked(i)
        assert dlg.add_button.isEnabled()
        assert "Angle" in dlg.selection_label.text()

    def test_display_four_atoms(self, make_dialog):
        """With four atoms selected the label shows 'Dihedral'."""
        dlg = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        assert dlg.add_button.isEnabled()
        assert "Dihedral" in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# add_constraint
# ---------------------------------------------------------------------------


class TestAddConstraint:
    def test_distance_constraint_added(self, make_dialog):
        """Adding two atoms creates a Distance constraint with a positive value and default force."""
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
        """add_constraint resets selected_atoms to an empty list after adding."""
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        assert dlg.selected_atoms == []

    def test_angle_constraint_added(self, make_dialog):
        """Adding three atoms creates an Angle constraint with value near tetrahedral angle."""
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

    def test_dihedral_constraint_added(self, make_dialog):
        """Adding four atoms creates a Dihedral constraint with correct atom indices."""
        dlg = make_dialog()
        # H2–C0–C1–H5: dihedral in ethane
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.on_atom_picked(5)
        dlg.add_constraint()

        assert len(dlg.constraints) == 1
        ctype, cidx, cval, _ = dlg.constraints[0]
        assert ctype == "Dihedral"
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
        """add_constraint inserts a row into constraint_table."""
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        assert dlg.constraint_table.rowCount() == 1

    def test_custom_force_constant_used(self, make_dialog):
        """Force constant entered by the user is stored in the constraint tuple."""
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
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox"):
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
        """remove_constraint deletes the selected row and its corresponding constraint."""
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
        """remove_all_constraints clears both the constraints list and the table."""
        dlg = make_dialog()
        self._add_distance(dlg, 0, 1)
        self._add_distance(dlg, 0, 2)
        dlg.remove_all_constraints()

        assert dlg.constraints == []
        assert dlg.constraint_table.rowCount() == 0

    def test_remove_all_on_empty_is_noop(self, make_dialog):
        """remove_all_constraints does not raise when called with no constraints."""
        dlg = make_dialog()
        dlg.remove_all_constraints()  # should not raise
        assert dlg.constraints == []

    def test_remove_button_disabled_after_remove_all(self, make_dialog):
        """Remove button is disabled after all constraints are cleared."""
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
        """Editing the value cell updates the constraint's stored value."""
        dlg = self._setup_one_constraint(make_dialog)
        original_type = dlg.constraints[0][0]

        self._set_cell_text(dlg, 0, 2, "2.000")
        dlg.on_cell_changed(0, 2)

        assert dlg.constraints[0][2] == pytest.approx(2.0)
        assert dlg.constraints[0][0] == original_type  # type unchanged

    def test_edit_force_column_updates_constraint(self, make_dialog):
        """Editing the force-constant cell updates the constraint's stored force value."""
        dlg = self._setup_one_constraint(make_dialog)
        self._set_cell_text(dlg, 0, 3, "2.5e4")
        dlg.on_cell_changed(0, 3)

        assert dlg.constraints[0][3] == pytest.approx(2.5e4)

    def test_invalid_value_reverts_to_original(self, make_dialog):
        """Entering a non-numeric cell value reverts the cell to the original constraint value."""
        dlg = self._setup_one_constraint(make_dialog)
        original_val = dlg.constraints[0][2]

        # Block signals so setText doesn't fire QMessageBox before patch is active
        self._set_cell_text(dlg, 0, 2, "xyz")
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox"):
            dlg.on_cell_changed(0, 2)

        # Internal value must be unchanged
        assert dlg.constraints[0][2] == pytest.approx(original_val)

    def test_non_value_column_ignored(self, make_dialog):
        """on_cell_changed does nothing when a non-editable column is touched."""
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
        """__init__ converts 4-element list constraints to tuples correctly."""
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
        """__init__ loads Distance, Angle, and Dihedral constraints simultaneously."""
        existing = [
            ["Distance", [0, 1], 1.54, 1.0e5],
            ["Angle", [2, 0, 1], 109.5, 5.0e4],
            ["Dihedral", [2, 0, 1, 5], 60.0, 1.0e5],
        ]
        dlg = make_dialog(constraints=existing)
        assert len(dlg.constraints) == 3
        assert dlg.constraint_table.rowCount() == 3

    def test_table_populated_on_load(self, make_dialog):
        """Preloaded constraints are shown in the table on dialog open."""
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
        """Force-field combo is pre-selected based on the current optimization method."""
        dlg = make_dialog(opt_method=opt_method)
        assert dlg.ff_combo.currentText() == expected_text

    def test_unknown_method_defaults_to_mmff94s(self, make_dialog):
        """An unrecognised optimization method defaults the force-field combo to MMFF94s."""
        dlg = make_dialog(opt_method="UNKNOWN_FORCE_FIELD")
        assert dlg.ff_combo.currentText() == "MMFF94s"


# ---------------------------------------------------------------------------
# reject — JSON-safe serialisation
# ---------------------------------------------------------------------------


class TestReject:
    def test_reject_converts_tuples_to_lists(self, make_dialog):
        """reject serialises constraints as JSON-safe lists, not tuples."""
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


# ---------------------------------------------------------------------------
# Optimization Execution (QThread)
# ---------------------------------------------------------------------------


class TestOptimizationExecution:
    @patch("moleditpy.ui.constrained_optimization_dialog.ConstrainedOptimizationThread")
    def test_apply_optimization_starts_thread(self, mock_thread_class, make_dialog):
        """apply_optimization disables the button and starts the worker thread."""
        dlg = make_dialog()
        mock_thread = MagicMock()
        mock_thread_class.return_value = mock_thread

        dlg.apply_optimization()

        # Check UI state
        assert not dlg.optimize_button.isEnabled()

        # Check thread was created and started
        mock_thread_class.assert_called_once()
        mock_thread.start.assert_called_once()

    def test_on_optimization_finished(self, make_dialog):
        """_on_optimization_finished re-enables the button and redraws the molecule."""
        dlg = make_dialog()
        dlg.optimize_button.setEnabled(False)

        # Mock conformer returning dummy coordinates
        mock_conf = MagicMock()
        pos_mock = MagicMock()
        pos_mock.x = 1.0
        pos_mock.y = 2.0
        pos_mock.z = 3.0
        mock_conf.GetAtomPosition.return_value = pos_mock

        # Setup mock cache
        dlg.main_window.view_3d_manager.atom_positions_3d = {
            i: [0, 0, 0] for i in range(dlg.mol.GetNumAtoms())
        }

        dlg._on_optimization_finished("MMFF94s", mock_conf)

        assert dlg.optimize_button.isEnabled()
        dlg.main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with(
            dlg.mol
        )
        dlg.main_window.edit_actions_manager.push_undo_state.assert_called_once()

        # Check coordinates were cached
        assert dlg.main_window.view_3d_manager.atom_positions_3d[0] == [1.0, 2.0, 3.0]

    def test_on_optimization_error(self, make_dialog):
        """_on_optimization_error re-enables the button and shows a critical dialog."""
        dlg = make_dialog()
        dlg.optimize_button.setEnabled(False)

        with patch(
            "moleditpy.ui.constrained_optimization_dialog.QMessageBox"
        ) as mock_mb:
            dlg._on_optimization_error("Test Error")

        assert dlg.optimize_button.isEnabled()
        mock_mb.critical.assert_called_once()
        assert "Test Error" in mock_mb.critical.call_args[0][2]


# ---------------------------------------------------------------------------
# Cancellation & undo integrity (dev-4.3.1 bug fixes)
# ---------------------------------------------------------------------------


class TestCancelAndUndoIntegrity:
    def _close_quietly(self, dlg):
        dlg.picking_enabled = False
        dlg.constraint_labels = []
        dlg.selection_labels = []
        dlg.reject()

    def test_finished_after_close_discards_result(self, make_dialog):
        """A result arriving after the dialog was closed must not touch the doc."""
        dlg = make_dialog()
        self._close_quietly(dlg)
        assert dlg._closed is True

        dlg.main_window.reset_mock()
        mock_conf = MagicMock()
        dlg._on_optimization_finished("MMFF94s", mock_conf)

        mock_conf.GetAtomPosition.assert_not_called()
        dlg.main_window.view_3d_manager.draw_molecule_3d.assert_not_called()
        dlg.main_window.edit_actions_manager.push_undo_state.assert_not_called()

    def test_error_after_close_is_suppressed(self, make_dialog):
        """An error arriving after close must not pop a message box."""
        dlg = make_dialog()
        self._close_quietly(dlg)

        with patch(
            "moleditpy.ui.constrained_optimization_dialog.QMessageBox"
        ) as mock_mb:
            dlg._on_optimization_error("late failure")

        mock_mb.critical.assert_not_called()

    def test_reject_pushes_undo_when_constraints_changed(self, make_dialog):
        """Closing the dialog after editing constraints must record an undo step."""
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = []
        self._close_quietly(dlg)

        dlg.main_window.edit_actions_manager.push_undo_state.assert_called_once()

    def test_reject_no_change_does_not_push_undo(self, make_dialog):
        """Closing without touching constraints must not create an undo entry."""
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = [
            ["Distance", [0, 1], 1.54, 1.0e5]
        ]
        self._close_quietly(dlg)

        dlg.main_window.edit_actions_manager.push_undo_state.assert_not_called()

    def test_finished_syncs_constraints_before_undo_push(self, make_dialog):
        """The undo snapshot must capture the constraints used for the run."""
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = []

        seen_at_push = []
        dlg.main_window.edit_actions_manager.push_undo_state.side_effect = (
            lambda: seen_at_push.append(
                list(dlg.main_window.edit_3d_manager.constraints_3d)
            )
        )

        mock_conf = MagicMock()
        pos_mock = MagicMock()
        pos_mock.x = pos_mock.y = pos_mock.z = 0.0
        mock_conf.GetAtomPosition.return_value = pos_mock
        dlg.main_window.view_3d_manager.atom_positions_3d = {
            i: [0, 0, 0] for i in range(dlg.mol.GetNumAtoms())
        }

        dlg._on_optimization_finished("MMFF94s", mock_conf)

        assert seen_at_push == [[["Distance", [0, 1], 1.54, 1.0e5]]]


# ---------------------------------------------------------------------------
# ConstrainedOptimizationThread
# ---------------------------------------------------------------------------


class TestOptimizationThread:
    def test_run_success_emits_finished(self, qapp):
        from moleditpy.ui.constrained_optimization_dialog import (
            ConstrainedOptimizationThread,
        )

        ff = MagicMock()
        thread = ConstrainedOptimizationThread(ff, max_iters=123)
        done = []
        thread.optimization_finished.connect(lambda: done.append(True))
        thread.run()

        ff.Minimize.assert_called_once_with(maxIts=123)
        assert done == [True]

    def test_run_error_emits_error(self, qapp):
        from moleditpy.ui.constrained_optimization_dialog import (
            ConstrainedOptimizationThread,
        )

        ff = MagicMock()
        ff.Minimize.side_effect = RuntimeError("boom")
        thread = ConstrainedOptimizationThread(ff)
        errs = []
        thread.error_occurred.connect(lambda m: errs.append(m))
        thread.run()

        assert errs and "boom" in errs[0]


# ---------------------------------------------------------------------------
# Force-field combo fallback branches
# ---------------------------------------------------------------------------


class TestForceFieldFallback:
    def test_legacy_uff_substring(self, make_dialog):
        dlg = make_dialog(opt_method="OLD_UFF_CONFIG")
        assert dlg.ff_combo.currentText() == "UFF"

    def test_legacy_mmff94s_substring(self, make_dialog):
        dlg = make_dialog(opt_method="LEGACY_MMFF94S_SET")
        assert dlg.ff_combo.currentText() == "MMFF94s"

    def test_legacy_mmff94_substring(self, make_dialog):
        dlg = make_dialog(opt_method="LEGACY_MMFF94_SET")
        assert dlg.ff_combo.currentText() == "MMFF94"

    def test_non_string_method_is_caught(self, make_dialog):
        # optimization_method as int -> .upper() raises -> except path
        dlg = make_dialog(opt_method=12345)
        assert dlg.ff_combo.currentText() in {"MMFF94s", "MMFF94", "UFF"}


# ---------------------------------------------------------------------------
# Guard clauses / early returns
# ---------------------------------------------------------------------------


class TestGuardClauses:
    def test_update_display_more_than_four(self, make_dialog):
        dlg = make_dialog()
        dlg.selected_atoms = [0, 1, 2, 3, 4]
        dlg.update_selection_display()
        assert "max 4" in dlg.selection_label.text()

    def test_add_constraint_with_one_atom_returns(self, make_dialog):
        dlg = make_dialog()
        dlg.selected_atoms = [0]
        dlg.add_constraint()
        assert dlg.constraints == []

    def test_remove_constraint_no_selection_returns(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        dlg.constraint_table.clearSelection()
        dlg.remove_constraint()
        assert len(dlg.constraints) == 1  # nothing removed

    def test_apply_optimization_no_conformer(self, make_dialog):
        mol = MagicMock()
        mol.GetNumConformers.return_value = 0
        dlg = make_dialog()
        dlg.mol = mol
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox") as mb:
            dlg.apply_optimization()
        mb.warning.assert_called_once()


# ---------------------------------------------------------------------------
# show_constraint_labels
# ---------------------------------------------------------------------------


class TestShowConstraintLabels:
    def _add_angle(self, dlg):
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()

    def test_no_selection_disables_remove(self, make_dialog):
        dlg = make_dialog()
        dlg.constraint_table.clearSelection()
        dlg.show_constraint_labels()
        assert not dlg.remove_button.isEnabled()

    def test_angle_labels_selected(self, make_dialog):
        dlg = make_dialog()
        dlg.main_window.view_3d_manager.atom_positions_3d = [
            [0.0, 0.0, 0.0] for _ in range(8)
        ]
        self._add_angle(dlg)
        dlg.constraint_table.selectRow(0)
        dlg.show_constraint_labels()
        assert dlg.remove_button.isEnabled()
        assert dlg.constraint_labels  # a label actor was added

    def test_dihedral_labels_and_3element_fallback(self, make_dialog):
        dlg = make_dialog()
        dlg.main_window.view_3d_manager.atom_positions_3d = [
            [0.0, 0.0, 0.0] for _ in range(8)
        ]
        # inject a legacy 3-element dihedral constraint
        dlg.constraints = [("Dihedral", (2, 0, 1, 5), 60.0)]
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.insertRow(0)
        dlg.constraint_table.setItem(0, 0, QTableWidgetItem("Dihedral"))
        dlg.constraint_table.blockSignals(False)
        dlg.constraint_table.selectRow(0)
        dlg.show_constraint_labels()
        assert dlg.constraint_labels

    def test_positions_none_returns_early(self, make_dialog):
        dlg = make_dialog()
        self._add_angle(dlg)
        dlg.main_window.view_3d_manager.atom_positions_3d = None
        dlg.constraint_labels = []
        dlg.constraint_table.selectRow(0)
        dlg.show_constraint_labels()
        assert dlg.constraint_labels == []

    def test_clear_labels_suppresses_remove_actor_error(self, make_dialog):
        dlg = make_dialog()
        dlg.constraint_labels = [MagicMock()]
        dlg.main_window.view_3d_manager.plotter.remove_actor.side_effect = RuntimeError(
            "no actor"
        )
        dlg.clear_constraint_labels()
        assert dlg.constraint_labels == []


# ---------------------------------------------------------------------------
# apply_optimization force-field + constraint plumbing
# ---------------------------------------------------------------------------


class TestApplyOptimizationPlumbing:
    def _add_all_types(self, dlg):
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.on_atom_picked(5)
        dlg.add_constraint()

    @patch("moleditpy.ui.constrained_optimization_dialog.ConstrainedOptimizationThread")
    def test_mmff_with_all_constraint_types(self, mock_thread, make_dialog):
        dlg = make_dialog()
        dlg.ff_combo.setCurrentText("MMFF94s")
        self._add_all_types(dlg)
        dlg.apply_optimization()
        mock_thread.assert_called_once()
        assert not dlg.optimize_button.isEnabled()

    @patch("moleditpy.ui.constrained_optimization_dialog.ConstrainedOptimizationThread")
    def test_uff_branch(self, mock_thread, make_dialog):
        dlg = make_dialog()
        dlg.ff_combo.setCurrentText("UFF")
        self._add_all_types(dlg)
        dlg.apply_optimization()
        mock_thread.assert_called_once()

    @patch("moleditpy.ui.constrained_optimization_dialog.ConstrainedOptimizationThread")
    def test_start_failure_reenables_button(self, mock_thread, make_dialog):
        dlg = make_dialog()
        mock_thread.return_value.start.side_effect = RuntimeError("cannot start")
        dlg.apply_optimization()
        assert dlg.optimize_button.isEnabled()


# ---------------------------------------------------------------------------
# _on_optimization_finished error/edge paths
# ---------------------------------------------------------------------------


class TestFinishedEdgePaths:
    def _finish(self, dlg):
        mock_conf = MagicMock()
        pos = MagicMock()
        pos.x = pos.y = pos.z = 0.0
        mock_conf.GetAtomPosition.return_value = pos
        dlg.main_window.view_3d_manager.atom_positions_3d = {
            i: [0, 0, 0] for i in range(dlg.mol.GetNumAtoms())
        }
        dlg._on_optimization_finished("MMFF94s", mock_conf)

    def test_3element_constraint_serialised(self, make_dialog):
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54)]  # legacy 3-element
        dlg.main_window.edit_3d_manager.constraints_3d = []
        self._finish(dlg)
        saved = dlg.main_window.edit_3d_manager.constraints_3d
        assert saved[0][3] == pytest.approx(1.0e5)

    def test_constraint_sync_error_is_caught(self, make_dialog):
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = []
        dlg.main_window.state_manager.update_window_title.side_effect = RuntimeError(
            "title fail"
        )
        # should not raise despite the sync error
        self._finish(dlg)
        dlg.main_window.view_3d_manager.draw_molecule_3d.assert_called_once()

    def test_outer_draw_error_is_caught(self, make_dialog):
        dlg = make_dialog()
        dlg.main_window.view_3d_manager.draw_molecule_3d.side_effect = RuntimeError(
            "draw fail"
        )
        # exception in the main body must be swallowed
        self._finish(dlg)


# ---------------------------------------------------------------------------
# reject error path + clear_selection + show_selection_labels
# ---------------------------------------------------------------------------


class TestRejectAndSelectionLabels:
    def test_reject_save_error_is_caught(self, make_dialog):
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.54, 1.0e5)]
        dlg.main_window.edit_3d_manager.constraints_3d = []
        dlg.main_window.state_manager.update_window_title.side_effect = RuntimeError(
            "boom"
        )
        dlg.picking_enabled = False
        dlg.constraint_labels = []
        dlg.selection_labels = []
        dlg.reject()  # must not raise
        assert dlg._closed is True

    def test_clear_selection(self, make_dialog):
        dlg = make_dialog()
        dlg.selected_atoms = [0, 1]
        dlg.clear_selection()
        assert dlg.selected_atoms == []
        assert "None" in dlg.selection_label.text()

    def test_show_selection_labels_with_real_positions(self, make_dialog):
        dlg = make_dialog()
        dlg.main_window.view_3d_manager.atom_positions_3d = [
            [0.0, 0.0, 0.0] for _ in range(8)
        ]
        dlg.selected_atoms = [0, 1]
        dlg.selection_labels = []
        dlg.show_selection_labels()
        assert dlg.selection_labels  # a label actor added

    def test_show_selection_labels_no_positions(self, make_dialog):
        dlg = make_dialog()
        dlg.main_window.view_3d_manager.atom_positions_3d = None
        dlg.selection_labels = []
        dlg.selected_atoms = [0]
        dlg.show_selection_labels()
        assert dlg.selection_labels == []


# ---------------------------------------------------------------------------
# on_cell_changed additional branches
# ---------------------------------------------------------------------------


class TestCellChangedBranches:
    def _one_angle(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(2)
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        return dlg

    def _set_text(self, dlg, row, col, text):
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.item(row, col).setText(text)
        dlg.constraint_table.blockSignals(False)

    def test_missing_item_returns(self, make_dialog):
        dlg = make_dialog()
        # No such row -> item() is None -> early return, no raise
        dlg.on_cell_changed(9, 2)

    def test_3element_force_column_edit(self, make_dialog):
        dlg = make_dialog()
        dlg.constraints = [("Distance", (0, 1), 1.5)]  # legacy 3-element
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.insertRow(0)
        dlg.constraint_table.setItem(0, 2, QTableWidgetItem("1.5"))
        dlg.constraint_table.setItem(0, 3, QTableWidgetItem("1.00e+05"))
        dlg.constraint_table.blockSignals(False)
        self._set_text(dlg, 0, 3, "3.0e4")
        dlg.on_cell_changed(0, 3)
        assert len(dlg.constraints[0]) == 4
        assert dlg.constraints[0][3] == pytest.approx(3.0e4)

    def test_invalid_angle_value_reverts_with_2dp(self, make_dialog):
        dlg = self._one_angle(make_dialog)
        orig = dlg.constraints[0][2]
        self._set_text(dlg, 0, 2, "notnum")
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox"):
            dlg.on_cell_changed(0, 2)
        assert dlg.constraints[0][2] == pytest.approx(orig)
        # value cell reverted to 2-decimal formatting for non-Distance
        assert dlg.constraint_table.item(0, 2).text() == f"{orig:.2f}"

    def test_invalid_force_value_reverts(self, make_dialog):
        dlg = self._one_angle(make_dialog)
        orig_force = dlg.constraints[0][3]
        self._set_text(dlg, 0, 3, "bogus")
        with patch("moleditpy.ui.constrained_optimization_dialog.QMessageBox"):
            dlg.on_cell_changed(0, 3)
        assert dlg.constraints[0][3] == pytest.approx(orig_force)
        assert dlg.constraint_table.item(0, 3).text() == f"{orig_force:.2e}"

    def test_index_error_when_no_constraint(self, make_dialog):
        dlg = make_dialog()
        # table row without a matching constraint entry -> IndexError path
        dlg.constraint_table.blockSignals(True)
        dlg.constraint_table.insertRow(0)
        dlg.constraint_table.setItem(0, 2, QTableWidgetItem("2.0"))
        dlg.constraint_table.blockSignals(False)
        dlg.on_cell_changed(0, 2)  # must not raise
        assert dlg.constraints == []


# ---------------------------------------------------------------------------
# keyPressEvent
# ---------------------------------------------------------------------------


class TestKeyPressEvent:
    def _key(self, key):
        return QKeyEvent(QEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier)

    def test_delete_removes_selected(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        dlg.constraint_table.selectRow(0)
        dlg.keyPressEvent(self._key(Qt.Key.Key_Delete))
        assert dlg.constraints == []

    def test_backspace_removes_selected(self, make_dialog):
        dlg = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(1)
        dlg.add_constraint()
        dlg.constraint_table.selectRow(0)
        dlg.keyPressEvent(self._key(Qt.Key.Key_Backspace))
        assert dlg.constraints == []

    @patch("moleditpy.ui.constrained_optimization_dialog.ConstrainedOptimizationThread")
    def test_enter_triggers_optimize(self, mock_thread, make_dialog):
        dlg = make_dialog()
        dlg.optimize_button.setEnabled(True)
        dlg.keyPressEvent(self._key(Qt.Key.Key_Return))
        mock_thread.assert_called_once()

    def test_other_key_falls_through(self, make_dialog):
        dlg = make_dialog()
        # An unrelated key must not raise and not remove anything
        dlg.keyPressEvent(self._key(Qt.Key.Key_A))
        assert dlg.constraints == []
