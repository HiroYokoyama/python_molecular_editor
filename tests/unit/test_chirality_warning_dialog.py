"""Unit tests for the chirality warning dialog and its post-conversion hook."""

from unittest.mock import MagicMock, patch

from PyQt6.QtCore import Qt

from moleditpy.core.stereo_check import ChiralityMismatch
from moleditpy.ui.chirality_warning_dialog import ChiralityWarningDialog
from moleditpy.ui.compute_logic import ComputeManager

MISMATCHES = [
    ChiralityMismatch(atom_id=3, symbol="C", rdkit_index=3, drawn="R", actual="S"),
    ChiralityMismatch(atom_id=13, symbol="C", rdkit_index=12, drawn="S", actual=None),
]


def _main_window(labels_on=False, menu_checked=False):
    mw = MagicMock()
    mw.view_3d_manager.show_chiral_labels = labels_on
    mw.view_3d_manager.current_mol = object()
    mw.init_manager.toggle_chiral_action.isChecked.return_value = menu_checked
    return mw


def test_dialog_is_non_modal_and_stays_on_top(app):
    """The warning never blocks the main window and stays above it."""
    dlg = ChiralityWarningDialog(_main_window(), MISMATCHES, 9, parent=None)
    try:
        assert dlg.isModal() is False
        assert dlg.windowFlags() & Qt.WindowType.WindowStaysOnTopHint
    finally:
        dlg.close()


def test_dialog_lists_every_mismatch(app):
    """One row per wrong center: atom, index, drawn and 3D configuration."""
    from PyQt6.QtWidgets import QLabel, QTableWidget

    dlg = ChiralityWarningDialog(_main_window(), MISMATCHES, 9, parent=None)
    try:
        table = dlg.findChild(QTableWidget)
        rows = [
            [table.item(r, c).text() for c in range(table.columnCount())]
            for r in range(table.rowCount())
        ]
        assert rows == [["C (ID 3)", "3", "R", "S"], ["C (ID 13)", "12", "S", "none"]]
        text = " ".join(lbl.text() for lbl in dlg.findChildren(QLabel))
        assert "2 of 9" in text
    finally:
        dlg.close()


def test_dialog_shows_chiral_labels_while_open(app):
    """Opening the warning turns 3D chiral labels on and redraws."""
    mw = _main_window(labels_on=False)
    dlg = ChiralityWarningDialog(mw, MISMATCHES, 9, parent=None)
    try:
        assert mw.view_3d_manager.show_chiral_labels is True
        mw.view_3d_manager.draw_molecule_3d.assert_called_once()
    finally:
        dlg.close()


def test_closing_restores_menu_setting(app):
    """Closing hands the labels back to the View menu (off here) and redraws."""
    mw = _main_window(labels_on=False, menu_checked=False)
    dlg = ChiralityWarningDialog(mw, MISMATCHES, 9, parent=None)
    mw.view_3d_manager.draw_molecule_3d.reset_mock()

    dlg.reject()

    assert mw.view_3d_manager.show_chiral_labels is False
    mw.view_3d_manager.draw_molecule_3d.assert_called_once()


def test_labels_already_on_are_left_alone(app):
    """If the user already shows chiral labels, neither open nor close redraws."""
    mw = _main_window(labels_on=True, menu_checked=True)
    dlg = ChiralityWarningDialog(mw, MISMATCHES, 9, parent=None)
    dlg.reject()

    assert mw.view_3d_manager.show_chiral_labels is True
    mw.view_3d_manager.draw_molecule_3d.assert_not_called()


def _compute(settings, mismatches):
    compute = ComputeManager.__new__(ComputeManager)
    compute.host = MagicMock()
    compute.host.init_manager.settings = settings
    compute.host.state_manager.data.atoms = {0: {"symbol": "C"}}
    compute._chirality_dialog = None
    patcher = patch(
        "moleditpy.ui.compute_logic.find_chirality_mismatches",
        return_value=mismatches,
    )
    return compute, patcher


def test_check_shows_dialog_on_mismatch(app):
    """A wrong 3D result opens the warning and says so in the status bar."""
    compute, patcher = _compute({}, MISMATCHES)
    with (
        patcher,
        patch("moleditpy.ui.compute_logic.drawn_chirality", return_value={1: "R"}),
        patch("moleditpy.ui.compute_logic.ChiralityWarningDialog") as dialog_cls,
    ):
        compute.check_chirality_against_2d()

    dialog_cls.assert_called_once()
    dialog_cls.return_value.show.assert_called_once()
    assert compute._chirality_dialog is dialog_cls.return_value
    assert "2 stereocenter(s)" in compute.host.update_status_message.call_args.args[0]


def test_check_is_quiet_when_structure_matches(app):
    """A correct 3D result shows nothing."""
    compute, patcher = _compute({}, [])
    with patcher, patch("moleditpy.ui.compute_logic.ChiralityWarningDialog") as dlg:
        compute.check_chirality_against_2d()
    dlg.assert_not_called()


def test_check_respects_setting(app):
    """With the option off, the structure is not even compared."""
    compute, patcher = _compute({"check_chirality_after_conversion": False}, MISMATCHES)
    with patcher as find, patch("moleditpy.ui.compute_logic.ChiralityWarningDialog"):
        compute.check_chirality_against_2d()
    find.assert_not_called()


def test_new_result_closes_previous_warning(app):
    """The warning about an older result is closed before the next check."""
    compute, patcher = _compute({}, [])
    old = MagicMock()
    compute._chirality_dialog = old
    with patcher:
        compute.check_chirality_against_2d()
    old.close.assert_called_once()
    assert compute._chirality_dialog is None


def test_conversion_finish_runs_the_check(app, mock_parser_host):
    """on_calculation_finished ends with the chirality check."""
    compute = ComputeManager(mock_parser_host)
    compute.check_chirality_against_2d = MagicMock()
    compute.on_calculation_finished(MagicMock())
    compute.check_chirality_against_2d.assert_called_once()
