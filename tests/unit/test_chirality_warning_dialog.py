"""Unit tests for the chirality warning dialog and its post-conversion hook."""

from unittest.mock import MagicMock, patch

import pytest

from PyQt6.QtCore import Qt

from moleditpy.core.stereo_check import ChiralityMismatch, EZMismatch
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


def test_labels_already_on_stay_on_but_lose_red_marks(app):
    """With labels already on, closing keeps them on and removes the red marks.

    Previously neither opening nor closing redrew in this case; the dialog now
    marks wrong centers, so both must redraw to add and then clear the marks.
    """
    mw = _main_window(labels_on=True, menu_checked=True)
    dlg = ChiralityWarningDialog(mw, MISMATCHES, 9, parent=None)
    assert mw.view_3d_manager.chirality_mismatches == {3: "S", 12: "?"}
    dlg.reject()

    assert mw.view_3d_manager.show_chiral_labels is True
    assert mw.view_3d_manager.chirality_mismatches == {}
    assert mw.view_3d_manager.draw_molecule_3d.call_count == 2


def _compute(settings, mismatches, ez_mismatches=()):
    compute = ComputeManager.__new__(ComputeManager)
    compute.host = MagicMock()
    compute.host.init_manager.settings = settings
    compute.host.state_manager.data.atoms = {0: {"symbol": "C"}}
    compute._chirality_dialog = None
    patcher = patch(
        "moleditpy.ui.compute_logic.find_chirality_mismatches",
        return_value=mismatches,
    )
    ez_patcher = patch(
        "moleditpy.ui.compute_logic.find_ez_mismatches",
        return_value=list(ez_mismatches),
    )
    ez_patcher.start()
    return compute, patcher


@pytest.fixture(autouse=True)
def _stop_patches():
    yield
    patch.stopall()


EZ = [
    EZMismatch(
        atom_ids=(4, 5),
        symbols=("C", "C"),
        rdkit_bond_index=4,
        rdkit_atom_indices=(4, 5),
        drawn="Z",
        actual="E",
    )
]


def test_check_shows_dialog_on_mismatch(app):
    """A wrong 3D result opens the warning and says so in the status bar."""
    compute, patcher = _compute({}, MISMATCHES, EZ)
    with (
        patcher,
        patch("moleditpy.ui.compute_logic.drawn_chirality", return_value={1: "R"}),
        patch("moleditpy.ui.compute_logic.drawn_ez", return_value={(4, 5): "Z"}),
        patch("moleditpy.ui.compute_logic.ChiralityWarningDialog") as dialog_cls,
    ):
        compute.check_chirality_against_2d()

    dialog_cls.assert_called_once()
    assert dialog_cls.call_args.kwargs["ez_mismatches"] == EZ
    assert dialog_cls.call_args.kwargs["total_double_bonds"] == 1
    dialog_cls.return_value.show.assert_called_once()
    assert compute._chirality_dialog is dialog_cls.return_value
    assert "3 stereo element(s)" in compute.host.update_status_message.call_args.args[0]


def test_ez_mismatch_alone_opens_the_warning(app):
    """A wrong double bond with every stereocenter right still warns."""
    compute, patcher = _compute({}, [], EZ)
    with (
        patcher,
        patch("moleditpy.ui.compute_logic.drawn_chirality", return_value={}),
        patch("moleditpy.ui.compute_logic.drawn_ez", return_value={(4, 5): "Z"}),
        patch("moleditpy.ui.compute_logic.ChiralityWarningDialog") as dialog_cls,
    ):
        compute.check_chirality_against_2d()
    dialog_cls.assert_called_once()


def test_check_is_quiet_when_structure_matches(app):
    """A correct 3D result shows nothing."""
    compute, patcher = _compute({}, [])
    with patcher, patch("moleditpy.ui.compute_logic.ChiralityWarningDialog") as dlg:
        compute.check_chirality_against_2d()
    dlg.assert_not_called()


def test_check_respects_setting(app):
    """With the conversion option off, the structure is not even compared."""
    compute, patcher = _compute({"check_stereo_after_conversion": False}, MISMATCHES)
    with patcher as find, patch("moleditpy.ui.compute_logic.ChiralityWarningDialog"):
        compute.check_chirality_against_2d()
    find.assert_not_called()


@pytest.mark.parametrize(
    "settings,runs",
    [
        ({}, False),  # after optimization: off by default
        ({"check_stereo_after_optimization": True}, True),
        # The conversion option does not switch on the optimization check.
        ({"check_stereo_after_conversion": True}, False),
    ],
)
def test_optimization_check_follows_its_own_setting(app, settings, runs):
    """After Optimize 3D the check runs only when its own option is on."""
    compute, patcher = _compute(settings, [])
    with patcher as find:
        compute.check_chirality_against_2d(after_optimization=True)
    assert find.called is runs


def test_new_result_closes_previous_warning(app):
    """The warning about an older result is closed before the next check."""
    compute, patcher = _compute({}, [])
    old = MagicMock()
    compute._chirality_dialog = old
    with patcher:
        compute.check_chirality_against_2d()
    old.close.assert_called_once()
    assert compute._chirality_dialog is None


def test_finish_without_worker_id_counts_as_conversion(app, mock_parser_host):
    """A result with no worker id (legacy/direct callers) is a conversion."""
    compute = ComputeManager(mock_parser_host)
    compute.check_chirality_against_2d = MagicMock()
    compute.on_calculation_finished(MagicMock())
    compute.check_chirality_against_2d.assert_called_once_with(after_optimization=False)


@pytest.mark.parametrize("is_conversion", [True, False])
def test_finish_tells_conversion_from_optimization(
    app, mock_parser_host, is_conversion
):
    """Conversion runs check as conversions; Optimize 3D runs as optimizations."""
    compute = ComputeManager(mock_parser_host)
    compute.check_chirality_against_2d = MagicMock()
    compute.active_worker_ids.add(7)
    if is_conversion:
        compute._conversion_run_ids.add(7)
    compute.on_calculation_finished((7, MagicMock()))

    compute.check_chirality_against_2d.assert_called_once_with(
        after_optimization=not is_conversion
    )
    assert 7 not in compute._conversion_run_ids


def test_halt_forgets_pending_conversions(app, mock_parser_host):
    """A halted conversion's id does not linger."""
    compute = ComputeManager(mock_parser_host)
    compute.active_worker_ids.add(3)
    compute._conversion_run_ids.add(3)
    compute.halt_conversion()
    assert compute._conversion_run_ids == set()


def test_plugin_optimizer_runs_the_optimization_check(app, mock_parser_host):
    """A plugin optimizer runs in-process, so it calls the check itself."""
    compute = ComputeManager(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = MagicMock()
    mock_parser_host.init_manager.optimization_method = "MY_OPT"
    mock_parser_host.init_manager.opt3d_method_labels = {}
    mock_parser_host.plugin_manager.optimization_methods = {"MY_OPT": {"label": "x"}}
    compute._run_plugin_optimization = MagicMock()
    compute.check_chirality_against_2d = MagicMock()

    compute.optimize_3d_structure()

    compute._run_plugin_optimization.assert_called_once()
    compute.check_chirality_against_2d.assert_called_once_with(after_optimization=True)


def test_dialog_lists_double_bonds(app):
    """E/Z rows read "C=C (ID a=b)", with RDKit atom indices joined by '='."""
    from PyQt6.QtWidgets import QTableWidget

    dlg = ChiralityWarningDialog(
        _main_window(), [], 0, parent=None, ez_mismatches=EZ, total_double_bonds=2
    )
    try:
        table = dlg.findChild(QTableWidget)
        row = [table.item(0, c).text() for c in range(table.columnCount())]
        assert row == ["C=C (ID 4=5)", "4=5", "Z", "E"]
        assert dlg.windowTitle() == "Stereochemistry Check"
    finally:
        dlg.close()


def test_dialog_marks_and_clears_double_bonds(app):
    """Open marks the wrong double bonds for the 3D view; close clears them."""
    mw = _main_window()
    dlg = ChiralityWarningDialog(mw, [], 0, parent=None, ez_mismatches=EZ)
    assert mw.view_3d_manager.ez_mismatches == {4: "E"}
    dlg.reject()
    assert mw.view_3d_manager.ez_mismatches == {}
