"""Unit tests for multi-record SDF reading and the record selector."""

from unittest.mock import MagicMock, patch

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QDialog
from rdkit import Chem, RDLogger

from moleditpy.ui.io_logic import IOManager, SdfSelectionCancelled
from moleditpy.ui.sdf_record_dialog import SdfRecordDialog
from moleditpy.utils.sdf_records import read_sdf_records

RDLogger.DisableLog("rdApp.*")

UNREADABLE = (
    "broken\n  x\n\n"
    "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0 0.0 0.0 Xx  0\n"
    "M  END\n$$$$\n"
)


def _record(smiles, name, tag=None):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_Name", name)
    block = Chem.MolToMolBlock(mol)
    if tag:
        block += f"> <SOURCE>\n{tag}\n\n"
    return block + "$$$$\n"


def _sdf(*records):
    return "".join(records)


# --- utils.sdf_records ---------------------------------------------------


def test_reads_every_record_in_order():
    """All records come back, in file order, with name, formula and data fields."""
    text = _sdf(
        _record("CCO", "ethanol", tag="a"),
        _record("c1ccccc1", "benzene", tag="b"),
        _record("CC(=O)O", "acetic acid", tag="c"),
    )
    records = read_sdf_records(text)

    assert [r.display_name() for r in records] == ["ethanol", "benzene", "acetic acid"]
    assert [r.formula for r in records] == ["C2H6O", "C6H6", "C2H4O2"]
    assert [r.mol.GetProp("SOURCE") for r in records] == ["a", "b", "c"]
    assert [r.index for r in records] == [0, 1, 2]


def test_unreadable_record_keeps_its_place_and_title():
    """A record RDKit rejects stays in the list, named from its title line."""
    text = _sdf(_record("CCO", "ethanol"), UNREADABLE, _record("C", "methane"))
    records = read_sdf_records(text)

    assert len(records) == 3
    assert records[1].readable is False
    assert records[1].display_name() == "broken"
    assert records[1].formula == ""
    assert records[1].num_atoms == 0
    assert records[2].display_name() == "methane"


def test_trailing_blank_lines_add_no_record():
    """Whitespace after the last "$$$$" is not an extra, empty record."""
    text = _sdf(_record("CCO", "ethanol"), _record("C", "methane")) + "\n\n  \n"
    assert len(read_sdf_records(text)) == 2


def test_last_record_without_separator_is_read():
    """A file that omits the final "$$$$" still yields its last record."""
    text = _record("CCO", "ethanol") + Chem.MolToMolBlock(Chem.MolFromSmiles("C"))
    records = read_sdf_records(text)
    assert len(records) == 2 and records[1].readable


def test_blank_title_gets_a_numbered_name():
    """A record without a title is listed as "Molecule N"."""
    records = read_sdf_records(_sdf(_record("CCO", ""), _record("C", "")))
    assert [r.display_name() for r in records] == ["Molecule 1", "Molecule 2"]


def test_explicit_hydrogens_are_kept():
    """Records are read with removeHs=False, as the single-record path was."""
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    mol.SetProp("_Name", "methane")
    records = read_sdf_records(Chem.MolToMolBlock(mol) + "$$$$\n")
    assert records[0].num_atoms == 5


# --- ui.sdf_record_dialog -------------------------------------------------


def _dialog(app, text):
    return SdfRecordDialog(read_sdf_records(text))


def test_dialog_header_matches_plugin(app):
    """Same title and header line as the Open Babel Conversion Tool selector."""
    from PyQt6.QtWidgets import QLabel

    dlg = _dialog(app, _sdf(_record("CCO", "ethanol"), _record("C", "methane")))
    assert dlg.windowTitle() == "Select Molecule"
    texts = [lbl.text() for lbl in dlg.findChildren(QLabel)]
    assert "Found 2 molecules. Please select one:" in texts


def test_dialog_lists_number_name_and_formula(app):
    """Each row reads "N: name (formula)"; unreadable rows say so and are disabled."""
    dlg = _dialog(app, _sdf(_record("CCO", "ethanol"), UNREADABLE))
    rows = [dlg.list_widget.item(i) for i in range(dlg.list_widget.count())]
    assert rows[0].text() == "1: ethanol (C2H6O)"
    assert rows[1].text() == "2: broken (unreadable)"
    assert not rows[1].flags() & Qt.ItemFlag.ItemIsEnabled


def test_dialog_preselects_first_readable_record(app):
    """The first readable record is selected, skipping unreadable ones."""
    dlg = _dialog(app, _sdf(UNREADABLE, _record("CCO", "ethanol")))
    assert dlg.selected_record().display_name() == "ethanol"


def test_dialog_returns_chosen_record(app):
    """selected_record follows the current row."""
    dlg = _dialog(app, _sdf(_record("CCO", "ethanol"), _record("C", "methane")))
    dlg.list_widget.setCurrentRow(1)
    assert dlg.selected_record().display_name() == "methane"


def test_dialog_never_returns_unreadable_record(app):
    """Forcing the current row onto an unreadable record yields None."""
    dlg = _dialog(app, _sdf(_record("CCO", "ethanol"), UNREADABLE))
    dlg.list_widget.setCurrentRow(1)
    assert dlg.selected_record() is None


def test_double_click_accepts(app):
    """Double-clicking a row accepts the dialog."""
    dlg = _dialog(app, _sdf(_record("CCO", "ethanol"), _record("C", "methane")))
    dlg.list_widget.itemDoubleClicked.emit(dlg.list_widget.item(1))
    assert dlg.result() == QDialog.DialogCode.Accepted


# --- IOManager integration -----------------------------------------------


def _io():
    host = MagicMock()
    return IOManager(host)


def _write(tmp_path, text, name="multi.sdf"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_single_record_sdf_opens_no_dialog(tmp_path):
    """One record loads directly, as before."""
    path = _write(tmp_path, _record("CCO", "ethanol"))
    with patch("moleditpy.ui.io_logic.SdfRecordDialog") as dialog_cls:
        mol = _io()._read_mol_or_sdf(path)
    dialog_cls.assert_not_called()
    assert mol.GetProp("_Name") == "ethanol"


def test_multi_record_sdf_loads_the_chosen_record(tmp_path):
    """With several records the chosen one is loaded, data fields included."""
    path = _write(
        tmp_path,
        _sdf(_record("CCO", "ethanol", tag="a"), _record("C", "methane", tag="b")),
    )
    with patch("moleditpy.ui.io_logic.SdfRecordDialog") as dialog_cls:
        dialog = dialog_cls.return_value
        dialog.exec.return_value = QDialog.DialogCode.Accepted
        dialog.selected_record.side_effect = lambda: dialog_cls.call_args.args[0][1]
        mol = _io()._read_mol_or_sdf(path)

    assert mol.GetProp("_Name") == "methane"
    assert mol.GetProp("SOURCE") == "b"


@pytest.mark.parametrize("accepted", [False, True])
def test_cancel_or_empty_choice_raises_cancelled(tmp_path, accepted):
    """Cancel, or OK with nothing usable selected, is a cancellation."""
    path = _write(tmp_path, _sdf(_record("CCO", "a"), _record("C", "b")))
    with patch("moleditpy.ui.io_logic.SdfRecordDialog") as dialog_cls:
        dialog = dialog_cls.return_value
        dialog.exec.return_value = (
            QDialog.DialogCode.Accepted if accepted else QDialog.DialogCode.Rejected
        )
        dialog.selected_record.return_value = None
        with pytest.raises(SdfSelectionCancelled):
            _io()._read_mol_or_sdf(path)


def test_file_with_no_readable_record_skips_the_dialog(tmp_path):
    """Nothing to choose from: no dialog, and the old fallback path runs."""
    path = _write(tmp_path, _sdf(UNREADABLE, UNREADABLE))
    with patch("moleditpy.ui.io_logic.SdfRecordDialog") as dialog_cls:
        _io()._read_mol_or_sdf(path)
    dialog_cls.assert_not_called()


def test_non_widget_host_is_not_used_as_parent(tmp_path):
    """A host that is not a QWidget (test fakes) is not passed as Qt parent."""
    path = _write(tmp_path, _sdf(_record("CCO", "a"), _record("C", "b")))
    with patch("moleditpy.ui.io_logic.SdfRecordDialog") as dialog_cls:
        dialog_cls.return_value.exec.return_value = QDialog.DialogCode.Rejected
        with pytest.raises(SdfSelectionCancelled):
            _io()._read_mol_or_sdf(path)
    assert dialog_cls.call_args.args[1] is None


@pytest.mark.parametrize(
    "loader,message",
    [
        ("load_mol_file", "Import cancelled."),
        ("load_mol_file_for_3d_viewing", "Open cancelled."),
    ],
)
def test_cancelling_is_quiet(tmp_path, loader, message):
    """Cancelling the selector shows a status message, not an error dialog."""
    path = _write(tmp_path, _sdf(_record("CCO", "a"), _record("C", "b")))
    io = _io()
    io.host.state_manager.check_unsaved_changes.return_value = True
    with (
        patch.object(IOManager, "_read_mol_or_sdf", side_effect=SdfSelectionCancelled),
        patch.object(IOManager, "_report_load_error") as report,
    ):
        getattr(io, loader)(path)

    report.assert_not_called()
    io.host.statusBar.return_value.showMessage.assert_called_with(message)
    io.host.edit_actions_manager.clear_all.assert_not_called()
