# -*- coding: utf-8 -*-
"""E2E robustness: corrupted input files must never raise unhandled
exceptions out of the loaders — every case must land in a handler that
reports the failure (status bar + dialog in GUI mode) and leaves the app
alive. Added with the dev-4.3.1 fix for the empty-.pmeraw EOFError crash."""

import pickle

import pytest
from PyQt6.QtWidgets import QMessageBox

CASES = {
    "empty": b"",
    "binary_garbage": bytes(range(256)) * 4,
    "text_garbage": b"this is not a molecule file at all\nrandom text\n",
}

VALID_XYZ = b"3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.469\nH 0.0 -0.757 -0.469\n"


def _xyz_cases():
    c = dict(CASES)
    c["bad_count"] = b"abc\ncomment\nO 0.0 0.0 0.0\n"
    c["count_too_large"] = b"999\ncomment\nO 0.0 0.0 0.0\n"
    c["letters_as_coords"] = b"1\ncomment\nO x y z\n"
    c["truncated"] = VALID_XYZ[: len(VALID_XYZ) // 2]
    c["negative_count"] = b"-5\ncomment\nO 0.0 0.0 0.0\n"
    return c


def _mol_cases():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.Compute2DCoords(mol)
    block = Chem.MolToMolBlock(mol).encode()
    c = dict(CASES)
    c["truncated"] = block[: len(block) // 2]
    c["bad_counts_line"] = block.replace(b"V2000", b"XXXXX", 1)
    return c


def _pmeprj_cases():
    c = dict(CASES)
    c["invalid_json"] = b"{ not json !!"
    c["wrong_format_field"] = b'{"format": "Something Else", "version": "1.0"}'
    c["right_format_wrong_types"] = (
        b'{"format": "PME Project", "version": "1.0", "atoms": "garbage",'
        b' "bonds": 42, "structure_2d": [1, 2, 3]}'
    )
    return c


def _pmeraw_cases():
    c = dict(CASES)
    c["pickle_of_list"] = pickle.dumps([1, 2, 3])
    c["truncated_pickle"] = pickle.dumps({"a": 1})[:-3]
    return c


def _batteries(window):
    return [
        ("xyz3d", ".xyz", _xyz_cases(), window.io_manager.load_xyz_for_3d_viewing),
        ("mol2d", ".mol", _mol_cases(), window.io_manager.load_mol_file),
        ("mol3d", ".mol", _mol_cases(), window.io_manager.load_mol_file_for_3d_viewing),
        ("sdf3d", ".sdf", _mol_cases(), window.io_manager.load_mol_file_for_3d_viewing),
        ("pmeprj", ".pmeprj", _pmeprj_cases(), window.io_manager.load_json_data),
        ("pmeraw", ".pmeraw", _pmeraw_cases(), window.io_manager.load_raw_data),
    ]


@pytest.mark.gui
def test_corrupted_files_never_crash_loaders(window, qtbot, tmp_path, monkeypatch):
    for name in ("critical", "warning", "information", "question"):
        monkeypatch.setattr(
            QMessageBox,
            name,
            staticmethod(lambda *a, **k: QMessageBox.StandardButton.No),
        )
    monkeypatch.setattr(
        window.state_manager, "check_unsaved_changes", lambda: True, raising=False
    )
    window.init_manager.settings["skip_chemistry_checks"] = True

    crashes = []
    for battery, ext, cases, loader in _batteries(window):
        for case_name, payload in cases.items():
            p = tmp_path / f"{battery}_{case_name}{ext}"
            p.write_bytes(payload)
            try:
                loader(str(p))
                qtbot.wait(30)
            except Exception as e:  # noqa: BLE001 - collecting every escape
                crashes.append(f"{battery}:{case_name} -> {type(e).__name__}: {e}")

    assert crashes == [], "Unhandled exceptions escaped the loaders:\n" + "\n".join(
        crashes
    )


@pytest.mark.gui
def test_corrupt_startup_file_does_not_block_window_construction(
    app, qtbot, tmp_path, monkeypatch
):
    """Regression: initial_file was loaded synchronously in the constructor,
    so the invalid-format dialog for a '{}' .pmeprj blocked MainWindow from
    opening until dismissed. The load is now deferred to the event loop."""
    import importlib

    warnings_shown = []
    monkeypatch.setattr(
        QMessageBox, "warning", lambda *a, **k: warnings_shown.append(a)
    )
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)

    bad = tmp_path / "bad.pmeprj"
    bad.write_text("{}", encoding="utf-8")

    MainWindow = importlib.import_module("moleditpy.ui.main_window").MainWindow
    win = MainWindow(initial_file=str(bad), safe_mode=True)
    qtbot.addWidget(win)

    assert warnings_shown == []  # constructor finished without raising the dialog

    qtbot.wait(150)
    assert warnings_shown, "deferred load must still report the invalid file"
    win.close()
