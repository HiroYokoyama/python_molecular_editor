import os
import pytest
import sys
from unittest.mock import MagicMock, patch

# --- ABSOLUTE AGGRESSIVE PATCHING ---
# 1. Force rdDetermineBonds import to succeed by injecting a mock into sys.modules and rdkit.Chem
import rdkit.Chem

mock_rd_mod = MagicMock()
sys.modules["rdkit.Chem.rdDetermineBonds"] = mock_rd_mod
rdkit.Chem.rdDetermineBonds = mock_rd_mod

# 2. Patch real PyQt6 classes directly
import PyQt6.QtWidgets

ACCEPTED_CODE = 1
REJECTED_CODE = 0

PyQt6.QtWidgets.QDialog.exec = lambda self: REJECTED_CODE
PyQt6.QtWidgets.QDialog.exec_ = lambda self: REJECTED_CODE
PyQt6.QtWidgets.QDialog.show = lambda self: None
PyQt6.QtWidgets.QInputDialog.getText = lambda *args, **kwargs: ("0", False)
PyQt6.QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: ("", "")
PyQt6.QtWidgets.QFileDialog.getSaveFileName = lambda *args, **kwargs: ("", "")
PyQt6.QtWidgets.QMessageBox.critical = lambda *args, **kwargs: None
PyQt6.QtWidgets.QMessageBox.warning = lambda *args, **kwargs: None
PyQt6.QtWidgets.QMessageBox.information = lambda *args, **kwargs: None
PyQt6.QtWidgets.QMessageBox.question = lambda *args, **kwargs: 65536  # No

# Patch target module namespace
import moleditpy.modules.main_window_molecular_parsers as mwm

mwm.QDialog = PyQt6.QtWidgets.QDialog
mwm.QInputDialog = PyQt6.QtWidgets.QInputDialog
mwm.QFileDialog = PyQt6.QtWidgets.QFileDialog
mwm.QMessageBox = PyQt6.QtWidgets.QMessageBox

from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_molecular_parsers import MainWindowMolecularParsers
from PyQt6.QtCore import QPointF, QTimer
from PyQt6.QtWidgets import QWidget, QApplication, QDialog


class DummyParser(QWidget, MainWindowMolecularParsers):
    def __init__(self, host):
        super().__init__()
        self._host = host
        self.data = host.data
        self.scene = host.scene
        # Settings as real dict to avoid MagicMock behavior
        self.settings = {"always_ask_charge": False, "skip_chemistry_checks": False}
        self.view_2d = host.view_2d
        self.plotter = host.plotter
        self.statusBar_mock = MagicMock()
        self.current_mol = None
        self.is_xyz_derived = False
        self.dragged_atom_info = None
        self.current_file_path = None

    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self):
        return self.statusBar_mock

    def check_unsaved_changes(self):
        return True

    def clear_2d_editor(self, push_to_undo=True):
        pass

    def restore_ui_for_editing(self):
        pass

    def reset_undo_stack(self):
        pass

    def fit_to_view(self):
        pass

    def _apply_chem_check_and_set_flags(self, mol, source_desc=""):
        pass

    def estimate_bonds_from_distances(self, mol):
        if self.settings.get("force_fallback_fail", False):
            raise RuntimeError("FallbackFailed")
        return super().estimate_bonds_from_distances(mol)


# --- Safety Timer ---
@pytest.fixture(autouse=True)
def kill_all_dialogs():
    timer = QTimer()

    def killer():
        for w in QApplication.topLevelWidgets():
            if isinstance(w, QDialog):
                try:
                    w.reject()
                except:
                    pass

    timer.timeout.connect(killer)
    timer.start(50)
    yield
    timer.stop()


# --- Tests ---


def test_load_mol_file_fallback_to_sd_supplier(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol1 = Chem.MolFromSmiles("C")
    sdf_path = tmp_path / "test.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol1)
    writer.close()
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.load_mol_file(str(sdf_path))
    assert len(parser.data.atoms) == 1


def test_load_xyz_always_ask_charge(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.settings["always_ask_charge"] = True
    xyz_path = tmp_path / "ask.xyz"
    xyz_path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = None
    mock_rd_mod.DetermineBonds.return_value = None
    with patch.object(mwm.QInputDialog, "getText", return_value=("1", True)):
        with patch.object(mwm, "QDialog", side_effect=Exception("Force Fallback")):
            mol = parser.load_xyz_file(str(xyz_path))
            assert mol is not None
            assert mol.GetIntProp("_xyz_charge") == 1


def test_load_xyz_charge_loop_cancel(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "cancel.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = RuntimeError("Fail")
    parser.settings["force_fallback_fail"] = True
    with patch.object(mwm.QInputDialog, "getText", return_value=("", False)):
        with patch.object(mwm, "QDialog", side_effect=Exception("Force Fallback")):
            result = parser.load_xyz_file(str(path))
            assert result is None


def test_load_xyz_unrecognized_symbol(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "unknown.xyz"
    xyz_path.write_text("1\nUnknown\nXx 0.0 0.0 0.0\n")
    parser.settings["skip_chemistry_checks"] = False
    with pytest.raises(ValueError, match="Unrecognized element symbol"):
        parser.load_xyz_file(str(xyz_path))


def test_save_as_xyz_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.data.add_atom("O", QPointF(0, 0))
    mol = Chem.MolFromSmiles("O")
    AllChem.Compute2DCoords(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "saved.xyz")
    with patch.object(
        mwm.QFileDialog, "getSaveFileName", return_value=(save_path, "*.xyz")
    ):
        parser.save_as_xyz()
    assert os.path.exists(save_path)


def test_load_mol_file_with_v2000_fix(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    block = Chem.MolToMolBlock(mol).replace("V2000", "     ")
    mol_path = tmp_path / "fixable.mol"
    mol_path.write_text(block)
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.load_mol_file(str(mol_path))
    assert len(parser.data.atoms) == 1


def test_load_xyz_recovery_loop_retries(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "retry.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = [RuntimeError("Fail"), None]
    parser.settings["force_fallback_fail"] = True
    with patch.object(mwm.QInputDialog, "getText", return_value=("1", True)):
        with patch.object(mwm, "QDialog", side_effect=Exception("Force Fallback")):
            mol = parser.load_xyz_file(str(path))
            assert mol is not None
            assert mol.GetIntProp("_xyz_charge") == 1


def test_load_mol_file_not_found(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    parser.load_mol_file("missing_parser_xyz_final.mol")
    parser.statusBar().showMessage.assert_any_call(
        "File not found: missing_parser_xyz_final.mol"
    )


def test_load_mol_file_invalid_format(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    bad_mol = tmp_path / "bad.mol"
    bad_mol.write_text("GARBAGE")
    parser.load_mol_file(str(bad_mol))
    msgs = [
        str(c.args[0]) for c in parser.statusBar().showMessage.call_args_list if c.args
    ]
    assert any(
        "Invalid MOL" in m or "Failed to read" in m or "Error loading" in m
        for m in msgs
    )


def test_save_as_xyz_charge_mult(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("[C+]")
    AllChem.EmbedMolecule(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "charge.xyz")
    with patch.object(
        mwm.QFileDialog, "getSaveFileName", return_value=(save_path, "*.xyz")
    ):
        parser.save_as_xyz()
    with open(save_path, "r") as f:
        assert "chrg = 1" in f.read()


def test_load_mol_file_malformed_counts(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    block = Chem.MolToMolBlock(mol).splitlines()
    block[3] = "  1  0  0  0  0  0  0  0  0  0999 V2000"
    path = tmp_path / "counts_m.mol"
    path.write_text("\n".join(block))
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.load_mol_file(str(path))
    assert len(parser.data.atoms) == 1


def test_load_xyz_complex_recovery_branches(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "complex.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = RuntimeError("Fail")
    parser.settings["force_fallback_fail"] = True
    with patch.object(
        mwm.QInputDialog, "getText", side_effect=[("1", True), ("", False)]
    ):
        with patch.object(mwm, "QDialog", side_effect=Exception("Force Fallback")):
            assert parser.load_xyz_file(str(path)) is None
            msgs = [
                str(c.args[0])
                for c in parser.statusBar().showMessage.call_args_list
                if c.args
            ]
            assert any("failed for that charge" in m for m in msgs)


def test_save_as_mol_no_current_path(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.current_file_path = None
    save_path = str(tmp_path / "new_save.mol")
    parser.data.add_atom("C", QPointF(0, 0))
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    parser.current_mol = mol
    with patch.object(
        mwm.QFileDialog, "getSaveFileName", return_value=(save_path, "*.mol")
    ):
        parser.save_as_mol()
    assert os.path.exists(save_path)


def test_load_xyz_skip_chemistry_via_button(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "c2.xyz"
    xyz_path.write_text("2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n")
    parser.settings["skip_chemistry_checks"] = True
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.HasProp("_xyz_skip_checks") or getattr(mol, "_xyz_skip_checks", False)


def test_fix_mol_block(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    invalid_counts = " 3 2  0  0  0  0  0  0  0  0\n"
    mol_block = "\n  Title\n\n" + invalid_counts + "  0.0 0.0 0.0 C\n"
    fixed = parser.fix_mol_block(mol_block)
    lines = fixed.splitlines()
    assert "V2000" in lines[3]
    assert len(lines[3]) >= 39
