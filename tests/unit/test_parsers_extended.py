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
import moleditpy.ui.io_logic as mwm

mwm.QDialog = PyQt6.QtWidgets.QDialog
mwm.QInputDialog = PyQt6.QtWidgets.QInputDialog
mwm.QFileDialog = PyQt6.QtWidgets.QFileDialog
mwm.QMessageBox = PyQt6.QtWidgets.QMessageBox

from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.io_logic import IOManager
from PyQt6.QtCore import QPointF, QTimer
from PyQt6.QtWidgets import QWidget, QApplication, QDialog


class DummyParser(IOManager):
    def __init__(self, host):
        self._host = host
        IOManager.__init__(self, host)
        
        # Required attributes (as properties to reflect host state)
        self.chem_check_failed = False
        self.chem_check_tried = False
        self.is_xyz_derived = False

    def __getattr__(self, name):
        # Prefer properties, then delegate to host
        return getattr(self._host, name)

    @property
    def state_manager(self): return self.host.state_manager
    @property
    def init_manager(self): return self.host.init_manager
    @property
    def view_3d_manager(self): return self.host.view_3d_manager
    @property
    def ui_manager(self): return self.host.ui_manager
    @property
    def edit_actions_manager(self): return self.host.edit_actions_manager
    @property
    def edit_3d_manager(self): return self.host.edit_3d_manager

    @property
    def data(self): return self.host.state_manager.data
    @property
    def scene(self): return self.host.init_manager.scene
    @property
    def settings(self): return self.host.init_manager.settings
    @property
    def view_2d(self): return self.host.init_manager.view_2d
    @property
    def plotter(self): return self.host.view_3d_manager.plotter

    @property
    def current_mol(self): return self.host.view_3d_manager.current_mol
    @current_mol.setter
    def current_mol(self, v): self.host.view_3d_manager.current_mol = v

    @property
    def current_file_path(self): return getattr(self._host, "current_file_path", None)
    @current_file_path.setter
    def current_file_path(self, v): self._host.current_file_path = v

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

    def _apply_chem_check_and_set_flags(self, mol, source_desc=None, force_skip=False):
        pass

    def prompt_for_charge(self):
        """Default: cancel (return None) so tests don't block on real dialog."""
        return None, False, False

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
    """Verify fallback to ForwardSDMolSupplier when standard MolBlock reading fails."""
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
    """Verify that charge dialog is shown when loading XYZ."""
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "ask.xyz"
    xyz_path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = None
    mock_rd_mod.DetermineBonds.return_value = None
    parser.settings["skip_chemistry_checks"] = False
    parser.host.prompt_for_charge = lambda: (1, True, False)
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.GetIntProp("_xyz_charge") == 1


def test_load_xyz_charge_loop_cancel(mock_parser_host, tmp_path):
    """Verify handling of user cancellation during the charge input loop."""
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "cancel.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    # Mock prompt to return cancel (ok=False)
    parser.host.prompt_for_charge.return_value = (0, False, False)
    result = parser.load_xyz_file(str(path))
    assert result is None


def test_load_xyz_unrecognized_symbol(mock_parser_host, tmp_path):
    """Test load_xyz_file raises ValueError for unrecognized element symbols."""
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "unknown.xyz"
    xyz_path.write_text("1\nUnknown\nXx 0.0 0.0 0.0\n")
    parser.settings["skip_chemistry_checks"] = False
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is None
    # Verify status bar message (approximate check)
    msgs = [str(c.args[0]) for c in parser.statusBar().showMessage.call_args_list if c.args]
    assert any("Unrecognized element symbol" in m for m in msgs)


def test_save_as_xyz_logic(mock_parser_host, tmp_path):
    """Verify saving a molecule as an XYZ file."""
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
    """Verify that malformed V2000 headers are fixed automatically during MOL load."""
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
    """Verify that the XYZ load recovery loop handles retries correctly."""
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "retry.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    # Ensure it fails first, then succeeds after user input
    parser.host.prompt_for_charge.side_effect = [(0, True, False), (0, True, False)]
    mol = parser.load_xyz_file(str(path))
    assert mol is not None
    assert mol.GetIntProp("_xyz_charge") == 1


def test_load_mol_file_not_found(mock_parser_host):
    """Verify error handling when a MOL file is not found."""
    parser = DummyParser(mock_parser_host)
    parser.load_mol_file("missing_parser_xyz_final.mol")
    parser.statusBar().showMessage.assert_any_call(
        "File not found: missing_parser_xyz_final.mol"
    )


def test_load_mol_file_invalid_format(mock_parser_host, tmp_path):
    """Verify error handling for invalid MOL file content."""
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
    """Verify that charge and multiplicity are written to XYZ file headers."""
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
    """Verify fixing of malformed counts lines in MOL files."""
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
    """Verify complex error recovery branches during XYZ loading."""
    parser = DummyParser(mock_parser_host)
    path = tmp_path / "complex.xyz"
    path.write_text("1\nC\nC 0 0 0\n")
    mock_rd_mod.DetermineBonds.side_effect = RuntimeError("Fail")
    # First call fails with charge 1, second call chooses skip chemistry
    calls = iter([(1, True, False), (0, True, True)])
    parser.host.prompt_for_charge = lambda: next(calls)
    mol = parser.load_xyz_file(str(path))
    assert mol is not None
    assert mol.HasProp("_xyz_skip_checks") or getattr(mol, "_xyz_skip_checks", False)


def test_save_as_mol_no_current_path(mock_parser_host, tmp_path):
    """Verify saving as MOL when no current file path is set."""
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
    """Verify skip chemistry check flag is set when user chooses to skip."""
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "c2.xyz"
    xyz_path.write_text("2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n")
    parser.settings["skip_chemistry_checks"] = True
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.HasProp("_xyz_skip_checks") or getattr(mol, "_xyz_skip_checks", False)



