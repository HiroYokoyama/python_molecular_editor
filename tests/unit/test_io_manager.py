"""
Unit tests for IOManager IO logic.
Covers:
  - prompt_for_charge
  - load_xyz_for_3d_viewing
  - save_3d_as_mol

Follows the same conventions as test_edit_actions_extended.py:
  - Real PyQt6 via a session-scoped QApplication
  - DummyHost with all required sub-managers
  - Module-level patch targets (no local imports in production code)
"""

import os
import sys
import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication, QDialog

# Make the local moleditpy package discoverable
workspace_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(workspace_src_path) and workspace_src_path not in sys.path:
    sys.path.insert(0, workspace_src_path)

from moleditpy.ui.io_logic import IOManager
from rdkit import Chem
from rdkit.Chem import AllChem


# ---------------------------------------------------------------------------
# Shared QApplication fixture (session-scoped — one instance for all tests)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def qapp():
    app = QApplication.instance() or QApplication([])
    return app


# ---------------------------------------------------------------------------
# DummyHost — mirrors DummyHost in test_edit_actions_extended.py
# ---------------------------------------------------------------------------

class DummyHost:
    def __init__(self):
        self.statusBar_mock = MagicMock()
        self.settings = {}
        self.is_xyz_derived = False
        self.initialization_complete = True
        self._is_restoring_state = False

        # Child managers
        self.init_manager = MagicMock()
        self.init_manager.settings = self.settings
        self.init_manager.current_file_path = ""
        self.state_manager = MagicMock()
        self.state_manager.check_unsaved_changes.return_value = True
        self.state_manager.has_unsaved_changes = False
        self.state_manager.data = MagicMock()
        self.state_manager.data.atoms = {}
        self.ui_manager = MagicMock()
        self.view_3d_manager = MagicMock()
        self.view_3d_manager.current_mol = None
        self.view_3d_manager.atom_id_to_rdkit_idx_map = {}
        self.view_3d_manager.plotter = MagicMock()
        self.edit_actions_manager = MagicMock()
        self.edit_actions_manager.clear_all.return_value = True
        self.plugin_manager = MagicMock()

    def statusBar(self):
        return self.statusBar_mock

    def push_undo_state(self):
        pass


# ===========================================================================
# Tests for prompt_for_charge
# ===========================================================================

# Key insight: QDialog is now a module-level import in io_logic.py, so we
# patch "moleditpy.ui.io_logic.QDialog".  The method checks:
#     if dialog.exec() != QDialog.DialogCode.Accepted
# With our patch, both QDialog (class) and the dialog instance come from
# the MockDialog; so we must return MockDialog.DialogCode.Accepted from exec()
# to pass the "Accepted" branch.  Returning anything else → Rejected path.

class TestPromptForCharge:
    """Tests for IOManager.prompt_for_charge()."""

    def _make_io(self):
        host = DummyHost()
        return IOManager(host), host

    # Shared patcher names
    _PATCHES = [
        "moleditpy.ui.io_logic.QVBoxLayout",
        "moleditpy.ui.io_logic.QHBoxLayout",
        "moleditpy.ui.io_logic.QLabel",
        "moleditpy.ui.io_logic.QDialogButtonBox",
    ]

    def test_accept_with_default_charge(self, qapp):
        """Accept with default text '0' → (0, True, False)."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "0"

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton"):

            dlg = MagicMock()
            dlg.exec.return_value = MockDlg.DialogCode.Accepted  # same sentinel
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge == 0
        assert ok is True
        assert skip is False

    def test_accept_with_positive_charge(self, qapp):
        """Text '2' → charge 2."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "2"

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton"):

            dlg = MagicMock()
            dlg.exec.return_value = MockDlg.DialogCode.Accepted
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge == 2
        assert ok is True
        assert skip is False

    def test_accept_with_negative_charge(self, qapp):
        """Text '-1' → charge -1."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "-1"

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton"):

            dlg = MagicMock()
            dlg.exec.return_value = MockDlg.DialogCode.Accepted
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge == -1
        assert ok is True
        assert skip is False

    def test_invalid_text_falls_back_to_zero(self, qapp):
        """Non-numeric text 'abc' → charge falls back to 0."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "abc"

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton"):

            dlg = MagicMock()
            dlg.exec.return_value = MockDlg.DialogCode.Accepted
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge == 0
        assert ok is True
        assert skip is False

    def test_cancel_returns_none_false_false(self, qapp):
        """Rejected exec → (None, False, False)."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "0"

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton"):

            dlg = MagicMock()
            # Return something that is != MockDlg.DialogCode.Accepted
            dlg.exec.return_value = MockDlg.DialogCode.Rejected
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge is None
        assert ok is False
        assert skip is False

    def test_skip_chemistry_returns_zero_true_true(self, qapp):
        """Skip button fires callback before exec returns Accepted → (0, True, True)."""
        io, _ = self._make_io()
        le = MagicMock()
        le.text.return_value = "0"

        captured_cb = []

        class FakeSkipBtn:
            def __init__(self, *a, **k):
                pass
            @property
            def clicked(self):
                class _Sig:
                    def connect(self_, cb):
                        captured_cb.append(cb)
                return _Sig()

        with patch("moleditpy.ui.io_logic.QDialog") as MockDlg, \
             patch("moleditpy.ui.io_logic.QLineEdit", return_value=le), \
             patch("moleditpy.ui.io_logic.QVBoxLayout"), \
             patch("moleditpy.ui.io_logic.QHBoxLayout"), \
             patch("moleditpy.ui.io_logic.QLabel"), \
             patch("moleditpy.ui.io_logic.QDialogButtonBox"), \
             patch("moleditpy.ui.io_logic.QPushButton", FakeSkipBtn):

            dlg = MagicMock()
            accepted = MockDlg.DialogCode.Accepted  # capture before exec fires

            def _exec():
                if captured_cb:
                    captured_cb[0]()          # fire skip callback
                return accepted               # return the same sentinel

            dlg.exec.side_effect = _exec
            MockDlg.return_value = dlg

            charge, ok, skip = io.prompt_for_charge()

        assert charge == 0
        assert ok is True
        assert skip is True


# ===========================================================================
# Tests for load_xyz_for_3d_viewing
# ===========================================================================

class TestLoadXYZFor3DViewing:
    """Tests for IOManager.load_xyz_for_3d_viewing()."""

    @staticmethod
    def _dummy_mol():
        mol = MagicMock()
        mol.HasProp.return_value = False
        mol.GetNumBonds.return_value = 0
        return mol

    def test_explicit_path_skips_dialog(self, qapp, tmp_path):
        """Passing file_path directly bypasses the file dialog."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        mol = self._dummy_mol()
        io.load_xyz_file = MagicMock(return_value=mol)

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        io.load_xyz_file.assert_called_once_with(str(xyz))
        host.edit_actions_manager.clear_all.assert_called_once_with(skip_check=True)
        assert host.view_3d_manager.current_mol is mol

    def test_dialog_provides_path(self, qapp, tmp_path):
        """Without explicit path, QFileDialog is shown and its result used."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        mol = self._dummy_mol()
        io.load_xyz_file = MagicMock(return_value=mol)

        with patch("moleditpy.ui.io_logic.QFileDialog.getOpenFileName",
                   return_value=(str(xyz), None)), \
             patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing()

        io.load_xyz_file.assert_called_once_with(str(xyz))
        assert host.view_3d_manager.current_mol is mol

    def test_dialog_cancelled_is_noop(self, qapp):
        """Empty string from dialog → early return, no state change."""
        host = DummyHost()
        io = IOManager(host)
        io.load_xyz_file = MagicMock()

        with patch("moleditpy.ui.io_logic.QFileDialog.getOpenFileName",
                   return_value=("", None)):
            io.load_xyz_for_3d_viewing()

        io.load_xyz_file.assert_not_called()
        assert host.view_3d_manager.current_mol is None

    def test_load_failure_shows_error(self, qapp, tmp_path):
        """load_xyz_file returning None shows error on status bar."""
        xyz = tmp_path / "bad.xyz"
        xyz.write_text("0\n\n")

        host = DummyHost()
        io = IOManager(host)
        io.load_xyz_file = MagicMock(return_value=None)

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        host.statusBar_mock.showMessage.assert_called()
        assert host.view_3d_manager.current_mol is None

    def test_ui_modes_enabled_on_success(self, qapp, tmp_path):
        """3D viewer UI mode methods are called after successful load."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        mol = self._dummy_mol()
        io.load_xyz_file = MagicMock(return_value=mol)

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        host.ui_manager._enter_3d_viewer_ui_mode.assert_called_once()
        host.ui_manager._enable_3d_features.assert_called_once_with(True)

    def test_is_xyz_derived_with_skip_prop(self, qapp, tmp_path):
        """is_xyz_derived=True when _xyz_skip_checks property is set on mol."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        mol = MagicMock()
        mol.HasProp.return_value = True
        mol.GetIntProp.return_value = 1
        mol.GetNumBonds.return_value = 0
        io.load_xyz_file = MagicMock(return_value=mol)

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        assert host.is_xyz_derived is True

    def test_is_xyz_derived_from_zero_bonds(self, qapp, tmp_path):
        """is_xyz_derived=True when molecule has 0 bonds (no skip flag)."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        mol = MagicMock()
        mol.HasProp.return_value = False
        mol.GetNumBonds.return_value = 0
        io.load_xyz_file = MagicMock(return_value=mol)

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        assert host.is_xyz_derived is True

    def test_current_file_path_updated(self, qapp, tmp_path):
        """init_manager.current_file_path is updated to the loaded file."""
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("2\ntest\nC  0.0  0.0  0.0\nO  1.2  0.0  0.0\n")

        host = DummyHost()
        io = IOManager(host)
        io.load_xyz_file = MagicMock(return_value=self._dummy_mol())

        with patch("moleditpy.ui.io_logic.QTimer"):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))

        assert host.init_manager.current_file_path == str(xyz)


# ===========================================================================
# Tests for save_3d_as_mol
# ===========================================================================

class TestSave3DAsMol:
    """Tests for IOManager.save_3d_as_mol()."""

    @staticmethod
    def _make_io_with_mol(smiles="C"):
        host = DummyHost()
        io = IOManager(host)
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        host.view_3d_manager.current_mol = mol
        return io, host

    def test_no_mol_shows_error_no_dialog(self, qapp):
        """current_mol is None → error message, dialog never opened."""
        host = DummyHost()
        io = IOManager(host)
        host.view_3d_manager.current_mol = None

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName") as mock_dlg:
            io.save_3d_as_mol()

        mock_dlg.assert_not_called()
        host.statusBar_mock.showMessage.assert_called_with(
            "Error: No 3D structure to save."
        )

    def test_dialog_cancelled_writes_nothing(self, qapp, tmp_path):
        """Empty path from dialog → no file written."""
        io, host = self._make_io_with_mol()
        out = tmp_path / "out.mol"

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                   return_value=("", None)):
            io.save_3d_as_mol()

        assert not out.exists()

    def test_mol_file_written(self, qapp, tmp_path):
        """MOL file is created and non-empty."""
        io, host = self._make_io_with_mol("CCO")
        out = tmp_path / "out.mol"

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                   return_value=(str(out), None)):
            io.save_3d_as_mol()

        assert out.exists()
        assert out.stat().st_size > 0

    def test_header_line_replaced(self, qapp, tmp_path):
        """Second line is replaced with 'MoleditPy Ver. ... 3D'."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "header.mol"

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                   return_value=(str(out), None)):
            io.save_3d_as_mol()

        lines = out.read_text(encoding="utf-8").splitlines()
        assert len(lines) > 1
        assert "MoleditPy Ver." in lines[1]
        assert "3D" in lines[1]

    def test_success_status_message(self, qapp, tmp_path):
        """Status bar shows success message after saving."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "success.mol"

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                   return_value=(str(out), None)):
            io.save_3d_as_mol()

        host.statusBar_mock.showMessage.assert_called()
        msg = host.statusBar_mock.showMessage.call_args[0][0]
        assert "saved" in msg.lower() or "3D" in msg

    def test_exception_shows_error_message(self, qapp, tmp_path):
        """RuntimeError during MolToMolBlock shows error on status bar."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "err.mol"

        with patch("moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                   return_value=(str(out), None)), \
             patch("moleditpy.ui.io_logic.Chem.MolToMolBlock",
                   side_effect=RuntimeError("boom")):
            io.save_3d_as_mol()

        host.statusBar_mock.showMessage.assert_called()
        msg = host.statusBar_mock.showMessage.call_args[0][0]
        assert "error" in msg.lower()
