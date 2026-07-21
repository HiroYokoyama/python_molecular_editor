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

from PyQt6.QtWidgets import QApplication

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
        self.is_restoring_state = False

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

    def set_current_molecule(self, mol):
        self.view_3d_manager.current_mol = mol

    def set_current_file_path(self, path):
        self.init_manager.current_file_path = path

    def set_has_unsaved_changes(self, value):
        self.state_manager.has_unsaved_changes = value

    def set_atom_id_to_rdkit_idx_map(self, mapping):
        self.view_3d_manager.atom_id_to_rdkit_idx_map = mapping

    def update_status_message(self, message, timeout=0):
        if timeout == 0:
            self.statusBar_mock.showMessage(message)
        else:
            self.statusBar_mock.showMessage(message, timeout)


# ===========================================================================
# Tests for prompt_for_charge
# ===========================================================================

# Key insight: QDialog is now a module-level import in io_logic.py, so we
# patch "moleditpy.ui.io_logic.QDialog".  The method checks:
#     if dialog.exec() != QDialog.DialogCode.Accepted
# With our patch, both QDialog (class) and the dialog instance come from
# the MockDialog; so we must return MockDialog.DialogCode.Accepted from exec()
# to pass the "Accepted" branch.  Returning anything else → Rejected path.


from contextlib import contextmanager


@contextmanager
def _make_charge_dialog(line_edit_text="0", accepted=True, skip_btn_class=None):
    """Patch all Qt widgets used by prompt_for_charge."""
    le = MagicMock()
    le.text.return_value = line_edit_text
    btn_patch = (
        patch("moleditpy.ui.io_logic.QPushButton", skip_btn_class)
        if skip_btn_class is not None
        else patch("moleditpy.ui.io_logic.QPushButton")
    )
    with (
        patch("moleditpy.ui.io_logic.QDialog") as MockDlg,
        patch("moleditpy.ui.io_logic.QLineEdit", return_value=le),
        patch("moleditpy.ui.io_logic.QVBoxLayout"),
        patch("moleditpy.ui.io_logic.QHBoxLayout"),
        patch("moleditpy.ui.io_logic.QLabel"),
        patch("moleditpy.ui.io_logic.QDialogButtonBox"),
        btn_patch,
    ):
        dlg = MagicMock()
        code = MockDlg.DialogCode.Accepted if accepted else MockDlg.DialogCode.Rejected
        dlg.exec.return_value = code
        MockDlg.return_value = dlg
        yield MockDlg, le


class TestPromptForCharge:
    def _make_io(self):
        return IOManager(DummyHost())

    @pytest.mark.parametrize(
        "text,expected_charge",
        [
            ("0", 0),
            ("2", 2),
            ("-1", -1),
            ("abc", 0),  # non-numeric falls back to 0
        ],
    )
    def test_accept_returns_charge(self, qapp, text, expected_charge):
        """Accepted dialog returns the entered charge, ok=True, skip=False."""
        io = self._make_io()
        with _make_charge_dialog(line_edit_text=text, accepted=True):
            charge, ok, skip = io.prompt_for_charge()
        assert charge == expected_charge
        assert ok is True
        assert skip is False

    def test_cancel_returns_none_false_false(self, qapp):
        """Cancelled dialog returns charge=None, ok=False, skip=False."""
        io = self._make_io()
        with _make_charge_dialog(accepted=False):
            charge, ok, skip = io.prompt_for_charge()
        assert charge is None
        assert ok is False
        assert skip is False

    def test_skip_chemistry_returns_zero_true_true(self, qapp):
        """Skip button returns charge=0, ok=True, skip=True."""
        io = self._make_io()
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

        with _make_charge_dialog(skip_btn_class=FakeSkipBtn) as (MockDlg, _):
            dlg = MockDlg.return_value
            accepted = MockDlg.DialogCode.Accepted

            def _exec():
                if captured_cb:
                    captured_cb[0]()
                return accepted

            dlg.exec.side_effect = _exec
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

        with (
            patch(
                "moleditpy.ui.io_logic.QFileDialog.getOpenFileName",
                return_value=(str(xyz), None),
            ),
            patch("moleditpy.ui.io_logic.QTimer"),
        ):
            io.load_xyz_for_3d_viewing()

        io.load_xyz_file.assert_called_once_with(str(xyz))
        assert host.view_3d_manager.current_mol is mol

    def test_dialog_cancelled_is_noop(self, qapp):
        """Empty string from dialog → early return, no state change."""
        host = DummyHost()
        io = IOManager(host)
        io.load_xyz_file = MagicMock()

        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getOpenFileName", return_value=("", None)
        ):
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

        host.ui_manager.enter_3d_viewer_mode.assert_called_once()
        host.ui_manager.enable_3d_features.assert_called_once_with(True)

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


class TestShowXYZData:
    """Tests for IOManager.show_xyz_data()."""

    def test_sets_current_mol_and_draws(self, qapp):
        """XYZ text display sets current_mol and enters 3D viewer mode."""
        host = DummyHost()
        io = IOManager(host)
        mol = MagicMock()
        mol.HasProp.return_value = True
        mol.GetIntProp.return_value = 1
        mol.GetNumBonds.return_value = 0
        io.load_xyz_block = MagicMock(return_value=mol)

        result = io.show_xyz_data("- 0.0 0.0 0.0\n", "plugin result")

        assert result is mol
        assert host.view_3d_manager.current_mol is mol
        host.view_3d_manager.draw_molecule_3d.assert_called_once_with(mol)
        host.ui_manager.enter_3d_viewer_mode.assert_called_once()
        host.ui_manager.enable_3d_features.assert_called_once_with(True)
        assert host.is_xyz_derived is True


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

        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getSaveFileName", return_value=("", None)
        ):
            io.save_3d_as_mol()

        assert not out.exists()

    def test_mol_file_written(self, qapp, tmp_path):
        """MOL file is created and non-empty."""
        io, host = self._make_io_with_mol("CCO")
        out = tmp_path / "out.mol"

        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
            return_value=(str(out), None),
        ):
            io.save_3d_as_mol()

        assert out.exists()
        assert out.stat().st_size > 0

    def test_header_line_replaced(self, qapp, tmp_path):
        """Second line is replaced with 'MoleditPy Ver. ... 3D'."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "header.mol"

        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
            return_value=(str(out), None),
        ):
            io.save_3d_as_mol()

        lines = out.read_text(encoding="utf-8").splitlines()
        assert len(lines) > 1
        assert "MoleditPy Ver." in lines[1]
        assert "3D" in lines[1]

    def test_success_status_message(self, qapp, tmp_path):
        """Status bar shows success message after saving."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "success.mol"

        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
            return_value=(str(out), None),
        ):
            io.save_3d_as_mol()

        host.statusBar_mock.showMessage.assert_called()
        msg = host.statusBar_mock.showMessage.call_args[0][0]
        assert "saved" in msg.lower() or "3D" in msg

    def test_exception_shows_error_message(self, qapp, tmp_path):
        """RuntimeError during MolToMolBlock shows error on status bar."""
        io, host = self._make_io_with_mol("C")
        out = tmp_path / "err.mol"

        with (
            patch(
                "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
                return_value=(str(out), None),
            ),
            patch(
                "moleditpy.ui.io_logic.Chem.MolToMolBlock",
                side_effect=RuntimeError("boom"),
            ),
        ):
            io.save_3d_as_mol()

        host.statusBar_mock.showMessage.assert_called()
        msg = host.statusBar_mock.showMessage.call_args[0][0]
        assert "error" in msg.lower()


# ===========================================================================
# Regression tests for code-review fixes
# ===========================================================================


class TestSave3DAsMolExtension:
    """save_3d_as_mol must append .mol when the user omits the extension."""

    def test_missing_extension_is_appended(self, qapp, tmp_path):
        host = DummyHost()
        io = IOManager(host)
        mol = Chem.MolFromSmiles("CCO")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        host.view_3d_manager.current_mol = mol

        out_no_ext = tmp_path / "bare_name"
        with patch(
            "moleditpy.ui.io_logic.QFileDialog.getSaveFileName",
            return_value=(str(out_no_ext), None),
        ):
            io.save_3d_as_mol()

        assert not out_no_ext.exists()
        assert (tmp_path / "bare_name.mol").exists()


class TestPromptForChargeInlineValidation:
    """Invalid input must keep the dialog open with an inline error message,
    not silently fall back to charge 0."""

    def test_invalid_charge_blocks_accept_and_shows_error(self, qapp):
        io = IOManager(DummyHost())
        le = MagicMock()
        le.text.return_value = "abc"
        with (
            patch("moleditpy.ui.io_logic.QDialog") as MockDlg,
            patch("moleditpy.ui.io_logic.QLineEdit", return_value=le),
            patch("moleditpy.ui.io_logic.QVBoxLayout"),
            patch("moleditpy.ui.io_logic.QHBoxLayout"),
            patch("moleditpy.ui.io_logic.QLabel") as MockLabel,
            patch("moleditpy.ui.io_logic.QDialogButtonBox") as MockBox,
            patch("moleditpy.ui.io_logic.QPushButton"),
        ):
            dlg = MockDlg.return_value
            dlg.exec.return_value = MockDlg.DialogCode.Rejected
            io.prompt_for_charge()

            # The OK handler was connected to the inline validator
            validator = MockBox.return_value.accepted.connect.call_args[0][0]

            dlg.accept.reset_mock()
            MockLabel.return_value.setText.reset_mock()

            # Invalid text: dialog stays open, error message set inline
            validator()
            dlg.accept.assert_not_called()
            msgs = [
                str(c[0][0])
                for c in MockLabel.return_value.setText.call_args_list
                if c[0]
            ]
            assert any("Invalid charge" in m for m in msgs)

            # Valid text: dialog accepts
            le.text.return_value = "-2"
            validator()
            dlg.accept.assert_called_once()


# ===========================================================================
# Tests for XYZ dummy/pseudo-atom handling (any XYZ must load cleanly)
# ===========================================================================


class TestXyzDummyAtoms:
    """A dummy/pseudo atom label must map to the RDKit wildcard (*) without
    the load appearing to fail. Regression: probing the periodic table with a
    non-element label used to make RDKit print an 'Element 'Xx' not found'
    violation to stderr, which reads as a load failure."""

    def _io(self):
        host = DummyHost()
        host.init_manager.settings = {
            "skip_chemistry_checks": False,
            "always_ask_charge": False,
        }
        return IOManager(host)

    @pytest.mark.parametrize("label", ["Xx", "XX", "X1", "Bq", "DA", "Q", "LP"])
    def test_dummy_label_normalizes_to_wildcard(self, qapp, label):
        io = self._io()
        symbol, is_dummy = io._normalize_xyz_symbol(label)
        assert symbol == "*"
        assert is_dummy is True

    @pytest.mark.parametrize(
        "label,expected", [("C", "C"), ("fe", "Fe"), ("og", "Og"), ("H", "H")]
    )
    def test_real_element_preserved(self, qapp, label, expected):
        io = self._io()
        symbol, is_dummy = io._normalize_xyz_symbol(label)
        assert symbol == expected
        assert is_dummy is False

    def test_dummy_atom_xyz_loads_to_mol(self, qapp):
        io = self._io()
        xyz = "3\n\nXx 0 0 0\nC 0 0 1.0\nO 0 1.2 0\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol is not None
        assert mol.GetNumAtoms() == 3
        assert mol.GetAtomWithIdx(0).GetSymbol() == "*"
        # binary round-trips (used for undo/redo and .pmeprj persistence)
        assert Chem.Mol(mol.ToBinary()).GetNumAtoms() == 3

    def test_normalize_does_not_emit_rdkit_violation(self, qapp, tmp_path):
        """Normalizing a non-element label must not print a C++ periodic-table
        violation to stderr — that stderr noise is the user-visible symptom."""
        io = self._io()
        # Capture at the file-descriptor level: RDKit's C++ logger writes to the
        # underlying stderr fd, not Python's sys.stderr.
        err_path = tmp_path / "stderr.txt"
        saved_fd = os.dup(2)
        with open(err_path, "w") as fh:
            os.dup2(fh.fileno(), 2)
            try:
                for label in ["Xx", "XX", "X1", "Bq", "Zz", "Qq"]:
                    io._normalize_xyz_symbol(label)
            finally:
                os.dup2(saved_fd, 2)
                os.close(saved_fd)
        captured = err_path.read_text()
        assert "not found" not in captured, captured
        assert "Violation" not in captured, captured


class TestXyzLoadRobustness:
    """Any structurally-parseable XYZ must yield a molecule (never None),
    even when bond perception fails or the atom-count header is wrong."""

    def _io(self, **settings):
        host = DummyHost()
        base = {"skip_chemistry_checks": False, "always_ask_charge": False}
        base.update(settings)
        host.init_manager.settings = base
        return IOManager(host)

    def test_bond_estimation_failure_still_loads_atoms(self, qapp):
        io = self._io()
        io.estimate_bonds_from_distances = MagicMock(
            side_effect=RuntimeError("bond perception blew up")
        )
        xyz = "3\n\nXx 0 0 0\nC 0 0 1.0\nO 0 1.2 0\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol is not None
        assert mol.GetNumAtoms() == 3
        assert mol.GetNumBonds() == 0

    def test_truncated_atom_count_loads_available_rows(self, qapp):
        io = self._io()
        # header claims 5 atoms, only 2 rows present
        xyz = "5\n\nXx 0 0 0\nC 0 0 1.0\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol is not None
        assert mol.GetNumAtoms() == 2

    def test_skip_checks_with_bond_failure_loads(self, qapp):
        io = self._io(skip_chemistry_checks=True)
        io.estimate_bonds_from_distances = MagicMock(
            side_effect=ValueError("bad geometry")
        )
        xyz = "2\n\nC 0 0 0\nO 0 0 1.2\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol is not None
        assert mol.GetNumAtoms() == 2


class TestXyzGhostAtomColumns:
    """Ghost-atom rows may carry a separator/label column, e.g. 'XX : x y z'.
    Coordinates must be read from the numeric columns, not blindly from
    parts[1..3] (which would hit the ':' and fail the whole file)."""

    def _io(self):
        host = DummyHost()
        host.init_manager.settings = {
            "skip_chemistry_checks": False,
            "always_ask_charge": False,
        }
        return IOManager(host)

    def test_extract_coords_skips_separator_column(self, qapp):
        io = self._io()
        assert io._extract_xyz_coords([":", "1", "2", "3"]) == (1.0, 2.0, 3.0)
        assert io._extract_xyz_coords(["1", "2", "3"]) == (1.0, 2.0, 3.0)
        # trailing extra column (e.g. a charge) is ignored
        assert io._extract_xyz_coords(["1", "2", "3", "0.5"]) == (1.0, 2.0, 3.0)
        # too few numbers -> None (caller raises a clear per-line error)
        assert io._extract_xyz_coords(["1", "2"]) is None

    def test_mixed_standard_and_ghost_columns_load(self, qapp):
        io = self._io()
        xyz = (
            "5\n\n"
            "C   -53.9 38.9 -0.9\n"
            "H   -51.6 42.7  1.1\n"
            "XX : -80.5 39.9 -3.0\n"
            "XX : -19.6 38.8 -3.2\n"
            "XX : -51.7 36.8 -1.6\n"
        )
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol is not None
        assert mol.GetNumAtoms() == 5
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(5)]
        assert symbols == ["C", "H", "*", "*", "*"]
        # the ghost atoms kept their real coordinates (not shifted by the ':')
        pos2 = mol.GetConformer().GetAtomPosition(2)
        assert round(pos2.x, 1) == -80.5
        assert round(pos2.z, 1) == -3.0


class TestXyzNonStandardWarning:
    """A non-standard column layout loads successfully but must be flagged: the
    parser records a count property and load_xyz_for_3d_viewing surfaces it on
    the status bar. Standard files stay silent."""

    def _io(self):
        host = DummyHost()
        host.init_manager.settings = {
            "skip_chemistry_checks": False,
            "always_ask_charge": False,
        }
        return IOManager(host), host

    def test_is_numeric_token(self, qapp):
        io, _ = self._io()
        assert io._is_numeric_token("1.5") is True
        assert io._is_numeric_token("-3.0e2") is True
        assert io._is_numeric_token(":") is False
        assert io._is_numeric_token("XX") is False

    def test_parser_records_nonstandard_row_count(self, qapp):
        io, _ = self._io()
        xyz = "3\n\nC 0 0 0\nXX : 1 0 0\nXX : 0 1 0\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert mol.HasProp("_xyz_nonstandard_rows")
        assert mol.GetIntProp("_xyz_nonstandard_rows") == 2

    def test_parser_no_property_for_standard_file(self, qapp):
        io, _ = self._io()
        xyz = "2\n\nC 0 0 0\nO 0 0 1.2\n"
        mol = io._mol_from_xyz_lines(xyz.splitlines())
        assert not mol.HasProp("_xyz_nonstandard_rows")

    def test_status_bar_warns_on_nonstandard(self, qapp, tmp_path):
        io, host = self._io()
        xyz = tmp_path / "ghost.xyz"
        xyz.write_text("3\n\nC 0 0 0\nXX : 1 0 0\nXX : 0 1 0\n")
        with patch("moleditpy.ui.io_logic.QTimer"), patch(
            "moleditpy.ui.io_logic.QMessageBox"
        ):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))
        msgs = [c[0][0] for c in host.statusBar_mock.showMessage.call_args_list if c[0]]
        assert any("non-standard columns in 2 row" in m for m in msgs), msgs

    def test_status_bar_silent_on_standard(self, qapp, tmp_path):
        io, host = self._io()
        xyz = tmp_path / "std.xyz"
        xyz.write_text("2\n\nC 0 0 0\nO 0 0 1.2\n")
        with patch("moleditpy.ui.io_logic.QTimer"), patch(
            "moleditpy.ui.io_logic.QMessageBox"
        ):
            io.load_xyz_for_3d_viewing(file_path=str(xyz))
        msgs = [c[0][0] for c in host.statusBar_mock.showMessage.call_args_list if c[0]]
        assert msgs
        assert all("non-standard" not in m for m in msgs), msgs
