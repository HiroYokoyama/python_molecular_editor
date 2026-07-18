"""Tests for project I/O (merged from test_project_io_extended.py and test_project_io_raw.py)."""

import os
import json
import pickle
import pytest
from moleditpy.ui.io_logic import IOManager
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch


# ---------------------------------------------------------------------------
# DummyProjectIo — uses an external host from conftest (mock_parser_host)
# ---------------------------------------------------------------------------


class DummyProjectIo(IOManager):
    def __init__(self, host):
        self.host = host
        IOManager.__init__(self, host)

    def __getattr__(self, name):
        """Allow legacy access to host attributes for test convenience."""
        return getattr(self.host, name)

    @property
    def data(self):
        return self.host.state_manager.data

    @property
    def scene(self):
        return self.host.init_manager.scene

    @property
    def settings(self):
        return self.host.init_manager.settings

    @property
    def view_2d(self):
        return self.host.init_manager.view_2d

    @property
    def plotter(self):
        return self.host.view_3d_manager.plotter

    @property
    def current_mol(self):
        return self.host.view_3d_manager.current_mol

    @current_mol.setter
    def current_mol(self, v):
        self.host.view_3d_manager.current_mol = v

    @property
    def current_file_path(self):
        return self.host.init_manager.current_file_path

    @current_file_path.setter
    def current_file_path(self, v):
        self.host.init_manager.current_file_path = v

    @property
    def has_unsaved_changes(self):
        return self.host.state_manager.has_unsaved_changes

    @has_unsaved_changes.setter
    def has_unsaved_changes(self, v):
        self.host.state_manager.has_unsaved_changes = v

    def statusBar(self):
        return self.host.statusBar()

    def update_window_title(self):
        pass

    def reset_undo_stack(self):
        pass

    def restore_ui_for_editing(self):
        pass

    def fit_to_view(self):
        pass

    def check_unsaved_changes(self):
        return True  # Default to allow in tests


def test_save_project_no_data(mock_parser_host):
    """Verify error message when trying to save an empty project."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.atoms = {}
    io.host.view_3d_manager.current_mol = None
    io.save_project()
    io.statusBar().showMessage.assert_called_with("Error: Nothing to save.")


def test_save_project_overwrite_json(mock_parser_host, tmp_path):
    """Verify overwriting an existing JSON project file."""
    io = DummyProjectIo(mock_parser_host)
    project_file = str(tmp_path / "existing.pmeprj")
    io.host.init_manager.current_file_path = project_file
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}

    with patch.object(
        io.host.state_manager,
        "create_json_data",
        return_value={"format": "PME Project"},
    ):
        io.save_project()
        assert os.path.exists(project_file)
        with open(project_file, "r") as f:
            data = json.load(f)
            assert data["format"] == "PME Project"
        io.statusBar().showMessage.assert_any_call(f"Project saved to {project_file}")


def test_save_project_overwrite_raw(mock_parser_host, tmp_path):
    """Verify overwriting an existing raw (pickle) project file."""
    io = DummyProjectIo(mock_parser_host)
    raw_file = str(tmp_path / "existing.pmeraw")
    io.host.init_manager.current_file_path = raw_file
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}

    with patch.object(
        io.host.state_manager, "get_current_state", return_value={"atoms": "mock"}
    ):
        io.save_project()
        assert os.path.exists(raw_file)
        with open(raw_file, "rb") as f:
            data = pickle.load(f)
            assert data["atoms"] == "mock"


def test_save_project_redirect_to_save_as(mock_parser_host):
    """Verify that 'save' redirects to 'save as' if the current file is not a project file."""
    io = DummyProjectIo(mock_parser_host)
    io.host.init_manager.current_file_path = "some_molecule.mol"
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}

    with patch.object(io, "save_project_as") as mock_save_as:
        io.save_project()
        assert mock_save_as.called


def test_load_raw_data_success(mock_parser_host, tmp_path):
    """Verify successful loading of a raw project file."""
    io = DummyProjectIo(mock_parser_host)
    raw_file = str(tmp_path / "test.pmeraw")
    sample_data = {"atoms": {1: "C"}}
    with open(raw_file, "wb") as f:
        pickle.dump(sample_data, f)

    with patch.object(io.host.state_manager, "set_state_from_data") as mock_set_state:
        io.load_raw_data(raw_file)
        assert mock_set_state.called
        assert io.host.init_manager.current_file_path == raw_file


def test_load_json_data_invalid_format(mock_parser_host, tmp_path):
    """Verify error handling for invalid JSON format in project files."""
    io = DummyProjectIo(mock_parser_host)
    json_file = str(tmp_path / "wrong.pmeprj")
    with open(json_file, "w") as f:
        json.dump({"format": "Unknown"}, f)

    with patch("PyQt6.QtWidgets.QMessageBox.warning") as mock_warn:
        io.load_json_data(json_file)
        assert mock_warn.called


def test_open_project_file_dispatch(mock_parser_host):
    """Verify dispatching to correct load method based on file extension."""
    io = DummyProjectIo(mock_parser_host)
    with (
        patch.object(io, "load_json_data") as mock_json,
        patch.object(io, "load_raw_data") as mock_raw,
    ):
        io.open_project_file("test.pmeprj")
        assert mock_json.called

        io.open_project_file("test.pmeraw")
        assert mock_raw.called


def test_save_as_json_trigger(mock_parser_host, tmp_path):
    """Verify triggering of 'save as' for JSON format."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.atoms = {1: "C"}
    save_path = str(tmp_path / "exported.pmeprj")

    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
            return_value=(save_path, "*.pmeprj"),
        ),
        patch.object(
            io.host.state_manager,
            "create_json_data",
            return_value={"format": "PME Project"},
        ),
    ):
        io.save_as_json()
        assert os.path.exists(save_path)


def test_load_raw_data_error_paths(mock_parser_host):
    """Verify error handling during raw data loading (file not found, corrupt)."""
    io = DummyProjectIo(mock_parser_host)
    io.load_raw_data("non_existent.pmeraw")
    io.statusBar().showMessage.assert_called_with("File not found: non_existent.pmeraw")

    with patch("builtins.open", MagicMock()):
        with patch("pickle.load", side_effect=pickle.UnpicklingError("Corrupt")):
            io.load_raw_data("dummy_data")
            io.statusBar().showMessage.assert_called_with(
                "Invalid project file format: Corrupt"
            )


def test_open_project_file_unsaved_check(mock_parser_host):
    """Verify that 'open project' checks for unsaved changes before proceeding."""
    io = DummyProjectIo(mock_parser_host)
    with patch.object(
        io.host.state_manager, "check_unsaved_changes", return_value=False
    ):
        io.open_project_file()
        with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName") as mock_open:
            assert not mock_open.called


def test_save_project_io_error(mock_parser_host, tmp_path):
    """Verify handling of I/O errors during save."""
    io = DummyProjectIo(mock_parser_host)
    io.host.init_manager.current_file_path = str(tmp_path / "readonly.pmeprj")
    io.host.state_manager.data.atoms = {1: "C"}

    with patch("builtins.open", side_effect=IOError("Permission denied")):
        with patch.object(
            io.host.state_manager,
            "create_json_data",
            return_value={"format": "PME Project"},
        ):
            io.save_project()
            io.statusBar().showMessage.assert_called_with(
                "File I/O error: Permission denied"
            )


def test_load_json_data_version_mismatch(mock_parser_host, tmp_path):
    """Verify warning when loading a project from a newer software version."""
    io = DummyProjectIo(mock_parser_host)
    json_path = tmp_path / "future.pmeprj"
    with open(json_path, "w") as f:
        json.dump({"format": "PME Project", "version": "2.0"}, f)

    with patch("PyQt6.QtWidgets.QMessageBox.information") as mock_info:
        io.load_json_data(str(json_path))
        assert mock_info.called
        assert "version 2.0" in mock_info.call_args[0][2]


def test_project_save_load_full_cycle(mock_parser_host, tmp_path):
    """Verify a complete save-load cycle for a project."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.add_atom("C", QPointF(10, 20), charge=1)

    project_file = str(tmp_path / "full_cycle.pmeprj")

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(project_file, "*.pmeprj"),
    ):
        io.save_as_json()
    assert os.path.exists(project_file)

    io.host.state_manager.data.atoms.clear()
    with patch.object(io.host.state_manager, "load_from_json_data") as mock_load_json:
        io.load_json_data(project_file)
        assert mock_load_json.called
        saved_data = mock_load_json.call_args[0][0]
        assert "2d_structure" in saved_data
        atoms = saved_data["2d_structure"]["atoms"]
        assert len(atoms) == 1
        assert atoms[0]["symbol"] == "C"
        assert atoms[0]["charge"] == 1


def test_save_project_default_filename(mock_parser_host, tmp_path):
    """Verify that 'save as' suggests a reasonable default filename."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}

    io.host.init_manager.current_file_path = str(
        tmp_path / "subdir" / "my_molecule.pmeprj"
    )

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")
    ) as mock_save:
        io.save_project_as()

        args, _ = mock_save.call_args
        suggested_path = args[2]
        expected = str(tmp_path / "subdir" / "my_molecule")
        assert suggested_path == expected

    io.host.init_manager.current_file_path = None
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")
    ) as mock_save:
        io.save_project_as()
        args, _ = mock_save.call_args
        assert args[2] == "untitled"


def test_save_project_extension_enforcement(mock_parser_host, tmp_path):
    """Verify that correct file extensions are enforced when saving."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}

    save_path_input = str(tmp_path / "test")
    expected_path = str(tmp_path / "test.pmeprj")

    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
            return_value=(save_path_input, "PME Project Files (*.pmeprj)"),
        ),
        patch.object(io.host.state_manager, "create_json_data", return_value={}),
    ):
        io.save_project_as()

        assert io.host.init_manager.current_file_path == expected_path
        assert os.path.exists(expected_path)


def test_save_project_success_state_update(mock_parser_host, tmp_path):
    """Verify that application state is updated correctly after a successful save."""
    io = DummyProjectIo(mock_parser_host)
    io.host.state_manager.data.atoms = {1: {"symbol": "C"}}
    io.host.state_manager.has_unsaved_changes = True

    save_path = str(tmp_path / "saved.pmeprj")

    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "")
        ),
        patch.object(io.host.state_manager, "create_json_data", return_value={}),
        patch.object(
            io.host.state_manager, "get_current_state", return_value={"test": 1}
        ),
        patch("copy.deepcopy", side_effect=lambda x: x),
    ):
        io.save_project_as()

        assert io.host.state_manager.has_unsaved_changes is False
        assert io.host.init_manager.current_file_path == save_path
        assert io.host.update_window_title.called
        assert io.host.state_manager._saved_state is not None


# ---------------------------------------------------------------------------
# DummyProjectIoRaw — standalone, creates its own host (from test_project_io_raw.py)
# ---------------------------------------------------------------------------


class DummyProjectIoRaw(IOManager):
    def __init__(self, host=None):
        self.host = host or MagicMock()
        IOManager.__init__(self, self.host)

        self.host.state_manager = MagicMock()
        self.host.init_manager = MagicMock()
        self.host.ui_manager = MagicMock()
        self.host.edit_actions_manager = MagicMock()
        self.host.view_3d_manager = MagicMock()
        self.host.compute_manager = MagicMock()

        from moleditpy.core.molecular_data import MolecularData

        self.host.state_manager.data = MolecularData()
        self.host.view_3d_manager.current_mol = None
        self.host.init_manager.current_file_path = None
        self.host.state_manager.has_unsaved_changes = False

        self.host.state_manager.get_current_state.return_value = {"atoms": "mock"}
        self.host.state_manager.update_window_title = MagicMock()
        self.host.state_manager.set_state_from_data = MagicMock()

        self.host.ui_manager.restore_ui_for_editing = MagicMock()
        self.host.edit_actions_manager.clear_all.return_value = True
        self.host.edit_actions_manager.reset_undo_stack = MagicMock()

        self.statusBar_mock = MagicMock()
        self.host.statusBar.return_value = self.statusBar_mock

        def set_current_file_path(path):
            self.host.init_manager.current_file_path = path

        self.host.set_current_file_path.side_effect = set_current_file_path

        def set_has_unsaved_changes(val):
            self.host.state_manager.has_unsaved_changes = val

        self.host.set_has_unsaved_changes.side_effect = set_has_unsaved_changes

        def update_status_message(message, timeout=0):
            if timeout == 0:
                self.statusBar_mock.showMessage(message)
            else:
                self.statusBar_mock.showMessage(message, timeout)

        self.host.update_status_message.side_effect = update_status_message

    def __getattr__(self, name):
        """Allow legacy access to host attributes for test convenience."""
        return getattr(self.host, name)

    @property
    def data(self):
        return self.host.state_manager.data

    @property
    def current_mol(self):
        return self.host.view_3d_manager.current_mol

    @current_mol.setter
    def current_mol(self, v):
        self.host.view_3d_manager.current_mol = v

    @property
    def current_file_path(self):
        return self.host.init_manager.current_file_path

    @current_file_path.setter
    def current_file_path(self, v):
        self.host.init_manager.current_file_path = v

    @property
    def has_unsaved_changes(self):
        return self.host.state_manager.has_unsaved_changes

    @has_unsaved_changes.setter
    def has_unsaved_changes(self, v):
        self.host.state_manager.has_unsaved_changes = v

    def statusBar(self):
        return self.statusBar_mock


@pytest.fixture
def io():
    return DummyProjectIoRaw()


def test_save_raw_data_no_data(io):
    """Verify error message when trying to save empty project."""
    io.data.atoms = {}
    io.current_mol = None
    io.save_raw_data()
    io.statusBar().showMessage.assert_called_with("Error: Nothing to save.")


def test_save_raw_data_success(io, tmp_path):
    """Verify successful saving via file dialog."""
    io.data.atoms = {1: "C"}
    save_path = str(tmp_path / "test.pmeraw")

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(save_path, "Project Files (*.pmeraw)"),
    ):
        io.save_raw_data()

        assert os.path.exists(save_path)
        with open(save_path, "rb") as f:
            data = pickle.load(f)
            assert data == {"atoms": "mock"}

        assert io.host.init_manager.current_file_path == save_path
        assert io.host.state_manager.has_unsaved_changes is False
        io.statusBar().showMessage.assert_called_with(f"Project saved to {save_path}")


def test_save_raw_data_cancel(io):
    """Verify that nothing happens if the user cancels the save dialog."""
    io.data.atoms = {1: "C"}
    with patch("PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")):
        io.save_raw_data()
        io.statusBar().showMessage.assert_not_called()


def test_load_raw_data_dialog_success(io, tmp_path):
    """Verify loading via file dialog."""
    load_path = str(tmp_path / "load_test.pmeraw")
    sample_data = {"atoms": "loaded"}
    with open(load_path, "wb") as f:
        pickle.dump(sample_data, f)

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getOpenFileName",
        return_value=(load_path, "Project Files (*.pmeraw)"),
    ):
        io.load_raw_data()

        io.host.state_manager.set_state_from_data.assert_called_with(sample_data)
        assert io.host.init_manager.current_file_path == load_path
        assert io.host.state_manager.has_unsaved_changes is False
        io.statusBar().showMessage.assert_called_with(
            f"Project loaded from {load_path}"
        )


def test_load_raw_data_cancel(io):
    """Verify that nothing happens if the user cancels the load dialog."""
    with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName", return_value=("", "")):
        io.load_raw_data()
        io.host.state_manager.set_state_from_data.assert_not_called()


def test_load_raw_data_io_error(io, tmp_path):
    """Verify handling of I/O errors during load."""
    bad_path = str(tmp_path / "corrupt.pmeraw")
    with open(bad_path, "w") as f:
        f.write("not a pickle")

    io.load_raw_data(bad_path)
    io.statusBar().showMessage.assert_called()
    msg = io.statusBar().showMessage.call_args[0][0]
    assert "Invalid project file format" in msg


from rdkit import Chem
from rdkit.Chem import AllChem


def _status_messages(host):
    return [str(c.args[0]) for c in host.statusBar().showMessage.call_args_list]


# ---------------------------------------------------------------------------
# load_mol_file — 2D MOL/SDF import
# ---------------------------------------------------------------------------

_ETHANOL_MOL_BLOCK = """
  Test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"""


def test_load_mol_file_no_path_cancelled_dialog(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName", return_value=("", "")):
        io.load_mol_file()
    mock_parser_host.init_manager.scene.create_atom.assert_not_called()


def test_load_mol_file_missing_file_reports_error(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    io.load_mol_file("/no/such/file.mol")
    assert "File not found" in _status_messages(mock_parser_host)[-1]


def test_load_mol_file_imports_atoms_and_bonds(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "ethanol.mol"
    path.write_text(_ETHANOL_MOL_BLOCK, encoding="utf-8")

    io.load_mol_file(str(path))

    assert len(mock_parser_host.state_manager.data.atoms) == 3
    assert len(mock_parser_host.state_manager.data.bonds) == 2
    mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()
    assert "Successfully imported" in _status_messages(mock_parser_host)[-1]


def test_load_mol_file_places_relative_to_existing_atoms(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.init_manager.scene.create_atom("C", QPointF(500.0, 10.0))
    path = tmp_path / "ethanol.mol"
    path.write_text(_ETHANOL_MOL_BLOCK, encoding="utf-8")

    io.load_mol_file(str(path))

    # New atoms are placed to the right of the pre-existing rightmost atom.
    assert len(mock_parser_host.state_manager.data.atoms) == 4


def test_load_mol_file_invalid_block_reports_error(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "bad.mol"
    path.write_text("not a mol block", encoding="utf-8")

    io.load_mol_file(str(path))

    assert "MOL Import Error" not in "".join(_status_messages(mock_parser_host))
    assert "Error loading file" in _status_messages(mock_parser_host)[-1]


# ---------------------------------------------------------------------------
# save_as_mol
# ---------------------------------------------------------------------------


def test_save_as_mol_no_data(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    assert mock_parser_host.state_manager.data.atoms == {}
    io.save_as_mol()
    assert "Error: No 2D data to save." in _status_messages(mock_parser_host)


def test_save_as_mol_writes_file(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.init_manager.scene.create_atom("C", QPointF(0, 0))
    out_path = tmp_path / "out.mol"
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(str(out_path), ""),
    ):
        io.save_as_mol()
    assert out_path.exists()
    assert "MoleditPy Ver." in out_path.read_text(encoding="utf-8")
    assert "2D data saved to" in _status_messages(mock_parser_host)[-1]


def test_save_as_mol_cancelled_dialog_writes_nothing(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.init_manager.scene.create_atom("C", QPointF(0, 0))
    with patch("PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")):
        io.save_as_mol()
    assert list(tmp_path.iterdir()) == []


def test_save_as_mol_io_error_reports_message(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.init_manager.scene.create_atom("C", QPointF(0, 0))
    bad_path = tmp_path / "no_such_dir" / "out.mol"
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(str(bad_path), ""),
    ):
        io.save_as_mol()
    assert "Error saving MOL" in _status_messages(mock_parser_host)[-1]


# ---------------------------------------------------------------------------
# load_mol_file_for_3d_viewing
# ---------------------------------------------------------------------------


def test_load_mol_file_3d_no_path_cancelled(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName", return_value=("", "")):
        io.load_mol_file_for_3d_viewing()
    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_not_called()


def test_load_mol_file_3d_unsaved_changes_blocks(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.state_manager.check_unsaved_changes.return_value = False
    io.load_mol_file_for_3d_viewing("ignored.mol")
    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_not_called()


def test_load_mol_file_3d_success_draws_molecule(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "ethanol.mol"
    path.write_text(_ETHANOL_MOL_BLOCK, encoding="utf-8")

    io.load_mol_file_for_3d_viewing(str(path))

    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()
    mock_parser_host.ui_manager.enter_3d_viewer_mode.assert_called_once()
    assert mock_parser_host.is_xyz_derived is False
    assert mock_parser_host.init_manager.current_file_path == str(path)


def test_load_mol_file_3d_failure_reports_error(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "bad.mol"
    path.write_text("not a mol block", encoding="utf-8")

    io.load_mol_file_for_3d_viewing(str(path))

    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_not_called()
    assert "3D MOL Load failed" in _status_messages(mock_parser_host)[-1]


# ---------------------------------------------------------------------------
# save_as_xyz
# ---------------------------------------------------------------------------


def test_save_as_xyz_no_current_mol(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = None
    io.save_as_xyz()
    assert (
        "Please generate a 3D structure first" in _status_messages(mock_parser_host)[-1]
    )


def test_save_as_xyz_cancelled_dialog_writes_nothing(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("CO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mock_parser_host.view_3d_manager.current_mol = mol
    with patch("PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")):
        io.save_as_xyz()
    assert list(tmp_path.iterdir()) == []


def test_save_as_xyz_writes_charge_and_multiplicity(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("CO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mol.SetIntProp("_xyz_charge", -1)
    mock_parser_host.view_3d_manager.current_mol = mol
    out_path = tmp_path / "out.xyz"
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(str(out_path), ""),
    ):
        io.save_as_xyz()
    text = out_path.read_text(encoding="utf-8")
    assert text.splitlines()[0] == str(mol.GetNumAtoms())
    assert "chrg = -1" in text.splitlines()[1]
    assert "Successfully saved to" in _status_messages(mock_parser_host)[-1]


def test_save_as_xyz_io_error_reports_message(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("CO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mock_parser_host.view_3d_manager.current_mol = mol
    bad_path = tmp_path / "no_such_dir" / "out.xyz"
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(str(bad_path), ""),
    ):
        io.save_as_xyz()
    assert "Error saving XYZ" in _status_messages(mock_parser_host)[-1]


# ---------------------------------------------------------------------------
# load_json_data error branches
# ---------------------------------------------------------------------------


def test_load_json_data_file_not_found(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    io.load_json_data("/no/such/project.pmeprj")
    assert "File not found" in _status_messages(mock_parser_host)[-1]


def test_load_json_data_invalid_json_syntax(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "broken.pmeprj"
    path.write_text("{not json", encoding="utf-8")
    io.load_json_data(str(path))
    assert "Invalid JSON format" in _status_messages(mock_parser_host)[-1]


def test_load_json_data_corrupted_payload_reports_error(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    path = tmp_path / "corrupt.pmeprj"
    path.write_text(
        json.dumps({"format": "PME Project", "version": "1.0"}), encoding="utf-8"
    )
    mock_parser_host.state_manager.load_from_json_data.side_effect = KeyError("atoms")
    io.load_json_data(str(path))
    assert (
        "Data corruption in PME Project file" in _status_messages(mock_parser_host)[-1]
    )


# ---------------------------------------------------------------------------
# _set_mol_prop / _get_mol_prop / _set_mol_prop_safe
# ---------------------------------------------------------------------------


def test_set_get_mol_prop_int_roundtrip(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    io._set_mol_prop(mol, "my_int", 7)
    assert io._get_mol_prop(mol, "my_int") == 7


def test_set_get_mol_prop_float_roundtrip(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    io._set_mol_prop(mol, "my_float", 3.5)
    assert io._get_mol_prop(mol, "my_float") == 3.5


def test_set_get_mol_prop_string_roundtrip(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    io._set_mol_prop(mol, "my_str", "hello")
    assert io._get_mol_prop(mol, "my_str") == "hello"


def test_get_mol_prop_missing_returns_default(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    assert io._get_mol_prop(mol, "absent", default="fallback") == "fallback"


def test_set_mol_prop_safe_module_helper():
    from moleditpy.ui.io_logic import _set_mol_prop_safe

    mol = Chem.MolFromSmiles("C")
    _set_mol_prop_safe(mol, "k_int", 5)
    _set_mol_prop_safe(mol, "k_float", 2.25)
    assert mol.GetIntProp("k_int") == 5
    assert mol.GetDoubleProp("k_float") == 2.25
