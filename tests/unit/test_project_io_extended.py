import os
import json
import pickle
import pytest
from rdkit import Chem
from moleditpy.ui.io_logic import IOManager
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QMessageBox
from unittest.mock import MagicMock, patch
import copy


class DummyProjectIo(IOManager):
    def __init__(self, host):
        self.host = host
        IOManager.__init__(self, host)

    def __getattr__(self, name):
        """Allow legacy access to host attributes for test convenience."""
        return getattr(self.host, name)

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
    def current_file_path(self): return self.host.init_manager.current_file_path
    @current_file_path.setter
    def current_file_path(self, v): self.host.init_manager.current_file_path = v

    @property
    def has_unsaved_changes(self): return self.host.state_manager.has_unsaved_changes
    @has_unsaved_changes.setter
    def has_unsaved_changes(self, v): self.host.state_manager.has_unsaved_changes = v

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

    with patch.object(io.host.state_manager, "create_json_data", return_value={"format": "PME Project"}):
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

    with patch.object(io.host.state_manager, "get_current_state", return_value={"atoms": "mock"}):
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
        patch.object(io.host.state_manager, "create_json_data", return_value={"format": "PME Project"}),
    ):
        io.save_as_json()
        assert os.path.exists(save_path)


def test_load_raw_data_error_paths(mock_parser_host):
    """Verify error handling during raw data loading (file not found, corrupt)."""
    io = DummyProjectIo(mock_parser_host)
    # File not found
    io.load_raw_data("non_existent.pmeraw")
    io.statusBar().showMessage.assert_called_with("File not found: non_existent.pmeraw")

    # Unpickling error
    with patch("builtins.open", MagicMock()):
        with patch("pickle.load", side_effect=pickle.UnpicklingError("Corrupt")):
            # Should handle exception and log/show error
            io.load_raw_data("dummy_data")
            # Verify error message
            io.statusBar().showMessage.assert_called_with(
                "Invalid project file format: Corrupt"
            )


def test_open_project_file_unsaved_check(mock_parser_host):
    """Verify that 'open project' checks for unsaved changes before proceeding."""
    io = DummyProjectIo(mock_parser_host)
    # Mock check_unsaved_changes to return False (cancel)
    with patch.object(io.host.state_manager, "check_unsaved_changes", return_value=False):
        io.open_project_file()
        # Should return early, not opening dialog
        with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName") as mock_open:
            assert not mock_open.called


def test_save_project_io_error(mock_parser_host, tmp_path):
    """Verify handling of I/O errors during save."""
    io = DummyProjectIo(mock_parser_host)
    io.host.init_manager.current_file_path = str(tmp_path / "readonly.pmeprj")
    io.host.state_manager.data.atoms = {1: "C"}

    with patch("builtins.open", side_effect=IOError("Permission denied")):
        with patch.object(
            io.host.state_manager, "create_json_data", return_value={"format": "PME Project"}
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
        # Should verify warning was shown
        assert mock_info.called
        assert "version 2.0" in mock_info.call_args[0][2]


def test_project_save_load_full_cycle(mock_parser_host, tmp_path):
    """Verify a complete save-load cycle for a project."""
    io = DummyProjectIo(mock_parser_host)
    # Populate some data
    io.host.state_manager.data.add_atom("C", QPointF(10, 20), charge=1)
    # Ensure item mock is present for coordinate extraction if needed
    io.host.state_manager.data.atoms[0]["item"] = MagicMock()
    io.host.state_manager.data.atoms[0]["item"].pos.return_value = QPointF(10, 20)

    project_file = str(tmp_path / "full_cycle.pmeprj")

    # Save
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(project_file, "*.pmeprj"),
    ):
        io.save_as_json()
    assert os.path.exists(project_file)

    # Load
    io.host.state_manager.data.atoms.clear()
    with patch.object(io.host.state_manager, "load_from_json_data") as mock_load_json:
        io.load_json_data(project_file)
        assert mock_load_json.called
        # Verify the data passed to load_from_json_data contains our atom
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

    # 1. With existing path
    io.host.init_manager.current_file_path = str(tmp_path / "subdir" / "my_molecule.pmeprj")

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=("", "")
    ) as mock_save:
        io.save_project_as()

        args, _ = mock_save.call_args
        suggested_path = args[2]
        # Should be in same dir, same name (without ext potentially, or used as is depending on logic)
        # Logic says: os.path.join(dirname, basename_no_ext)
        expected = str(tmp_path / "subdir" / "my_molecule")
        assert suggested_path == expected

    # 2. No existing path
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

    # User types "test" without extension
    save_path_input = str(tmp_path / "test")
    # Expected output
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
    io.host.state_manager.update_window_title = MagicMock()

    save_path = str(tmp_path / "saved.pmeprj")

    # Mock deepcopy to just return the object (or a copy), avoiding issues with MagicMocks
    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "")
        ),
        patch.object(io.host.state_manager, "create_json_data", return_value={}),
        patch.object(io.host.state_manager, "get_current_state", return_value={"test": 1}),
        patch("copy.deepcopy", side_effect=lambda x: x),
    ):
        io.save_project_as()

        assert io.host.state_manager.has_unsaved_changes is False
        assert io.host.init_manager.current_file_path == save_path
        assert io.host.state_manager.update_window_title.called
        assert io.host.state_manager._saved_state is not None
