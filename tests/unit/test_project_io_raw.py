import os
import pickle
import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import QPointF
from moleditpy.ui.project_io import MainWindowProjectIo

class DummyProjectIo(MainWindowProjectIo):
    def __init__(self):
        self.data = MagicMock()
        self.data.atoms = {}
        self.current_mol = None
        self.current_file_path = None
        self.has_unsaved_changes = False
        self._saved_state = None
        self.statusBar_mock = MagicMock()
        self.update_window_title = MagicMock()
        self.reset_undo_stack = MagicMock()
        self.restore_ui_for_editing = MagicMock()
        self.set_state_from_data = MagicMock()
        self.clear_all = MagicMock(return_value=True)
        self.fit_to_view = MagicMock()
        self.get_current_state = MagicMock(return_value={"atoms": "mock"})

    def statusBar(self):
        return self.statusBar_mock

@pytest.fixture
def io():
    return DummyProjectIo()

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
    
    with patch("PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "Project Files (*.pmeraw)")):
        io.save_raw_data()
        
        assert os.path.exists(save_path)
        with open(save_path, "rb") as f:
            data = pickle.load(f)
            assert data == {"atoms": "mock"}
        
        assert io.current_file_path == save_path
        assert io.has_unsaved_changes is False
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
    
    with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName", return_value=(load_path, "Project Files (*.pmeraw)")):
        io.load_raw_data()
        
        io.set_state_from_data.assert_called_with(sample_data)
        assert io.current_file_path == load_path
        assert io.has_unsaved_changes is False
        io.statusBar().showMessage.assert_called_with(f"Project loaded from {load_path}")

def test_load_raw_data_cancel(io):
    """Verify that nothing happens if the user cancels the load dialog."""
    with patch("PyQt6.QtWidgets.QFileDialog.getOpenFileName", return_value=("", "")):
        io.load_raw_data()
        io.set_state_from_data.assert_not_called()

def test_load_raw_data_io_error(io, tmp_path):
    """Verify handling of I/O errors during load."""
    bad_path = str(tmp_path / "corrupt.pmeraw")
    with open(bad_path, "w") as f:
        f.write("not a pickle")
    
    io.load_raw_data(bad_path)
    io.statusBar().showMessage.assert_called()
    # Should show invalid format error
    msg = io.statusBar().showMessage.call_args[0][0]
    assert "Invalid project file format" in msg or "Error loading project file" in msg
