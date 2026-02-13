import os
import json
import pickle
import pytest
from rdkit import Chem
from moleditpy.modules.main_window_project_io import MainWindowProjectIo
from PyQt6.QtWidgets import QMessageBox
from unittest.mock import MagicMock, patch

class DummyProjectIo(MainWindowProjectIo):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.current_file_path = None
        self.has_unsaved_changes = False
        self.current_mol = None
        self._saved_state = None
        self.statusBar_mock = MagicMock()
    
    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self): return self.statusBar_mock
    def update_window_title(self): pass
    def reset_undo_stack(self): pass
    def restore_ui_for_editing(self): pass
    def fit_to_view(self): pass
    def check_unsaved_changes(self): return True # Default to allow in tests

def test_save_project_no_data(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    io.data.atoms = {}
    io.current_mol = None
    io.save_project()
    io.statusBar().showMessage.assert_called_with("Error: Nothing to save.")

def test_save_project_overwrite_json(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    project_file = str(tmp_path / "existing.pmeprj")
    io.current_file_path = project_file
    io.data.atoms = {1: {'symbol': 'C'}}
    
    with patch.object(io, 'create_json_data', return_value={"format": "PME Project"}):
        io.save_project()
        assert os.path.exists(project_file)
        with open(project_file, 'r') as f:
            data = json.load(f)
            assert data["format"] == "PME Project"
        io.statusBar().showMessage.assert_any_call(f"Project saved to {project_file}")

def test_save_project_overwrite_raw(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    raw_file = str(tmp_path / "existing.pmeraw")
    io.current_file_path = raw_file
    io.data.atoms = {1: {'symbol': 'C'}}
    
    with patch.object(io, 'get_current_state', return_value={"atoms": "mock"}):
        io.save_project()
        assert os.path.exists(raw_file)
        with open(raw_file, 'rb') as f:
            data = pickle.load(f)
            assert data["atoms"] == "mock"

def test_save_project_redirect_to_save_as(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    io.current_file_path = "some_molecule.mol"
    io.data.atoms = {1: {'symbol': 'C'}}
    
    with patch.object(io, 'save_project_as') as mock_save_as:
        io.save_project()
        assert mock_save_as.called

def test_load_raw_data_success(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    raw_file = str(tmp_path / "test.pmeraw")
    sample_data = {"atoms": {1: "C"}}
    with open(raw_file, 'wb') as f:
        pickle.dump(sample_data, f)
    
    with patch.object(io, 'set_state_from_data') as mock_set_state:
        io.load_raw_data(raw_file)
        assert mock_set_state.called
        assert io.current_file_path == raw_file

def test_load_json_data_invalid_format(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    json_file = str(tmp_path / "wrong.pmeprj")
    with open(json_file, 'w') as f:
        json.dump({"format": "Unknown"}, f)
    
    with patch('PyQt6.QtWidgets.QMessageBox.warning') as mock_warn:
        io.load_json_data(json_file)
        assert mock_warn.called

def test_open_project_file_dispatch(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    with patch.object(io, 'load_json_data') as mock_json, \
         patch.object(io, 'load_raw_data') as mock_raw:
        
        io.open_project_file("test.pmeprj")
        assert mock_json.called
        
        io.open_project_file("test.pmeraw")
        assert mock_raw.called

def test_save_as_json_trigger(mock_parser_host, tmp_path):
    io = DummyProjectIo(mock_parser_host)
    io.data.atoms = {1: "C"}
    save_path = str(tmp_path / "exported.pmeprj")
    
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.pmeprj")), \
         patch.object(io, 'create_json_data', return_value={"format": "PME Project"}):
        io.save_as_json()
        assert os.path.exists(save_path)

def test_load_raw_data_error_paths(mock_parser_host):
    io = DummyProjectIo(mock_parser_host)
    # File not found
    io.load_raw_data("non_existent.pmeraw")
    io.statusBar().showMessage.assert_called_with("File not found: non_existent.pmeraw")
    
    # Unpickling error
    with patch('builtins.open', MagicMock()):
        with patch('pickle.load', side_effect=pickle.UnpicklingError("Corrupt")):
            io.load_raw_data("dummy.pmeraw")
            io.statusBar().showMessage.assert_called_with("Invalid project file format: Corrupt")
