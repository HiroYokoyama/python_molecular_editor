import pytest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_edit_3d import MainWindowEdit3d
from PyQt6.QtCore import Qt, QPointF
from unittest.mock import MagicMock, patch
import numpy as np


class DummyEdit3d(MainWindowEdit3d):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.view_2d = host.view_2d
        self.view_3d = host.view_3d
        self.settings = host.settings
        self.current_mol = host.current_mol
        self.plotter = host.plotter
        self.selected_atoms_for_measurement = []
        self.measurement_mode = False
        self.atom_positions_3d = []
        self.measurement_actor = None
        self.measurement_labels = []
        self.active_dialogs = []
        self.edit_3d_action = MagicMock()
        self.is_3d_edit_mode = False

    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self):
        return self._host.statusBar()

    def update_measurement_labels_display(self):
        pass

    def update_2d_measurement_labels(self):
        pass

    def add_measurement_label(self, idx, num):
        pass

    def close_all_3d_edit_dialogs(self):
        pass

    def toggle_3d_edit_mode(self, mode):
        pass


def test_calculate_distance_logic(mock_parser_host):
    """Verify calculation of distance between 3D points."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.atom_positions_3d = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
    dist = edit3d.calculate_distance(0, 1)
    assert dist == pytest.approx(1.5)


def test_calculate_angle_logic(mock_parser_host):
    """Verify calculation of angle between three 3D points."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.atom_positions_3d = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    angle = edit3d.calculate_angle(0, 1, 2)
    assert angle == pytest.approx(90.0)


def test_calculate_dihedral_logic(mock_parser_host):
    """Verify calculation of dihedral angle between four 3D points."""
    edit3d = DummyEdit3d(mock_parser_host)
    # Ensure 3D vectors (non-collinear for cross product)
    edit3d.atom_positions_3d = [
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 1.0, 1.0],
    ]
    dihedral = edit3d.calculate_dihedral(0, 1, 2, 3)
    # Dihedral should be non-zero
    assert abs(dihedral) > 0


def test_handle_measurement_atom_selection(mock_parser_host):
    """Verify handling of atom selection for measurements."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.selected_atoms_for_measurement = []
    edit3d.handle_measurement_atom_selection(10)
    assert 10 in edit3d.selected_atoms_for_measurement

    # Second time should return (currently handle_measurement_atom_selection doesn't pop)
    edit3d.handle_measurement_atom_selection(10)
    assert len(edit3d.selected_atoms_for_measurement) == 1


def test_clear_3d_selection(mock_parser_host):
    """Verify clearing of 3D selection highlights."""
    edit3d = DummyEdit3d(mock_parser_host)
    mock_actor = MagicMock()
    mock_actor.name = "atom_10_highlight"
    edit3d.plotter.renderer.actors = {"atom_10_highlight": mock_actor}

    edit3d.clear_3d_selection()
    assert edit3d.plotter.remove_actor.called


def test_toggle_measurement_mode(mock_parser_host):
    """Verify toggling of measurement mode."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.toggle_measurement_mode(True)
    assert edit3d.measurement_mode is True
    assert edit3d.statusBar().showMessage.called

    edit3d.toggle_measurement_mode(False)
    assert edit3d.measurement_mode is False


def test_calculate_and_display_measurements_trigger(mock_parser_host):
    """Verify triggering of measurement calculation and display."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.selected_atoms_for_measurement = [0, 1]
    edit3d.atom_positions_3d = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]

    with patch.object(edit3d, "display_measurement_text") as mock_display:
        edit3d.calculate_and_display_measurements()
        assert mock_display.called
        # It passes a list of strings
        assert "Distance" in mock_display.call_args[0][0][0]
