import pytest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.edit_3d_logic import Edit3DManager
from PyQt6.QtCore import Qt, QPointF, QRectF
from PyQt6.QtWidgets import QGraphicsTextItem
from unittest.mock import MagicMock, patch
import numpy as np


class DummyEdit3d(Edit3DManager):
    def __init__(self, host):
        super().__init__(host)
        self._host = host
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene
        self.view_2d = host.init_manager.view_2d
        self.view_3d = host.view_3d_manager.view_3d
        self.settings = host.init_manager.settings
        self.current_mol = host.view_3d_manager.current_mol
        self.plotter = host.view_3d_manager.plotter
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
    positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
    edit3d.atom_positions_3d = positions
    mock_parser_host.view_3d_manager.atom_positions_3d = positions
    dist = edit3d.calculate_distance(0, 1)
    assert dist == pytest.approx(1.5)


def test_calculate_angle_logic(mock_parser_host):
    """Verify calculation of angle between three 3D points."""
    edit3d = DummyEdit3d(mock_parser_host)
    positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    edit3d.atom_positions_3d = positions
    mock_parser_host.view_3d_manager.atom_positions_3d = positions
    angle = edit3d.calculate_angle(0, 1, 2)
    assert angle == pytest.approx(90.0)


def test_calculate_dihedral_logic(mock_parser_host):
    """Verify calculation of dihedral angle between four 3D points."""
    edit3d = DummyEdit3d(mock_parser_host)
    # Ensure 3D vectors (non-collinear for cross product)
    positions = [
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 1.0, 1.0],
    ]
    edit3d.atom_positions_3d = positions
    mock_parser_host.view_3d_manager.atom_positions_3d = positions
    dihedral = edit3d.calculate_dihedral(0, 1, 2, 3)
    # For these coordinates, the dihedral is exactly 45.0°
    assert abs(dihedral - 45.0) < 0.01, f"Expected 45.0°, got {dihedral}°"


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


def test_add_2d_measurement_label_adds_item(mock_parser_host):
    """Verify add_2d_measurement_label creates a label and registers it."""
    edit3d = DummyEdit3d(mock_parser_host)

    atom_item = MagicMock()
    atom_item.pos.return_value = QPointF(10.0, 20.0)
    atom_item.boundingRect.return_value = QRectF(0, 0, 16.0, 16.0)

    edit3d.add_2d_measurement_label(atom_item, "1.234 Å")

    # Item must be registered
    assert len(edit3d.measurement_label_items_2d) == 1
    label = edit3d.measurement_label_items_2d[0]
    assert isinstance(label, QGraphicsTextItem)

    # Correct text
    assert label.toPlainText() == "1.234 Å"

    # Z-value must be 2000 for top-most display
    assert label.zValue() == 2000

    # Must have been added to the scene
    mock_parser_host.init_manager.scene.addItem.assert_called_once_with(label)

    # Position offset: x = 10 + 16/4 + 2 = 16, y = 20 - 16/4 - 8 = 8
    pos = label.pos()
    assert pos.x() == pytest.approx(16.0)
    assert pos.y() == pytest.approx(8.0)


def test_add_2d_measurement_label_accumulates(mock_parser_host):
    """Verify multiple labels are each appended to measurement_label_items_2d."""
    edit3d = DummyEdit3d(mock_parser_host)

    for i in range(3):
        atom_item = MagicMock()
        atom_item.pos.return_value = QPointF(float(i * 10), 0.0)
        atom_item.boundingRect.return_value = QRectF(0, 0, 8.0, 8.0)
        edit3d.add_2d_measurement_label(atom_item, f"label_{i}")

    assert len(edit3d.measurement_label_items_2d) == 3
    assert mock_parser_host.init_manager.scene.addItem.call_count == 3


def test_clear_2d_measurement_labels_removes_items(mock_parser_host):
    """Verify clear_2d_measurement_labels removes all items from the scene."""
    edit3d = DummyEdit3d(mock_parser_host)

    mock_label = MagicMock(spec=QGraphicsTextItem)
    mock_label.scene.return_value = mock_parser_host.init_manager.scene
    edit3d.measurement_label_items_2d = [mock_label]

    with patch("moleditpy.ui.edit_3d_logic.sip_isdeleted_safe", return_value=False):
        edit3d.clear_2d_measurement_labels()

    mock_parser_host.init_manager.scene.removeItem.assert_called_once_with(mock_label)
    assert edit3d.measurement_label_items_2d == []


def test_clear_2d_measurement_labels_skips_deleted(mock_parser_host):
    """Verify clear_2d_measurement_labels skips sip-deleted items."""
    edit3d = DummyEdit3d(mock_parser_host)

    mock_label = MagicMock(spec=QGraphicsTextItem)
    edit3d.measurement_label_items_2d = [mock_label]

    with patch("moleditpy.ui.edit_3d_logic.sip_isdeleted_safe", return_value=True):
        edit3d.clear_2d_measurement_labels()

    mock_parser_host.init_manager.scene.removeItem.assert_not_called()
    assert edit3d.measurement_label_items_2d == []


def test_clear_2d_measurement_labels_skips_unscened(mock_parser_host):
    """Verify that items not attached to a scene are not passed to removeItem."""
    edit3d = DummyEdit3d(mock_parser_host)

    mock_label = MagicMock(spec=QGraphicsTextItem)
    mock_label.scene.return_value = None  # item has no scene
    edit3d.measurement_label_items_2d = [mock_label]

    with patch("moleditpy.ui.edit_3d_logic.sip_isdeleted_safe", return_value=False):
        edit3d.clear_2d_measurement_labels()

    mock_parser_host.init_manager.scene.removeItem.assert_not_called()
    assert edit3d.measurement_label_items_2d == []


def test_find_rdkit_atom_index_no_mol(mock_parser_host):
    """Verify find_rdkit_atom_index returns None when no molecule is loaded."""
    edit3d = DummyEdit3d(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = None

    atom_item = MagicMock()
    assert edit3d.find_rdkit_atom_index(atom_item) is None


def test_find_rdkit_atom_index_no_item(mock_parser_host):
    """Verify find_rdkit_atom_index returns None for a None atom_item."""
    edit3d = DummyEdit3d(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = MagicMock()

    assert edit3d.find_rdkit_atom_index(None) is None


def test_find_rdkit_atom_index_with_map(mock_parser_host):
    """Verify find_rdkit_atom_index returns the mapped RDKit index."""
    edit3d = DummyEdit3d(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = MagicMock()
    mock_parser_host.atom_id_to_rdkit_idx_map = {5: 2}

    atom_item = MagicMock()
    atom_item.atom_id = 5

    assert edit3d.find_rdkit_atom_index(atom_item) == 2


def test_find_rdkit_atom_index_missing_map_entry(mock_parser_host):
    """Verify find_rdkit_atom_index returns None when atom_id not in map."""
    edit3d = DummyEdit3d(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = MagicMock()
    mock_parser_host.atom_id_to_rdkit_idx_map = {5: 2}

    atom_item = MagicMock()
    atom_item.atom_id = 99  # not in map

    assert edit3d.find_rdkit_atom_index(atom_item) is None


def test_display_measurement_text_empty_clears_actor(mock_parser_host):
    """Verify display_measurement_text with empty list removes existing actor."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.measurement_text_actor = MagicMock()

    edit3d.display_measurement_text([])

    mock_parser_host.view_3d_manager.plotter.remove_actor.assert_called_once()
    assert edit3d.measurement_text_actor is None


def test_display_measurement_text_adds_actor(mock_parser_host):
    """Verify display_measurement_text calls add_text with joined lines."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.measurement_text_actor = None
    mock_parser_host.init_manager.settings["background_color"] = "#000000"  # dark

    edit3d.display_measurement_text(["Distance 1-2: 1.540 Å", "Angle: 109.5°"])

    call_kwargs = mock_parser_host.view_3d_manager.plotter.add_text.call_args
    text_arg = call_kwargs[0][0]
    assert "Distance 1-2: 1.540 Å" in text_arg
    assert "Angle: 109.5°" in text_arg
    assert call_kwargs[1].get("color") == "white"


def test_display_measurement_text_light_background(mock_parser_host):
    """Verify display_measurement_text uses black text on a light background."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.measurement_text_actor = None
    mock_parser_host.init_manager.settings["background_color"] = "#FFFFFF"  # light

    edit3d.display_measurement_text(["Distance: 1.0 Å"])

    call_kwargs = mock_parser_host.view_3d_manager.plotter.add_text.call_args
    assert call_kwargs[1].get("color") == "black"


def test_toggle_atom_selection_3d_add(mock_parser_host):
    """Verify toggle_atom_selection_3d adds an atom not yet selected."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.selected_atoms_3d = set()

    with patch.object(edit3d, "update_3d_selection_display") as mock_update:
        edit3d.toggle_atom_selection_3d(3)

    assert 3 in edit3d.selected_atoms_3d
    mock_update.assert_called_once()


def test_toggle_atom_selection_3d_remove(mock_parser_host):
    """Verify toggle_atom_selection_3d removes an already-selected atom."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.selected_atoms_3d = {3}

    with patch.object(edit3d, "update_3d_selection_display") as mock_update:
        edit3d.toggle_atom_selection_3d(3)

    assert 3 not in edit3d.selected_atoms_3d
    mock_update.assert_called_once()


def test_toggle_atom_selection_3d_idempotent_add(mock_parser_host):
    """Verify toggling different atoms accumulates them independently."""
    edit3d = DummyEdit3d(mock_parser_host)
    edit3d.selected_atoms_3d = set()

    with patch.object(edit3d, "update_3d_selection_display"):
        edit3d.toggle_atom_selection_3d(1)
        edit3d.toggle_atom_selection_3d(2)

    assert edit3d.selected_atoms_3d == {1, 2}
