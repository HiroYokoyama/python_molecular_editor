import pytest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_export import MainWindowExport
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QMessageBox
from unittest.mock import MagicMock, patch
import numpy as np


class DummyExport(MainWindowExport):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.plotter = host.plotter
        self.view_2d = host.view_2d
        self.statusBar_mock = MagicMock()
        self.current_file_path = None
        self.current_mol = host.current_mol

    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self):
        return self.statusBar_mock


def test_create_multi_material_obj_advanced(mock_parser_host, tmp_path):
    """Verify that create_multi_material_obj creates both .obj and .mtl files."""
    exporter = DummyExport(mock_parser_host)
    mesh1 = MagicMock()
    mesh1.points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    mesh1.n_points = 3
    mesh1.n_cells = 1
    cell1 = MagicMock()
    cell1.type = 5  # VTK_TRIANGLE
    cell1.point_ids = [0, 1, 2]
    mesh1.get_cell.return_value = cell1

    meshes_with_colors = [{"mesh": mesh1, "color": [255, 0, 0], "name": "Atom1"}]
    obj_path = str(tmp_path / "test.obj")
    mtl_path = str(tmp_path / "test.mtl")

    exporter.create_multi_material_obj(meshes_with_colors, obj_path, mtl_path)
    assert os.path.exists(obj_path)
    assert os.path.exists(mtl_path)
    assert os.path.getsize(obj_path) > 0
    assert os.path.getsize(mtl_path) > 0
    with open(obj_path, 'r') as f:
        content = f.read()
        assert 'mtllib' in content
        assert 'v ' in content


def test_export_2d_png_basic_trigger(mock_parser_host, tmp_path):
    """Verify that export_2d_png successfully creates a PNG file after user confirmation."""
    exporter = DummyExport(mock_parser_host)
    exporter.data.add_atom("C", QPointF(100, 100))
    save_path = str(tmp_path / "export.png")

    item = MagicMock(spec=AtomItem)
    item.__class__ = AtomItem
    item.isVisible.return_value = True
    item.sceneBoundingRect.return_value = QRectF(90, 90, 20, 20)

    exporter.scene.items.return_value = [item]
    exporter.scene.backgroundBrush.return_value = QColor(255, 255, 255)

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "*.png")
    ):
        with patch(
            "PyQt6.QtWidgets.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ):
            exporter.export_2d_png()
            # Verify file exists and is not empty
            assert os.path.exists(save_path)
            assert os.path.getsize(save_path) > 0


def test_export_2d_svg_trigger(mock_parser_host, tmp_path):
    """Verify that export_2d_svg successfully creates an SVG file."""
    exporter = DummyExport(mock_parser_host)
    exporter.data.add_atom("C", QPointF(100, 100))
    save_path = str(tmp_path / "export.svg")
    item = MagicMock(spec=AtomItem)
    item.__class__ = AtomItem
    item.isVisible.return_value = True
    item.sceneBoundingRect.return_value = QRectF(0, 0, 100, 100)
    exporter.scene.items.return_value = [item]

    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "*.svg")
    ):
        with patch(
            "PyQt6.QtWidgets.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ):
            exporter.export_2d_svg()
    assert os.path.exists(save_path)
    assert os.path.getsize(save_path) > 0
    with open(save_path, 'r') as f:
        assert '<svg' in f.read().lower()


def test_export_stl_error_no_mol(mock_parser_host):
    """Verify that export_stl shows an error message when no molecule is present."""
    exporter = DummyExport(mock_parser_host)
    exporter.current_mol = None
    exporter.export_stl()
    exporter.statusBar().showMessage.assert_any_call(
        "Error: Please generate a 3D structure first."
    )


def test_export_obj_mtl_error_no_mol(mock_parser_host):
    """Verify that export_obj_mtl shows an error message when no molecule is present."""
    exporter = DummyExport(mock_parser_host)
    exporter.current_mol = None
    exporter.export_obj_mtl()
    exporter.statusBar().showMessage.assert_any_call(
        "Error: Please generate a 3D structure first."
    )


def test_export_stl_success_trigger(mock_parser_host, tmp_path):
    """Verify that export_stl triggers the mesh save logic for a valid 3D molecule."""
    exporter = DummyExport(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    exporter.current_mol = mol
    save_path = str(tmp_path / "test.stl")
    mesh = MagicMock()
    mesh.n_points = 100  # Ensure it has points

    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
            return_value=(save_path, "*.stl"),
        ),
        patch.object(exporter, "export_from_3d_view_no_color", return_value=mesh),
    ):
        exporter.export_stl()
        assert mesh.save.called
    # Verify save called with correct filename and binary mode
    args, kwargs = mesh.save.call_args
    assert args[0] == save_path
    assert kwargs.get("binary") is True


def test_export_obj_mtl_success_trigger(mock_parser_host, tmp_path):
    """Verify that export_obj_mtl triggers the multi-material OBJ creation logic."""
    exporter = DummyExport(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    exporter.current_mol = mol
    save_path = str(tmp_path / "test.obj")
    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
            return_value=(save_path, "*.obj"),
        ) as mock_file_dialog,
        patch.object(
            exporter,
            "export_from_3d_view_with_colors",
            return_value=[{"mesh": MagicMock(), "color": [1, 1, 1], "name": "A"}],
        ),
        patch.object(exporter, "create_multi_material_obj") as mock_create,
    ):
        exporter.export_obj_mtl()
        assert mock_file_dialog.called # Explicit assert for file dialog call
        assert mock_create.called


def test_export_3d_png_logic(mock_parser_host, tmp_path):
    """Verify that export_3d_png triggers the plotter screenshot logic."""
    exporter = DummyExport(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    exporter.current_mol = mol
    save_path = str(tmp_path / "test3d.png")
    mock_plotter = MagicMock()
    exporter.plotter = mock_plotter
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "*.png")
    ):
        with patch(
            "PyQt6.QtWidgets.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ):
            exporter.export_3d_png()
            assert mock_plotter.screenshot.called


def test_export_color_stl_logic(mock_parser_host, tmp_path):
    """Verify that export_color_stl triggers the status bar message (success indicator)."""
    exporter = DummyExport(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    exporter.current_mol = mol
    save_path = str(tmp_path / "color.stl")
    mesh = MagicMock()
    mesh.n_points = 100
    with (
        patch(
            "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
            return_value=(save_path, "*.stl"),
        ),
        patch.object(
            exporter,
            "export_from_3d_view_with_colors",
            return_value=[{"mesh": mesh, "color": [1, 1, 1], "name": "A"}],
        ),
    ):
        exporter.export_color_stl()
        assert exporter.statusBar().showMessage.called


def test_export_from_3d_view_with_colors_complex_splitting(mock_parser_host):
    """Test the complex logic of splitting a mesh by per-vertex colors."""
    exporter = DummyExport(mock_parser_host)

    # Create a mock actor
    actor = MagicMock()
    actor_name = "test_actor"

    # Create a mock PyVista/VTK mesh with point data colors
    mesh = MagicMock()
    mesh.n_points = 4

    # Simulate PolyData behavior
    mesh.copy.return_value = mesh
    mesh.point_data = {
        # 4 vertices, 2 red, 2 blue
        "red": np.array([255, 255, 0, 0]),
        "green": np.array([0, 0, 0, 0]),
        "blue": np.array([0, 0, 255, 255]),
    }

    # extraction result simulation
    def extract_points_side_effect(point_inds, adjacent_cells=False):
        sub = MagicMock()
        sub.n_points = len(point_inds)
        return sub

    mesh.extract_points.side_effect = extract_points_side_effect

    # Setup renderer actors
    exporter.plotter.renderer.actors = {actor_name: actor}

    # Setup mapper to return our mesh
    actor.mapper.input = mesh

    with patch(
        "moleditpy.modules.main_window_export.pv.PolyData", MagicMock
    ):  # used for type check
        res = exporter.export_from_3d_view_with_colors()

        # Expect 2 submeshes: one Red, one Blue
        assert len(res) == 2
        # Check colors
        c1 = res[0]["color"]
        c2 = res[1]["color"]

        # One should be [255,0,0], other [0,0,255]
        colors = sorted([tuple(c1), tuple(c2)])
        assert colors[0] == (0, 0, 255)
        assert colors[1] == (255, 0, 0)


def test_export_2d_png_hides_items(mock_parser_host, tmp_path):
    """Test that export_2d_png hides non-atom items and restores them."""
    exporter = DummyExport(mock_parser_host)
    save_path = str(tmp_path / "hide_test.png")

    # 1. Add atom (should stay visible)
    atom = MagicMock(spec=AtomItem)
    atom.isVisible.return_value = True
    atom.sceneBoundingRect.return_value = QRectF(0, 0, 10, 10)

    # 2. Add random item (should be hidden)
    # e.g. a Selection Box or temporary line
    other_item = MagicMock()
    other_item.isVisible.return_value = True
    # Not instance of AtomItem/BondItem

    exporter.scene.items.return_value = [atom, other_item]

    # Ensure guard clause (if not self.data.atoms) passes
    exporter.data.atoms = {0: "dummy"}

    # We mock QFileDialog and render to avoid actual GUI dependency
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(save_path, "*.png")
    ):
        with patch(
            "PyQt6.QtWidgets.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ):
            with patch("PyQt6.QtGui.QPainter"):  # Mock painter
                with patch("PyQt6.QtGui.QImage"):  # Mock image
                    exporter.export_2d_png()

    # Verify other_item was process (hidden then restored)
    # The code calls item.hide() then item.setVisible(was_visible)
    assert other_item.hide.called
    other_item.setVisible.assert_called_with(True)
