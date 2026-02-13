import pytest
import os
from rdkit import Chem
from moleditpy.modules.main_window_export import MainWindowExport
from moleditpy.modules.atom_item import AtomItem
from moleditpy.modules.bond_item import BondItem
from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QColor
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
    
    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self): return self.statusBar_mock

def test_create_multi_material_obj_advanced(mock_parser_host, tmp_path):
    exporter = DummyExport(mock_parser_host)
    mesh1 = MagicMock()
    mesh1.points = np.array([[0,0,0], [1,0,0], [0,1,0]], dtype=float)
    mesh1.n_points = 3
    mesh1.n_cells = 1
    cell1 = MagicMock()
    cell1.type = 5 # VTK_TRIANGLE
    cell1.point_ids = [0, 1, 2]
    mesh1.get_cell.return_value = cell1
    
    meshes_with_colors = [{'mesh': mesh1, 'color': [255, 0, 0], 'name': 'Atom1'}]
    obj_path = str(tmp_path / "test.obj")
    mtl_path = str(tmp_path / "test.mtl")
    
    exporter.create_multi_material_obj(meshes_with_colors, obj_path, mtl_path)
    assert os.path.exists(obj_path)
    assert os.path.exists(mtl_path)

def test_export_2d_png_basic_trigger(mock_parser_host, tmp_path):
    exporter = DummyExport(mock_parser_host)
    exporter.data.add_atom("C", QPointF(100, 100))
    save_path = str(tmp_path / "export.png")
    
    item = MagicMock(spec=AtomItem)
    item.__class__ = AtomItem
    item.isVisible.return_value = True
    item.sceneBoundingRect.return_value = QRectF(90, 90, 20, 20)
    
    exporter.scene.items.return_value = [item]
    exporter.scene.backgroundBrush.return_value = QColor(255, 255, 255)
    
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.png")):
        with patch('PyQt6.QtWidgets.QMessageBox.question', return_value=MagicMock(value=16384)):
            exporter.export_2d_png()
    # verify it finishes

def test_export_2d_svg_trigger(mock_parser_host, tmp_path):
    exporter = DummyExport(mock_parser_host)
    exporter.data.add_atom("C", QPointF(100, 100))
    save_path = str(tmp_path / "export.svg")
    item = MagicMock(spec=AtomItem)
    item.__class__ = AtomItem
    item.isVisible.return_value = True
    item.sceneBoundingRect.return_value = QRectF(0, 0, 100, 100)
    exporter.scene.items.return_value = [item]
    
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.svg")):
        with patch('PyQt6.QtWidgets.QMessageBox.question', return_value=MagicMock(value=16384)): # Yes
            exporter.export_2d_svg()
    assert os.path.exists(save_path)

def test_export_stl_error_no_mol(mock_parser_host):
    exporter = DummyExport(mock_parser_host)
    exporter.current_mol = None 
    exporter.export_stl()
    exporter.statusBar().showMessage.assert_called_with("Error: Please generate a 3D structure first.")

def test_export_obj_mtl_error_no_mol(mock_parser_host):
    exporter = DummyExport(mock_parser_host)
    exporter.current_mol = None
    exporter.export_obj_mtl()
    exporter.statusBar().showMessage.assert_called_with("Error: Please generate a 3D structure first.")
