import os
import json
from rdkit import Chem
from moleditpy.ui.io_logic import IOManager
from moleditpy.core.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch


class DummyProjectIo(IOManager):  # IOManager is same class
    def __init__(self, host):
        self._host = host
        self.host = self  # self-referential so host-attr writes land on io
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene

    def __getattr__(self, name):
        return getattr(self._host, name)

    def check_unsaved_changes(self):
        return self._host.check_unsaved_changes()

    def clear_2d_editor(self, push_to_undo=True):
        self._host.clear_2d_editor(push_to_undo)

    def create_json_data(self):
        return self._host.create_json_data()


def test_project_save_load_logic(mock_parser_host, tmp_path):
    """Verify project save and load logic for .pmeprj files."""
    io_handler = DummyProjectIo(mock_parser_host)
    # Use scene.create_atom to ensure atoms have 'item' mocks (required for create_json_data)
    io_handler.scene.create_atom("C", QPointF(10.5, 20.7), charge=1)

    project_file = str(tmp_path / "test.pmeprj")
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName",
        return_value=(project_file, "*.pmeprj"),
    ):
        io_handler.save_project_as()

    assert os.path.exists(project_file)

    # Reload
    io_handler.data.atoms.clear()
    io_handler.data.bonds.clear()
    io_handler.load_json_data(project_file)

    assert len(io_handler.data.atoms) == 1
    assert next(iter(io_handler.data.atoms.values()))["charge"] == 1


def test_xyz_export_logic(mock_parser_host, tmp_path):
    """Verify XYZ export logic."""
    io_handler = DummyProjectIo(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    io_handler.current_mol = mol

    xyz_file = str(tmp_path / "test.xyz")
    with patch(
        "PyQt6.QtWidgets.QFileDialog.getSaveFileName", return_value=(xyz_file, "*.xyz")
    ):
        io_handler.save_as_xyz()

    assert os.path.exists(xyz_file)
    xyz_mol = Chem.MolFromXYZFile(xyz_file)
    assert xyz_mol.GetNumAtoms() == mol.GetNumAtoms()
