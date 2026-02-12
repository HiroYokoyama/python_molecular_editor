import os
from rdkit import Chem
from moleditpy.modules.main_window_molecular_parsers import MainWindowMolecularParsers
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch

class DummyParser(MainWindowMolecularParsers):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
    
    def __getattr__(self, name):
        # Delegate to host for any status bar, actions, or other main window methods
        return getattr(self._host, name)

    def check_unsaved_changes(self): return self._host.check_unsaved_changes()
    def clear_2d_editor(self, push_to_undo=True): self._host.clear_2d_editor(push_to_undo)

def test_fix_mol_block(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    invalid_counts = " 3 2  0  0  0  0  0  0  0  0\n"
    mol_block = "\n  Title\n\n" + invalid_counts + "  0.0 0.0 0.0 C\n"
    fixed = parser.fix_mol_block(mol_block)
    lines = fixed.splitlines()
    assert "V2000" in lines[3]
    assert len(lines[3]) >= 39

def test_load_mol_file_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    ref_mol = Chem.MolFromSmiles("CC")
    mol_content = Chem.MolToMolBlock(ref_mol)
    mol_file = tmp_path / "test.mol"
    mol_file.write_text(mol_content)
    
    # Mock view properties
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    parser.load_mol_file(str(mol_file))
    
    assert len(parser.data.atoms) == 2
    assert any(a['symbol'] == 'C' for a in parser.data.atoms.values())

def test_xyz_parsing_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    # NO LEADING SPACES in XYZ
    xyz_content = "3\nWater\nO 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n"
    xyz_file = tmp_path / "water.xyz"
    xyz_file.write_text(xyz_content)
    
    mol = parser.load_xyz_file(str(xyz_file))
    assert mol is not None
    assert mol.GetNumAtoms() == 3
    assert any(a.GetSymbol() == 'O' for a in mol.GetAtoms())

def test_load_xyz_file_with_estimation(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "2\nEthane\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n"
    xyz_file = tmp_path / "ethane.xyz"
    xyz_file.write_text(xyz_content)
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    mol = parser.load_xyz_file(str(xyz_file))
    assert mol is not None
    assert mol.GetNumAtoms() == 2
    assert mol.GetNumBonds() == 1

def test_save_as_mol_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.data.add_atom("C", QPointF(0, 0))
    save_path = str(tmp_path / "saved.mol")
    
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.mol")):
        parser.save_as_mol()
    
    assert os.path.exists(save_path)
    content = open(save_path).read()
    assert "C  " in content
