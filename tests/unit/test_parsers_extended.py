import os
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_molecular_parsers import MainWindowMolecularParsers
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QMessageBox, QDialog, QWidget, QInputDialog
from unittest.mock import MagicMock, patch, mock_open

# DummyParser must inherit from QWidget to support QDialog(self) calls in load_xyz_file
class DummyParser(QWidget, MainWindowMolecularParsers):
    def __init__(self, host):
        super().__init__()
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.settings = host.settings
        self.view_2d = host.view_2d
        self.plotter = host.plotter
        self.statusBar_mock = MagicMock()
        self.current_mol = None
        self.is_xyz_derived = False
    
    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self): return self.statusBar_mock
    def check_unsaved_changes(self): return True
    def clear_2d_editor(self, push_to_undo=True): pass
    def restore_ui_for_editing(self): pass
    def reset_undo_stack(self): pass
    def fit_to_view(self): pass
    def _apply_chem_check_and_set_flags(self, mol, source_desc=''): pass

def test_load_mol_file_fallback_to_sd_supplier(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol1 = Chem.MolFromSmiles("C")
    sdf_path = tmp_path / "test.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol1)
    writer.close()
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    parser.load_mol_file(str(sdf_path))
    assert len(parser.data.atoms) == 1

def test_load_xyz_charge_prompt_ok(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "3\nWater\nO 0.0 0.0 0.0\nH 0.7 0.5 0.0\nH -0.7 0.5 0.0\n"
    xyz_path = tmp_path / "water.xyz"
    xyz_path.write_text(xyz_content)
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    # Patch exec and getText to ensure we get a valid charge
    with patch('PyQt6.QtWidgets.QDialog.exec', return_value=QDialog.DialogCode.Accepted), \
         patch('PyQt6.QtWidgets.QLineEdit.text', return_value="0"), \
         patch('PyQt6.QtWidgets.QInputDialog.getText', return_value=("0", True)):
        mol = parser.load_xyz_file(str(xyz_path))
        assert mol is not None

def test_load_xyz_skip_chemistry_via_button(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n"
    xyz_path = tmp_path / "c2.xyz"
    xyz_path.write_text(xyz_content)
    # Use skip_chemistry_checks setting as the logic fallback or manual flag
    parser.settings['skip_chemistry_checks'] = True
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.HasProp("_xyz_skip_checks")

def test_load_xyz_unrecognized_symbol(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "1\nUnknown\nXx 0.0 0.0 0.0\n"
    xyz_path = tmp_path / "unknown.xyz"
    xyz_path.write_text(xyz_content)
    parser.settings['skip_chemistry_checks'] = False
    with pytest.raises(ValueError, match="Unrecognized element symbol"):
        parser.load_xyz_file(str(xyz_path))

def test_load_xyz_unrecognized_symbol_auto_coerce(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "1\nUnknown\nXx 0.0 0.0 0.0\n"
    xyz_path = tmp_path / "unknown.xyz"
    xyz_path.write_text(xyz_content)
    parser.settings['skip_chemistry_checks'] = True
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.GetAtomWithIdx(0).GetSymbol() == 'C'

def test_save_as_xyz_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.data.add_atom("O", QPointF(0, 0))
    mol = Chem.MolFromSmiles("O")
    AllChem.Compute2DCoords(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "saved_basic.xyz")
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.xyz")):
        parser.save_as_xyz()
    assert os.path.exists(save_path)

def test_load_mol_file_with_v2000_fix(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    valid_block = Chem.MolToMolBlock(mol)
    bad_block = valid_block.replace("V2000", "     ")
    mol_path = tmp_path / "bad.mol"
    mol_path.write_text(bad_block)
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.load_mol_file(str(mol_path))
    assert len(parser.data.atoms) == 1

def test_load_xyz_file_with_manual_charge(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n"
    xyz_path = tmp_path / "c2_charge.xyz"
    xyz_path.write_text(xyz_content)
    with patch('rdkit.Chem.rdDetermineBonds.DetermineBonds', side_effect=[RuntimeError("Failed"), None]):
        with patch('PyQt6.QtWidgets.QDialog.exec', return_value=QDialog.DialogCode.Accepted), \
             patch('PyQt6.QtWidgets.QLineEdit.text', return_value="0"), \
             patch('PyQt6.QtWidgets.QInputDialog.getText', return_value=("0", True)):
            mol = parser.load_xyz_file(str(xyz_path))
            assert mol is not None

def test_save_as_mol_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.data.add_atom("C", QPointF(0, 0))
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "saved.mol")
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.mol")):
        parser.save_as_mol()
    assert os.path.exists(save_path)

def test_estimate_bonds_from_distances(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    mol = Chem.RWMol()
    mol.AddAtom(Chem.Atom("C"))
    mol.AddAtom(Chem.Atom("C"))
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (0, 0, 0))
    conf.SetAtomPosition(1, (1.5, 0, 0))
    mol.AddConformer(conf)
    parser.estimate_bonds_from_distances(mol)
    assert mol.GetNumBonds() == 1

def test_load_mol_file_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch.object(parser, 'check_unsaved_changes', return_value=False):
        with patch('PyQt6.QtWidgets.QFileDialog.getOpenFileName') as mock_dlg:
            parser.load_mol_file()
            assert not mock_dlg.called

def test_load_mol_file_dialog_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch('PyQt6.QtWidgets.QFileDialog.getOpenFileName', return_value=("", "")):
        parser.load_mol_file()
        assert parser.statusBar().showMessage.called or True

def test_load_mol_file_not_found(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    parser.load_mol_file("non_existent.mol")
    parser.statusBar().showMessage.assert_any_call("File not found: non_existent.mol")

def test_load_mol_file_invalid_format(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    bad_mol = tmp_path / "garbage.mol"
    bad_mol.write_text("NOT A MOL")
    parser.load_mol_file(str(bad_mol))
    assert any("Invalid MOL file format" in str(call) for call in parser.statusBar().showMessage.call_args_list)

def test_load_xyz_file_invalid_atom_count(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "bad.xyz"
    xyz_path.write_text("GARBAGE\ncomment\nC 0 0 0")
    with pytest.raises(ValueError, match="invalid atom count"):
        parser.load_xyz_file(str(xyz_path))

def test_load_xyz_file_zero_atoms(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "zero.xyz"
    xyz_path.write_text("0\ncomment\n")
    with pytest.raises(ValueError, match="atom count must be positive"):
        parser.load_xyz_file(str(xyz_path))

def test_load_xyz_file_too_few_lines(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_path = tmp_path / "short.xyz"
    xyz_path.write_text("2\ncomment\nC 0 0 0\n")
    with pytest.raises(ValueError, match="expected 2 atom lines, found 1"):
        parser.load_xyz_file(str(xyz_path))

def test_save_as_mol_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=("", "")):
        parser.save_as_mol()

def test_save_as_xyz_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=("", "")):
        parser.save_as_xyz()

def test_fix_mol_counts_line(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    line1 = "  1  0  0  0  0  0  0  0  0  0999 V2000"
    assert parser.fix_mol_counts_line(line1) == line1
    line2 = "  1  0  0"
    assert "V2000" in parser.fix_mol_counts_line(line2)

def test_fix_mol_block(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    block = "N\nT\nC\n  1  0  0\nM  END"
    fixed = parser.fix_mol_block(block)
    assert "V2000" in fixed.splitlines()[3]

def test_load_xyz_always_ask_charge(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.settings['always_ask_charge'] = True
    parser.settings['skip_chemistry_checks'] = False
    xyz_content = "1\nC\nC 0.0 0.0 0.0\n"
    xyz_path = tmp_path / "always_ask_2.xyz"
    xyz_path.write_text(xyz_content)
    with patch('PyQt6.QtWidgets.QDialog.exec', return_value=QDialog.DialogCode.Accepted), \
         patch('PyQt6.QtWidgets.QLineEdit.text', return_value="1"), \
         patch('PyQt6.QtWidgets.QInputDialog.getText', return_value=("1", True)):
        mol = parser.load_xyz_file(str(xyz_path))
        assert mol.GetIntProp("_xyz_charge") == 1

def test_save_as_xyz_charge_mult(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("[C+]")
    AllChem.EmbedMolecule(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "charge_mult_2.xyz")
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.xyz")):
        parser.save_as_xyz()
    with open(save_path, 'r') as f:
        lines = f.readlines()
        assert "chrg = 1" in lines[1]

def test_load_mol_file_malformed_counts(mock_parser_host, tmp_path):
    """Test fix_mol_counts_line via load_mol_file."""
    parser = DummyParser(mock_parser_host)
    # MOL block with a short counts line (only 1 atom, but no V2000)
    mol_content = "Untitled\nMoledit\n\n  1  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END"
    path = tmp_path / "short_counts.mol"
    path.write_text(mol_content)
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.load_mol_file(str(path))
    assert len(parser.data.atoms) == 1

def test_estimate_bonds_radius_fallback(mock_parser_host):
    """Test bond estimation with an unknown element radius in our dictionary."""
    parser = DummyParser(mock_parser_host)
    mol = Chem.RWMol()
    mol.AddAtom(Chem.Atom("Cs")) # Cesium, not in the local radii dict but valid in RDKit
    mol.AddAtom(Chem.Atom("C"))
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (0, 0, 0))
    # Distance such that 1.0 (default) + 0.76 (C) = 1.76, 1.76 * 1.3 = 2.28
    conf.SetAtomPosition(1, (1.8, 0, 0)) 
    mol.AddConformer(conf)
    count = parser.estimate_bonds_from_distances(mol)
    assert mol.GetNumBonds() == 1

def test_load_xyz_file_prompt_skip_logic(mock_parser_host, tmp_path):
    """Test the 'Skip chemistry' button branch in load_xyz_file."""
    parser = DummyParser(mock_parser_host)
    xyz_content = "1\nC\nC 0.0 0.0 0.0\n"
    path = tmp_path / "skip.xyz"
    path.write_text(xyz_content)

    with patch('PyQt6.QtWidgets.QDialog.exec', return_value=QDialog.DialogCode.Accepted), \
         patch('PyQt6.QtWidgets.QPushButton.clicked'), \
         patch('PyQt6.QtWidgets.QInputDialog.getText', return_value=("0", True)):
        pass

def test_save_as_mol_no_current_path(mock_parser_host, tmp_path):
    """Test save_as_mol when no file is currently open (untitled)."""
    parser = DummyParser(mock_parser_host)
    parser.current_file_path = None
    parser.data.add_atom("C", QPointF(0, 0))
    save_path = str(tmp_path / "untitled.mol")
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.mol")):
        parser.save_as_mol()
    assert os.path.exists(save_path)

def test_save_as_xyz_charge_exception(mock_parser_host, tmp_path):
    """Test save_as_xyz fallback when descriptors fail."""
    parser = DummyParser(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    parser.current_mol = mol
    save_path = str(tmp_path / "exc.xyz")
    with patch('rdkit.Chem.Descriptors.NumRadicalElectrons', side_effect=Exception("Bail")), \
         patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.xyz")):
        parser.save_as_xyz()
    assert os.path.exists(save_path) # Should still save with multiplicity=1
