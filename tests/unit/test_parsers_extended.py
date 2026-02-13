import os
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_molecular_parsers import MainWindowMolecularParsers
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch, mock_open

class DummyParser(MainWindowMolecularParsers):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
    
    def __getattr__(self, name):
        return getattr(self._host, name)

    def check_unsaved_changes(self): return self._host.check_unsaved_changes()
    def clear_2d_editor(self, push_to_undo=True): self._host.clear_2d_editor(push_to_undo)

def test_load_mol_file_fallback_to_sd_supplier(mock_parser_host, tmp_path):
    """Test that non-.mol files (like .sdf) use SDMolSupplier and load the first record."""
    parser = DummyParser(mock_parser_host)
    mol1 = Chem.MolFromSmiles("C")
    mol2 = Chem.MolFromSmiles("O")
    sdf_path = tmp_path / "test.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol1)
    writer.write(mol2)
    writer.close()
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    parser.load_mol_file(str(sdf_path))
    assert len(parser.data.atoms) == 1
    assert list(parser.data.atoms.values())[0]['symbol'] == 'C'

def test_load_xyz_charge_prompt_ok(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "3\nWater\nO 0.0 0.0 0.0\nH 0.7 0.5 0.0\nH -0.7 0.5 0.0\n"
    xyz_path = tmp_path / "water.xyz"
    xyz_path.write_text(xyz_content)
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    parser.settings['skip_chemistry_checks'] = False
    
    with patch('PyQt6.QtWidgets.QDialog.exec', return_value=1):
        with patch('PyQt6.QtWidgets.QLineEdit.text', return_value="0"):
            parser.scene.find_bond_between.return_value = None
            mol = parser.load_xyz_file(str(xyz_path))
            assert mol is not None

def test_load_xyz_skip_chemistry_in_dialog(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n"
    xyz_path = tmp_path / "c2.xyz"
    xyz_path.write_text(xyz_content)
    parser.settings['skip_chemistry_checks'] = True
    mol = parser.load_xyz_file(str(xyz_path))
    assert mol is not None
    assert mol.GetProp("_xyz_skip_checks") == "1"

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
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
    parser._host.current_mol = mol
    save_path = str(tmp_path / "saved.xyz")
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=(save_path, "*.xyz")):
        parser.save_as_xyz()
    assert os.path.exists(save_path)

def test_load_mol_file_with_v2000_fix(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    # Generate a valid MOL block and strip V2000
    mol = Chem.MolFromSmiles("C")
    AllChem.Compute2DCoords(mol)
    valid_block = Chem.MolToMolBlock(mol)
    bad_block = valid_block.replace("V2000", "     ") # Replace with spaces to simulate missing tag
    
    mol_path = tmp_path / "bad.mol"
    mol_path.write_text(bad_block)
    
    parser.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    parser.view_2d.mapToScene.return_value = QPointF(0, 0)
    
    # This should trigger fix_mol_counts_line which adds V2000 back
    parser.load_mol_file(str(mol_path))
    assert len(parser.data.atoms) == 1

def test_load_xyz_file_with_manual_charge(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    xyz_content = "2\nC2\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\n"
    xyz_path = tmp_path / "c2.xyz"
    xyz_path.write_text(xyz_content)
    with patch('rdkit.Chem.rdDetermineBonds.DetermineBonds', side_effect=[RuntimeError("Failed"), None]):
        with patch('PyQt6.QtWidgets.QDialog.exec', return_value=1):
            with patch('PyQt6.QtWidgets.QLineEdit.text', return_value="0"):
                mol = parser.load_xyz_file(str(xyz_path))
                assert mol is not None

def test_save_as_mol_logic(mock_parser_host, tmp_path):
    parser = DummyParser(mock_parser_host)
    parser.data.add_atom("C", QPointF(0, 0))
    mol = Chem.MolFromSmiles("C")
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
    parser._host.current_mol = mol
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
    mock_parser_host.check_unsaved_changes.return_value = False
    with patch('PyQt6.QtWidgets.QFileDialog.getOpenFileName') as mock_dlg:
        parser.load_mol_file()
        assert not mock_dlg.called

def test_load_mol_file_dialog_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    mock_parser_host.check_unsaved_changes.return_value = True
    with patch('PyQt6.QtWidgets.QFileDialog.getOpenFileName', return_value=("", "")):
        parser.load_mol_file()
        assert not hasattr(parser, 'current_file_path') or parser.current_file_path != ""

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
        pass

def test_save_as_xyz_cancel(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch('PyQt6.QtWidgets.QFileDialog.getSaveFileName', return_value=("", "")):
        parser.save_as_xyz()
        pass

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

def test_load_mol_file_sdf_path(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    with patch('rdkit.Chem.SDMolSupplier') as mock_suppl:
        mock_suppl.return_value = [Chem.MolFromSmiles("C")]
        parser.load_mol_file("test.sdf")
        assert mock_suppl.called

def test_load_xyz_file_symbol_capitalization(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    content = "1\n\nca 0.0 0.0 0.0" 
    with patch("builtins.open", mock_open(read_data=content)):
        try:
            parser.load_xyz_file("test.xyz")
        except Exception:
            pass

def test_load_xyz_file_skip_checks_setting(mock_parser_host):
    parser = DummyParser(mock_parser_host)
    parser.settings['skip_chemistry_checks'] = True
    content = "1\n\nC 0.0 0.0 0.0"
    with patch("builtins.open", mock_open(read_data=content)), \
         patch.object(parser, 'estimate_bonds_from_distances') as mock_est:
        parser.load_xyz_file("test.xyz")
        assert mock_est.called
