import pytest
from moleditpy.ui.angle_dialog import AngleDialog
from moleditpy.ui.dihedral_dialog import DihedralDialog
from rdkit import Chem
from unittest.mock import MagicMock

def test_angle_dialog_wrapping(qtbot):
    mol = Chem.MolFromSmiles("CCO")
    main_window = MagicMock()
    dialog = AngleDialog(mol, main_window)
    qtbot.addWidget(dialog)

    # Set selection
    dialog.atom1_idx = 0
    dialog.atom2_idx = 1
    dialog.atom3_idx = 2

    # Set text and apply
    dialog.angle_input.setText("190")
    dialog.adjust_angle = MagicMock()
    
    dialog.apply_changes()
    
    assert dialog.angle_input.text() == "-170.00"
    dialog.adjust_angle.assert_called_once_with(-170.0)

def test_dihedral_dialog_wrapping(qtbot):
    mol = Chem.MolFromSmiles("CCCC")
    main_window = MagicMock()
    dialog = DihedralDialog(mol, main_window)
    qtbot.addWidget(dialog)

    # Set selection
    dialog.atom1_idx = 0
    dialog.atom2_idx = 1
    dialog.atom3_idx = 2
    dialog.atom4_idx = 3

    # Set text and apply
    dialog.dihedral_input.setText("-200")
    dialog.adjust_dihedral = MagicMock()
    
    dialog.apply_changes()
    
    assert dialog.dihedral_input.text() == "160.00"
    dialog.adjust_dihedral.assert_called_once_with(160.0)
