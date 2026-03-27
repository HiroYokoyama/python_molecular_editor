import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.bond_length_dialog import BondLengthDialog
from moleditpy.ui.angle_dialog import AngleDialog
from moleditpy.ui.dihedral_dialog import DihedralDialog
from moleditpy.ui.alignment_dialog import AlignmentDialog
from moleditpy.ui.translation_dialog import TranslationDialog
from moleditpy.ui.move_group_dialog import MoveGroupDialog

@pytest.fixture
def mol():
    m = Chem.MolFromSmiles("CCCCCOC")
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    return m

def test_bond_length_dialog_launch(window, qtbot, mol):
    dialog = BondLengthDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    # In headless mode, isVisible() may return False. Check windowTitle as proxy for UI init.
    assert dialog.windowTitle() == "Adjust Bond Length"
    dialog.close()

def test_angle_dialog_launch(window, qtbot, mol):
    dialog = AngleDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Adjust Angle"
    dialog.close()

def test_dihedral_dialog_launch(window, qtbot, mol):
    dialog = DihedralDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Adjust Dihedral Angle"
    dialog.close()

def test_alignment_dialog_launch(window, qtbot, mol):
    dialog = AlignmentDialog(mol, window, axis="x")
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Align to X-axis"
    dialog.close()

def test_translation_dialog_launch(window, qtbot, mol):
    dialog = TranslationDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Translation"
    dialog.close()

def test_move_group_dialog_launch(window, qtbot, mol):
    dialog = MoveGroupDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Move Group"
    dialog.close()
