import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.planarize_dialog import PlanarizeDialog
from moleditpy.modules.mirror_dialog import MirrorDialog
from moleditpy.modules.align_plane_dialog import AlignPlaneDialog
from moleditpy.modules.constrained_optimization_dialog import ConstrainedOptimizationDialog
from moleditpy.modules.periodic_table_dialog import PeriodicTableDialog
# from moleditpy.modules.color_settings_dialog import ColorSettingsDialog # Color settings might be in a different module

@pytest.fixture
def mol():
    m = Chem.MolFromSmiles("CCCCCOC")
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    return m

def test_planarize_dialog_launch(window, qtbot, mol):
    dialog = PlanarizeDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Planarize"
    dialog.close()

def test_mirror_dialog_launch(window, qtbot, mol):
    dialog = MirrorDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Mirror Molecule"
    dialog.close()

@pytest.mark.parametrize("plane", ["xy", "xz", "yz"])
def test_align_plane_dialog_launch(window, qtbot, mol, plane):
    plane_names = {"xy": "XY", "xz": "XZ", "yz": "YZ"}
    dialog = AlignPlaneDialog(mol, window, plane)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == f"Align to {plane_names[plane]} Plane"
    dialog.close()

def test_constrained_optimization_dialog_launch(window, qtbot, mol):
    dialog = ConstrainedOptimizationDialog(mol, window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Constrained Optimization"
    dialog.close()

def test_periodic_table_dialog_launch(window, qtbot):
    dialog = PeriodicTableDialog(window)
    qtbot.add_widget(dialog)
    dialog.show()
    assert dialog.windowTitle() == "Select an Element"
    dialog.close()
