"""Axis alignment and replacement with real Qt and VTK rendering."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from moleditpy.ui.alignment_dialog import AlignmentDialog

pytestmark = pytest.mark.full_gui


def test_axis_alignment_and_replacement_with_native_renderer(
    full_window, qtbot, tmp_path
):
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    full_window.set_current_molecule(mol)
    full_window.view_3d_manager.draw_molecule_3d(mol)
    full_window.view_3d_manager.plotter.reset_camera()
    dialog = AlignmentDialog(mol, full_window, "x", [0, 1], parent=full_window)
    full_window.dialog_manager._open_3d_edit_dialog(dialog)
    qtbot.waitExposed(dialog)

    dialog.apply_alignment()

    positions = mol.GetConformer().GetPositions()
    np.testing.assert_allclose((positions[1] - positions[0])[1:], 0, atol=1e-6)
    assert full_window.grab().save(str(tmp_path / "axis-alignment.png"))
    assert dialog.grab().save(str(tmp_path / "axis-dialog.png"))
    replacement = Chem.Mol(mol)
    full_window.set_current_molecule(replacement)
    full_window.view_3d_manager.draw_molecule_3d(replacement)
    assert dialog._invalidated
    assert not dialog.isVisible()
    assert not full_window.edit_3d_manager.active_3d_dialogs
