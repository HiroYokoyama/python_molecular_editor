"""Viewer imports must protect the old document and provide an undo baseline."""

from unittest.mock import MagicMock

import numpy as np
import pytest
from PyQt6.QtCore import QPointF
from rdkit import Chem
from rdkit.Chem import AllChem


@pytest.mark.parametrize("source", ["xyz_file", "xyz_text", "mol_file"])
def test_cancel_viewer_import_preserves_document(window, monkeypatch, source):
    """Cancel leaves 2D content, dirty state, history, and molecule untouched."""
    window.scene.create_atom("N", QPointF(0, 0))
    window.edit_actions_manager.push_undo_state()
    before = window.state_manager.get_current_state()
    history = list(window.edit_actions_manager.undo_stack)
    monkeypatch.setattr(window.state_manager, "check_unsaved_changes", lambda: False)
    parser = MagicMock()
    for name in ("load_xyz_file", "load_xyz_block", "_read_mol_or_sdf"):
        monkeypatch.setattr(window.io_manager, name, parser)
    if source == "xyz_file":
        window.io_manager.load_xyz_for_3d_viewing("cancel.xyz")
    elif source == "xyz_text":
        assert window.io_manager.show_xyz_data("cancel") is None
    else:
        window.io_manager.load_mol_file_for_3d_viewing("cancel.mol")
    parser.assert_not_called()
    assert window.state_manager.get_current_state() == before
    assert window.edit_actions_manager.undo_stack == history
    assert window.state_manager.has_unsaved_changes


@pytest.mark.parametrize("source", ["xyz_file", "xyz_text", "mol_file"])
def test_first_edit_undo_restores_imported_geometry(window, monkeypatch, source):
    """Undo/redo after import preserves the molecule, coordinates and viewer mode."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    original = mol.GetConformer().GetPositions().copy()
    monkeypatch.setattr(window.state_manager, "check_unsaved_changes", lambda: True)
    for name in ("load_xyz_file", "load_xyz_block", "_read_mol_or_sdf"):
        monkeypatch.setattr(window.io_manager, name, lambda _: mol)
    monkeypatch.setattr(window.view_3d_manager, "draw_molecule_3d", MagicMock())
    monkeypatch.setattr(window.io_manager, "_plotter_view_isometric", lambda: None)
    monkeypatch.setattr(window.io_manager, "_plotter_render", lambda: None)
    if source == "xyz_file":
        window.io_manager.load_xyz_for_3d_viewing("import.xyz")
    elif source == "xyz_text":
        assert window.io_manager.show_xyz_data("data") is mol
    else:
        window.io_manager.load_mol_file_for_3d_viewing("import.mol")
    assert not window.state_manager.has_unsaved_changes
    assert len(window.edit_actions_manager.undo_stack) == 1
    assert not window.ui_manager.is_2d_editable
    changed = original + [1.0, 2.0, 3.0]
    for i, pos in enumerate(changed):
        mol.GetConformer().SetAtomPosition(i, pos)
    window.edit_actions_manager.push_undo_state()
    window.edit_actions_manager.undo()
    np.testing.assert_allclose(
        window.current_mol.GetConformer().GetPositions(), original
    )
    assert not window.ui_manager.is_2d_editable
    window.edit_actions_manager.redo()
    np.testing.assert_allclose(
        window.current_mol.GetConformer().GetPositions(), changed
    )
