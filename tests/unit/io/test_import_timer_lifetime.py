"""Viewer redraws must not outlive their native window."""

from unittest.mock import MagicMock

import pytest
from PySide6.QtCore import QCoreApplication, QEvent
from PySide6.QtWidgets import QMainWindow
from rdkit import Chem

from moleditpy.ui.io_logic import IOManager


@pytest.mark.parametrize("destroy_window", [False, True])
def test_import_redraw_only_runs_while_window_exists(qtbot, destroy_window):
    host = QMainWindow()
    qtbot.addWidget(host)
    for name in (
        "state_manager",
        "edit_actions_manager",
        "view_3d_manager",
        "ui_manager",
        "set_current_molecule",
        "set_atom_id_to_rdkit_idx_map",
        "update_status_message",
        "set_current_file_path",
        "set_has_unsaved_changes",
    ):
        setattr(host, name, MagicMock())
    mol = Chem.MolFromSmiles("C")
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    manager = IOManager(host)
    manager._read_mol_or_sdf = lambda _: mol
    manager._plotter_view_isometric = MagicMock()
    manager._plotter_render = MagicMock()

    manager.load_mol_file_for_3d_viewing("fixture.mol")
    host.view_3d_manager.draw_molecule_3d.assert_called_once_with(mol)
    if destroy_window:
        host.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    qtbot.wait(160)

    expected = 0 if destroy_window else 1
    assert manager._plotter_view_isometric.call_count == expected
    assert manager._plotter_render.call_count == expected
