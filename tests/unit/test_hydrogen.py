import pytest
from moleditpy.ui.edit_actions_logic import EditActionsManager
from moleditpy.core.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch


class DummyEditActions(EditActionsManager):
    def __init__(self, host):
        super().__init__(host)
        self._host = host
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene
        self.view_2d = host.init_manager.view_2d
        self.plotter = host.view_3d_manager.plotter
        self.settings = host.init_manager.settings
        self.current_mol = host.view_3d_manager.current_mol
        self.host.init_manager.current_file_path = host.init_manager.current_file_path
        self.host.state_manager.has_unsaved_changes = host.state_manager.has_unsaved_changes
        self.main_window_edit_actions = self

    def statusBar(self):
        return self._host.statusBar()

    def push_undo_state(self):
        self._host.push_undo_state()


def test_add_hydrogen_atoms_app_logic(mock_parser_host):
    """Verify add_hydrogen_atoms creates items in the scene based on RDKit logic."""
    actions = DummyEditActions(mock_parser_host)
    # Add a Carbon atom via scene to trigger item mock setup
    c_id = actions.scene.create_atom("C", QPointF(100, 100))

    # Run the app logic
    actions.add_hydrogen_atoms()

    # scene.create_atom should have been called 4 times for 'H'
    h_calls = [
        c
        for c in actions.scene.create_atom.call_args_list
        if c.args[0] == "H" or c.kwargs.get("symbol") == "H"
    ]
    assert len(h_calls) == 4
    # scene.create_bond should have 4 calls
    assert actions.scene.create_bond.call_count == 4


def test_remove_hydrogen_atoms_app_logic(mock_parser_host):
    """Verify remove_hydrogen_atoms finds and deletes H items using app logic."""
    actions = DummyEditActions(mock_parser_host)

    # Add C and H via scene
    c_id = actions.scene.create_atom("C", QPointF(100, 100))
    h_id = actions.scene.create_atom("H", QPointF(150, 100))

    # Create bond using items (mocked items should exist now)
    actions.scene.create_bond(
        actions.data.atoms[c_id]["item"], actions.data.atoms[h_id]["item"]
    )

    # We need to ensure sip_isdeleted_safe doesn't block deletion in test
    with patch("moleditpy.ui.edit_actions_logic.sip_isdeleted_safe", return_value=False):
        actions.remove_hydrogen_atoms()

    # Should call scene.delete_items with a set containing the H item
    actions.scene.delete_items.assert_called()
    deleted_set = actions.scene.delete_items.call_args[0][0]
    # In mock_parser_host, host.state_manager.data.atoms[h_id]['item'] is a MagicMock
    h_item = actions.data.atoms[h_id]["item"]
    assert h_item in deleted_set
    assert actions.data.atoms[c_id]["item"] not in deleted_set
