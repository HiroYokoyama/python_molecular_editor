import pytest
from PyQt6.QtCore import Qt, QPointF, QTimer
from PyQt6.QtWidgets import QMainWindow
from unittest.mock import MagicMock, patch

from moleditpy.modules.main_window_edit_actions import MainWindowEditActions
from moleditpy.modules.molecular_data import MolecularData
from moleditpy.modules.atom_item import AtomItem


class DummyEditActions(MainWindowEditActions):
    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.view_2d = host.view_2d
        self.undo_stack = MagicMock()
        self.is_xyz_derived = False
        self._ih_update_counter = 0

    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self):
        return self._host.statusBar()

    def update_window_title(self):
        pass

    def clear_2d_editor(self, push_to_undo=True):
        pass

    def update_edit_menu_actions(self):
        pass

    def create_atom_id_mapping(self):
        return {}

    def push_undo_state(self):
        pass


def test_rotate_molecule_2d_basic(mock_parser_host):
    """Test 2D rotation of the entire molecule."""
    editor = DummyEditActions(mock_parser_host)
    aid1 = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    aid2 = mock_parser_host.scene.create_atom("C", QPointF(10, 0))
    a1 = mock_parser_host.data.atoms[aid1]["item"]
    a2 = mock_parser_host.data.atoms[aid2]["item"]

    editor.scene.selectedItems.return_value = []

    with patch(
        "moleditpy.modules.main_window_edit_actions.sip_isdeleted_safe",
        return_value=False,
    ):
        editor.rotate_molecule_2d(180)

    assert a1.setPos.called
    assert a2.setPos.called


def test_resolve_overlapping_groups_basic(mock_parser_host):
    """Test overlapping group resolution."""
    editor = DummyEditActions(mock_parser_host)

    aid1 = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    a1 = mock_parser_host.data.atoms[aid1]["item"]

    aid2 = mock_parser_host.scene.create_atom(
        "O", QPointF(0, 0.1)
    )  # Within 0.5 threshold
    a2 = mock_parser_host.data.atoms[aid2]["item"]

    editor.scene.items.return_value = [a1, a2]

    with patch(
        "moleditpy.modules.main_window_edit_actions.sip_isdeleted_safe",
        return_value=False,
    ):
        editor.resolve_overlapping_groups()

    assert a2.setPos.called


def test_update_implicit_hydrogens_main_logic(mock_parser_host):
    """Test calculation of implicit hydrogens for display."""
    editor = DummyEditActions(mock_parser_host)
    aid = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    a_item = mock_parser_host.data.atoms[aid]["item"]
    # Ensure item.scene() returns something so it's not skipped
    a_item.scene.return_value = MagicMock()

    from rdkit import Chem

    mol = Chem.MolFromSmiles("C")
    for atom in mol.GetAtoms():
        atom.SetIntProp("_original_atom_id", aid)

    # Mock QTimer.singleShot to call the function immediately
    def mock_single_shot(ms, func):
        func()

    with patch.object(editor.data, "to_rdkit_mol", return_value=mol):
        with patch("PyQt6.QtCore.QTimer.singleShot", side_effect=mock_single_shot):
            editor.update_implicit_hydrogens()

    assert a_item.implicit_h_count == 4


def test_clipboard_copy_serialization(mock_parser_host):
    """Test copy selection MimeData generation."""
    editor = DummyEditActions(mock_parser_host)

    aid = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    a_item = mock_parser_host.data.atoms[aid]["item"]
    editor.scene.selectedItems.return_value = [a_item]

    mock_clipboard = MagicMock()
    with patch("PyQt6.QtWidgets.QApplication.clipboard", return_value=mock_clipboard):
        editor.copy_selection()
        assert mock_clipboard.setMimeData.called
