"""Unit tests for EditActionsManager 2D edit operations."""

import os
import sys
import pytest
import numpy as np
from PyQt6.QtWidgets import QDialog
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch

from moleditpy.ui.edit_actions_logic import EditActionsManager
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem
from rdkit import Chem

# Ensure local moleditpy is discoverable
_workspace_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_workspace_src) and _workspace_src not in sys.path:
    sys.path.insert(0, _workspace_src)


class DummyEditActions(EditActionsManager):
    def __init__(self, host):
        super().__init__(host)
        self._host = host
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene
        self.view_2d = host.init_manager.view_2d
        self.edit_actions_manager.undo_stack = MagicMock()
        self.is_xyz_derived = False
        self.main_window_edit_actions = self
        self.ih_update_counter = 0

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
    a1 = mock_parser_host.scene.atom_items[aid1]
    a2 = mock_parser_host.scene.atom_items[aid2]

    editor.scene.selectedItems.return_value = []

    with patch(
        "moleditpy.ui.edit_actions_logic.sip_isdeleted_safe",
        return_value=False,
    ):
        editor.rotate_molecule_2d(180)

    assert a1.setPos.called
    assert a2.setPos.called


def test_resolve_overlapping_groups_basic(mock_parser_host):
    """Test overlapping group resolution."""
    editor = DummyEditActions(mock_parser_host)

    aid1 = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    a1 = mock_parser_host.scene.atom_items[aid1]

    aid2 = mock_parser_host.scene.create_atom(
        "O", QPointF(0, 0.1)
    )  # Within 0.5 threshold
    a2 = mock_parser_host.scene.atom_items[aid2]

    editor.scene.items.return_value = [a1, a2]

    with patch(
        "moleditpy.ui.edit_actions_logic.sip_isdeleted_safe",
        return_value=False,
    ):
        editor.resolve_overlapping_groups()

    assert a2.setPos.called


def test_update_implicit_hydrogens_main_logic(mock_parser_host):
    """Test calculation of implicit hydrogens for display."""
    editor = DummyEditActions(mock_parser_host)
    aid = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    a_item = mock_parser_host.scene.atom_items[aid]
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
    a_item = mock_parser_host.scene.atom_items[aid]
    editor.scene.selectedItems.return_value = [a_item]

    mock_clipboard = MagicMock()
    with patch("PyQt6.QtWidgets.QApplication.clipboard", return_value=mock_clipboard):
        editor.copy_selection()
        assert mock_clipboard.setMimeData.called


# =============================================================================
# Extended tests (merged from test_edit_actions_extended.py)
# =============================================================================


class DummyHost:
    def __init__(self):
        self.statusBar_mock = MagicMock()
        self.settings = {}
        self.is_xyz_derived = False
        self.chem_check_tried = False
        self.chem_check_failed = False
        self.initialization_complete = True
        self._is_restoring_state = False

        self.init_manager = MagicMock()
        self.init_manager.settings = self.settings
        self.state_manager = MagicMock()
        self.ui_manager = MagicMock()
        self.view_3d_manager = MagicMock()
        self.edit_3d_manager = MagicMock()
        self.plugin_manager = MagicMock()
        self.edit_actions_manager = None

        self.init_manager.scene = MagicMock()
        self.state_manager.data = MagicMock()
        self.init_manager.scene.data = self.state_manager.data
        self.init_manager.view_2d = MagicMock()
        self.view_3d_manager.plotter = MagicMock()

        self.init_manager.scene.atom_items = {}
        self.init_manager.scene.bond_items = {}

        self.init_manager.cut_action = MagicMock()
        self.init_manager.copy_action = MagicMock()
        self.init_manager.paste_action = MagicMock()
        self.init_manager.undo_action = MagicMock()
        self.init_manager.redo_action = MagicMock()
        self.init_manager.measurement_action = MagicMock()
        self.init_manager.edit_3d_action = MagicMock()
        self.init_manager.optimize_3d_button = MagicMock()

    def statusBar(self):
        return self.statusBar_mock

    def push_undo_state(self):
        pass

    def update_window_title(self):
        pass

    def check_unsaved_changes(self):
        return True

    def restore_ui_for_editing(self):
        pass

    def clear_3d_selection(self):
        pass

    def clear_2d_editor(self, push_to_undo=True):
        pass

    def reset_undo_stack(self):
        pass

    def reset_zoom(self):
        pass

    def update_atom_id_menu_text(self):
        pass

    def update_atom_id_menu_state(self):
        pass

    def toggle_measurement_mode(self, state):
        pass

    def toggle_3d_edit_mode(self, state):
        pass

    def _enable_3d_features(self, state):
        pass

    def clear_2d_measurement_labels(self):
        pass

    def update_2d_measurement_labels(self):
        pass

    def set_current_molecule(self, mol):
        self.view_3d_manager.current_mol = mol

    def clear_3d_view(self):
        self.view_3d_manager.current_mol = None

    def set_constraints_3d(self, constraints):
        self.edit_3d_manager.constraints_3d = constraints

    def get_constraints_3d(self):
        return self.edit_3d_manager.constraints_3d

    def set_has_unsaved_changes(self, value):
        self.state_manager.has_unsaved_changes = value

    def set_current_file_path(self, path):
        self.init_manager.current_file_path = path

    def get_current_file_path(self):
        return self.init_manager.current_file_path

    def set_molecule_data(self, data):
        self.state_manager.data = data
        if self.init_manager.scene:
            self.init_manager.scene.data = data

    def get_molecule_data(self):
        return self.state_manager.data

    def update_status_message(self, message, timeout=0):
        if timeout == 0:
            self.statusBar_mock.showMessage(message)
        else:
            self.statusBar_mock.showMessage(message, timeout)

    def get_settings(self):
        return self.init_manager.settings


class TestEditActionsExtended:
    @pytest.fixture
    def manager(self):
        host = DummyHost()
        mgr = EditActionsManager(host)
        host.edit_actions_manager = mgr

        if not hasattr(sys.modules[mgr.__module__], "CLIPBOARD_MIME_TYPE"):
            setattr(
                sys.modules[mgr.__module__],
                "CLIPBOARD_MIME_TYPE",
                "application/x-moleditpy-pme",
            )

        return mgr

    def test_apply_chem_check_force_skip(self, manager):
        """force_skip=True bypasses the chemistry check without setting any flags."""
        mol = Chem.MolFromSmiles("C")
        manager.apply_chem_check_and_set_flags(mol, force_skip=True)
        assert manager.host.chem_check_tried is False
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_settings_skip(self, manager):
        """skip_chemistry_checks=True in settings prevents the check from running."""
        manager.host.init_manager.settings["skip_chemistry_checks"] = True
        mol = Chem.MolFromSmiles("C")
        manager.apply_chem_check_and_set_flags(mol)
        assert manager.host.chem_check_tried is False
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_success(self, manager):
        """A valid molecule sets chem_check_tried=True and chem_check_failed=False."""
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")
        manager.apply_chem_check_and_set_flags(mol)
        assert manager.host.chem_check_tried is True
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_failure(self, manager):
        """A sanitization error sets chem_check_failed=True and disables the 3D optimize button."""
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")
        with patch(
            "rdkit.Chem.SanitizeMol", side_effect=ValueError("Invalid molecule")
        ):
            manager.apply_chem_check_and_set_flags(mol, source_desc="Test")

        assert manager.host.chem_check_tried is True
        assert manager.host.chem_check_failed is True
        assert manager.host.init_manager.optimize_3d_button.setEnabled.called
        manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(
            False
        )

    def test_clear_xyz_flags_with_mol_arg(self, manager):
        """_clear_xyz_flags removes XYZ-specific props and clears is_xyz_derived."""
        mol = Chem.MolFromSmiles("C")
        mol.SetProp("_xyz_skip_checks", "1")
        mol._xyz_skip_checks = True
        mol.xyz_atom_data = {}
        manager.host.is_xyz_derived = True
        manager.host.chem_check_failed = False

        manager._clear_xyz_flags(mol)

        assert not mol.HasProp("_xyz_skip_checks")
        assert not hasattr(mol, "_xyz_skip_checks")
        assert not hasattr(mol, "xyz_atom_data")
        assert manager.host.is_xyz_derived is False
        manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(True)

    def test_update_edit_menu_actions(self, manager):
        """update_edit_menu_actions enables cut/copy/paste when items are selected."""
        manager.host.init_manager.scene.selectedItems.return_value = [MagicMock()]

        mock_clipboard = MagicMock()
        mock_mime = MagicMock()
        mock_mime.hasFormat.return_value = True
        mock_clipboard.mimeData.return_value = mock_mime

        with patch(
            "PyQt6.QtWidgets.QApplication.clipboard", return_value=mock_clipboard
        ):
            manager.update_edit_menu_actions()

        manager.host.init_manager.cut_action.setEnabled.assert_called_with(True)
        manager.host.init_manager.copy_action.setEnabled.assert_called_with(True)
        manager.host.init_manager.paste_action.setEnabled.assert_called_with(True)

    def test_open_rotate_2d_dialog(self, manager):
        """open_rotate_2d_dialog calls rotate_molecule_2d with the user's chosen angle."""
        manager.rotate_molecule_2d = MagicMock()
        mock_dialog = MagicMock()
        mock_dialog.exec.return_value = QDialog.DialogCode.Accepted
        mock_dialog.get_angle.return_value = 45.0

        with patch(
            "moleditpy.ui.edit_actions_logic.Rotate2DDialog", return_value=mock_dialog
        ):
            manager.open_rotate_2d_dialog()

        manager.rotate_molecule_2d.assert_called_with(45.0)
        assert manager.last_rotation_angle == 45.0

    def test_rotate_molecule_2d_full(self, manager):
        """rotate_molecule_2d updates atom positions and stores them in molecule data."""
        manager.host.init_manager.scene.selectedItems.return_value = []
        atom1 = MagicMock(spec=AtomItem)
        atom1.atom_id = 1
        atom1.pos.return_value = QPointF(0, 0)
        manager.host.state_manager.data.atoms = {
            1: {"symbol": "C", "pos": QPointF(0, 0)}
        }
        manager.host.init_manager.scene.atom_items = {1: atom1}

        with patch(
            "moleditpy.core.mol_geometry.rotate_2d_points", return_value={1: (10, 10)}
        ):
            with patch(
                "moleditpy.ui.edit_actions_logic.sip_isdeleted_safe", return_value=False
            ):
                manager.rotate_molecule_2d(90.0)

        atom1.setPos.assert_called()
        manager.host.state_manager.data.set_atom_pos.assert_called()

    def test_select_all(self, manager):
        """select_all selects every AtomItem and BondItem in the scene."""
        atom = MagicMock(spec=AtomItem)
        bond = MagicMock(spec=BondItem)
        manager.host.init_manager.scene.items.return_value = [atom, bond, MagicMock()]

        manager.select_all()

        atom.setSelected.assert_called_with(True)
        bond.setSelected.assert_called_with(True)

    def test_clear_all(self, manager):
        """clear_all(skip_check=True) returns True and shows the cleared status message."""
        manager.host.edit_3d_manager.measurement_mode = True
        manager.host.edit_3d_manager.is_3d_edit_mode = True

        result = manager.clear_all(skip_check=True)

        assert result is True
        manager.host.statusBar().showMessage.assert_called_with("Cleared all data.")

    def test_cut_selection(self, manager):
        """cut_selection copies then deletes the selected items."""
        manager.copy_selection = MagicMock()
        item = MagicMock()
        manager.host.init_manager.scene.selectedItems.return_value = [item]
        manager.host.init_manager.scene.delete_items.return_value = True

        manager.cut_selection()

        manager.copy_selection.assert_called()
        manager.host.init_manager.scene.delete_items.assert_called()
        manager.host.statusBar().showMessage.assert_called_with("Cut selection.", 2000)

    def test_cut_selection_no_selection(self, manager):
        """cut_selection is a no-op when nothing is selected."""
        manager.copy_selection = MagicMock()
        manager.host.init_manager.scene.selectedItems.return_value = []

        manager.cut_selection()

        manager.copy_selection.assert_not_called()

    def test_adjust_molecule_positions_no_collision(self, manager):
        """Fragments already far apart are not moved by collision avoidance."""
        mol = Chem.MolFromSmiles("C.C")
        from rdkit.Chem import AllChem

        AllChem.EmbedMultipleConfs(mol, numConfs=1)
        conf = mol.GetConformer()

        conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        conf.SetAtomPosition(1, (10.0, 0.0, 0.0))

        frags = [[0], [1]]
        pos0_before = list(conf.GetAtomPosition(0))
        pos1_before = list(conf.GetAtomPosition(1))

        manager.adjust_molecule_positions_to_avoid_collisions(mol, frags)

        assert list(conf.GetAtomPosition(0)) == pos0_before
        assert list(conf.GetAtomPosition(1)) == pos1_before

    def test_adjust_molecule_positions_with_collision(self, manager):
        """Colliding fragments are separated by collision avoidance."""
        mol = Chem.MolFromSmiles("C.C")
        from rdkit.Chem import AllChem

        AllChem.EmbedMultipleConfs(mol, numConfs=1)
        conf = mol.GetConformer()

        conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        conf.SetAtomPosition(1, (0.5, 0.0, 0.0))

        frags = [[0], [1]]
        pos0_before = np.array(list(conf.GetAtomPosition(0)))
        pos1_before = np.array(list(conf.GetAtomPosition(1)))

        manager.adjust_molecule_positions_to_avoid_collisions(mol, frags)

        pos0_after = np.array(list(conf.GetAtomPosition(0)))
        pos1_after = np.array(list(conf.GetAtomPosition(1)))

        assert not np.array_equal(pos0_before, pos0_after)
        assert not np.array_equal(pos1_before, pos1_after)
        assert np.linalg.norm(pos0_after - pos1_after) > np.linalg.norm(
            pos0_before - pos1_before
        )

    def test_adjust_molecule_positions_single_fragment(self, manager):
        """A single fragment is not moved by collision avoidance."""
        mol = Chem.MolFromSmiles("C")
        from rdkit.Chem import AllChem

        AllChem.EmbedMolecule(mol)
        conf = mol.GetConformer()

        frags = [[0]]
        pos_before = list(conf.GetAtomPosition(0))

        manager.adjust_molecule_positions_to_avoid_collisions(mol, frags)

        assert list(conf.GetAtomPosition(0)) == pos_before

    def test_apply_chem_check_missing_button(self, manager):
        """apply_chem_check_and_set_flags does not raise when optimize_3d_button is absent."""
        del manager.host.init_manager.optimize_3d_button
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")

        with patch("rdkit.Chem.SanitizeMol", side_effect=ValueError("Invalid module")):
            manager.apply_chem_check_and_set_flags(mol)

    def test_clear_xyz_flags_current_mol(self, manager):
        """_clear_xyz_flags(mol=None) uses the current_mol from the view manager."""
        mol = Chem.MolFromSmiles("C")
        mol.SetProp("_xyz_skip_checks", "1")
        manager.host.view_3d_manager.current_mol = mol
        manager.host.is_xyz_derived = True

        manager._clear_xyz_flags(mol=None)

        assert not mol.HasProp("_xyz_skip_checks")
        assert manager.host.is_xyz_derived is False

    def test_clear_xyz_flags_missing_zoom(self, manager):
        """_clear_xyz_flags does not raise when reset_zoom is absent."""
        del manager.host.view_3d_manager.reset_zoom

        mol = Chem.MolFromSmiles("C")
        manager._clear_xyz_flags(mol)


def test_push_undo_state_caps_stack_depth(monkeypatch):
    """The undo stack drops its oldest entries beyond UNDO_STACK_MAX_DEPTH."""
    from moleditpy.utils.constants import UNDO_STACK_MAX_DEPTH

    host = MagicMock()
    host.is_restoring_state = False
    host.view_3d_manager.current_mol = None
    host.state_manager.data.atoms = {}
    host.state_manager.data.bonds = {}

    mgr = EditActionsManager(host)
    monkeypatch.setattr(mgr, "update_implicit_hydrogens", lambda: None)
    monkeypatch.setattr(mgr, "update_undo_redo_actions", lambda: None)

    for i in range(UNDO_STACK_MAX_DEPTH + 25):
        host.state_manager.data.next_atom_id = i
        host.state_manager.get_current_state.return_value = {
            "atoms": {},
            "bonds": {},
            "_next_atom_id": i,
        }
        mgr.push_undo_state()

    assert len(mgr.undo_stack) == UNDO_STACK_MAX_DEPTH
    # Newest state kept, oldest dropped
    assert mgr.undo_stack[-1]["_next_atom_id"] == UNDO_STACK_MAX_DEPTH + 24
    assert mgr.undo_stack[0]["_next_atom_id"] == 25


def test_push_undo_state_detects_constraint_change(monkeypatch):
    """Two states identical except for constraints_3d must both be pushed.

    Regression: the dedup comparison ignored constraints_3d, so a
    constraints-only edit (e.g. added in the Constrained Optimization
    dialog and saved on close) was silently deduplicated away.
    """
    host = MagicMock()
    host.is_restoring_state = False
    host.view_3d_manager.current_mol = None
    host.state_manager.data.atoms = {}
    host.state_manager.data.bonds = {}
    host.state_manager.data.next_atom_id = 0
    host.edit_3d_manager.constraints_3d = []

    mgr = EditActionsManager(host)
    monkeypatch.setattr(mgr, "update_implicit_hydrogens", lambda: None)
    monkeypatch.setattr(mgr, "update_undo_redo_actions", lambda: None)

    def _state():
        return {
            "atoms": {},
            "bonds": {},
            "_next_atom_id": 0,
            "constraints_3d": [
                list(c) for c in host.edit_3d_manager.constraints_3d
            ],
        }

    host.state_manager.get_current_state.side_effect = _state

    mgr.push_undo_state()
    assert len(mgr.undo_stack) == 1

    # Identical state: deduplicated
    mgr.push_undo_state()
    assert len(mgr.undo_stack) == 1

    # Same atoms/bonds, new constraint: must create a new entry
    host.edit_3d_manager.constraints_3d = [["Distance", [0, 1], 1.54, 1.0e5]]
    mgr.push_undo_state()
    assert len(mgr.undo_stack) == 2
