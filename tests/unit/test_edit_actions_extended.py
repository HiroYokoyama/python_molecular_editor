import os
import sys
import pytest
import math
import numpy as np
from unittest.mock import MagicMock, patch
from PyQt6.QtWidgets import QApplication, QDialog
from PyQt6.QtCore import QPointF

# Ensure local moleditpy is discoverable
workspace_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(workspace_src_path) and workspace_src_path not in sys.path:
    sys.path.insert(0, workspace_src_path)

from moleditpy.ui.edit_actions_logic import EditActionsManager
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem
from rdkit import Chem

class DummyHost:
    def __init__(self):
        self.statusBar_mock = MagicMock()
        self.settings = {}
        self.is_xyz_derived = False
        self.chem_check_tried = False
        self.chem_check_failed = False
        self.initialization_complete = True
        self._is_restoring_state = False
        
        # Child Managers
        self.init_manager = MagicMock()
        self.init_manager.settings = self.settings
        self.state_manager = MagicMock()
        self.ui_manager = MagicMock()
        self.view_3d_manager = MagicMock()
        self.edit_3d_manager = MagicMock()
        self.plugin_manager = MagicMock()
        self.edit_actions_manager = None # Set by fixture
        
        # Shortcuts for frequently accessed items on managers
        self.init_manager.scene = MagicMock()
        self.state_manager.data = MagicMock()
        self.init_manager.view_2d = MagicMock()
        self.view_3d_manager.plotter = MagicMock()
        
        # Action mocks on init_manager
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

class TestEditActionsExtended:
    @pytest.fixture
    def manager(self):
        host = DummyHost()
        mgr = EditActionsManager(host)
        host.edit_actions_manager = mgr
        
        # CLIPBOARD_MIME_TYPE must be available
        if not hasattr(sys.modules[mgr.__module__], "CLIPBOARD_MIME_TYPE"):
            setattr(sys.modules[mgr.__module__], "CLIPBOARD_MIME_TYPE", "application/x-moleditpy-pme")
        
        return mgr

    def test_apply_chem_check_force_skip(self, manager):
        mol = Chem.MolFromSmiles("C")
        manager._apply_chem_check_and_set_flags(mol, force_skip=True)
        assert manager.host.chem_check_tried is False
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_settings_skip(self, manager):
        manager.host.init_manager.settings["skip_chemistry_checks"] = True
        mol = Chem.MolFromSmiles("C")
        manager._apply_chem_check_and_set_flags(mol)
        assert manager.host.chem_check_tried is False
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_success(self, manager):
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")
        manager._apply_chem_check_and_set_flags(mol)
        assert manager.host.chem_check_tried is True
        assert manager.host.chem_check_failed is False

    def test_apply_chem_check_failure(self, manager):
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")
        with patch("rdkit.Chem.SanitizeMol", side_effect=ValueError("Invalid molecule")):
            manager._apply_chem_check_and_set_flags(mol, source_desc="Test")
        
        assert manager.host.chem_check_tried is True
        assert manager.host.chem_check_failed is True
        # Verify optimize_3d_button is disabled
        assert manager.host.init_manager.optimize_3d_button.setEnabled.called
        manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(False)

    def test_clear_xyz_flags_with_mol_arg(self, manager):
        mol = Chem.MolFromSmiles("C")
        mol.SetProp("_xyz_skip_checks", "1")
        mol._xyz_skip_checks = True
        mol._xyz_atom_data = {}
        manager.host.is_xyz_derived = True
        manager.host.chem_check_failed = False
        
        manager._clear_xyz_flags(mol)
        
        assert not mol.HasProp("_xyz_skip_checks")
        assert not hasattr(mol, "_xyz_skip_checks")
        assert manager.host.is_xyz_derived is False
        manager.host.init_manager.optimize_3d_button.setEnabled.assert_called_with(True)

    def test_update_edit_menu_actions(self, manager):
        manager.host.init_manager.scene.selectedItems.return_value = [MagicMock()]
        
        mock_clipboard = MagicMock()
        mock_mime = MagicMock()
        mock_mime.hasFormat.return_value = True
        mock_clipboard.mimeData.return_value = mock_mime
        
        with patch("PyQt6.QtWidgets.QApplication.clipboard", return_value=mock_clipboard):
            manager.update_edit_menu_actions()
            
        manager.host.init_manager.cut_action.setEnabled.assert_called_with(True)
        manager.host.init_manager.copy_action.setEnabled.assert_called_with(True)
        manager.host.init_manager.paste_action.setEnabled.assert_called_with(True)

    def test_open_rotate_2d_dialog(self, manager):
        manager.rotate_molecule_2d = MagicMock()
        mock_dialog = MagicMock()
        mock_dialog.exec.return_value = QDialog.DialogCode.Accepted
        mock_dialog.get_angle.return_value = 45.0
        
        with patch("moleditpy.ui.edit_actions_logic.Rotate2DDialog", return_value=mock_dialog):
            manager.open_rotate_2d_dialog()
            
        manager.rotate_molecule_2d.assert_called_with(45.0)
        assert manager.last_rotation_angle == 45.0

    def test_rotate_molecule_2d_full(self, manager):
        manager.host.init_manager.scene.selectedItems.return_value = []
        atom1 = MagicMock(spec=AtomItem)
        atom1.atom_id = 1
        atom1.pos.return_value = QPointF(0, 0)
        manager.host.state_manager.data.atoms = {1: {"item": atom1, "symbol": "C", "pos": QPointF(0, 0)}}
        
        with patch("moleditpy.core.mol_geometry.rotate_2d_points", return_value={1: (10, 10)}):
            with patch("moleditpy.ui.edit_actions_logic.sip_isdeleted_safe", return_value=False):
                manager.rotate_molecule_2d(90.0)
                
        atom1.setPos.assert_called()
        manager.host.state_manager.data.set_atom_pos.assert_called()

    def test_select_all(self, manager):
        atom = MagicMock(spec=AtomItem)
        bond = MagicMock(spec=BondItem)
        manager.host.init_manager.scene.items.return_value = [atom, bond, MagicMock()]
        
        manager.select_all()
        
        atom.setSelected.assert_called_with(True)
        bond.setSelected.assert_called_with(True)

    def test_clear_all(self, manager):
        manager.host.edit_3d_manager.measurement_mode = True
        manager.host.edit_3d_manager.is_3d_edit_mode = True
        
        result = manager.clear_all(skip_check=True)
        
        assert result is True
        manager.host.statusBar().showMessage.assert_called_with("Cleared all data.")

    def test_cut_selection(self, manager):
        manager.copy_selection = MagicMock()
        item = MagicMock()
        manager.host.init_manager.scene.selectedItems.return_value = [item]
        manager.host.init_manager.scene.delete_items.return_value = True
        
        manager.cut_selection()
        
        manager.copy_selection.assert_called()
        manager.host.init_manager.scene.delete_items.assert_called()
        manager.host.statusBar().showMessage.assert_called_with("Cut selection.", 2000)

    def test_cut_selection_no_selection(self, manager):
        manager.copy_selection = MagicMock()
        manager.host.init_manager.scene.selectedItems.return_value = []
        
        manager.cut_selection()
        
        manager.copy_selection.assert_not_called()

    def test_adjust_molecule_positions_no_collision(self, manager):
        mol = Chem.MolFromSmiles("C.C")
        from rdkit.Chem import AllChem
        AllChem.EmbedMultipleConfs(mol, numConfs=1)
        conf = mol.GetConformer()
        
        # Set them far apart
        conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        conf.SetAtomPosition(1, (10.0, 0.0, 0.0))
        
        frags = [[0], [1]]
        pos0_before = list(conf.GetAtomPosition(0))
        pos1_before = list(conf.GetAtomPosition(1))
        
        manager.adjust_molecule_positions_to_avoid_collisions(mol, frags)
        
        assert list(conf.GetAtomPosition(0)) == pos0_before
        assert list(conf.GetAtomPosition(1)) == pos1_before

    def test_adjust_molecule_positions_with_collision(self, manager):
        mol = Chem.MolFromSmiles("C.C")
        from rdkit.Chem import AllChem
        AllChem.EmbedMultipleConfs(mol, numConfs=1)
        conf = mol.GetConformer()
        
        # Set them very close (colliding)
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
        assert np.linalg.norm(pos0_after - pos1_after) > np.linalg.norm(pos0_before - pos1_before)

    def test_adjust_molecule_positions_single_fragment(self, manager):
        mol = Chem.MolFromSmiles("C")
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol)
        conf = mol.GetConformer()
        
        frags = [[0]]
        pos_before = list(conf.GetAtomPosition(0))
        
        manager.adjust_molecule_positions_to_avoid_collisions(mol, frags)
        
        assert list(conf.GetAtomPosition(0)) == pos_before

    def test_apply_chem_check_missing_button(self, manager):
        # Remove optimize_3d_button from mock
        del manager.host.init_manager.optimize_3d_button
        manager.host.init_manager.settings["skip_chemistry_checks"] = False
        mol = Chem.MolFromSmiles("C")
        
        # Should log error but not crash when sanitization fails
        with patch("rdkit.Chem.SanitizeMol", side_effect=ValueError("Invalid module")):
            with patch("moleditpy.ui.edit_actions_logic.logging.error") as mock_log:
                manager._apply_chem_check_and_set_flags(mol)
                mock_log.assert_called()

    def test_clear_xyz_flags_current_mol(self, manager):
        mol = Chem.MolFromSmiles("C")
        mol.SetProp("_xyz_skip_checks", "1")
        # Ensure view_3d_manager is mocked and current_mol is set
        manager.host.view_3d_manager.current_mol = mol
        manager.host.is_xyz_derived = True
        
        manager._clear_xyz_flags(mol=None) # Use current_mol
        
        assert not mol.HasProp("_xyz_skip_checks")
        assert manager.host.is_xyz_derived is False

    def test_clear_xyz_flags_missing_zoom(self, manager):
        # Remove reset_zoom from view_3d_manager mock
        del manager.host.view_3d_manager.reset_zoom
        
        mol = Chem.MolFromSmiles("C")
        with patch("moleditpy.ui.edit_actions_logic.logging.error") as mock_log:
            manager._clear_xyz_flags(mol)
            mock_log.assert_called()
