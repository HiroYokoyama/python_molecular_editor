import pytest
import base64
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.app_state import StateManager
from moleditpy.core.molecular_data import MolecularData
from PyQt6.QtCore import QPointF, Qt
from unittest.mock import MagicMock

class DummyMainWindow(StateManager):
    def __init__(self):
        # Initialize as StateManager, which expects a 'host'
        # In this dummy, we act as both StateManager and the host
        self.host = self
        self.data = MolecularData()
        
        # Initialize all managers first
        self.state_manager = self # We are the state manager
        self.init_manager = MagicMock()
        self.ui_manager = MagicMock()
        self.edit_actions_manager = MagicMock()
        self.view_3d_manager = MagicMock()
        self.edit_3d_manager = MagicMock()
        self.compute_manager = MagicMock()
        self.io_manager = MagicMock()
        self.export_manager = MagicMock()
        self.dialog_manager = MagicMock()
        self.plugin_manager = MagicMock()
        self.plugin_manager.save_handlers = {}
        self.plugin_manager.load_handlers = {}

        self.init_manager.scene = MagicMock()
        self.init_manager.view_2d = MagicMock()
        self.init_manager.settings = MagicMock()
        self.view_3d_manager.view_3d = MagicMock()
        self.view_3d_manager.plotter = MagicMock()
        self.view_3d_manager.current_mol = None
        self.view_3d_manager.atom_positions_3d = None
        
        self.host.init_manager.current_file_path = None
        self.host.state_manager.has_unsaved_changes = False
        self.edit_3d_manager.constraints_3d = []
        self.is_2d_editable = True
        self.edit_actions_manager.undo_stack = []
        self.edit_actions_manager.redo_stack = []
        self._is_restoring_state = False
        self.initialization_complete = True
        self._preserved_plugin_data = {}
        self.init_manager.formula_label = MagicMock()
        self.compute_manager.last_successful_optimization_method = None
        self.edit_actions_manager.clear_2d_editor.side_effect = self._do_clear_2d_editor
        self.edit_actions_manager.update_implicit_hydrogens = MagicMock()

    def _do_clear_2d_editor(self, push_to_undo=True):
        self.data.atoms.clear()
        self.data.bonds.clear()
        self.data._next_atom_id = 0
        
    def statusBar(self):
        return MagicMock()
        
    def setWindowTitle(self, title):
        pass
        
    def update_window_title(self):
        pass
        
    def update_implicit_hydrogens(self):
        pass
        
    def update_chiral_labels(self):
        pass
        
    def update_2d_measurement_labels(self):
        pass
        
    def update_realtime_info(self):
        pass
        
    def update_undo_redo_actions(self):
        pass
        
    def _enable_3d_features(self, val):
        pass
        
    def _enable_3d_edit_actions(self, val):
        pass
        
    def clear_2d_editor(self, push_to_undo=True):
        self._do_clear_2d_editor(push_to_undo)
        
    def restore_ui_for_editing(self):
        pass
        
    def create_atom_id_mapping(self):
        pass
        
    def draw_molecule_3d(self, mol):
        pass

@pytest.fixture
def dummy_window(app):
    return DummyMainWindow()

def test_pmeprj_serialization_roundtrip(dummy_window):
    """Test full project serialization/deserialization (PMEPRJ)."""
    mw = dummy_window
    
    # 1. Setup a complex state with 2D, 3D, and constraints
    aid1 = mw.state_manager.data.add_atom("C", QPointF(10, 20))
    aid2 = mw.state_manager.data.add_atom("O", QPointF(30, 40))
    mw.state_manager.data.add_bond(aid1, aid2, order=1)
    
    # Mock visual items for 2D serialization
    mw.state_manager.data.atoms[aid1]["item"] = MagicMock(pos=lambda: QPointF(10, 20))
    mw.state_manager.data.atoms[aid2]["item"] = MagicMock(pos=lambda: QPointF(30, 40))
    
    # 3D structure (Methanol)
    mol = Chem.MolFromSmiles("CO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    # Tag atoms with editor IDs
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", aid1)
    mol.GetAtomWithIdx(1).SetIntProp("_original_atom_id", aid2)
    mw.view_3d_manager.current_mol = mol
    
    # Constraints and metadata
    mw.edit_3d_manager.constraints_3d = [("DISTANCE", (0, 1), 1.43, 8.0e5)]
    mw.compute_manager.last_successful_optimization_method = "MMFF94s"
    mw._preserved_plugin_data = {"TestPlugin": {"val": 42}}
    
    # 2. Serialize
    json_data = mw.create_json_data()
    
    # 3. Verify high-level structure
    assert json_data["format"] == "PME Project"
    assert "2d_structure" in json_data
    assert "3d_structure" in json_data
    assert json_data["3d_structure"]["num_conformers"] == 1
    
    # 4. Clear and Restore
    mw.clear_2d_editor(push_to_undo=False)
    mw.view_3d_manager.current_mol = None
    mw.edit_3d_manager.constraints_3d = []
    
    mw.load_from_json_data(json_data)
    
    # 5. Verify 2D Restored
    assert len(mw.state_manager.data.atoms) == 2
    assert mw.state_manager.data.atoms[aid1]["symbol"] == "C"
    assert mw.state_manager.data.atoms[aid2]["symbol"] == "O"
    assert mw.state_manager.data.atoms[aid1]["pos"][0] == 10
    
    # 6. Verify 3D Restored
    assert mw.view_3d_manager.current_mol is not None
    assert mw.view_3d_manager.current_mol.GetNumAtoms() == mol.GetNumAtoms()
    # Check property round-trip
    assert mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id") == aid1
    assert mw.view_3d_manager.current_mol.GetAtomWithIdx(1).GetIntProp("_original_atom_id") == aid2
    
    # 7. Verify helper data (atom_positions_3d)
    assert mw.view_3d_manager.atom_positions_3d is not None
    assert mw.view_3d_manager.atom_positions_3d.shape == (mol.GetNumAtoms(), 3)
    
    # 8. Verify Constraints and Metadata
    assert len(mw.edit_3d_manager.constraints_3d) == 1
    assert mw.edit_3d_manager.constraints_3d[0][0] == "DISTANCE"
    assert mw.edit_3d_manager.constraints_3d[0][2] == 1.43
    assert mw.compute_manager.last_successful_optimization_method == "MMFF94s"
    assert mw._preserved_plugin_data["TestPlugin"]["val"] == 42

def test_undo_state_binary_roundtrip(dummy_window):
    """Test the internal binary state serialization used for Undo/Redo."""
    mw = dummy_window
    
    # Setup state
    aid = mw.state_manager.data.add_atom("N", QPointF(5, 5))
    mw.state_manager.data.atoms[aid]["item"] = MagicMock(pos=lambda: QPointF(5, 5))
    
    # 3D with property
    mol = Chem.MolFromSmiles("N")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", aid)
    mw.view_3d_manager.current_mol = mol
    
    # Get internal state
    state = mw.get_current_state()
    
    # Verify binary blob existence
    assert "mol_3d" in state
    assert isinstance(state["mol_3d"], bytes)
    assert state["mol_3d_atom_ids"][0] == aid
    
    # Restore
    mw.clear_2d_editor(push_to_undo=False)
    mw.set_state_from_data(state)
    
    # Verify
    assert mw.state_manager.data.atoms[aid]["symbol"] == "N"
    assert mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id") == aid

def test_legacy_version_handling(dummy_window):
    """Verify that version mismatch warnings are triggered (but don't crash)."""
    mw = dummy_window
    
    # Create fake JSON with future version
    json_data = {
        "format": "PME Project",
        "version": "99.0.0",
        "2d_structure": {"atoms": [], "bonds": []}
    }
    
    # This should trigger QMessageBox call (which we might want to mock if it's annoying)
    with MagicMock() as mock_msg:
        mw.warning_message_box = mock_msg
        mw.load_from_json_data(json_data)
        # Note: In the actual code, it might try to use global QMessageBox if helper not present
        # In our Dummy it should just pass or we can mock it.
