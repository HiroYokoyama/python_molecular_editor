import os
import sys
import pytest
from PyQt6.QtWidgets import QApplication

# Minimal path setup to ensure moleditpy is discoverable
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'moleditpy', 'src'))
if os.path.isdir(src_path) and src_path not in sys.path:
    sys.path.insert(0, src_path)

@pytest.fixture(scope="session")
def app():
    """QApplication session-wide instance for unit tests that might need it (e.g. for parsers using QFileDialog mocks)."""
    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)
    return q_app

@pytest.fixture
def mock_parser_host(app):
    """A dummy host for MainWindow helper classes to satisfy their interface."""
    from unittest.mock import MagicMock
    from moleditpy.modules.molecular_data import MolecularData
    host = MagicMock()
    # Data and Scene
    host.data = MolecularData()
    host.scene = MagicMock()
    host.scene.data = host.data
    
    # UI and Settings
    host.view_2d = MagicMock()
    host.plotter = MagicMock()
    host.settings = {'skip_chemistry_checks': True}
    host.is_2d_editable = True
    host.current_mol = None
    host.current_file_path = None
    host.has_unsaved_changes = False
    host.constraints_3d = []
    host.is_2d_editable = True
    host.is_xyz_derived = False
    host.last_successful_optimization_method = None
    host._preserved_plugin_data = {}
    host.plugin_manager = MagicMock()
    host.plugin_manager.save_handlers = {}
    host.initialization_complete = True
    host.formula_label = MagicMock()
    host.statusBar = MagicMock()
    host.statusBar.return_value = MagicMock()
    
    # Logic methods typically on MainWindow (AppState)
    from moleditpy.modules.main_window_app_state import MainWindowAppState
    host.create_json_data.side_effect = lambda: MainWindowAppState.create_json_data(host)
    host.load_from_json_data.side_effect = lambda data: MainWindowAppState.load_from_json_data(host, data)
    
    host.push_undo_state = MagicMock()
    host.reset_undo_stack = MagicMock()
    host.update_window_title = MagicMock()
    host.clear_2d_editor = MagicMock()
    host.fit_to_view = MagicMock()
    host.activate_select_mode = MagicMock()
    host.analysis_action = MagicMock()
    host.optimize_3d_button = MagicMock()
    host._enable_3d_edit_actions = MagicMock()
    host._enable_3d_features = MagicMock()
    host.update_implicit_hydrogens = MagicMock()
    host.update_chiral_labels = MagicMock()
    host.restore_ui_for_editing = MagicMock()
    host.update_2d_measurement_labels = MagicMock()
    host._is_restoring_state = False
    host.optimize_3d_button = MagicMock()

    # Ensure scene.create_atom/bond actually update host.data if not mocked per-test
    def default_create_atom(symbol, pos, charge=0, radical=0):
        from moleditpy.modules.atom_item import AtomItem
        aid = host.data.add_atom(symbol, pos, charge=charge, radical=radical)
        # We MUST set 'item' because app logic often accesses host.data.atoms[aid]['item']
        # Using spec=AtomItem to pass isinstance(item, AtomItem) checks
        item = MagicMock(spec=AtomItem)
        item.pos.return_value = pos
        item.atom_id = aid  # Expected by some logic
        host.data.atoms[aid]['item'] = item
        return aid

    def default_create_bond(a1_item, a2_item, bond_order=1, bond_stereo=0):
        # Extract IDs if items are passed (app logic style)
        id1 = a1_item.atom_id if hasattr(a1_item, 'atom_id') else a1_item
        id2 = a2_item.atom_id if hasattr(a2_item, 'atom_id') else a2_item
        return host.data.add_bond(id1, id2, order=bond_order, stereo=bond_stereo)

    host.scene.create_atom.side_effect = default_create_atom
    host.scene.create_bond.side_effect = default_create_bond
    
    # Debug helper: print status bar messages to console
    def print_status(msg, timeout=0):
        print(f"STATUS: {msg}", flush=True)
    host.statusBar().showMessage.side_effect = print_status

    return host
