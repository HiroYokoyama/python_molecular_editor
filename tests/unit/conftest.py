import os
import sys
import pytest
import json
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QPointF, QDateTime, Qt

# Minimal path setup to ensure moleditpy is discoverable
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'moleditpy', 'src'))
if os.path.isdir(src_path) and src_path not in sys.path:
    sys.path.insert(0, src_path)

@pytest.fixture(scope="session")
def app():
    """QApplication session-wide instance for unit tests that might need it."""
    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)
    return q_app

@pytest.fixture
def mock_parser_host(app):
    """A dummy host for MainWindow helper classes to satisfy their interface."""
    from unittest.mock import MagicMock
    from moleditpy.modules.molecular_data import MolecularData
    from moleditpy.modules.atom_item import AtomItem
    
    host = MagicMock()
    # Data and Scene
    host.data = MolecularData()
    host.scene = MagicMock()
    host.scene.data = host.data
    
    # UI and Settings
    host.view_2d = MagicMock()
    host.view_3d = MagicMock()
    host.plotter = MagicMock()
    
    # Robust settings mock (dict-like)
    class MockSettings(dict):
        def value(self, key, default=None):
            return self.get(key, default)
        def setValue(self, key, value):
            self[key] = value
        def set(self, key, value):
            self[key] = value
            
    host.settings = MockSettings({
        'skip_chemistry_checks': True, 
        'optimization_method': 'UFF_RDKIT',
        'conversion_mode': 'rdkit'
    })

    host.is_2d_editable = True
    host.current_mol = None
    host.current_file_path = None
    host.has_unsaved_changes = False
    host.constraints_3d = []
    host.is_xyz_derived = False
    host.active_worker_ids = set()
    host.halt_ids = set()
    host._ih_update_counter = 0

    # Status Bar
    host.statusBar_mock = MagicMock()
    host.statusBar.return_value = host.statusBar_mock
    
    # Simplified Logic methods for IO (avoiding MainWindowAppState dependencies)
    def create_json_data():
        json_data = {
            "format": "PME Project",
            "version": "1.0",
            "2d_structure": {
                "atoms": [],
                "bonds": [],
                "next_atom_id": host.data._next_atom_id
            }
        }
        for aid, data in host.data.atoms.items():
            pos = data['item'].pos()
            json_data["2d_structure"]["atoms"].append({
                "id": aid, "symbol": data['symbol'], 
                "x": pos.x(), "y": pos.y(),
                "charge": data.get('charge', 0),
                "radical": data.get('radical', 0)
            })
        for (id1, id2), bdata in host.data.bonds.items():
            json_data["2d_structure"]["bonds"].append({
                "atom1": id1, "atom2": id2, "order": bdata['order'], "stereo": bdata.get('stereo', 0)
            })
        return json_data

    def load_from_json_data(data):
        host.data.atoms.clear()
        host.data.bonds.clear()
        s2d = data.get("2d_structure", {})
        for adata in s2d.get("atoms", []):
            host.scene.create_atom(adata["symbol"], QPointF(adata["x"], adata["y"]), 
                                 charge=adata.get("charge", 0), radical=adata.get("radical", 0))
        for bdata in s2d.get("bonds", []):
            a1 = host.data.atoms[bdata["atom1"]]['item']
            a2 = host.data.atoms[bdata["atom2"]]['item']
            host.scene.create_bond(a1, a2, bond_order=bdata["order"], bond_stereo=bdata.get("stereo", 0))
        host.data._next_atom_id = s2d.get("next_atom_id", 0)

    host.create_json_data.side_effect = create_json_data
    host.load_from_json_data.side_effect = load_from_json_data
    
    host.push_undo_state = MagicMock()
    host.reset_undo_stack = MagicMock()
    host.update_window_title = MagicMock()
    host.clear_2d_editor = MagicMock()
    host.restore_ui_for_editing = MagicMock()
    host.activate_select_mode = MagicMock()
    host._enable_3d_features = MagicMock()
    host._enable_3d_edit_actions = MagicMock()
    host.check_unsaved_changes = MagicMock(return_value=True)

    # Scene helpers
    host.scene.find_bond_between.return_value = None

    def default_create_atom(symbol, pos, charge=0, radical=0):
        aid = host.data.add_atom(symbol, pos, charge=charge, radical=radical)
        item = MagicMock(spec=AtomItem)
        item.pos.return_value = pos
        item.atom_id = aid
        item.symbol = symbol
        item.charge = charge
        item.radical = radical
        item.has_problem = False
        item.__class__ = AtomItem
        item.scene.return_value = host.scene
        host.data.atoms[aid]['item'] = item
        return aid

    def default_create_bond(a1_item, a2_item, bond_order=1, bond_stereo=0):
        id1 = a1_item.atom_id if hasattr(a1_item, 'atom_id') else a1_item
        id2 = a2_item.atom_id if hasattr(a2_item, 'atom_id') else a2_item
        return host.data.add_bond(id1, id2, order=bond_order, stereo=bond_stereo)

    host.scene.create_atom.side_effect = default_create_atom
    host.scene.create_bond.side_effect = default_create_bond
    
    return host
