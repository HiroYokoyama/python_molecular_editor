import os
import sys
import pytest
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QPointF

# Minimal path setup to ensure moleditpy is discoverable
src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(src_path) and src_path not in sys.path:
    sys.path.insert(0, src_path)


@pytest.fixture(scope="session")
def app():
    """QApplication session-wide instance with platform-aware teardown."""
    is_headless = os.environ.get("MOLEDITPY_HEADLESS", "0") == "1"
    is_offscreen = os.environ.get("QT_QPA_PLATFORM") == "offscreen"

    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)

    yield q_app

    # Platform-aware teardown: prevent 0x80010108 on Windows but avoid CI segfaults
    if not (is_headless or is_offscreen):
        try:
            # Aggressively neutralize blocking dialogs/filters before bulk closure
            # This prevents 0xc0000374 (heap corruption) on Windows during rapid teardown.
            for widget in QApplication.topLevelWidgets():
                try:
                    # 1. Force has_unsaved_changes to False to skip "Save?" prompts
                    if hasattr(widget, "state_manager"):
                        try:
                            widget.state_manager.has_unsaved_changes = False
                        except Exception:
                            pass
                    elif hasattr(widget, "has_unsaved_changes"):
                        widget.has_unsaved_changes = False

                    # 2. Detach UIManager filters to prevent recursive/blocking logic
                    if hasattr(widget, "ui_manager"):
                        try:
                            widget.removeEventFilter(widget.ui_manager)
                        except Exception:
                            pass

                    # 3. Fast-path closeEvent
                    widget.closeEvent = lambda event: event.accept()
                except Exception:
                    pass

            q_app.closeAllWindows()
            for _ in range(10):
                q_app.processEvents()

            import gc

            gc.collect()
            q_app.processEvents()

            try:
                import colorama

                colorama.deinit()
            except (ImportError, Exception):
                pass
        except Exception:
            pass

    q_app.quit()


@pytest.fixture
def mock_parser_host(app):
    """A dummy host for MainWindow helper classes to satisfy their interface."""
    from unittest.mock import MagicMock
    from moleditpy.core.molecular_data import MolecularData
    from moleditpy.ui.atom_item import AtomItem

    host = MagicMock()

    # 0. Initialize Managers first to avoid overwriting attributes later
    host.state_manager = MagicMock()
    host.init_manager = MagicMock()
    host.ui_manager = MagicMock()
    host.edit_actions_manager = MagicMock()
    host.view_3d_manager = MagicMock()
    host.edit_3d_manager = MagicMock()
    host.compute_manager = MagicMock()
    host.io_manager = MagicMock()
    host.export_manager = MagicMock()
    host.dialog_manager = MagicMock()
    host.plugin_manager = MagicMock()

    # 1. Data and Scene
    host.state_manager.data = MolecularData()
    host.init_manager.scene = MagicMock()
    # Ensure scene.data is the same object
    host.init_manager.scene.data = host.state_manager.data

    # 2. UI and Settings
    host.init_manager.view_2d = MagicMock()
    host.view_3d_manager.view_3d = MagicMock()
    host.view_3d_manager.plotter = MagicMock()
    host.view_3d_manager.current_mol = None

    # Robust settings mock (dict-like)
    class MockSettings(dict):
        def value(self, key, default=None):
            return self.get(key, default)

        def setValue(self, key, value):
            self[key] = value

        def set(self, key, value):
            self[key] = value

    host.init_manager.settings = MockSettings(
        {
            "skip_chemistry_checks": True,
            "optimization_method": "UFF_RDKIT",
            "conversion_mode": "rdkit",
        }
    )

    # 3. Manager proxies on host for convenience
    # These are allowed in test files to keep tests readable and backward compatible.
    type(host).settings = property(lambda self: self.init_manager.settings)
    type(host).data = property(lambda self: self.state_manager.data)
    type(host).scene = property(lambda self: self.init_manager.scene)
    type(host).view_2d = property(lambda self: self.init_manager.view_2d)
    type(host).view_3d = property(lambda self: self.view_3d_manager.view_3d)
    type(host).plotter = property(lambda self: self.view_3d_manager.plotter)

    # 4. State variables
    host.ui_manager.is_2d_editable = True
    host.init_manager.current_file_path = None
    host.state_manager.has_unsaved_changes = False
    host.edit_3d_manager.constraints_3d = []
    host.is_xyz_derived = False
    host.compute_manager.active_worker_ids = set()
    host.compute_manager.halt_ids = set()
    host.ih_update_counter = 0

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
                "next_atom_id": host.state_manager.data.next_atom_id,
            },
        }
        for aid, data in host.state_manager.data.atoms.items():
            item = host.init_manager.scene.atom_items.get(aid)
            pos = item.pos() if item else QPointF(0, 0)
            json_data["2d_structure"]["atoms"].append(
                {
                    "id": aid,
                    "symbol": data["symbol"],
                    "x": pos.x(),
                    "y": pos.y(),
                    "charge": data.get("charge", 0),
                    "radical": data.get("radical", 0),
                }
            )
        for (id1, id2), bdata in host.state_manager.data.bonds.items():
            json_data["2d_structure"]["bonds"].append(
                {
                    "atom1": id1,
                    "atom2": id2,
                    "order": bdata["order"],
                    "stereo": bdata.get("stereo", 0),
                }
            )
        return json_data

    def load_from_json_data(data):
        host.state_manager.data.atoms.clear()
        host.state_manager.data.bonds.clear()
        s2d = data.get("2d_structure", {})
        for adata in s2d.get("atoms", []):
            host.init_manager.scene.create_atom(
                adata["symbol"],
                QPointF(adata["x"], adata["y"]),
                charge=adata.get("charge", 0),
                radical=adata.get("radical", 0),
            )
        for bdata in s2d.get("bonds", []):
            a1 = host.init_manager.scene.atom_items[bdata["atom1"]]
            a2 = host.init_manager.scene.atom_items[bdata["atom2"]]
            host.init_manager.scene.create_bond(
                a1, a2, bond_order=bdata["order"], bond_stereo=bdata.get("stereo", 0)
            )
        host.state_manager.data._next_atom_id = s2d.get("next_atom_id", 0)

    # Apply side effects to both host and state_manager for maximum compatibility
    host.create_json_data.side_effect = create_json_data
    host.load_from_json_data.side_effect = load_from_json_data
    host.state_manager.create_json_data.side_effect = create_json_data
    host.state_manager.load_from_json_data.side_effect = load_from_json_data

    # Common mock behaviors
    host.state_manager.check_unsaved_changes.return_value = True
    host.state_manager.update_window_title = MagicMock()
    host.edit_actions_manager.push_undo_state = MagicMock()
    host.edit_actions_manager.clear_2d_editor = MagicMock()
    host.edit_actions_manager.reset_undo_stack = MagicMock()
    host.ui_manager.restore_ui_for_editing = MagicMock()
    host.ui_manager.activate_select_mode = MagicMock()
    host.ui_manager._enable_3d_features = MagicMock()
    host.ui_manager._enable_3d_edit_actions = MagicMock()
    host.check_unsaved_changes.return_value = True

    # Scene helpers
    host.init_manager.scene.atom_items = {}
    host.init_manager.scene.bond_items = {}
    host.init_manager.scene.find_bond_between.return_value = None

    def default_create_atom(symbol, pos, charge=0, radical=0):
        aid = host.state_manager.data.add_atom(
            symbol, pos, charge=charge, radical=radical
        )
        item = MagicMock(spec=AtomItem)
        item.pos.return_value = pos
        item.atom_id = aid
        item.symbol = symbol
        item.charge = charge
        item.radical = radical
        item.has_problem = False
        item.__class__ = AtomItem
        item.scene.return_value = host.init_manager.scene
        host.init_manager.scene.atom_items[aid] = item
        return aid

    def default_create_bond(a1_item, a2_item, bond_order=1, bond_stereo=0):
        id1 = a1_item.atom_id if hasattr(a1_item, "atom_id") else a1_item
        id2 = a2_item.atom_id if hasattr(a2_item, "atom_id") else a2_item
        return host.state_manager.data.add_bond(
            id1, id2, order=bond_order, stereo=bond_stereo
        )

    # 5. Mediator methods
    def set_current_molecule(mol):
        host.view_3d_manager.current_mol = mol

    def set_3d_atom_positions(positions):
        host.view_3d_manager.atom_positions_3d = positions

    def clear_3d_view():
        host.view_3d_manager.current_mol = None
        if host.view_3d_manager.plotter:
            host.view_3d_manager.plotter.clear()

    def set_plotter_camera_position(camera_position):
        if host.view_3d_manager.plotter:
            host.view_3d_manager.plotter.camera_position = camera_position

    def set_plotter_picker(picker):
        if host.view_3d_manager.plotter:
            host.view_3d_manager.plotter.picker = picker

    def set_constraints_3d(constraints):
        host.edit_3d_manager.constraints_3d = constraints

    def get_constraints_3d():
        return host.edit_3d_manager.constraints_3d

    def set_has_unsaved_changes(value):
        host.state_manager.has_unsaved_changes = value

    def set_current_file_path(path):
        host.init_manager.current_file_path = path

    def get_current_file_path():
        return host.init_manager.current_file_path

    def get_molecule_data():
        return host.state_manager.data

    def set_molecule_data(data):
        host.state_manager.data = data
        if host.init_manager.scene:
            host.init_manager.scene.data = data

    def set_settings_dirty(value):
        host.init_manager.settings_dirty = value

    def set_is_2d_editable(value):
        host.ui_manager.is_2d_editable = value

    def set_optimization_method(method):
        host.init_manager.optimization_method = method

    def set_3d_edit_mode(enabled):
        host.edit_3d_manager.is_3d_edit_mode = enabled

    def is_3d_measurement_mode():
        return bool(host.edit_3d_manager.measurement_mode)

    def set_scene_mode(mode):
        if host.init_manager.scene:
            host.init_manager.scene.mode = mode

    def set_scene_atom_symbol(symbol):
        if host.init_manager.scene:
            host.init_manager.scene.current_atom_symbol = symbol

    def set_scene_bond_properties(order, stereo=0):
        if host.init_manager.scene:
            host.init_manager.scene.bond_order = order
            host.init_manager.scene.bond_stereo = stereo

    def set_scene_user_template_data(data):
        if host.init_manager.scene:
            host.init_manager.scene.user_template_data = data

    def update_status_message(message, timeout=0):
        if host.statusBar():
            if timeout == 0:
                host.statusBar().showMessage(message)
            else:
                host.statusBar().showMessage(message, timeout)

    def get_settings():
        return host.init_manager.settings

    host.set_current_molecule.side_effect = set_current_molecule
    host.set_3d_atom_positions.side_effect = set_3d_atom_positions
    host.clear_3d_view.side_effect = clear_3d_view
    host.set_plotter_camera_position.side_effect = set_plotter_camera_position
    host.set_plotter_picker.side_effect = set_plotter_picker
    host.set_constraints_3d.side_effect = set_constraints_3d
    host.get_constraints_3d.side_effect = get_constraints_3d
    host.set_has_unsaved_changes.side_effect = set_has_unsaved_changes
    host.set_current_file_path.side_effect = set_current_file_path
    host.get_current_file_path.side_effect = get_current_file_path
    host.get_molecule_data.side_effect = get_molecule_data
    host.set_molecule_data.side_effect = set_molecule_data
    host.set_settings_dirty.side_effect = set_settings_dirty
    host.set_is_2d_editable.side_effect = set_is_2d_editable
    host.set_optimization_method.side_effect = set_optimization_method
    host.set_3d_edit_mode.side_effect = set_3d_edit_mode
    host.is_3d_measurement_mode.side_effect = is_3d_measurement_mode
    host.set_scene_mode.side_effect = set_scene_mode
    host.set_scene_atom_symbol.side_effect = set_scene_atom_symbol
    host.set_scene_bond_properties.side_effect = set_scene_bond_properties
    host.set_scene_user_template_data.side_effect = set_scene_user_template_data
    host.update_status_message.side_effect = update_status_message
    host.get_settings.side_effect = get_settings
    host.init_manager.scene.create_atom.side_effect = default_create_atom
    host.init_manager.scene.create_bond.side_effect = default_create_bond

    return host
