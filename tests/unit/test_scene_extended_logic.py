"""Unit tests for MoleculeScene extended atom and bond logic."""

import pytest
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QTransform
from moleditpy.ui.molecule_scene import MoleculeScene
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem
from unittest.mock import MagicMock, patch


@pytest.fixture
def scene_setup(app):
    mock_window = MagicMock()
    mock_window.is_2d_editable = True

    # Mock settings as a dict with value() method
    class MockSettings(dict):
        def value(self, key, default=None):
            return self.get(key, default)

    mock_window.settings = MockSettings(
        {"atom_label_font_size": 10, "bond_width": 2.0, "atom_color_C": "#000000"}
    )

    # Use real dicts for data
    data = MagicMock()
    data.atoms = {}
    data.bonds = {}

    scene = MoleculeScene(data, mock_window)
    mock_view = MagicMock()
    mock_view.transform.return_value = None
    scene.views = MagicMock(return_value=[mock_view])

    # Use real QPointF for all coordinates
    scene.press_pos = None

    yield scene, data, mock_window


def create_mock_event(pos, button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.button.return_value = button
    event.pos.return_value = pos
    event.scenePos.return_value = pos
    return event


def test_scene_ez_toggle_logic(scene_setup):
    """Test E/Z stereo toggling logic in mouseReleaseEvent."""
    scene, data, win = scene_setup
    scene.mode = "bond_2_5"

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=2, stereo=0)
    data.bonds[(0, 1)] = {"order": 2, "stereo": 0, "item": bond}

    pos = QPointF(25, 0)
    event = create_mock_event(pos)

    with patch.object(MoleculeScene, "itemAt", return_value=bond):
        scene.press_pos = pos  # Manually set to ensure is_click
        scene.mouseReleaseEvent(event)

    assert bond.stereo == 3
    assert data.bonds[(0, 1)]["stereo"] == 3


def test_scene_item_deletion_path(scene_setup):
    """Test item deletion logic in mouseReleaseEvent."""
    scene, data, win = scene_setup
    scene.mode = "delete"

    atom = AtomItem(0, "C", QPointF(0, 0))
    scene.addItem(atom)
    scene.delete_items = MagicMock(return_value=True)

    pos = QPointF(10, 10)
    event = create_mock_event(pos)

    scene.press_pos = pos  # Manually set to ensure is_click

    with patch.object(MoleculeScene, "itemAt", return_value=atom):
        scene.mouseReleaseEvent(event)

    assert scene.delete_items.called


def test_scene_bond_inversion(scene_setup):
    """Test bond inversion logic in mouseReleaseEvent."""
    scene, data, win = scene_setup
    scene.mode = "bond_1_1"
    scene.bond_order = 1
    scene.bond_stereo = 1

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=1, stereo=1)
    data.bonds[(0, 1)] = {"order": 1, "stereo": 1, "item": bond}

    def mock_add_bond(id1, id2, order=1, stereo=0):
        bid = (id1, id2)
        data.bonds[bid] = {"order": order, "stereo": stereo}
        return (bid, True)

    data.add_bond.side_effect = mock_add_bond

    pos = QPointF(25, 0)
    event = create_mock_event(pos)

    with patch.object(MoleculeScene, "itemAt", return_value=bond):
        scene.press_pos = pos  # Manually set to ensure is_click
        scene.mouseReleaseEvent(event)

    assert data.remove_bond.called
    assert (1, 0) in data.bonds


def test_update_user_template_preview(scene_setup):
    """Test Template preview logic."""
    scene, data, win = scene_setup
    scene.user_template_data = {
        "atoms": [
            {"x": 0, "y": 0, "symbol": "C", "id": 0},
            {"x": 50, "y": 0, "symbol": "C", "id": 1},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 1}],
    }
    scene.template_preview = MagicMock()
    scene.update_user_template_preview(QPointF(100, 100))
    assert "points" in scene.template_context


def test_benzene_template_rotation_logic(scene_setup):
    """Test benzene template rotation alignment."""
    scene, data, win = scene_setup

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    data.atoms = {0: {"item": a1, "symbol": "C"}, 1: {"item": a2, "symbol": "C"}}

    # Mock items() to include a bond
    bond = BondItem(a1, a2, order=1)
    scene.find_bond_between = MagicMock(return_value=bond)
    scene.items = MagicMock(return_value=[bond])

    bonds_info = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
    points = [
        QPointF(0, 0),
        QPointF(50, 0),
        QPointF(75, 43),
        QPointF(50, 86),
        QPointF(0, 86),
        QPointF(-25, 43),
    ]
    atoms_data = [{"symbol": "C", "id": i} for i in range(6)]
    context = {
        "points": points,
        "bonds_info": bonds_info,
        "atoms_data": atoms_data,
        "attachment_atom": None,
    }

    def mock_add_atom(symbol, pos, charge=0, radical=0):
        new_id = len(data.atoms) + 100
        data.atoms[new_id] = {"symbol": symbol, "item": MagicMock(spec=AtomItem)}
        return new_id

    data.add_atom.side_effect = mock_add_atom

    def mock_add_bond(id1, id2, order=1, stereo=0):
        bid = (id1, id2)
        data.bonds[bid] = {"order": order, "stereo": stereo}
        return (bid, True)

    data.add_bond.side_effect = mock_add_bond

    scene.add_user_template_fragment(context)
    assert data.add_atom.called


def test_wedge_flip_via_key_w_pushes_undo(scene_setup):
    """Regression: pressing W on an already-wedge bond flips its direction
    (swaps the data-model key and the bond's atoms) but left order/stereo
    unchanged, so any_bond_changed stayed False and no undo checkpoint was
    recorded for a real data mutation."""
    scene, data, win = scene_setup

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=1, stereo=1)  # already a wedge
    data.bonds[(0, 1)] = {"order": 1, "stereo": 1, "item": bond}

    event = MagicMock()
    event.key.return_value = Qt.Key.Key_W
    event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

    win.ui_manager.is_2d_editable = True
    win.edit_actions_manager.push_undo_state.reset_mock()

    with (
        patch.object(MoleculeScene, "itemAt", return_value=bond),
        patch.object(scene, "find_atom_near", return_value=None),
        patch.object(scene, "update_all_items"),
    ):
        scene.keyPressEvent(event)

    # Direction flipped: key swapped in the data model, atoms swapped on item
    assert (1, 0) in data.bonds and (0, 1) not in data.bonds
    assert bond.atom1 is a2 and bond.atom2 is a1
    # The mutation must be recorded in undo history
    win.edit_actions_manager.push_undo_state.assert_called_once()


# ---------------------------------------------------------------------------
# keyPressEvent — hover Z/E stereo toggle, 1/2/3 bond creation, delete, space
# ---------------------------------------------------------------------------


def _key_event(key, modifiers=Qt.KeyboardModifier.NoModifier):
    event = MagicMock()
    event.key.return_value = key
    event.modifiers.return_value = modifiers
    return event


def _wire_add_atom_bond(data):
    def mock_add_atom(symbol, pos, charge=0, radical=0):
        new_id = len(data.atoms) + 100
        data.atoms[new_id] = {"symbol": symbol, "pos": pos}
        return new_id

    data.add_atom.side_effect = mock_add_atom

    def mock_add_bond(id1, id2, order=1, stereo=0):
        key = (id1, id2)
        data.bonds[key] = {"order": order, "stereo": stereo}
        return (key, "created")

    data.add_bond.side_effect = mock_add_bond


def test_key_z_toggles_hovered_double_bond_to_z_isomer(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=2, stereo=0)
    a1.bonds.append(bond)
    a2.bonds.append(bond)
    data.bonds[(0, 1)] = {"order": 2, "stereo": 0}
    scene.hovered_item = bond

    event = _key_event(Qt.Key.Key_Z)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "update_all_items"),
    ):
        scene.keyPressEvent(event)

    assert data.bonds[(0, 1)]["stereo"] == 3
    assert bond.stereo == 3
    win.edit_actions_manager.push_undo_state.assert_called_once()
    event.accept.assert_called_once()


def test_key_e_toggles_hovered_double_bond_to_e_isomer(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=2, stereo=0)
    a1.bonds.append(bond)
    a2.bonds.append(bond)
    data.bonds[(0, 1)] = {"order": 2, "stereo": 0}
    scene.hovered_item = bond

    event = _key_event(Qt.Key.Key_E)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "update_all_items"),
    ):
        scene.keyPressEvent(event)

    assert data.bonds[(0, 1)]["stereo"] == 4
    win.edit_actions_manager.push_undo_state.assert_called_once()


def test_key_1_from_selected_atom_creates_new_atom_and_bond(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True
    _wire_add_atom_bond(data)

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.atom_id = 0
    data.atoms[0] = {"symbol": "C", "pos": QPointF(0, 0)}
    scene.atom_items[0] = a1
    scene.addItem(a1)
    a1.setSelected(True)

    before_bonds = dict(data.bonds)
    event = _key_event(Qt.Key.Key_1)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "find_atom_near", return_value=None),
    ):
        scene.keyPressEvent(event)

    assert len(data.atoms) == 2  # original + newly created
    assert len(data.bonds) == len(before_bonds) + 1
    win.edit_actions_manager.push_undo_state.assert_called_once()


def test_key_3_bonds_to_existing_nearby_atom_without_new_atom(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True
    _wire_add_atom_bond(data)

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(40, 0))
    a1.atom_id, a2.atom_id = 0, 1
    data.atoms[0] = {"symbol": "C", "pos": QPointF(0, 0)}
    data.atoms[1] = {"symbol": "C", "pos": QPointF(40, 0)}
    scene.atom_items[0] = a1
    scene.atom_items[1] = a2
    scene.addItem(a1)
    scene.addItem(a2)
    a1.setSelected(True)

    event = _key_event(Qt.Key.Key_3)
    # First find_atom_near call resolves the cursor (no atom there); the second,
    # inside the Key_3 handler, resolves the snap target to a2.
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "find_atom_near", side_effect=[None, a2]),
    ):
        scene.keyPressEvent(event)

    assert len(data.atoms) == 2  # no new atom created
    assert (0, 1) in data.bonds
    assert data.bonds[(0, 1)]["order"] == 3
    win.edit_actions_manager.push_undo_state.assert_called_once()


def test_key_delete_aborts_active_temp_line_placement(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.atom_id = 0
    scene.addItem(a1)
    line = scene.addLine(0, 0, 10, 10)
    scene.temp_line = line
    scene.start_atom = a1
    scene.start_pos = QPointF(0, 0)
    scene.initial_positions_in_event = {0: QPointF(0, 0)}

    event = _key_event(Qt.Key.Key_Delete)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "find_atom_near", return_value=None),
    ):
        scene.keyPressEvent(event)

    assert scene.temp_line is None
    assert scene.start_atom is None
    assert scene.start_pos is None
    assert scene.initial_positions_in_event == {}
    event.accept.assert_called_once()


def test_key_delete_removes_selected_atom_and_pushes_undo(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    data.atoms[0] = {"symbol": "C"}
    data.atoms[1] = {"symbol": "C"}
    scene.atom_items[0] = a1
    scene.atom_items[1] = a2
    scene.addItem(a1)
    scene.addItem(a2)
    a2.setSelected(True)
    data.remove_atom = MagicMock()
    data.remove_bond = MagicMock(return_value=False)

    event = _key_event(Qt.Key.Key_Delete)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "find_atom_near", return_value=None),
    ):
        scene.keyPressEvent(event)

    assert 1 not in scene.atom_items
    data.remove_atom.assert_called_once_with(1)
    win.edit_actions_manager.push_undo_state.assert_called_once()
    win.statusBar().showMessage.assert_any_call("Deleted selected items.")


def test_key_delete_clears_scene_when_no_atoms_remain(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True

    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.atom_id = 0
    scene.atom_items[0] = a1
    scene.addItem(a1)
    a1.setSelected(True)
    data.remove_atom = MagicMock()
    data.remove_bond = MagicMock(return_value=False)
    # Simulate that the data model is now empty (real MolecularData would be)
    data.atoms = {}

    event = _key_event(Qt.Key.Key_Delete)
    with (
        patch.object(MoleculeScene, "itemAt", return_value=None),
        patch.object(scene, "find_atom_near", return_value=None),
    ):
        scene.keyPressEvent(event)

    assert scene.temp_line is None
    assert scene.start_atom is None
    event.accept.assert_called_once()


def test_key_space_switches_to_select_mode_when_not_already(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True
    scene.mode = "atom_C"

    event = _key_event(Qt.Key.Key_Space)
    with patch.object(MoleculeScene, "itemAt", return_value=None):
        scene.keyPressEvent(event)

    win.ui_manager.activate_select_mode.assert_called_once()
    win.edit_actions_manager.select_all.assert_not_called()
    event.accept.assert_called_once()


def test_key_space_selects_all_when_already_in_select_mode(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = True
    scene.mode = "select"

    event = _key_event(Qt.Key.Key_Space)
    with patch.object(MoleculeScene, "itemAt", return_value=None):
        scene.keyPressEvent(event)

    win.edit_actions_manager.select_all.assert_called_once()
    win.ui_manager.activate_select_mode.assert_not_called()


def test_key_press_noop_when_not_2d_editable(scene_setup):
    scene, data, win = scene_setup
    win.ui_manager.is_2d_editable = False
    scene.mode = "select"

    event = _key_event(Qt.Key.Key_Space)
    with patch.object(MoleculeScene, "itemAt", return_value=None):
        scene.keyPressEvent(event)

    win.edit_actions_manager.select_all.assert_not_called()


# ---------------------------------------------------------------------------
# SceneQueryMixin — find_atom_near / update_bond_stereo direct tests
# ---------------------------------------------------------------------------


def test_find_atom_near_returns_none_for_none_pos(scene_setup):
    scene, data, win = scene_setup
    assert scene.find_atom_near(None) is None


def test_find_atom_near_finds_close_atom(scene_setup):
    scene, data, win = scene_setup
    # AtomItem.shape() reads the view transform; give it an identity transform.
    scene.views()[0].transform.return_value = QTransform()
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.atom_id = 0
    scene.addItem(a1)
    found = scene.find_atom_near(QPointF(2, 2), tol=14.0)
    assert found is a1


def test_find_atom_near_returns_none_when_out_of_tolerance(scene_setup):
    scene, data, win = scene_setup
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a1.atom_id = 0
    scene.addItem(a1)
    assert scene.find_atom_near(QPointF(500, 500), tol=14.0) is None


def test_update_bond_stereo_noop_for_none_bond(scene_setup):
    scene, data, win = scene_setup
    scene.update_bond_stereo(None, 3)  # must not raise


def test_update_bond_stereo_noop_for_single_bond(scene_setup):
    scene, data, win = scene_setup
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=1, stereo=0)
    scene.update_bond_stereo(bond, 3)
    assert bond.stereo == 0  # unchanged: not a double bond


def test_update_bond_stereo_warns_when_bond_missing_from_model(scene_setup):
    scene, data, win = scene_setup
    a1 = AtomItem(0, "C", QPointF(0, 0))
    a2 = AtomItem(1, "C", QPointF(50, 0))
    a1.atom_id, a2.atom_id = 0, 1
    bond = BondItem(a1, a2, order=2, stereo=0)
    data.bonds = {}  # bond not registered

    scene.update_bond_stereo(bond, 3)

    win.statusBar().showMessage.assert_called_once()
    assert "not found in model" in win.statusBar().showMessage.call_args.args[0]
