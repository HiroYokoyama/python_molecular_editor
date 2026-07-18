"""
Additional coverage tests for ui/molecule_scene.py.

Fills gaps left by existing test_scene_interactions.py / test_scene_advanced.py.

Covers:
  - get_setting: safe gateway to settings with fallback
  - update_connected_bonds: calls update_position on all connected bonds
  - clear_all_problem_flags: resets has_problem, returns True only when flags were set
  - purge_deleted_items: clears _deleted_items list
  - set_hovered_item: stores reference
  - E/Z stereo cycling (bond_2_5 mode): None->Z(3)->E(4)->None(0)
  - Bond direction inversion when same stereo clicked again
  - mouseDoubleClickEvent in select mode: BFS connected-component selection
  - mouseDoubleClickEvent in charge/radical modes
"""

import os
import sys
from unittest.mock import MagicMock, patch

from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QTransform

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from moleditpy.ui.molecule_scene import MoleculeScene
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _scene(host):
    scene = MoleculeScene(host.state_manager.data, host)
    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()
    mock_view.viewport.return_value = MagicMock()
    scene.views = MagicMock(return_value=[mock_view])
    return scene


def _event(pos=QPointF(100, 100), button=Qt.MouseButton.LeftButton):
    ev = MagicMock()
    ev.scenePos.return_value = pos
    ev.button.return_value = button
    return ev


# ---------------------------------------------------------------------------
# get_setting
# ---------------------------------------------------------------------------


class TestGetSetting:
    def test_returns_value_from_settings(self, mock_parser_host):
        """get_setting returns the value stored in init_manager.settings."""
        mock_parser_host.init_manager.settings["my_key"] = "my_value"
        scene = _scene(mock_parser_host)
        assert scene.get_setting("my_key") == "my_value"

    def test_returns_default_when_key_missing(self, mock_parser_host):
        """get_setting returns the default value when the key is absent."""
        scene = _scene(mock_parser_host)
        assert scene.get_setting("no_such_key", "fallback") == "fallback"

    def test_returns_default_when_window_is_none(self, mock_parser_host):
        """get_setting returns the default value when scene.window is None."""
        scene = _scene(mock_parser_host)
        scene.window = None
        assert scene.get_setting("any_key", 99) == 99


# ---------------------------------------------------------------------------
# update_connected_bonds
# ---------------------------------------------------------------------------


class TestUpdateConnectedBonds:
    def test_calls_update_position_on_bond(self, mock_parser_host):
        """update_connected_bonds calls update_position on each live bond."""
        scene = _scene(mock_parser_host)
        bond = MagicMock(spec=BondItem)
        bond.update_position = MagicMock()
        atom = MagicMock(spec=AtomItem)
        atom.bonds = [bond]

        with patch(
            "moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False
        ):
            scene.update_connected_bonds([atom])

        bond.update_position.assert_called_once()

    def test_skips_sip_deleted_bonds(self, mock_parser_host):
        """update_connected_bonds skips bonds that have been sip-deleted."""
        scene = _scene(mock_parser_host)
        bond = MagicMock(spec=BondItem)
        bond.update_position = MagicMock()
        atom = MagicMock(spec=AtomItem)
        atom.bonds = [bond]

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=True):
            scene.update_connected_bonds([atom])

        bond.update_position.assert_not_called()

    def test_deduplicates_shared_bond(self, mock_parser_host):
        """Bond shared by two atoms should only be updated once."""
        scene = _scene(mock_parser_host)
        bond = MagicMock(spec=BondItem)
        bond.update_position = MagicMock()
        a1 = MagicMock(spec=AtomItem)
        a1.bonds = [bond]
        a2 = MagicMock(spec=AtomItem)
        a2.bonds = [bond]

        with patch(
            "moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False
        ):
            scene.update_connected_bonds([a1, a2])

        bond.update_position.assert_called_once()

    def test_handles_atom_without_bonds_attribute(self, mock_parser_host):
        """update_connected_bonds does not raise when an atom has no bonds attribute."""
        scene = _scene(mock_parser_host)
        atom = object()  # no 'bonds' attribute
        scene.update_connected_bonds([atom])  # must not raise


# ---------------------------------------------------------------------------
# clear_all_problem_flags
# ---------------------------------------------------------------------------


class TestClearAllProblemFlags:
    def test_returns_true_when_flags_were_set(self, mock_parser_host):
        """clear_all_problem_flags returns True when at least one atom had has_problem set."""
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = scene.atom_items[aid]
        item.has_problem = True
        item.update = MagicMock()

        assert scene.clear_all_problem_flags() is True

    def test_returns_false_when_no_flags_set(self, mock_parser_host):
        """clear_all_problem_flags returns False when no atom has has_problem set."""
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = scene.atom_items[aid]
        item.has_problem = False

        assert scene.clear_all_problem_flags() is False

    def test_flag_is_reset_to_false(self, mock_parser_host):
        """clear_all_problem_flags resets has_problem to False on each atom."""
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = scene.atom_items[aid]
        item.has_problem = True
        item.update = MagicMock()

        scene.clear_all_problem_flags()
        assert item.has_problem is False


# ---------------------------------------------------------------------------
# purge_deleted_items
# ---------------------------------------------------------------------------


class TestPurgeDeletedItems:
    def test_noop_on_empty_list(self, mock_parser_host):
        """purge_deleted_items is a no-op when the deleted items list is empty."""
        scene = _scene(mock_parser_host)
        scene._deleted_items = []
        scene.purge_deleted_items()  # must not raise

    def test_clears_the_list(self, mock_parser_host):
        """purge_deleted_items empties the _deleted_items list."""
        scene = _scene(mock_parser_host)
        obj = MagicMock()
        obj.hide = MagicMock()
        obj.bonds = None
        scene._deleted_items = [obj]

        with patch(
            "moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False
        ):
            scene.purge_deleted_items()

        assert scene._deleted_items == []

    def test_calls_hide_on_valid_objects(self, mock_parser_host):
        """purge_deleted_items calls hide() on each live item."""
        scene = _scene(mock_parser_host)
        obj = MagicMock()
        obj.hide = MagicMock()
        obj.bonds = None
        scene._deleted_items = [obj]

        with patch(
            "moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False
        ):
            scene.purge_deleted_items()

        obj.hide.assert_called_once()

    def test_skips_already_sip_deleted(self, mock_parser_host):
        """purge_deleted_items skips items that have already been sip-deleted."""
        scene = _scene(mock_parser_host)
        obj = MagicMock()
        obj.hide = MagicMock()
        scene._deleted_items = [obj]

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=True):
            scene.purge_deleted_items()

        obj.hide.assert_not_called()


# ---------------------------------------------------------------------------
# set_hovered_item
# ---------------------------------------------------------------------------


class TestSetHoveredItem:
    def test_stores_item(self, mock_parser_host):
        """set_hovered_item stores the provided item reference."""
        scene = _scene(mock_parser_host)
        item = MagicMock()
        scene.set_hovered_item(item)
        assert scene.hovered_item is item

    def test_accepts_none(self, mock_parser_host):
        """set_hovered_item accepts None to clear the hovered item."""
        scene = _scene(mock_parser_host)
        scene.set_hovered_item(None)
        assert scene.hovered_item is None


# ---------------------------------------------------------------------------
# E/Z stereo cycling (bond_2_5 mode)
# ---------------------------------------------------------------------------


class TestEZStereoCycling:
    def _setup(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "bond_2_5"

        a1_id = scene.create_atom("C", QPointF(0, 0))
        a2_id = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(
            scene.atom_items[a1_id],
            scene.atom_items[a2_id],
            bond_order=2,
        )
        bond_key = (a1_id, a2_id)
        bond_item = scene.bond_items[bond_key]
        bond_item.order = 2
        bond_item.stereo = 0
        bond_item.set_stereo = MagicMock()
        return scene, bond_item

    def test_none_to_z_on_first_click(self, mock_parser_host):
        """Clicking a double bond with no stereo sets it to Z (stereo=3)."""
        scene, bond_item = self._setup(mock_parser_host)

        with (
            patch(
                "moleditpy.ui.molecule_scene.QApplication.startDragDistance",
                return_value=1000,
            ),
            patch.object(MoleculeScene, "itemAt", return_value=bond_item),
            patch.object(MoleculeScene, "update_bond_stereo") as mock_stereo,
            patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent"),
            patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent"),
        ):
            ev = _event()
            scene.mousePressEvent(ev)
            scene.mouseReleaseEvent(ev)

        mock_stereo.assert_called_once_with(bond_item, 3)

    def test_z_to_e_on_click(self, mock_parser_host):
        """Clicking a Z double bond advances it to E (stereo=4)."""
        scene, bond_item = self._setup(mock_parser_host)
        bond_item.stereo = 3  # already Z

        with (
            patch(
                "moleditpy.ui.molecule_scene.QApplication.startDragDistance",
                return_value=1000,
            ),
            patch.object(MoleculeScene, "itemAt", return_value=bond_item),
            patch.object(MoleculeScene, "update_bond_stereo") as mock_stereo,
            patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent"),
            patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent"),
        ):
            ev = _event()
            scene.mousePressEvent(ev)
            scene.mouseReleaseEvent(ev)

        mock_stereo.assert_called_once_with(bond_item, 4)

    def test_e_to_none_on_click(self, mock_parser_host):
        """Clicking an E double bond cycles it back to no stereo (stereo=0)."""
        scene, bond_item = self._setup(mock_parser_host)
        bond_item.stereo = 4  # already E

        with (
            patch(
                "moleditpy.ui.molecule_scene.QApplication.startDragDistance",
                return_value=1000,
            ),
            patch.object(MoleculeScene, "itemAt", return_value=bond_item),
            patch.object(MoleculeScene, "update_bond_stereo") as mock_stereo,
            patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent"),
            patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent"),
        ):
            ev = _event()
            scene.mousePressEvent(ev)
            scene.mouseReleaseEvent(ev)

        mock_stereo.assert_called_once_with(bond_item, 0)


# ---------------------------------------------------------------------------
# Bond direction inversion
# ---------------------------------------------------------------------------


class TestBondDirectionInversion:
    def test_stereo_bond_click_inverts_atom_order(self, mock_parser_host):
        """Clicking a stereo bond with matching stereo type inverts the atom order."""
        scene = _scene(mock_parser_host)
        scene.mode = "bond_1"
        scene.bond_order = 1
        scene.bond_stereo = 1  # Wedge

        a1_id = scene.create_atom("C", QPointF(0, 0))
        a2_id = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(
            scene.atom_items[a1_id],
            scene.atom_items[a2_id],
            bond_order=1,
            bond_stereo=1,
        )
        bond_item = scene.bond_items[(a1_id, a2_id)]
        bond_item.order = 1
        bond_item.stereo = 1
        bond_item.atom1 = scene.atom_items[a1_id]
        bond_item.atom2 = scene.atom_items[a2_id]
        bond_item.update_position = MagicMock()

        with (
            patch(
                "moleditpy.ui.molecule_scene.QApplication.startDragDistance",
                return_value=1000,
            ),
            patch.object(MoleculeScene, "itemAt", return_value=bond_item),
            patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent"),
            patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent"),
        ):
            ev = _event()
            scene.mousePressEvent(ev)
            scene.mouseReleaseEvent(ev)

        # Reversed bond key should now exist
        assert (a2_id, a1_id) in mock_parser_host.state_manager.data.bonds


# ---------------------------------------------------------------------------
# mouseDoubleClickEvent — select mode BFS
# ---------------------------------------------------------------------------


class TestDoubleClickSelectMode:
    def test_selects_connected_atom(self, mock_parser_host):
        """Double-click in select mode selects the connected component via BFS."""
        scene = _scene(mock_parser_host)
        scene.mode = "select"

        a1_id = scene.create_atom("C", QPointF(0, 0))
        a2_id = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(
            scene.atom_items[a1_id],
            scene.atom_items[a2_id],
        )
        atom_item = scene.atom_items[a1_id]
        atom_item.setSelected = MagicMock()

        with (
            patch.object(MoleculeScene, "itemAt", return_value=atom_item),
            patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False),
        ):
            scene.mouseDoubleClickEvent(_event())

        atom_item.setSelected.assert_called_with(True)

    def test_bond_2_5_mode_accepts_event(self, mock_parser_host):
        """Double-click in bond_2_5 mode accepts the event and does not raise."""
        scene = _scene(mock_parser_host)
        scene.mode = "bond_2_5"
        ev = _event()
        ev.accept = MagicMock()
        with patch.object(MoleculeScene, "itemAt", return_value=None):
            scene.mouseDoubleClickEvent(ev)
        ev.accept.assert_called()


# ---------------------------------------------------------------------------
# mouseDoubleClickEvent — charge / radical modes
# ---------------------------------------------------------------------------


class TestDoubleClickChargeRadical:
    def test_radical_increments_on_double_click(self, mock_parser_host):
        """Double-clicking an atom in radical mode increments its radical count."""
        scene = _scene(mock_parser_host)
        scene.mode = "radical"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = scene.atom_items[aid]
        item.radical = 0
        item.prepareGeometryChange = MagicMock()
        item.update_style = MagicMock()

        with (
            patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
            patch.object(MoleculeScene, "itemAt", return_value=item),
        ):
            scene.mouseDoubleClickEvent(_event(QPointF(100, 100)))

        assert item.radical == 1

    def test_charge_plus_increments_on_double_click(self, mock_parser_host):
        """Double-clicking an atom in charge_plus mode increments its charge."""
        scene = _scene(mock_parser_host)
        scene.mode = "charge_plus"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = scene.atom_items[aid]
        item.charge = 0
        item.prepareGeometryChange = MagicMock()
        item.update_style = MagicMock()

        with (
            patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
            patch.object(MoleculeScene, "itemAt", return_value=item),
        ):
            scene.mouseDoubleClickEvent(_event(QPointF(100, 100)))

        assert item.charge == 1

    def test_charge_minus_decrements_on_double_click(self, mock_parser_host):
        """Double-clicking an atom in charge_minus mode decrements its charge."""
        scene = _scene(mock_parser_host)
        scene.mode = "charge_minus"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = scene.atom_items[aid]
        item.charge = 0
        item.prepareGeometryChange = MagicMock()
        item.update_style = MagicMock()

        with (
            patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
            patch.object(MoleculeScene, "itemAt", return_value=item),
        ):
            scene.mouseDoubleClickEvent(_event(QPointF(100, 100)))

        assert item.charge == -1


# ---------------------------------------------------------------------------
# restore_atoms_and_bonds — undo/redo dict-of-dicts restoration
# ---------------------------------------------------------------------------


class TestRestoreAtomsAndBonds:
    def test_restores_atoms_and_bonds_from_state(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        raw_atoms = {
            1: {"symbol": "C", "pos": (0.0, 0.0), "charge": 0, "radical": 0},
            2: {"symbol": "O", "pos": (50.0, 0.0), "charge": -1, "radical": 1},
        }
        raw_bonds = {(1, 2): {"order": 2, "stereo": 0}}

        scene.restore_atoms_and_bonds(raw_atoms, raw_bonds)

        assert set(scene.atom_items.keys()) == {1, 2}
        assert (1, 2) in scene.bond_items
        assert scene.data.atoms[2]["charge"] == -1
        assert scene.data.atoms[2]["radical"] == 1
        assert scene.data.bonds[(1, 2)]["order"] == 2
        # Bond registered on both endpoint items
        assert len(scene.atom_items[1].bonds) == 1
        assert len(scene.atom_items[2].bonds) == 1

    def test_skips_bond_when_endpoint_atom_missing(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        raw_atoms = {1: {"symbol": "C", "pos": (0.0, 0.0)}}
        raw_bonds = {(1, 99): {"order": 1, "stereo": 0}}

        scene.restore_atoms_and_bonds(raw_atoms, raw_bonds)

        assert scene.bond_items == {}
        assert 1 in scene.atom_items


# ---------------------------------------------------------------------------
# restore_atoms_and_bonds_from_json — PMEPRJ list-of-dicts restoration
# ---------------------------------------------------------------------------


class TestRestoreFromJson:
    def test_restores_from_json_lists(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        atoms_2d = [
            {"id": 1, "symbol": "C", "x": 0.0, "y": 0.0},
            {"id": 2, "symbol": "N", "x": 50.0, "y": 0.0, "charge": 1},
        ]
        bonds_2d = [{"atom1": 1, "atom2": 2, "order": 1, "stereo": 2}]

        scene.restore_atoms_and_bonds_from_json(atoms_2d, bonds_2d)

        assert set(scene.atom_items.keys()) == {1, 2}
        assert scene.data.atoms[2]["charge"] == 1
        assert scene.data.bonds[(1, 2)]["stereo"] == 2
        assert (1, 2) in scene.bond_items

    def test_json_skips_bond_with_unknown_atom(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        atoms_2d = [{"id": 1, "symbol": "C", "x": 0.0, "y": 0.0}]
        bonds_2d = [{"atom1": 1, "atom2": 7, "order": 1}]

        scene.restore_atoms_and_bonds_from_json(atoms_2d, bonds_2d)

        assert scene.bond_items == {}


# ---------------------------------------------------------------------------
# update_ring_info_2d
# ---------------------------------------------------------------------------


class TestUpdateRingInfo2D:
    def test_noop_when_no_atoms(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.update_ring_info_2d()  # must not raise

    def test_marks_ring_bonds_and_sets_center(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        # Cyclopropane: three carbons forming a triangle
        raw_atoms = {
            1: {"symbol": "C", "pos": (0.0, 0.0)},
            2: {"symbol": "C", "pos": (50.0, 0.0)},
            3: {"symbol": "C", "pos": (25.0, 43.0)},
        }
        raw_bonds = {
            (1, 2): {"order": 1, "stereo": 0},
            (2, 3): {"order": 1, "stereo": 0},
            (1, 3): {"order": 1, "stereo": 0},
        }
        scene.restore_atoms_and_bonds(raw_atoms, raw_bonds)

        scene.update_ring_info_2d()

        assert all(b.is_in_ring for b in scene.bond_items.values())
        assert all(b.ring_center is not None for b in scene.bond_items.values())

    def test_acyclic_bonds_not_marked(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        raw_atoms = {
            1: {"symbol": "C", "pos": (0.0, 0.0)},
            2: {"symbol": "C", "pos": (50.0, 0.0)},
        }
        raw_bonds = {(1, 2): {"order": 1, "stereo": 0}}
        scene.restore_atoms_and_bonds(raw_atoms, raw_bonds)

        scene.update_ring_info_2d()

        assert not scene.bond_items[(1, 2)].is_in_ring


# ---------------------------------------------------------------------------
# refresh_mode_state / leaveEvent
# ---------------------------------------------------------------------------


class TestRefreshModeState:
    def test_updates_template_preview_when_cursor_in_viewport(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "template_benzene"
        scene.update_template_preview = MagicMock()
        view = scene.views()[0]
        view.isVisible.return_value = True
        view.viewport().rect().contains.return_value = True
        view.mapToScene.return_value = QPointF(10, 20)

        scene.refresh_mode_state()

        scene.update_template_preview.assert_called_once_with(QPointF(10, 20))

    def test_no_preview_when_cursor_outside_viewport(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "template_benzene"
        scene.update_template_preview = MagicMock()
        view = scene.views()[0]
        view.isVisible.return_value = True
        view.viewport().rect().contains.return_value = False

        scene.refresh_mode_state()

        scene.update_template_preview.assert_not_called()

    def test_no_preview_when_not_template_mode(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "select"
        scene.update_template_preview = MagicMock()
        view = scene.views()[0]
        view.isVisible.return_value = True
        view.viewport().rect().contains.return_value = True

        scene.refresh_mode_state()

        scene.update_template_preview.assert_not_called()


class TestLeaveEvent:
    def test_hides_template_preview(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.template_preview = MagicMock()
        scene.leaveEvent(MagicMock())
        scene.template_preview.hide.assert_called_once()


# ---------------------------------------------------------------------------
# mousePressEvent — right-button mode handlers
# ---------------------------------------------------------------------------


def _rclick(pos=QPointF(10, 10)):
    return _event(pos=pos, button=Qt.MouseButton.RightButton)


class TestRightClickPress:
    def test_radical_mode_resets_radical(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "radical"
        aid = scene.create_atom("C", QPointF(10, 10))
        item = scene.atom_items[aid]
        item.radical = 2
        scene.data.atoms[aid]["radical"] = 2

        ev = _rclick()
        with patch.object(MoleculeScene, "itemAt", return_value=item):
            scene.mousePressEvent(ev)

        assert item.radical == 0
        assert scene.data.atoms[aid]["radical"] == 0
        ev.accept.assert_called_once()

    def test_charge_mode_resets_charge(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "charge_plus"
        aid = scene.create_atom("N", QPointF(10, 10))
        item = scene.atom_items[aid]
        item.charge = 1
        scene.data.atoms[aid]["charge"] = 1

        ev = _rclick()
        with patch.object(MoleculeScene, "itemAt", return_value=item):
            scene.mousePressEvent(ev)

        assert item.charge == 0
        assert scene.data.atoms[aid]["charge"] == 0

    def test_bond25_right_click_clears_ez_label(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "bond_2_5"
        a1 = scene.create_atom("C", QPointF(0, 0))
        a2 = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(scene.atom_items[a1], scene.atom_items[a2], bond_order=2)
        bond = scene.bond_items[(a1, a2)]
        bond.order = 2
        bond.stereo = 3  # currently Z
        bond.set_stereo = MagicMock()
        scene.data.bonds[(a1, a2)]["stereo"] = 3

        ev = _rclick()
        with patch.object(MoleculeScene, "itemAt", return_value=bond):
            scene.mousePressEvent(ev)

        bond.set_stereo.assert_called_once_with(0)
        assert scene.data.bonds[(a1, a2)]["stereo"] == 0
        mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()

    def test_delete_mode_right_click_removes_atom(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "delete"
        aid = scene.create_atom("C", QPointF(10, 10))
        item = scene.atom_items[aid]

        ev = _rclick()
        with (
            patch.object(MoleculeScene, "itemAt", return_value=item),
            patch.object(scene, "delete_items", return_value=True) as del_items,
        ):
            scene.mousePressEvent(ev)

        del_items.assert_called_once()
        assert item in del_items.call_args.args[0]
        mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()

    def test_right_click_multi_selection_deletes_all(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "delete"
        a1 = scene.create_atom("C", QPointF(0, 0))
        a2 = scene.create_atom("C", QPointF(50, 0))
        it1, it2 = scene.atom_items[a1], scene.atom_items[a2]
        scene.selectedItems = MagicMock(return_value=[it1, it2])

        ev = _rclick()
        with (
            patch.object(MoleculeScene, "itemAt", return_value=it1),
            patch.object(scene, "delete_items", return_value=True) as del_items,
        ):
            scene.mousePressEvent(ev)

        deleted = del_items.call_args.args[0]
        assert it1 in deleted and it2 in deleted
        ev.accept.assert_called_once()

    def test_right_click_empty_space_is_noop(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.mode = "delete"

        ev = _rclick()
        with patch.object(MoleculeScene, "itemAt", return_value=None):
            scene.mousePressEvent(ev)

        ev.accept.assert_not_called()
        mock_parser_host.edit_actions_manager.push_undo_state.assert_not_called()
