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
import pytest
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
        mock_parser_host.init_manager.settings["my_key"] = "my_value"
        scene = _scene(mock_parser_host)
        assert scene.get_setting("my_key") == "my_value"

    def test_returns_default_when_key_missing(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        assert scene.get_setting("no_such_key", "fallback") == "fallback"

    def test_returns_default_when_window_is_none(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene.window = None
        assert scene.get_setting("any_key", 99) == 99


# ---------------------------------------------------------------------------
# update_connected_bonds
# ---------------------------------------------------------------------------


class TestUpdateConnectedBonds:
    def test_calls_update_position_on_bond(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        bond = MagicMock(spec=BondItem)
        bond.update_position = MagicMock()
        atom = MagicMock(spec=AtomItem)
        atom.bonds = [bond]

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False):
            scene.update_connected_bonds([atom])

        bond.update_position.assert_called_once()

    def test_skips_sip_deleted_bonds(self, mock_parser_host):
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

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False):
            scene.update_connected_bonds([a1, a2])

        bond.update_position.assert_called_once()

    def test_handles_atom_without_bonds_attribute(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        atom = object()  # no 'bonds' attribute
        scene.update_connected_bonds([atom])  # must not raise


# ---------------------------------------------------------------------------
# clear_all_problem_flags
# ---------------------------------------------------------------------------


class TestClearAllProblemFlags:
    def test_returns_true_when_flags_were_set(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
        item.has_problem = True
        item.update = MagicMock()

        assert scene.clear_all_problem_flags() is True

    def test_returns_false_when_no_flags_set(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
        item.has_problem = False

        assert scene.clear_all_problem_flags() is False

    def test_flag_is_reset_to_false(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        aid = scene.create_atom("C", QPointF(0, 0))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
        item.has_problem = True
        item.update = MagicMock()

        scene.clear_all_problem_flags()
        assert item.has_problem is False


# ---------------------------------------------------------------------------
# purge_deleted_items
# ---------------------------------------------------------------------------


class TestPurgeDeletedItems:
    def test_noop_on_empty_list(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        scene._deleted_items = []
        scene.purge_deleted_items()  # must not raise

    def test_clears_the_list(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        obj = MagicMock()
        obj.hide = MagicMock()
        obj.bonds = None
        scene._deleted_items = [obj]

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False):
            scene.purge_deleted_items()

        assert scene._deleted_items == []

    def test_calls_hide_on_valid_objects(self, mock_parser_host):
        scene = _scene(mock_parser_host)
        obj = MagicMock()
        obj.hide = MagicMock()
        obj.bonds = None
        scene._deleted_items = [obj]

        with patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False):
            scene.purge_deleted_items()

        obj.hide.assert_called_once()

    def test_skips_already_sip_deleted(self, mock_parser_host):
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
        scene = _scene(mock_parser_host)
        item = MagicMock()
        scene.set_hovered_item(item)
        assert scene.hovered_item is item

    def test_accepts_none(self, mock_parser_host):
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
            mock_parser_host.state_manager.data.atoms[a1_id]["item"],
            mock_parser_host.state_manager.data.atoms[a2_id]["item"],
            bond_order=2,
        )
        bond_key = (a1_id, a2_id)
        bond_item = mock_parser_host.state_manager.data.bonds[bond_key]["item"]
        bond_item.order = 2
        bond_item.stereo = 0
        bond_item.set_stereo = MagicMock()
        return scene, bond_item

    def test_none_to_z_on_first_click(self, mock_parser_host):
        scene, bond_item = self._setup(mock_parser_host)

        with (
            patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000),
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
        scene, bond_item = self._setup(mock_parser_host)
        bond_item.stereo = 3  # already Z

        with (
            patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000),
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
        scene, bond_item = self._setup(mock_parser_host)
        bond_item.stereo = 4  # already E

        with (
            patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000),
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
        scene = _scene(mock_parser_host)
        scene.mode = "bond_1"
        scene.bond_order = 1
        scene.bond_stereo = 1  # Wedge

        a1_id = scene.create_atom("C", QPointF(0, 0))
        a2_id = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(
            mock_parser_host.state_manager.data.atoms[a1_id]["item"],
            mock_parser_host.state_manager.data.atoms[a2_id]["item"],
            bond_order=1,
            bond_stereo=1,
        )
        bond_item = mock_parser_host.state_manager.data.bonds[(a1_id, a2_id)]["item"]
        bond_item.order = 1
        bond_item.stereo = 1
        bond_item.atom1 = mock_parser_host.state_manager.data.atoms[a1_id]["item"]
        bond_item.atom2 = mock_parser_host.state_manager.data.atoms[a2_id]["item"]
        bond_item.update_position = MagicMock()

        with (
            patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000),
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
        scene = _scene(mock_parser_host)
        scene.mode = "select"

        a1_id = scene.create_atom("C", QPointF(0, 0))
        a2_id = scene.create_atom("C", QPointF(50, 0))
        scene.create_bond(
            mock_parser_host.state_manager.data.atoms[a1_id]["item"],
            mock_parser_host.state_manager.data.atoms[a2_id]["item"],
        )
        atom_item = mock_parser_host.state_manager.data.atoms[a1_id]["item"]
        atom_item.setSelected = MagicMock()

        with (
            patch.object(MoleculeScene, "itemAt", return_value=atom_item),
            patch("moleditpy.ui.molecule_scene.sip_isdeleted_safe", return_value=False),
        ):
            scene.mouseDoubleClickEvent(_event())

        atom_item.setSelected.assert_called_with(True)

    def test_bond_2_5_mode_accepts_event(self, mock_parser_host):
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
        scene = _scene(mock_parser_host)
        scene.mode = "radical"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
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
        scene = _scene(mock_parser_host)
        scene.mode = "charge_plus"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
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
        scene = _scene(mock_parser_host)
        scene.mode = "charge_minus"

        aid = scene.create_atom("C", QPointF(100, 100))
        item = mock_parser_host.state_manager.data.atoms[aid]["item"]
        item.charge = 0
        item.prepareGeometryChange = MagicMock()
        item.update_style = MagicMock()

        with (
            patch("moleditpy.ui.molecule_scene.isinstance", return_value=True),
            patch.object(MoleculeScene, "itemAt", return_value=item),
        ):
            scene.mouseDoubleClickEvent(_event(QPointF(100, 100)))

        assert item.charge == -1
