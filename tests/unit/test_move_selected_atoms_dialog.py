"""
Unit tests for MoveSelectedAtomsDialog (ui/move_selected_atoms_dialog.py).

Covers:
  - on_atom_picked: toggles clicked atom individually (no BFS)
  - apply_translation: no-selection guard, invalid input guard, correct coordinate update for selected atoms only
  - apply_rotation: no-selection guard, invalid input guard, rotation math (90° around centroid) for selected atoms only
  - reset_translation_inputs / reset_rotation_inputs
  - clear_selection: clears all state
"""

import os
import sys

import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

# ---------------------------------------------------------------------------
# Molecules
# ---------------------------------------------------------------------------


def _ethane():
    """Ethane with 3D coords: 2 C atoms (idx 0,1) + 6 H atoms (idx 2-7)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _make_main_window(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None):
        from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_main_window(_mol)
        # Patch show_atom_labels and clear_atom_labels to avoid pyvista calls
        with (
            patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
            patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
        ):
            dlg = MoveSelectedAtomsDialog(_mol, mw)
        entry = (dlg, _mol, mw)
        created.append(entry)
        return entry

    yield _factory

    for dlg, *_ in created:
        try:
            dlg.picking_enabled = False
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# on_atom_picked
# ---------------------------------------------------------------------------


class TestOnAtomPicked:
    def test_picks_atom_individually_no_bfs(self, make_dialog):
        """Picking an atom should select ONLY that atom (no BFS connected components)."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
        assert dlg.selected_atoms == {0}

    def test_picking_again_toggles_deselect(self, make_dialog):
        """Re-picking the same atom deselects it."""
        dlg, mol, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
            assert 0 in dlg.selected_atoms
            dlg.on_atom_picked(0)
            assert 0 not in dlg.selected_atoms


# ---------------------------------------------------------------------------
# apply_translation
# ---------------------------------------------------------------------------


class TestApplyTranslation:
    def test_no_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_trans_input.setText("bad")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_translation_updates_only_selected_atoms(self, make_dialog):
        dlg, mol, _ = make_dialog()
        dlg.selected_atoms.add(0)  # select carbon 0 only

        before = mol.GetConformer().GetPositions().copy()
        dlg.x_trans_input.setText("1.0")
        dlg.y_trans_input.setText("2.0")
        dlg.z_trans_input.setText("3.0")

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_translation()

        after = mol.GetConformer().GetPositions()
        # Atom 0 should have shifted
        assert after[0] == pytest.approx(before[0] + [1.0, 2.0, 3.0], abs=1e-4)
        # Atom 1 (unselected carbon) should NOT have shifted
        assert after[1] == pytest.approx(before[1], abs=1e-4)


# ---------------------------------------------------------------------------
# apply_rotation
# ---------------------------------------------------------------------------


class TestApplyRotation:
    def test_no_atoms_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_rot_input.setText("not_a_number")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_rotation_updates_only_selected_atoms_around_centroid(self, make_dialog):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import Geometry

        # Build a 2-atom molecule at known positions: [1,0,0] and [-1,0,0], and another at [5,5,5]
        mol = Chem.MolFromSmiles("[H][H].[H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(1.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Geometry.Point3D(-1.0, 0.0, 0.0))
        conf.SetAtomPosition(2, Geometry.Point3D(5.0, 5.0, 5.0))

        dlg, mol, mw = make_dialog(mol=mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [5.0, 5.0, 5.0]]
        )

        dlg.selected_atoms.update([0, 1])  # select atoms 0 and 1
        dlg.z_rot_input.setText("90.0")

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()

        after = mol.GetConformer().GetPositions()
        # Centroid of 0 and 1 is [0,0,0]. 90° rotation around Z moves them to [0,1,0] and [0,-1,0]
        assert after[0] == pytest.approx([0.0, 1.0, 0.0], abs=1e-5)
        assert after[1] == pytest.approx([0.0, -1.0, 0.0], abs=1e-5)
        # Atom 2 (unselected) should NOT rotate
        assert after[2] == pytest.approx([5.0, 5.0, 5.0], abs=1e-5)


# ---------------------------------------------------------------------------
# Reset inputs
# ---------------------------------------------------------------------------


class TestResetInputs:
    def test_reset_translation_inputs(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.x_trans_input.setText("5.5")
        dlg.y_trans_input.setText("-3.0")
        dlg.z_trans_input.setText("1.2")
        dlg.reset_translation_inputs()
        assert dlg.x_trans_input.text() == "0.0"
        assert dlg.y_trans_input.text() == "0.0"
        assert dlg.z_trans_input.text() == "0.0"

    def test_reset_rotation_inputs(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.x_rot_input.setText("45.0")
        dlg.y_rot_input.setText("90.0")
        dlg.z_rot_input.setText("180.0")
        dlg.reset_rotation_inputs()
        assert dlg.x_rot_input.text() == "0.0"
        assert dlg.y_rot_input.text() == "0.0"
        assert dlg.z_rot_input.text() == "0.0"


# ---------------------------------------------------------------------------
# clear_selection
# ---------------------------------------------------------------------------


class TestClearSelection:
    def test_clear_removes_selected_atoms(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.update([0, 1])
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0


# ---------------------------------------------------------------------------
# Extra tests for Properties, Highlights, and Event Filter
# ---------------------------------------------------------------------------


def test_group_atoms_property(make_dialog):
    dlg, _, _ = make_dialog()
    dlg.group_atoms = {0, 2}
    assert dlg.selected_atoms == {0, 2}
    assert dlg.group_atoms == {0, 2}


def test_preselected_atoms_init(qapp):
    from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

    mol = _ethane()
    mw = _make_main_window(mol)
    with (
        patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
        patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
    ):
        dlg = MoveSelectedAtomsDialog(mol, mw, preselected_atoms={0, 1})
    assert dlg.selected_atoms == {0, 1}
    dlg.close()


def test_show_atom_labels_none_positions(make_dialog):
    dlg, mol, mw = make_dialog()
    mw.view_3d_manager.atom_positions_3d = None
    with patch("moleditpy.ui.move_selected_atoms_dialog.logging.error") as mock_log:
        dlg.selected_atoms.add(0)
        dlg.show_atom_labels()
        mock_log.assert_called_once_with(
            "atom_positions_3d is None in update_atom_labels"
        )


def test_show_and_clear_atom_labels(make_dialog):
    dlg, mol, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.selected_atoms.update([0, 1])

    with patch("moleditpy.ui.move_selected_atoms_dialog.pv") as mock_pv:
        mock_pd = MagicMock()
        mock_pv.PolyData.return_value = mock_pd
        mock_pd.glyph.return_value = MagicMock()
        mock_pv.Sphere.return_value = MagicMock()
        dlg.show_atom_labels()
        plotter.add_mesh.assert_called_once()
        plotter.render.assert_called()

    dlg.clear_atom_labels()
    plotter.remove_actor.assert_called()


def test_event_filter_unrelated_obj(make_dialog):
    from PyQt6.QtWidgets import QWidget

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    event = MagicMock()
    assert dlg.eventFilter(QWidget(), event) is False


def test_event_filter_double_click(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    event = QMouseEvent(
        QEvent.Type.MouseButtonDblClick,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    dlg.drag_state["potential_drag"] = True
    assert dlg.eventFilter(plotter.interactor, event) is False
    assert dlg.drag_state["potential_drag"] is False


def test_event_filter_mouse_press_with_selection(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.selected_atoms.add(0)
    event = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is False


def test_event_filter_mouse_press_selects_atom(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    plotter.interactor.GetEventPosition.return_value = (100, 100)

    event = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    with (
        patch(
            "moleditpy.ui.move_selected_atoms_dialog.pick_atom_index_from_screen",
            return_value=0,
        ),
        patch.object(dlg, "on_atom_picked") as mock_pick,
    ):
        assert dlg.eventFilter(plotter.interactor, event) is True
        mock_pick.assert_called_once_with(0)
        assert dlg.drag_state["consume_next_left_release"] is True


def test_event_filter_mouse_press_empty_space(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    plotter.interactor.GetEventPosition.return_value = (100, 100)

    event = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    with patch(
        "moleditpy.ui.move_selected_atoms_dialog.pick_atom_index_from_screen",
        return_value=None,
    ):
        assert dlg.eventFilter(plotter.interactor, event) is False


def test_event_filter_mouse_move_hover_cursor(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    plotter.interactor.GetEventPosition.return_value = (100, 100)

    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(100.0, 100.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )

    dlg.selected_atoms.add(0)

    # Cursor over selected atom
    with patch(
        "moleditpy.ui.move_selected_atoms_dialog.pick_atom_index_from_screen",
        return_value=0,
    ):
        dlg.eventFilter(plotter.interactor, event)
        plotter.setCursor.assert_called_with(Qt.CursorShape.OpenHandCursor)

    # Cursor not over selected atom
    with patch(
        "moleditpy.ui.move_selected_atoms_dialog.pick_atom_index_from_screen",
        return_value=1,
    ):
        dlg.eventFilter(plotter.interactor, event)
        plotter.setCursor.assert_called_with(Qt.CursorShape.ArrowCursor)


def test_event_filter_mouse_release_consume(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    event = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    dlg.drag_state["consume_next_left_release"] = True
    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg.drag_state["consume_next_left_release"] is False


def test_mouse_move_potential_drag_to_actual_drag(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.drag_state["potential_drag"] = True
    dlg.drag_state["drag_start_pos"] = (100, 100)

    plotter.interactor.GetEventPosition.return_value = (110, 110)

    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(110.0, 110.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg.drag_state["is_dragging_group"] is True
    assert dlg.drag_state["potential_drag"] is False
    plotter.setCursor.assert_called_with(Qt.CursorShape.ClosedHandCursor)


def test_mouse_move_during_actual_drag(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.drag_state["is_dragging_group"] = True
    dlg.drag_state["drag_start_pos"] = (100, 100)
    dlg.drag_state["mouse_moved_during_drag"] = False

    plotter.interactor.GetEventPosition.return_value = (105, 105)
    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(105.0, 105.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg.drag_state["mouse_moved_during_drag"] is True


def test_mouse_release_no_movement_toggles_atom(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.drag_state["potential_drag"] = True
    dlg.drag_state["mouse_moved_during_drag"] = False
    dlg.clicked_atom_for_toggle = 0

    event = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    with patch.object(dlg, "on_atom_picked") as mock_pick:
        assert dlg.eventFilter(plotter.interactor, event) is True
        mock_pick.assert_called_once_with(0)
        assert dlg.drag_state["potential_drag"] is False


def test_mouse_release_with_movement_resets_drag_state(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.drag_state["is_dragging_group"] = True
    dlg.drag_state["mouse_moved_during_drag"] = True
    dlg.clicked_atom_for_toggle = 0

    event = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(110.0, 110.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    with patch.object(dlg, "on_atom_picked") as mock_pick:
        assert dlg.eventFilter(plotter.interactor, event) is True
        mock_pick.assert_not_called()
        assert dlg.drag_state["is_dragging_group"] is False


def test_handle_mouse_press_exceptions(make_dialog):
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    plotter.interactor.GetEventPosition.side_effect = AttributeError("Mock error")

    event = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100.0, 100.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is False


def test_on_atom_picked_during_drag_ignored(make_dialog):
    dlg, _, _ = make_dialog()
    dlg.drag_state["is_dragging_group"] = True
    dlg.selected_atoms.clear()
    dlg.on_atom_picked(0)
    assert 0 not in dlg.selected_atoms


def test_update_display_many_atoms(make_dialog):
    dlg, _, _ = make_dialog()
    dlg.selected_atoms.update([0, 1, 2, 3, 4, 5, 6])
    dlg.update_display()
    text = dlg.widgets["selection_label"].text()
    assert "..." in text
    assert "Selected: 7 atoms" in text
