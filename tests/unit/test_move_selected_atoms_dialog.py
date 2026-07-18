from PyQt6.QtCore import QEvent, QPointF, Qt
from PyQt6.QtGui import QMouseEvent

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
        """apply_translation shows a warning when no atoms are selected."""
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        """apply_translation shows a warning when a translation input is non-numeric."""
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_trans_input.setText("bad")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_translation()
        mb.warning.assert_called_once()

    def test_translation_updates_only_selected_atoms(self, make_dialog):
        """apply_translation moves only selected atoms and leaves others unchanged."""
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
        """apply_rotation shows a warning when no atoms are selected."""
        dlg, _, _ = make_dialog()
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_invalid_input_shows_warning(self, make_dialog):
        """apply_rotation shows a warning when a rotation input is non-numeric."""
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.x_rot_input.setText("not_a_number")
        with patch("moleditpy.ui.move_selected_atoms_dialog.QMessageBox") as mb:
            dlg.apply_rotation()
        mb.warning.assert_called_once()

    def test_rotation_updates_only_selected_atoms_around_centroid(self, make_dialog):
        """apply_rotation rotates only selected atoms around their centroid."""
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
        """reset_translation_inputs sets all three translation fields to '0.0'."""
        dlg, _, _ = make_dialog()
        dlg.x_trans_input.setText("5.5")
        dlg.y_trans_input.setText("-3.0")
        dlg.z_trans_input.setText("1.2")
        dlg.reset_translation_inputs()
        assert dlg.x_trans_input.text() == "0.0"
        assert dlg.y_trans_input.text() == "0.0"
        assert dlg.z_trans_input.text() == "0.0"

    def test_reset_rotation_inputs(self, make_dialog):
        """reset_rotation_inputs sets all three rotation fields to '0.0'."""
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
        """clear_selection empties the selected_atoms set."""
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.update([0, 1])
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0


# ---------------------------------------------------------------------------
# Extra tests for Properties, Highlights, and Event Filter
# ---------------------------------------------------------------------------


def test_group_atoms_property(make_dialog):
    """group_atoms property aliases selected_atoms for compatibility with the mixin."""
    dlg, _, _ = make_dialog()
    dlg.group_atoms = {0, 2}
    assert dlg.selected_atoms == {0, 2}
    assert dlg.group_atoms == {0, 2}


def test_preselected_atoms_init(qapp):
    """MoveSelectedAtomsDialog pre-populates selected_atoms from preselected_atoms."""
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
    """show_atom_labels logs a warning and does nothing when atom_positions_3d is None."""
    dlg, mol, mw = make_dialog()
    mw.view_3d_manager.atom_positions_3d = None
    with patch("moleditpy.ui.move_selected_atoms_dialog.logging.warning") as mock_log:
        dlg.selected_atoms.add(0)
        dlg.show_atom_labels()
        mock_log.assert_called_once_with(
            "atom_positions_3d is None in update_atom_labels"
        )


def test_show_and_clear_atom_labels(make_dialog):
    """show_atom_labels calls add_mesh; clear_atom_labels removes the actor."""
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
    """eventFilter returns False for events on objects other than the VTK interactor."""
    from PyQt6.QtWidgets import QWidget

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter
    event = MagicMock()
    assert dlg.eventFilter(QWidget(), event) is False


def test_event_filter_double_click(make_dialog):
    """A double-click event resets potential_drag and returns False."""
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

    dlg.potential_drag = True
    assert dlg.eventFilter(plotter.interactor, event) is False
    assert dlg.potential_drag is False


def test_event_filter_mouse_press_with_selection(make_dialog):
    """A left-press with atoms already selected delegates to VTK (returns False)."""
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
    """A left-press that hits an atom picks it and returns True."""
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
        assert dlg._consume_next_left_release is True


def test_event_filter_mouse_press_empty_space(make_dialog):
    """A left-press that misses all atoms returns False."""
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
    """Mouse-move over a selected atom changes the cursor to OpenHand; otherwise to Arrow."""
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
    """A mouse-release is consumed and clears _consume_next_left_release when the flag is set."""
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

    dlg._consume_next_left_release = True
    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg._consume_next_left_release is False


def test_mouse_move_potential_drag_to_actual_drag(make_dialog):
    """A mouse-move beyond the threshold converts a potential drag to an actual drag."""
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.potential_drag = True
    dlg.drag_start_pos = (100, 100)

    plotter.interactor.GetEventPosition.return_value = (115, 115)

    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(110.0, 110.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg.is_dragging_group is True
    assert dlg.potential_drag is False
    plotter.setCursor.assert_called_with(Qt.CursorShape.ClosedHandCursor)


def test_mouse_move_during_actual_drag(make_dialog):
    """Mouse-move while is_dragging_group is True sets mouse_moved_during_drag."""
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.is_dragging_group = True
    dlg.drag_start_pos = (100, 100)
    dlg.mouse_moved_during_drag = False

    plotter.interactor.GetEventPosition.return_value = (110, 110)
    event = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(105.0, 105.0),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )

    assert dlg.eventFilter(plotter.interactor, event) is True
    assert dlg.mouse_moved_during_drag is True


def test_mouse_release_no_movement_toggles_atom(make_dialog):
    """A release with no movement calls on_atom_picked to toggle the atom selection."""
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.potential_drag = True
    dlg.mouse_moved_during_drag = False
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
        assert dlg.potential_drag is False


def test_mouse_release_with_movement_resets_drag_state(make_dialog):
    """A release after dragging resets drag state without calling on_atom_picked."""
    from PyQt6.QtCore import QEvent, Qt, QPointF
    from PyQt6.QtGui import QMouseEvent

    dlg, _, mw = make_dialog()
    plotter = MagicMock()
    mw.view_3d_manager.plotter = plotter

    dlg.is_dragging_group = True
    dlg.drag_start_pos = (100, 100)
    dlg.mouse_moved_during_drag = True
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
        assert dlg.is_dragging_group is False


def test_handle_mouse_press_exceptions(make_dialog):
    """eventFilter returns False and does not raise when GetEventPosition throws."""
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
    """on_atom_picked is a no-op while is_dragging_group is True."""
    dlg, _, _ = make_dialog()
    dlg.is_dragging_group = True
    dlg.selected_atoms.clear()
    dlg.on_atom_picked(0)
    assert 0 not in dlg.selected_atoms


def test_update_display_many_atoms(make_dialog):
    """update_display shows '...' and total count when more than 5 atoms are selected."""
    dlg, _, _ = make_dialog()
    dlg.selected_atoms.update([0, 1, 2, 3, 4, 5, 6])
    dlg.update_display()
    text = dlg.widgets["selection_label"].text()
    assert "..." in text
    assert "Selected: 7 atoms" in text


# ---------------------------------------------------------------------------
# click-to-deselect / eventFilter
# ---------------------------------------------------------------------------


class TestClickToDeselect:
    def test_eventfilter_delegates_to_vtk_when_atoms_selected(self, make_dialog):
        """eventFilter returns False to let VTK handle the press when atoms are selected."""
        dlg, mol, mw = make_dialog()
        dlg.selected_atoms.add(0)

        ev = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(100.0, 100.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )

        # eventFilter returns False to delegate to VTK
        result = dlg.eventFilter(mw.view_3d_manager.plotter.interactor, ev)
        assert result is False

    def test_eventfilter_handles_press_when_no_atoms_selected(self, make_dialog):
        """A left-press with no prior selection picks the clicked atom and returns True."""
        dlg, mol, mw = make_dialog()
        dlg.selected_atoms.clear()

        ev = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(100.0, 100.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        mw.view_3d_manager.plotter.interactor.GetEventPosition.return_value = (100, 100)

        with (
            patch(
                "moleditpy.ui.move_selected_atoms_dialog.pick_atom_index_from_screen",
                return_value=0,
            ) as mock_pick,
            patch.object(dlg, "show_atom_labels"),
            patch.object(dlg, "clear_atom_labels"),
        ):
            try:
                result = dlg.eventFilter(mw.view_3d_manager.plotter.interactor, ev)
                print("TEST DEBUG - result:", result)
                print("TEST DEBUG - call count:", mock_pick.call_count)
            except Exception as e:
                print("TEST DEBUG - exception:", e)
                import traceback

                traceback.print_exc()
                raise e
            assert result is True
            assert 0 in dlg.selected_atoms

    def test_eventfilter_release_click_only_deselects_atom(self, make_dialog):
        """A release after a pure click (no drag) calls on_atom_picked to toggle the atom."""
        dlg, mol, mw = make_dialog()
        dlg.selected_atoms.add(0)
        dlg.clicked_atom_for_toggle = 0
        dlg.is_dragging_group = False
        dlg.potential_drag = True
        dlg.drag_start_pos = (100, 100)
        dlg.mouse_moved_during_drag = False

        ev = QMouseEvent(
            QEvent.Type.MouseButtonRelease,
            QPointF(100.0, 100.0),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )

        with patch.object(dlg, "on_atom_picked") as mock_on_pick:
            result = dlg.eventFilter(mw.view_3d_manager.plotter.interactor, ev)
            assert result is True
            mock_on_pick.assert_called_once_with(0)


def test_show_atom_labels_camera_restore(qapp):
    """show_atom_labels restores the camera position after adding highlight meshes."""
    from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

    mol = _ethane()
    mw = _make_main_window(mol)
    plotter = MagicMock()
    plotter.camera_position = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
    mw.view_3d_manager.plotter = plotter

    with (
        patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
        patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
    ):
        dlg = MoveSelectedAtomsDialog(mol, mw)

    dlg.selected_atoms.add(0)

    with (
        patch("moleditpy.ui.move_selected_atoms_dialog.pv") as mock_pv,
        patch.object(dlg, "clear_atom_labels") as mock_clear,
    ):
        mock_pd = MagicMock()
        mock_pv.PolyData.return_value = mock_pd
        mock_pd.glyph.return_value = MagicMock()

        MoveSelectedAtomsDialog.show_atom_labels(dlg)

        mock_clear.assert_called_once()

        args, kwargs = plotter.add_mesh.call_args
        assert kwargs.get("reset_camera") is False

        assert plotter.camera_position == [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
        plotter.render.assert_called()


def test_box_selection_toggle_button(qapp):
    """Test that the Box Selection toggle enables/disables rectangle picking."""
    from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

    mol = _ethane()
    mw = _make_main_window(mol)

    with (
        patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
        patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
    ):
        dlg = MoveSelectedAtomsDialog(mol, mw)
        plotter = mw.view_3d_manager.plotter
        plotter.enable_rectangle_picking = MagicMock()
        plotter.disable_picking = MagicMock()
        plotter.interactor = MagicMock()
        plotter.interactor.GetInteractorStyle.return_value = "mock_style"

        btn = dlg.widgets["box_select_btn"]

        # Turn ON
        btn.setChecked(True)
        dlg.toggle_box_selection(True)
        plotter.enable_rectangle_picking.assert_called_once()
        assert dlg.original_style == "mock_style"

        # Turn OFF
        btn.setChecked(False)
        dlg.toggle_box_selection(False)
        plotter.disable_picking.assert_called_once()


def test_box_selection_click_to_reset(qapp):
    """Test that a single click clears the selection when Box Selection is ON."""
    from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog
    from PyQt6.QtCore import QEvent, QPoint, Qt
    from PyQt6.QtGui import QMouseEvent

    mol = _ethane()
    mw = _make_main_window(mol)

    with (
        patch.object(MoveSelectedAtomsDialog, "show_atom_labels"),
        patch.object(MoveSelectedAtomsDialog, "clear_atom_labels"),
    ):
        dlg = MoveSelectedAtomsDialog(mol, mw)
        dlg.selected_atoms.add(0)

        btn = dlg.widgets["box_select_btn"]
        btn.setChecked(True)
        # We don't actually need to call toggle_box_selection(True)
        # since we are manually passing the event to eventFilter.

        class MockEvent:
            def __init__(self, t, p):
                self._t = t
                self._p = p

            def type(self):
                return self._t

            def pos(self):
                return self._p

        press_ev = MockEvent(QEvent.Type.MouseButtonPress, QPoint(100, 100))
        release_ev = MockEvent(QEvent.Type.MouseButtonRelease, QPoint(100, 100))

        # interactor must match the one checked inside eventFilter
        interactor = mw.view_3d_manager.plotter.interactor

        with patch.object(dlg, "clear_selection") as mock_clear:
            dlg.eventFilter(interactor, press_ev)
            dlg.eventFilter(interactor, release_ev)
            mock_clear.assert_called_once()


# ---------------------------------------------------------------------------
# Box selection style restore (fast-click freeze root cause)
# ---------------------------------------------------------------------------


def test_box_selection_off_restores_style_via_pyvista(make_dialog):
    """Restore must go through iren.style so pyvista bookkeeping tracks it."""
    dlg, _mol, mw = make_dialog()
    sentinel = object()
    dlg.original_style = sentinel
    dlg.widgets["box_select_btn"] = MagicMock()
    plotter = mw.view_3d_manager.plotter

    dlg.toggle_box_selection(False)

    plotter.disable_picking.assert_called_once()
    assert plotter.iren.style is sentinel
    plotter.interactor.SetInteractorStyle.assert_not_called()


# ---------------------------------------------------------------------------
# on_rectangle_picked — PyVista box-selection callback
# ---------------------------------------------------------------------------


class _Selection:
    def __init__(self, viewport):
        self.viewport = viewport


def test_on_rectangle_picked_ignores_object_without_viewport(make_dialog):
    dlg, _mol, _mw = make_dialog()
    dlg.selected_atoms.add(3)
    dlg.on_rectangle_picked(object())  # no .viewport attribute
    assert dlg.selected_atoms == {3}  # untouched


def test_on_rectangle_picked_small_box_clears_selection(make_dialog):
    dlg, _mol, _mw = make_dialog()
    dlg.selected_atoms.update({0, 1})
    with patch.object(type(dlg), "clear_selection") as clear:
        # 10x10 box is below the 15px click threshold
        dlg.on_rectangle_picked(_Selection((100, 100, 110, 110)))
    clear.assert_called_once()


def test_on_rectangle_picked_selects_atoms_inside_box(make_dialog):
    dlg, _mol, mw = make_dialog()
    renderer = mw.view_3d_manager.plotter.renderer
    n_atoms = len(mw.view_3d_manager.atom_positions_3d)
    # Atom 0 maps inside the box; all others map far outside it.
    display_points = [(50.0, 50.0, 0.0)] + [(9999.0, 9999.0, 0.0)] * (n_atoms - 1)
    renderer.GetDisplayPoint.side_effect = display_points

    with (
        patch.object(type(dlg), "show_atom_labels") as show,
        patch.object(type(dlg), "update_display") as upd,
    ):
        dlg.on_rectangle_picked(_Selection((0, 0, 100, 100)))

    assert dlg.selected_atoms == {0}
    show.assert_called_once()
    upd.assert_called_once()


def test_on_rectangle_picked_no_atoms_inside_box_leaves_selection(make_dialog):
    dlg, _mol, mw = make_dialog()
    renderer = mw.view_3d_manager.plotter.renderer
    n_atoms = len(mw.view_3d_manager.atom_positions_3d)
    renderer.GetDisplayPoint.side_effect = [(9999.0, 9999.0, 0.0)] * n_atoms

    with (
        patch.object(type(dlg), "show_atom_labels") as show,
        patch.object(type(dlg), "update_display"),
    ):
        dlg.on_rectangle_picked(_Selection((0, 0, 100, 100)))

    assert dlg.selected_atoms == set()
    show.assert_not_called()  # nothing added -> no relabel


def test_on_rectangle_picked_none_plotter_guard(make_dialog):
    dlg, _mol, mw = make_dialog()
    mw.view_3d_manager.plotter = None
    dlg.on_rectangle_picked(_Selection((0, 0, 100, 100)))  # must not raise
    assert dlg.selected_atoms == set()


# ---------------------------------------------------------------------------
# reject — restore interactor when closing mid box-selection
# ---------------------------------------------------------------------------


def test_reject_turns_off_active_box_selection(make_dialog):
    dlg, _mol, _mw = make_dialog()
    btn = MagicMock()
    btn.isChecked.return_value = True
    dlg.widgets["box_select_btn"] = btn

    with patch.object(type(dlg), "toggle_box_selection") as toggle:
        dlg.reject()

    btn.setChecked.assert_called_once_with(False)
    toggle.assert_called_once_with(False)


def test_reject_when_box_selection_inactive(make_dialog):
    dlg, _mol, _mw = make_dialog()
    btn = MagicMock()
    btn.isChecked.return_value = False
    dlg.widgets["box_select_btn"] = btn

    with patch.object(type(dlg), "toggle_box_selection") as toggle:
        dlg.reject()

    toggle.assert_not_called()
