"""
Tests for Dialog3DPickingMixin (ui/dialog_3d_picking_mixin.py).

Covers:
  - __init__ default state
  - eventFilter: None event → False; plotter None → False; mol None → False
  - eventFilter: non-interactor object → False
  - eventFilter: left-click on atom → calls on_atom_picked, returns True
  - eventFilter: left-click miss (no atom actor) → returns False
  - eventFilter: mouse-move tracking (sets _mouse_moved)
  - eventFilter: mouse-release on background (no drag) → calls clear_selection
  - eventFilter: mouse-release after drag → no clear_selection
  - enable_picking / disable_picking: installs and removes event filter
  - clear_atom_labels: removes actors and empties list
  - add_selection_label: calls plotter.add_point_labels, appends to list
  - show_atom_labels_for: clears then adds multiple labels
"""

import numpy as np
import pytest
from unittest.mock import MagicMock, call, patch
from PyQt6.QtCore import QEvent, Qt, QPoint
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import QApplication, QDialog

from moleditpy.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin


# ---------------------------------------------------------------------------
# Minimal concrete class
# ---------------------------------------------------------------------------


class _PickingDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, main_window):
        QDialog.__init__(self)
        Dialog3DPickingMixin.__init__(self)
        self.main_window = main_window
        self.mol = MagicMock()
        self.mol.GetNumAtoms.return_value = 4
        self._picked = []

    def on_atom_picked(self, atom_idx: int) -> None:
        self._picked.append(atom_idx)


def _make_mw():
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.5, 0.0, 0.0],
        ]
    )
    plotter = MagicMock()
    plotter.interactor = MagicMock()
    mw.view_3d_manager.plotter = plotter
    return mw


def _make_dlg(app, mw=None):
    if mw is None:
        mw = _make_mw()
    return _PickingDialog(mw), mw


def _left_press_event(pos=None):
    ev = MagicMock(spec=QMouseEvent)
    ev.type.return_value = QEvent.Type.MouseButtonPress
    ev.button.return_value = Qt.MouseButton.LeftButton
    ev.pos.return_value = QPoint(10, 10) if pos is None else pos
    return ev


def _move_event(pos):
    ev = MagicMock(spec=QMouseEvent)
    ev.type.return_value = QEvent.Type.MouseMove
    ev.button.return_value = Qt.MouseButton.NoButton
    ev.pos.return_value = pos
    return ev


def _release_event():
    ev = MagicMock(spec=QMouseEvent)
    ev.type.return_value = QEvent.Type.MouseButtonRelease
    ev.button.return_value = Qt.MouseButton.LeftButton
    ev.pos.return_value = QPoint(10, 10)
    return ev


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------


def test_init_defaults(app):
    dlg, _ = _make_dlg(app)
    assert dlg.picking_enabled is False
    assert dlg._mouse_press_pos is None
    assert dlg._mouse_moved is False
    assert dlg.selection_labels == []


# ---------------------------------------------------------------------------
# eventFilter — early-return paths
# ---------------------------------------------------------------------------


def test_eventfilter_none_event_returns_false(app):
    dlg, _ = _make_dlg(app)
    assert dlg.eventFilter(MagicMock(), None) is False


def test_eventfilter_plotter_none_returns_false(app):
    dlg, mw = _make_dlg(app)
    mw.view_3d_manager.plotter = None
    assert dlg.eventFilter(MagicMock(), _left_press_event()) is False


def test_eventfilter_mol_none_returns_false(app):
    dlg, _ = _make_dlg(app)
    dlg.mol = None
    assert dlg.eventFilter(MagicMock(), _left_press_event()) is False


def test_eventfilter_non_interactor_returns_false(app):
    dlg, mw = _make_dlg(app)
    other_obj = MagicMock()
    assert dlg.eventFilter(other_obj, _left_press_event()) is False


# ---------------------------------------------------------------------------
# eventFilter — atom click hit
# ---------------------------------------------------------------------------


def test_eventfilter_atom_click_calls_on_atom_picked(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter

    # Atom actor matches
    mw.view_3d_manager.atom_actor = MagicMock()
    plotter.picker.GetActor.return_value = mw.view_3d_manager.atom_actor
    plotter.picker.GetPickPosition.return_value = (0.0, 0.0, 0.0)
    plotter.interactor.GetEventPosition.return_value = (10, 10)

    # atom 0 at [0,0,0] — distance = 0 < threshold
    atom = MagicMock()
    atom.GetAtomicNum.return_value = 6
    dlg.mol.GetAtomWithIdx.return_value = atom

    result = dlg.eventFilter(plotter.interactor, _left_press_event())
    assert result is True
    assert 0 in dlg._picked


def test_eventfilter_atom_click_miss_returns_false(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter

    # Actor does NOT match atom_actor → background click
    mw.view_3d_manager.atom_actor = MagicMock()
    plotter.picker.GetActor.return_value = MagicMock()  # different object

    result = dlg.eventFilter(plotter.interactor, _left_press_event())
    assert result is False
    assert dlg._picked == []


# ---------------------------------------------------------------------------
# eventFilter — mouse move tracking
# ---------------------------------------------------------------------------


def test_eventfilter_mouse_move_sets_moved_flag(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg._mouse_press_pos = QPoint(0, 0)
    dlg._mouse_moved = False

    # Move far enough (manhattan length > 3)
    ev = _move_event(QPoint(10, 10))
    dlg.eventFilter(plotter.interactor, ev)
    assert dlg._mouse_moved is True


def test_eventfilter_mouse_move_small_does_not_set_flag(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg._mouse_press_pos = QPoint(0, 0)
    dlg._mouse_moved = False

    ev = _move_event(QPoint(1, 1))  # manhattan = 2, below threshold of 3
    dlg.eventFilter(plotter.interactor, ev)
    assert dlg._mouse_moved is False


# ---------------------------------------------------------------------------
# eventFilter — mouse release clears selection on pure click
# ---------------------------------------------------------------------------


def test_eventfilter_release_pure_click_calls_clear_selection(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg._mouse_press_pos = QPoint(0, 0)
    dlg._mouse_moved = False
    dlg.clear_selection = MagicMock()

    dlg.eventFilter(plotter.interactor, _release_event())
    dlg.clear_selection.assert_called_once()
    assert dlg._mouse_press_pos is None


def test_eventfilter_release_after_drag_no_clear_selection(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg._mouse_press_pos = QPoint(0, 0)
    dlg._mouse_moved = True
    dlg.clear_selection = MagicMock()

    dlg.eventFilter(plotter.interactor, _release_event())
    dlg.clear_selection.assert_not_called()
    assert dlg._mouse_press_pos is None


# ---------------------------------------------------------------------------
# enable_picking / disable_picking
# ---------------------------------------------------------------------------


def test_enable_picking_installs_event_filter(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg.enable_picking()
    plotter.interactor.installEventFilter.assert_called_once_with(dlg)
    assert dlg.picking_enabled is True


def test_enable_picking_none_plotter_no_crash(app):
    dlg, mw = _make_dlg(app)
    mw.view_3d_manager.plotter = None
    dlg.enable_picking()  # should not raise
    assert dlg.picking_enabled is False


def test_disable_picking_removes_event_filter(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg.picking_enabled = True
    dlg.disable_picking()
    plotter.interactor.removeEventFilter.assert_called_once_with(dlg)
    assert dlg.picking_enabled is False


def test_disable_picking_when_not_enabled_is_noop(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg.picking_enabled = False
    dlg.disable_picking()
    plotter.interactor.removeEventFilter.assert_not_called()


# ---------------------------------------------------------------------------
# clear_atom_labels
# ---------------------------------------------------------------------------


def test_clear_atom_labels_removes_actors(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    actor1, actor2 = MagicMock(), MagicMock()
    dlg.selection_labels = [actor1, actor2]

    dlg.clear_atom_labels()

    plotter.remove_actor.assert_any_call(actor1)
    plotter.remove_actor.assert_any_call(actor2)
    assert dlg.selection_labels == []


def test_clear_atom_labels_none_plotter_empties_list(app):
    dlg, mw = _make_dlg(app)
    mw.view_3d_manager.plotter = None
    dlg.selection_labels = [MagicMock()]
    dlg.clear_atom_labels()
    assert dlg.selection_labels == []


# ---------------------------------------------------------------------------
# add_selection_label
# ---------------------------------------------------------------------------


def test_add_selection_label_calls_add_point_labels(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    plotter.add_point_labels.return_value = MagicMock()

    dlg.add_selection_label(0, "A1")

    plotter.add_point_labels.assert_called_once()
    assert len(dlg.selection_labels) == 1


def test_add_selection_label_none_plotter_no_crash(app):
    dlg, mw = _make_dlg(app)
    mw.view_3d_manager.plotter = None
    dlg.add_selection_label(0, "A1")  # should not raise
    assert dlg.selection_labels == []


def test_add_selection_label_none_positions_no_crash(app):
    dlg, mw = _make_dlg(app)
    mw.view_3d_manager.atom_positions_3d = None
    dlg.add_selection_label(0, "A1")
    assert dlg.selection_labels == []


# ---------------------------------------------------------------------------
# show_atom_labels_for
# ---------------------------------------------------------------------------


def test_show_atom_labels_for_clears_then_adds(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    plotter.add_point_labels.return_value = MagicMock()

    # Seed existing labels
    dlg.selection_labels = [MagicMock()]

    dlg.show_atom_labels_for([(0, "A"), (1, "B")])

    # Should have cleared then added 2 new ones
    plotter.remove_actor.assert_called()
    assert len(dlg.selection_labels) == 2


def test_show_atom_labels_for_empty_list_clears_all(app):
    dlg, mw = _make_dlg(app)
    plotter = mw.view_3d_manager.plotter
    dlg.selection_labels = [MagicMock()]

    dlg.show_atom_labels_for([])
    assert dlg.selection_labels == []
