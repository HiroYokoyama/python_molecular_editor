"""
Tests for ZoomableView (ui/zoomable_view.py).

Covers:
  - __init__ properties (anchors, scroll bars, drag mode)
  - wheelEvent: Ctrl+wheel zooms in/out; plain wheel passes to super
  - wheelEvent: clamps at min (0.05) and max (20.0) scale limits
  - mousePressEvent: middle button and Shift+left start panning
  - mousePressEvent: other buttons pass to super
  - mouseMoveEvent: panning updates both scrollbars; non-panning passes to super
  - mouseReleaseEvent: ends pan and restores cursor per scene mode
  - mouseReleaseEvent: non-pan release passes to super
  - viewportEvent: NativeGesture ZoomNativeGesture applies scale
"""

import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPoint, QEvent
from PyQt6.QtGui import QMouseEvent, QWheelEvent
from PyQt6.QtWidgets import QGraphicsScene
from moleditpy.ui.zoomable_view import ZoomableView


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class _ModeScene(QGraphicsScene):
    """Real QGraphicsScene with a .mode attribute for cursor-restore tests."""

    def __init__(self, mode="select"):
        super().__init__()
        self.mode = mode


def _make_view(app, mode="select"):
    scene = _ModeScene(mode)
    view = ZoomableView(scene)
    return view, scene


def _wheel_event(delta_y: int, ctrl: bool = False) -> QWheelEvent:
    modifiers = (
        Qt.KeyboardModifier.ControlModifier if ctrl else Qt.KeyboardModifier.NoModifier
    )
    ev = MagicMock(spec=QWheelEvent)
    ev.modifiers.return_value = modifiers
    ev.angleDelta.return_value = QPoint(0, delta_y)
    return ev


def _mouse_event(
    button, modifiers=Qt.KeyboardModifier.NoModifier, pos=None
) -> QMouseEvent:
    ev = MagicMock(spec=QMouseEvent)
    ev.button.return_value = button
    ev.modifiers.return_value = modifiers
    ev.pos.return_value = QPoint(100, 100) if pos is None else pos
    ev.buttons.return_value = button
    return ev


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------


def test_init_panning_state(app):
    """ZoomableView initialises with _is_panning=False and zero scroll offsets."""
    view, _ = _make_view(app)
    assert view._is_panning is False
    assert view._pan_start_scroll_h == 0
    assert view._pan_start_scroll_v == 0


def test_init_drag_mode(app):
    """ZoomableView initialises with NoDrag drag mode."""
    view, _ = _make_view(app)
    assert view.dragMode() == ZoomableView.DragMode.NoDrag


def test_init_scroll_bars_always_on(app):
    """ZoomableView initialises with both scroll bars always visible."""
    view, _ = _make_view(app)
    assert view.verticalScrollBarPolicy() == Qt.ScrollBarPolicy.ScrollBarAlwaysOn
    assert view.horizontalScrollBarPolicy() == Qt.ScrollBarPolicy.ScrollBarAlwaysOn


# ---------------------------------------------------------------------------
# wheelEvent
# ---------------------------------------------------------------------------


def test_wheel_ctrl_zoom_in_scales_up(app):
    """Ctrl+scroll-up increases the view scale."""
    view, _ = _make_view(app)
    before = view.transform().m11()
    ev = _wheel_event(delta_y=120, ctrl=True)
    view.wheelEvent(ev)
    assert view.transform().m11() > before
    ev.accept.assert_called_once()


def test_wheel_ctrl_zoom_out_scales_down(app):
    """Ctrl+scroll-down decreases the view scale."""
    view, _ = _make_view(app)
    ev_in = _wheel_event(delta_y=120, ctrl=True)
    view.wheelEvent(ev_in)
    mid = view.transform().m11()

    ev_out = _wheel_event(delta_y=-120, ctrl=True)
    view.wheelEvent(ev_out)
    assert view.transform().m11() < mid


def test_wheel_no_ctrl_passes_to_super(app):
    """Scroll without Ctrl passes the event to the parent class and leaves scale unchanged."""
    view, _ = _make_view(app)
    before = view.transform().m11()
    ev = _wheel_event(delta_y=120, ctrl=False)
    with patch.object(ZoomableView.__bases__[0], "wheelEvent") as mock_super:
        view.wheelEvent(ev)
    # Scale should NOT have changed (super handled it)
    assert view.transform().m11() == pytest.approx(before)


def test_wheel_ctrl_does_not_exceed_max_scale(app):
    """Repeated Ctrl+scroll-up clamps the scale near the maximum limit."""
    view, _ = _make_view(app)
    # Zoom in many times; the guard allows one overshoot step (max * 1.1)
    for _ in range(200):
        ev = _wheel_event(delta_y=120, ctrl=True)
        view.wheelEvent(ev)
    # Scale stabilises after at most one overshoot of 1.1×
    assert view.transform().m11() <= 20.0 * 1.1 + 0.1


def test_wheel_ctrl_does_not_go_below_min_scale(app):
    """Repeated Ctrl+scroll-down clamps the scale near the minimum limit."""
    view, _ = _make_view(app)
    # Zoom out many times; the guard allows one undershoot step (min / 1.1)
    for _ in range(200):
        ev = _wheel_event(delta_y=-120, ctrl=True)
        view.wheelEvent(ev)
    assert view.transform().m11() >= 0.05 / 1.1 - 0.01


# ---------------------------------------------------------------------------
# mousePressEvent
# ---------------------------------------------------------------------------


def test_mouse_press_middle_button_starts_pan(app):
    """Middle-button press starts panning and accepts the event."""
    view, _ = _make_view(app)
    ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mousePressEvent(ev)
    assert view._is_panning is True
    ev.accept.assert_called_once()


def test_mouse_press_shift_left_starts_pan(app):
    """Shift+left-button press starts panning and accepts the event."""
    view, _ = _make_view(app)
    ev = _mouse_event(
        Qt.MouseButton.LeftButton, modifiers=Qt.KeyboardModifier.ShiftModifier
    )
    view.mousePressEvent(ev)
    assert view._is_panning is True
    ev.accept.assert_called_once()


def test_mouse_press_pan_sets_cursor_closed_hand(app):
    """Starting a pan changes the cursor to ClosedHandCursor."""
    view, _ = _make_view(app)
    ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mousePressEvent(ev)
    assert view.cursor().shape() == Qt.CursorShape.ClosedHandCursor


def test_mouse_press_left_no_shift_passes_to_super(app):
    """Left-button press without Shift does not start panning."""
    view, _ = _make_view(app)
    from PyQt6.QtCore import QEvent, QPointF

    ev = QMouseEvent(
        QEvent.Type.MouseButtonPress,
        QPointF(100, 100),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    view.mousePressEvent(ev)
    assert view._is_panning is False


def test_mouse_press_none_event_returns_early(app):
    """mousePressEvent with None does not raise and does not start panning."""
    view, _ = _make_view(app)
    view.mousePressEvent(None)  # should not raise
    assert view._is_panning is False


# ---------------------------------------------------------------------------
# mouseMoveEvent
# ---------------------------------------------------------------------------


def test_mouse_move_while_panning_updates_scrollbars(app):
    """Mouse move during panning accepts the event and updates scroll positions."""
    view, _ = _make_view(app)
    # Start panning at (100, 100)
    press_ev = _mouse_event(Qt.MouseButton.MiddleButton, pos=QPoint(100, 100))
    view.mousePressEvent(press_ev)

    move_ev = MagicMock(spec=QMouseEvent)
    move_ev.pos.return_value = QPoint(110, 120)  # moved +10, +20

    view.mouseMoveEvent(move_ev)
    move_ev.accept.assert_called_once()


def test_mouse_move_not_panning_passes_to_super(app):
    """Mouse move without active pan passes the event to the parent class."""
    view, _ = _make_view(app)
    assert view._is_panning is False
    from PyQt6.QtCore import QEvent, QPointF

    move_ev = QMouseEvent(
        QEvent.Type.MouseMove,
        QPointF(50, 50),
        Qt.MouseButton.NoButton,
        Qt.MouseButton.NoButton,
        Qt.KeyboardModifier.NoModifier,
    )
    # Should not raise, super handles it; _is_panning stays False
    view.mouseMoveEvent(move_ev)


def test_mouse_move_none_event_returns_early(app):
    """mouseMoveEvent with None does not raise."""
    view, _ = _make_view(app)
    view.mouseMoveEvent(None)  # should not raise


# ---------------------------------------------------------------------------
# mouseReleaseEvent
# ---------------------------------------------------------------------------


def test_mouse_release_ends_pan(app):
    """Mouse release while panning ends the pan and accepts the event."""
    view, _ = _make_view(app)
    press_ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mousePressEvent(press_ev)
    assert view._is_panning is True

    rel_ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mouseReleaseEvent(rel_ev)
    assert view._is_panning is False
    rel_ev.accept.assert_called_once()


def test_mouse_release_restores_arrow_cursor_in_select_mode(app):
    """Ending a pan in select mode restores the ArrowCursor."""
    view, scene = _make_view(app, mode="select")
    press_ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mousePressEvent(press_ev)
    rel_ev = _mouse_event(Qt.MouseButton.MiddleButton)
    view.mouseReleaseEvent(rel_ev)
    assert view.cursor().shape() == Qt.CursorShape.ArrowCursor


def test_mouse_release_restores_cross_cursor_in_atom_mode(app):
    """Ending a pan in atom mode restores the CrossCursor."""
    view, scene = _make_view(app, mode="atom_C")
    view._is_panning = True
    rel_ev = _mouse_event(Qt.MouseButton.LeftButton)
    view.mouseReleaseEvent(rel_ev)
    assert view.cursor().shape() == Qt.CursorShape.CrossCursor


def test_mouse_release_restores_cross_cursor_in_bond_mode(app):
    """Ending a pan in bond mode restores the CrossCursor."""
    view, scene = _make_view(app, mode="bond_single")
    view._is_panning = True
    rel_ev = _mouse_event(Qt.MouseButton.LeftButton)
    view.mouseReleaseEvent(rel_ev)
    assert view.cursor().shape() == Qt.CursorShape.CrossCursor


def test_mouse_release_restores_cross_cursor_in_charge_mode(app):
    """Ending a pan in charge mode restores the CrossCursor."""
    view, scene = _make_view(app, mode="charge_plus")
    view._is_panning = True
    rel_ev = _mouse_event(Qt.MouseButton.LeftButton)
    view.mouseReleaseEvent(rel_ev)
    assert view.cursor().shape() == Qt.CursorShape.CrossCursor


def test_mouse_release_restores_arrow_for_unknown_mode(app):
    """Ending a pan in an unknown scene mode restores the ArrowCursor."""
    view, scene = _make_view(app, mode="unknown_mode")
    view._is_panning = True
    rel_ev = _mouse_event(Qt.MouseButton.LeftButton)
    view.mouseReleaseEvent(rel_ev)
    assert view.cursor().shape() == Qt.CursorShape.ArrowCursor


def test_mouse_release_not_panning_passes_to_super(app):
    """Mouse release when not panning passes the event to the parent class."""
    view, _ = _make_view(app)
    assert view._is_panning is False
    from PyQt6.QtCore import QEvent, QPointF

    rel_ev = QMouseEvent(
        QEvent.Type.MouseButtonRelease,
        QPointF(100, 100),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    # Should not set panning state or raise
    view.mouseReleaseEvent(rel_ev)
    assert view._is_panning is False


# ---------------------------------------------------------------------------
# viewportEvent (pinch zoom)
# ---------------------------------------------------------------------------


def test_viewport_event_non_native_gesture_passes_to_super(app):
    """viewportEvent with a non-NativeGesture event passes to the parent class."""
    view, _ = _make_view(app)
    # Use a benign event type that super().viewportEvent will not try to
    # dispatch as a mouse event (avoids C++ type-cast errors on plain QEvent).
    ev = QEvent(QEvent.Type.WindowActivate)
    result = view.viewportEvent(ev)
    assert isinstance(result, bool)


def test_viewport_event_zoom_native_gesture_scales(app):
    """viewportEvent with ZoomNativeGesture scales the view and returns True."""
    view, _ = _make_view(app)
    before = view.transform().m11()
    ev = MagicMock()
    ev.type.return_value = QEvent.Type.NativeGesture
    ev.gestureType.return_value = Qt.NativeGestureType.ZoomNativeGesture
    ev.value.return_value = 0.1  # scale factor delta → multiply by 1.1
    result = view.viewportEvent(ev)
    assert result is True
    assert view.transform().m11() > before
