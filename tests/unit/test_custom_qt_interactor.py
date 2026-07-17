"""Unit tests for CustomQtInteractor fast-click handling."""

from unittest.mock import MagicMock, patch

from PyQt6.QtCore import QEvent, QPointF, Qt
from PyQt6.QtGui import QMouseEvent

from moleditpy.ui.custom_qt_interactor import CustomQtInteractor

_BASE = CustomQtInteractor.__mro__[1]


def _make_widget():
    return CustomQtInteractor(main_window=None)


def _double_click_event():
    return QMouseEvent(
        QEvent.Type.MouseButtonDblClick,
        QPointF(10.0, 10.0),
        QPointF(10.0, 10.0),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )


def test_double_click_redispatched_as_plain_press(app):
    widget = _make_widget()
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        widget.mouseDoubleClickEvent(_double_click_event())

        assert mock_press.call_count == 1
        forwarded = mock_press.call_args[0][0]
        assert forwarded.type() == QEvent.Type.MouseButtonPress
        assert forwarded.button() == Qt.MouseButton.LeftButton


def test_fast_consecutive_presses_all_pass_through(app):
    widget = _make_widget()
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        for _ in range(4):
            widget.mousePressEvent(MagicMock())

        assert mock_press.call_count == 4


def test_release_passes_through_unconditionally(app):
    widget = _make_widget()
    with patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release:
        widget.mouseReleaseEvent(MagicMock())
        widget.mouseReleaseEvent(MagicMock())

        assert mock_release.call_count == 2

def _mw_with_view2d():
    mw = MagicMock()
    return mw, mw.init_manager.view_2d


def test_init_stores_main_window(app):
    mw = MagicMock()
    widget = CustomQtInteractor(main_window=mw)
    assert widget.main_window is mw


def test_init_default_main_window_is_none(app):
    assert _make_widget().main_window is None


def test_wheel_zooms_then_returns_focus_to_2d_view(app):
    mw, view_2d = _mw_with_view2d()
    widget = CustomQtInteractor(main_window=mw)
    with patch.object(_BASE, "wheelEvent", create=True) as mock_wheel:
        event = MagicMock()
        widget.wheelEvent(event)

        mock_wheel.assert_called_once_with(event)
        view_2d.setFocus.assert_called_once_with()


def test_wheel_without_main_window_still_zooms(app):
    widget = _make_widget()
    with patch.object(_BASE, "wheelEvent", create=True) as mock_wheel:
        widget.wheelEvent(MagicMock())

        assert mock_wheel.call_count == 1


def test_wheel_without_view_2d_does_not_crash(app):
    mw = MagicMock()
    mw.init_manager = object()  # no view_2d attribute
    widget = CustomQtInteractor(main_window=mw)
    with patch.object(_BASE, "wheelEvent", create=True):
        widget.wheelEvent(MagicMock())


def test_release_returns_focus_to_2d_view(app):
    mw, view_2d = _mw_with_view2d()
    widget = CustomQtInteractor(main_window=mw)
    with patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release:
        event = MagicMock()
        widget.mouseReleaseEvent(event)

        mock_release.assert_called_once_with(event)
        view_2d.setFocus.assert_called_once_with()


def test_release_without_main_window_still_forwards(app):
    widget = _make_widget()
    with patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release:
        widget.mouseReleaseEvent(MagicMock())

        assert mock_release.call_count == 1


def test_release_without_view_2d_does_not_crash(app):
    mw = MagicMock()
    mw.init_manager = object()
    widget = CustomQtInteractor(main_window=mw)
    with patch.object(_BASE, "mouseReleaseEvent", create=True):
        widget.mouseReleaseEvent(MagicMock())


def test_synthetic_press_preserves_button_position_modifiers(app):
    widget = _make_widget()
    event = QMouseEvent(
        QEvent.Type.MouseButtonDblClick,
        QPointF(12.5, 20.25),
        QPointF(112.5, 120.25),
        Qt.MouseButton.RightButton,
        Qt.MouseButton.RightButton,
        Qt.KeyboardModifier.ShiftModifier,
    )
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        widget.mouseDoubleClickEvent(event)

        forwarded = mock_press.call_args[0][0]
        assert forwarded.type() == QEvent.Type.MouseButtonPress
        assert forwarded.position() == QPointF(12.5, 20.25)
        assert forwarded.globalPosition() == QPointF(112.5, 120.25)
        assert forwarded.button() == Qt.MouseButton.RightButton
        assert forwarded.buttons() == Qt.MouseButton.RightButton
        assert forwarded.modifiers() == Qt.KeyboardModifier.ShiftModifier


def test_double_click_accepts_event(app):
    widget = _make_widget()
    event = _double_click_event()
    with patch.object(_BASE, "mousePressEvent", create=True):
        widget.mouseDoubleClickEvent(event)

    assert event.isAccepted()


def test_double_click_broken_event_suppressed(app):
    widget = _make_widget()
    event = MagicMock()
    event.position.side_effect = AttributeError("deleted C++ object")
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        widget.mouseDoubleClickEvent(event)  # must not raise

        mock_press.assert_not_called()
        event.accept.assert_not_called()


def test_double_click_bad_event_types_suppressed(app):
    widget = _make_widget()
    # MagicMock return values make the QMouseEvent constructor raise TypeError
    event = MagicMock()
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        widget.mouseDoubleClickEvent(event)  # must not raise

        mock_press.assert_not_called()
        event.accept.assert_not_called()
