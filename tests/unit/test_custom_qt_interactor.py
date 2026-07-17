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
