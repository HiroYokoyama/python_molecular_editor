"""Unit tests for CustomQtInteractor fast-click swallow pairing."""

from unittest.mock import MagicMock, patch

from PyQt6.QtCore import Qt

from moleditpy.ui.custom_qt_interactor import CustomQtInteractor

_BASE = CustomQtInteractor.__mro__[1]


def _make_widget():
    widget = CustomQtInteractor(main_window=None)
    widget._click_count = 0
    widget._last_click_time = 0.0
    widget._swallowed_buttons = set()
    return widget


def _event(button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.button.return_value = button
    return event


def test_third_fast_left_click_and_its_release_are_swallowed(app):
    widget = _make_widget()
    with (
        patch.object(_BASE, "mousePressEvent", create=True) as mock_press,
        patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release,
    ):
        widget.mousePressEvent(_event())
        widget.mousePressEvent(_event())
        third = _event()
        widget.mousePressEvent(third)

        assert mock_press.call_count == 2
        third.accept.assert_called_once()
        assert Qt.MouseButton.LeftButton in widget._swallowed_buttons

        # The paired release is swallowed exactly once
        widget.mouseReleaseEvent(_event())
        mock_release.assert_not_called()
        assert Qt.MouseButton.LeftButton not in widget._swallowed_buttons

        # The next release passes through again
        widget.mouseReleaseEvent(_event())
        assert mock_release.call_count == 1


def test_double_click_swallows_press_release_pair(app):
    widget = _make_widget()
    with patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release:
        event = _event()
        widget.mouseDoubleClickEvent(event)

        event.accept.assert_called_once()
        assert widget._click_count == 2
        assert Qt.MouseButton.LeftButton in widget._swallowed_buttons

        widget.mouseReleaseEvent(_event())
        mock_release.assert_not_called()
        assert not widget._swallowed_buttons


def test_right_click_not_counted_toward_triple_click_filter(app):
    widget = _make_widget()
    with patch.object(_BASE, "mousePressEvent", create=True) as mock_press:
        widget.mousePressEvent(_event())
        widget.mousePressEvent(_event())
        widget.mousePressEvent(_event(Qt.MouseButton.RightButton))

        assert mock_press.call_count == 3
        assert not widget._swallowed_buttons


def test_release_of_other_button_passes_while_left_latched(app):
    widget = _make_widget()
    widget._swallowed_buttons = {Qt.MouseButton.LeftButton}
    with patch.object(_BASE, "mouseReleaseEvent", create=True) as mock_release:
        widget.mouseReleaseEvent(_event(Qt.MouseButton.RightButton))

        assert mock_release.call_count == 1
        assert Qt.MouseButton.LeftButton in widget._swallowed_buttons


def test_leave_event_clears_stale_latches_when_no_button_held(app):
    widget = _make_widget()
    widget._swallowed_buttons = {Qt.MouseButton.LeftButton}
    with (
        patch.object(_BASE, "leaveEvent", create=True) as mock_leave,
        patch("moleditpy.ui.custom_qt_interactor.QApplication") as mock_qapp,
    ):
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.NoButton
        widget.leaveEvent(MagicMock())

        assert not widget._swallowed_buttons
        mock_leave.assert_called_once()


def test_leave_event_keeps_latches_while_button_held(app):
    widget = _make_widget()
    widget._swallowed_buttons = {Qt.MouseButton.LeftButton}
    with (
        patch.object(_BASE, "leaveEvent", create=True),
        patch("moleditpy.ui.custom_qt_interactor.QApplication") as mock_qapp,
    ):
        mock_qapp.mouseButtons.return_value = Qt.MouseButton.LeftButton
        widget.leaveEvent(MagicMock())

        assert Qt.MouseButton.LeftButton in widget._swallowed_buttons
