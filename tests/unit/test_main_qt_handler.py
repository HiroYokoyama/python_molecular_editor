"""Unit tests for the Qt message handler (_qt_message_handler) in moleditpy.main."""

import logging
import sys
import os
from unittest.mock import MagicMock, patch

import pytest

_SRC = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from PyQt6.QtCore import QtMsgType

from moleditpy.main import _DOWNGRADED_QT_PATTERNS, _qt_message_handler


@pytest.mark.parametrize("pattern", _DOWNGRADED_QT_PATTERNS)
def test_downgraded_pattern_routes_to_debug(pattern):
    msg = f"Qt: {pattern} details"
    with patch("logging.debug") as mock_debug, patch("logging.log") as mock_log:
        _qt_message_handler(QtMsgType.QtWarningMsg, MagicMock(), msg)
    mock_debug.assert_called_once_with("Qt: %s", msg)
    mock_log.assert_not_called()


@pytest.mark.parametrize(
    "mode, expected_level",
    [
        (QtMsgType.QtDebugMsg, logging.DEBUG),
        (QtMsgType.QtInfoMsg, logging.INFO),
        (QtMsgType.QtWarningMsg, logging.WARNING),
        (QtMsgType.QtCriticalMsg, logging.ERROR),
        (QtMsgType.QtFatalMsg, logging.CRITICAL),
    ],
)
def test_known_mode_maps_to_correct_level(mode, expected_level):
    msg = "some Qt message"
    with patch("logging.log") as mock_log, patch("logging.debug") as mock_debug:
        _qt_message_handler(mode, MagicMock(), msg)
    mock_log.assert_called_once_with(expected_level, "Qt: %s", msg)
    mock_debug.assert_not_called()


def test_unknown_mode_falls_back_to_warning():
    msg = "unknown mode message"
    with patch("logging.log") as mock_log:
        _qt_message_handler(999, MagicMock(), msg)
    mock_log.assert_called_once_with(logging.WARNING, "Qt: %s", msg)
