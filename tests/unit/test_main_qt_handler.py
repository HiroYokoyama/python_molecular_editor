"""Unit tests for the Qt message handler and startup logging helpers in moleditpy.main."""

import json
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

from moleditpy.main import (
    _DOWNGRADED_QT_PATTERNS,
    _qt_message_handler,
    _read_startup_log_settings,
)


@pytest.mark.parametrize("pattern", _DOWNGRADED_QT_PATTERNS)
def test_downgraded_pattern_routes_to_debug(pattern):
    """Messages matching a downgraded pattern are logged at DEBUG, not WARNING."""
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
    """Each Qt message type is routed to the corresponding Python logging level."""
    msg = "some Qt message"
    with patch("logging.log") as mock_log, patch("logging.debug") as mock_debug:
        _qt_message_handler(mode, MagicMock(), msg)
    mock_log.assert_called_once_with(expected_level, "Qt: %s", msg)
    mock_debug.assert_not_called()


def test_unknown_mode_falls_back_to_warning():
    """Unknown Qt message type defaults to logging.WARNING."""
    msg = "unknown mode message"
    with patch("logging.log") as mock_log:
        _qt_message_handler(999, MagicMock(), msg)
    mock_log.assert_called_once_with(logging.WARNING, "Qt: %s", msg)


# ---------------------------------------------------------------------------
# _read_startup_log_settings
# ---------------------------------------------------------------------------


def test_read_startup_log_settings_both_true(tmp_path):
    """Returns (True, True) when both keys are set in settings.json."""
    settings_dir = tmp_path / ".moleditpy"
    settings_dir.mkdir()
    (settings_dir / "settings.json").write_text(
        json.dumps({"log_to_file": True, "log_level_debug": True}), encoding="utf-8"
    )
    with patch("os.path.expanduser", return_value=str(tmp_path)):
        result = _read_startup_log_settings()
    assert result == (True, True)


def test_read_startup_log_settings_defaults_false(tmp_path):
    """Returns (False, False) when keys are absent from settings.json."""
    settings_dir = tmp_path / ".moleditpy"
    settings_dir.mkdir()
    (settings_dir / "settings.json").write_text(json.dumps({}), encoding="utf-8")
    with patch("os.path.expanduser", return_value=str(tmp_path)):
        result = _read_startup_log_settings()
    assert result == (False, False)


def test_read_startup_log_settings_missing_file(tmp_path):
    """Returns (False, False) when settings.json does not exist."""
    with patch("os.path.expanduser", return_value=str(tmp_path)):
        result = _read_startup_log_settings()
    assert result == (False, False)


def test_read_startup_log_settings_invalid_json(tmp_path):
    """Returns (False, False) on malformed JSON."""
    settings_dir = tmp_path / ".moleditpy"
    settings_dir.mkdir()
    (settings_dir / "settings.json").write_text("not json{{", encoding="utf-8")
    with patch("os.path.expanduser", return_value=str(tmp_path)):
        result = _read_startup_log_settings()
    assert result == (False, False)
