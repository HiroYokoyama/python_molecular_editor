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

import moleditpy.main as main_mod
from moleditpy.main import (
    _DOWNGRADED_QT_PATTERNS,
    _ErrorDialogHandler,
    _qt_message_handler,
    _read_startup_log_settings,
)


def _record(message="boom", level=logging.ERROR, exc_info=None, **extra):
    """Build a logging.LogRecord as the handler would receive it."""
    rec = logging.LogRecord(
        name="test", level=level, pathname=__file__, lineno=1,
        msg=message, args=(), exc_info=exc_info,
    )
    for key, value in extra.items():
        setattr(rec, key, value)
    return rec


# ---------------------------------------------------------------------------
# _ErrorDialogHandler.emit — guards, dedup, deferral (state is per-instance)
# ---------------------------------------------------------------------------


def test_handler_skips_without_qapplication():
    """No GUI event loop → nothing is queued (the log record is the surface)."""
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QTimer"
    ) as qtimer:
        qapp.instance.return_value = None
        _ErrorDialogHandler().emit(_record())
    qtimer.singleShot.assert_not_called()


def test_handler_skips_off_gui_thread():
    """A worker-thread record must not queue a dialog (QMessageBox unsafe)."""
    app = MagicMock()
    app.thread.return_value = "gui"
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QThread"
    ) as qthread, patch.object(main_mod, "QTimer") as qtimer:
        qapp.instance.return_value = app
        qthread.currentThread.return_value = "worker"
        _ErrorDialogHandler().emit(_record())
    qtimer.singleShot.assert_not_called()


def test_handler_honors_no_dialog_extra():
    """extra={'no_dialog': True} opts a record out entirely."""
    with patch.object(main_mod, "QTimer") as qtimer:
        _ErrorDialogHandler().emit(_record(no_dialog=True))
    qtimer.singleShot.assert_not_called()


def test_handler_queues_once_then_dedups():
    """Identical message+traceback queues one dialog; the repeat is suppressed."""
    app = MagicMock()
    app.thread.return_value = "gui"
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QThread"
    ) as qthread, patch.object(main_mod, "QTimer") as qtimer:
        qapp.instance.return_value = app
        qthread.currentThread.return_value = "gui"
        handler = _ErrorDialogHandler()
        handler.emit(_record("same error"))
        handler.emit(_record("same error"))
    assert qtimer.singleShot.call_count == 1


def test_handler_emit_never_raises():
    """emit must swallow internal failures (logging contract)."""
    app = MagicMock()
    app.thread.return_value = "gui"
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QThread"
    ) as qthread, patch.object(main_mod, "QTimer") as qtimer:
        qapp.instance.return_value = app
        qthread.currentThread.return_value = "gui"
        qtimer.singleShot.side_effect = RuntimeError("boom")
        handler = _ErrorDialogHandler()
        with patch.object(handler, "handleError") as handle_err:
            handler.emit(_record())  # must not raise
        handle_err.assert_called_once()


# ---------------------------------------------------------------------------
# _ErrorDialogHandler._show — duplicate suppression + safety
# ---------------------------------------------------------------------------


def test_show_yields_to_existing_modal():
    """If a site/plugin QMessageBox.critical is already up, stay quiet."""
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QMessageBox"
    ) as qmb:
        qapp.activeModalWidget.return_value = MagicMock()  # a modal is showing
        _ErrorDialogHandler()._show("boom", "")
    qmb.assert_not_called()


def test_show_displays_non_blocking_when_no_modal():
    """With no other modal up, the dialog is shown non-blocking (show, not exec)
    and held alive so it is not garbage-collected away."""
    handler = _ErrorDialogHandler()
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QMessageBox"
    ) as qmb:
        qapp.activeModalWidget.return_value = None
        handler._show("boom", "Traceback (most recent call last): ...")
    qmb.return_value.setDetailedText.assert_called_once()
    qmb.return_value.show.assert_called_once()
    qmb.return_value.exec.assert_not_called()  # never blocks the event loop
    assert qmb.return_value in handler._open_boxes


def test_show_swallows_its_own_errors():
    """A failure while building the dialog must not propagate."""
    handler = _ErrorDialogHandler()
    with patch.object(main_mod, "QApplication") as qapp, patch.object(
        main_mod, "QMessageBox"
    ) as qmb:
        qapp.activeModalWidget.return_value = None
        qmb.side_effect = RuntimeError("Qt is gone")
        handler._show("boom", "")  # must not raise


def test_details_names_log_file_when_file_logging_on():
    """The informative text points at the log file only when one is active."""
    with_file = _ErrorDialogHandler(log_path="/home/u/.moleditpy/moleditpy.log")
    assert "/home/u/.moleditpy/moleditpy.log" in with_file._details()

    without_file = _ErrorDialogHandler()
    hint = without_file._details()
    assert "terminal" in hint.lower()
    assert "moleditpy.log" not in hint


# ---------------------------------------------------------------------------
# setup_logging — the dialog handler is a GUI feature, skipped when headless
# ---------------------------------------------------------------------------


def _dialog_handlers():
    return [
        h for h in logging.getLogger().handlers
        if isinstance(h, _ErrorDialogHandler)
    ]


def test_dialog_handler_not_attached_when_headless(monkeypatch, tmp_path):
    """MOLEDITPY_HEADLESS must suppress the modal dialog handler entirely,
    otherwise an offscreen test run hangs on box.exec()."""
    monkeypatch.setenv("MOLEDITPY_HEADLESS", "1")
    monkeypatch.setattr("os.path.expanduser", lambda _p: str(tmp_path))
    existing = _dialog_handlers()
    for handler in existing:
        logging.getLogger().removeHandler(handler)
    try:
        main_mod.setup_logging()
        assert _dialog_handlers() == []
    finally:
        for handler in _dialog_handlers():
            logging.getLogger().removeHandler(handler)


def test_dialog_handler_attached_when_not_headless(monkeypatch, tmp_path):
    """A real (non-headless) run installs exactly one dialog handler."""
    monkeypatch.delenv("MOLEDITPY_HEADLESS", raising=False)
    monkeypatch.setattr("os.path.expanduser", lambda _p: str(tmp_path))
    for handler in _dialog_handlers():
        logging.getLogger().removeHandler(handler)
    try:
        main_mod.setup_logging()
        assert len(_dialog_handlers()) == 1
    finally:
        for handler in _dialog_handlers():
            logging.getLogger().removeHandler(handler)


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
