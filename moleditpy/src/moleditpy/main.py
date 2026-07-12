#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import ctypes
import hashlib
import json
import logging
import logging.handlers
import os
import sys
import argparse
import time
import traceback
from typing import Any, Dict, Optional

from .utils.constants import VERSION

# VERSION is resolved above (before Qt) so --version works without launching the app.

from PyQt6.QtCore import QtMsgType, QThread, QTimer, qInstallMessageHandler
from PyQt6.QtWidgets import QApplication, QMessageBox

from .ui.main_window import MainWindow

_QT_LOG_LEVEL = {
    QtMsgType.QtDebugMsg: logging.DEBUG,
    QtMsgType.QtInfoMsg: logging.INFO,
    QtMsgType.QtWarningMsg: logging.WARNING,
    QtMsgType.QtCriticalMsg: logging.ERROR,
    QtMsgType.QtFatalMsg: logging.CRITICAL,
}

_DOWNGRADED_QT_PATTERNS = ("Retrying to obtain clipboard",)


def _qt_message_handler(mode: QtMsgType, _context: Any, message: Optional[str]) -> None:
    """Route Qt log messages to Python logging, downgrading known noisy warnings."""
    message = message or ""
    for pattern in _DOWNGRADED_QT_PATTERNS:
        if pattern in message:
            logging.debug("Qt: %s", message)
            return
    logging.log(_QT_LOG_LEVEL.get(mode, logging.WARNING), "Qt: %s", message)


class _ErrorDialogHandler(logging.Handler):
    """Surface ERROR/CRITICAL log records in a dialog, once per unique message.

    The uncaught-exception hook routes through here too (it logs at ERROR with
    ``exc_info``), so every genuine failure — logged explicitly or bubbling out
    of Qt — gets one dialog. Warnings/info never reach it (handler level is
    ERROR). All state is per-instance; the one instance is created and attached
    in :func:`setup_logging` (only for a real GUI run, never headless).

    Guards: GUI thread + live QApplication only (a worker-thread record, where a
    QMessageBox is unsafe, stays log-only); a time-windowed dedup so a fast
    repeating error (e.g. from a QTimer slot) collapses to one dialog while a
    genuine retry seconds later surfaces again; and ``extra={"no_dialog": True}``
    to opt a record out. ``emit`` never raises.

    The dialog is shown non-blocking (``show``, not ``exec``): it never spins a
    nested event loop, so it cannot freeze the app — or a GUI-enabled test run —
    while it is on screen.
    """

    # Identical errors within this many seconds share one dialog; a repeat
    # after it surfaces again.
    _DEDUP_WINDOW_S = 10.0

    def __init__(self, log_path: Optional[str] = None) -> None:
        super().__init__(level=logging.ERROR)
        # signature -> monotonic timestamp of its last shown dialog.
        self._last_shown: Dict[str, float] = {}
        # Held so a non-blocking box is not garbage-collected before it closes.
        self._open_boxes: set = set()
        self._log_path = log_path

    def emit(self, record: logging.LogRecord) -> None:
        try:
            if getattr(record, "no_dialog", False):
                return
            app = QApplication.instance()
            if app is None or QThread.currentThread() is not app.thread():
                return
            if app.property("moleditpy_shutting_down"):
                return

            message = record.getMessage()
            # Source location (as the terminal log shows), then any traceback.
            detail = f"{record.pathname}:{record.lineno}"
            if record.exc_info:
                detail += "\n\n" + "".join(traceback.format_exception(*record.exc_info))

            signature = hashlib.sha1(
                (message + detail).encode("utf-8", "replace")
            ).hexdigest()
            now = time.monotonic()
            last = self._last_shown.get(signature)
            if last is not None and now - last < self._DEDUP_WINDOW_S:
                return
            self._last_shown[signature] = now

            # Defer so a site's own QMessageBox.critical (raised right after the
            # log call) is already modal when _show checks activeModalWidget.
            QTimer.singleShot(0, lambda: self._show(message, detail))
        except Exception:  # pragma: no cover - logging must never raise
            self.handleError(record)

    def _details(self) -> str:
        """Where to find the full error text, honest about file logging."""
        if self._log_path:
            return (
                "Full details are in the terminal output and the log file:\n"
                f"{self._log_path}"
            )
        return (
            "Full details are in the terminal output. Enable "
            "Settings ▸ Other ▸ 'Save log to file' to also keep them "
            "on disk."
        )

    def _show(self, message: str, detail: str) -> None:
        """Show one non-blocking error dialog, unless another modal is up.

        Sites (and plugins) log the error and then immediately raise their own
        ``QMessageBox.critical``; by the time this deferred call fires, that
        modal is up, so we detect it and stay quiet — surfacing the generic
        dialog only for errors nobody else reports.
        """
        if QApplication.activeModalWidget() is not None:
            return
        try:
            box = QMessageBox()
            box.setIcon(QMessageBox.Icon.Critical)
            box.setWindowTitle("MoleditPy — Error")
            box.setText(message or "An unexpected error occurred.")
            box.setInformativeText(self._details())
            if detail:
                box.setDetailedText(detail)
            box.setStandardButtons(QMessageBox.StandardButton.Ok)
            box.setModal(True)
            self._open_boxes.add(box)
            box.finished.connect(lambda _result, b=box: self._open_boxes.discard(b))
            box.show()
        except Exception:
            logging.debug("Failed to show error dialog", exc_info=True)


def _read_startup_log_settings() -> tuple[bool, bool]:
    """Read log_to_file and log_level_debug from settings.json before Qt starts.

    Returns (log_to_file, log_level_debug). Falls back to (False, False) on any error.
    """
    settings_path = os.path.join(os.path.expanduser("~"), ".moleditpy", "settings.json")
    try:
        with open(settings_path, encoding="utf-8") as f:
            data = json.load(f)
        return bool(data.get("log_to_file", False)), bool(
            data.get("log_level_debug", False)
        )
    except (OSError, json.JSONDecodeError, ValueError):
        return False, False


def setup_logging() -> None:
    """Configure root logger and install a global unhandled-exception handler."""
    log_to_file, log_level_debug = _read_startup_log_settings()
    level = logging.DEBUG if log_level_debug else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(name)s (%(pathname)s:%(lineno)d): %(message)s"

    logging.basicConfig(
        level=level,
        format=fmt,
        stream=sys.stdout,
        force=True,
    )

    active_log_path: Optional[str] = None
    if log_to_file:
        log_dir = os.path.join(os.path.expanduser("~"), ".moleditpy")
        try:
            os.makedirs(log_dir, exist_ok=True)
            log_path = os.path.join(log_dir, "moleditpy.log")
            fh = logging.handlers.RotatingFileHandler(
                log_path, maxBytes=1_048_576, backupCount=3, encoding="utf-8"
            )
            fh.setLevel(level)
            fh.setFormatter(logging.Formatter(fmt))
            logging.getLogger().addHandler(fh)
            logging.info("File logging enabled: %s", log_path)
            active_log_path = log_path
        except OSError as e:
            logging.warning("Could not open log file: %s", e)

    # Surface every ERROR/CRITICAL record in a dialog (see _ErrorDialogHandler);
    # uncaught exceptions flow through here too via handle_exception's log call.
    # It is a GUI feature, so it is skipped in headless mode (MOLEDITPY_HEADLESS),
    # where a blocking modal would hang the run.
    if not os.environ.get("MOLEDITPY_HEADLESS"):
        logging.getLogger().addHandler(_ErrorDialogHandler(log_path=active_log_path))

    def handle_exception(exc_type: Any, exc_value: Any, exc_traceback: Any) -> None:
        """Log unhandled exceptions; the ERROR dialog handler surfaces them.

        The dialog fires only for genuine failures — this log call and explicit
        logging.error/critical calls — never for warnings/info.
        """
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        logging.error(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )

    sys.excepthook = handle_exception


def main() -> None:
    """Parse CLI arguments, configure logging, and launch the GUI."""
    setup_logging()

    if sys.platform == "win32":
        # Taskbar grouping ID follows the major version
        major = VERSION.split(".")[0] if VERSION and VERSION != "Unknown" else "0"
        myappid = f"hyoko.moleditpy.{major}"
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    parser = argparse.ArgumentParser(
        prog="moleditpy",
        description="MoleditPy — A Python-based molecular editing software",
    )
    parser.add_argument("file", nargs="?", default=None, help="File to open on startup")
    _variant = " (Linux)" if "moleditpy_linux" in (__file__ or "") else ""
    parser.add_argument(
        "--version", action="version", version=f"MoleditPy{_variant} {VERSION}"
    )
    parser.add_argument(
        "--safe",
        action="store_true",
        default=False,
        help="Start in safe mode: skip loading all plugins",
    )
    parser.add_argument(
        "--install-plugin",
        metavar="PATH",
        help="Install a plugin from a .py file, .zip, or folder (Headless)",
    )
    # parse_known_args so Qt's own argv flags (e.g. -platform) are passed through
    args, remaining = parser.parse_known_args()

    if args.install_plugin:
        plugin_path = os.path.abspath(args.install_plugin)
        if not os.path.exists(plugin_path):
            print(f"Error: Plugin path not found: {plugin_path}")
            sys.exit(1)

        try:
            from moleditpy.plugins.plugin_manager import PluginManager
        except ImportError:
            from .plugins.plugin_manager import PluginManager

        pm = PluginManager()
        sha256 = pm.compute_sha256(plugin_path)

        metadata_file = plugin_path
        if os.path.isdir(plugin_path):
            init_py = os.path.join(plugin_path, "__init__.py")
            if os.path.exists(init_py):
                metadata_file = init_py

        info = (
            pm.get_plugin_info_safe(metadata_file)
            if metadata_file.endswith(".py")
            else {}
        )

        print("\n" + "=" * 40)
        print("      PLUGIN INSTALLATION (HEADLESS)")
        print("=" * 40)
        print(f" Name:        {info.get('name', os.path.basename(plugin_path))}")
        print(f" Author:      {info.get('author', 'Unknown')}")
        print(f" Version:     {info.get('version', 'Unknown')}")
        print(f" Description: {info.get('description', 'No description')}")
        print("-" * 40)
        print(f" Path:        {plugin_path}")
        print(f" SHA-256:     {sha256}")
        print("=" * 40)

        confirm = (
            input("\nDo you want to proceed with installation? (y/N): ").strip().lower()
        )
        if confirm == "y":
            success, msg = pm.install_plugin(plugin_path)
            if success:
                print(f"Success: {msg}")
                sys.exit(0)
            else:
                print(f"Error: {msg}")
                sys.exit(1)
        else:
            print("Installation aborted.")
            sys.exit(0)

    app = QApplication([sys.argv[0]] + remaining)
    qInstallMessageHandler(_qt_message_handler)
    window = MainWindow(initial_file=args.file, safe_mode=args.safe)
    window.show()

    if sys.platform == "win32":
        try:
            from PyQt6.QtCore import QTimer

            QTimer.singleShot(100, lambda: window.setWindowIcon(window.windowIcon()))
        except (
            Exception
        ):  # [COSMETIC] Icon refresh is best-effort; Qt timing errors are non-fatal.
            logging.debug("Suppressed non-critical error", exc_info=True)

    sys.exit(app.exec())
