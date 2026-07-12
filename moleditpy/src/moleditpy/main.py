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
import traceback
from typing import Any, Optional

from .utils.constants import VERSION

# VERSION is resolved above (before Qt) so --version works without launching the app.

from PyQt6.QtCore import QtMsgType, QThread, qInstallMessageHandler
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


# Traceback signatures already shown this session, so a repeating exception
# (e.g. one thrown from a QTimer slot or a paint event) surfaces a dialog once
# instead of storming the user with an unclosable stack of modals.
_shown_exception_signatures: set[str] = set()


def _show_exception_dialog(
    exc_type: Any, exc_value: Any, exc_traceback: Any
) -> None:
    """Surface an uncaught exception in a modal dialog, at most once per bug.

    Only runs when a GUI event loop is up and we are on the GUI thread — an
    exception raised before/without the QApplication, or from a worker thread
    (where a QMessageBox is unsafe), falls back to the already-emitted log
    record. Anything that goes wrong here is swallowed to DEBUG so the
    excepthook can never recurse or mask the original error.
    """
    try:
        app = QApplication.instance()
        if app is None or QThread.currentThread() is not app.thread():
            return

        tb_text = "".join(
            traceback.format_exception(exc_type, exc_value, exc_traceback)
        )
        signature = hashlib.sha1(tb_text.encode("utf-8", "replace")).hexdigest()
        if signature in _shown_exception_signatures:
            return
        _shown_exception_signatures.add(signature)

        log_path = os.path.join(
            os.path.expanduser("~"), ".moleditpy", "moleditpy.log"
        )
        box = QMessageBox()
        box.setIcon(QMessageBox.Icon.Critical)
        box.setWindowTitle("MoleditPy — Unexpected Error")
        box.setText(
            "An unexpected error occurred. The application may be unstable; "
            "saving your work and restarting is recommended."
        )
        box.setInformativeText(
            f"{exc_type.__name__}: {exc_value}\n\n"
            f"Details are in the log:\n{log_path}"
        )
        box.setDetailedText(tb_text)
        box.setStandardButtons(QMessageBox.StandardButton.Ok)
        box.exec()
    except Exception:
        logging.debug("Failed to show exception dialog", exc_info=True)


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
        except OSError as e:
            logging.warning("Could not open log file: %s", e)

    def handle_exception(exc_type: Any, exc_value: Any, exc_traceback: Any) -> None:
        """Log unhandled exceptions and, if a GUI is up, surface them once.

        This fires only for genuinely uncaught exceptions (real bugs) — never
        for logged warnings/info, which do not pass through sys.excepthook.
        """
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        logging.error(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )
        _show_exception_dialog(exc_type, exc_value, exc_traceback)

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
