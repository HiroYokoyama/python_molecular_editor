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
import sys
import argparse
import logging
import os
from typing import Any

from .utils.constants import VERSION

# VERSION is resolved above (before Qt) so --version works without launching the app.

from PyQt6.QtCore import QtMsgType, qInstallMessageHandler
from PyQt6.QtWidgets import QApplication

from .ui.main_window import MainWindow

_QT_LOG_LEVEL = {
    QtMsgType.QtDebugMsg: logging.DEBUG,
    QtMsgType.QtInfoMsg: logging.INFO,
    QtMsgType.QtWarningMsg: logging.WARNING,
    QtMsgType.QtCriticalMsg: logging.ERROR,
    QtMsgType.QtFatalMsg: logging.CRITICAL,
}

_DOWNGRADED_QT_PATTERNS = ("Retrying to obtain clipboard",)


def _qt_message_handler(mode: QtMsgType, _context: Any, message: str) -> None:
    """Route Qt log messages to Python logging, downgrading known noisy warnings."""
    for pattern in _DOWNGRADED_QT_PATTERNS:
        if pattern in message:
            logging.debug("Qt: %s", message)
            return
    logging.log(_QT_LOG_LEVEL.get(mode, logging.WARNING), "Qt: %s", message)


def setup_logging() -> None:
    """Configure root logger and install a global unhandled-exception handler."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s (%(pathname)s:%(lineno)d): %(message)s",
        stream=sys.stdout,
        force=True,
    )

    def handle_exception(exc_type: Any, exc_value: Any, exc_traceback: Any) -> None:
        """Log unhandled exceptions using the configured logging system."""
        if issubclass(exc_type, KeyboardInterrupt):
            # Allow keyboard interrupt to exit normally
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
        myappid = "hyoko.moleditpy_linux.1.0"
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
            from moleditpy_linux.plugins.plugin_manager import PluginManager
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
            pass

    sys.exit(app.exec())
