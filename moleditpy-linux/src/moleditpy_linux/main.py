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

from PyQt6.QtWidgets import QApplication

try:
    from .ui.main_window import MainWindow
except ImportError:
    # Add the parent directory (src) to sys.path so 'moleditpy.*' imports work
    src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    from moleditpy.ui.main_window import MainWindow


def setup_logging():
    """Configure global logging to standard output."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stdout,
    )


def main():
    # Setup logging as early as possible
    setup_logging()

    # --- Additional handling for Windows taskbar icon ---
    if sys.platform == "win32":
        myappid = "hyoko.moleditpy.1.0"  # Application-specific ID (arbitrary)
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    parser = argparse.ArgumentParser(
        prog="moleditpy", description="MoleditPy molecular editor"
    )
    parser.add_argument("file", nargs="?", default=None, help="File to open on startup")
    parser.add_argument(
        "--safe",
        action="store_true",
        default=False,
        help="Start in safe mode: skip loading all plugins",
    )
    # parse_known_args so Qt's own argv flags (e.g. -platform) are passed through
    args, remaining = parser.parse_known_args()

    app = QApplication([sys.argv[0]] + remaining)
    window = MainWindow(initial_file=args.file, safe_mode=args.safe)
    window.show()
    sys.exit(app.exec())
