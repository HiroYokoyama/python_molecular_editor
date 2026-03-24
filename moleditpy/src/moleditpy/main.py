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

from PyQt6.QtWidgets import QApplication

try:
    from .modules.main_window import MainWindow
except ImportError:
    from modules.main_window import MainWindow


def main():
    # --- Additional handling for Windows taskbar icon ---
    if sys.platform == "win32":
        myappid = "hyoko.moleditpy.1.0"  # Application-specific ID (arbitrary)
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication(sys.argv)
    file_path = sys.argv[1] if len(sys.argv) > 1 else None
    window = MainWindow(initial_file=file_path)
    window.show()
    sys.exit(app.exec())
