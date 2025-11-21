#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""


import sys
import ctypes

from PyQt6.QtWidgets import QApplication

try:
    # When executed as part of the package (python -m moleditpy)
    from .modules.main_window import MainWindow
except Exception:
    # When executed as a standalone script (python main.py) the package-relative
    # import won't work; fall back to absolute import that works with sys.path
    from modules.main_window import MainWindow

def main():
    # --- Windows タスクバーアイコンのための追加処理 ---
    if sys.platform == 'win32':
        myappid = 'hyoko.moleditpy.1.0' # アプリケーション固有のID（任意）
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication(sys.argv)
    file_path = sys.argv[1] if len(sys.argv) > 1 else None
    window = MainWindow(initial_file=file_path)
    window.show()
    sys.exit(app.exec())
