#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

# PyQt6 Modules
from PyQt6.QtCore import pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import QMainWindow

try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .main_window_app_state import MainWindowAppState
    from .main_window_compute import MainWindowCompute
    from .main_window_dialog_manager import MainWindowDialogManager
    from .main_window_edit_3d import MainWindowEdit3d
    from .main_window_edit_actions import MainWindowEditActions
    from .main_window_export import MainWindowExport
    from .main_window_main_init import MainWindowMainInit
    from .main_window_molecular_parsers import MainWindowMolecularParsers
    from .main_window_project_io import MainWindowProjectIo
    from .main_window_string_importers import MainWindowStringImporters
    from .main_window_ui_manager import MainWindowUiManager
    from .main_window_view_3d import MainWindowView3d
    from .main_window_view_loaders import MainWindowViewLoaders
except (AttributeError, RuntimeError, TypeError):
    # Fallback to absolute imports for script-style execution
    from modules.main_window_app_state import MainWindowAppState
    from modules.main_window_compute import MainWindowCompute
    from modules.main_window_dialog_manager import MainWindowDialogManager
    from modules.main_window_edit_3d import MainWindowEdit3d
    from modules.main_window_edit_actions import MainWindowEditActions
    from modules.main_window_export import MainWindowExport
    from modules.main_window_main_init import MainWindowMainInit
    from modules.main_window_molecular_parsers import MainWindowMolecularParsers
    from modules.main_window_project_io import MainWindowProjectIo
    from modules.main_window_string_importers import MainWindowStringImporters
    from modules.main_window_ui_manager import MainWindowUiManager
    from modules.main_window_view_3d import MainWindowView3d
    from modules.main_window_view_loaders import MainWindowViewLoaders


class MainWindow(
    MainWindowAppState,
    MainWindowCompute,
    MainWindowDialogManager,
    MainWindowEdit3d,
    MainWindowEditActions,
    MainWindowExport,
    MainWindowMainInit,
    MainWindowMolecularParsers,
    MainWindowProjectIo,
    MainWindowStringImporters,
    MainWindowUiManager,
    MainWindowView3d,
    MainWindowViewLoaders,
    QMainWindow,
):
    # start_calculation carries the MOL block and an options object (second arg)
    start_calculation = pyqtSignal(str, object)

    def __init__(self, initial_file=None):
        QMainWindow.__init__(self)
        
        # Initialize properties
        self._is_restoring_state = False
        
        # Initialize features via Mixins
        # MainWindowMainInit handles the bulk of the UI and data setup
        MainWindowMainInit.__init__(self, initial_file)
        # MainWindowAppState handles undo/redo stack and app-wide state tracking
        MainWindowAppState.__init__(self)

        # Backwards compatibility for legacy plugins using delegation attributes
        for attr in [
            "main_window_app_state", "main_window_compute", "main_window_dialog_manager",
            "main_window_edit_3d", "main_window_edit_actions", "main_window_export",
            "main_window_main_init", "main_window_molecular_parsers", "main_window_project_io",
            "main_window_string_importers", "main_window_ui_manager", "main_window_view_3d",
            "main_window_view_loaders"
        ]:
            setattr(self, attr, self)

    def set_atom_from_periodic_table(self, symbol):
        self.set_mode(f"atom_{symbol}")
