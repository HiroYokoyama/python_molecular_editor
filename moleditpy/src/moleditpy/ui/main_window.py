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
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QMainWindow

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .app_state import StateManager
    from .compute_logic import ComputeManager
    from .dialog_logic import DialogManager
    from .edit_3d_logic import Edit3DManager
    from .edit_actions_logic import EditActionsManager
    from .io_logic import IOManager
    from .export_logic import ExportManager
    from .main_window_init import MainInitManager
    from .string_importers import StringImporterManager
    from .ui_manager import UIManager
    from .view_3d_logic import View3DManager
except (AttributeError, RuntimeError, TypeError):
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.app_state import StateManager
    from moleditpy.ui.edit_3d_logic import Edit3DManager
    from moleditpy.ui.edit_actions_logic import EditActionsManager
    from moleditpy.ui.io_logic import IOManager
    from moleditpy.ui.export_logic import ExportManager
    from moleditpy.ui.main_window_init import MainInitManager
    from moleditpy.ui.string_importers import StringImporterManager
    from moleditpy.ui.ui_manager import UIManager
    from moleditpy.ui.view_3d_logic import View3DManager


class MainWindow(QMainWindow):
    # start_calculation carries the MOL block and an options object (second arg)
    start_calculation = pyqtSignal(str, object)

    def __init__(self, initial_file=None, safe_mode=False):
        QMainWindow.__init__(self)

        # Initialize properties
        self._is_restoring_state = False

        # Initialize Managers (Composition)
        self.export_manager = ExportManager(self)
        self.view_3d_manager = View3DManager(self)
        self.edit_3d_manager = Edit3DManager(self)
        self.edit_actions_manager = EditActionsManager(self)
        self.compute_manager = ComputeManager(self)
        self.dialog_manager = DialogManager(self)
        self.io_manager = IOManager(self)
        self.state_manager = StateManager(self)
        self.init_manager = MainInitManager(self)
        self.string_importer_manager = StringImporterManager(self)
        self.ui_manager = UIManager(self)

        # Handle all logic via Event Filter (including Close events)
        self.installEventFilter(self.ui_manager)

        # Initialize features via Mixins
        # MainWindowMainInit handles the bulk of the UI and data setup
        self.init_manager.__init__(self, initial_file, safe_mode=safe_mode)
        # MainWindowAppState handles undo/redo stack and app-wide state tracking
        self.state_manager.__init__(self)


