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
    from .app_state import MainWindowAppState
    from .compute_logic import ComputeManager
    from .dialog_logic import DialogManager
    from .edit_3d_logic import Edit3DManager
    from .edit_actions_logic import EditActionsManager
    from .export_logic import ExportManager
    from .main_window_init import MainWindowMainInit
    from .molecular_parsers import MainWindowMolecularParsers
    from .project_io import MainWindowProjectIo
    from .string_importers import MainWindowStringImporters
    from .ui_manager import MainWindowUiManager
    from .view_3d_logic import View3DManager
    from .view_loaders import MainWindowViewLoaders
except (AttributeError, RuntimeError, TypeError):
    # Fallback to absolute imports for script-style execution
    from moleditpy_linux.ui.app_state import MainWindowAppState
    from moleditpy_linux.ui.edit_3d_logic import Edit3DManager
    from moleditpy_linux.ui.edit_actions_logic import EditActionsManager
    from moleditpy_linux.ui.export_logic import ExportManager
    from moleditpy_linux.ui.main_window_init import MainWindowMainInit
    from moleditpy_linux.ui.molecular_parsers import MainWindowMolecularParsers
    from moleditpy_linux.ui.project_io import MainWindowProjectIo
    from moleditpy_linux.ui.string_importers import MainWindowStringImporters
    from moleditpy_linux.ui.ui_manager import MainWindowUiManager
    from moleditpy_linux.ui.view_3d_logic import View3DManager
    from moleditpy_linux.ui.view_loaders import MainWindowViewLoaders


class MainWindow(
    MainWindowAppState,
    MainWindowMainInit,
    MainWindowMolecularParsers,
    MainWindowProjectIo,
    MainWindowStringImporters,
    MainWindowUiManager,
    MainWindowViewLoaders,
    QMainWindow,
):
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

        self._delegate_manager_methods(self.export_manager)
        self._delegate_manager_methods(self.view_3d_manager)
        self._delegate_manager_methods(self.edit_3d_manager)
        self._delegate_manager_methods(self.edit_actions_manager)
        self._delegate_manager_methods(self.compute_manager)
        self._delegate_manager_methods(self.dialog_manager)

        # Initialize features via Mixins
        # MainWindowMainInit handles the bulk of the UI and data setup
        MainWindowMainInit.__init__(self, initial_file, safe_mode=safe_mode)
        # MainWindowAppState handles undo/redo stack and app-wide state tracking
        MainWindowAppState.__init__(self)

        # Backwards compatibility for legacy plugins using delegation attributes
        delegations = {
            "main_window_app_state": self,
            "main_window_compute": self,
            "main_window_dialog_manager": self,
            "main_window_edit_3d": self.edit_3d_manager,
            "main_window_edit_actions": self.edit_actions_manager,
            "main_window_export": self.export_manager,
            "main_window_main_init": self,
            "main_window_molecular_parsers": self,
            "main_window_project_io": self,
            "main_window_string_importers": self,
            "main_window_ui_manager": self,
            "main_window_view_3d": self.view_3d_manager,
            "main_window_view_loaders": self,
        }
        for attr, target in delegations.items():
            setattr(self, attr, target)

    def _delegate_manager_methods(self, manager):
        """Bind public methods from manager to this instance for backwards compatibility."""
        for name in dir(manager):
            # Also delegate _compute_h_counts explicitly for internal use by other mixins
            if name == "_compute_h_counts":
                setattr(self, name, getattr(manager, name))
            elif not name.startswith("_") and callable(getattr(manager, name)):
                if not hasattr(self, name):
                    setattr(self, name, getattr(manager, name))

    # Proxy properties for View3DManager state synchronization
    @property
    def current_3d_style(self):
        return self.view_3d_manager.current_3d_style

    @current_3d_style.setter
    def current_3d_style(self, value):
        self.view_3d_manager.current_3d_style = value

    @property
    def atom_info_display_mode(self):
        return self.view_3d_manager.atom_info_display_mode

    @atom_info_display_mode.setter
    def atom_info_display_mode(self, value):
        self.view_3d_manager.atom_info_display_mode = value

    @property
    def atom_positions_3d(self):
        return self.view_3d_manager.atom_positions_3d

    @atom_positions_3d.setter
    def atom_positions_3d(self, value):
        self.view_3d_manager.atom_positions_3d = value

    @property
    def show_chiral_labels(self):
        return self.view_3d_manager.show_chiral_labels

    @show_chiral_labels.setter
    def show_chiral_labels(self, value):
        self.view_3d_manager.show_chiral_labels = value

    @property
    def is_3d_edit_mode(self):
        return self.edit_3d_manager.is_3d_edit_mode

    @is_3d_edit_mode.setter
    def is_3d_edit_mode(self, value):
        self.edit_3d_manager.is_3d_edit_mode = value

    def set_atom_from_periodic_table(self, symbol):
        self.set_mode(f"atom_{symbol}")

    def __getattr__(self, name):
        """Dynamic delegation to managers for backwards compatibility."""
        # Prevent recursion for manager attributes themselves
        if name in [
            "export_manager",
            "view_3d_manager",
            "edit_3d_manager",
            "edit_actions_manager",
        ]:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

        # Core attributes that MUST be on host (skip list)
        core_attrs = {
            "settings",
            "plotter",
            "statusBar",
            "current_mol",
            "scene",
            "view_2d",
            "data",
            "plugin_manager",
            "current_file_path",
        }
        if name in core_attrs or name.startswith("_"):
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

        # Try to find the attribute in any of the injected managers
        for mgr_attr in [
            "export_manager",
            "view_3d_manager",
            "edit_3d_manager",
            "edit_actions_manager",
        ]:
            manager = getattr(self, mgr_attr, None)
            if manager:
                # Avoid hasattr() as it triggers mgr.__getattr__
                # Check if name is in manager's class hierarchy or instance dict
                mgr_type = type(manager)
                if (
                    any(name in c.__dict__ for c in mgr_type.mro())
                    or name in manager.__dict__
                ):
                    return getattr(manager, name)

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )
