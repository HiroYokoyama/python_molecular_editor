#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import contextlib
import json
import math
import os
import sys

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem

# PyQt6 Modules
from PyQt6.QtCore import QLineF, QPointF, QRectF, Qt, QTimer, QUrl
from PyQt6.QtGui import (
    QAction,
    QActionGroup,
    QBrush,
    QColor,
    QDesktopServices,
    QFont,
    QIcon,
    QKeySequence,
    QPainter,
    QPen,
    QPixmap,
    QPolygonF,
)
from PyQt6.QtWidgets import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QToolBar,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

try:
    from .. import OBABEL_AVAILABLE
    from ..utils.default_settings import DEFAULT_SETTINGS
    from ..plugins.plugin_manager import PluginManager
    from ..utils.system_utils import detect_system_dark_mode
except ImportError:
    from moleditpy import OBABEL_AVAILABLE
    from moleditpy.utils.default_settings import DEFAULT_SETTINGS
    from moleditpy.plugins.plugin_manager import PluginManager
    from moleditpy.utils.system_utils import detect_system_dark_mode

try:
    import winreg
except ImportError:
    winreg = None


try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except (AttributeError, RuntimeError, TypeError):
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .color_settings_dialog import ColorSettingsDialog
    from ..utils.constants import DEFAULT_CPK_COLORS, NUM_DASHES, VERSION
    from .custom_qt_interactor import CustomQtInteractor
    from ..core.molecular_data import MolecularData
    from .molecule_scene import MoleculeScene
    from .settings_dialog import SettingsDialog
    from .zoomable_view import ZoomableView
except (AttributeError, RuntimeError, TypeError):
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.color_settings_dialog import ColorSettingsDialog
    from moleditpy.utils.constants import NUM_DASHES, VERSION
    from moleditpy.ui.custom_qt_interactor import CustomQtInteractor
    from moleditpy.core.molecular_data import MolecularData
    from moleditpy.ui.molecule_scene import MoleculeScene
    from moleditpy.ui.settings_dialog import SettingsDialog
    from moleditpy.ui.zoomable_view import ZoomableView


class MainInitManager:
    """Feature class separated from main_window.py"""

    def __init__(self, host, initial_file=None, safe_mode=False):
        self.host = host
        self.host.init_manager = self
        # This helper is not used as a mixin in this project; initialization
        # happens on the `QMainWindow` base class in
        # `MainWindow.__init__` directly.
        self.host.setAcceptDrops(True)
        self.settings_dir = os.path.join(os.path.expanduser("~"), ".moleditpy")
        self.settings_file = os.path.join(self.settings_dir, "settings.json")
        self.settings = {}  # Will be populated by load_settings
        self.load_settings()
        self.host.initial_settings = self.settings.copy()
        self.host.setWindowTitle("MoleditPy Ver. " + VERSION)
        self.host.setGeometry(100, 100, 1400, 800)
        self.host.state_manager.data = MolecularData()
        self.host.view_3d_manager.current_mol = None
        self.host.ui_manager.is_2d_editable = True
        self.host.is_xyz_derived = (
            False  # Flag indicating if the molecule is derived from XYZ
        )
        # Chemical check flags: whether a chemical/sanitization check was attempted and whether it failed
        self.host.chem_check_tried = False
        self.host.chem_check_failed = False
        self.host._template_dialog = None
        self.host._picking_consumed = False
        self.host.init_manager.mode_actions = {}

        # Variable tracking the saved state
        self.host.state_manager.has_unsaved_changes = False
        self.settings_dirty = True
        self.host.init_manager.current_file_path = None
        self.host.initialization_complete = False
        self.host._ih_update_counter = 0

        # Initialization of the plugin manager
        if safe_mode:
            print("Safe mode: plugins disabled.")
            self.host.plugin_manager = None
        else:
            try:
                self.host.plugin_manager = PluginManager()
            except (AttributeError, RuntimeError, ValueError) as e:
                print(f"Failed to initialize PluginManager: {e}")
                self.host.plugin_manager = None

        # Dictionary holding data for plugins that haven't been loaded
        self._preserved_plugin_data = {}

        self.init_ui()
        self.init_worker_thread()
        self.host.ui_manager._setup_3d_picker()
        # Install event filter to capture window close events (handled in UIManager)
        self.host.installEventFilter(self.host.ui_manager)

        # RDKit Warm-up (initial execution cost)
        try:
            # Create a molecule with a variety of common atoms to ensure
            # the valence/H-count machinery is fully initialized.
            warmup_smiles = "OC(N)C(S)P"
            warmup_mol = Chem.MolFromSmiles(warmup_smiles)
            if warmup_mol:
                for atom in warmup_mol.GetAtoms():
                    atom.GetNumImplicitHs()
        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"RDKit warm-up failed: {e}")

        self.host.state_manager.reset_undo_stack()
        self.host.init_manager.scene.selectionChanged.connect(self.host.edit_actions_manager.update_edit_menu_actions)
        QApplication.clipboard().dataChanged.connect(self.host.edit_actions_manager.update_edit_menu_actions)

        self.host.edit_actions_manager.update_edit_menu_actions()

        if initial_file:
            self.load_command_line_file(initial_file)

        QTimer.singleShot(0, self.apply_initial_settings)
        # Camera initialization flag (permits reset only during the first draw)
        self._camera_initialized = False

        # Set initial menu text and state
        self.host.view_3d_manager.update_atom_id_menu_text()
        self.host.view_3d_manager.update_atom_id_menu_state()

        # Set initialization complete
        self.host.initialization_complete = True
        self.host.state_manager.update_window_title()  # Update title after initialization is complete
        # Ensure initial keyboard/mouse focus is placed on the 2D view
        # when opening a file or starting the application. This avoids
        # accidental focus landing on toolbar/buttons (e.g. Optimize 2D).
        try:
            QTimer.singleShot(0, self.host.init_manager.view_2d.setFocus)
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

    def init_ui(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if sys.platform == "win32":
            icon_path = os.path.join(script_dir, "..", "assets", "icon.ico")
        else:
            icon_path = os.path.join(script_dir, "..", "assets", "icon.png")

        # Create a QIcon object directly from the file path
        if os.path.exists(icon_path):
            app_icon = QIcon(icon_path)

            # Set the icon for both the window and the application
            self.host.setWindowIcon(app_icon)
            QApplication.instance().setWindowIcon(app_icon)
        else:
            print(f"Warning: Icon file not found: {icon_path}")

        self._init_main_layout()
        self._init_toolbars()
        self.init_menu_bar()

        self._setup_action_groups(self.toolbar, self.toolbar_bottom)

    def init_menu_bar(self):
        menu_bar = self.host.menuBar()

        self._init_file_menu(menu_bar)
        self._init_edit_menu(menu_bar)
        self._init_view_menu(menu_bar)
        self._init_analysis_menu(menu_bar)
        self._init_edit_3d_menu(menu_bar)
        self._init_plugin_menu(menu_bar)
        self._init_settings_menu(menu_bar)
        self._init_help_menu(menu_bar)

        # Consistently set initial state for 3D-related features
        self.host.ui_manager._enable_3d_features(False)

        # Finally, populate plugins now that all menus are created
        self.update_plugin_menu(self.plugin_menu)

    def init_worker_thread(self):
        # Initialize shared state for calculation runs.
        self.halt_ids = set()
        self.next_conversion_id = 1
        self.active_worker_ids = set()
        # Track active threads for diagnostics/cleanup (weak references ok)
        try:
            self.host.compute_manager._active_calc_threads = []
        except (AttributeError, RuntimeError, ValueError, TypeError):
            self.host.compute_manager._active_calc_threads = []

    def load_command_line_file(self, file_path):
        """Open file specified by command-line argument"""
        if not file_path or not os.path.exists(file_path):
            return

        # Helper for extension
        _, ext_with_dot = os.path.splitext(file_path)
        ext_with_dot = ext_with_dot.lower()
        # Legacy variable name (no dot)
        file_ext = ext_with_dot.lstrip(".")

        # 1. Custom Plugin Openers
        if ext_with_dot in self.host.plugin_manager.file_openers:
            openers = self.host.plugin_manager.file_openers[ext_with_dot]
            # Iterate through openers (already sorted by priority)
            for opener_info in openers:
                try:
                    callback = opener_info["callback"]
                    # Try to call the opener
                    callback(file_path)

                    self.host.init_manager.current_file_path = file_path
                    self.host.state_manager.update_window_title()
                    return  # Success
                except (AttributeError, RuntimeError, ValueError) as e:
                    print(
                        f"Plugin opener failed for '{opener_info.get('plugin', 'Unknown')}': {e}"
                    )
                    # If this opener fails, try the next one or fall through to default
                    continue

        if file_ext in ["mol", "sdf"]:
            self.host.io_manager.load_mol_file_for_3d_viewing(file_path)
        elif file_ext == "xyz":
            self.host.io_manager.load_xyz_for_3d_viewing(file_path)
        elif file_ext in ["pmeraw", "pmeprj"]:
            self.host.io_manager.open_project_file(file_path=file_path)
        else:
            self.host.statusBar().showMessage(f"Unsupported file type: {file_ext}")

    def apply_initial_settings(self):
        """Apply saved settings to the 3D view after UI initialization is complete."""

        try:
            self.update_cpk_colors_from_settings()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu/settings sync errors

        if self.host.view_3d_manager.plotter and self.host.view_3d_manager.plotter.renderer:
            bg_color = self.host.init_manager.settings.get("background_color", "#919191")
            self.host.view_3d_manager.plotter.set_background(bg_color)
            self.host.view_3d_manager.apply_3d_settings()
            
            # Redraw if molecule exists
            mol = getattr(self.host.view_3d_manager, "current_mol", None)
            if mol:
                self.host.view_3d_manager.draw_molecule_3d(mol)

        try:
            if hasattr(self.host.init_manager, 'scene') and self.host.init_manager.scene:
                # Apply 2D background color
                bg_color_2d = self.host.init_manager.settings.get("background_color_2d", "#FFFFFF")
                self.host.init_manager.scene.setBackgroundBrush(QBrush(QColor(bg_color_2d)))

                for it in list(self.host.init_manager.scene.items()):
                    if hasattr(it, "update_style"):  # [SAFE]
                        it.update_style()
                self.host.init_manager.scene.update()
                for v in list(self.host.init_manager.scene.views()):
                    v.viewport().update()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu/settings sync errors

    def update_cpk_colors_from_settings(self):
        """Update global CPK_COLORS and CPK_COLORS_PV from saved settings overrides.

        This modifies the in-memory CPK_COLORS mapping (not persisted until settings are saved).
        Only keys present in self.host.init_manager.settings['cpk_colors'] are changed; other elements keep the defaults.
        """
        try:
            overrides = self.host.init_manager.settings.get("cpk_colors", {}) or {}

            modules_to_update = []
            for name, mod in sys.modules.items():
                if (name.endswith("utils.constants") or name == "constants" or "moleditpy.utils.constants" in name):
                    if hasattr(mod, "CPK_COLORS") and isinstance(mod.CPK_COLORS, dict):
                        modules_to_update.append(mod)
            
            if not modules_to_update:
                try:
                    from . import constants as constants_mod
                    modules_to_update.append(constants_mod)
                except ImportError:
                    import moleditpy.utils.constants as constants_mod
                    modules_to_update.append(constants_mod)

            for constants_mod in modules_to_update:
                constants_mod.CPK_COLORS.clear()
                for k, v in DEFAULT_CPK_COLORS.items():
                    constants_mod.CPK_COLORS[k] = (
                        QColor(v) if not isinstance(v, QColor) else v
                    )
                
                # Apply overrides from settings
                for k, hexv in overrides.items():
                    if isinstance(hexv, str) and hexv:
                        constants_mod.CPK_COLORS[k] = QColor(hexv)

                # Rebuild the PV representation in-place too
                if hasattr(constants_mod, "CPK_COLORS_PV"):
                    constants_mod.CPK_COLORS_PV.clear()
                    for k, c in constants_mod.CPK_COLORS.items():
                        constants_mod.CPK_COLORS_PV[k] = [c.redF(), c.greenF(), c.blueF()]
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error(f"REPORT ERROR: Missing attribute 'CPK_COLORS_PV' on constants_mod")
        except (AttributeError, RuntimeError, TypeError, ValueError) as e:
            print(f"Failed to update CPK colors from settings: {e}")

    def open_settings_dialog(self):
        dialog = SettingsDialog(self.host.init_manager.settings, self.host)
        # Settings application and 3D view updates are handled by the accept() method.
        dialog.exec()

    def reset_all_settings_menu(self):
        """Prompt the user and reset all application settings to defaults."""
        if not self._confirm_settings_reset():
            return

        try:
            self._perform_settings_reset()
            self._refresh_ui_after_reset()
            QMessageBox.information(
                self.host, "Reset Complete", "All settings have been reset to defaults."
            )
        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.warning(self, "Reset Failed", f"Could not reset settings: {e}")

    def _confirm_settings_reset(self):
        """Show a confirmation dialog for resetting settings."""
        return (
            QMessageBox.question(
                self.host,
                "Reset All Settings",
                "Are you sure you want to reset all settings to defaults?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            == QMessageBox.StandardButton.Yes
        )

    def _perform_settings_reset(self):
        """Delete the settings file and reload defaults."""
        if os.path.exists(self.host.init_manager.settings_file):
            os.remove(self.host.init_manager.settings_file)
        self.load_settings()
        self.host.init_manager.settings_dirty = True

    def _refresh_ui_after_reset(self):
        """Update all UI components to reflect the reset settings."""
        # 1. Refresh related dialogs if open
        for w in QApplication.topLevelWidgets():
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                if isinstance(w, ColorSettingsDialog):
                    w.refresh_ui()
                if isinstance(w, SettingsDialog):
                    w.update_ui_from_settings(self.host.init_manager.settings)

        # 2. Update internal state and sync CPK colors
        self.optimization_method = self.host.init_manager.settings.get(
            "optimization_method", "MMFF_RDKIT"
        )
        try:
            self.update_cpk_colors_from_settings()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(f"Suppressed exception: {e}")

        # 3. Synchronize menu actions with new settings
        self._sync_settings_to_menu_actions()

        # 4. Refresh 2D and 3D views
        self._refresh_views_after_reset()

    def _sync_settings_to_menu_actions(self):
        """Synchronize menu action checked states with current settings."""
        with contextlib.suppress(AttributeError, RuntimeError, TypeError):
            # Optimization actions
            if hasattr(self, "opt3d_actions"):
                current_method = (getattr(self, "optimization_method", "") or "").upper()
                for key, action in self.opt3d_actions.items():
                    # Use case-insensitive comparison for robustness
                    action.setChecked(key.upper() == current_method)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(f"REPORT ERROR: Missing attribute 'opt3d_actions' on self")

            # Conversion actions
            if hasattr(self, "conv_actions"):
                mode = (
                    self.host.init_manager.settings.get("3d_conversion_mode", "fallback") or ""
                ).lower()
                for key, action in self.conv_actions.items():
                    action.setChecked(key.lower() == mode)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(f"REPORT ERROR: Missing attribute 'conv_actions' on self")

            # Intermolecular interaction
            if hasattr(self, "intermolecular_rdkit_action"):
                self.intermolecular_rdkit_action.setChecked(
                    self.host.init_manager.settings.get("optimize_intermolecular_interaction_rdkit", True)
                )
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(f"REPORT ERROR: Missing attribute 'intermolecular_rdkit_action' on self")

    def _refresh_views_after_reset(self):
        """Refresh 2D and 3D views after settings reset."""
        # Refresh 3D View
        try:
            self.host.view_3d_manager.apply_3d_settings()
            if hasattr(self.host, "current_mol") and self.host.view_3d_manager.current_mol:
                self.host.view_3d_manager.draw_molecule_3d(self.host.view_3d_manager.current_mol)
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(f"Suppressed exception: {e}")

        # Refresh 2D View
        if hasattr(self.host.init_manager, 'scene') and self.host.init_manager.scene:
            try:
                bg_c = self.host.init_manager.settings.get("background_color_2d", "#FFFFFF")
                self.host.init_manager.scene.setBackgroundBrush(QBrush(QColor(bg_c)))
                for item in self.host.init_manager.scene.items():
                    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                        if hasattr(item, "update_style"):
                            item.update_style()
                        else:  # [REPORT ERROR MISSING ATTRIBUTE]
                            logging.error(f"REPORT ERROR: Missing attribute 'update_style' on item")
                self.host.init_manager.scene.update()
                for v in self.host.init_manager.scene.views():
                    v.viewport().update()
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                logging.debug(f"Suppressed exception: {e}")

    def load_settings(self):
        """Load settings from a JSON file, or use defaults if the file is missing."""
        # 1. Start with default settings
        self.host.init_manager.settings = self._get_default_settings()

        # 2. Try to load from user's settings file
        try:
            if os.path.exists(self.host.init_manager.settings_file):
                with open(self.host.init_manager.settings_file, "r", encoding="utf-8") as f:
                    loaded_settings = json.load(f)

                    # 3. Handle legacy settings migration
                    self._migrate_legacy_settings(loaded_settings)

                    # 4. Update settings with loaded values
                    self.host.init_manager.settings.update(loaded_settings)
        except (AttributeError, RuntimeError, ValueError, TypeError, IOError):
            # Use defaults on any error
            pass

        # 4.5 Save initial settings copy for change detection
        self.host.initial_settings = self.host.init_manager.settings.copy()

        # 5. Apply loaded settings to application state
        self.host.view_3d_manager.show_chiral_labels = self.host.init_manager.settings.get("show_chiral_labels", False)
        # Apply optimization method
        if "optimization_method" in self.host.init_manager.settings:
            self.host.init_manager.optimization_method = self.host.init_manager.settings["optimization_method"]

    def save_settings(self):
        try:
            if not os.path.exists(self.host.init_manager.settings_dir):
                os.makedirs(self.host.init_manager.settings_dir)
            with open(self.host.init_manager.settings_file, "w", encoding="utf-8") as f:
                json.dump(self.host.init_manager.settings, f, indent=4)
            self.host.init_manager.settings_dirty = False
        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error saving settings: {e}")

    def update_plugin_menu(self, plugin_menu):
        """Discovers plugins and updates the plugin menu actions."""
        if not self.host.plugin_manager:
            return

        # Cleanup
        self._clear_all_plugin_actions(plugin_menu)

        # Re-add Manager
        manage_plugins_action = QAction("Plugin Manager...", self.host)
        manage_plugins_action.triggered.connect(
            lambda: self._show_plugin_manager(plugin_menu)
        )
        plugin_menu.addAction(manage_plugins_action)
        plugin_menu.addSeparator()

        # Discover
        plugins = self.host.plugin_manager.discover_plugins(self.host)

        # Integrate
        self._update_style_menu_with_plugins()
        self._add_registered_plugin_actions()
        self._add_plugin_toolbar_actions()
        self._add_legacy_plugin_actions(plugin_menu, plugins)
        self._integrate_plugin_export_actions()
        self._integrate_plugin_file_openers()
        self._integrate_plugin_analysis_tools()

    def _show_plugin_manager(self, plugin_menu):
        """Displays the plugin manager window and refreshes the menu."""
        if not self.host.plugin_manager:
            QMessageBox.information(
                self.host, "Safe Mode", "Plugins are disabled (safe mode)."
            )
            return
        from ..plugins.plugin_manager_window import PluginManagerWindow

        dlg = PluginManagerWindow(self.host.plugin_manager, self.host)
        dlg.exec()
        self.update_plugin_menu(plugin_menu)

    def _clear_all_plugin_actions(self, plugin_menu):
        """Clear all plugin actions from the plugin menu and other UI components."""
        PLUGIN_ACTION_TAG = "plugin_managed"

        def clear_menu(menu):
            if not menu:
                return
            for act in list(menu.actions()):
                if act.data() == PLUGIN_ACTION_TAG:
                    menu.removeAction(act)
                elif act.menu():
                    clear_menu(act.menu())

        plugin_menu.clear()
        for top_action in self.host.menuBar().actions():
            if top_action.menu():
                clear_menu(top_action.menu())

        if hasattr(self.host.init_manager, 'export_button') and self.host.init_manager.export_button.menu():
            clear_menu(self.host.init_manager.export_button.menu())

    def _update_style_menu_with_plugins(self):
        """Update the 3D style menu with custom styles from plugins."""
        if not hasattr(self.host, "style_button") or not self.host.init_manager.style_button.menu():
            return

        style_menu = self.host.init_manager.style_button.menu()
        style_group = next(
            (a.actionGroup() for a in style_menu.actions() if a.actionGroup()), None
        )

        if style_group and self.host.plugin_manager.custom_3d_styles:
            if not style_menu.actions()[-1].isSeparator():
                style_menu.addSeparator()

            for style_name in self.host.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self.host, checkable=True)
                    action.triggered.connect(
                        lambda checked=False, s=style_name: self.host.view_3d_manager.set_3d_style(s)
                    )
                    style_menu.addAction(action)
                    style_group.addAction(action)

    def _add_registered_plugin_actions(self):
        """Add actions that have been explicitly registered via the plugin manager."""
        PLUGIN_ACTION_TAG = "plugin_managed"
        if not self.host.plugin_manager.menu_actions:
            return

        for action_def in self.host.plugin_manager.menu_actions:
            path = action_def["path"]
            callback = action_def["callback"]
            text = action_def["text"]

            parts = path.split("/")
            top_level_title = parts[0]
            current_menu = next(
                (
                    a.menu()
                    for a in self.host.menuBar().actions()
                    if a.menu() and a.text().replace("&", "") == top_level_title
                ),
                None,
            )

            if not current_menu:
                current_menu = self.host.menuBar().addMenu(top_level_title)

            for part in parts[1:-1]:
                sub = next(
                    (
                        a.menu()
                        for a in current_menu.actions()
                        if a.menu() and a.text().replace("&", "") == part
                    ),
                    None,
                )
                current_menu = sub if sub else current_menu.addMenu(part)

            actions = current_menu.actions()
            if (
                actions
                and not actions[-1].isSeparator()
                and actions[-1].data() != PLUGIN_ACTION_TAG
            ):
                sep = current_menu.addSeparator()
                sep.setData(PLUGIN_ACTION_TAG)

            action = QAction(text or parts[-1], self.host)
            action.triggered.connect(callback)
            if action_def.get("shortcut"):
                action.setShortcut(QKeySequence(action_def["shortcut"]))
            action.setData(PLUGIN_ACTION_TAG)
            current_menu.addAction(action)

    def _add_plugin_toolbar_actions(self):
        """Add toolbar actions registered by plugins."""
        if not hasattr(self, "plugin_toolbar"):
            return

        self.plugin_toolbar.clear()
        if self.host.plugin_manager.toolbar_actions:
            self.plugin_toolbar.show()
            for action_def in self.host.plugin_manager.toolbar_actions:
                action = QAction(action_def["text"], self.host)
                action.triggered.connect(action_def["callback"])
                if action_def["icon"] and os.path.exists(action_def["icon"]):
                    action.setIcon(QIcon(action_def["icon"]))
                if action_def["tooltip"]:
                    action.setToolTip(action_def["tooltip"])
                self.plugin_toolbar.addAction(action)
        else:
            self.plugin_toolbar.hide()

    def _add_legacy_plugin_actions(self, plugin_menu, plugins):
        """Add folder-based legacy plugin actions to the plugin menu."""
        if not plugins:
            no_plugin = QAction("(No plugins found)", self.host)
            no_plugin.setEnabled(False)
            plugin_menu.addAction(no_plugin)
            return

        # Categorize
        categorized = {}
        root = []
        for p in plugins:
            if hasattr(p["module"], "run"):
                cat = p.get("category", p.get("rel_folder", "")).strip()
                if cat:
                    categorized.setdefault(cat, []).append(p)
                else:
                    root.append(p)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(f"REPORT ERROR: Missing attribute 'run' on object")

        # Build categorized menus
        for cat in sorted(categorized.keys()):
            current_parent = plugin_menu
            for part in cat.split(os.sep):
                sub = next(
                    (
                        a.menu()
                        for a in current_parent.actions()
                        if a.menu() and a.text().replace("&", "") == part
                    ),
                    None,
                )
                current_parent = sub if sub else current_parent.addMenu(part)

            for p in sorted(categorized[cat], key=lambda x: x["name"]):
                a = QAction(p["name"], self.host)
                a.triggered.connect(
                    lambda checked, mod=p["module"]: self.host.plugin_manager.run_plugin(
                        mod, self
                    )
                )
                current_parent.addAction(a)

        # Add root items
        for p in sorted(root, key=lambda x: x["name"]):
            a = QAction(p["name"], self.host)
            a.triggered.connect(
                lambda checked, mod=p["module"]: self.host.plugin_manager.run_plugin(
                    mod, self
                )
            )
            plugin_menu.addAction(a)

    def _integrate_plugin_export_actions(self):
        """Add plugin-provided export actions to the export menu."""
        if not self.host.plugin_manager.export_actions:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        main_export_menu = None
        for top_action in self.host.menuBar().actions():
            if top_action.menu() and top_action.text().replace("&", "") == "File":
                main_export_menu = next(
                    (
                        a.menu()
                        for a in top_action.menu().actions()
                        if a.menu() and a.text().replace("&", "") == "Export"
                    ),
                    None,
                )
                if main_export_menu:
                    break

        targets = []
        if hasattr(self.host.init_manager, 'export_button') and self.host.init_manager.export_button.menu():
            targets.append(self.host.init_manager.export_button.menu())
        if main_export_menu:
            targets.append(main_export_menu)

        for menu in targets:
            sep = menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for exp in self.host.plugin_manager.export_actions:
                a = QAction(exp["label"], self.host)
                a.triggered.connect(exp["callback"])
                a.setData(PLUGIN_ACTION_TAG)
                menu.addAction(a)

    def _integrate_plugin_file_openers(self):
        """Add plugin-provided file openers to the import menu."""
        if not hasattr(self.host, "import_menu") or not self.host.plugin_manager.file_openers:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        sep = self.import_menu.addSeparator()
        sep.setData(PLUGIN_ACTION_TAG)

        plugin_map = {}
        for ext, openers in self.host.plugin_manager.file_openers.items():
            for info in openers:
                plugin_map.setdefault(info.get("plugin", "Plugin"), {})[ext] = info[
                    "callback"
                ]

        for p_name, ext_map in sorted(plugin_map.items()):
            exts = sorted(ext_map.keys())
            ext_str = "/".join(exts)
            if len(ext_str) > 30:
                cutoff = ext_str.rfind("/", 0, 30)
                ext_str = (ext_str[:cutoff] if cutoff != -1 else ext_str[:30]) + "/..."

            filter_str = (
                f"{p_name} Files ({' '.join(['*' + e for e in exts])});;All Files (*)"
            )

            def make_cb(m, f, n):
                def _cb():
                    fpath, _ = QFileDialog.getOpenFileName(
                        self.host, f"Import {n} Files", "", f
                    )
                    if fpath:
                        ext = os.path.splitext(fpath)[1].lower()
                        if ext in m:
                            m[ext](fpath)
                            self.host.init_manager.current_file_path = fpath
                            self.host.state_manager.update_window_title()

                return _cb

            a = QAction(f"Import {ext_str} ({p_name})...", self.host)
            a.triggered.connect(make_cb(ext_map, filter_str, p_name))
            a.setData(PLUGIN_ACTION_TAG)
            self.import_menu.addAction(a)

    def _integrate_plugin_analysis_tools(self):
        """Add plugin-provided analysis tools to the analysis menu."""
        analysis_menu = next(
            (
                a.menu()
                for a in self.host.menuBar().actions()
                if a.text().replace("&", "") == "Analysis"
            ),
            None,
        )
        if analysis_menu and self.host.plugin_manager.analysis_tools:
            PLUGIN_ACTION_TAG = "plugin_managed"
            sep = analysis_menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for tool in self.host.plugin_manager.analysis_tools:
                a = QAction(f"{tool['label']} ({tool.get('plugin', 'Plugin')})", self.host)
                a.triggered.connect(tool["callback"])
                a.setData(PLUGIN_ACTION_TAG)
                analysis_menu.addAction(a)

    # --- UI Initialization Helpers ---
    def _init_main_layout(self):
        """Initialize the main layout with splitter and panels."""
        self.host.init_manager.splitter = QSplitter(Qt.Orientation.Horizontal)
        # Make splitter handle thicker for visibility
        self.host.init_manager.splitter.setHandleWidth(8)
        # Improve splitter handle style
        self.host.init_manager.splitter.setStyleSheet("""
            QSplitter::handle {
                background-color: #ccc;
                border: 1px solid #999;
                border-radius: 4px;
                margin: 2px;
            }
            QSplitter::handle:hover {
                background-color: #aaa;
            }
            QSplitter::handle:pressed {
                background-color: #888;
            }
        """)
        self.host.setCentralWidget(self.host.init_manager.splitter)

        left_pane = QWidget()
        left_pane.setAcceptDrops(True)
        left_layout = QVBoxLayout(left_pane)
        self._init_left_panel(left_layout)
        self.host.init_manager.splitter.addWidget(left_pane)

        right_pane = QWidget()
        right_layout = QVBoxLayout(right_pane)
        self._init_right_panel(right_layout)
        self.host.init_manager.splitter.addWidget(right_pane)

        # Monitor splitter movement
        self.host.init_manager.splitter.splitterMoved.connect(self.host.ui_manager.on_splitter_moved)
        self.host.init_manager.splitter.setSizes([600, 600])

        # Set tooltip for splitter handle
        QTimer.singleShot(100, self.host.ui_manager.setup_splitter_tooltip)

        # Settings to separate status bar segments
        self.status_bar = self.host.statusBar()
        self.host.init_manager.formula_label = QLabel("")  # Create label to be displayed on the right
        # Add margin to the right end for better appearance
        self.host.init_manager.formula_label.setStyleSheet("padding-right: 8px;")
        # Add label as a permanent widget on the right
        self.status_bar.addPermanentWidget(self.host.init_manager.formula_label)

    # --- Settings and Plugin Helpers ---
    def _get_default_settings(self):
        """Return a dictionary of default application settings."""
        return DEFAULT_SETTINGS.copy()

    def _migrate_legacy_settings(self, loaded_settings):
        """Handle migration of legacy settings keys."""
        # Migrate old key 'clean_up_2d_on_conversion'
        if "clean_up_2d_on_conversion" in loaded_settings:
            if loaded_settings["clean_up_2d_on_conversion"]:
                loaded_settings["3d_conversion_mode"] = "rdkit"
            del loaded_settings["clean_up_2d_on_conversion"]

        # Migrate old key 'use_obabel_optimization'
        if "use_obabel_optimization" in loaded_settings:
            if loaded_settings["use_obabel_optimization"]:
                loaded_settings["optimization_method"] = "MMFF94_OBABEL"
            del loaded_settings["use_obabel_optimization"]

    def _clear_plugin_ui_elements(self, plugin_menu):
        """Clean up tagged plugin actions from menus and toolbars."""
        # 1. Clear plugin-specific actions from Plugin menu (excluding the Manager)
        for action in plugin_menu.actions():
            if action.data() == "plugin_action":
                plugin_menu.removeAction(action)

        # 2. Clear plugin-specific toolbars or buttons
        if hasattr(self, "plugin_toolbar"):
            self.plugin_toolbar.clear()
            self.plugin_toolbar.hide()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(f"REPORT ERROR: Missing attribute 'plugin_toolbar' on self")

    def _init_left_panel(self, left_layout):
        """Initialize the left panel (2D view and buttons)."""
        self.host.init_manager.scene = MoleculeScene(self.host.state_manager.data, self.host)
        self.host.init_manager.scene.setSceneRect(-4000, -4000, 4000, 4000)
        self.host.init_manager.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.host.init_manager.view_2d = ZoomableView(self.host.init_manager.scene, self.host)
        self.host.init_manager.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.host.init_manager.view_2d.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        left_layout.addWidget(self.host.init_manager.view_2d, 1)

        self.host.init_manager.view_2d.scale(0.75, 0.75)

        # --- Left panel button layout ---
        left_buttons_layout = QHBoxLayout()
        self.host.init_manager.cleanup_button = QPushButton("Clean Up 2D")
        self.host.init_manager.cleanup_button.clicked.connect(self.host.edit_actions_manager.clean_up_2d_structure)
        left_buttons_layout.addWidget(self.host.init_manager.cleanup_button)

        self.host.init_manager.convert_button = QPushButton("Convert 2D to 3D")
        self.host.init_manager.convert_button.clicked.connect(self.host.compute_manager.trigger_conversion)
        # Allow right-click to open a temporary conversion-mode menu
        try:
            self.host.init_manager.convert_button.setContextMenuPolicy(
                Qt.ContextMenuPolicy.CustomContextMenu
            )
            self.host.init_manager.convert_button.customContextMenuRequested.connect(
                self.host.compute_manager.show_convert_menu
            )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

        left_buttons_layout.addWidget(self.host.init_manager.convert_button)
        left_layout.addLayout(left_buttons_layout)

    def _init_right_panel(self, right_layout):
        """Initialize the right panel (3D view and buttons)."""
        self.host.view_3d_manager.plotter = CustomQtInteractor(
            right_layout.parentWidget(), main_window=self.host, lighting="none"
        )
        self.host.view_3d_manager.plotter.setAcceptDrops(False)
        self.host.view_3d_manager.plotter.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.host.view_3d_manager.plotter.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)

        # 2. Add 3D view to layout
        right_layout.addWidget(self.host.view_3d_manager.plotter, 1)
        # self.host.view_3d_manager.plotter.installEventFilter(self)
        # 3. Create horizontal layout for buttons
        right_buttons_layout = QHBoxLayout()

        # 3D Optimize button
        self.host.init_manager.optimize_3d_button = QPushButton("Optimize 3D")
        self.host.init_manager.optimize_3d_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.host.init_manager.optimize_3d_button.clicked.connect(self.host.compute_manager.optimize_3d_structure)
        self.host.init_manager.optimize_3d_button.setEnabled(False)
        # Initialized via _enable_3d_features(False)
        # Allow right-click to open a temporary optimization-method menu
        try:
            self.host.init_manager.optimize_3d_button.setContextMenuPolicy(
                Qt.ContextMenuPolicy.CustomContextMenu
            )
            self.host.init_manager.optimize_3d_button.customContextMenuRequested.connect(
                self.host.compute_manager.show_optimize_menu
            )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

        right_buttons_layout.addWidget(self.host.init_manager.optimize_3d_button)

        # Export button with menu
        self.host.init_manager.export_button = QToolButton()
        self.host.init_manager.export_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.host.init_manager.export_button.setText("Export 3D")
        self.host.init_manager.export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.host.init_manager.export_button.setEnabled(False)  # Initially disabled

        export_menu = QMenu(self.host)
        export_mol_action = QAction("Export as MOL...", self.host)
        export_mol_action.triggered.connect(self.host.io_manager.save_3d_as_mol)
        export_menu.addAction(export_mol_action)

        export_xyz_action = QAction("Export as XYZ...", self.host)
        export_xyz_action.triggered.connect(self.host.io_manager.save_as_xyz)
        export_menu.addAction(export_xyz_action)

        export_png_action = QAction("Export as PNG...", self.host)
        export_png_action.triggered.connect(self.host.export_manager.export_3d_png)
        export_menu.addAction(export_png_action)

        self.host.init_manager.export_button.setMenu(export_menu)
        right_buttons_layout.addWidget(self.host.init_manager.export_button)

        # 4. Add horizontal layout to vertical layout
        right_layout.addLayout(right_buttons_layout)

    def _init_toolbars(self):
        """Initialize and organize toolbars."""
        # Row 1: Main Toolbar
        self.toolbar = QToolBar("Main Toolbar")
        self.host.addToolBar(self.toolbar)

        # Row 2: Templates
        with contextlib.suppress(AttributeError):
            self.host.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        self.toolbar_bottom = QToolBar("Templates Toolbar")
        self.host.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.toolbar_bottom)

        # Row 3: Plugins
        with contextlib.suppress(AttributeError):
            self.host.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        self.plugin_toolbar = QToolBar("Plugin Toolbar")
        self.host.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.plugin_toolbar)
        self.plugin_toolbar.hide()

    def _setup_action_groups(self, toolbar, toolbar_bottom):
        """Set up action groups and tool actions."""
        self.host.init_manager.tool_group = QActionGroup(self.host)
        self.host.init_manager.tool_group.setExclusive(True)

        self._add_atom_actions(toolbar)
        self._add_bond_actions(toolbar)
        self._add_charge_radical_actions(toolbar)
        self._add_template_actions(toolbar_bottom)
        self._add_3d_edit_actions(toolbar)

        # Set default tool
        self.host.ui_manager.set_mode("atom_C")
        if "atom_C" in self.host.init_manager.mode_actions:
            self.host.init_manager.mode_actions["atom_C"].setChecked(True)

    def _add_atom_actions(self, toolbar):
        """Add standard atom selection actions to the toolbar."""
        actions_data = [
            ("Select", "select", "Space"),
            ("C", "atom_C", "c"),
            ("H", "atom_H", "h"),
            ("B", "atom_B", "b"),
            ("N", "atom_N", "n"),
            ("O", "atom_O", "o"),
            ("S", "atom_S", "s"),
            ("Si", "atom_Si", "Shift+S"),
            ("P", "atom_P", "p"),
            ("F", "atom_F", "f"),
            ("Cl", "atom_Cl", "Shift+C"),
            ("Br", "atom_Br", "Shift+B"),
            ("I", "atom_I", "i"),
            ("Other...", "atom_other", ""),
        ]

        for text, mode, shortcut in actions_data:
            if text == "C":
                toolbar.addSeparator()
            action = QAction(text, self.host, checkable=(mode != "atom_other"))
            if shortcut:
                action.setToolTip(f"{text} ({shortcut})")

            if mode == "atom_other":
                action.triggered.connect(self.host.dialog_manager.open_periodic_table_dialog)
                self.host.init_manager.other_atom_action = action
            else:
                action.triggered.connect(lambda c, m=mode: self.host.ui_manager.set_mode(m))
                self.host.init_manager.mode_actions[mode] = action
                self.host.init_manager.tool_group.addAction(action)
            toolbar.addAction(action)
        toolbar.addSeparator()

    def _add_bond_actions(self, toolbar):
        """Add bond type selection actions with custom icons."""
        bonds = [
            ("Single Bond", "bond_1_0", "1", "single"),
            ("Double Bond", "bond_2_0", "2", "double"),
            ("Triple Bond", "bond_3_0", "3", "triple"),
            ("Wedge Bond", "bond_1_1", "W", "wedge"),
            ("Dash Bond", "bond_1_2", "D", "dash"),
            ("Toggle E/Z", "bond_2_5", "E/Z", "ez_toggle"),
        ]

        for text, mode, shortcut, icon_type in bonds:
            action = QAction(self.host)
            action.setIcon(self._create_bond_icon(icon_type))
            action.setToolTip(f"{text} ({shortcut})")
            action.setCheckable(True)
            action.triggered.connect(lambda checked, m=mode: self.host.ui_manager.set_mode(m))
            self.host.init_manager.mode_actions[mode] = action
            toolbar.addAction(action)
            self.host.init_manager.tool_group.addAction(action)
        toolbar.addSeparator()

    def _add_charge_radical_actions(self, toolbar):
        """Add charge and radical modification actions."""
        ops = [
            ("+ Charge", "charge_plus", "Increase Atom Charge (+)"),
            ("- Charge", "charge_minus", "Decrease Atom Charge (-)"),
            ("Radical", "radical", "Toggle Radical (0/1/2) (.)"),
        ]

        for text, mode, tooltip in ops:
            action = QAction(text, self.host, checkable=True)
            action.setToolTip(tooltip)
            action.triggered.connect(lambda c, m=mode: self.host.ui_manager.set_mode(m))
            self.host.init_manager.mode_actions[mode] = action
            toolbar.addAction(action)
            self.host.init_manager.tool_group.addAction(action)

    def _add_template_actions(self, toolbar_bottom):
        """Add structural template actions (rings, etc.) to the bottom toolbar."""
        toolbar_bottom.addWidget(QLabel(" Templates:"))
        templates = [("Benzene", "template_benzene", 6)] + [
            (f"{i}-Ring", f"template_{i}", i) for i in range(3, 10)
        ]

        for text, mode, n in templates:
            action = QAction(self.host)
            action.setCheckable(True)
            is_benzene = text == "Benzene"
            action.setIcon(self._create_template_icon(n, is_benzene=is_benzene))
            action.setToolTip(
                f"{text} Template (4)" if is_benzene else f"{text} Template"
            )
            action.triggered.connect(lambda c, m=mode: self.host.ui_manager.set_mode(m))
            self.host.init_manager.mode_actions[mode] = action
            toolbar_bottom.addAction(action)
            self.host.init_manager.tool_group.addAction(action)

        user_action = QAction("USER", self.host, checkable=True)
        user_action.setToolTip("Open User Templates Dialog")
        user_action.triggered.connect(self.host.dialog_manager.open_template_dialog_and_activate)
        self.host.init_manager.mode_actions["template_user"] = user_action
        toolbar_bottom.addAction(user_action)
        self.host.init_manager.tool_group.addAction(user_action)

    def _add_3d_edit_actions(self, toolbar):
        """Add 3D-specific selection and manipulation actions."""
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        toolbar.addWidget(spacer)

        self.host.init_manager.measurement_action = QAction("3D Select", self.host, checkable=True)
        self.host.init_manager.measurement_action.setToolTip(
            "Enable distance, angle, and dihedral measurement in 3D view"
        )
        self.host.init_manager.measurement_action.triggered.connect(self.host.edit_3d_manager.toggle_measurement_mode)
        toolbar.addAction(self.host.init_manager.measurement_action)

        self.host.init_manager.edit_3d_action = QAction("3D Drag", self.host, checkable=True)
        self.host.init_manager.edit_3d_action.setToolTip(
            "Toggle 3D atom dragging mode (Hold Alt for temporary mode)"
        )
        self.host.init_manager.edit_3d_action.toggled.connect(self.host.ui_manager.toggle_3d_edit_mode)
        toolbar.addAction(self.host.init_manager.edit_3d_action)

        self.host.init_manager.style_button = QToolButton()
        self.host.init_manager.style_button.setText("3D Style")
        self.host.init_manager.style_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(self.host.init_manager.style_button)

        style_menu = QMenu(self.host)
        self.host.init_manager.style_button.setMenu(style_menu)
        style_group = QActionGroup(self.host)
        style_group.setExclusive(True)

        for name, key in [
            ("Ball & Stick", "ball_and_stick"),
            ("CPK (Space-filling)", "cpk"),
            ("Wireframe", "wireframe"),
            ("Stick", "stick"),
        ]:
            action = QAction(name, self.host, checkable=True)
            if key == "ball_and_stick":
                action.setChecked(True)
            action.triggered.connect(
                lambda checked=False, k=key: (
                    self.host.view_3d_manager.set_3d_style(k),
                    self.host.view_3d_manager.draw_molecule_3d(self.host.view_3d_manager.current_mol)
                    if getattr(self.host.view_3d_manager, "current_mol", None)
                    else None,
                )
            )
            style_menu.addAction(action)
            style_group.addAction(action)

        if self.host.plugin_manager and self.host.plugin_manager.custom_3d_styles:
            style_menu.addSeparator()
            for style_name in self.host.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self.host, checkable=True)
                    action.triggered.connect(
                        lambda checked=False, s=style_name: (
                            self.host.view_3d_manager.set_3d_style(s),
                            self.host.view_3d_manager.draw_molecule_3d(self.host.view_3d_manager.current_mol)
                            if getattr(self.host.view_3d_manager, "current_mol", None)
                            else None,
                        )
                    )
                    style_menu.addAction(action)
                    style_group.addAction(action)

    def _create_bond_icon(self, bond_type, size=32):
        """Generate a QIcon for a specific bond type."""
        fg = self._get_icon_foreground_color()
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        p1, p2 = QPointF(6, size / 2), QPointF(size - 6, size / 2)
        line = QLineF(p1, p2)
        painter.setPen(QPen(fg, 2))
        painter.setBrush(QBrush(fg))

        if bond_type == "single":
            painter.drawLine(line)
        elif bond_type == "double":
            v = line.unitVector().normalVector()
            offset = QPointF(v.dx(), v.dy()) * 2.5
            painter.drawLine(line.translated(offset))
            painter.drawLine(line.translated(-offset))
        elif bond_type == "triple":
            v = line.unitVector().normalVector()
            offset = QPointF(v.dx(), v.dy()) * 3.0
            painter.drawLine(line)
            painter.drawLine(line.translated(offset))
            painter.drawLine(line.translated(-offset))
        elif bond_type == "wedge":
            normal = line.unitVector().normalVector()
            offset = QPointF(normal.dx(), normal.dy()) * 5.0
            painter.drawPolygon(QPolygonF([p1, p2 + offset, p2 - offset]))
        elif bond_type == "dash":
            normal = line.unitVector().normalVector()
            for i in range(NUM_DASHES + 1):
                t = i / NUM_DASHES
                start = p1 * (1 - t) + p2 * t
                offset = QPointF(normal.dx(), normal.dy()) * (10 * t) / 2.0
                painter.setPen(QPen(fg, 1.5))
                painter.drawLine(start - offset, start + offset)
        elif bond_type == "ez_toggle":
            p1z, p2z = QPointF(6, size * 0.75), QPointF(size - 6, size * 0.75)
            linez = QLineF(p1z, p2z)
            v = linez.unitVector().normalVector()
            offset = QPointF(v.dx(), v.dy()) * 2.0
            painter.drawLine(linez.translated(offset))
            painter.drawLine(linez.translated(-offset))
            font = painter.font()
            font.setPointSize(10)
            font.setBold(True)
            painter.setFont(font)
            painter.drawText(
                QRectF(0, 0, size, size * 0.6), Qt.AlignmentFlag.AlignCenter, "Z⇌E"
            )

        painter.end()
        return QIcon(pixmap)

    def _create_template_icon(self, n, is_benzene=False, size=32):
        """Generate a QIcon for a structural template (ring)."""
        fg = self._get_icon_foreground_color()
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(fg, 2))
        center = QPointF(size / 2, size / 2)
        radius = size / 2 - 4
        points = []
        angle_step = 2 * math.pi / n
        start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2
        for i in range(n):
            angle = start_angle + i * angle_step
            points.append(
                QPointF(
                    center.x() + radius * math.cos(angle),
                    center.y() + radius * math.sin(angle),
                )
            )
        painter.drawPolygon(QPolygonF(points))
        if is_benzene:
            painter.drawEllipse(center, radius * 0.6, radius * 0.6)
        if n in [7, 8, 9]:
            painter.setFont(QFont("Arial", 10, QFont.Weight.Bold))
            painter.setPen(QPen(fg, 1))
            painter.drawText(
                QRectF(0, 0, size, size), Qt.AlignmentFlag.AlignCenter, str(n)
            )
        painter.end()
        return QIcon(pixmap)

    def _get_icon_foreground_color(self):
        """Determine appropriate icon foreground color based on theme/settings."""
        with contextlib.suppress(Exception):
            fg = self.host.init_manager.settings.get("icon_foreground")
            if fg and QColor(fg).isValid():
                return QColor(fg)

        with contextlib.suppress(Exception):
            os_pref = detect_system_dark_mode()
            if os_pref is not None:
                return QColor("#FFFFFF") if os_pref else QColor("#000000")

        with contextlib.suppress(Exception):
            bg = QColor(self.host.init_manager.settings.get("background_color", "#919191"))
            if bg.isValid():
                lum = 0.2126 * bg.redF() + 0.7152 * bg.greenF() + 0.0722 * bg.blueF()
                return QColor("#FFFFFF") if lum < 0.5 else QColor("#000000")
        return QColor("#000000")

    # --- Menu Bar Initialization Helpers ---
    def _init_file_menu(self, menu_bar):
        """Initialize the File menu."""
        file_menu = menu_bar.addMenu("&File")

        # === Project Operations ===
        new_action = QAction("&New", self.host)
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.host.edit_actions_manager.clear_all)
        file_menu.addAction(new_action)

        load_project_action = QAction("&Open Project...", self.host)
        load_project_action.setShortcut("Ctrl+O")
        load_project_action.triggered.connect(self.host.io_manager.open_project)
        file_menu.addAction(load_project_action)

        save_action = QAction("&Save Project", self.host)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.host.io_manager.save_project)
        file_menu.addAction(save_action)

        save_as_action = QAction("Save Project &As...", self.host)
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.host.io_manager.save_project_as)
        file_menu.addAction(save_as_action)

        save_template_action = QAction("Save 2D as Template...", self.host)
        save_template_action.triggered.connect(self.host.dialog_manager.save_2d_as_template)
        file_menu.addAction(save_template_action)

        file_menu.addSeparator()

        # === Import ===
        self.import_menu = file_menu.addMenu("Import")
        load_mol_action = QAction("MOL/SDF File...", self.host)
        load_mol_action.triggered.connect(self.host.io_manager.load_mol_file)
        self.import_menu.addAction(load_mol_action)

        import_smiles_action = QAction("SMILES...", self.host)
        import_smiles_action.triggered.connect(self.host.string_importer_manager.import_smiles_dialog)
        self.import_menu.addAction(import_smiles_action)

        import_inchi_action = QAction("InChI...", self.host)
        import_inchi_action.triggered.connect(self.host.string_importer_manager.import_inchi_dialog)
        self.import_menu.addAction(import_inchi_action)

        self.import_menu.addSeparator()
        load_3d_mol_action = QAction("3D MOL/SDF (3D View Only)...", self.host)
        load_3d_mol_action.triggered.connect(self.host.io_manager.load_mol_file_for_3d_viewing)
        self.import_menu.addAction(load_3d_mol_action)

        load_3d_xyz_action = QAction("3D XYZ (3D View Only)...", self.host)
        load_3d_xyz_action.triggered.connect(self.host.io_manager.load_xyz_for_3d_viewing)
        self.import_menu.addAction(load_3d_xyz_action)

        # === Export ===
        export_menu = file_menu.addMenu("Export")
        export_pmeraw_action = QAction("PME Raw Format...", self.host)
        export_pmeraw_action.triggered.connect(self.host.io_manager.save_raw_data)
        export_menu.addAction(export_pmeraw_action)
        export_menu.addSeparator()

        export_2d_menu = export_menu.addMenu("2D Formats")
        save_mol_action = QAction("MOL File...", self.host)
        save_mol_action.triggered.connect(self.host.io_manager.save_as_mol)
        export_2d_menu.addAction(save_mol_action)

        export_2d_png_action = QAction("PNG Image...", self.host)
        export_2d_png_action.triggered.connect(self.host.export_manager.export_2d_png)
        export_2d_menu.addAction(export_2d_png_action)

        export_2d_svg_action = QAction("SVG Image...", self.host)
        export_2d_svg_action.triggered.connect(self.host.export_manager.export_2d_svg)
        export_2d_menu.addAction(export_2d_svg_action)

        export_3d_menu = export_menu.addMenu("3D Formats")
        save_3d_mol_action = QAction("MOL File...", self.host)
        save_3d_mol_action.triggered.connect(self.host.io_manager.save_3d_as_mol)
        export_3d_menu.addAction(save_3d_mol_action)

        save_xyz_action = QAction("XYZ File...", self.host)
        save_xyz_action.triggered.connect(self.host.io_manager.save_as_xyz)
        export_3d_menu.addAction(save_xyz_action)

        export_3d_png_action = QAction("PNG Image...", self.host)
        export_3d_png_action.triggered.connect(self.host.export_manager.export_3d_png)
        export_3d_menu.addAction(export_3d_png_action)
        export_3d_menu.addSeparator()

        export_stl_action = QAction("STL File...", self.host)
        export_stl_action.triggered.connect(self.host.export_manager.export_stl)
        export_3d_menu.addAction(export_stl_action)

        export_obj_action = QAction("OBJ/MTL (with colors)...", self.host)
        export_obj_action.triggered.connect(self.host.export_manager.export_obj_mtl)
        export_3d_menu.addAction(export_obj_action)

        file_menu.addSeparator()
        quit_action = QAction("Quit", self.host)
        quit_action.triggered.connect(self.host.close)
        file_menu.addAction(quit_action)

    def _init_edit_menu(self, menu_bar):
        """Initialize the Edit menu."""
        edit_menu = menu_bar.addMenu("&Edit")
        self.host.init_manager.undo_action = QAction("Undo", self.host)
        self.host.init_manager.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.host.init_manager.undo_action.triggered.connect(self.host.edit_actions_manager.undo)
        edit_menu.addAction(self.host.init_manager.undo_action)

        self.host.init_manager.redo_action = QAction("Redo", self.host)
        self.host.init_manager.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.host.init_manager.redo_action.triggered.connect(self.host.edit_actions_manager.redo)
        edit_menu.addAction(self.host.init_manager.redo_action)

        edit_menu.addSeparator()
        self.host.init_manager.cut_action = QAction("Cut", self.host)
        self.host.init_manager.cut_action.setShortcut(QKeySequence.StandardKey.Cut)
        self.host.init_manager.cut_action.triggered.connect(self.host.edit_actions_manager.cut_selection)
        edit_menu.addAction(self.host.init_manager.cut_action)

        self.host.init_manager.copy_action = QAction("Copy", self.host)
        self.host.init_manager.copy_action.setShortcut(QKeySequence.StandardKey.Copy)
        self.host.init_manager.copy_action.triggered.connect(self.host.edit_actions_manager.copy_selection)
        edit_menu.addAction(self.host.init_manager.copy_action)

        self.host.init_manager.paste_action = QAction("Paste", self.host)
        self.host.init_manager.paste_action.setShortcut(QKeySequence.StandardKey.Paste)
        self.host.init_manager.paste_action.triggered.connect(self.host.edit_actions_manager.paste_from_clipboard)
        edit_menu.addAction(self.host.init_manager.paste_action)

        edit_menu.addSeparator()
        add_hydrogen_action = QAction("Add Hydrogens", self.host)
        add_hydrogen_action.setToolTip(
            "Add explicit hydrogens based on RDKit implicit counts"
        )
        add_hydrogen_action.triggered.connect(self.host.edit_actions_manager.add_hydrogen_atoms)
        edit_menu.addAction(add_hydrogen_action)

        remove_hydrogen_action = QAction("Remove Hydrogens", self.host)
        remove_hydrogen_action.triggered.connect(self.host.edit_actions_manager.remove_hydrogen_atoms)
        edit_menu.addAction(remove_hydrogen_action)

        edit_menu.addSeparator()
        rotate_2d_action = QAction("Rotate 2D...", self.host)
        rotate_2d_action.setShortcut(QKeySequence("Ctrl+R"))
        rotate_2d_action.triggered.connect(self.host.edit_actions_manager.open_rotate_2d_dialog)
        edit_menu.addAction(rotate_2d_action)

        edit_menu.addSeparator()
        optimize_2d_action = QAction("Clean Up 2D", self.host)
        optimize_2d_action.setShortcut(QKeySequence("Ctrl+J"))
        optimize_2d_action.triggered.connect(self.host.edit_actions_manager.clean_up_2d_structure)
        edit_menu.addAction(optimize_2d_action)

        convert_3d_action = QAction("Convert 2D to 3D", self.host)
        convert_3d_action.setShortcut(QKeySequence("Ctrl+K"))
        convert_3d_action.triggered.connect(self.host.compute_manager.trigger_conversion)
        edit_menu.addAction(convert_3d_action)

        optimize_3d_action = QAction("Optimize 3D", self.host)
        optimize_3d_action.setShortcut(QKeySequence("Ctrl+L"))
        optimize_3d_action.triggered.connect(self.host.compute_manager.optimize_3d_structure)
        edit_menu.addAction(optimize_3d_action)

        edit_menu.addSeparator()
        select_all_action = QAction("Select All", self.host)
        select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.host.edit_actions_manager.select_all)
        edit_menu.addAction(select_all_action)

        clear_all_action = QAction("Clear All", self.host)
        clear_all_action.setShortcut(QKeySequence("Ctrl+Shift+C"))
        clear_all_action.triggered.connect(self.host.edit_actions_manager.clear_all)
        edit_menu.addAction(clear_all_action)

    def _init_view_menu(self, menu_bar):
        """Initialize the View menu."""
        view_menu = menu_bar.addMenu("&View")
        zoom_in_action = QAction("Zoom In", self.host)
        zoom_in_action.setShortcut(QKeySequence.StandardKey.ZoomIn)
        zoom_in_action.triggered.connect(self.host.view_3d_manager.zoom_in)
        view_menu.addAction(zoom_in_action)

        zoom_out_action = QAction("Zoom Out", self.host)
        zoom_out_action.setShortcut(QKeySequence.StandardKey.ZoomOut)
        zoom_out_action.triggered.connect(self.host.view_3d_manager.zoom_out)
        view_menu.addAction(zoom_out_action)

        reset_zoom_action = QAction("Reset Zoom", self.host)
        reset_zoom_action.setShortcut(QKeySequence("Ctrl+0"))
        reset_zoom_action.triggered.connect(self.host.view_3d_manager.reset_zoom)
        view_menu.addAction(reset_zoom_action)

        fit_action = QAction("Fit to View", self.host)
        fit_action.setShortcut(QKeySequence("Ctrl+9"))
        fit_action.triggered.connect(self.host.view_3d_manager.fit_to_view)
        view_menu.addAction(fit_action)

        view_menu.addSeparator()
        reset_3d_view_action = QAction("Reset 3D View", self.host)
        reset_3d_view_action.triggered.connect(
            lambda: self.host.view_3d_manager.plotter.reset_camera() if hasattr(self.host.view_3d_manager, 'plotter') else None
        )
        reset_3d_view_action.setShortcut(QKeySequence("Ctrl+Shift+R"))
        view_menu.addAction(reset_3d_view_action)

        redraw_menu_action = QAction("Redraw 3D Molecule", self.host)
        redraw_menu_action.triggered.connect(self.host.edit_actions_manager.redraw_molecule_3d)
        self.host.redraw_menu_action = redraw_menu_action
        view_menu.addAction(self.host.redraw_menu_action)

        view_menu.addSeparator()
        layout_menu = view_menu.addMenu("Panel Layout")
        equal_panels_action = QAction("Equal Panels (50:50)", self.host)
        equal_panels_action.setShortcut(QKeySequence("Ctrl+1"))
        equal_panels_action.triggered.connect(lambda: self.host.ui_manager.set_panel_layout(50, 50))
        layout_menu.addAction(equal_panels_action)

        layout_2d_focus_action = QAction("2D Focus (70:30)", self.host)
        layout_2d_focus_action.setShortcut(QKeySequence("Ctrl+2"))
        layout_2d_focus_action.triggered.connect(lambda: self.host.ui_manager.set_panel_layout(70, 30))
        layout_menu.addAction(layout_2d_focus_action)

        layout_3d_focus_action = QAction("3D Focus (30:70)", self.host)
        layout_3d_focus_action.setShortcut(QKeySequence("Ctrl+3"))
        layout_3d_focus_action.triggered.connect(lambda: self.host.ui_manager.set_panel_layout(30, 70))
        layout_menu.addAction(layout_3d_focus_action)

        layout_menu.addSeparator()
        toggle_2d_panel_action = QAction("Toggle 2D Panel", self.host)
        toggle_2d_panel_action.setShortcut(QKeySequence("Ctrl+H"))
        toggle_2d_panel_action.triggered.connect(self.host.ui_manager.toggle_2d_panel)
        layout_menu.addAction(toggle_2d_panel_action)

        view_menu.addSeparator()
        self.toggle_chiral_action = QAction("Show Chiral Labels", self.host, checkable=True)
        self.toggle_chiral_action.setChecked(self.host.view_3d_manager.show_chiral_labels)
        self.toggle_chiral_action.triggered.connect(self.host.view_3d_manager.toggle_chiral_labels_display)
        view_menu.addAction(self.toggle_chiral_action)

        view_menu.addSeparator()
        atom_info_menu = view_menu.addMenu("3D Atom Info Display")
        self.show_atom_id_action = QAction(
            "Show Original ID / Index", self.host, checkable=True
        )
        self.show_atom_id_action.triggered.connect(
            lambda: self.host.view_3d_manager.toggle_atom_info_display("id")
        )
        atom_info_menu.addAction(self.show_atom_id_action)

        self.show_rdkit_id_action = QAction("Show RDKit Index", self.host, checkable=True)
        self.show_rdkit_id_action.triggered.connect(
            lambda: self.host.view_3d_manager.toggle_atom_info_display("rdkit_id")
        )
        atom_info_menu.addAction(self.show_rdkit_id_action)

        self.show_atom_coords_action = QAction(
            "Show Coordinates (X,Y,Z)", self.host, checkable=True
        )
        self.show_atom_coords_action.triggered.connect(
            lambda: self.host.view_3d_manager.toggle_atom_info_display("coords")
        )
        atom_info_menu.addAction(self.show_atom_coords_action)

        self.show_atom_symbol_action = QAction(
            "Show Element Symbol", self.host, checkable=True
        )
        self.show_atom_symbol_action.triggered.connect(
            lambda: self.host.view_3d_manager.toggle_atom_info_display("symbol")
        )
        atom_info_menu.addAction(self.show_atom_symbol_action)

    def _init_analysis_menu(self, menu_bar):
        """Initialize the Analysis menu."""
        analysis_menu = menu_bar.addMenu("&Analysis")
        self.host.init_manager.analysis_action = QAction("Show Analysis...", self.host)
        self.host.init_manager.analysis_action.triggered.connect(self.host.dialog_manager.open_analysis_window)
        self.host.init_manager.analysis_action.setEnabled(False)
        analysis_menu.addAction(self.host.init_manager.analysis_action)

    def _init_edit_3d_menu(self, menu_bar):
        """Initialize the 3D Edit menu."""
        edit_3d_menu = menu_bar.addMenu("3D &Edit")
        translation_action = QAction("Translation...", self.host)
        translation_action.triggered.connect(self.host.dialog_manager.open_translation_dialog)
        translation_action.setEnabled(False)
        edit_3d_menu.addAction(translation_action)
        self.host.translation_action = translation_action

        move_group_action = QAction("Move Group...", self.host)
        move_group_action.triggered.connect(self.host.dialog_manager.open_move_group_dialog)
        move_group_action.setEnabled(False)
        edit_3d_menu.addAction(move_group_action)
        self.host.move_group_action = move_group_action

        edit_3d_menu.addSeparator()
        align_menu = edit_3d_menu.addMenu("Align to")
        align_menu.setEnabled(False)
        self.host.align_menu = align_menu

        axis_align_menu = align_menu.addMenu("Axis")
        align_x_action = QAction("X-axis", self.host)
        align_x_action.triggered.connect(lambda: self.host.dialog_manager.open_alignment_dialog("x"))
        align_x_action.setEnabled(False)
        axis_align_menu.addAction(align_x_action)
        self.host.align_x_action = align_x_action

        align_y_action = QAction("Y-axis", self.host)
        align_y_action.triggered.connect(lambda: self.host.dialog_manager.open_alignment_dialog("y"))
        align_y_action.setEnabled(False)
        axis_align_menu.addAction(align_y_action)
        self.host.align_y_action = align_y_action

        align_z_action = QAction("Z-axis", self.host)
        align_z_action.triggered.connect(lambda: self.host.dialog_manager.open_alignment_dialog("z"))
        align_z_action.setEnabled(False)
        axis_align_menu.addAction(align_z_action)
        self.host.align_z_action = align_z_action

        plane_align_menu = align_menu.addMenu("Plane")
        alignplane_xy_action = QAction("XY-plane", self.host)
        alignplane_xy_action.triggered.connect(
            lambda: self.host.dialog_manager.open_align_plane_dialog("xy")
        )
        alignplane_xy_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xy_action)
        self.host.alignplane_xy_action = alignplane_xy_action

        alignplane_xz_action = QAction("XZ-plane", self.host)
        alignplane_xz_action.triggered.connect(
            lambda: self.host.dialog_manager.open_align_plane_dialog("xz")
        )
        alignplane_xz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xz_action)
        self.host.alignplane_xz_action = alignplane_xz_action

        alignplane_yz_action = QAction("YZ-plane", self.host)
        alignplane_yz_action.triggered.connect(
            lambda: self.host.dialog_manager.open_align_plane_dialog("yz")
        )
        alignplane_yz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_yz_action)
        self.host.alignplane_yz_action = alignplane_yz_action

        edit_3d_menu.addSeparator()
        mirror_action = QAction("Mirror...", self.host)
        mirror_action.triggered.connect(self.host.dialog_manager.open_mirror_dialog)
        mirror_action.setEnabled(False)
        edit_3d_menu.addAction(mirror_action)
        self.host.mirror_action = mirror_action

        edit_3d_menu.addSeparator()
        planarize_action = QAction("Planarize...", self.host)
        planarize_action.triggered.connect(lambda: self.host.dialog_manager.open_planarize_dialog(None))
        planarize_action.setEnabled(False)
        edit_3d_menu.addAction(planarize_action)
        self.host.planarize_action = planarize_action

        edit_3d_menu.addSeparator()
        bond_length_action = QAction("Adjust Bond Length...", self.host)
        bond_length_action.triggered.connect(self.host.dialog_manager.open_bond_length_dialog)
        bond_length_action.setEnabled(False)
        edit_3d_menu.addAction(bond_length_action)
        self.host.bond_length_action = bond_length_action

        angle_action = QAction("Adjust Angle...", self.host)
        angle_action.triggered.connect(self.host.dialog_manager.open_angle_dialog)
        angle_action.setEnabled(False)
        edit_3d_menu.addAction(angle_action)
        self.host.angle_action = angle_action

        dihedral_action = QAction("Adjust Dihedral Angle...", self.host)
        dihedral_action.triggered.connect(self.host.dialog_manager.open_dihedral_dialog)
        dihedral_action.setEnabled(False)
        edit_3d_menu.addAction(dihedral_action)
        self.host.dihedral_action = dihedral_action

        edit_3d_menu.addSeparator()
        constrained_opt_action = QAction("Constrained Optimization...", self.host)
        constrained_opt_action.triggered.connect(
            self.host.dialog_manager.open_constrained_optimization_dialog
        )
        constrained_opt_action.setEnabled(False)
        edit_3d_menu.addAction(constrained_opt_action)
        self.host.constrained_opt_action = constrained_opt_action

    def _init_plugin_menu(self, menu_bar):
        """Initialize the Plugin menu."""
        self.plugin_menu = menu_bar.addMenu("&Plugin")
        manage_plugins_action = QAction("Plugin Manager...", self.host)

        def show_plugin_manager():
            if not self.host.plugin_manager:
                QMessageBox.information(
                    self.host, "Safe Mode", "Plugins are disabled (safe mode)."
                )
                return
            from ..plugins.plugin_manager_window import PluginManagerWindow

            dlg = PluginManagerWindow(self.host.plugin_manager, self.host)
            dlg.exec()
            self.update_plugin_menu(self.plugin_menu)

        manage_plugins_action.triggered.connect(show_plugin_manager)
        self.plugin_menu.addAction(manage_plugins_action)
        self.plugin_menu.addSeparator()

    def _init_settings_menu(self, menu_bar):
        """Initialize the Settings menu."""
        settings_menu = menu_bar.addMenu("&Settings")
        view_settings_action = QAction("Settings...", self.host)
        view_settings_action.triggered.connect(self.host.dialog_manager.open_settings_dialog)
        settings_menu.addAction(view_settings_action)

        color_action = QAction("CPK Colors...", self.host)
        color_action.triggered.connect(self.host.dialog_manager.open_color_settings_dialog)
        settings_menu.addAction(color_action)

        conversion_menu = settings_menu.addMenu("3D Conversion")
        conv_group = QActionGroup(self.host)
        conv_group.setExclusive(True)

        def _set_conv_mode(mode):
            try:
                self.host.init_manager.settings["3d_conversion_mode"] = mode
                self.host.init_manager.settings_dirty = True
                self.host.statusBar().showMessage(f"3D conversion mode set to: {mode}")
            except Exception as e:
                logging.debug(f"Suppressed exception: {e}")

        conv_options = [
            ("Fallback", "fallback"),
            ("RDKit only", "rdkit"),
            ("Open Babel only", "obabel"),
            ("Direct (use 2D coords + add H)", "direct"),
        ]
        self.conv_actions = {}
        for label, key in conv_options:
            a = QAction(label, self.host)
            a.setCheckable(True)
            if key == "obabel" and not OBABEL_AVAILABLE:
                a.setEnabled(False)
            a.triggered.connect(lambda checked, m=key: _set_conv_mode(m))
            conversion_menu.addAction(a)
            conv_group.addAction(a)
            self.conv_actions[key] = a

        saved_conv = self.host.init_manager.settings.get("3d_conversion_mode", "fallback")
        if (
            saved_conv not in self.conv_actions
            or not self.conv_actions[saved_conv].isEnabled()
        ):
            saved_conv = (
                "rdkit"
                if self.conv_actions.get("rdkit", QAction(self.host)).isEnabled()
                else "fallback"
            )

        if saved_conv in self.conv_actions:
            self.conv_actions[saved_conv].setChecked(True)
        self.host.init_manager.settings["3d_conversion_mode"] = saved_conv

        optimization_menu = settings_menu.addMenu("3D Optimization Settings")
        opt_methods = [
            ("MMFF94s (RDKit)", "MMFF_RDKIT"),
            ("MMFF94 (RDKit)", "MMFF94_RDKIT"),
            ("UFF (RDKit)", "UFF_RDKIT"),
            ("MMFF94s (Open Babel)", "MMFF94s_OBABEL"),
            ("MMFF94 (Open Babel)", "MMFF94_OBABEL"),
            ("UFF (Open Babel)", "UFF_OBABEL"),
            ("GAFF (Open Babel)", "GAFF_OBABEL"),
            ("Ghemical (Open Babel)", "GHEMICAL_OBABEL"),
        ]
        self.opt3d_method_labels = {key.upper(): label for (label, key) in opt_methods}
        opt_group = QActionGroup(self.host)
        opt_group.setExclusive(True)
        self.opt3d_actions = {}
        for label, key in opt_methods:
            action = QAction(label, self.host)
            action.setCheckable(True)
            if key.endswith("_OBABEL") and not OBABEL_AVAILABLE:
                action.setEnabled(False)
            action.triggered.connect(
                lambda checked, m=key: self.host.compute_manager.set_optimization_method(m)
            )
            optimization_menu.addAction(action)
            opt_group.addAction(action)
            self.opt3d_actions[key] = action

        optimization_menu.addSeparator()
        self.host.intermolecular_rdkit_action = QAction(
            "Consider Intermolecular Interaction for RDKit", self.host
        )
        self.host.intermolecular_rdkit_action.setCheckable(True)
        self.host.intermolecular_rdkit_action.setChecked(
            self.host.init_manager.settings.get("optimize_intermolecular_interaction_rdkit", True)
        )
        self.host.intermolecular_rdkit_action.triggered.connect(
            self.host.compute_manager.toggle_intermolecular_interaction_rdkit
        )
        optimization_menu.addAction(self.host.intermolecular_rdkit_action)

        saved_opt = (self.host.init_manager.settings.get("optimization_method") or "MMFF_RDKIT").upper()
        if (
            saved_opt in self.opt3d_actions
            and self.opt3d_actions[saved_opt].isEnabled()
        ):
            self.opt3d_actions[saved_opt].setChecked(True)
            self.optimization_method = saved_opt
        else:
            self.opt3d_actions["MMFF_RDKIT"].setChecked(True)
            self.optimization_method = "MMFF_RDKIT"

        settings_menu.addSeparator()
        reset_settings_action = QAction("Reset All Settings", self.host)
        reset_settings_action.triggered.connect(self.reset_all_settings_menu)
        settings_menu.addAction(reset_settings_action)

    def _init_help_menu(self, menu_bar):
        """Initialize the Help menu."""
        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("About", self.host)
        about_action.triggered.connect(self.host.dialog_manager.show_about_dialog)
        help_menu.addAction(about_action)

        github_action = QAction("GitHub", self.host)
        github_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl("https://github.com/HiroYokoyama/python_molecular_editor")
            )
        )
        help_menu.addAction(github_action)

        github_wiki_action = QAction("GitHub Wiki", self.host)
        github_wiki_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl("https://github.com/HiroYokoyama/python_molecular_editor/wiki")
            )
        )
        help_menu.addAction(github_wiki_action)

        manual_action = QAction("User Manual", self.host)
        manual_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl(
                    "https://hiroyokoyama.github.io/python_molecular_editor/manual/manual"
                )
            )
        )
        help_menu.addAction(manual_action)
