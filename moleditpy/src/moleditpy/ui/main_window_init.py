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
except ImportError:
    from moleditpy import OBABEL_AVAILABLE


try:
    import winreg
except ImportError:
    winreg = None

try:
    from ..plugins.plugin_manager import PluginManager
except ImportError:
    from moleditpy.plugins.plugin_manager import PluginManager


try:
    from ..utils.system_utils import detect_system_dark_mode
except ImportError:
    from moleditpy.utils.system_utils import detect_system_dark_mode


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


# --- Class Definition ---
class MainWindowMainInit:
    """Feature class separated from main_window.py"""

    def __init__(self, initial_file=None, safe_mode=False):
        # This helper is not used as a mixin in this project; initialization
        # happens on the `QMainWindow` base class in
        # `MainWindow.__init__` directly.
        self.setAcceptDrops(True)
        self.settings_dir = os.path.join(os.path.expanduser("~"), ".moleditpy")
        self.settings_file = os.path.join(self.settings_dir, "settings.json")
        self.settings = {}  # Will be populated by load_settings
        self.load_settings()
        self.initial_settings = self.settings.copy()
        self.setWindowTitle("MoleditPy Ver. " + VERSION)
        self.setGeometry(100, 100, 1400, 800)
        self.data = MolecularData()
        self.current_mol = None
        # self.current_3d_style is now in View3DManager
        # self.atom_info_display_mode is now in View3DManager
        # self.current_atom_info_labels is now in View3DManager
        # self.is_3d_edit_mode is now in Edit3DManager
        # self.dragged_atom_info is now in Edit3DManager
        # self.atom_actor is now in View3DManager
        self.is_2d_editable = True
        self.is_xyz_derived = (
            False  # Flag indicating if the molecule is derived from XYZ
        )
        # Chemical check flags: whether a chemical/sanitization check was attempted and whether it failed
        self.chem_check_tried = False
        self.chem_check_failed = False
        # self.axes_actor and self.axes_widget are now in View3DManager
        self._template_dialog = None  # Reference to the template dialog
        self.undo_stack = []
        self.redo_stack = []
        # self.constraints_3d is now in Edit3DManager
        self.mode_actions = {}

        # Variable tracking the saved state
        self.has_unsaved_changes = False
        # Flag to delay disk writes of the settings file
        # If set to True, settings are updated in memory and saved together upon application exit.
        self.settings_dirty = True
        self.current_file_path = None  # Path of the currently open file
        self.initialization_complete = False  # Initialization completion flag
        # Token to invalidate pending implicit-hydrogen UI updates
        self._ih_update_counter = 0

        # Variables for measurement functionality and 3D selection are now managed by Edit3DManager
        # self.measurement_mode is now in Edit3DManager
        # self.selected_atoms_for_measurement is now in Edit3DManager
        # self.measurement_labels is now in Edit3DManager
        # self.measurement_text_actor is now in Edit3DManager
        # self.measurement_label_items_2d is now in Edit3DManager
        # self.selected_atoms_3d is now in Edit3DManager
        # self.active_3d_dialogs is now in Edit3DManager

        # Initialization of the plugin manager
        if safe_mode:
            print("Safe mode: plugins disabled.")
            self.plugin_manager = None
        else:
            try:
                self.plugin_manager = PluginManager()
            except (AttributeError, RuntimeError, ValueError) as e:
                print(f"Failed to initialize PluginManager: {e}")
                self.plugin_manager = None

        # Dictionary holding data for plugins that haven't been loaded
        self._preserved_plugin_data = {}

        self.init_ui()
        self.init_worker_thread()
        self._setup_3d_picker()

        # --- RDKit Warm-up (initial execution cost) ---
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

        self.reset_undo_stack()
        self.scene.selectionChanged.connect(self.update_edit_menu_actions)
        QApplication.clipboard().dataChanged.connect(self.update_edit_menu_actions)

        self.update_edit_menu_actions()

        if initial_file:
            self.load_command_line_file(initial_file)

        QTimer.singleShot(0, self.apply_initial_settings)
        # Camera initialization flag (permits reset only during the first draw)
        self._camera_initialized = False

        # Set initial menu text and state
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

        # Set initialization complete
        self.initialization_complete = True
        self.update_window_title()  # Update title after initialization is complete
        # Ensure initial keyboard/mouse focus is placed on the 2D view
        # when opening a file or starting the application. This avoids
        # accidental focus landing on toolbar/buttons (e.g. Optimize 2D).
        try:
            QTimer.singleShot(0, self.view_2d.setFocus)
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

    def init_ui(self):
        # 1. Get the path to the directory where the current script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # 2. Build the full path to the icon file in the 'assets' folder
        # Windows taskbar prefers .ico; fall back to .png on other platforms
        if sys.platform == "win32":
            icon_path = os.path.join(script_dir, "..", "assets", "icon.ico")
        else:
            icon_path = os.path.join(script_dir, "..", "assets", "icon.png")

        # 3. Create a QIcon object directly from the file path
        if os.path.exists(icon_path):  # Check if file exists
            app_icon = QIcon(icon_path)

            # 4. Set the icon for both the window and the application
            self.setWindowIcon(app_icon)
            QApplication.instance().setWindowIcon(app_icon)
        else:
            print(f"Warning: Icon file not found: {icon_path}")

        self._init_main_layout()
        self._init_toolbars()
        self.init_menu_bar()

        self._setup_action_groups(self.toolbar, self.toolbar_bottom)

    def init_menu_bar(self):
        menu_bar = self.menuBar()

        self._init_file_menu(menu_bar)
        self._init_edit_menu(menu_bar)
        self._init_view_menu(menu_bar)
        self._init_analysis_menu(menu_bar)
        self._init_edit_3d_menu(menu_bar)
        self._init_plugin_menu(menu_bar)
        self._init_settings_menu(menu_bar)
        self._init_help_menu(menu_bar)

        # Consistently set initial state for 3D-related features
        self._enable_3d_features(False)

        # Finally, populate plugins now that all menus are created
        self.update_plugin_menu(self.plugin_menu)

    def init_worker_thread(self):
        # Initialize shared state for calculation runs.
        self.halt_ids = set()
        self.next_conversion_id = 1
        self.active_worker_ids = set()
        # Track active threads for diagnostics/cleanup (weak references ok)
        try:
            self._active_calc_threads = []
        except (AttributeError, RuntimeError, ValueError, TypeError):
            self._active_calc_threads = []

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
        if ext_with_dot in self.plugin_manager.file_openers:
            openers = self.plugin_manager.file_openers[ext_with_dot]
            # Iterate through openers (already sorted by priority)
            for opener_info in openers:
                try:
                    callback = opener_info["callback"]
                    # Try to call the opener
                    callback(file_path)

                    self.current_file_path = file_path
                    self.update_window_title()
                    return  # Success
                except (AttributeError, RuntimeError, ValueError) as e:
                    print(
                        f"Plugin opener failed for '{opener_info.get('plugin', 'Unknown')}': {e}"
                    )
                    # If this opener fails, try the next one or fall through to default
                    continue

        if file_ext in ["mol", "sdf"]:
            self.load_mol_file_for_3d_viewing(file_path)
        elif file_ext == "xyz":
            self.load_xyz_for_3d_viewing(file_path)
        elif file_ext in ["pmeraw", "pmeprj"]:
            self.open_project_file(file_path=file_path)
        else:
            self.statusBar().showMessage(f"Unsupported file type: {file_ext}")

    def apply_initial_settings(self):
        """Apply saved settings to the 3D view after UI initialization is complete."""

        try:
            self.update_cpk_colors_from_settings()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu/settings sync errors

        if self.plotter and self.plotter.renderer:
            bg_color = self.settings.get("background_color", "#919191")
            self.plotter.set_background(bg_color)
            self.apply_3d_settings()

        try:
            if hasattr(self, "scene") and self.scene:
                # Apply 2D background color
                bg_color_2d = self.settings.get("background_color_2d", "#FFFFFF")
                self.scene.setBackgroundBrush(QBrush(QColor(bg_color_2d)))

                for it in list(self.scene.items()):
                    if hasattr(it, "update_style"):
                        it.update_style()
                self.scene.update()
                for v in list(self.scene.views()):
                    v.viewport().update()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu/settings sync errors

    def update_cpk_colors_from_settings(self):
        """Update global CPK_COLORS and CPK_COLORS_PV from saved settings overrides.

        This modifies the in-memory CPK_COLORS mapping (not persisted until settings are saved).
        Only keys present in self.settings['cpk_colors'] are changed; other elements keep the defaults.
        """
        try:
            # Overridden CPK settings are stored in self.settings['cpk_colors'].
            # To ensure that 2D modules (e.g., atom_item.py) which imported the
            # `CPK_COLORS` mapping from `moleditpy.utils.constants` at import time see
            # updates, mutate the mapping in-place on the constants module
            # instead of rebinding a new local variable here.
            overrides = self.settings.get("cpk_colors", {}) or {}

            # Import the constants module so we can update mappings directly
            try:
                from . import constants as constants_mod
            except ImportError:
                import moleditpy.utils.constants as constants_mod

            # Reset constants.CPK_COLORS to defaults but keep the same dict
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
            constants_mod.CPK_COLORS_PV.clear()
            for k, c in constants_mod.CPK_COLORS.items():
                constants_mod.CPK_COLORS_PV[k] = [c.redF(), c.greenF(), c.blueF()]
        except (AttributeError, RuntimeError, TypeError, ValueError) as e:
            print(f"Failed to update CPK colors from settings: {e}")

    def open_settings_dialog(self):
        dialog = SettingsDialog(self.settings, self)
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
                self, "Reset Complete", "All settings have been reset to defaults."
            )
        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.warning(self, "Reset Failed", f"Could not reset settings: {e}")

    def _confirm_settings_reset(self):
        """Show a confirmation dialog for resetting settings."""
        return (
            QMessageBox.question(
                self,
                "Reset All Settings",
                "Are you sure you want to reset all settings to defaults?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            == QMessageBox.StandardButton.Yes
        )

    def _perform_settings_reset(self):
        """Delete the settings file and reload defaults."""
        if os.path.exists(self.settings_file):
            os.remove(self.settings_file)
        self.load_settings()
        self.settings_dirty = True

    def _refresh_ui_after_reset(self):
        """Update all UI components to reflect the reset settings."""
        # 1. Refresh related dialogs if open
        for w in QApplication.topLevelWidgets():
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                if isinstance(w, ColorSettingsDialog):
                    w.refresh_ui()
                if isinstance(w, SettingsDialog):
                    w.update_ui_from_settings(self.settings)

        # 2. Update internal state and sync CPK colors
        self.optimization_method = self.settings.get(
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
                current_method = (self.optimization_method or "").upper()
                for key, action in self.opt3d_actions.items():
                    # Use case-insensitive comparison for robustness
                    action.setChecked(key.upper() == current_method)

            # Conversion actions
            if hasattr(self, "conv_actions"):
                mode = (
                    self.settings.get("3d_conversion_mode", "fallback") or ""
                ).lower()
                for key, action in self.conv_actions.items():
                    action.setChecked(key.lower() == mode)

            # Intermolecular interaction
            if hasattr(self, "intermolecular_rdkit_action"):
                self.intermolecular_rdkit_action.setChecked(
                    self.settings.get("optimize_intermolecular_interaction_rdkit", True)
                )

    def _refresh_views_after_reset(self):
        """Refresh 2D and 3D views after settings reset."""
        # Refresh 3D View
        try:
            self.apply_3d_settings()
            if hasattr(self, "current_mol") and self.current_mol:
                self.draw_molecule_3d(self.current_mol)
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(f"Suppressed exception: {e}")

        # Refresh 2D View
        if hasattr(self, "scene") and self.scene:
            try:
                bg_c = self.settings.get("background_color_2d", "#FFFFFF")
                self.scene.setBackgroundBrush(QBrush(QColor(bg_c)))
                for item in self.scene.items():
                    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                        if hasattr(item, "update_style"):
                            item.update_style()
                self.scene.update()
                for v in self.scene.views():
                    v.viewport().update()
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                logging.debug(f"Suppressed exception: {e}")

    def load_settings(self):
        """Load settings from a JSON file, or use defaults if the file is missing."""
        # 1. Start with default settings
        self.settings = self._get_default_settings()

        # 2. Try to load from user's settings file
        try:
            if hasattr(self, "settings_file") and os.path.exists(self.settings_file):
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    loaded_settings = json.load(f)

                    # 3. Handle legacy settings migration
                    self._migrate_legacy_settings(loaded_settings)

                    # 4. Update settings with loaded values
                    self.settings.update(loaded_settings)
        except (AttributeError, RuntimeError, ValueError, TypeError, IOError):
            # Use defaults on any error
            pass

        # 5. Apply loaded settings to application state
        self.show_chiral_labels = self.settings.get("show_chiral_labels", False)
        # Apply optimization method
        if "optimization_method" in self.settings:
            self.optimization_method = self.settings["optimization_method"]

    def save_settings(self):
        try:
            if not os.path.exists(self.settings_dir):
                os.makedirs(self.settings_dir)
            with open(self.settings_file, "w", encoding="utf-8") as f:
                json.dump(self.settings, f, indent=4)
            self.settings_dirty = False
        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error saving settings: {e}")

    def update_plugin_menu(self, plugin_menu):
        """Discovers plugins and updates the plugin menu actions."""
        if not self.plugin_manager:
            return

        # 1. Cleanup
        self._clear_all_plugin_actions(plugin_menu)

        # 2. Re-add Manager
        manage_plugins_action = QAction("Plugin Manager...", self)
        manage_plugins_action.triggered.connect(
            lambda: self._show_plugin_manager(plugin_menu)
        )
        plugin_menu.addAction(manage_plugins_action)
        plugin_menu.addSeparator()

        # 3. Discover
        plugins = self.plugin_manager.discover_plugins(self)

        # 4. Integrate
        self._update_style_menu_with_plugins()
        self._add_registered_plugin_actions()
        self._add_plugin_toolbar_actions()
        self._add_legacy_plugin_actions(plugin_menu, plugins)
        self._integrate_plugin_export_actions()
        self._integrate_plugin_file_openers()
        self._integrate_plugin_analysis_tools()

    def _show_plugin_manager(self, plugin_menu):
        """Displays the plugin manager window and refreshes the menu."""
        if not self.plugin_manager:
            QMessageBox.information(
                self, "Safe Mode", "Plugins are disabled (safe mode)."
            )
            return
        from ..plugins.plugin_manager_window import PluginManagerWindow

        dlg = PluginManagerWindow(self.plugin_manager, self)
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
        for top_action in self.menuBar().actions():
            if top_action.menu():
                clear_menu(top_action.menu())

        if hasattr(self, "export_button") and self.export_button.menu():
            clear_menu(self.export_button.menu())

    def _update_style_menu_with_plugins(self):
        """Update the 3D style menu with custom styles from plugins."""
        if not hasattr(self, "style_button") or not self.style_button.menu():
            return

        style_menu = self.style_button.menu()
        style_group = next(
            (a.actionGroup() for a in style_menu.actions() if a.actionGroup()), None
        )

        if style_group and self.plugin_manager.custom_3d_styles:
            if not style_menu.actions()[-1].isSeparator():
                style_menu.addSeparator()

            for style_name in self.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self, checkable=True)
                    action.triggered.connect(
                        lambda checked=False, s=style_name: self.set_3d_style(s)
                    )
                    style_menu.addAction(action)
                    style_group.addAction(action)

    def _add_registered_plugin_actions(self):
        """Add actions that have been explicitly registered via the plugin manager."""
        PLUGIN_ACTION_TAG = "plugin_managed"
        if not self.plugin_manager.menu_actions:
            return

        for action_def in self.plugin_manager.menu_actions:
            path = action_def["path"]
            callback = action_def["callback"]
            text = action_def["text"]

            parts = path.split("/")
            top_level_title = parts[0]
            current_menu = next(
                (
                    a.menu()
                    for a in self.menuBar().actions()
                    if a.menu() and a.text().replace("&", "") == top_level_title
                ),
                None,
            )

            if not current_menu:
                current_menu = self.menuBar().addMenu(top_level_title)

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

            action = QAction(text or parts[-1], self)
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
        if self.plugin_manager.toolbar_actions:
            self.plugin_toolbar.show()
            for action_def in self.plugin_manager.toolbar_actions:
                action = QAction(action_def["text"], self)
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
            no_plugin = QAction("(No plugins found)", self)
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
                a = QAction(p["name"], self)
                a.triggered.connect(
                    lambda checked, mod=p["module"]: self.plugin_manager.run_plugin(
                        mod, self
                    )
                )
                current_parent.addAction(a)

        # Add root items
        for p in sorted(root, key=lambda x: x["name"]):
            a = QAction(p["name"], self)
            a.triggered.connect(
                lambda checked, mod=p["module"]: self.plugin_manager.run_plugin(
                    mod, self
                )
            )
            plugin_menu.addAction(a)

    def _integrate_plugin_export_actions(self):
        """Add plugin-provided export actions to the export menu."""
        if not self.plugin_manager.export_actions:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        main_export_menu = None
        for top_action in self.menuBar().actions():
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
        if hasattr(self, "export_button") and self.export_button.menu():
            targets.append(self.export_button.menu())
        if main_export_menu:
            targets.append(main_export_menu)

        for menu in targets:
            sep = menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for exp in self.plugin_manager.export_actions:
                a = QAction(exp["label"], self)
                a.triggered.connect(exp["callback"])
                a.setData(PLUGIN_ACTION_TAG)
                menu.addAction(a)

    def _integrate_plugin_file_openers(self):
        """Add plugin-provided file openers to the import menu."""
        if not hasattr(self, "import_menu") or not self.plugin_manager.file_openers:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        sep = self.import_menu.addSeparator()
        sep.setData(PLUGIN_ACTION_TAG)

        plugin_map = {}
        for ext, openers in self.plugin_manager.file_openers.items():
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
                        self, f"Import {n} Files", "", f
                    )
                    if fpath:
                        ext = os.path.splitext(fpath)[1].lower()
                        if ext in m:
                            m[ext](fpath)
                            self.current_file_path = fpath
                            self.update_window_title()

                return _cb

            a = QAction(f"Import {ext_str} ({p_name})...", self)
            a.triggered.connect(make_cb(ext_map, filter_str, p_name))
            a.setData(PLUGIN_ACTION_TAG)
            self.import_menu.addAction(a)

    def _integrate_plugin_analysis_tools(self):
        """Add plugin-provided analysis tools to the analysis menu."""
        analysis_menu = next(
            (
                a.menu()
                for a in self.menuBar().actions()
                if a.text().replace("&", "") == "Analysis"
            ),
            None,
        )
        if analysis_menu and self.plugin_manager.analysis_tools:
            PLUGIN_ACTION_TAG = "plugin_managed"
            sep = analysis_menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for tool in self.plugin_manager.analysis_tools:
                a = QAction(f"{tool['label']} ({tool.get('plugin', 'Plugin')})", self)
                a.triggered.connect(tool["callback"])
                a.setData(PLUGIN_ACTION_TAG)
                analysis_menu.addAction(a)

    # --- UI Initialization Helpers ---
    def _init_main_layout(self):
        """Initialize the main layout with splitter and panels."""
        self.splitter = QSplitter(Qt.Orientation.Horizontal)
        # Make splitter handle thicker for visibility
        self.splitter.setHandleWidth(8)
        # Improve splitter handle style
        self.splitter.setStyleSheet("""
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
        self.setCentralWidget(self.splitter)

        left_pane = QWidget()
        left_pane.setAcceptDrops(True)
        left_layout = QVBoxLayout(left_pane)
        self._init_left_panel(left_layout)
        self.splitter.addWidget(left_pane)

        right_pane = QWidget()
        right_layout = QVBoxLayout(right_pane)
        self._init_right_panel(right_layout)
        self.splitter.addWidget(right_pane)

        # Monitor splitter movement
        self.splitter.splitterMoved.connect(self.on_splitter_moved)
        self.splitter.setSizes([600, 600])

        # Set tooltip for splitter handle
        QTimer.singleShot(100, self.setup_splitter_tooltip)

        # Settings to separate status bar segments
        self.status_bar = self.statusBar()
        self.formula_label = QLabel("")  # Create label to be displayed on the right
        # Add margin to the right end for better appearance
        self.formula_label.setStyleSheet("padding-right: 8px;")
        # Add label as a permanent widget on the right
        self.status_bar.addPermanentWidget(self.formula_label)

    # --- Settings and Plugin Helpers ---
    def _get_default_settings(self):
        """Return a dictionary of default application settings."""
        return {
            "background_color": "#FFFFFF",
            "atom_label_color": "#000000",
            "bond_color": "#808080",
            "show_chiral_labels": False,
            "theme": "light",
            "window_size": [1200, 800],
            "window_position": [100, 100],
            "splitter_sizes": [600, 600],
            "last_dir": "",
            "3d_conversion_mode": "rdkit",
            "optimization_method": "MMFF_RDKIT",
            "optimize_intermolecular_interaction_rdkit": True,
            "use_high_fidelity_selection": True,
            "selection_color": "#FFD700",
            "high_quality_meshing": True,
        }

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

    def _init_left_panel(self, left_layout):
        """Initialize the left panel (2D view and buttons)."""
        self.scene = MoleculeScene(self.data, self)
        self.scene.setSceneRect(-4000, -4000, 4000, 4000)
        self.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.view_2d = ZoomableView(self.scene, self)
        self.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view_2d.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        left_layout.addWidget(self.view_2d, 1)

        self.view_2d.scale(0.75, 0.75)

        # --- Left panel button layout ---
        left_buttons_layout = QHBoxLayout()
        self.cleanup_button = QPushButton("Clean Up 2D")
        self.cleanup_button.clicked.connect(self.clean_up_2d_structure)
        left_buttons_layout.addWidget(self.cleanup_button)

        self.convert_button = QPushButton("Convert 2D to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        # Allow right-click to open a temporary conversion-mode menu
        try:
            self.convert_button.setContextMenuPolicy(
                Qt.ContextMenuPolicy.CustomContextMenu
            )
            self.convert_button.customContextMenuRequested.connect(
                self.show_convert_menu
            )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

        left_buttons_layout.addWidget(self.convert_button)
        left_layout.addLayout(left_buttons_layout)

    def _init_right_panel(self, right_layout):
        """Initialize the right panel (3D view and buttons)."""
        self.plotter = CustomQtInteractor(
            right_layout.parentWidget(), main_window=self, lighting="none"
        )
        self.plotter.setAcceptDrops(False)
        self.plotter.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.plotter.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)

        # 2. Add 3D view to layout
        right_layout.addWidget(self.plotter, 1)
        # self.plotter.installEventFilter(self)
        # 3. Create horizontal layout for buttons
        right_buttons_layout = QHBoxLayout()

        # 3D Optimize button
        self.optimize_3d_button = QPushButton("Optimize 3D")
        self.optimize_3d_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.optimize_3d_button.clicked.connect(self.optimize_3d_structure)
        self.optimize_3d_button.setEnabled(False)
        # Initialized via _enable_3d_features(False)
        # Allow right-click to open a temporary optimization-method menu
        try:
            self.optimize_3d_button.setContextMenuPolicy(
                Qt.ContextMenuPolicy.CustomContextMenu
            )
            self.optimize_3d_button.customContextMenuRequested.connect(
                self.show_optimize_menu
            )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical UI/menu initialization errors

        right_buttons_layout.addWidget(self.optimize_3d_button)

        # Export button with menu
        self.export_button = QToolButton()
        self.export_button.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.export_button.setText("Export 3D")
        self.export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.export_button.setEnabled(False)  # Initially disabled

        export_menu = QMenu(self)
        export_mol_action = QAction("Export as MOL...", self)
        export_mol_action.triggered.connect(self.save_3d_as_mol)
        export_menu.addAction(export_mol_action)

        export_xyz_action = QAction("Export as XYZ...", self)
        export_xyz_action.triggered.connect(self.save_as_xyz)
        export_menu.addAction(export_xyz_action)

        export_png_action = QAction("Export as PNG...", self)
        export_png_action.triggered.connect(self.export_3d_png)
        export_menu.addAction(export_png_action)

        self.export_button.setMenu(export_menu)
        right_buttons_layout.addWidget(self.export_button)

        # 4. Add horizontal layout to vertical layout
        right_layout.addLayout(right_buttons_layout)

    def _init_toolbars(self):
        """Initialize and organize toolbars."""
        # Row 1: Main Toolbar
        self.toolbar = QToolBar("Main Toolbar")
        self.addToolBar(self.toolbar)

        # Row 2: Templates
        with contextlib.suppress(AttributeError):
            self.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        self.toolbar_bottom = QToolBar("Templates Toolbar")
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.toolbar_bottom)

        # Row 3: Plugins
        with contextlib.suppress(AttributeError):
            self.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        self.plugin_toolbar = QToolBar("Plugin Toolbar")
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.plugin_toolbar)
        self.plugin_toolbar.hide()

    def _setup_action_groups(self, toolbar, toolbar_bottom):
        """Set up action groups and tool actions."""
        self.tool_group = QActionGroup(self)
        self.tool_group.setExclusive(True)

        self._add_atom_actions(toolbar)
        self._add_bond_actions(toolbar)
        self._add_charge_radical_actions(toolbar)
        self._add_template_actions(toolbar_bottom)
        self._add_3d_edit_actions(toolbar)

        # Set default tool
        self.set_mode("atom_C")
        if "atom_C" in self.mode_actions:
            self.mode_actions["atom_C"].setChecked(True)

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
            action = QAction(text, self, checkable=(mode != "atom_other"))
            if shortcut:
                action.setToolTip(f"{text} ({shortcut})")

            if mode == "atom_other":
                action.triggered.connect(self.open_periodic_table_dialog)
                self.other_atom_action = action
            else:
                action.triggered.connect(lambda c, m=mode: self.set_mode(m))
                self.mode_actions[mode] = action
                self.tool_group.addAction(action)
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
            action = QAction(self)
            action.setIcon(self._create_bond_icon(icon_type))
            action.setToolTip(f"{text} ({shortcut})")
            action.setCheckable(True)
            action.triggered.connect(lambda checked, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            toolbar.addAction(action)
            self.tool_group.addAction(action)
        toolbar.addSeparator()

    def _add_charge_radical_actions(self, toolbar):
        """Add charge and radical modification actions."""
        ops = [
            ("+ Charge", "charge_plus", "Increase Atom Charge (+)"),
            ("- Charge", "charge_minus", "Decrease Atom Charge (-)"),
            ("Radical", "radical", "Toggle Radical (0/1/2) (.)"),
        ]

        for text, mode, tooltip in ops:
            action = QAction(text, self, checkable=True)
            action.setToolTip(tooltip)
            action.triggered.connect(lambda c, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            toolbar.addAction(action)
            self.tool_group.addAction(action)

    def _add_template_actions(self, toolbar_bottom):
        """Add structural template actions (rings, etc.) to the bottom toolbar."""
        toolbar_bottom.addWidget(QLabel(" Templates:"))
        templates = [("Benzene", "template_benzene", 6)] + [
            (f"{i}-Ring", f"template_{i}", i) for i in range(3, 10)
        ]

        for text, mode, n in templates:
            action = QAction(self)
            action.setCheckable(True)
            is_benzene = text == "Benzene"
            action.setIcon(self._create_template_icon(n, is_benzene=is_benzene))
            action.setToolTip(
                f"{text} Template (4)" if is_benzene else f"{text} Template"
            )
            action.triggered.connect(lambda c, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            toolbar_bottom.addAction(action)
            self.tool_group.addAction(action)

        user_action = QAction("USER", self, checkable=True)
        user_action.setToolTip("Open User Templates Dialog")
        user_action.triggered.connect(self.open_template_dialog_and_activate)
        self.mode_actions["template_user"] = user_action
        toolbar_bottom.addAction(user_action)
        self.tool_group.addAction(user_action)

    def _add_3d_edit_actions(self, toolbar):
        """Add 3D-specific selection and manipulation actions."""
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        toolbar.addWidget(spacer)

        self.measurement_action = QAction("3D Select", self, checkable=True)
        self.measurement_action.setToolTip(
            "Enable distance, angle, and dihedral measurement in 3D view"
        )
        self.measurement_action.triggered.connect(self.toggle_measurement_mode)
        toolbar.addAction(self.measurement_action)

        self.edit_3d_action = QAction("3D Drag", self, checkable=True)
        self.edit_3d_action.setToolTip(
            "Toggle 3D atom dragging mode (Hold Alt for temporary mode)"
        )
        self.edit_3d_action.toggled.connect(self.toggle_3d_edit_mode)
        toolbar.addAction(self.edit_3d_action)

        self.style_button = QToolButton()
        self.style_button.setText("3D Style")
        self.style_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(self.style_button)

        style_menu = QMenu(self)
        self.style_button.setMenu(style_menu)
        style_group = QActionGroup(self)
        style_group.setExclusive(True)

        for name, key in [
            ("Ball & Stick", "ball_and_stick"),
            ("CPK (Space-filling)", "cpk"),
            ("Wireframe", "wireframe"),
            ("Stick", "stick"),
        ]:
            action = QAction(name, self, checkable=True)
            if key == "ball_and_stick":
                action.setChecked(True)
            action.triggered.connect(
                lambda checked=False, k=key: (
                    self.set_3d_style(k),
                    self.draw_molecule_3d(self.current_mol)
                    if getattr(self, "current_mol", None)
                    else None,
                )
            )
            style_menu.addAction(action)
            style_group.addAction(action)

        if self.plugin_manager and self.plugin_manager.custom_3d_styles:
            style_menu.addSeparator()
            for style_name in self.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self, checkable=True)
                    action.triggered.connect(
                        lambda checked=False, s=style_name: (
                            self.set_3d_style(s),
                            self.draw_molecule_3d(self.current_mol)
                            if getattr(self, "current_mol", None)
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
            fg = self.settings.get("icon_foreground")
            if fg and QColor(fg).isValid():
                return QColor(fg)

        with contextlib.suppress(Exception):
            os_pref = detect_system_dark_mode()
            if os_pref is not None:
                return QColor("#FFFFFF") if os_pref else QColor("#000000")

        with contextlib.suppress(Exception):
            bg = QColor(self.settings.get("background_color", "#FFFFFF"))
            if bg.isValid():
                lum = 0.2126 * bg.redF() + 0.7152 * bg.greenF() + 0.0722 * bg.blueF()
                return QColor("#FFFFFF") if lum < 0.5 else QColor("#000000")
        return QColor("#000000")

    # --- Menu Bar Initialization Helpers ---
    def _init_file_menu(self, menu_bar):
        """Initialize the File menu."""
        file_menu = menu_bar.addMenu("&File")

        # === Project Operations ===
        new_action = QAction("&New", self)
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.clear_all)
        file_menu.addAction(new_action)

        load_project_action = QAction("&Open Project...", self)
        load_project_action.setShortcut("Ctrl+O")
        load_project_action.triggered.connect(self.open_project_file)
        file_menu.addAction(load_project_action)

        save_action = QAction("&Save Project", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_project)
        file_menu.addAction(save_action)

        save_as_action = QAction("Save Project &As...", self)
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.save_project_as)
        file_menu.addAction(save_as_action)

        save_template_action = QAction("Save 2D as Template...", self)
        save_template_action.triggered.connect(self.save_2d_as_template)
        file_menu.addAction(save_template_action)

        file_menu.addSeparator()

        # === Import ===
        self.import_menu = file_menu.addMenu("Import")
        load_mol_action = QAction("MOL/SDF File...", self)
        load_mol_action.triggered.connect(self.load_mol_file)
        self.import_menu.addAction(load_mol_action)

        import_smiles_action = QAction("SMILES...", self)
        import_smiles_action.triggered.connect(self.import_smiles_dialog)
        self.import_menu.addAction(import_smiles_action)

        import_inchi_action = QAction("InChI...", self)
        import_inchi_action.triggered.connect(self.import_inchi_dialog)
        self.import_menu.addAction(import_inchi_action)

        self.import_menu.addSeparator()
        load_3d_mol_action = QAction("3D MOL/SDF (3D View Only)...", self)
        load_3d_mol_action.triggered.connect(self.load_mol_file_for_3d_viewing)
        self.import_menu.addAction(load_3d_mol_action)

        load_3d_xyz_action = QAction("3D XYZ (3D View Only)...", self)
        load_3d_xyz_action.triggered.connect(self.load_xyz_for_3d_viewing)
        self.import_menu.addAction(load_3d_xyz_action)

        # === Export ===
        export_menu = file_menu.addMenu("Export")
        export_pmeraw_action = QAction("PME Raw Format...", self)
        export_pmeraw_action.triggered.connect(self.save_raw_data)
        export_menu.addAction(export_pmeraw_action)
        export_menu.addSeparator()

        export_2d_menu = export_menu.addMenu("2D Formats")
        save_mol_action = QAction("MOL File...", self)
        save_mol_action.triggered.connect(self.save_as_mol)
        export_2d_menu.addAction(save_mol_action)

        export_2d_png_action = QAction("PNG Image...", self)
        export_2d_png_action.triggered.connect(self.export_2d_png)
        export_2d_menu.addAction(export_2d_png_action)

        export_2d_svg_action = QAction("SVG Image...", self)
        export_2d_svg_action.triggered.connect(self.export_2d_svg)
        export_2d_menu.addAction(export_2d_svg_action)

        export_3d_menu = export_menu.addMenu("3D Formats")
        save_3d_mol_action = QAction("MOL File...", self)
        save_3d_mol_action.triggered.connect(self.save_3d_as_mol)
        export_3d_menu.addAction(save_3d_mol_action)

        save_xyz_action = QAction("XYZ File...", self)
        save_xyz_action.triggered.connect(self.save_as_xyz)
        export_3d_menu.addAction(save_xyz_action)

        export_3d_png_action = QAction("PNG Image...", self)
        export_3d_png_action.triggered.connect(self.export_3d_png)
        export_3d_menu.addAction(export_3d_png_action)
        export_3d_menu.addSeparator()

        export_stl_action = QAction("STL File...", self)
        export_stl_action.triggered.connect(self.export_stl)
        export_3d_menu.addAction(export_stl_action)

        export_obj_action = QAction("OBJ/MTL (with colors)...", self)
        export_obj_action.triggered.connect(self.export_obj_mtl)
        export_3d_menu.addAction(export_obj_action)

        file_menu.addSeparator()
        quit_action = QAction("Quit", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    def _init_edit_menu(self, menu_bar):
        """Initialize the Edit menu."""
        edit_menu = menu_bar.addMenu("&Edit")
        self.undo_action = QAction("Undo", self)
        self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.triggered.connect(self.undo)
        edit_menu.addAction(self.undo_action)

        self.redo_action = QAction("Redo", self)
        self.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.redo_action.triggered.connect(self.redo)
        edit_menu.addAction(self.redo_action)

        edit_menu.addSeparator()
        self.cut_action = QAction("Cut", self)
        self.cut_action.setShortcut(QKeySequence.StandardKey.Cut)
        self.cut_action.triggered.connect(self.cut_selection)
        edit_menu.addAction(self.cut_action)

        self.copy_action = QAction("Copy", self)
        self.copy_action.setShortcut(QKeySequence.StandardKey.Copy)
        self.copy_action.triggered.connect(self.copy_selection)
        edit_menu.addAction(self.copy_action)

        self.paste_action = QAction("Paste", self)
        self.paste_action.setShortcut(QKeySequence.StandardKey.Paste)
        self.paste_action.triggered.connect(self.paste_from_clipboard)
        edit_menu.addAction(self.paste_action)

        edit_menu.addSeparator()
        add_hydrogen_action = QAction("Add Hydrogens", self)
        add_hydrogen_action.setToolTip(
            "Add explicit hydrogens based on RDKit implicit counts"
        )
        add_hydrogen_action.triggered.connect(self.add_hydrogen_atoms)
        edit_menu.addAction(add_hydrogen_action)

        remove_hydrogen_action = QAction("Remove Hydrogens", self)
        remove_hydrogen_action.triggered.connect(self.remove_hydrogen_atoms)
        edit_menu.addAction(remove_hydrogen_action)

        edit_menu.addSeparator()
        rotate_2d_action = QAction("Rotate 2D...", self)
        rotate_2d_action.setShortcut(QKeySequence("Ctrl+R"))
        rotate_2d_action.triggered.connect(self.open_rotate_2d_dialog)
        edit_menu.addAction(rotate_2d_action)

        edit_menu.addSeparator()
        optimize_2d_action = QAction("Clean Up 2D", self)
        optimize_2d_action.setShortcut(QKeySequence("Ctrl+J"))
        optimize_2d_action.triggered.connect(self.clean_up_2d_structure)
        edit_menu.addAction(optimize_2d_action)

        convert_3d_action = QAction("Convert 2D to 3D", self)
        convert_3d_action.setShortcut(QKeySequence("Ctrl+K"))
        convert_3d_action.triggered.connect(self.trigger_conversion)
        edit_menu.addAction(convert_3d_action)

        optimize_3d_action = QAction("Optimize 3D", self)
        optimize_3d_action.setShortcut(QKeySequence("Ctrl+L"))
        optimize_3d_action.triggered.connect(self.optimize_3d_structure)
        edit_menu.addAction(optimize_3d_action)

        edit_menu.addSeparator()
        select_all_action = QAction("Select All", self)
        select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.select_all)
        edit_menu.addAction(select_all_action)

        clear_all_action = QAction("Clear All", self)
        clear_all_action.setShortcut(QKeySequence("Ctrl+Shift+C"))
        clear_all_action.triggered.connect(self.clear_all)
        edit_menu.addAction(clear_all_action)

    def _init_view_menu(self, menu_bar):
        """Initialize the View menu."""
        view_menu = menu_bar.addMenu("&View")
        zoom_in_action = QAction("Zoom In", self)
        zoom_in_action.setShortcut(QKeySequence.StandardKey.ZoomIn)
        zoom_in_action.triggered.connect(self.zoom_in)
        view_menu.addAction(zoom_in_action)

        zoom_out_action = QAction("Zoom Out", self)
        zoom_out_action.setShortcut(QKeySequence.StandardKey.ZoomOut)
        zoom_out_action.triggered.connect(self.zoom_out)
        view_menu.addAction(zoom_out_action)

        reset_zoom_action = QAction("Reset Zoom", self)
        reset_zoom_action.setShortcut(QKeySequence("Ctrl+0"))
        reset_zoom_action.triggered.connect(self.reset_zoom)
        view_menu.addAction(reset_zoom_action)

        fit_action = QAction("Fit to View", self)
        fit_action.setShortcut(QKeySequence("Ctrl+9"))
        fit_action.triggered.connect(self.fit_to_view)
        view_menu.addAction(fit_action)

        view_menu.addSeparator()
        reset_3d_view_action = QAction("Reset 3D View", self)
        reset_3d_view_action.triggered.connect(
            lambda: self.plotter.reset_camera() if hasattr(self, "plotter") else None
        )
        reset_3d_view_action.setShortcut(QKeySequence("Ctrl+Shift+R"))
        view_menu.addAction(reset_3d_view_action)

        self.redraw_menu_action = QAction("Redraw 3D Molecule", self)
        self.redraw_menu_action.triggered.connect(self.redraw_molecule_3d)
        view_menu.addAction(self.redraw_menu_action)

        view_menu.addSeparator()
        layout_menu = view_menu.addMenu("Panel Layout")
        equal_panels_action = QAction("Equal Panels (50:50)", self)
        equal_panels_action.setShortcut(QKeySequence("Ctrl+1"))
        equal_panels_action.triggered.connect(lambda: self.set_panel_layout(50, 50))
        layout_menu.addAction(equal_panels_action)

        layout_2d_focus_action = QAction("2D Focus (70:30)", self)
        layout_2d_focus_action.setShortcut(QKeySequence("Ctrl+2"))
        layout_2d_focus_action.triggered.connect(lambda: self.set_panel_layout(70, 30))
        layout_menu.addAction(layout_2d_focus_action)

        layout_3d_focus_action = QAction("3D Focus (30:70)", self)
        layout_3d_focus_action.setShortcut(QKeySequence("Ctrl+3"))
        layout_3d_focus_action.triggered.connect(lambda: self.set_panel_layout(30, 70))
        layout_menu.addAction(layout_3d_focus_action)

        layout_menu.addSeparator()
        toggle_2d_panel_action = QAction("Toggle 2D Panel", self)
        toggle_2d_panel_action.setShortcut(QKeySequence("Ctrl+H"))
        toggle_2d_panel_action.triggered.connect(self.toggle_2d_panel)
        layout_menu.addAction(toggle_2d_panel_action)

        view_menu.addSeparator()
        self.toggle_chiral_action = QAction("Show Chiral Labels", self, checkable=True)
        self.toggle_chiral_action.setChecked(self.show_chiral_labels)
        self.toggle_chiral_action.triggered.connect(self.toggle_chiral_labels_display)
        view_menu.addAction(self.toggle_chiral_action)

        view_menu.addSeparator()
        atom_info_menu = view_menu.addMenu("3D Atom Info Display")
        self.show_atom_id_action = QAction(
            "Show Original ID / Index", self, checkable=True
        )
        self.show_atom_id_action.triggered.connect(
            lambda: self.toggle_atom_info_display("id")
        )
        atom_info_menu.addAction(self.show_atom_id_action)

        self.show_rdkit_id_action = QAction("Show RDKit Index", self, checkable=True)
        self.show_rdkit_id_action.triggered.connect(
            lambda: self.toggle_atom_info_display("rdkit_id")
        )
        atom_info_menu.addAction(self.show_rdkit_id_action)

        self.show_atom_coords_action = QAction(
            "Show Coordinates (X,Y,Z)", self, checkable=True
        )
        self.show_atom_coords_action.triggered.connect(
            lambda: self.toggle_atom_info_display("coords")
        )
        atom_info_menu.addAction(self.show_atom_coords_action)

        self.show_atom_symbol_action = QAction(
            "Show Element Symbol", self, checkable=True
        )
        self.show_atom_symbol_action.triggered.connect(
            lambda: self.toggle_atom_info_display("symbol")
        )
        atom_info_menu.addAction(self.show_atom_symbol_action)

    def _init_analysis_menu(self, menu_bar):
        """Initialize the Analysis menu."""
        analysis_menu = menu_bar.addMenu("&Analysis")
        self.analysis_action = QAction("Show Analysis...", self)
        self.analysis_action.triggered.connect(self.open_analysis_window)
        self.analysis_action.setEnabled(False)
        analysis_menu.addAction(self.analysis_action)

    def _init_edit_3d_menu(self, menu_bar):
        """Initialize the 3D Edit menu."""
        edit_3d_menu = menu_bar.addMenu("3D &Edit")
        translation_action = QAction("Translation...", self)
        translation_action.triggered.connect(self.open_translation_dialog)
        translation_action.setEnabled(False)
        edit_3d_menu.addAction(translation_action)
        self.translation_action = translation_action

        move_group_action = QAction("Move Group...", self)
        move_group_action.triggered.connect(self.open_move_group_dialog)
        move_group_action.setEnabled(False)
        edit_3d_menu.addAction(move_group_action)
        self.move_group_action = move_group_action

        edit_3d_menu.addSeparator()
        align_menu = edit_3d_menu.addMenu("Align to")
        align_menu.setEnabled(False)
        self.align_menu = align_menu

        axis_align_menu = align_menu.addMenu("Axis")
        align_x_action = QAction("X-axis", self)
        align_x_action.triggered.connect(lambda: self.open_alignment_dialog("x"))
        align_x_action.setEnabled(False)
        axis_align_menu.addAction(align_x_action)
        self.align_x_action = align_x_action

        align_y_action = QAction("Y-axis", self)
        align_y_action.triggered.connect(lambda: self.open_alignment_dialog("y"))
        align_y_action.setEnabled(False)
        axis_align_menu.addAction(align_y_action)
        self.align_y_action = align_y_action

        align_z_action = QAction("Z-axis", self)
        align_z_action.triggered.connect(lambda: self.open_alignment_dialog("z"))
        align_z_action.setEnabled(False)
        axis_align_menu.addAction(align_z_action)
        self.align_z_action = align_z_action

        plane_align_menu = align_menu.addMenu("Plane")
        alignplane_xy_action = QAction("XY-plane", self)
        alignplane_xy_action.triggered.connect(
            lambda: self.open_align_plane_dialog("xy")
        )
        alignplane_xy_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xy_action)
        self.alignplane_xy_action = alignplane_xy_action

        alignplane_xz_action = QAction("XZ-plane", self)
        alignplane_xz_action.triggered.connect(
            lambda: self.open_align_plane_dialog("xz")
        )
        alignplane_xz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xz_action)
        self.alignplane_xz_action = alignplane_xz_action

        alignplane_yz_action = QAction("YZ-plane", self)
        alignplane_yz_action.triggered.connect(
            lambda: self.open_align_plane_dialog("yz")
        )
        alignplane_yz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_yz_action)
        self.alignplane_yz_action = alignplane_yz_action

        edit_3d_menu.addSeparator()
        mirror_action = QAction("Mirror...", self)
        mirror_action.triggered.connect(self.open_mirror_dialog)
        mirror_action.setEnabled(False)
        edit_3d_menu.addAction(mirror_action)
        self.mirror_action = mirror_action

        edit_3d_menu.addSeparator()
        planarize_action = QAction("Planarize...", self)
        planarize_action.triggered.connect(lambda: self.open_planarize_dialog(None))
        planarize_action.setEnabled(False)
        edit_3d_menu.addAction(planarize_action)
        self.planarize_action = planarize_action

        edit_3d_menu.addSeparator()
        bond_length_action = QAction("Adjust Bond Length...", self)
        bond_length_action.triggered.connect(self.open_bond_length_dialog)
        bond_length_action.setEnabled(False)
        edit_3d_menu.addAction(bond_length_action)
        self.bond_length_action = bond_length_action

        angle_action = QAction("Adjust Angle...", self)
        angle_action.triggered.connect(self.open_angle_dialog)
        angle_action.setEnabled(False)
        edit_3d_menu.addAction(angle_action)
        self.angle_action = angle_action

        dihedral_action = QAction("Adjust Dihedral Angle...", self)
        dihedral_action.triggered.connect(self.open_dihedral_dialog)
        dihedral_action.setEnabled(False)
        edit_3d_menu.addAction(dihedral_action)
        self.dihedral_action = dihedral_action

        edit_3d_menu.addSeparator()
        constrained_opt_action = QAction("Constrained Optimization...", self)
        constrained_opt_action.triggered.connect(
            self.open_constrained_optimization_dialog
        )
        constrained_opt_action.setEnabled(False)
        edit_3d_menu.addAction(constrained_opt_action)
        self.constrained_opt_action = constrained_opt_action

    def _init_plugin_menu(self, menu_bar):
        """Initialize the Plugin menu."""
        self.plugin_menu = menu_bar.addMenu("&Plugin")
        manage_plugins_action = QAction("Plugin Manager...", self)

        def show_plugin_manager():
            if not self.plugin_manager:
                QMessageBox.information(
                    self, "Safe Mode", "Plugins are disabled (safe mode)."
                )
                return
            from ..plugins.plugin_manager_window import PluginManagerWindow

            dlg = PluginManagerWindow(self.plugin_manager, self)
            dlg.exec()
            self.update_plugin_menu(self.plugin_menu)

        manage_plugins_action.triggered.connect(show_plugin_manager)
        self.plugin_menu.addAction(manage_plugins_action)
        self.plugin_menu.addSeparator()

    def _init_settings_menu(self, menu_bar):
        """Initialize the Settings menu."""
        settings_menu = menu_bar.addMenu("&Settings")
        view_settings_action = QAction("Settings...", self)
        view_settings_action.triggered.connect(self.open_settings_dialog)
        settings_menu.addAction(view_settings_action)

        color_action = QAction("CPK Colors...", self)
        color_action.triggered.connect(
            lambda: ColorSettingsDialog(self.settings, parent=self).exec_()
        )
        settings_menu.addAction(color_action)

        conversion_menu = settings_menu.addMenu("3D Conversion")
        conv_group = QActionGroup(self)
        conv_group.setExclusive(True)

        def _set_conv_mode(mode):
            try:
                self.settings["3d_conversion_mode"] = mode
                self.settings_dirty = True
                self.statusBar().showMessage(f"3D conversion mode set to: {mode}")
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
            a = QAction(label, self)
            a.setCheckable(True)
            if key == "obabel" and not OBABEL_AVAILABLE:
                a.setEnabled(False)
            a.triggered.connect(lambda checked, m=key: _set_conv_mode(m))
            conversion_menu.addAction(a)
            conv_group.addAction(a)
            self.conv_actions[key] = a

        saved_conv = self.settings.get("3d_conversion_mode", "fallback")
        if (
            saved_conv not in self.conv_actions
            or not self.conv_actions[saved_conv].isEnabled()
        ):
            saved_conv = (
                "rdkit"
                if self.conv_actions.get("rdkit", QAction(self)).isEnabled()
                else "fallback"
            )

        if saved_conv in self.conv_actions:
            self.conv_actions[saved_conv].setChecked(True)
        self.settings["3d_conversion_mode"] = saved_conv

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
        opt_group = QActionGroup(self)
        opt_group.setExclusive(True)
        self.opt3d_actions = {}
        for label, key in opt_methods:
            action = QAction(label, self)
            action.setCheckable(True)
            if key.endswith("_OBABEL") and not OBABEL_AVAILABLE:
                action.setEnabled(False)
            action.triggered.connect(
                lambda checked, m=key: self.set_optimization_method(m)
            )
            optimization_menu.addAction(action)
            opt_group.addAction(action)
            self.opt3d_actions[key] = action

        optimization_menu.addSeparator()
        self.intermolecular_rdkit_action = QAction(
            "Consider Intermolecular Interaction for RDKit", self
        )
        self.intermolecular_rdkit_action.setCheckable(True)
        self.intermolecular_rdkit_action.setChecked(
            self.settings.get("optimize_intermolecular_interaction_rdkit", True)
        )
        self.intermolecular_rdkit_action.triggered.connect(
            self.toggle_intermolecular_interaction_rdkit
        )
        optimization_menu.addAction(self.intermolecular_rdkit_action)

        saved_opt = (self.settings.get("optimization_method") or "MMFF_RDKIT").upper()
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
        reset_settings_action = QAction("Reset All Settings", self)
        reset_settings_action.triggered.connect(self.reset_all_settings_menu)
        settings_menu.addAction(reset_settings_action)

    def _init_help_menu(self, menu_bar):
        """Initialize the Help menu."""
        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(about_action)

        github_action = QAction("GitHub", self)
        github_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl("https://github.com/HiroYokoyama/python_molecular_editor")
            )
        )
        help_menu.addAction(github_action)

        github_wiki_action = QAction("GitHub Wiki", self)
        github_wiki_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl("https://github.com/HiroYokoyama/python_molecular_editor/wiki")
            )
        )
        help_menu.addAction(github_wiki_action)

        manual_action = QAction("User Manual", self)
        manual_action.triggered.connect(
            lambda: QDesktopServices.openUrl(
                QUrl(
                    "https://hiroyokoyama.github.io/python_molecular_editor/manual/manual"
                )
            )
        )
        help_menu.addAction(manual_action)
