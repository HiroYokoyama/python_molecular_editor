#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Any, Callable, List, Optional, Union


class PluginContext:
    """
    PluginContext provides a safe interface for plugins to interact with the application.
    It is passed to the `initialize(context)` function of the plugin.
    """

    def __init__(self, manager, plugin_name: str):
        self._manager = manager
        self._plugin_name = plugin_name

    def add_menu_action(
        self,
        path: str,
        callback: Callable,
        text: Optional[str] = None,
        icon: Optional[str] = None,
        shortcut: Optional[str] = None,
    ):
        """
        Register a menu action.

        Args:
            path: Menu path, e.g., "File/Import", "Edit", or "MyPlugin" (top level).
            callback: Function to call when triggered.
            text: Label for the action (defaults to last part of path if None).
            icon: Path to icon or icon name (optional).
            shortcut: Keyboard shortcut (optional).
        """
        self._manager.register_menu_action(
            self._plugin_name, path, callback, text, icon, shortcut
        )

    def register_menu_action(
        self,
        path: str,
        text_or_callback: Union[str, Callable],
        callback: Optional[Callable] = None,
        icon: Optional[str] = None,
        shortcut: Optional[str] = None,
    ):
        """Backward-compatible alias for add_menu_action.
        Supports old 3-arg style: register_menu_action(path, text, callback).
        """
        if callable(text_or_callback):
            # New style: (path, callback, ...)
            self.add_menu_action(path, text_or_callback, None, icon, shortcut)
        else:
            # Old style: (path, text, callback)
            self.add_menu_action(path, callback, text_or_callback, icon, shortcut)

    def add_plugin_menu(
        self,
        path: str,
        callback: Callable,
        text: Optional[str] = None,
        icon: Optional[str] = None,
        shortcut: Optional[str] = None,
    ):
        """
        Register an action nested inside the Plugin menu.

        Equivalent to add_menu_action("Plugin/<path>", ...).
        Use this instead of add_menu_action when you want your plugin to appear
        as a nested folder inside the Plugin menu rather than as a top-level menu.

        Args:
            path: Sub-path within the Plugin menu, e.g. "Utility/My Tool"
                  results in Plugin > Utility > My Tool.
            callback: Function to call when triggered.
            text: Label override (defaults to last part of path).
            icon: Path to icon (optional).
            shortcut: Keyboard shortcut (optional).
        """
        full_path = f"Plugin/{path.lstrip('/')}"
        self.add_menu_action(full_path, callback, text, icon, shortcut)

    def add_toolbar_action(
        self,
        callback: Callable,
        text: str,
        icon: Optional[str] = None,
        tooltip: Optional[str] = None,
    ):
        """
        Register a toolbar action.
        """
        self._manager.register_toolbar_action(
            self._plugin_name, callback, text, icon, tooltip
        )

    def register_drop_handler(self, callback: Callable[[str], bool], priority: int = 0):
        """
        Register a handler for file drops.

        Args:
            callback: Function taking (file_path) -> bool. Returns True if handled.
            priority: Higher priority handlers are tried first.
        """
        self._manager.register_drop_handler(self._plugin_name, callback, priority)

    def get_3d_controller(self) -> "Plugin3DController":
        """
        Returns a controller to manipulate the 3D scene (e.g. colors).
        """
        return Plugin3DController(self._manager.get_main_window())

    def show_status_message(self, message: str, timeout: int = 3000) -> None:
        """
        Display a message in the application status bar.
        """
        self._manager.show_status_message(message, timeout)

    def push_undo_checkpoint(self) -> None:
        """
        Create an undo checkpoint for the current state.
        Call this AFTER making modifications to the molecule to ensure the
        new state is saved to the undo history.
        """
        self._manager.push_undo_checkpoint()

    def refresh_3d_view(self) -> None:
        """
        Force a refresh (re-render) of the 3D scene.
        """
        self._manager.refresh_3d_view()

    def reset_3d_camera(self) -> None:
        """
        Resets the 3D camera to fit the current molecule.
        """
        self._manager.reset_3d_camera()

    def get_selected_atom_indices(self) -> List[int]:
        """
        Returns a list of RDKit atom indices currently selected in the 2D or 3D view.
        Note: RDKit indices are returned, which map to the current_mol.
        """
        return self._manager.get_selected_atom_indices()

    def register_window(self, window_id: str, window: Any) -> None:
        """
        Register a custom plugin window/dialog with the application.
        This allows the application to manage the window lifecycle.
        Windows are namespaced by the plugin name automatically.
        """
        self._manager.register_window(self._plugin_name, window_id, window)

    def get_window(self, window_id: str) -> Optional[Any]:
        """
        Retrieve a previously registered window by its ID.
        """
        return self._manager.get_window(self._plugin_name, window_id)

    def get_main_window(self) -> Any:
        """
        Returns the raw MainWindow instance.
        Use with caution; prefer specific methods when available.
        """
        return self._manager.get_main_window()

    @property
    def current_mol(self) -> Any:
        """
        Get or set the current molecule (RDKit Mol object). Shortcut for current_molecule.
        """
        mw = self.get_main_window()
        return mw.current_mol if mw else None

    @current_mol.setter
    def current_mol(self, mol: Any):
        mw = self.get_main_window()
        if mw:
            mw.current_mol = mol
            if hasattr(mw.view_3d_manager, "draw_molecule_3d"):
                mw.view_3d_manager.draw_molecule_3d(mol)

    @property
    def current_molecule(self) -> Any:
        """Alias for current_mol for backward compatibility."""
        return self.current_mol

    @current_molecule.setter
    def current_molecule(self, mol: Any):
        self.current_mol = mol

    @property
    def plotter(self) -> Any:
        """
        Returns the PyVista plotter from the MainWindow.
        """
        mw = self.get_main_window()
        return mw.plotter if mw else None

    @property
    def scene(self) -> Any:
        """
        Returns the 2D MoleculeScene from the MainWindow.
        """
        mw = self.get_main_window()
        return mw.scene if mw else None

    def add_export_action(self, label: str, callback: Callable):
        """
        Register a custom export action.

        Args:
            label: Text to display in the Export menu (e.g., "Export as MyFormat...").
            callback: Function to call when triggered.
        """
        self._manager.register_export_action(self._plugin_name, label, callback)

    def register_optimization_method(
        self, method_name: str, callback: Callable[[Any], bool]
    ):
        """
        Register a custom 3D optimization method.

        Args:
            method_name: Name of the method to display in 3D Optimization menu.
            callback: Function taking (rdkit_mol) -> bool (success).
                      Modifies the molecule in-place.
        """
        self._manager.register_optimization_method(
            self._plugin_name, method_name, callback
        )

    def register_file_opener(
        self, extension: str, callback: Callable[[str], None], priority: int = 0
    ):
        """
        Register a handler for opening a specific file extension.

        Args:
            extension: File extension including dot, e.g. ".xyz".
            callback: Function taking (file_path) -> None.
                      Should load the file into the main window.
            priority: Higher priority handlers are tried first (default 0).
        """
        self._manager.register_file_opener(
            self._plugin_name, extension, callback, priority
        )

    def add_analysis_tool(self, label: str, callback: Callable):
        """
        Register a tool in the Analysis menu.

        Args:
            label: Text to display in the menu.
            callback: Function to contact when triggered.
        """
        self._manager.register_analysis_tool(self._plugin_name, label, callback)

    def register_save_handler(self, callback: Callable[[], dict]):
        """
        Register a callback to save state into the project file.

        Args:
            callback: Function returning a dict of serializable data.
        """
        self._manager.register_save_handler(self._plugin_name, callback)

    def register_load_handler(self, callback: Callable[[dict], None]):
        """
        Register a callback to restore state from the project file.

        Args:
            callback: Function receiving the dict of saved data.
        """
        self._manager.register_load_handler(self._plugin_name, callback)

    def register_3d_context_menu(self, callback: Callable, label: str):
        """Deprecated: This method does nothing. Kept for backward compatibility."""
        print(
            f"Warning: Plugin '{self._plugin_name}' uses deprecated 'register_3d_context_menu'. This API has been removed."
        )

    def register_3d_style(self, style_name: str, callback: Callable[[Any, Any], None]):
        """
        Register a custom 3D rendering style.

        Args:
            style_name: Name of the style (must be unique).
            callback: Function taking (main_window, mol) -> None.
                      Should fully handle drawing the molecule in the 3D view.
        """
        self._manager.register_3d_style(self._plugin_name, style_name, callback)

    def register_document_reset_handler(self, callback: Callable[[], None]):
        """
        Register a callback to be called when a new document is created (File→New).

        Args:
            callback: Function with no arguments that resets plugin state.
        """
        self._manager.register_document_reset_handler(self._plugin_name, callback)

    def get_setting(self, key: str, default: Any = None) -> Any:
        """
        Get a plugin-specific persistent setting.

        Settings are stored in the app settings dict under 'plugin.<plugin_name>.<key>'
        and are persisted across sessions with the application settings.

        Args:
            key: Setting key name.
            default: Value to return if the setting is not found.
        """
        mw = self.get_main_window()
        if mw and hasattr(mw, "init_manager") and hasattr(mw.init_manager, "settings"):
            namespaced = f"plugin.{self._plugin_name}.{key}"
            return mw.init_manager.settings.get(namespaced, default)
        return default

    def set_setting(self, key: str, value: Any) -> None:
        """
        Save a plugin-specific persistent setting.

        Settings are stored in the app settings dict under 'plugin.<plugin_name>.<key>'
        and are saved when the application saves its settings.

        Args:
            key: Setting key name.
            value: Value to store (must be JSON-serializable).
        """
        mw = self.get_main_window()
        if mw and hasattr(mw, "init_manager") and hasattr(mw.init_manager, "settings"):
            namespaced = f"plugin.{self._plugin_name}.{key}"
            mw.init_manager.settings[namespaced] = value
            if hasattr(mw.init_manager, "settings_dirty"):
                mw.init_manager.settings_dirty = True


class Plugin3DController:
    """Helper to manipulate the 3D scene."""

    def __init__(self, main_window):
        self._mw = main_window

    def _get_v3d(self):
        """Helper to get the 3D manager."""
        return getattr(self._mw, "view_3d_manager", None)

    def set_atom_color(self, atom_index: int, color_hex: str):
        """
        Set the color of a specific atom in the 3D view.
        Args:
            atom_index: RDKit atom index.
            color_hex: Hex string e.g., "#FF0000".
        """
        v3d = self._get_v3d()
        if v3d:
            v3d.update_atom_color_override(atom_index, color_hex)
            if hasattr(self._mw, "plotter") and self._mw.plotter:
                self._mw.plotter.render()

    def set_bond_color(self, bond_index: int, color_hex: str):
        """
        Set the color of a specific bond in the 3D view.

        Args:
             bond_index: RDKit bond index.
             color_hex: Hex string e.g., "#00FF00".
        """
        v3d = self._get_v3d()
        if v3d:
            v3d.update_bond_color_override(bond_index, color_hex)
            if hasattr(self._mw, "plotter") and self._mw.plotter:
                self._mw.plotter.render()

    def set_bond_color_by_atoms(self, atom_idx1: int, atom_idx2: int, color_hex: str):
        """
        Set the color of the bond between two atoms.

        Args:
            atom_idx1: First RDKit atom index.
            atom_idx2: Second RDKit atom index.
            color_hex: Hex string e.g., "#00FF00".
        """
        mol = getattr(self._mw, "current_mol", None)
        if not mol:
            return

        bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
        if bond:
            self.set_bond_color(bond.GetIdx(), color_hex)
