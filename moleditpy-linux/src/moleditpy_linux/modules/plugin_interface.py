#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Callable, Optional, Any

class PluginContext:
    """
    PluginContext provides a safe interface for plugins to interact with the application.
    It is passed to the `initialize(context)` function of the plugin.
    """
    def __init__(self, manager, plugin_name: str):
        self._manager = manager
        self._plugin_name = plugin_name

    def add_menu_action(self, path: str, callback: Callable, text: Optional[str] = None, icon: Optional[str] = None, shortcut: Optional[str] = None):
        """
        Register a menu action.
        
        Args:
            path: Menu path, e.g., "File/Import", "Edit", or "MyPlugin" (top level).
            callback: Function to call when triggered.
            text: Label for the action (defaults to last part of path if None).
            icon: Path to icon or icon name (optional).
            shortcut: Keyboard shortcut (optional).
        """
        self._manager.register_menu_action(self._plugin_name, path, callback, text, icon, shortcut)

    def add_toolbar_action(self, callback: Callable, text: str, icon: Optional[str] = None, tooltip: Optional[str] = None):
        """
        Register a toolbar action.
        """
        self._manager.register_toolbar_action(self._plugin_name, callback, text, icon, tooltip)

    def register_drop_handler(self, callback: Callable[[str], bool], priority: int = 0):
        """
        Register a handler for file drops.
        
        Args:
            callback: Function taking (file_path) -> bool. Returns True if handled.
            priority: Higher priority handlers are tried first.
        """
        self._manager.register_drop_handler(self._plugin_name, callback, priority)



    def get_3d_controller(self) -> 'Plugin3DController':
        """
        Returns a controller to manipulate the 3D scene (e.g. colors).
        """
        return Plugin3DController(self._manager.get_main_window())

    def get_main_window(self) -> Any:
        """
        Returns the raw MainWindow instance. 
        Use with caution; prefer specific methods when available.
        """
        return self._manager.get_main_window()

    def add_export_action(self, label: str, callback: Callable):
        """
        Register a custom export action.
        
        Args:
            label: Text to display in the Export menu (e.g., "Export as MyFormat...").
            callback: Function to call when triggered.
        """
        self._manager.register_export_action(self._plugin_name, label, callback)

    def register_optimization_method(self, method_name: str, callback: Callable[[Any], bool]):
        """
        Register a custom 3D optimization method.
        
        Args:
            method_name: Name of the method to display in 3D Optimization menu.
            callback: Function taking (rdkit_mol) -> bool (success).
                      Modifies the molecule in-place.
        """
        self._manager.register_optimization_method(self._plugin_name, method_name, callback)

    def register_file_opener(self, extension: str, callback: Callable[[str], None]):
        """
        Register a handler for opening a specific file extension.
        
        Args:
            extension: File extension including dot, e.g. ".xyz".
            callback: Function taking (file_path) -> None.
                      Should load the file into the main window.
        """
        self._manager.register_file_opener(self._plugin_name, extension, callback)

    def register_file_opener(self, extension: str, callback: Callable[[str], None]):
        """
        Register a handler for opening a specific file extension.
        
        Args:
            extension: File extension including dot, e.g. ".xyz".
            callback: Function taking (file_path) -> None.
                      Should load the file into the main window.
        """
        self._manager.register_file_opener(self._plugin_name, extension, callback)

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
    def register_load_handler(self, callback: Callable[[dict], None]):
        """
        Register a callback to restore state from the project file.
        
        Args:
            callback: Function receiving the dict of saved data.
        """
        self._manager.register_load_handler(self._plugin_name, callback)

    def register_3d_context_menu(self, callback: Callable, label: str):
        """Deprecated: This method does nothing. Kept for backward compatibility."""
        print(f"Warning: Plugin '{self._plugin_name}' uses deprecated 'register_3d_context_menu'. This API has been removed.")

    def register_3d_style(self, style_name: str, callback: Callable[[Any, Any], None]):
        """
        Register a custom 3D rendering style.
        
        Args:
            style_name: Name of the style (must be unique).
            callback: Function taking (main_window, mol) -> None.
                      Should fully handle drawing the molecule in the 3D view.
        """
        self._manager.register_3d_style(self._plugin_name, style_name, callback)






class Plugin3DController:
    """Helper to manipulate the 3D scene."""
    def __init__(self, main_window):
        self._mw = main_window

    def set_atom_color(self, atom_index: int, color_hex: str):
        """
        Set the color of a specific atom in the 3D view.
        Args:
            atom_index: RDKit atom index.
            color_hex: Hex string e.g., "#FF0000".
        """
        # This will need to hook into the actual 3D view logic
        if hasattr(self._mw, 'main_window_view_3d'):
             # Logic to update color map and trigger redraw
             # For now we can assume we might need to expose a method in view_3d
             self._mw.main_window_view_3d.update_atom_color_override(atom_index, color_hex)
             self._mw.plotter.render()

    def set_bond_color(self, bond_index: int, color_hex: str):
        """
        Set the color of a specific bond in the 3D view.

        Args:
             bond_index: RDKit bond index.
             color_hex: Hex string e.g., "#00FF00".
        """
        if hasattr(self._mw, 'main_window_view_3d'):
            self._mw.main_window_view_3d.update_bond_color_override(bond_index, color_hex)
            self._mw.plotter.render()
