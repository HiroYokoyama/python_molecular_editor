#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plugin_manager.py
Manages discovery, loading, and execution of external plugins.
"""

import os
import sys
import importlib.util
import traceback
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtCore import QUrl
from PyQt6.QtWidgets import QMessageBox

class PluginManager:
    def __init__(self):
        self.plugin_dir = os.path.join(os.path.expanduser('~'), '.moleditpy', 'plugins')
        self.plugins = [] # List of {"name": str, "module": module_obj}

    def ensure_plugin_dir(self):
        """Creates the plugin directory if it creates doesn't exist."""
        if not os.path.exists(self.plugin_dir):
            try:
                os.makedirs(self.plugin_dir)
            except OSError as e:
                print(f"Error creating plugin directory: {e}")

    def open_plugin_folder(self):
        """Opens the plugin directory in the OS file explorer."""
        self.ensure_plugin_dir()
        QDesktopServices.openUrl(QUrl.fromLocalFile(self.plugin_dir))

    def discover_plugins(self, parent=None):
        """
        Scans the plugin directory for .py files and attempts to import them.
        Returns a list of valid loaded plugins.
        """
        self.ensure_plugin_dir()
        self.plugins = []
        
        if not os.path.exists(self.plugin_dir):
            return []

        for filename in os.listdir(self.plugin_dir):
            if filename.endswith(".py") and not filename.startswith("__"):
                filepath = os.path.join(self.plugin_dir, filename)
                try:
                    # Dynamically import the module
                    spec = importlib.util.spec_from_file_location(filename[:-3], filepath)
                    if spec and spec.loader:
                        module = importlib.util.module_from_spec(spec)
                        sys.modules[spec.name] = module # helper for relative imports if needed
                        spec.loader.exec_module(module)

                        # Check for required attributes
                        plugin_name = getattr(module, 'PLUGIN_NAME', filename[:-3])
                        
                        # Validate that it has a run function
                        if hasattr(module, 'run') and callable(module.run):
                            self.plugins.append({
                                'name': plugin_name,
                                'module': module
                            })
                        else:
                            print(f"Plugin {filename} skipped: Missing 'run(main_window)' function.")
                except Exception as e:
                    # Robust error handling with user notification
                    msg = f"Failed to load plugin {filename}:\n{e}"
                    print(msg)
                    traceback.print_exc()
                    if parent:
                        QMessageBox.warning(parent, "Plugin Load Error", msg)
        
        return self.plugins

    def run_plugin(self, module, main_window):
        """Executes the plugin's run method."""
        try:
            module.run(main_window)
        except Exception as e:
            QMessageBox.critical(main_window, "Plugin Error", f"Error running plugin '{getattr(module, 'PLUGIN_NAME', 'Unknown')}':\n{e}")
            traceback.print_exc()

