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
        Recursively scans the plugin directory for .py files and attempts to import them.
        Ignores __pycache__ and other directories starting with "__".
        Returns a list of valid loaded plugins.
        """
        self.ensure_plugin_dir()
        self.plugins = []
        
        if not os.path.exists(self.plugin_dir):
            return []

        for root, dirs, files in os.walk(self.plugin_dir):
            # Modify dirs in-place to skip hidden directories and __pycache__
            dirs[:] = [d for d in dirs if not d.startswith('__') and d != '__pycache__']
            
            for filename in files:
                if filename.endswith(".py") and not filename.startswith("__"):
                    filepath = os.path.join(root, filename)
                    
                    # Calculate relative folder path for menu structure
                    # equivalent to: rel_path = os.path.relpath(root, self.plugin_dir)
                    # if root is plugin_dir, rel_path is '.'
                    rel_folder = os.path.relpath(root, self.plugin_dir)
                    if rel_folder == '.':
                        rel_folder = ""
                        
                    try:
                        # Unique module name based on file path to avoid conflicts
                        # e.g. plugins.subdir.myplugin
                        module_name = os.path.splitext(os.path.relpath(filepath, self.plugin_dir))[0].replace(os.sep, '.')
                        
                        spec = importlib.util.spec_from_file_location(module_name, filepath)
                        if spec and spec.loader:
                            module = importlib.util.module_from_spec(spec)
                            sys.modules[spec.name] = module 
                            spec.loader.exec_module(module)

                            # Check for required attributes
                            plugin_name = getattr(module, 'PLUGIN_NAME', filename[:-3])
                            
                            # Valid plugin if it has 'run' OR 'autorun'
                            has_run = hasattr(module, 'run') and callable(module.run)
                            has_autorun = hasattr(module, 'autorun') and callable(module.autorun)
                             
                            if has_run:
                                self.plugins.append({
                                    'name': plugin_name,
                                    'module': module,
                                    'rel_folder': rel_folder
                                })
                            
                            if has_autorun:
                                try:
                                    if parent:
                                        module.autorun(parent)
                                    else:
                                        print(f"Skipping autorun for {plugin_name}: parent not provided.")
                                except Exception as e:
                                    print(f"Error executing autorun for {filename}: {e}")
                                    traceback.print_exc()

                            if not has_run and not has_autorun:
                                print(f"Plugin {filename} skipped: Missing 'run(main_window)' or 'autorun(main_window)' function.")

                    except Exception as e:
                        # Robust error handling
                        msg = f"Failed to load plugin {filename}:\n{e}"
                        print(msg)
                        traceback.print_exc()
                        if parent:
                            # Use print/status bar instead of popups for non-critical failures during bulk load?
                            # For now, keep it visible but maybe less intrusive if many fail?
                            # sticking to original logic just in catch block
                            pass 
        
        return self.plugins

    def run_plugin(self, module, main_window):
        """Executes the plugin's run method."""
        try:
            module.run(main_window)
        except Exception as e:
            QMessageBox.critical(main_window, "Plugin Error", f"Error running plugin '{getattr(module, 'PLUGIN_NAME', 'Unknown')}':\n{e}")
            traceback.print_exc()

