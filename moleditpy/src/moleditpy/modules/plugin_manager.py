#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
plugin_manager.py
Manages discovery, loading, and execution of external plugins.
"""

import os
import sys
import shutil
import importlib.util
import traceback
import ast
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtCore import QUrl
from PyQt6.QtWidgets import QMessageBox

try:
    from .plugin_interface import PluginContext
except ImportError:
    # Fallback if running as script
    from modules.plugin_interface import PluginContext

class PluginManager:
    def __init__(self, main_window=None):
        self.plugin_dir = os.path.join(os.path.expanduser('~'), '.moleditpy', 'plugins')
        self.plugins = [] # List of dicts
        self.main_window = main_window
        
        # Registries for actions
        self.menu_actions = [] # List of (plugin_name, path, callback, text, icon, shortcut)
        self.toolbar_actions = [] 
        self.drop_handlers = [] # List of (priority, plugin_name, callback)
        
        # Extended Registries (Added to prevent lazy initialization "monkey patching")
        self.export_actions = [] 
        self.optimization_methods = {}
        self.file_openers = {}
        self.analysis_tools = []
        self.save_handlers = {}
        self.load_handlers = {}
        self.custom_3d_styles = {} # style_name -> {'plugin': name, 'callback': func}

    def get_main_window(self):
        return self.main_window

    def set_main_window(self, mw):
        self.main_window = mw

    def ensure_plugin_dir(self):
        """Creates the plugin directory if it doesn't exist."""
        if not os.path.exists(self.plugin_dir):
            try:
                os.makedirs(self.plugin_dir)
            except OSError as e:
                print(f"Error creating plugin directory: {e}")

    def open_plugin_folder(self):
        """Opens the plugin directory in the OS file explorer."""
        self.ensure_plugin_dir()
        QDesktopServices.openUrl(QUrl.fromLocalFile(self.plugin_dir))

    def install_plugin(self, file_path):
        """Copies a plugin file to the plugin directory."""
        self.ensure_plugin_dir()
        try:
            filename = os.path.basename(file_path)
            dest_path = os.path.join(self.plugin_dir, filename)
            shutil.copy2(file_path, dest_path)
            # Reload plugins after install
            if self.main_window:
                self.discover_plugins(self.main_window)
            return True, f"Installed {filename}"
        except Exception as e:
            return False, str(e)

    def discover_plugins(self, parent=None):
        """
        Recursively scans the plugin directory.
        Supports both legacy autorun(parent) and new initialize(context).
        """
        if parent:
            self.main_window = parent
            
        self.ensure_plugin_dir()
        self.plugins = []
        self.menu_actions = []
        self.toolbar_actions = []
        self.drop_handlers = []
        
        # Clear extended registries
        self.export_actions = [] 
        self.optimization_methods = {}
        self.file_openers = {}
        self.analysis_tools = []
        self.save_handlers = {}
        self.load_handlers = {}
        self.custom_3d_styles = {}
        
        if not os.path.exists(self.plugin_dir):
            return []

        for root, dirs, files in os.walk(self.plugin_dir):
            # Modify dirs in-place to skip hidden directories and __pycache__
            dirs[:] = [d for d in dirs if not d.startswith('__') and d != '__pycache__']
            
            for filename in files:
                if filename.endswith(".py") and not filename.startswith("__"):
                    filepath = os.path.join(root, filename)
                    rel_folder = os.path.relpath(root, self.plugin_dir)
                    if rel_folder == '.':
                        rel_folder = ""
                        
                    try:
                        module_name = os.path.splitext(os.path.relpath(filepath, self.plugin_dir))[0].replace(os.sep, '.')
                        
                        spec = importlib.util.spec_from_file_location(module_name, filepath)
                        if spec and spec.loader:
                            module = importlib.util.module_from_spec(spec)
                            sys.modules[spec.name] = module 
                            spec.loader.exec_module(module)

                            # --- Metadata Extraction ---
                            plugin_name = getattr(module, 'PLUGIN_NAME', filename[:-3])
                            plugin_version = getattr(module, 'PLUGIN_VERSION', getattr(module, '__version__', 'Unknown'))
                            plugin_author = getattr(module, 'PLUGIN_AUTHOR', getattr(module, '__author__', 'Unknown'))
                            plugin_desc = getattr(module, 'PLUGIN_DESCRIPTION', getattr(module, '__doc__', ''))
                            
                            # Clean up docstring if used as description
                            if plugin_desc:
                                plugin_desc = plugin_desc.strip().split('\n')[0]

                            # check for interface compliance
                            has_run = hasattr(module, 'run') and callable(module.run)
                            has_autorun = hasattr(module, 'autorun') and callable(module.autorun)
                            has_init = hasattr(module, 'initialize') and callable(module.initialize)
                            
                            status = "Loaded"
                            
                            # Execute loading logic
                            if has_init:
                                context = PluginContext(self, plugin_name)
                                try:
                                    module.initialize(context)
                                except Exception as e:
                                    status = f"Error (Init): {e}"
                                    print(f"Plugin {plugin_name} initialize error: {e}")
                                    traceback.print_exc()
                            elif has_autorun:
                                try:
                                    if self.main_window:
                                        module.autorun(self.main_window)
                                    else:
                                        status = "Skipped (No MW)"
                                except Exception as e:
                                    status = f"Error (Autorun): {e}"
                                    print(f"Plugin {plugin_name} autorun error: {e}")
                                    traceback.print_exc()
                            elif not has_run:
                                status = "No Entry Point"

                            self.plugins.append({
                                'name': plugin_name,
                                'version': plugin_version,
                                'author': plugin_author,
                                'description': plugin_desc,
                                'module': module,
                                'rel_folder': rel_folder,
                                'status': status,
                                'filepath': filepath,
                                'has_run': has_run # for menu manual run
                            })
                            
                    except Exception as e:
                        print(f"Failed to load plugin {filename}: {e}")
                        traceback.print_exc()
        
        return self.plugins

    def run_plugin(self, module, main_window):
        """Executes the plugin's run method (Legacy manual trigger)."""
        try:
            module.run(main_window)
        except Exception as e:
            QMessageBox.critical(main_window, "Plugin Error", f"Error running plugin '{getattr(module, 'PLUGIN_NAME', 'Unknown')}':\n{e}")
            traceback.print_exc()

    # --- Registration Callbacks ---
    def register_menu_action(self, plugin_name, path, callback, text, icon, shortcut):
        self.menu_actions.append({
            'plugin': plugin_name, 'path': path, 'callback': callback,
            'text': text, 'icon': icon, 'shortcut': shortcut
        })
    
    def register_toolbar_action(self, plugin_name, callback, text, icon, tooltip):
        self.toolbar_actions.append({
            'plugin': plugin_name, 'callback': callback, 
            'text': text, 'icon': icon, 'tooltip': tooltip
        })


        
    def register_drop_handler(self, plugin_name, callback, priority):
        self.drop_handlers.append({
            'priority': priority, 'plugin': plugin_name, 'callback': callback
        })
        # Sort by priority desc
        self.drop_handlers.sort(key=lambda x: x['priority'], reverse=True)

    def register_export_action(self, plugin_name, label, callback):
        self.export_actions.append({
            'plugin': plugin_name, 'label': label, 'callback': callback
        })

    def register_optimization_method(self, plugin_name, method_name, callback):
        # Key by upper-case method name for consistency
        self.optimization_methods[method_name.upper()] = {
            'plugin': plugin_name, 'callback': callback, 'label': method_name
        }

    def register_file_opener(self, plugin_name, extension, callback):
        # Normalize extension to lowercase
        ext = extension.lower()
        if not ext.startswith('.'):
            ext = '.' + ext
        self.file_openers[ext] = {
            'plugin': plugin_name, 'callback': callback
        }

    # Analysis Tools registration
    def register_analysis_tool(self, plugin_name, label, callback):
        self.analysis_tools.append({'plugin': plugin_name, 'label': label, 'callback': callback})

    # State Persistence registration
    def register_save_handler(self, plugin_name, callback):
        self.save_handlers[plugin_name] = callback

    def register_load_handler(self, plugin_name, callback):
        self.load_handlers[plugin_name] = callback

    def register_3d_style(self, plugin_name, style_name, callback):
        self.custom_3d_styles[style_name] = {
            'plugin': plugin_name, 'callback': callback
        }

    def get_plugin_info_safe(self, file_path):
        """Extracts plugin metadata using AST parsing (safe, no execution)."""
        info = {
            'name': os.path.basename(file_path),
            'version': 'Unknown',
            'author': 'Unknown',
            'description': ''
        }
        try:
             with open(file_path, "r", encoding="utf-8") as f:
                 tree = ast.parse(f.read())
             
             for node in tree.body:
                 if isinstance(node, ast.Assign):
                     for target in node.targets:
                         if isinstance(target, ast.Name):
                             # Helper to extract value
                             val = None
                             if isinstance(node.value, ast.Constant): # Py3.8+
                                 val = node.value.value
                             elif hasattr(ast, 'Str') and isinstance(node.value, ast.Str): # Py3.7 and below
                                 val = node.value.s
                             
                             if val is not None:
                                 if target.id == 'PLUGIN_NAME':
                                     info['name'] = val
                                 elif target.id == 'PLUGIN_VERSION':
                                     info['version'] = val
                                 elif target.id == 'PLUGIN_AUTHOR':
                                     info['author'] = val
                                 elif target.id == 'PLUGIN_DESCRIPTION':
                                     info['description'] = val
                                 elif target.id == '__version__' and info['version'] == 'Unknown':
                                     info['version'] = val
                                 elif target.id == '__author__' and info['author'] == 'Unknown':
                                     info['author'] = val
                 
                 # Docstring extraction
                 if isinstance(node, ast.Expr):
                     val = None
                     if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                          val = node.value.value
                     elif hasattr(ast, 'Str') and isinstance(node.value, ast.Str):
                          val = node.value.s
                          
                     if val and not info['description']:
                          info['description'] = val.strip().split('\n')[0]

        except Exception as e:
            print(f"Error parsing plugin info: {e}")
        return info



