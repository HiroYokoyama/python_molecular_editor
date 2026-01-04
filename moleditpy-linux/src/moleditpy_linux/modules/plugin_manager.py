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
import zipfile
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
        """Copies a plugin file to the plugin directory. Supports .py and .zip."""
        self.ensure_plugin_dir()
        try:
            # Handle trailing slash and normalize path
            file_path = os.path.normpath(file_path)
            filename = os.path.basename(file_path)
            
            if os.path.isdir(file_path):
                # Copy entire directory
                dest_path = os.path.join(self.plugin_dir, filename)
                if os.path.exists(dest_path):
                     # Option 1: Overwrite (remove then copy) - safer for clean install
                     if os.path.isdir(dest_path):
                         shutil.rmtree(dest_path)
                     else:
                         os.remove(dest_path)
                
                # Copy directory, ignoring cache files
                shutil.copytree(file_path, dest_path, ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '.git'))
                msg = f"Installed package {filename}"
            elif filename.lower().endswith('.zip'):
                # Extract ZIP contents
                with zipfile.ZipFile(file_path, 'r') as zf:
                    # Smart Extraction: Check if ZIP has a single top-level folder
                    # Fix for paths with backslashes on Windows if zip was created on Windows
                    roots = set()
                    for name in zf.namelist():
                        # Normalize path separators to forward slash for consistent check
                        name = name.replace('\\', '/')
                        parts = name.split('/')
                        if parts[0]:
                            roots.add(parts[0])
                    
                    is_nested = (len(roots) == 1)
                    
                    if is_nested:
                        # Case A: ZIP contains a single folder (e.g. MyPlugin/init.py)
                        top_folder = list(roots)[0]
                        
                        # Guard: If the single item is __init__.py, we MUST create a wrapper folder
                        # otherwise we pollute the plugin_dir root.
                        if top_folder == "__init__.py":
                            is_nested = False
                    
                    if is_nested:
                        # Case A (Confirmed): Extract directly
                        dest_path = os.path.join(self.plugin_dir, top_folder)
                        
                        # Clean Install: Remove existing folder to prevent stale files
                        if os.path.exists(dest_path):
                             if os.path.isdir(dest_path):
                                 shutil.rmtree(dest_path)
                             else:
                                 os.remove(dest_path)
                        
                        zf.extractall(self.plugin_dir)
                        msg = f"Installed package {top_folder} (from ZIP)"
                    else:
                        # Case B: ZIP is flat (e.g. file1.py, file2.py or just __init__.py)
                        # Extract into a new folder named after the ZIP file
                        folder_name = os.path.splitext(filename)[0]
                        dest_path = os.path.join(self.plugin_dir, folder_name)
                        
                        if os.path.exists(dest_path):
                            if os.path.isdir(dest_path):
                                shutil.rmtree(dest_path)
                            else:
                                os.remove(dest_path)
                        
                        os.makedirs(dest_path)
                        zf.extractall(dest_path)
                        msg = f"Installed package {folder_name} (from Flat ZIP)"
            else:
                # Standard file copy
                dest_path = os.path.join(self.plugin_dir, filename)
                if os.path.exists(dest_path):
                    if os.path.isdir(dest_path):
                        shutil.rmtree(dest_path)
                shutil.copy2(file_path, dest_path)
                msg = f"Installed {filename}"

            # Reload plugins after install
            if self.main_window:
                self.discover_plugins(self.main_window)
            return True, msg
        except Exception as e:
            return False, str(e)

    def discover_plugins(self, parent=None):
        """
        Hybrid discovery:
        - Folders with '__init__.py' -> Treated as single package plugin.
        - Folders without '__init__.py' -> Treated as category folders (scan for .py inside).
        """
        if parent:
            self.main_window = parent
            
        self.ensure_plugin_dir()
        self.plugins = []
        # Clear registries
        self.menu_actions = []
        self.toolbar_actions = []
        self.drop_handlers = []
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
            # Exclude hidden directories
            dirs[:] = [d for d in dirs if not d.startswith('__') and d != '__pycache__']
            
            # [Check] Is current dir a package (plugin body)?
            if "__init__.py" in files:
                # === Case 1: Package Plugin (Folder is the plugin) ===
                
                # Stop recursion into this folder
                dirs[:] = []
                
                entry_point = os.path.join(root, "__init__.py")
                # Category is relative path to parent folder
                rel_path = os.path.relpath(os.path.dirname(root), self.plugin_dir)
                category = rel_path if rel_path != "." else ""
                
                # Module name is the folder name
                module_name = os.path.basename(root)
                
                self._load_single_plugin(entry_point, module_name, category)
                
            else:
                # === Case 2: Category Folder (Load individual .py files) ===
                
                # Category is relative path to current folder
                rel_path = os.path.relpath(root, self.plugin_dir)
                category = rel_path if rel_path != "." else ""

                for filename in files:
                    if filename.endswith(".py") and not filename.startswith("__"):
                        entry_point = os.path.join(root, filename)
                        module_name = os.path.splitext(filename)[0]
                        
                        self._load_single_plugin(entry_point, module_name, category)
        
        return self.plugins

    def _load_single_plugin(self, filepath, module_name, category):
        """Common loading logic for both single-file and package plugins."""
        try:
            # Ensure unique module name by including category path
            # e.g. Analysis.Docking
            unique_module_name = f"{category.replace(os.sep, '.')}.{module_name}" if category else module_name
            unique_module_name = unique_module_name.strip(".")

            spec = importlib.util.spec_from_file_location(unique_module_name, filepath)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                sys.modules[spec.name] = module
                
                # Inject category info
                module.PLUGIN_CATEGORY = category 
                
                spec.loader.exec_module(module)

                # --- Metadata Extraction ---
                # Metadata
                # Priority: PLUGIN_XXX > __xxx__ > Fallback
                plugin_name = getattr(module, 'PLUGIN_NAME', module_name)
                plugin_version = getattr(module, 'PLUGIN_VERSION', getattr(module, '__version__', 'Unknown'))
                plugin_author = getattr(module, 'PLUGIN_AUTHOR', getattr(module, '__author__', 'Unknown'))
                plugin_desc = getattr(module, 'PLUGIN_DESCRIPTION', getattr(module, '__doc__', ''))
                plugin_category = getattr(module, 'PLUGIN_CATEGORY', category)
                
                # Additional cleanup for docstring (strip whitespace)
                if plugin_desc is None: plugin_desc = ""
                plugin_desc = str(plugin_desc).strip()
                
                 # Handle version tuple
                if isinstance(plugin_version, tuple):
                    plugin_version = ".".join(map(str, plugin_version))
                    
                # Interface compliance
                has_run = hasattr(module, 'run') and callable(module.run)
                has_autorun = hasattr(module, 'autorun') and callable(module.autorun)
                has_init = hasattr(module, 'initialize') and callable(module.initialize)
                
                status = "Loaded"
                
                # Execute initialization
                if has_init:
                    context = PluginContext(self, plugin_name)
                    # Pass category to context if needed, currently not storing it in context directly but could be useful
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
                    'category': plugin_category, # Store category
                    'status': status,
                    'filepath': filepath,
                    'has_run': has_run
                })
                
        except Exception as e:
            print(f"Failed to load plugin {module_name}: {e}")
            traceback.print_exc()

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
                 targets = []
                 if isinstance(node, ast.Assign):
                     targets = node.targets
                 elif isinstance(node, ast.AnnAssign):
                     targets = [node.target]
                 
                 for target in targets:
                     if isinstance(target, ast.Name):
                         # Helper to extract value
                         val = None
                         if node.value: # AnnAssign might presumably not have value? (though usually does for module globals)
                             if isinstance(node.value, ast.Constant): # Py3.8+
                                 val = node.value.value
                             elif hasattr(ast, 'Str') and isinstance(node.value, ast.Str): # Py3.7 and below
                                 val = node.value.s
                             elif isinstance(node.value, ast.Tuple):
                                 # Handle version tuples e.g. (1, 0, 0)
                                 try:
                                     # Extract simple constants from tuple
                                     elts = []
                                     for elt in node.value.elts:
                                         if isinstance(elt, ast.Constant):
                                             elts.append(elt.value)
                                         elif hasattr(ast, 'Num') and isinstance(elt, ast.Num):
                                             elts.append(elt.n)
                                     val = ".".join(map(str, elts))
                                 except:
                                     pass
                         
                         if val is not None:
                             if target.id == 'PLUGIN_NAME':
                                 info['name'] = val
                             elif target.id == 'PLUGIN_VERSION':
                                 info['version'] = val
                             elif target.id == 'PLUGIN_AUTHOR':
                                 info['author'] = val
                             elif target.id == 'PLUGIN_DESCRIPTION':
                                 info['description'] = val
                             elif target.id == 'PLUGIN_CATEGORY':
                                 info['category'] = val
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



