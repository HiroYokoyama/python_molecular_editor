# Plugin System Internals

This document explains the internal architecture of the MoleditPy plugin system. For a guide on *how to write* plugins, please refer to the [Plugin Development Manual (V4)](PLUGIN_DEVELOPMENT_MANUAL_V4.md) (or the archived [V3](old/PLUGIN_DEVELOPMENT_MANUAL_V3.md)).

## Overview

The MoleditPy plugin system allows developers to extend the application's functionality without modifying the core codebase. It provides a managed API surface to ensure stability across architectural refactors.

| Component | Responsibility | Key Class | Location |
| :--- | :--- | :--- | :--- |
| **`PluginManager`** | The backend engine. Discovers, loads, and manages the lifecycle of plugins. | `PluginManager` | `plugins/plugin_manager.py` |
| **`PluginContext`** | The API surface exposed to plugins. Scoped interface for developers. | `PluginContext` | `plugins/plugin_interface.py` |
| **`PluginManagerWindow`** | The frontend UI. Allows users to view, install, and remove plugins. | `PluginManagerWindow` | `plugins/plugin_manager_window.py` |

---

## 1. Discovery and Loading

### Location
Plugins are discovered from the `~/.moleditpy/plugins` directory in the user's home folder.

### Discovery Strategy
The `PluginManager` scans the directory for two types of extensions:
1.  **Package Plugins**: Folders containing an `__init__.py` file. The folder name serves as the module name.
2.  **Script Plugins**: Individual `.py` files within "Category Folders" (e.g., `plugins/Analysis/calc.py`).

### Metadata Parsing
The manager uses Python's `ast` (Abstract Syntax Tree) module to safely extract metadata (e.g., `PLUGIN_NAME`, `PLUGIN_VERSION`) without executing the plugin's code. This allows the Plugin Manager UI to display information even for plugins that might have syntax errors or missing dependencies.

---

## 2. Extension Points

Plugins interact with the application through the `PluginContext`. Key extension points include:

- **UI Customization**: Adding actions to the main menu or toolbars.
- **File I/O**: Registering custom file openers or exporters for specialized chemical formats.
- **Analysis**: Adding new tools to the Analysis window or the 3D viewer.
- **State Persistence**: Registering custom save/load handlers for plugin-specific data within the `.pmeprj` project file.

---

## 3. Plugin Management UI

The `PluginManagerWindow` (`plugins/plugin_manager_window.py`) provides a user-friendly interface for:
- **Status Monitoring**: Viewing which plugins are successfully loaded and which encountered errors.
- **Installation**: Drag-and-drop support for `.py` or `.zip` files, including smart unzipping logic.
- **Hot Reloading**: Removing a plugin or adding a new one triggers an immediate refresh of the application menus.
