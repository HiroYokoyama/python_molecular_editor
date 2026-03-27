# Plugin System Internals

This document explains the internal architecture of the MoleditPy plugin system. For a guide on *how to write* plugins, please refer to the **Plugin Development Manual**.

## Core Components

The plugin system is composed of three main classes:

| Class | File | Location | Role |
| :--- | :--- | :--- | :--- |
| **`PluginManager`** | `plugin_manager.py` | `plugins/` | The backend engine. Discovers, loads, and manages the lifecycle of plugins. |
| **`PluginContext`** | `plugin_interface.py` | `plugins/` | The API surface exposed to plugins. Scoped interface for developers. |
| **`PluginManagerWindow`** | `plugin_manager_window.py` | `ui/` | The frontend UI. Allows users to view, install, and remove plugins. |

---

## 1. Discovery and Loading (`PluginManager`)

### Location
Plugins are loaded from `~/.moleditpy/plugins` (User Home Directory).

### Discovery Process (`discover_plugins`)
The manager scans the plugin directory recursively:
1.  **Packages**: If a folder contains `__init__.py`, it is treated as a single "Package Plugin". The folder name becomes the module name.
2.  **Scripts**: If a folder does *not* contain `__init__.py`, it is treated as a "Category Folder". The manager scans it for individual `.py` files.
    - Example: `plugins/Analysis/docking.py` -> Category: `Analysis`, Module: `docking`.

### Metadata Extraction
Before importing, the manager uses Python's `ast` (Abstract Syntax Tree) module to safely parse the file without executing code. It looks for global variables:
- `PLUGIN_NAME`
- `PLUGIN_VERSION`
- `PLUGIN_AUTHOR`
- `PLUGIN_DESCRIPTION`
- `PLUGIN_CATEGORY`

### Initialization
1.  **`initialize(context)`**: The preferred entry point. The manager passes a `PluginContext` instance tailored to that plugin.
2.  **`autorun(main_window)`**: (Legacy) Called with the raw MainWindow.
3.  **`run(main_window)`**: (Legacy) Manually triggered via menu.

---

## 2. API Surface (`PluginContext`)

The `PluginContext` acts as a facade. When a plugin calls `context.add_menu_action(...)`, the context forwards this to `PluginManager.register_menu_action(...)`, automatically injecting the plugin's name.

### Extension Points
Plugins can register into specific registries in the `PluginManager`:
- **Menus & Toolbars**: `add_menu_action`, `add_toolbar_action`.
- **File I/O**: `register_file_opener` (for custom extensions), `add_export_action`.
- **Analysis**: `add_analysis_tool`.
- **3D Compute**: `register_optimization_method`.
- **State Persistence**: `register_save_handler` / `register_load_handler` (saved into the `.molproj` JSON).
- **Drag & Drop**: `register_drop_handler`.

---

## 3. Management UI (`PluginManagerWindow`)

- **Table View**: Lists all discovered plugins with their status (Loaded, Error, etc.).
- **Installation**:
    - Accepts Drag & Drop of `.py` files or `.zip` archives.
    - **Smart Unzipping**: Auto-detects if a ZIP contains a single root folder or flat files and handles extraction accordingly to keep the plugin folder clean.
- **Removal**: Deletes the plugin file or folder and triggers a hot reload.
