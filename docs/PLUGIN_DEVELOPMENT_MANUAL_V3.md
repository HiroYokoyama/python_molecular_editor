# MoleditPy Plugin Development Manual (Version 3)

This manual describes the current plugin API for the composition-based MoleditPy architecture.
It replaces the old direct-`mw.<method>` style with the safer `PluginContext` interface.

> [!IMPORTANT]
> New plugins should use `initialize(context)` and the `PluginContext` methods.
> Avoid relying on raw `MainWindow` attributes unless you are doing advanced Qt or rendering work.

## 1. Plugin Basics

A plugin is either a single Python file or a package folder with an `__init__.py` entry point.

### Required metadata

Define these at module scope:

| Variable | Purpose |
| --- | --- |
| `PLUGIN_NAME` | Human-readable plugin name shown in the UI. |
| `PLUGIN_VERSION` | Plugin version string. |
| `PLUGIN_AUTHOR` | Author or team name. |
| `PLUGIN_DESCRIPTION` | Short description shown in the plugin manager. |

### Entry point

Use `initialize(context)` for new plugins.

```python
PLUGIN_NAME = "My Plugin"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Your Name"
PLUGIN_DESCRIPTION = "Does something useful."

def initialize(context):
    context.add_menu_action("Tools/My Action", lambda: run_action(context))
```

### Package layout

Folder-based plugins are supported and are the preferred format for anything non-trivial.

```text
plugins/
  MyPlugin/
    __init__.py
    logic.py
    assets/
```

Use relative imports inside the package:

```python
from .logic import run_action

PLUGIN_NAME = "My Plugin"
PLUGIN_VERSION = "1.0.0"


def initialize(context):
    context.add_plugin_menu("Utility/My Action", lambda: run_action(context))
```

## 2. Current Architecture

MoleditPy now uses a composition-based `MainWindow`. Most behavior lives in manager objects instead of direct `MainWindow` methods.

Common managers on the main window include:

- `state_manager`
- `view_3d_manager`
- `edit_3d_manager`
- `edit_actions_manager`
- `compute_manager`
- `dialog_manager`
- `io_manager`
- `export_manager`
- `ui_manager`
- `init_manager`
- `string_importer_manager`

For plugins, the important point is simple:

- Prefer `PluginContext` methods first.
- Use `context.get_main_window()` only when you need a Qt-level operation that the context does not expose.
- Do not assume old `mw.<method>` helpers still exist; many of them were moved into managers.

## 3. PluginContext Reference

`context` is the supported plugin API surface.

### 3.1 Menus and toolbar actions

#### `add_menu_action(path, callback, text=None, icon=None, shortcut=None)`
Add an action to the menu bar.

Examples:

- `File/Import/My Loader`
- `Tools/My Tool`
- `MyPlugin/Action`

#### `register_menu_action(path, text_or_callback, callback=None, icon=None, shortcut=None)`
Compatibility alias for `add_menu_action`.

Preferred style:

```python
context.add_menu_action("Tools/Run", lambda: run_tool(context))
```

Legacy style still works:

```python
context.register_menu_action("Tools/Run", "Run", run_tool)
```

#### `add_plugin_menu(path, callback, text=None, icon=None, shortcut=None)`
Adds an action under the shared `Plugin` menu.
This is the preferred shortcut when you want your plugin to live in a grouped plugin namespace.

```python
context.add_plugin_menu("Utility/Inspect", lambda: open_inspector(context))
```

#### `add_toolbar_action(callback, text, icon=None, tooltip=None)`
Add a button to the plugin toolbar.

#### `add_analysis_tool(label, callback)`
Add a tool to the top-level `Analysis` menu.
Use this for read-only data inspection or reporting tools.

#### `add_export_action(label, callback)`
Add a tool to the `Export` menu.
Use this for file writers, summaries, or export helpers.

### 3.2 Messages and feedback

#### `show_status_message(message, timeout=3000)`
Shows a message in the status bar.

### 3.3 Files and drag-and-drop

#### `register_file_opener(extension, callback, priority=0)`
Register a loader for a specific file extension.
The callback receives the file path as a string.

#### `register_drop_handler(callback, priority=0)`
Register a drop handler for files dropped onto the editor window.
The callback must return `True` if it handled the file.

### 3.4 Molecule state and undo

#### `current_mol` / `current_molecule`
Get or set the active RDKit molecule.
These are aliases.
Setting the molecule updates the 3D view through the current rendering pipeline.

#### `get_selected_atom_indices()`
Returns the selected atom indices from the current 2D/3D selection.

#### `push_undo_checkpoint()`
Pushes a state snapshot for undo history.
Call this after you finish making changes.

### 3.5 View control

#### `refresh_3d_view()`
Forces a redraw of the 3D view.
Use this after in-place visual updates.

#### `reset_3d_camera()`
Resets and fits the 3D camera to the current molecule.

#### `get_3d_controller()`
Returns a `Plugin3DController` for temporary visual overrides.

### 3.6 Windows and persistent settings

#### `register_window(window_id, window)`
Register a Qt window or dialog so the application can keep track of it.
Window IDs are namespaced by plugin name automatically.

#### `get_window(window_id)`
Fetch a previously registered window.
Returns `None` if it does not exist.

#### `get_setting(key, default=None)`
Read a persistent plugin setting.

#### `set_setting(key, value)`
Store a persistent plugin setting.
Values must be JSON-serializable.

### 3.7 Advanced access

#### `get_main_window()`
Returns the raw `MainWindow` instance.
Use this for Qt-specific work or for advanced rendering hooks that need the underlying window.

#### `plotter`
The current PyVista plotter, exposed for advanced 3D customization.

#### `scene`
The current 2D molecule scene.

## 4. Legacy Support

Some older plugin entry points still work, but they are compatibility features.

### `run(mw)`
Legacy entry point for menu-driven plugins.
Use `initialize(context)` for new work.

### `autorun(mw)`
Legacy startup hook.
Use `initialize(context)` and explicit registration instead whenever possible.

### `register_3d_context_menu(...)`
This API is deprecated and does nothing.
Do not build new plugins around it.

## 5. Working With The New Composition Model

The biggest plugin migration change is that many old `mw.<method>` calls were moved into managers.

Use these rules:

- Prefer `PluginContext` methods first.
- If you need the raw window, call `context.get_main_window()`.
- If you need 3D actions, prefer the context helpers rather than calling view internals directly.
- If you need app state, use the state/undo helpers instead of reaching into `mw` attributes.

### Useful direct manager access on `MainWindow`

If you do need the raw window, these are the current manager paths that replace many old methods:

- `mw.view_3d_manager.draw_molecule_3d(mol)`
- `mw.view_3d_manager.current_mol`
- `mw.edit_actions_manager.push_undo_state()`
- `mw.state_manager.update_window_title()`
- `mw.state_manager.update_realtime_info()`
- `mw.init_manager.current_file_path`
- `mw.init_manager.settings`
- `mw.ui_manager`
- `mw.compute_manager`

## 6. Migration Cheat Sheet

| Old pattern | New pattern |
| --- | --- |
| `mw.statusBar().showMessage(msg)` | `context.show_status_message(msg)` |
| `mw.push_undo_state()` | `context.push_undo_checkpoint()` |
| `mw.current_mol` | `context.current_mol` or `context.current_molecule` |
| `mw.draw_molecule_3d(mol)` | `context.current_mol = mol` |
| `mw.plotter.reset_camera()` | `context.reset_3d_camera()` |
| `mw.my_dialog = dlg` | `context.register_window("my_dialog", dlg)` |
| `context.register_menu_action(path, text, cb)` | `context.add_menu_action(path, cb, text)` |
| `add_menu_action("Plugin/X/...", cb)` | `context.add_plugin_menu("X/...", cb)` |
| `mw.register_menu_action(...)` | use `context.add_menu_action(...)` from plugins |
| `mw.get_context(...)` | use `PluginContext` from `initialize(context)` |

## 7. Recommended Patterns

### 7.1 Menu callback

```python
def initialize(context):
    context.add_menu_action("Tools/Do Thing", lambda: do_thing(context))
```

### 7.2 Read and modify the current molecule

```python
from rdkit import Chem

def do_thing(context):
    mol = context.current_mol
    if mol is None:
        context.show_status_message("No molecule loaded.", 2000)
        return

    mol = Chem.AddHs(mol)
    context.current_mol = mol
    context.push_undo_checkpoint()
    context.show_status_message("Molecule updated.", 2000)
```

### 7.3 3D highlighting

```python
def highlight_atoms(context):
    ctrl = context.get_3d_controller()
    for idx in context.get_selected_atom_indices():
        ctrl.set_atom_color(idx, "#00FF00")
    context.refresh_3d_view()
```

### 7.4 Singleton window

```python
from PyQt6.QtWidgets import QDialog


def toggle_panel(context):
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        return

    win = QDialog(context.get_main_window())
    context.register_window("main_panel", win)
    win.show()
```

## 8. Testing and Validation

You can test plugin logic with a mocked `PluginContext` in unit tests.

```python
from unittest.mock import MagicMock

context = MagicMock()
context.current_mol = None
```

For repo-level validation, these scripts are useful:

- `python scripts\\scan_plugin_disconnections.py`
- `python scripts\\scan_disconnections_v2.py`

Use them after larger API changes to catch stale plugin references.

## 9. Best Practices

- Keep callbacks small and fast.
- Wrap risky logic in `try/except` and report errors with `show_status_message`.
- Register persistent UI with `register_window`.
- Use `context.get_main_window()` only when the context API does not cover the task.
- Prefer menu callbacks that accept `context`, not bare access to internal globals.

## 10. Example Plugin

```python
from rdkit.Chem import Descriptors
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Molecule Stats"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Your Name"
PLUGIN_DESCRIPTION = "Shows basic molecule statistics."


def initialize(context):
    context.add_analysis_tool("Show Molecular Weight", lambda: show_weight(context))


def show_weight(context):
    mol = context.current_mol
    if mol is None:
        context.show_status_message("Load a molecule first.", 2000)
        return

    weight = Descriptors.MolWt(mol)
    QMessageBox.information(
        context.get_main_window(),
        "Molecule Stats",
        f"Molecular weight: {weight:.2f}"
    )
```

## 11. Final Notes

This manual tracks the current plugin architecture, where the application is organized around managers rather than a monolithic `MainWindow` API.
If you are updating an old plugin, start by replacing direct `mw.<method>` calls with `PluginContext` helpers, then only fall back to `get_main_window()` for the few cases that truly need it.
