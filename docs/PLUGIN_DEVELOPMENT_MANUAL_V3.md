# MoleditPy Plugin Development Manual (Version 3.0)

Welcome to the **Version 3.0** of the MoleditPy Plugin API. This version introduces a fully decoupled, namespaced architecture designed for high stability, clean memory management, and long-term maintainability.

> [!IMPORTANT]
> **API PHILOSOPHY**: In Version 3.0, it is strongly recommended that you avoid accessing the `MainWindow` directly via monkey-patching or unverified attributes. Instead, you should use the stable methods provided by the `PluginContext`. 
> 
> While we aim to maintain internal stability as much as possible, the application's core structure (e.g., manager names, attribute paths) may change significantly between major versions. By using the `PluginContext`, your plugin is protected by a stable abstraction layer that ensures long-term compatibility even when the core application is refactored.

---

## 1. Getting Started

A MoleditPy plugin is either a **single `.py` file** or a **folder package**. Both are placed in your **user plugin directory** and discovered automatically at startup.

**Plugin directory location:**
- **Windows**: `C:\Users\<YourName>\.moleditpy\plugins\`
- **Linux / macOS**: `~/.moleditpy/plugins/`

Create the folder if it does not exist, drop your `.py` file in, and restart the app (or use **Plugins > Reload All Plugins**).

### 1.1 Quick Start: The Minimal Plugin

The absolute smallest working plugin. Save as `~/.moleditpy/plugins/hello_world.py`:

```python
PLUGIN_NAME = "Hello World"

def initialize(context):
    context.add_menu_action("Tools/Say Hello", lambda: context.show_status_message("Hello!"))
```

Two lines. No imports needed for basic status messages. The plugin appears under **Tools > Say Hello** after the next reload.

### 1.2 Single-File vs. Folder Plugins

| | Single file (`my_plugin.py`) | Folder package (`MyPlugin/`) |
|---|---|---|
| **Best for** | Simple tools, one-file scripts | Complex tools with sub-modules or assets |
| **Entry point** | `my_plugin.py` | `MyPlugin/__init__.py` |
| **Imports** | Top-level only | Use relative imports (`from .logic import …`) |

Start with a single file. Promote to a folder package only when the file needs splitting.

### 1.3 Plugin Metadata
Define these at the top of your script. They are used for the UI and the internal registry.

| Variable | Description |
| :--- | :--- |
| `PLUGIN_NAME` | **(Required)** The human-readable name of the plugin. |
| `PLUGIN_VERSION` | Version string (e.g., `"1.0.2"` or `"2026.03.31"`). |
| `PLUGIN_AUTHOR` | Name of the developer. |
| `PLUGIN_DESCRIPTION` | A short summary shown in the Plugin Manager. |
| `PLUGIN_CATEGORY` | Optional category (e.g., `"Analysis"`, `"Visualization"`). |

### 1.4 Folder-based Plugins (Packages)
For complex plugins, use a folder structure. MoleditPy will treat the folder as a single plugin if it contains an `__init__.py`.

**Structure:**
```text
plugins/
  └── HighValuePlugin/      <-- Folder Name
       ├── __init__.py      <-- Entry point & Metadata (REQUIRED)
       ├── logic.py         <-- Sub-module
       └── assets/          <-- Icons or data
```

In your `__init__.py`, use **relative imports**:
```python
from .logic import process_data

PLUGIN_NAME = "High Value Tool"
PLUGIN_VERSION = "1.0"

def initialize(context):
    context.add_menu_action("Tools/Process", lambda: process_data(context))
```

---

## 2. Core API: The PluginContext Reference

The `context` object passed to `initialize(context)` is your safe proxy to the application's core logic.

### Quick Reference

| Category | Method / Property | Description |
|---|---|---|
| **UI** | `add_menu_action(path, cb, ...)` | Add item to any main menu |
| **UI** | `add_plugin_menu(path, cb, ...)` | Add item inside the Plugins menu |
| **UI** | `add_analysis_tool(label, cb)` | Add tool to the Analysis menu |
| **UI** | `add_export_action(label, cb)` | Add option to the Export menu |
| **UI** | `add_toolbar_action(cb, text, ...)` | Add button to the Plugin Toolbar |
| **UI** | `show_status_message(msg, ms)` | Temporary message in status bar |
| **Files** | `register_file_opener(ext, cb, priority)` | Handle a file extension (Import + CLI) |
| **Files** | `register_drop_handler(cb, priority)` | Handle drag-and-drop onto the window |
| **Molecule** | `current_molecule` | Get / set the active RDKit mol |
| **Molecule** | `get_selected_atom_indices()` | Indices of user-selected atoms |
| **Molecule** | `push_undo_checkpoint()` | Snapshot state to undo stack |
| **3D** | `refresh_3d_view()` | Lightweight 3D redraw |
| **3D** | `draw_molecule_3d(mol)` | Full 3D scene rebuild |
| **3D** | `reset_3d_camera()` | Fit camera to molecule |
| **3D** | `get_3d_controller()` | Atom/bond color overrides |
| **3D** | `register_3d_style(name, cb)` | Custom visualization mode |
| **3D** | `register_optimization_method(name, cb)` | Custom geometry optimizer |
| **3D** | `plotter` | Direct PyVista plotter access |
| **Project** | `register_save_handler(cb)` | Save plugin data to `.pmeprj` |
| **Project** | `register_load_handler(cb)` | Restore plugin data from `.pmeprj` |
| **Project** | `register_document_reset_handler(cb)` | Reset on File > New |
| **Settings** | `get_setting(key, default)` | Read a persisted setting |
| **Settings** | `set_setting(key, value)` | Write a persisted setting |
| **Windows** | `register_window(id, win)` | Keep a Qt window alive |
| **Windows** | `get_window(id)` | Retrieve a registered window |
| **Access** | `get_main_window()` | Raw `MainWindow` (use sparingly) |
| **Access** | `scene` | Direct 2D `MoleculeScene` access |

---

### 2.1 UI & Feedback

The following methods allow your plugin to communicate with the user via the main interface.

#### `show_status_message(message, timeout=3000)`
Display a temporary message in the application's bottom status bar.
- **message** (`str`): The text to display.
- **timeout** (`int`): Duration in milliseconds before the message disappears. Default is 3000ms.

#### `add_menu_action(path, callback, text=None, icon=None, shortcut=None)`
Register a custom item in the main menu. MoleditPy will automatically create any sub-menus defined in the path.
- **path** (`str`): The full menu path.
    - Use `File/My Action` to add to existing menus.
    - Use `MyPlugin/Action` to create a new top-level menu.
- **callback** (`Callable`): Function to execute when the action is triggered.
- **text** (`str`, optional): The label of the menu item. Defaults to the last part of `path`.
- **icon** (`str`, optional): Path to an image file or a standard icon name.
- **shortcut** (`str`, optional): Keyboard shortcut (e.g., `"Ctrl+Shift+X"`).

> [!NOTE]
> `register_menu_action` is a deprecated alias kept for V2 compatibility. New plugins should always use `add_menu_action`.

#### `add_plugin_menu(path, callback, text=None, icon=None, shortcut=None)`
Register an action nested inside the **Plugins** menu. This is the preferred way to keep the main menu bar clean if your plugin has many tools.
- **path** (`str`): The sub-path within the Plugin menu (e.g., `"Utils/My Tool"`).
- **callback** (`Callable`): Function to execute.
- **text** (`str`, optional): Label for the action.
- **icon** (`str`, optional): Icon path.
- **shortcut** (`str`, optional): Keyboard shortcut.

#### `add_toolbar_action(callback, text, icon=None, tooltip=None)`
Add a button to the dedicated **Plugin Toolbar**. The toolbar is hidden until at least one plugin registers an action — it appears automatically the first time this method is called.
- **callback** (`Callable`): Function to execute when clicked.
- **text** (`str`): Label for the button.
- **icon** (`str`, optional): Icon path.
- **tooltip** (`str`, optional): Hover text.

#### `add_analysis_tool(label, callback)`
Register a tool in the top-level **Analysis** menu. This is the preferred location for non-modifying data processing tools.
- **label** (`str`): Text to display in the menu.
- **callback** (`Callable`): Function to execute.

#### `add_export_action(label, callback)`
Register an action in the **Export** menu. Use this for custom file formats or data summaries.
- **label** (`str`): Text to display (e.g., `"Export as MyFormat..."`).
- **callback** (`Callable`): Function to execute.

---

### 2.2 Files & Interoperability

These methods allow your plugin to handle external files and drag-and-drop events.

#### `register_file_opener(extension, callback, priority=0)`
Register a handler for opening a specific file type. Used for the **File > Import** menu and **Command Line** startup. If multiple plugins register the same extension, all appear in the Import menu as separate entries (e.g., "Import .xyz (Plugin A)"). For CLI/open commands, the plugin with the highest priority wins.
- **extension** (`str`): File extension including the dot (e.g., `".xyz"`, `".cub"`).
- **callback** (`Callable[[str], None]`): Function that receives the absolute file path and handles the loading logic.
- **priority** (`int`): Higher values run first. Use **negative values** (e.g., `-1`) for fallback handlers that only run when no other plugin claims the file.

```python
context.register_file_opener(".xyz", open_xyz)          # standard
context.register_file_opener(".xyz", my_opener, priority=100)   # takes precedence
context.register_file_opener(".xyz", fallback, priority=-1)     # last resort
```

#### `register_drop_handler(callback, priority=0)`
Register a handler for files dropped onto the main 2D/3D editor window. **This is distinct from `register_file_opener`** — drop handlers are only triggered by drag-and-drop, not by the Import menu or command-line file arguments. Register both if you need both paths.
- **callback** (`Callable[[str], bool]`): Function that receives the dropped file path. Must return `True` if it successfully handled the file, `False` to pass to the next handler.
- **priority** (`int`): Handlers with higher priority are checked first. Use negative values (e.g., `-1`) for fallback handlers that run only when nothing else claims the file.

---

### 2.3 Molecular State & Undo

Methods for interacting with the active molecule and managing the undo stack.

#### `current_molecule` (Property)
Get or set the active RDKit molecule object.
- **Getter**: Returns the `rdkit.Chem.Mol` currently loaded in the editor.
- **Setter**: Replaces the active molecule. Automatically triggers 2D/3D redrawing of the scene.

> [!NOTE]
> `current_mol` is a shorthand alias for `current_molecule` — both are available.

#### `get_selected_atom_indices()`
Returns the RDKit indices of atoms currently selected by the user in either the 2D canvas or the 3D viewer.

#### `push_undo_checkpoint()`
Snapshots the current application state and adds it to the Undo history.
> [!IMPORTANT]
> **Usage Timing**: You should call this method **AFTER** you have finished modifying the molecule. The system only pushes a new state if it detects a difference from the previous one.

---

### 2.4 Lifecycle & Project Data

Manage how your plugin interacts with the application's overall lifecycle and the `.pmeprj` project format.

#### `register_save_handler(callback)`
Register a callback to save custom state into the application's project file.
- **callback** (`Callable[[], dict]`): Must return a dictionary of JSON-serializable data.

#### `register_load_handler(callback)`
Register a callback to restore custom state from the project file.
- **callback** (`Callable[[dict], None]`): Receives the dictionary previously saved by your save handler.

#### `register_document_reset_handler(callback)`
Register a callback to be invoked when a new document is created (**File > New**). Use this to reset your plugin's internal state when the user clears all data.
- **callback** (`Callable[[], None]`): Function with no arguments that resets plugin state.
- **Note**: This handler is called AFTER the main application has cleared all molecular data but BEFORE the "Cleared" message.

---

### 2.5 3D Visualization & Engine

These methods and properties allow your plugin to extend the core rendering and computational capabilities of MoleditPy.

#### `refresh_3d_view()`
Force the 3D window to redraw. Use this after performing minor visual changes (like color overrides) or manual coordinate updates.

#### `draw_molecule_3d(mol)`
Directly trigger a full redraw of the 3D scene using a specific RDKit molecule. This is more intensive than `refresh_3d_view()` as it rebuilds all 3D actors.
- **mol** (`rdkit.Chem.Mol`): The molecule to render.

#### `reset_3d_camera()`
Zoom in and re-center the 3D viewport to perfectly fit the current molecule.

#### `get_3d_controller()`
Returns a `Plugin3DController` instance (see Section 3). Use this for high-level visual overrides (atom/bond colors).

#### `register_3d_style(style_name, callback)`
Register a fully custom 3D visualization mode. This allows you to completely bypass the standard engine's rendering for specific research needs.
- **style_name** (`str`): Unique name for the style (e.g., `"vdw_surface"`).
- **callback** (`Callable[[MainWindow, rdkit.Chem.Mol], None]`): Function responsible for the entire drawing process. Access the PyVista plotter via `mw.plotter`.

#### `register_optimization_method(method_name, callback)`
Add a custom geometry optimizer to the **Compute > Optimize Geometry** menu.
- **method_name** (`str`): Name as it appears in the menu.
- **callback** (`Callable[[rdkit.Chem.Mol], bool]`): Function that receives the RDKit molecule. It should modify the coordinates in-place and return `True` on success.

#### `plotter` (Property)
Direct access to the `pyvista.Plotter` instance. Use for adding custom actors, text, or shapes to the 3D scene (e.g., `context.plotter.add_mesh(...)`).

---

### 2.6 Window & State Management (Namespaced)

In V2, plugins often attached windows directly to `mw`. In V3, use the **Namespaced Registry**. This keeps your windows alive in memory and prevents ID collisions with other plugins.

| Method | Description |
| :--- | :--- |
| `register_window(id, window)` | Stores a Qt window/dialog. The ID is automatically prefixed with your `PLUGIN_NAME`. |
| `get_window(id)` | Retrieves your registered window by ID. Returns `None` if not found or if the window was deleted. |

---

### 2.7 Settings & Persistence

Plugins can store persistent settings that are saved in the global application configuration.

#### `get_setting(key, default=None)`
Retrieve a plugin-specific setting. The key is automatically namespaced by your plugin name.
- **key** (`str`): The setting name.
- **default** (`Any`): Value to return if the setting is not found.

#### `set_setting(key, value)`
Save a plugin-specific setting. These are saved to the user's disk when the application closes.
- **key** (`str`): The setting name.
- **value** (`Any`): The value to store (must be JSON-serializable).

> [!CAUTION]
> **Persistence Limit**: Settings saved via `set_setting` reside in the application's global `settings.json`. If the user triggers **"Reset All Settings"** via the main menu, these settings will be **REMOVED**.

#### Isolated Storage (Companion JSON)
If your plugin needs to persist data that must survive an application-wide reset, or if you have complex data structures, use a dedicated JSON file instead of `set_setting`.

The correct file name and location depend on your plugin type:

| Plugin type | File to use | Why it survives updates |
|---|---|---|
| **Single-file** (`my_plugin.py`) | `my_plugin.json` — same name, same folder | Installer only overwrites the `.py` file |
| **Folder** (`MyPlugin/__init__.py`) | `MyPlugin/settings.json` | Plugin Installer backs up and restores `settings.json` before/after replacing the folder |

> [!IMPORTANT]
> For folder plugins, the file **must** be named exactly `settings.json`. The Plugin Installer only preserves this specific filename during updates — any other name will be wiped.

**Single-file plugin** — name your JSON after the plugin file:

```python
import os, json

def _config_path():
    base = os.path.splitext(os.path.abspath(__file__))[0]
    return base + ".json"          # e.g. my_plugin.json next to my_plugin.py

def load_config():
    path = _config_path()
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"my_prop": 42}         # defaults

def save_config(data):
    with open(_config_path(), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
```

**Folder plugin** — always use `settings.json` in the package directory:

```python
import os, json

def _config_path():
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(plugin_dir, "settings.json")   # MyPlugin/settings.json

def load_config():
    path = _config_path()
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"my_prop": 42}         # defaults

def save_config(data):
    with open(_config_path(), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
```

---

## 3. Visual Overrides: Plugin3DController

Obtain a controller via `ctrl = context.get_3d_controller()`. These overrides are temporary visual changes that do not modify the actual RDKit molecule coordinates or properties.

- `set_atom_color(atom_index: int, color_hex: str)`: Override the color of a specific atom (e.g., `"#FF0000"` for red).
- `set_bond_color(bond_index: int, color_hex: str)`: Override the color of a specific RDKit bond by its index.
- `set_bond_color_by_atoms(idx1: int, idx2: int, color_hex: str)`: Helper to find and color the bond connecting two atom indices.

---

## 4. Modern Workflow Examples

These patterns come up repeatedly in real plugins. Each snippet shows the helper function and how it is wired into `initialize`.

### 4.1 The "Singleton Dialog" Pattern
Always check if your window is already open before creating a new one.

```python
PLUGIN_NAME = "My Tool"

def initialize(context):
    context.add_menu_action("Tools/Open Viewer", lambda: toggle_viewer(context))

def toggle_viewer(context):
    win = context.get_window("main_panel")
    if win:
        win.show()
        win.raise_()
        return

    win = MyCustomDialog(context.get_main_window())
    context.register_window("main_panel", win)
    win.show()
```

### 4.2 The "Modify-Then-Push" Undo Pattern
Correct ordering is essential for the Undo/Redo system to function.

```python
from rdkit.Chem import AllChem

PLUGIN_NAME = "Geometry Tools"

def initialize(context):
    context.add_menu_action("Tools/Center Molecule", lambda: center_molecule(context))

def center_molecule(context):
    mol = context.current_molecule
    if not mol:
        return

    # 1. Modify coordinates
    AllChem.ComputeCanonicalTransform(mol)

    # 2. Push back — this triggers the redraw
    context.current_molecule = mol

    # 3. Checkpoint AFTER setting (only saves if state changed)
    context.push_undo_checkpoint()

    context.show_status_message("Molecule centered.", 2000)
```

### 4.3 Highlighting Selection
Color atoms the user has selected.

```python
PLUGIN_NAME = "Selection Highlighter"

def initialize(context):
    context.add_menu_action("Tools/Highlight Selection", lambda: highlight_selection(context))

def highlight_selection(context):
    indices = context.get_selected_atom_indices()
    ctrl = context.get_3d_controller()

    for idx in indices:
        ctrl.set_atom_color(idx, "#00FF00")  # green

    context.refresh_3d_view()
```

---

## 5. Migration Guide (V2 to V3)

| Legacy Pattern (Direct Access) | Modern Pattern (Best Practice) |
| :--- | :--- |
| `mw.statusBar().showMessage(m)` | `context.show_status_message(m)` |
| `mw.push_undo_state()` | `context.push_undo_checkpoint()` |
| `mw.current_mol` | `context.current_molecule` |
| `mw.plotter.reset_camera()` | `context.reset_3d_camera()` |
| `mw.my_tool = win` | `context.register_window("my_tool", win)` |
| `mw.draw_molecule_3d(...)` | `context.refresh_3d_view()` or `context.draw_molecule_3d(mol)` |
| `context.register_menu_action(path, text, cb)` | `context.add_menu_action(path, cb, text)` |
| `mw.settings.get(key)` | `context.get_setting(key)` |

> [!NOTE]
> `mw.trigger_conversion()` is still accessible via `context.get_main_window()` and delegates to the internal compute manager. Prefer `context.add_menu_action` or `context.add_analysis_tool` for triggering computation workflows from menus.

---

## 6. Advanced: MainWindow Internals (Low-Level)

While V3 encourages using `context`, sometimes you need direct access to the `MainWindow` (`mw`) for specialized Qt or PyVista operations. Obtain it via `mw = context.get_main_window()`.

### 6.1 Core Proxy Properties (Convenience)
V3 maintains several proxy properties on `mw` for backward compatibility and convenience.

| Attribute / Method | Description |
| :--- | :--- |
| `mw.plotter` | Proxy for `mw.view_3d_manager.plotter`. Direct PyVista plotter access. |
| `mw.scene` | Proxy for `mw.init_manager.scene`. Direct 2D graphics scene access. |
| `mw.current_mol` | Proxy for `mw.view_3d_manager.current_mol`. The active RDKit molecule. |
| `mw.draw_molecule_3d(mol)` | Proxy method to trigger a full 3D redraw of the scene. |

### 6.2 The Managed Architecture (V3)
In Version 3.0, most core logic is separated into specialized **Managers**. If a feature is not available as a proxy on `mw`, you should look in the corresponding manager.

| Attribute Path | Type | Description |
| :--- | :--- | :--- |
| `mw.init_manager.settings` | `dict` | Global application settings (stored in `settings.json`). |
| `mw.init_manager.splitter` | `QSplitter` | The UI divider between the 2D and 3D editor panels. |
| `mw.ui_manager` | `UIManager` | Handles status bar messages, editor modes, and UI state updates. |
| `mw.io_manager` | `IOManager` | Handles loading/saving molecules and project files. |
| `mw.compute_manager` | `ComputeManager` | Manages 3D coordinate conversion and background calculations. |
| `mw.view_3d_manager` | `View3DManager` | Manages the PyVista engine, styles, and 3D visualization. |
| `mw.state_manager` | `StateManager` | Manages the undo/redo stack and molecular data lifecycle. |

#### Example: Setting the Editor Mode
```python
mw = context.get_main_window()
# Correct V3 way to switch to Carbon draw mode
mw.ui_manager.set_mode("atom_C")
```

#### Example: Triggering Coordinate Conversion (2D to 3D)
```python
mw = context.get_main_window()
# Correct V3 way to trigger the internal conversion logic (if molecule is already loaded)
mw.compute_manager.trigger_conversion()
```

### 6.3 Common Tasks (Quick Reference)

For operations not yet wrapped by `PluginContext`, access `mw` directly:

| Task | Code |
|---|---|
| Load a MOL/SDF file | `mw.io_manager.load_mol_file("path/to/file.mol")` |
| Load XYZ for 3D viewing only | `mw.io_manager.load_xyz_for_3d_viewing("path/to/file.xyz")` |
| Clear the 3D scene | `mw.plotter.clear()` |
| Switch the active 3D style | `mw.view_3d_manager.set_3d_style("ball_and_stick")` |
| Center the 3D camera | `mw.plotter.reset_camera()` |
| Show a status message | `mw.statusBar().showMessage("msg", 3000)` |
| Switch 2D editor mode | `mw.ui_manager.set_mode("atom_C")` |
| Trigger 2D→3D conversion | `mw.compute_manager.trigger_conversion()` |

---

## 7. Legacy Support (Compatibility Mode)

MoleditPy 3.0 remains compatible with older plugin patterns, though they are considered deprecated.

### 7.1 The `run(mw)` Function

The simplest possible plugin with a `run` function. MoleditPy automatically adds an entry under the **Plugins** menu using `PLUGIN_NAME`. No `initialize` needed.

```python
# plugins/simple_counter.py
PLUGIN_NAME = "Atom Counter"
PLUGIN_VERSION = "1.0"

def run(main_window):
    """Called when the user clicks Plugins > Atom Counter."""
    mol = main_window.current_mol
    if mol is None:
        main_window.statusBar().showMessage("No molecule loaded.", 3000)
        return
    n = mol.GetNumAtoms()
    main_window.statusBar().showMessage(f"Atom count: {n}", 4000)
```

This is the fastest way to ship a one-action tool. The trade-off is that you get `main_window` directly (V2 style) instead of the stable `PluginContext` proxy. For new plugins, prefer `initialize(context)` so you benefit from the V3 API.

> [!TIP]
> You can mix both: define `run(mw)` for the quick Plugins-menu entry **and** `initialize(context)` for additional menu actions or file handlers in the same file.

```python
# Mixing run() and initialize() in the same file
PLUGIN_NAME = "Bond Inspector"
PLUGIN_VERSION = "1.0"

def run(main_window):
    """Quick Plugins-menu action."""
    mol = main_window.current_mol
    if mol:
        main_window.statusBar().showMessage(f"Bonds: {mol.GetNumBonds()}", 3000)

def initialize(context):
    """Register extra menu entries via the V3 API."""
    context.add_analysis_tool("Bond Inspector", lambda: _show_bonds(context))

def _show_bonds(context):
    mol = context.current_molecule
    if mol:
        context.show_status_message(f"Bonds: {mol.GetNumBonds()}", 3000)
```

### 7.2 The `autorun(mw)` Function

Executed immediately at startup before the main window is shown. Use sparingly — slow `autorun` code delays application launch.

```python
PLUGIN_NAME = "Startup Logger"

def autorun(main_window):
    """Runs once at startup."""
    import logging
    logging.info("MoleditPy started — Startup Logger active.")
```

Prefer `initialize(context)` for all modern plugins. `autorun` has no equivalent in the V3 API on purpose; startup side-effects belong inside `initialize`.

---

## 8. Best Practices & Troubleshooting

### 8.1 Thread Safety
**CRITICAL**: All plugin callbacks run in the Main UI Thread.
- **Fast operations** (< 100ms) are fine.
- **Slow operations** (heavy QM, large loops) will **freeze the entire application**.
- **Solution**: Use `QThread` or Python's `threading` for heavy work and return results via signals.

### 8.2 Error Handling
Always wrap your logic in `try...except` blocks. An uncaught exception in a plugin can cause the entire application to hang or crash silently.
```python
def my_callback(context):
    try:
        # unsafe logic
    except Exception as e:
        context.show_status_message(f"Plugin Error: {e}", 5000)
```

### 8.3 Memory Safety
Use `context.register_window()` for all persistent UI elements (dialogs, panels). If you don't register them, Python's garbage collector might delete the window while it's still being used by the user.

### 8.4 Atom Indices
MoleditPy uses **0-based RDKit indices**. These may differ from atom numbers in some chemical file formats. Always verify your mappings.

### 8.5 Hot Reloading
Use **Plugins > Reload All Plugins** to pick up code changes without restarting. V3 plugins survive reloads correctly as long as all persistent windows are registered via `context.register_window()`.

---

## 9. Complete Integrated Example

Combines menu registration, file import, atom selection, undo, and a singleton dialog.

```python
from rdkit import Chem
from PyQt6.QtWidgets import QMessageBox, QDockWidget, QLabel, QVBoxLayout, QWidget
from PyQt6.QtCore import Qt

PLUGIN_NAME = "Research Assistant"
PLUGIN_VERSION = "3.1.0"

def initialize(context):
    context.add_menu_action("Research/Analyze Selection", lambda: analyze(context))
    context.add_menu_action("Research/Show Panel", lambda: show_panel(context))
    context.register_file_opener(".res", lambda path: load_res(path, context))

# --- Actions ---

def analyze(context):
    indices = context.get_selected_atom_indices()
    if not indices:
        context.show_status_message("Select atoms first.", 2000)
        return
    mol = context.current_molecule
    symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in indices]
    context.show_status_message(f"Selected: {', '.join(symbols)}", 4000)

def show_panel(context):
    dock = context.get_window("panel")
    if dock:
        dock.show()
        dock.raise_()
        return

    mw = context.get_main_window()
    dock = QDockWidget("Research Assistant", mw)
    content = QWidget()
    layout = QVBoxLayout(content)
    layout.addWidget(QLabel("Research Assistant active."))
    dock.setWidget(content)
    mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    context.register_window("panel", dock)
    dock.show()

def load_res(path, context):
    try:
        with open(path, "r") as f:
            smiles = f.read().strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Could not parse SMILES in file.")
        context.current_molecule = mol
        context.push_undo_checkpoint()
        context.reset_3d_camera()
        context.show_status_message(f"Loaded {path}")
    except Exception as e:
        QMessageBox.critical(context.get_main_window(), "Load Error", str(e))
```

---

## 10. Cookbook & Examples

### 10.1 Analysis Tool: Molecular Weight
A tool that calculates properties and shows them in a popup. Use `context.current_molecule` and `context.show_status_message`.

```python
from rdkit.Chem import Descriptors
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Weight Calculator"

def initialize(context):
    def run_calc():
        mol = context.current_molecule
        if not mol:
            context.show_status_message("Error: No molecule loaded!", 2000)
            return

        weight = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()

        QMessageBox.information(
            context.get_main_window(),
            "Stats",
            f"Molecular Weight: {weight:.2f}\nAtoms: {num_atoms}"
        )

    context.add_analysis_tool("Show Molecular Weight", run_calc)
```

### 10.2 Custom 3D Component: Add Sphere
Use the `plotter` object to add persistent custom 3D geometries.

```python
import pyvista as pv

PLUGIN_NAME = "3D Decorator"

def initialize(context):
    context.add_menu_action("Visuals/Add Red Sphere", lambda: add_sphere(context))

def add_sphere(context):
    mw = context.get_main_window()
    plotter = mw.plotter

    # Create and add a sphere mesh
    sphere = pv.Sphere(radius=2.0)
    plotter.add_mesh(sphere, color="red", opacity=0.3, name="custom_plugin_sphere")

    context.refresh_3d_view()
    context.show_status_message("Sphere added to scene.")
```

### 10.3 Persistent Dock Panel (Singleton)
V3 uses `register_window` to keep panels alive and manage their unique identity across the entire session.

```python
from PyQt6.QtWidgets import QDockWidget, QLabel, QVBoxLayout, QWidget
from PyQt6.QtCore import Qt

PLUGIN_NAME = "Info Panel"

def initialize(context):
    context.add_menu_action("Tools/Show Info Panel", lambda: toggle_panel(context))

def toggle_panel(context):
    # 1. Check if the panel already exists in the manager
    dock = context.get_window("main_dock")
    if dock:
        dock.show()
        dock.raise_()
        return

    # 2. Create new panel if not found
    mw = context.get_main_window()
    dock = QDockWidget("Plugin Information", mw)
    content = QWidget()
    layout = QVBoxLayout(content)
    layout.addWidget(QLabel("Telemetry Active..."))
    dock.setWidget(content)

    mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)

    # 3. Register it so context.get_window() works next time
    context.register_window("main_dock", dock)
    dock.show()
```

### 10.4 Custom File Importer
Refactor your loaders to use `current_molecule` and `push_undo_checkpoint()` to ensure the user can revert the import.

```python
from rdkit import Chem

PLUGIN_NAME = "Simple SMILES Loader"

def initialize(context):
    context.register_file_opener(".smiles", lambda path: load_smiles(path, context))

def load_smiles(path, context):
    try:
        with open(path, 'r') as f:
            smiles = f.read().strip()

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)

            # 1. Update the app
            context.current_molecule = mol

            # 2. Save to Undo history (MANDATORY for importers)
            context.push_undo_checkpoint()

            # 3. Adjust View
            context.reset_3d_camera()
            context.show_status_message(f"Imported {path} successfully.")
    except Exception as e:
        context.show_status_message(f"Import Failed: {e}")
```

### 10.5 Custom 3D Style (Dynamic Overrides)
Custom styles in V3 should combine standard drawing with custom elements.

```python
PLUGIN_NAME = "High-Contrast Style"

def initialize(context):
    context.register_3d_style("High Contrast", draw_hc_style)

def draw_hc_style(mw, mol):
    # 1. Get the plotter from the window passed by the engine
    plotter = mw.plotter

    # 2. Add a background gradient or text
    plotter.set_background("black", top="gray")
    plotter.add_text("HIGH CONTRAST MODE", position='upper_left', color='yellow')

    # 3. Trigger standard draw with a style override
    # Note: Accessing view_3d_manager on mw is allowed in style callbacks
    mw.view_3d_manager.draw_standard_3d_style(mol, style_override='stick')
```

### 10.6 Custom Optimization Method
Optimization callbacks must modify the molecule **in-place**.

```python
from rdkit.Chem import AllChem

PLUGIN_NAME = "Fast UFF Optimizer"

def initialize(context):
    context.register_optimization_method("Quick UFF", run_uff)

def run_uff(mol):
    """
    Called by the 'Optimize' logic.
    Modify 'mol' in-place. Return True on success.
    """
    if mol.GetNumAtoms() == 0:
        return False

    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return True
```

### 10.7 Persistent User Preferences
Use `get_setting` and `set_setting` to remember simple user choices across sessions. For preferences that must survive an application-wide settings reset or plugin updates, use the companion JSON pattern from [Section 2.7](#27-settings--persistence) instead.

```python
PLUGIN_NAME = "Smart Labels"

def initialize(context):
    # Load user preference (default to True)
    show_labels = context.get_setting("show_labels", True)
    
    context.add_menu_action("Labels/Toggle", lambda: toggle_labels(context))
    context.show_status_message(f"Labels are {'on' if show_labels else 'off'}")

def toggle_labels(context):
    current = context.get_setting("show_labels", True)
    new_state = not current
    context.set_setting("show_labels", new_state)
    context.show_status_message(f"Labels {'on' if new_state else 'off'}")
    context.refresh_3d_view()
```

---

## 11. UI & UX Style Guide

To ensure your plugin feels like a native part of MoleditPy, follow these design principles:

### 11.1 Color Palette
MoleditPy uses a "Sleek Dark" theme. If you build custom widgets, use these QColor hints:
- **Background**: `QColor(30, 30, 35)`
- **Accent**: `QColor(0, 120, 215)`
- **Secondary Text**: `QColor(180, 180, 180)`

### 11.2 Icons
Prefer using the system's built-in icons where possible. If you provide your own, use **24x24 pixel PNGs** with transparency for toolbars.

---

## 12. Testing Your Plugins

You can test your plugin logic without running the full MoleditPy GUI by mocking the `PluginContext`.

```python
import unittest
from unittest.mock import MagicMock
from rdkit import Chem

class TestMyPlugin(unittest.TestCase):
    def test_logic(self):
        # 1. Create a mock context
        context = MagicMock()
        context.current_molecule = Chem.MolFromSmiles("CCO")

        # 2. Run your logic
        my_plugin_function(context)

        # 3. Verify side effects
        context.push_undo_checkpoint.assert_called_once()
```

---

## 13. Conclusion & Support

Thank you for contributing to the MoleditPy ecosystem! If you encounter any bugs or need new API features, please reach out via GitHub Issues.

**Happy Coding!**
