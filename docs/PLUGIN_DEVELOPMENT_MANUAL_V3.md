# MoleditPy Plugin Development Manual (Version 3.0)

Welcome to the **Version 3.0** of the MoleditPy Plugin API. This version introduces a fully decoupled, namespaced architecture designed for high stability, clean memory management, and long-term maintainability.

> [!IMPORTANT]
> **API PHILOSOPHY**: In V3, you MUST avoid accessing the `MainWindow` directly via monkey-patching or unverified attributes. Instead, use the methods provided by the `PluginContext`. This ensures your plugin remains functional even if the main application's internal structure is refactored.

---

## 1. Getting Started

A MoleditPy plugin is a Python script (or folder package) that defines metadata and an `initialize(context)` entry point.

### 1.1 Plugin Metadata
Define these at the top of your script. They are used for the UI and the internal registry.

| Variable | Description |
| :--- | :--- |
| `PLUGIN_NAME` | **(Required)** The human-readable name of the plugin. |
| `PLUGIN_VERSION` | Version string (e.g., `"1.0.2"` or `"2026.03.31"`). |
| `PLUGIN_AUTHOR` | Name of the developer. |
| `PLUGIN_DESCRIPTION` | A short summary shown in the Plugin Manager. |

### 1.2 Folder-based Plugins (Packages)
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

```python
PLUGIN_NAME = "System Diagnostic Tool"
PLUGIN_VERSION = "3.0.1"
PLUGIN_AUTHOR = "Research Team"
PLUGIN_DESCRIPTION = "Analyzes molecular symmetry and logs diagnostics."

def initialize(context):
    """
    Entry point called by the PluginManager during application startup.
    :param context: PluginContext - Your interface to the application.
    """
    # Add a menu item
    context.add_menu_action("Tools/Run Diagnostic", lambda: run_diag(context))

    # Add a toolbar button
    context.add_toolbar_action(
        callback=lambda: run_diag(context),
        text="Diag",
        tooltip="Run System Diagnostics"
    )
```

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

#### `register_menu_action(path, text_or_callback, callback=None, icon=None, shortcut=None)`
Backward-compatible alias for `add_menu_action`. Supports two calling styles:
- **New style** (preferred): `register_menu_action(path, callback)` — equivalent to `add_menu_action`.
- **Old style** (V2 compat): `register_menu_action(path, text, callback)` — text and callback positions are swapped.

> [!NOTE]
> New plugins should use `add_menu_action` instead of `register_menu_action`.

#### `add_toolbar_action(callback, text, icon=None, tooltip=None)`
Add a button to the dedicated **Plugin Toolbar**.
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
Register a handler for opening a specific file type. Used for the **File > Import** menu and **Command Line** startup.
- **extension** (`str`): File extension including the dot (e.g., `".xyz"`, `".cub"`).
- **callback** (`Callable[[str], None]`): Function that receives the absolute file path and handles the loading logic.
- **priority** (`int`): Handlers with higher priority are checked first. Use **negative values** (e.g., -1) to register a fallback handler.

#### `register_drop_handler(callback, priority=0)`
Register a handler for files dropped onto the main 2D/3D editor window.
- **callback** (`Callable[[str], bool]`): Function that receives the dropped file path. Must return `True` if it successfully handled the file, `False` otherwise.
- **priority** (`int`): Handlers with higher priority are checked first.

---

### 2.3 Molecular State & Undo

Methods for interacting with the active molecule and managing the undo stack.

#### `current_molecule` (Property)
Get or set the active RDKit molecule object.
- **Getter**: Returns the `rdkit.Chem.Mol` currently loaded in the editor.
- **Setter**: Replaces the active molecule. Automatically triggers 2D/3D redrawing of the scene.

> [!NOTE]
> `current_mol` is an alias for `current_molecule` — both are available.

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

#### `scene` (Property)
Direct access to the 2D `MoleculeScene`.

#### `get_main_window()`
Returns the raw `MainWindow` instance. **Use with caution** — prefer specific context methods when available.

---

### 2.6 Window & State Management (Namespaced)

In V2, plugins often attached windows directly to `mw`. In V3, use the **Namespaced Registry**. This keeps your windows alive in memory and prevents ID collisions with other plugins.

| Method | Description |
| :--- | :--- |
| `register_window(id, window)` | Stores a Qt window/dialog. The ID is automatically prefixed with your `PLUGIN_NAME`. |
| `get_window(id)` | Retrieves your registered window by ID. Returns `None` if not found or if the window was deleted. |

---

## 3. Visual Overrides: Plugin3DController

Obtain a controller via `ctrl = context.get_3d_controller()`. These overrides are temporary visual changes that do not modify the actual RDKit molecule coordinates or properties.

- `set_atom_color(atom_index: int, color_hex: str)`: Override the color of a specific atom (e.g., `"#FF0000"` for red).
- `set_bond_color(bond_index: int, color_hex: str)`: Override the color of a specific RDKit bond by its index.
- `set_bond_color_by_atoms(idx1: int, idx2: int, color_hex: str)`: Helper to find and color the bond connecting two atom indices.

---

## 4. Modern Workflow Examples

### 4.1 The "Singleton Dialog" Pattern
Always check if your window is already open before creating a new one.

```python
def toggle_viewer(context):
    # 1. Try to get existing window
    win = context.get_window("main_panel")

    if win:
        win.show()
        win.raise_()
        return

    # 2. Create and Register if it's the first time
    win = MyCustomDialog(context.get_main_window())
    context.register_window("main_panel", win)
    win.show()
```

### 4.2 The "Modify-Then-Push" Undo Pattern
Correct ordering is essential for the Undo/Redo system to function.

```python
def center_molecule(context):
    mol = context.current_molecule
    if not mol: return

    # 1. Modify the molecule coordinates
    AllChem.ComputeCanonicalTransform(mol)

    # 2. Update the application state
    context.current_molecule = mol

    # 3. Save to History (Must come after setting the mol)
    context.push_undo_checkpoint()

    # 4. Success feedback
    context.show_status_message("Molecule centered.", 2000)
```

### 4.3 Highlighting Selection
Color atoms that the user has selected.

```python
def highlight_selection(context):
    indices = context.get_selected_atom_indices()
    ctrl = context.get_3d_controller()

    for idx in indices:
        ctrl.set_atom_color(idx, "#00FF00")  # Highlight green

    context.refresh_3d_view()
```

---

## 5. Migration Guide (V2 to V3)

| Legacy Pattern (Avoid) | Modern Pattern (Use This) |
| :--- | :--- |
| `mw.statusBar().showMessage(m)` | `context.show_status_message(m)` |
| `mw.push_undo_state()` | `context.push_undo_checkpoint()` |
| `mw.current_mol` | `context.current_molecule` |
| `mw.plotter.reset_camera()` | `context.reset_3d_camera()` |
| `mw.my_tool = win` | `context.register_window("my_tool", win)` |
| `mw.draw_molecule_3d(...)` | `context.refresh_3d_view()` |
| `context.register_menu_action(path, text, cb)` | `context.add_menu_action(path, cb, text)` |

> [!NOTE]
> `mw.trigger_conversion()` is still accessible via `context.get_main_window()` and delegates to the internal compute manager. Prefer `context.add_menu_action` or `context.add_analysis_tool` for triggering computation workflows from menus.

---

## 6. Advanced: MainWindow Internals (Low-Level)

While V3 encourages using `context`, sometimes you need direct access to the `MainWindow` (`mw`) for specialized Qt or PyVista operations. Obtain it via `mw = context.get_main_window()`.

| Attribute / Method | Type | Description |
| :--- | :--- | :--- |
| `mw.settings` | `dict` | Global application settings (colors, defaults). |
| `mw.plotter` | `pv.Plotter` | The raw PyVista plotter instance. |
| `mw.scene` | `MoleculeScene` | The 2D editor graphics scene. |
| `mw.current_mol` | `rdkit.Chem.Mol` | The active molecule (same as `context.current_molecule`). |
| `mw.splitter` | `QSplitter` | The UI divider between 2D and 3D views. |
| `mw.load_mol_file(path)` | — | Standard file loader. |
| `mw.set_mode(mode_str)` | — | Switch editor mode (e.g., `'select'`, `'erase'`, `'atom_C'`). |
| `mw.trigger_conversion()` | — | Trigger the built-in structure conversion (delegates to `compute_manager`). |

---

## 7. Legacy Support (Compatibility Mode)

MoleditPy 3.0 remains compatible with older plugin patterns, though they are considered deprecated.

### 7.1 The `run(mw)` Function
If you define a function `run(main_window)`, MoleditPy will add your plugin to the standard "Plugins" menu. This is useful for simple scripts that don't need the advanced V3 ecosystem.

### 7.2 The `autorun(mw)` Function
Will be executed immediately on startup. Use `initialize(context)` for modern equivalents.

---

## 8. Debugging & Advanced Logistics

### 8.1 File Extension Conflicts
- **Import Menu**: If multiple plugins accept the same extension, they will appear in the Import menu as separate entries (e.g., "Import .xyz (Plugin A)" and "Import .xyz (Plugin B)").
- **Command Line**: If the user opens a file via the command line, the plugin with the **highest priority** wins. If priorities are equal, the first registered plugin is used.

### 8.2 Hot Reloading
MoleditPy scans for plugins at startup. To see code changes, simply use the **Plugins > Reload All Plugins** menu item. V3 plugins are designed to survive reloads if their windows are correctly registered via `context.register_window()`.

### 8.3 Threading & Freezing
All plugin callbacks run in the Main UI Thread. Performance-heavy tasks will freeze the GUI. Use `QThread` for calculations and `context.show_status_message` for progress updates.

---

## 9. Best Practices & Troubleshooting

### 9.1 Thread Safety
**CRITICAL**: All plugin callbacks run in the Main UI Thread.
- **Fast operations** (< 100ms) are fine.
- **Slow operations** (heavy QM, large loops) will **freeze the entire application**.
- **Solution**: Use `QThread` or Python's `threading` for heavy work and return results via signals.

### 9.2 Error Handling
Always wrap your logic in `try...except` blocks. An uncaught exception in a plugin can cause the entire application to hang or crash silently.
```python
def my_callback(context):
    try:
        # unsafe logic
    except Exception as e:
        context.show_status_message(f"Plugin Error: {e}", 5000)
```

### 9.3 Memory Safety
Use `context.register_window()` for all persistent UI elements (dialogs, panels). If you don't register them, Python's garbage collector might delete the window while it's still being used by the user.

### 9.4 Atom Indices
Remember that MoleditPy uses **0-based RDKit indices**. These may differ from "Atom Numbers" in some chemical file formats. Always verify your mappings.

---

## 10. Complete Integrated Example

The following example combines multiple V3 features: namespaced windows, undo history, and status bar updates.

```python
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Research Assistant"
PLUGIN_VERSION = "3.1.0"

def initialize(context):
    # 1. Add tool to menu
    context.add_menu_action("Research/Analyze Selection", lambda: process(context))

    # 2. Register for specific file types
    context.register_file_opener(".res", lambda path: load_research(path, context))

def process(context):
    indices = context.get_selected_atom_indices()
    if not indices:
        context.show_status_message("Select atoms first!", 2000)
        return

    context.show_status_message(f"Analyzing {len(indices)} atoms...")
    # Logic goes here...

def load_research(path, context):
    try:
        # Load logic...
        context.push_undo_checkpoint()
        context.show_status_message("Loaded research data.")
    except Exception as e:
        QMessageBox.critical(context.get_main_window(), "Error", str(e))
```

---

## 11. Cookbook & Examples

### 11.1 Analysis Tool: Molecular Weight
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

### 11.2 Custom 3D Component: Add Sphere
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

### 11.3 Persistent Dock Panel (Singleton)
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

### 11.4 Custom File Importer
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

### 11.5 Custom 3D Style (Dynamic Overrides)
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

### 11.6 Custom Optimization Method
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

---

## 12. UI & UX Style Guide

To ensure your plugin feels like a native part of MoleditPy, follow these design principles:

### 12.1 Color Palette
MoleditPy uses a "Sleek Dark" theme. If you build custom widgets, use these QColor hints:
- **Background**: `QColor(30, 30, 35)`
- **Accent**: `QColor(0, 120, 215)`
- **Secondary Text**: `QColor(180, 180, 180)`

### 12.2 Icons
Prefer using the system's built-in icons where possible. If you provide your own, use **24x24 pixel PNGs** with transparency for toolbars.

---

## 13. Testing Your Plugins

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

## 14. Distribution & Sharing

### 14.1 The `.mpp` Format (Future)
We are working on a unified plugin package format (`.mpp`). For now, simply zip your plugin folder and share it with other users.

### 14.2 Repository Registry
If you want your plugin to be listed in the official **Community Plugin Browser**, please submit a Pull Request to the main MoleditPy repository with your plugin added to the `contrib/` directory.

---

## 15. Conclusion & Support

Thank you for contributing to the MoleditPy ecosystem! If you encounter any bugs or need new API features, please reach out via GitHub Issues.

**Happy Coding!**
