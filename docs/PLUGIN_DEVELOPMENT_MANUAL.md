# MoleditPy Plugin Development Manual

Welcome to the **MoleditPy Plugin Development Manual**! MoleditPy features a robust plugin system that allows you to extend the application with new tools, menus, visualizations, and file formats.

## 1. Introduction

MoleditPy plugins are Python scripts (`.py`) placed in your user plugin directory:
- **Windows**: `C:\Users\<YourName>\.moleditpy\plugins\`
- **Linux/macOS**: `~/.moleditpy/plugins/`

The application automatically discovers and loads these plugins on startup.

## 2. Plugin Structure

A modern MoleditPy plugin consists of **Metadata** and an **Initialization Hook**.

### Metadata Variables
Define these at the top level of your script:
```python
PLUGIN_NAME = "My Super Plugin"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "Jane Doe"
PLUGIN_DESCRIPTION = "Adds super powers to MoleditPy."
```

### The `initialize` Function
This is the entry point. It receives a `context` object (of type `PluginContext`) that you use to register your extensions.

```python
def initialize(context):
    # Register your hooks here
    context.add_menu_action("My Plugin/Action", my_callback)
```

## 2.1 Folder Plugins (Packages)

For larger plugins, you can organize your code into a folder (a Python package) instead of a single `.py` file.

**Structure:**
```text
plugins/
  └── MyComplexPlugin/      <-- Folder Name
       ├── __init__.py      <-- Entry point & Metadata (REQUIRED)
       ├── utils.py         <-- Other modules
       └── assets/          <-- Images, data files
            └── icon.png
```

**Requirements:**
1.  **`__init__.py` is mandatory**: The folder MUST contain an `__init__.py` file to be recognized as a single plugin.
2.  **Metadata & Entry Point**: Place your `PLUGIN_NAME`, `PLUGIN_VERSION`, and the `initialize(context)` function inside `__init__.py`.
3.  **Relative Imports**: You can import other modules in your folder using standard relative imports (e.g., `from .utils import my_helper`).

**Example `__init__.py`:**
```python
# __init__.py
from .utils import helper_function

PLUGIN_NAME = "My Complex Plugin"
PLUGIN_VERSION = "1.0"

def initialize(context):
    helper_function(context)
```

> [!NOTE]
> **Update Behavior:** The **Plugin Installer** plugin completely replaces the plugin folder when updating a folder-based plugin. It specifically preserves **only** `settings.json` (by backing it up before the update and restoring it afterwards). Any other user data or configuration files stored inside the plugin folder will be lost.

---

## 3. The PluginContext API

The `context` object (Type: `PluginContext`) passed to your `initialize` function provides the following methods for interacting with the application.

### Quick Summary

| Category | Method | Description |
| :--- | :--- | :--- |
| **UI** | `add_menu_action` | Add item to main menu. |
| **UI** | `add_toolbar_action` | Add button to plugin toolbar. |
| **UI** | `add_analysis_tool` | Add tool to Analysis menu. |
| **UI** | `add_export_action` | Add option to Export menu. |
| **UI** | `show_status_message`| Show feedback in the status bar. |
| **Files** | `register_file_opener` | Handle specific file extensions. |
| **Files** | `register_drop_handler` | Handle drag-and-drop files. |
| **State** | `register_save_handler` | Save data to `.pmeprj`. |
| **State** | `register_load_handler` | Load data from `.pmeprj`. |
| **State** | `register_document_reset_handler` | Reset plugin state on File→New. |
| **State** | `push_undo_checkpoint` | Save current state to undo history. |
| **3D** | `register_optimization_method`| Add geometry optimization method. |
| **3D** | `register_3d_style` | Add custom 3D visualization style. |
| **3D** | `refresh_3d_view` | Redraw the 3D scene. |
| **3D** | `reset_3d_camera` | Fit camera to current molecule. |
| **Windows** | `register_window` | Store and manage plugin windows. |
| **Windows** | `get_window` | Retrieve a registered window. |
| **Access** | `current_molecule` | Get/Set the active RDKit molecule. |
| **Access** | `get_selected_atom_indices` | Get indices of selected atoms. |
| **Access** | `get_3d_controller` | Get helper for 3D manipulation. |
| **Access** | `get_main_window` | Get raw MainWindow instance. |


### 3.1 UI & Menus

#### `add_menu_action(path, callback, text=None, icon=None, shortcut=None)`
Register a custom item in the main menu.

- **path** (`str`): The full menu path. 
    - Use `File/My Action` to add to existing menus.
    - Use `MyPlugin/Action` to create a new top-level menu.
- **callback** (`Callable`): Function to execute when the action is triggered.
- **text** (`str`, optional): The label of the menu item. Defaults to the last part of `path`.
- **icon** (`str`, optional): specialized icon name or path.
- **shortcut** (`str`, optional): Keyboard shortcut (e.g., "Ctrl+Shift+X").

**Example Usage:**
```python
context.add_menu_action("Edit/My Action", my_func, shortcut="Ctrl+M")
```

#### `add_toolbar_action(callback, text, icon=None, tooltip=None)`
Add a button to the dedicated **Plugin Toolbar**. This toolbar is only visible if at least one plugin registers an action.

- **callback** (`Callable`): Function to execute when clicked.
- **text** (`str`): Label for the button.
- **icon** (`str`, optional): Icon name or path.
- **tooltip** (`str`, optional): Hover text.

**Example Usage:**
```python
context.add_toolbar_action(my_func, "My Tool", icon="icon.png", tooltip="Run My Tool")
```

#### `add_analysis_tool(label, callback)`
Register a tool in the top-level **Analysis** menu.

- **label** (`str`): Text to display in the menu.
- **callback** (`Callable`): Function to execute.

**Example Usage:**
```python
context.add_analysis_tool("Calculate Mass", run_mass_calc)
```

#### `add_export_action(label, callback)`
Register an action in the **Export** menu (accessed via File > Export or the status bar export button).

- **label** (`str`): Text to display (e.g., "Export as MyFormat...").
- **callback** (`Callable`): Function to execute.

**Example Usage:**
```python
context.add_export_action("Export as JSON...", export_json_func)
```

#### `show_status_message(message, timeout=3000)`
Display a temporary message in the application's bottom status bar.

- **message** (`str`): The text to display.
- **timeout** (`int`): Duration in milliseconds before the message disappears.

**Example Usage:**
```python
context.show_status_message("Calculation complete!", 2000)
```

### 3.2 File Handling & Project State

#### `register_file_opener(extension, callback, priority=0)`
Register a handler for opening a specific file type. This handler is used for both the **Import** menu and **Command Line** file opening (e.g., `moleditpy myfile.ext`).

- **extension** (`str`): File extension (e.g., `.xyz`, `.cub`).
- **callback** (`Callable[[str], None]`): Function that receives the absolute file path and handles loading.
- **priority** (`int`): Handlers with higher priority are checked first. Default is 0. You can use **negative values** (e.g., -1) to register a fallback handler that only runs if no other plugin handles the file.

**Example Usage:**
```python
# Standard registration
context.register_file_opener(".mysim", open_simulation_file)

# High priority opener (overrides default or other plugins)
context.register_file_opener(".common", my_opener, priority=100)

# Fallback opener (runs only if nobody else claims it)
context.register_file_opener(".common", fallback_opener, priority=-1)
```

#### `register_drop_handler(callback, priority=0)`
Register a handler for files dropped onto the main window. **Note**: This is distinct from `register_file_opener`; handlers here are ONLY triggered by drag-and-drop operations, not by the Import menu or Command Line.

- **callback** (`Callable[[str], bool]`): Function that receives the dropped file path. Must return `True` if it handled the file, `False` otherwise.
- **priority** (`int`): Handlers with higher priority are checked first. Default is 0. Negative values (e.g., -1) can be used for fallback handlers.

**Example Usage:**
```python
def handle_drop(path):
    if path.endswith(".log"):
        parse_log(path)
        return True
    return False

# High priority
context.register_drop_handler(handle_drop, priority=10)

# Fallback (runs last)
context.register_drop_handler(my_fallback_func, priority=-1)
```

#### `register_save_handler(callback)`
Register a callback to save custom state into the application's project file (`.pmeprj`).

- **callback** (`Callable[[], dict]`): Function that must return a dictionary of JSON-serializable data.

**Example Usage:**
```python
context.register_save_handler(lambda: {"my_plugin_version": "1.0", "active": True})
```

#### `register_load_handler(callback)`
Register a callback to restore custom state from the project file.

- **callback** (`Callable[[dict], None]`): Function that receives the dictionary previously saved by your save handler.

**Example Usage:**
```python
def restore_state(data):
    if "my_plugin_version" in data:
        print(f"Restored plugin state: {data['active']}")
context.register_load_handler(restore_state)
```

#### `register_document_reset_handler(callback)`
Register a callback to be invoked when a new document is created (File→New). Use this to reset your plugin's internal state when the user clears all data.

- **callback** (`Callable[[], None]`): Function with no arguments that resets plugin state.

**Example Usage:**
```python
# Plugin state
calculation_cache = {}
last_result = None

def on_document_reset():
    """Reset plugin data when user creates a new document."""
    global calculation_cache, last_result
    calculation_cache.clear()
    last_result = None
    print("Plugin state has been reset")

context.register_document_reset_handler(on_document_reset)
```

**Use Cases:**
- Clear cached calculation results
- Reset UI state (close dialogs, clear selections)
- Release resources associated with the previous document
- Reset counters or temporary data

**Note:** This handler is called AFTER the main application has cleared all molecular data but BEFORE the "Cleared all data" status message is displayed.

#### `push_undo_checkpoint()`
Snapshots the current application state and adds it to the Undo history.

> [!IMPORTANT]
> **Usage Timing**: You should call this method **AFTER** you have finished modifying the molecule. The system only pushes a new state if it detects a difference from the previous one.

**Example Usage:**
```python
# 1. Modify the molecule
context.current_molecule = my_new_mol
# 2. Push checkpoint so user can Undo the change
context.push_undo_checkpoint()
```

### 3.3 Computation & 3D Visualization

#### `register_optimization_method(method_name, callback)`
Add a new method to the **Compute > Optimize Geometry** menu.

- **method_name** (`str`): Name of the method (e.g., "My Forcefield").
- **callback** (`Callable[[rdkit.Chem.Mol], bool]`): Function that receives the RDKit molecule object. It should modify the molecule's geometry in-place and return `True` on success.

**Example Usage:**
(See **Example 6** for a full implementation)
```python
context.register_optimization_method("My Opt", my_optimization_func)
```

#### `register_3d_style(style_name, callback)`
Register a fully custom 3D visualization mode.

- **style_name** (`str`): Unique name for the style.
- **callback** (`Callable[[MainWindow, rdkit.Chem.Mol], None]`): Function responsible for the entire drawing process.

**Example Usage:**
```python
def draw_red_spheres(mw, mol):
    mw.plotter.add_mesh(pv.Sphere(radius=1.0), color='red')
context.register_3d_style("Red Spheres", draw_red_spheres)
```

#### `refresh_3d_view()`
Force the 3D scene to re-render. Call this after changing custom meshes or manually updating atom colors via the controller.

#### `reset_3d_camera()`
Resets the 3D camera to fit the current molecule in the view.

#### `get_selected_atom_indices()`
Returns a list of RDKit atom indices that are currently selected in either the 2D editor or the 3D viewer.

**Example Usage:**
```python
selected = context.get_selected_atom_indices()
if selected:
    context.show_status_message(f"Working on {len(selected)} atoms")
```

#### `get_3d_controller()`
Returns a `Plugin3DController` instance for simplified manipulation of the 3D scene.

**Example Usage:**
```python
ctrl = context.get_3d_controller()
ctrl.set_atom_color(0, "#00FF00") # Set first atom to Green
```

#### `get_main_window()`
Returns the raw `MainWindow` instance. This gives you unrestricted access to the application internals (PyQt widgets, RDKit mol, PyVista plotter, etc.).

**Example Usage:**
```python
mw = context.get_main_window()
mw.statusBar().showMessage("Plugin accessed Main Window!")
plotter = mw.plotter
mw.statusBar().showMessage("Plugin accessed Main Window!")
plotter = mw.plotter
scene = mw.scene
```

#### `current_molecule`
Property to get or set the active RDKit molecule (`mw.current_mol`).

- **Get**: Returns the current `rdkit.Chem.Mol` object (or `None`).
- **Set**: Updates the main window's molecule and triggers a 3D redraw.

**Example Usage:**
```python
# Get
mol = context.current_molecule
if mol:
    print(mol.GetNumAtoms())

# Set
from rdkit import Chem
context.current_molecule = Chem.MolFromSmiles("C")
```

### 3.4 Window & Dialog Management

Plugins can register their custom windows with the application to manage their lifecycle and avoid garbage collection issues.

#### `register_window(window_id, window)`
Store a reference to a window or dialog. 

> [!TIP]
> **Namespacing**: Window IDs are automatically namespaced by your plugin name. If your plugin is "MyTool" and you register `"main"`, its internal ID becomes `MyTool:main`. This prevents collisions with other plugins.

#### `get_window(window_id)`
Retrieve a previously registered window by its local ID. Returns `None` if not found.

**Example Usage (Singleton Pattern):**
```python
def show_my_tool(context):
    win = context.get_window("my_dialog")
    if win:
        win.show()
        win.raise_()
    else:
        win = MyPluginDialog(context.get_main_window())
        context.register_window("my_dialog", win)
        win.show()
```


---

## 3b. Helper Classes

### Plugin3DController

Obtained via `context.get_3d_controller()`.

#### `set_atom_color(atom_index, color_hex)`
Override the color of a specific atom in the 3D view.

- **atom_index** (`int`): The 0-based index of the atom.
- **color_hex** (`str`): Hex color code (e.g., "#FF0000").

#### `set_bond_color(bond_index, color_hex)`
Override the color of a specific bond.

- **bond_index** (`int`): The 0-based index of the bond.
- **color_hex** (`str`): Hex color code.

#### `set_bond_color_by_atoms(atom_idx1, atom_idx2, color_hex)`
Override the color of the bond defined by two atom indices. This is often more convenient than finding the bond index manually.

- **atom_idx1** (`int`): RDKit index of the first atom.
- **atom_idx2** (`int`): RDKit index of the second atom.
- **color_hex** (`str`): Hex color code.


## 4. Advanced: MainWindow Internals

The `context.get_main_window()` method returns the raw `MainWindow` instance (Type: `PyQt6.QtWidgets.QMainWindow`). This object holds the entire application state.

### Key Attributes

| Attribute | Type | Description |
| :--- | :--- | :--- |
| `mw.current_mol` | `rdkit.Chem.Mol` | The active molecule object. |
| `mw.plotter` | `pyvista.Plotter` | The 3D view controller. Use to add meshes, actors, or change camera. |
| `mw.scene` | `MoleculeScene` | The 2D editor scene (inherits `QGraphicsScene`). |
| `mw.settings` | `dict` | Global application settings (e.g., colors, optimizer config). |
| `mw.splitter` | `QSplitter` | The container dividing 2D (index 0) and 3D (index 1) views. |

### Useful Internal Methods

| Method | Description |
| :--- | :--- |
| `mw.push_undo_state()` | Snapshots the current molecule and adds it to the Undo stack. Call this **after** modifying `mw.current_mol`. |
| `mw.draw_molecule_3d(mol)` | Renders the given molecule in the 3D window. |
| `mw.set_mode(mode_str)` | Switches the 2D editor mode (e.g., `'select'`, `'atom_C'`, `'bond_1'`). |
| `mw.statusBar().showMessage(msg)` | Displays text in the bottom status bar. |

### Common Tasks

| Task | Code Snippet |
| :--- | :--- |
| **Load MOL File** | `mw.load_mol_file("path/to/file.mol")` |
| **Load XYZ (view only)** | `mw.load_xyz_for_3d_viewing("path/to/file.xyz")` |
| **Clear 3D View** | `mw.plotter.clear()` |
| **Update 3D Style** | `mw.main_window_view_3d.update_3d_style()` |
| **Center Camera** | `mw.plotter.reset_camera()` |

> [!TIP]
> **Helper Access**: Many core functions are split into helper modules. You can often access their methods directly via `mw.main_window_helper_name.method()`, although direct calls are safer where `mw` exposes them (like `mw.load_mol_file`).

---

## 5. Debugging & Best Practices

### Thread Safety Warning
All plugin callbacks run in the **Main UI Thread**. 
- **Fast operations** (< 100ms) are fine.
- **Slow operations** (Heavy QM calculations, large loops) will **freeze the entire application**.
- **Solution**: Use Python's `threading` module or `QThread` for heavy work.

### Error Handling
If your plugin callback raises an unhandled exception, it may crash the application or simply do nothing (silently fail), depending on where it hook is.
**Recommendation**: Wrap your callback logic in a `try...except` block and use `QMessageBox` to show errors to the user.

```python
def my_callback():
    try:
        # Dangerous code
        risky_operation()
    except Exception as e:
        from PyQt6.QtWidgets import QMessageBox
        # Safe access to MW might fail if MW is not stored, so catch broadly
        print(f"Plugin Error: {e}")
```

### Plugin Naming & Menu Grouping
- **Unique Names**: Ensure your `PLUGIN_NAME` is unique. Plugins with identical names may be grouped together in the UI, which can lead to confusion.
- **Menu Sorting**: The **Import** menu is sorted alphabetically by Plugin Name.

### File Extension Conflicts
- **Import Menu**: If multiple plugins accept the same extension, they will all appear in the Import menu as separate entries (e.g., "Import .xyz (Plugin A)" and "Import .xyz (Plugin B)").
- **Command Line / Open Command**: If the user opens a file directly (e.g. via CLI), the plugin with the **highest priority** wins. If priorities are equal, the behavior is undefined (first registered wins).

### Hot Reloading
MoleditPy loads plugins at startup. To see changes in your code, you must currently **restart the application** or use the **"Reload Plugins"** feature if available (Settings dependent).

---

## 6. Legacy Support (Pre-2.2)

Older plugins used a different structure which is still supported but less capable.

### The `run` Function
If you define `run(main_window)`, your plugin will automatically appear in the "Plugins" menu. When clicked, `run` is executed with the raw `MainWindow` object.

```python
# Legacy Style
def run(main_window):
    mol = main_window.current_mol
    # ... direct access to implementation details ...
```

### The `autorun` Function
Legacy plugins often used `autorun(main_window)` to execute code immediately on startup. Use `initialize(context)` for this purpose in modern plugins.

```python
def autorun(main_window):
    print("I run immediately!")
```

**Note**: You can mix both! Define `initialize(context)` for startup hooks and `run(main_window)` to show a manual action in the main "Plugins" menu.

---

## 7. Full Example

```python
import os
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Complete Example"
PLUGIN_VERSION = "2.0"
PLUGIN_AUTHOR = "MoleditPy Team"

def initialize(context):
    print("Initializing Example Plugin...")

    # 1. Add a menu item
    def show_info():
        # You can access the main window via context if needed
        mw = context.get_main_window()
        QMessageBox.information(mw, "Info", "Plugin is active!")
    
    context.add_menu_action("Example Plugin/Show Info", show_info)

    # 2. Add an analysis tool
    context.add_analysis_tool("Count Atoms", lambda: print("Counting..."))

    # 3. Handle a custom file type
    context.register_file_opener(".ex", lambda path: print(f"Opening {path}"))

# Optional: Add a 'Run' entry to the Plugins menu
def run(mw):
    QMessageBox.information(mw, "Manual Run", "You clicked me in the Plugins menu!")
```

## 8. Cookbook / Examples

Here are common example references, updated for the modern API.
Note: For advanced logic involving the 2D scene or 3D plotter, you will often use `context.get_main_window()` (`mw`) to access `mw.scene`, `mw.data`, or `mw.plotter`.

### Example 1: Analysis Tool (Molecular Weight)
Add a tool to the **Analysis** menu that calculates properties.

```python
from rdkit.Chem import Descriptors
from PyQt6.QtWidgets import QMessageBox

PLUGIN_NAME = "Calculate Stats"

def initialize(context):
    def run_calc():
        # Access the main application to get the current molecule
        mw = context.get_main_window()
        mol = mw.current_mol
        
        if not mol:
            QMessageBox.warning(mw, "Error", "No molecule loaded!")
            return

        weight = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        
        QMessageBox.information(mw, "Stats", f"Molecular Weight: {weight:.2f}\nAtoms: {num_atoms}")

    context.add_analysis_tool("Show Molecular Weight", run_calc)
```

### Example 2: Custom 3D Visualization (Add Sphere)
Use the `plotter` object (PyVista) to add custom 3D geometries.

```python
import pyvista as pv

PLUGIN_NAME = "Add 3D Sphere"

def initialize(context):
    context.add_menu_action("Visuals/Add Sphere", lambda: add_sphere(context))

def add_sphere(context):
    mw = context.get_main_window()
    plotter = mw.plotter
    
    # Create and add a sphere
    sphere = pv.Sphere(radius=2.0)
    plotter.add_mesh(sphere, color="red", opacity=0.3, name="custom_sphere")
    plotter.render()
```

### Example 3: Custom Dock Panel (Singleton Pattern)
Add a permanent side panel UI using the modern window registration API to prevent duplicates and ensure cleanup.

```python
from PyQt6.QtWidgets import QDockWidget, QLabel, QVBoxLayout, QWidget, QPushButton
from PyQt6.QtCore import Qt

PLUGIN_NAME = "My Panel"

def initialize(context):
    context.add_menu_action("Tools/Show My Panel", lambda: toggle_panel(context))

def toggle_panel(context):
    mw = context.get_main_window()
    
    # Check if already registered
    dock = context.get_window("main_dock")
    if dock:
        dock.show()
        return

    # Create new dock
    dock = QDockWidget("My Tools", mw)
    content = QWidget()
    layout = QVBoxLayout(content)
    layout.addWidget(QLabel("Hello from Plugin!"))
    dock.setWidget(content)
    
    mw.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    
    # Register to manage lifecycle
    context.register_window("main_dock", dock)
    dock.show()
```

### Example 4: Custom File Importer
Register a handler for a custom file format (e.g. `.simple`).

```python
PLUGIN_NAME = "Simple Importer"

def initialize(context):
    context.register_file_opener(".simple", lambda path: load_simple(path, context))

def load_simple(path, context):
    mw = context.get_main_window()
    print(f"Parsing {path}...")
    
    # 1. Parse your file -> create RDKit molecule
    # (Here we create a dummy molecule for demonstration)
    from rdkit import Chem
    mol = Chem.MolFromSmiles("C1CCCCC1") # Cyclohexane

    if mol:
        # 2. Update the main window's current molecule
        # You can use the convenience property:
        context.current_molecule = mol
        
        # 3. Commit this state to the Undo Stack mechanism
        # (This ensures the user can Undo this import)
        # Note: 'context.current_molecule =' triggers a redraw but does not 
        # auto-push to Undo stack, so we still access mw for that.
        mw.push_undo_state()
        
        # 4. Update the 3D View (handled by setter, but camera reset is manual)
        mw.plotter.reset_camera()

        mw.plotter.reset_camera()
        
        # 5. Optional: Update UI mode (if switching from 2D)
        # mw._enter_3d_viewer_ui_mode()

```

### Example 5: Custom 3D Style (Native Registration)
Register a new visualization mode that appears in the 3D Style menu.

```python
import pyvista as pv

PLUGIN_NAME = "My Style Plugin"

def draw_custom_style(mw, mol):
    # 1. Base drawing (optional helper to draw standard atoms)
    # mw.main_window_view_3d.draw_standard_3d_style(mol, style_override='ball_and_stick')
    
    # 2. Add custom visualization
    mw.plotter.add_text("Custom Style Active", position='upper_left')
    mw.plotter.add_mesh(pv.Sphere(radius=5), color='blue', opacity=0.2)

def initialize(context):
    context.register_3d_style("My Blue Sphere", draw_custom_style)
```

### Example 6: Custom Optimization Method
Register an optimization algorithm for the 3D structure. The callback receives the RDKit molecule object and must modify it **in-place**.

> [!IMPORTANT]
> The optimization callback runs in the **Main UI Thread**. Long-running calculations will freeze the application interface.
> - For quick optimizations (e.g., standard RDKit force fields), this is acceptable.
> - For complex/slow calculations, you should manage your own threading or use `QApplication.processEvents()` with caution.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

PLUGIN_NAME = "Naive Optimization"

def optimize_naive(mol):
    """
    Callback for 3D optimization.
    
    Args:
        mol (rdkit.Chem.Mol): The molecule object to optimize.
                              You must modify this object IN-PLACE.
                              Ensure it has a conformer (mol.GetNumConformers() > 0).
    
    Returns:
        bool: True if optimization succeeded, False otherwise.
              Returning False will display a "Optimization method returned failure" message.
    """
    try:
        # Check if molecule has atoms
        if mol.GetNumAtoms() == 0:
            return False
            
        print("Running Custom Optimization...")

        # 1. Modify the molecule structure in-place
        # Example: Simple UFF optimization using RDKit
        # (maxIters=200 is small enough to run in the main thread without major freezing)
        res = AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        
        # 2. Return success status
        # AllChem.UFFOptimizeMolecule returns 0 if converged, 1 if more iterations needed.
        # Both can be considered "success" in that the coordinates were updated.
        return True # Generally return True unless a critical error occurred

    except Exception as e:
        print(f"Optimization failed: {e}")
        return False

def initialize(context):
    # Register the method. 
    # It will appear in the Right-Click Context Menu of the "Optimize 3D" button
    # in the right-hand panel.
    # Note: Plugin methods cannot currently be set as the default global optimizer 
    # in "Settings" > "3D Optimization Settings".
    context.register_optimization_method("Naive UFF", optimize_naive)
```





