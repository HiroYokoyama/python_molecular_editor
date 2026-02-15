# MoleditPy Architecture Documentation

## Overview
MoleditPy is a Python-based molecular editing software built with **PyQt6** for the GUI and **RDKit** for chemical informatics features (SMILES parsing, 3D conformation generation, etc.).

The application follows a modular design where the main window functionality is split across multiple mixin classes (traits pattern) to keep the `MainWindow` class manageable.

## Developer Guide / Feature Map

Use this map to find the code responsible for specific features.

| Feature Area | Key Modules / Classes | Primary Functions/Methods |
| :--- | :--- | :--- |
| **Startup / Init** | `main_window_main_init.py` | `init_ui`, `init_menu_bar` |
| **File I/O (Project)** | `main_window_project_io.py` | `save_project`, `load_raw_data` |
| **MOL/XYZ I/O** | `main_window_molecular_parsers.py` | `load_mol_file`, `save_as_mol`, `save_as_xyz` |
| **SMILES/InChI Import** | `main_window_string_importers.py` | `load_from_smiles`, `load_from_inchi` |
| **View Loading** | `main_window_view_loaders.py` | `load_xyz_for_3d_viewing`, `save_3d_as_mol` |
| **3D Rendering** | `main_window_view_3d.py` | `draw_standard_3d_style`, `fit_to_view` |
| **3D Optimization** | `main_window_compute.py`, `calculation_worker.py` | `optimize_3d_structure`, `run_optimization` |
| **2D Edit Operations** | `main_window_edit_actions.py` | `copy_selection`, `clean_up_2d_structure` |
| **3D Edit Dialogs** | `main_window_dialog_manager.py` | `open_planarize_dialog`, `open_alignment_dialog` |
| **3D Measurements** | `main_window_edit_3d.py` | `calculate_distance`, `toggle_measurement_mode` |
| **Media/Mesh Export** | `main_window_export.py` | `export_2d_png`, `export_3d_png` |
| **UI Management** | `main_window_ui_manager.py` | `set_mode`, `toggle_2d_panel` |
| **Dialog Handling** | `main_window_dialog_manager.py` | `open_settings_dialog`, `open_analysis_window` |
| **Undo/Redo** | `main_window_app_state.py` | `push_undo_state`, `undo`, `redo` |

## Architecture Breakdown

For in-depth details on specific subsystems, strictly refer to these dedicated documents:

| Subsystem | Document | Description |
| :--- | :--- | :--- |
| **Data & Scene** | [CORE_DATA_STRUCTURES.md](CORE_DATA_STRUCTURES.md) | Data model (`MolecularData`), Scene Graph (`MoleculeScene`), and Visual Items. |
| **Computations** | [ANALYSIS_AND_CALCULATIONS.md](ANALYSIS_AND_CALCULATIONS.md) | 3D generation (`CalculationWorker`), Optimization, and Property Analysis. |
| **UI Components** | [DIALOGS_AND_UI.md](DIALOGS_AND_UI.md) | Catalog of all Dialogs, Windows, and Custom Widgets. |
| **Plugins** | [PLUGIN_SYSTEM_INTERNAL.md](PLUGIN_SYSTEM_INTERNAL.md) | Internal architecture of the plugin system. |
| **Plugin Dev** | [Plugin Development Manual](https://github.com/HiroYokoyama/moleditpy-plugins/blob/main/PLUGIN_DEVELOPMENT_MANUAL.md) | External guide for creating plugins. |

### Complete Architecture Details

#### Core Components

##### Entry Point
- **`main.py`**: The application entry point. It handles:
    - initializing the `QApplication`.
    - Setting the AppUserModelID on Windows (for taskbar icon grouping).
    - Parsing command-line arguments (for opening files).
    - Instantiating and showing the `MainWindow`.

##### Data Model
The application separates the chemical data model from the visual representation.

###### `molecular_data.py` available as `MolecularData`
Holds the pure chemical data state.
- **Attributes**:
    - `atoms`: Dictionary mapping `atom_id` to atom data (symbol, position, charge, radical).
    - `bonds`: Dictionary mapping `(id1, id2)` tuples to bond data (order, stereo).
    - `adjacency_list`: Graph representation for efficient traversal.
- **Key Methods**:
    - `add_atom()`, `remove_atom()`: Manage atoms.
    - `add_bond()`, `remove_bond()`: Manage bonds.
    - `to_rdkit_mol()`: Converts the internal data model to an `RDKit` molecule object. This handles complex tasks like stereochemistry assignment (E/Z, R/S) and 3D coordinate generation.
    - `to_mol_block()`: Generates a V2000 MOL block string for saving.

##### Scene Graph (Visual Layer)
The editor uses `QGraphicsScene` for the 2D editing canvas.

###### `molecule_scene.py` available as `MoleculeScene`
Inherits from `QGraphicsScene`. It manages user interactions for editing molecules in 2D.
- **Modes**: The scene operates in various modes (e.g., `select`, `bond_1_0` (single bond), `atom_C` (Carbon atom), `charge_plus`).
- **Events**: Handles `mousePress`, `mouseMove`, `mouseRelease` for:
    - Drawing atoms and bonds.
    - Selecting items.
    - Dragging atoms (which updates connected bonds).
    - Context menus (right-click deletion).
- **Template System**: Logic for previewing and adding ring templates (Benzene, etc.).

###### `atom_item.py` available as `AtomItem`
Inherits from `QGraphicsItem`. Represents a single atom in the 2D scene.
- **Rendering**: Draws the atom symbol, implicit hydrogens, charge, and radical markers.
- **Styling**: Adapts to global settings (font size, font family, color).
- **Hit Detection**: `shape()` is defined for precise clicking.

###### `bond_item.py` available as `BondItem`
Inherits from `QGraphicsItem`. Represents a bond between two `AtomItem`s.
- **Connection**: Holds references to two `AtomItem` instances (`atom1`, `atom2`).
- **Rendering**: Draws lines based on bond order (Single, Double, Triple).
    - **Stereochemistry**: Draws Wedges and Dashes for 3D cues.
    - **Ring Visualization**: Detects if the bond is part of a ring (via RDKit analysis on the fly) to draw the "inside" double bond line correctly.
    - **E/Z Labels**: Renders "E" or "Z" labels for stereochemically defined double bonds.
    - **E/Z Labels**: Renders "E" or "Z" labels for stereochemically defined double bonds.

##### Background Processing
###### `calculation_worker.py` available as `CalculationWorker`
Handles heavy computational tasks to prevent GUI freezing.
- **Threading**: Run in a dedicated `QThread`.
- **Responsibilities**:
    - **3D Embedding**: Generates 3D coordinates from 2D graphs using RDKit (with Open Babel fallback).
    - **Optimization**: Performs force-field minimization (MMFF94/UFF).
    - **Stereo**: Ensures 2D stereochemistry (E/Z, R/S) is correctly preserved in 3D.

#### Main Window Architecture
The `MainWindow` class uses multiple inheritance (mixins) to separate concerns.

##### `main_window_app_state.py` available as `MainWindowAppState`
Manages the application's runtime state and undo/redo history.
- **State Capture**:
    - `get_current_state()`: Serializes atoms, bonds, and 3D data.
    - `set_state_from_data()`: Restores state, rebuilding the scene.
- **Undo/Redo**: Implements a stack-based undo system. Smartly diffs states to avoid redundant history entries.
- **Radical Support**: Correctly persists radical electrons in the serialized state.

##### `main_window_project_io.py` available as `MainWindowProjectIo`
Manages project persistence.
- **Formats**:
    - `.pmeprj` (JSON): The default, human-readable project format.
    - `.pmeraw` (Pickle): Binary format for full Python object fidelity.
- **Functionality**: Saving and loading project files, preserving full application state including UI settings.

##### `main_window_molecular_parsers.py` available as `MainWindowMolecularParsers`
Dedicated to loading and saving molecular file formats.
- **MOL/SDF**: `load_mol_file()`/`save_as_mol()` handle .mol/.sdf files, fixing headers and generating 2D coordinates if missing.
- **XYZ**: `load_xyz_file()`/`save_as_xyz()` handle .xyz files with prompts for charge and optional chemistry checks (bond perception).

##### `main_window_view_loaders.py` available as `MainWindowViewLoaders`
Specialized loaders for the 3D viewer.
- **Purpose**: Loading files specifically for visualization (read-only 3D view usually).
- **3D Logic**: Auto-generates 3D coordinates if missing (e.g. from 2D MOL files) to ensure something is shown in the 3D viewer.
- **Integration**: Switches UI to 3D mode upon successful load.

##### `main_window_ui_manager.py` available as `MainWindowUiManager`
Handles the core UI setup and event management.
- **UI Logic**: updates status bar, handles drag-and-drop.
- **Modes**: toggles between editing modes (Select, Atom, Bond, Template, Charge, Radical).
- **Layout**: manages 2D/3D splitter and visibility of panels.
##### `main_window_compute.py`
**Purpose**: Handles computational tasks and chemical conversions.
**Key Responsibilities**:
- **2D to 3D Conversion**: Triggers RDKit-based embedding to generate 3D coordinates from 2D structures (`trigger_conversion`).
- **3D Optimization**: interface for MMFF/UFF optimization of 3D structures (`optimize_3d_structure`).
- **Chemistry Checks**: Performs validation and sanitization of molecules before calculations.

##### `main_window_edit_actions.py`
**Purpose**: Manages core editing operations in the 2D editor.
**Key Responsibilities**:
- **Clipboard Operations**: Copy, Cut, Paste functionality for molecular structures.
- **Hydrogen Management**: Adding and removing explicit hydrogen atoms (`add_hydrogen_atoms`, `remove_hydrogen_atoms`).
- **Structure Cleanup**: Logic for "Clean Up 2D" to re-generate 2D coordinates.
- **Selection Tools**: Select All, Clear All, and deletion logic.

##### `main_window_edit_3d.py`
**Purpose**: Provides tools for interacting with and measuring the 3D scene.
**Key Responsibilities**:
- **Measurements**: Logic for measuring distances, angles, and dihedrals between atoms in 3D (`toggle_measurement_mode`).
- **3D Selection**: Handling selection of atoms in the 3D utility for operations like alignment and measurements.

##### `main_window_view_3d.py`
**Purpose**: Handles the rendering and visualization capabilities of the 3D viewer.
**Key Responsibilities**:
- **Drawing Styles**: Implements 3D representation styles: Ball & Stick, Stick, Wireframe, and CPK (`draw_standard_3d_style`).
- **Visualization Logic**: Uses `pyvista` and `vtk` to render atoms and bonds based on the current style settings.
- **Camera & Zoom**: Manages camera controls, zooming, and view fitting (`fit_to_view`).

##### `main_window_export.py`
**Purpose**: Detailed export capabilities for 3D and 2D content (Mesh/Media).
**Key Responsibilities**:
- **3D Model Export**: Exports the 3D scene to STL (color support limited) and OBJ/MTL formats for external modeling software (`export_obj_mtl`).
- **Image Export**: Renders and saves high-quality PNG images of the 2D or 3D views (`export_3d_png`, `export_2d_png`).

##### `main_window_string_importers.py`
**Purpose**: Facilitates importing molecules from text representations.
**Key Responsibilities**:
- **SMILES/InChI**: Dialogs and parsing logic to load molecules from SMILES or InChI strings (`load_from_smiles`, `load_from_inchi`).

##### `main_window_main_init.py`
**Purpose**: Handles the initialization phase of the main window.
**Key Responsibilities**:
- **UI Setup**: Constructs the main UI layout, toolbars, and menus (`init_ui`, `init_menu_bar`).
- **Settings Management**: Loads and applies user preferences from `settings.json`.
- **Startup Logic**: Performs RDKit warm-up and initial state configuration.

##### `main_window_dialog_manager.py`
**Purpose**: Centralized manager for application dialogs.
**Key Responsibilities**:
- **Dialog Lifecycle**: opens and tracks instances of various dialogs (Settings, Analysis, Periodic Table, etc.) to prevent duplicate windows.
- **3D Edit Dialogs**: Manages dialogs for geometric operations like `AlignPlaneDialog`, `PlanarizeDialog`, and `AlignmentDialog`.

## Complete Module Inventory

A comprehensive list of all source files in `moleditpy/modules/` and their purposes.

### Core & Data
| File | Description |
| :--- | :--- |
| `__init__.py` | Package initialization and export control. |
| `constants.py` | Global constants (colors, fonts, physical constants). |
| `molecular_data.py` | Core `MolecularData` class (atoms, bonds, graph). |
| `mol_geometry.py` | Pure math helpers for geometry checks (graph traversal, angles). |
| `sip_isdeleted_safe.py` | Safety helper for checking deleted C++ PyQt objects. |

### Main Window (Mixins)
| File | Description |
| :--- | :--- |
| `main_window.py` | Main class definition (inherits from mixins). |
| `main_window_main_init.py` | Initialization parsing and setup. |
| `main_window_ui_manager.py` | UI state, layout, and mode switching. |
| `main_window_app_state.py` | Undo/Redo and state serialization. |
| `main_window_project_io.py` | Project file (.pmeprj) loading/saving. |
| `main_window_molecular_parsers.py` | MOL/XYZ file I/O logic. |
| `main_window_string_importers.py` | SMILES/InChI import logic. |
| `main_window_view_loaders.py` | Loaders specific to the 3D viewer. |
| `main_window_view_3d.py` | 3D rendering and view management. |
| `main_window_edit_3d.py` | 3D measurements and interaction tools. |
| `main_window_edit_actions.py` | 2D editing actions (Cut/Copy/Paste, Hydrogen handling). |
| `main_window_compute.py` | Computational task management (3D conversion/opt). |
| `main_window_export.py` | Image and mesh export functionality. |
| `main_window_dialog_manager.py` | Centralized dialog opening/tracking. |

### Visualization (2D/3D)
| File | Description |
| :--- | :--- |
| `molecule_scene.py` | `QGraphicsScene` for 2D editing. |
| `zoomable_view.py` | `QGraphicsView` wrapper for 2D zooming/panning. |
| `atom_item.py` | 2D Atom graphics item. |
| `bond_item.py` | 2D Bond graphics item. |
| `template_preview_item.py` | Item for previewing templates before placement. |
| `template_preview_view.py` | Widget for viewing templates in dialogs. |
| `custom_qt_interactor.py` | `QVTKRenderWindowInteractor` wrapper. |
| `custom_interactor_style.py` | Custom VTK interaction style (picking, trackball). |

### Dialogs & Windows
| File | Description |
| :--- | :--- |
| `about_dialog.py` | Application "About" dialog. |
| `analysis_window.py` | Physicochemical properties report window. |
| `settings_dialog.py` | Main application preferences. |
| `periodic_table_dialog.py` | Element selection grid. |
| `color_settings_dialog.py` | Custom CPK color configuration. |
| `user_template_dialog.py` | User template management. |
| `translation_dialog.py` | 3D translation tool. |
| `move_group_dialog.py` | 3D substructure movement tool. |
| `align_plane_dialog.py` | Align atoms to Cartesian planes. |
| `alignment_dialog.py` | Align atoms to axes. |
| `mirror_dialog.py` | Reflect structure across a plane. |
| `planarize_dialog.py` | Flatten atoms to a best-fit plane. |
| `bond_length_dialog.py` | Adjust bond distances. |
| `angle_dialog.py` | Adjust bond angles. |
| `dihedral_dialog.py` | Adjust torsion angles. |
| `constrained_optimization_dialog.py` | Setup for constrained FF minimization. |
| `dialog3_d_picking_mixin.py` | Helper for dialogs that need 3D atom picking. |

### Background & Plugins
| File | Description |
| :--- | :--- |
| `calculation_worker.py` | Background thread implementation for 3D tasks. |
| `plugin_manager.py` | Plugin discovery and loading logic. |
| `plugin_interface.py` | Base class/interface for plugins. |
| `plugin_manager_window.py` | UI for managing installed plugins. |
