# MoleditPy Architecture Documentation

## Overview
MoleditPy is a Python-based molecular editing software built with **PyQt6** for the GUI and **RDKit** for chemical informatics features (SMILES parsing, 3D conformation generation, etc.).

The application follows a modular, decoupled design. The codebase is organized into specialized packages, and complex UI components are decomposed into view-specific code and independent logic modules.

## Developer Guide / Feature Map

Use this map to find the code responsible for specific features. The project is located in `moleditpy/src/moleditpy/`.

| Feature Area | Core Logic (`core/`) | UI Component (`ui/`) | Logic Helper (`ui/*_logic.py`) |
| :--- | :--- | :--- | :--- |
| **Data Model** | `molecular_data.py` | `molecule_scene.py` | |
| **3D Conversion** | `compute_engine.py` | `main_window.py` | `compute_logic.py` |
| **Optimization** | `calculation_worker.py` | `main_window.py` | `compute_logic.py` |
| **File I/O** | `project_io.py` | `main_window.py` | |
| **MOL/XYZ/SMILES** | `molecular_parsers.py` | `main_window.py` | |
| **3D Rendering** | | `view_3d.py` | `view_3d_logic.py` |
| **2D Editing** | | `edit_actions.py` | `edit_actions_logic.py` |
| **3D Editing** | | `edit_3d.py` | `edit_3d_logic.py` |
| **Undo/Redo** | `app_state.py` | `main_window.py` | |
| **Dialogs** | | `dialog_manager.py` | `dialog_logic.py` |

---

## Package Structure

### 1. `moleditpy.core`
Contains all "pure" chemical and application logic, independent of the GUI. This package has a **strict >80% test coverage** requirement.
- **`molecular_data.py`**: The central `MolecularData` class. Manages the chemical graph (atoms, bonds).
- **`compute_engine.py`**: Orchestrates 3D embedding and coordination between RDKit and Open Babel.
- **`calculation_worker.py`**: Thread-safe worker for background calculations.
- **`molecular_parsers.py`**: Robust parsers for MOL, SDF, and XYZ formats.
- **`mol_geometry.py`**: Geometric calculations and chemical integrity checks.
- **`app_state.py`**: State serialization for projects and Undo/Redo operations.
- **`project_io.py`**: Persistence logic for `.pmeprj` and `.pmeraw` files.

### 2. `moleditpy.ui`
Contains all PyQt6/VTK-based interface components.
- **`main_window.py`**: The top-level window, composed of multiple functional mixins.
- **`molecule_scene.py`**: Interactive 2D drawing canvas based on `QGraphicsScene`.
- **`view_3d.py`**: High-performance 3D visualization using `PyVista`.
- **`*_logic.py`**: Special modules (e.g., `compute_logic.py`, `edit_actions_logic.py`) that extract complex state management and decision logic from the UI classes to improve testability.
- **Dialogs**: All application dialogs (Settings, Periodic Table, Geometry controls) are located here.

### 3. `moleditpy.utils`
Common utilities used across the project.
- **`constants.py`**: Shared colors, styles, and physical constants.
- **`sip_isdeleted_safe.py`**: Utilities for safe interaction with underlying C++ objects in PyQt.

### 4. `moleditpy.plugins`
The extensibility layer.
- **`plugin_manager.py`**: Handles discovery and dynamic loading of Python-based plugins.
- **`plugin_interface.py`**: Defines the API for third-party extensions.

---

## Architectural Patterns

### 1. Separation of Concerns (Logic-UI Decoupling)
The project actively migrates business logic out of UI classes. 
- **Pattern**: A UI class (e.g., `MainWindow`) inherits from a mixin (e.g., `MainWindowCompute`), which in turn delegates complex operations to a logic-only module (e.g., `ui/compute_logic.py` or `core/compute_engine.py`).
- **Benefit**: This allows the core chemical logic to be tested without instantiating heavy GUI components.

### 2. Defensive GUI Management
To avoid crashes in the hybrid Python/C++ environment:
- **Rule**: Always check if a C++ object still exists using `sip_isdeleted_safe` before access.
- **Rule**: Wrap UI-triggered computations in broad but logged exception handlers to ensure application stability even if a specific operation fails.

### 3. Trait-based MainWindow
The `MainWindow` class is a thin wrapper that inherits from over a dozen mixins (located in `moleditpy/src/moleditpy/ui/`). This "traits" pattern keeps feature areas (Export, Edit, Compute, IO) isolated.

---

## Specialized Subsystems

For in-depth technical details on specific subsystems, refer to:

- [CORE_DATA_STRUCTURES.md](CORE_DATA_STRUCTURES.md): Data model and Scene Graph details.
- [ANALYSIS_AND_CALCULATIONS.md](ANALYSIS_AND_CALCULATIONS.md): 3D generation and optimization pipelines.
- [DIALOGS_AND_UI.md](DIALOGS_AND_UI.md): UI component inventory and custom widgets.
- [PLUGIN_SYSTEM_INTERNAL.md](PLUGIN_SYSTEM_INTERNAL.md): Technical details of the extension API.
