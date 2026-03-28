# MoleditPy Architecture Documentation

## Overview
MoleditPy is a Python-based molecular editing platform designed for the rapid preparation and visualization of chemical structures for DFT calculations. 

The application has undergone a significant architectural modernization, moving from a monolithic inheritance-based design to a **Composition-over-Inheritance** model. This ensures that core chemical logic, UI state management, and visual components are decoupled, highly testable, and robust against crashes in the hybrid Python/C++ environment (PyQt6/VTK).

---

## Technical Strategy

### 1. Composition-Based MainWindow
The `MainWindow` (the application's heartbeat) is no longer a single, massive class. It is now a **Composition Hub** that orchestrates specialized **Managers** and **Mixins**:

- **Managers (Delegation)**: These are independent objects (e.g., `View3DManager`, `ComputeManager`) that own specific functional areas. They are instantiated within `MainWindow` and handle complex state transitions and external library coordination (e.g., VTK and RDKit).
- **Mixins (Traits)**: These are behavioral classes (e.g., `MainWindowAppState`, `MainWindowMolecularParsers`) that provide standardized features like Undo/Redo or File I/O. The `MainWindow` inherits from these "traits" to gain their capabilities.

### 2. Backward Compatibility Layer (`moleditpy.modules`)
To prevent the major directory refactoring from breaking existing plugins, the application includes a **dynamic proxy layer**. Any plugin importing from legacy paths (e.g., `moleditpy.molecular_data`) is transparently redirected to the new location (e.g., `moleditpy.core.molecular_data`) via the `moleditpy.modules` system.

### 3. Application Lifecycle and Safety
- **Centralized Entry Point**: `main.py` serves as the unified bootstrapper.
- **Global Error Capture**: A custom `sys.excepthook` catches unhandled exceptions globally, logging them as `ERROR` entries with full tracebacks before a crash.
- **Enhanced Traceability**: Logging is configured with absolute system paths (`%(pathname)s`) for one-click navigation to failing code during development.

---

## Package Structure

| Package | Role | Key Files |
| :--- | :--- | :--- |
| **`moleditpy.core`** | **Pure Scientific Logic**: Stateless chemical algorithms and the fundamental data model. | `molecular_data.py`, `mol_geometry.py` |
| **`moleditpy.ui`** | **GUI & UI Logic**: PyQt6/VTK components, Managers, and I/O logic coupled with the interface. | `main_window.py`, `compute_engine.py`, `project_io.py`, `view_3d_logic.py` |
| **`moleditpy.utils`** | **Common Utilities**: Threading helpers, system constants, and safe SIP interaction. | `constants.py`, `sip_isdeleted_safe.py`, `system_utils.py` |
| **`moleditpy.plugins`**| **Extension Layer**: Discovery and dynamic loading of Python scripts. | `plugin_manager.py`, `plugin_interface.py` |
| **`moleditpy.modules`**| **Compatibility Proxy**: Redirects legacy imports to the refactored structure. | `__init__.py` |
| **`moleditpy.assets`** | **Media Content**: Application icons and graphical assets. | `icon.png`, `icon.ico` |

---

## Developer Feature Map

| Feature Area | Core Logic (`core/`) | UI Manager (`ui/`) | Logic Helper (`ui/*_logic.py`) |
| :--- | :--- | :--- | :--- |
| **Data Model** | `molecular_data.py` | `molecule_scene.py` | |
| **3D Conversion** | | `compute_engine.py` | `compute_logic.py` |
| **Optimization** | | `calculation_worker.py` | `compute_logic.py` |
| **File I/O** | | `project_io.py` | |
| **MOL/XYZ/SMILES** | | `molecular_parsers.py` | |
| **3D Rendering** | | `view_3d.py` | `view_3d_logic.py` |
| **2D Editing** | | `edit_actions.py` | `edit_actions_logic.py` |
| **3D Editing** | | `edit_3d.py` | `edit_3d_logic.py` |
| **Undo/Redo** | | `app_state.py` | (Integrated into `MainWindow`) |
| **Dialogs** | | `dialog_manager.py` | `dialog_logic.py` |

---

## Architectural Patterns

### 1. Separation of Concerns (Logic-UI Decoupling)
The project actively migrates business logic out of UI classes. 
- **Pattern**: A UI class (e.g., `MainWindow`) inherits from a mixin (e.g., `MainWindowCompute`), which in turn delegates complex operations to a logic-only module (e.g., `ui/compute_logic.py` or `core/compute_engine.py`).
- **Benefit**: This allows the core chemical logic to be tested without instantiating heavy GUI components.

### 2. Defensive GUI Management
To avoid crashes in the hybrid Python/C++ environment:
- **Rule**: Always check if a C++ object still exists using `sip_isdeleted_safe` before access.
- **Rule**: Wrap UI-triggered computations in logged exception handlers. The application also benefits from a **global exception hook** in `main.py` that captures crashes in GUI slots or background threads as a secondary safety net.

### 3. Trait-based MainWindow
The `MainWindow` class is a thin wrapper that inherits from over a dozen mixins (located in `moleditpy/src/moleditpy/ui/`). This "traits" pattern keeps feature areas (Export, Edit, Compute, IO) isolated while providing a unified API for the user.

---

## Specialized Subsystems

For in-depth technical details on specific subsystems, refer to:

- [CORE_DATA_STRUCTURES.md](CORE_DATA_STRUCTURES.md): Data model and Scene Graph details.
- [ANALYSIS_AND_CALCULATIONS.md](ANALYSIS_AND_CALCULATIONS.md): 3D generation and optimization pipelines.
- [DIALOGS_AND_UI.md](DIALOGS_AND_UI.md): UI component inventory and custom widgets.
- [PLUGIN_SYSTEM_INTERNAL.md](PLUGIN_SYSTEM_INTERNAL.md): Technical details of the extension API.
