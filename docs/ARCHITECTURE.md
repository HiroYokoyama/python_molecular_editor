# MoleditPy Architecture Documentation

## Overview
MoleditPy is a Python-based molecular editing platform designed for the rapid preparation and visualization of chemical structures for DFT calculations. 

The application has undergone a significant architectural modernization, moving from a monolithic inheritance-based design to a **Composition-over-Inheritance** model. This ensures that core chemical logic, UI state management, and visual components are decoupled, highly testable, and robust against crashes in the hybrid Python/C++ environment (PyQt6/VTK).

---

## Technical Strategy

### 1. Composition-Based Manager Architecture
The `MainWindow` (the application's heartbeat) is no longer a collection of inheritance-heavy mixins. It is now a **Composition Hub** that delegates functional responsibilities to specialized **Managers**:

- **Composition-over-Inheritance**: Functional logic (e.g., File I/O, 3D Rendering, Scientific Editing) is encapsulated within independent Manager classes (`IOManager`, `View3DManager`, `StateManager`, etc.).
- **Host-Centric Delegation**: Managers reference the `MainWindow` via a `self.host` property, ensuring clear boundaries and preventing the "massive object" anti-pattern.
- **Unified State**: Application state (undo stacks, current selection, mode configuration) is centralized in `StateManager`, ensuring absolute synchronization across decoupled components.

### 2. "Never Hide Errors" Diagnostic Strategy
To ensure stability in the hybrid Python/C++ environment (PyQt6/VTK), the application employs a **Diagnostic Hardening** strategy:

- **Visible Failures**: We prioritize "loud" failures over silent inconsistencies. Every `hasattr()` check or `try-except` block used for architectural routing must include a diagnostic `logging.error()` or re-raise.
- **Active Hardening**: The codebase is instrumented with defensive `else: logging.error()` branches for missing attributes, providing high-fidelity tracebacks directly from the architectural surface.
- **Global Error Capture**: A custom `sys.excepthook` in `main.py` catches unhandled exceptions, ensuring that even unpredicted GUI crashes leave a clear audit trail.
- **Path-Aware Logging**: Logs include absolute system paths (`%(pathname)s`) for immediate navigation to failing architectural disconnections.

---

## Package Structure

| Package | Role | Key Files |
| :--- | :--- | :--- |
| **`moleditpy.core`** | **Pure Scientific Logic**: Stateless chemical algorithms and the fundamental data model. | `molecular_data.py`, `mol_geometry.py` |
| **`moleditpy.ui`** | **GUI & UI Logic**: PyQt6/VTK components, Managers, and I/O logic coupled with the interface. | `main_window.py`, `compute_logic.py`, `io_logic.py`, `view_3d_logic.py` |
| **`moleditpy.utils`** | **Common Utilities**: Threading helpers, system constants, and safe SIP interaction. | `constants.py`, `sip_isdeleted_safe.py`, `system_utils.py` |
| **`moleditpy.plugins`**| **Extension Layer**: Discovery and dynamic loading of Python scripts. | `plugin_manager.py`, `plugin_interface.py` |
| **`moleditpy.assets`** | **Media Content**: Application icons and graphical assets. | `icon.png`, `icon.ico` |

---

## Developer Feature Map

| Feature Area | Core Logic (`core/`) | UI Manager (`ui/`) | Logic Helper (`ui/*_logic.py`) |
| :--- | :--- | :--- | :--- |
| **Data Model** | `molecular_data.py` | `molecule_scene.py` | |
| **3D Conversion** | | `calculation_worker.py` | `compute_logic.py` |
| **Optimization** | | `calculation_worker.py` | `compute_logic.py` |
| **File I/O** | | `io_logic.py` | |
| **MOL/XYZ/SMILES** | | `string_importers.py` | |
| **3D Rendering** | | | `view_3d_logic.py` |
| **2D Editing** | | `edit_actions_logic.py` | |
| **3D Editing** | | `edit_3d_logic.py` | |
| **Undo/Redo** | | `app_state.py` | (Integrated into `MainWindow`) |
| **Dialogs** | | `dialog_logic.py` | |

---

## Architectural Patterns

### 1. Separation of Concerns (Logic-UI Decoupling)
The project identifies a strict boundary between UI presentation and scientific logic:
- **Pattern**: A UI Manager (e.g., `ComputeManager`) aggregates complex scientific algorithms and interacts with the GUI host via the `self.host` property.
- **Benefit**: This allows the core chemical logic to be compartmentalized, preventing the intermingling of VTK/RDKit operations with PyQt6 boilerplate.

### 2. State Management (Single Source of Truth)
Instead of distributing application states (undo stacks, mode flags, molecule data) across multiple mixins, the application uses a unified **StateManager**:
- **Consistency**: All cross-component interactions (e.g., updating the 3D view after a 2D edit) are mediated by the `StateManager` to ensure functional synchronization.
- **Undo/Redo Integrity**: All state transitions are captured in a centralized undo stack, ensuring that the application can always return to a known-good configuration.

### 3. Manager-Host Delegation (The 'Host' Pattern)
The `MainWindow` acts as a pure container for specialized Managers.
- **Explicit Access**: Managers communicate with the host explicitly through `self.host`, avoiding the ambiguity of legacy method inheritance.
- **Lifecycle Control**: Managers are initialized by the `MainInitManager`, ensuring a clear startup and shutdown sequence.

---

## Specialized Subsystems

For in-depth technical details on specific subsystems, refer to:

- [CORE_DATA_STRUCTURES.md](CORE_DATA_STRUCTURES.md): Data model and Scene Graph details.
- [ANALYSIS_AND_CALCULATIONS.md](ANALYSIS_AND_CALCULATIONS.md): 3D generation and optimization pipelines.
- [DIALOGS_AND_UI.md](DIALOGS_AND_UI.md): UI component inventory and custom widgets.
- [PLUGIN_SYSTEM_INTERNAL.md](PLUGIN_SYSTEM_INTERNAL.md): Technical details of the extension API.
