# MoleditPy Test Suite

## 1\. Overview

This directory contains the automated test suite for the **MoleditPy** molecular editor. These tests are designed to verify the application's core logic (unit tests) and its graphical user interface interactions (GUI tests).

The suite uses the `pytest` framework, leveraging `pytest-qt` for GUI testing and `pytest-mock` (or the built-in `unittest.mock`) for isolating components.

## 2\. Technologies & Frameworks

  * **Testing Framework**: `pytest`
  * **GUI Testing**: `pytest-qt`
  * **Mocking**: `pytest-mock` (and `unittest.mock`)
  * **Core Application**: `PyQt6`
  * **Mocked Dependencies**: `RDKit` (for 3D generation), `PyVista` (for 3D rendering), and various `PyQt6` dialogs (`QMessageBox`, `QFileDialog`) to ensure non-blocking test execution.

## 3\. Test Structure

The suite is organized into the following key files:

### `pytest.ini`

This file configures the `pytest` runner:

  * **Markers**: Defines custom markers to separate test types:
      * `unit`: For unit tests that do not require a GUI.
      * `gui`: For tests that interact with the PyQt6 GUI.
  * **Configuration**: Sets options like `norecursedirs` to avoid external packages and `addopts` to enforce strict marker usage.

### `conftest.py`

This fixture file sets up the test environment for `pytest`.

  * **`app` fixture**: Provides a session-scoped `QApplication` instance.
  * **`window` fixture**: The primary setup fixture. It runs for each test, creating a fresh `MainWindow` instance.
  * **Mocking**: This fixture is responsible for heavily mocking external or slow components to ensure stable and fast tests:
      * **3D Conversion**: Mocks the `CalculationWorker` thread. When a 2D-to-3D conversion is triggered, it immediately calls the `on_calculation_finished` slot with a dummy RDKit molecule, bypassing the actual computation.
      * **PyVista/VTK**: Replaces the `CustomQtInteractor` (the 3D viewport) with a `DummyPlotter` (a simple `QWidget`) to prevent an actual 3D window from launching.
      * **Dialogs**: Mocks all blocking dialogs (`QMessageBox`, `QFileDialog`, `QInputDialog`, `QDialog.exec`) to return default "success" values (e.g., "Yes" to questions, a fake path for file dialogs).

### `test_plugin_manager.py`

This file verifies the **integration aspects** of the plugin system, specifically file system interactions and safe metadata parsing.

  * **Tests**:
      * `test_install_and_discover_single_file`: Verifies the end-to-end installation of a single `.py` plugin file.
      * `test_install_zip`: Verifies the extraction and installation of zipped plugins.
      * `test_plugin_registration`: Checks if a plugin can correctly register UI actions (menus/toolbars) in the `PluginContext`.
      * `test_ast_metadata_parsing`: Ensures safe, non-executing extraction of plugin metadata (Author, Version) using AST.

### `test_main_app.py`

This is the primary integration test suite, executing user workflows in a headless environment.

| Test Function | Feature Under Test | Verification Scope |
| :--- | :--- | :--- |
| **`test_app_launch`** | **Startup & Initialization** | Verifies that the application launches correctly, ensuring all UI components, database connections, and signal slots are wired without runtime errors. |
| **`test_mode_change_atom`** | **Toolbar Logic** | Verifies the `ModeManager` state transitions, ensuring that clicking toolbar buttons correctly updates the global application mode and cursor state. |
| **`test_draw_atom_on_click`** | **Scene Events (Mouse)** | Verifies the entire event chain for atom creation: `Mouse Press` -> `Scene Event Filter` -> `Coordinate Transform` -> `Data Model Update`. |
| **`test_draw_bond_on_drag`** | **Complex Interactions** | Verifies the state machine for bond creation, testing the press-drag-release sequence and satisfying valid valency checks. |
| **`test_2d_to_3d_conversion`** | **Scientific Workflow** | Verifies the critical handoff between the 2D sketcher and the 3D embedding engine, confirming that 2D topology is correctly interpreted by the RDKit conformer generator. |
| **`test_undo_redo`** | **Command Pattern** | Validates the Undo/Redo stack, ensuring that state snapshots are correctly captured before destructive actions and restored accurately. |
| **`test_copy_paste`** | **Clipboard & MIME** | Verifies custom serialization logic, ensuring that internal molecular data can be copied and pasted within the app (and potentially to external apps favoring text formats). |
| **`test_file_import_smiles`** | **Data Ingestion** | Tests the integration of the SMILES parser with the UI thread, ensuring valid strings generate the correct 2D graph structure. |
| **`test_save_project_as`** | **Project Persistence** | Verifies the full JSON serialization pipeline, checking that all workspace state (including viewport camera positions) is saved to disk. |
| **`test_toggle_3d_atom_info`** | **View Options** | Verifies that toggling UI view options triggers the correct update in the 3D render window (e.g., generating/removing text actors). |
| **`test_draw_bond_to_existing_atom`** | **Graph Connectivity** | Verifies the snapping logic, ensuring that drawing a bond to an existing atom correctly merges the graph nodes rather than creating overlaps. |
| **`test_delete_atom_on_right_click`** | **Data Integrity** | Verifies cascading deletion logic: ensuring that removing an atom also cleanly removes its incident bonds and updates neighbor lists. |
| **`test_key_press_change_atom`** | **Shortcuts & Hotkeys** | Verifies that keyboard inputs are correctly intercepted and mapped to element changes or tool switches. |
| **`test_add_remove_hydrogens`** | **Algorithmic Utilities** | Verifies that the 'Add/Remove Hydrogens' UI commands correctly invoke the backend valency algorithms and update the visualization. |
| **`test_additional_dialogs_launch.py`** | **Secondary Dialogs** | Smoke tests for less common dialogs (e.g., Mirror, Align Plane). |
| **`test_dialog_launch.py`** | **Dialog Smoke Tests** | Ensures all primary dialogs can be instantiated and shown without errors. |
| **`test_main_window_settings.py`** | **Settings Integration** | Verifies that settings applied in the dialog correctly propagate to the main application. |
| **`test_main_window_ui_integration.py`** | **UI/Data Sync** | Verifies consistent state between the data model and multiple UI views. |
| **`test_molecule_scene_events.py`** | **Scene Events** | Detailed verification of mouse and keyboard event filtering in the 2D scene. |
| **`test_plugin_manager_redundant.py`** | **Plugin Edge Cases** | Tests for redundant plugin registration and error handling in plugin loading. |
| **`test_edit_actions_gui_extended.py`** | **Edit Actions GUI** | Verifies the `Rotate2DDialog` GUI in a full PyQt6 context.<br>• **Angle Sync**: Confirms the spinbox and slider remain synchronized during initialization and value changes. |

## 4\. Setup & Running

### Installation

Before running the tests, install the required `pytest` plugins:

```bash
pip install pytest pytest-qt pytest-mock
```

You will also need the application's main dependencies, such as `PyQt6`.
