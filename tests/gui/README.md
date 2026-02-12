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

This file contains tests for the `PluginManager` module, ensuring robust plugin handling.

  * **Tests**:
      * `test_init`: Verifies plugin manager initialization.
      * `test_install_and_discover_single_file`: Tests installing and loading single-file plugins.
      * `test_install_zip`: Tests installing and extracting ZIP plugins.
      * `test_plugin_registration`: Verifies `PluginContext` action registration.
      * `test_ast_metadata_parsing`: Tests safe metadata extraction from plugin files.

### `test_main_app.py`

This is the main test file, containing all unit and GUI test cases.

  * **Test Helpers**:

      * `get_action(toolbar, tooltip_text)`: Finds a `QAction` on a toolbar.
      * `click_scene(...)`: Simulates a mouse click at a specific `QPointF` in the 2D scene.
      * `drag_scene(...)`: Simulates a mouse drag between two points in the 2D scene.

  * **Test Cases List**:

      * **Unit Tests (`@pytest.mark.unit`)**:

          * See `tests/unit/README.md` for details on unit test coverage.

      * **GUI Tests (`@pytest.mark.gui`)**:

          | Test Function | Verified Application Behavior | Meaningful Logic Verified |
          | :--- | :--- | :--- |
          | `test_app_launch` | Window initialization | • Verifies critical startup wiring (signals, default modes). |
          | `test_mode_change_atom` | Toolbar interaction | • Verifies global `ModeManager` state updates safely. |
          | `test_draw_atom_on_click` | Scene interaction | • Verifies `QGraphicsScene` event handling logic, coordinate transformation, and data model sync. |
          | `test_draw_bond_on_drag` | Drag & Drop logic | • Verifies complex state machine transitions (press -> drag -> release) for bond creation. |
          | `test_2d_to_3d_conversion` | Workflow Integration | • Verifies the entire pipeline: 2D Data -> RDKit embedding -> 3D Actor generation -> View update. |
          | `test_undo_redo` | State Persistence | • Verifies the `Command` pattern implementation (drawing -> undo -> restore state). |
          | `test_copy_paste` | Clipboard Logic | • Verifies custom MIME type serialization/deserialization logic. |
          | `test_file_import_smiles` | I/O Logic | • Verifies parsing logic integration with the UI thread. |
          | `test_save_project_as` | Persistence | • Verifies JSON serialization structure and file writing logic. |
          | `test_toggle_3d_atom_info` | View Logic | • Verifies that `show_atom_info` correctly iterates actors and updates labels based on internal state. |
          | `test_draw_bond_to_existing_atom` | Connectivity | • Verifies logic for snapping to existing atoms and updating adjacency lists. |
          | `test_delete_atom_on_right_click` | Data Integrity | • Verifies cascading deletes (removing usage in bonds/adjacency) when an atom is removed. |
          | `test_key_press_change_atom` | Keyboard Shortcuts | • Verifies that key events correctly map to specific scientific commands. |
          | `test_add_remove_hydrogens` | Algorithmic Actions | • Verifies the UI trigger correctly invokes the backend valency calculation logic. |

## 4\. Setup & Running

### Installation

Before running the tests, install the required `pytest` plugins:

```bash
pip install pytest pytest-qt pytest-mock
```

You will also need the application's main dependencies, such as `PyQt6`.
