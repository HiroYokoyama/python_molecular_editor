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

### `test_main_app.py`

This is the main test file, containing all unit and GUI test cases.

  * **Test Helpers**:

      * `get_action(toolbar, tooltip_text)`: Finds a `QAction` on a toolbar.
      * `click_scene(...)`: Simulates a mouse click at a specific `QPointF` in the 2D scene.
      * `drag_scene(...)`: Simulates a mouse drag between two points in the 2D scene.

  * **Unit Tests (`@pytest.mark.unit`)**:

      * Focus on the `MolecularData` class.
      * Test adding/removing atoms and bonds.
      * Verify correct handling of bond sorting (normal vs. stereo bonds).
      * Test RDKit conversion logic, including 2D stereochemistry (Wedge/Dash) and E/Z bond stereo.

  * **GUI Tests (`@pytest.mark.gui`)**:

      * Cover all major application features by simulating user interaction via `qtbot`.
      * **2D Editor**: Drawing atoms, bonds, templates (Benzene); changing modes (charge, radical); Undo/Redo; Copy/Paste; Clear All; 2D Cleanup.
      * **3D Conversion**: Clicking the "Convert 2D -\> 3D" button and verifying the UI state change.
      * **3D View**: Toggling 3D styles (CPK), measurement mode, 3D drag mode, and atom info labels.
      * **File I/O**: Importing SMILES, saving/loading projects (`.pmeprj`), loading `.mol` files in viewer mode, and drag-and-drop file handling.
      * **Keyboard Shortcuts**: Testing key presses ('O' for Oxygen, '2' for double bond, 'Delete' for deletion).

## 4\. Setup & Running

### Installation

Before running the tests, install the required `pytest` plugins:

```bash
pip install pytest pytest-qt pytest-mock
```

You will also need the application's main dependencies, such as `PyQt6`.

### Running Tests

To run the test suite, execute `pytest` from the root directory:

```bash
# Run all tests
pytest

# Run only unit tests
pytest -m unit

# Run only GUI tests
pytest -m gui
```
