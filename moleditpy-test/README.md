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

  * **Test Cases List**:

      * **Unit Tests (`@pytest.mark.unit`)**:

          * `test_molecular_data_add_atom`: Tests adding an atom to the data model.
          * `test_molecular_data_add_bond`: Tests adding a standard bond.
          * `test_molecular_data_add_stereo_bond`: Tests that stereo bonds (Wedge/Dash) are stored correctly without sorting.
          * `test_molecular_data_remove_atom`: Tests that removing an atom also removes its associated bonds.
          * `test_molecular_data_remove_bond`: Tests removing a bond.
          * `test_to_rdkit_mol_stereo`: Tests conversion of Wedge/Dash 2D stereo to RDKit.
          * `test_to_rdkit_mol_ez_stereo`: Tests conversion of E/Z bond stereo to RDKit.

      * **GUI Tests (`@pytest.mark.gui`)**:

          * `test_app_launch`: Checks if the application window launches correctly.
          * `test_mode_change_atom`: Tests changing the drawing mode via the atom toolbar (e.g., to 'N').
          * `test_mode_change_bond`: Tests changing the drawing mode via the bond toolbar (e.g., to 'Double Bond').
          * `test_draw_atom_on_click`: Tests creating an atom by clicking on the scene.
          * `test_draw_bond_on_drag`: Tests creating two atoms and a bond by dragging on the scene.
          * `test_draw_bond_to_existing_atom`: Tests drawing a bond between two existing atoms.
          * `test_change_atom_symbol_on_click`: Tests changing an existing atom's element by clicking on it.
          * `test_change_bond_order_on_click`: Tests changing an existing bond's order by clicking on it.
          * `test_delete_atom_on_right_click`: Tests deleting an atom via right-click.
          * `test_charge_mode_click`: Tests applying positive and negative charges to an atom.
          * `test_2d_to_3d_conversion`: Tests the 2D to 3D conversion button and UI state change.
          * `test_optimize_3d`: Tests the 3D optimization button.
          * `test_change_3d_style`: Tests changing the 3D view style (e.g., to 'CPK').
          * `test_undo_redo`: Tests the Undo and Redo actions.
          * `test_clear_all`: Tests the "Clear All" action.
          * `test_copy_paste`: Tests copying and pasting a selection.
          * `test_file_import_smiles`: Tests importing a structure from a SMILES string.
          * `test_key_press_change_atom`: Tests changing an atom's element via keyboard shortcut (e.g., 'o' for Oxygen).
          * `test_key_press_change_bond`: Tests changing a bond's order via keyboard shortcut (e.g., '2' for double bond).
          * `test_radical_mode_toggle`: Tests toggling radical states (0, 1, 2) on an atom.
          * `test_delete_key_selection`: Tests deleting selected items using the 'Delete' key.
          * `test_draw_benzene_template`: Tests placing the benzene template.
          * `test_open_settings_dialog`: Tests opening the 3D View Settings dialog.
          * `test_toggle_measurement_mode`: Tests toggling the 3D measurement/selection mode.
          * `test_toggle_3d_edit_mode`: Tests toggling the 3D atom drag mode.
          * `test_add_remove_hydrogens`: Tests the "Add Hydrogens" and "Remove Hydrogens" menu actions.
          * `test_2d_cleanup`: Tests the "Optimize 2D" (cleanup) button.
          * `test_3d_viewer_mode_mol`: Tests loading a .mol file directly into 3D viewer mode, disabling 2D editing.
          * `test_open_3d_edit_dialogs`: Tests if 3D edit dialogs (Translate, Planarize) can be opened.
          * `test_save_project_as`: Tests saving the current state as a .pmeprj project file.
          * `test_open_project`: Tests loading a .pmeprj project file.
          * `test_toggle_3d_atom_info`: Tests toggling 3D atom info labels (ID, Coords, Symbol).
          * `test_user_template_dialog_save_and_use`: Tests saving the current 2D structure as a user template and using it.
          * `test_implicit_hydrogens_update`: Tests that implicit hydrogen counts update automatically after drawing.
          * `test_drag_drop_mol_file_on_3d_view`: Tests drag-and-drop of a .mol file onto the 3D view (for 3D viewer mode).
          * `test_drag_drop_mol_file_on_2d_view`: Tests drag-and-drop of a .mol file onto the 2D view (for 2D editing).

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
