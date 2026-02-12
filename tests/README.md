# MoleditPy Test Suite

This directory contains the testing infrastructure for MoleditPy.

## Detailed Test Logic Documentation

This section describes exactly what logic each test suite verifies, ensuring that tests cover meaningful application behavior rather than just implementation details.

### Unit Tests (`tests/unit/`)

| Test File | Verified Logic | Key Scenarios | App Logic vs Library |
| :--- | :--- | :--- | :--- |
| **`test_app_logic.py`** | **Async Workflow & Safety** | • **Worker Halt**: Verifies that calculation threads respect the halt flag.<br>• **E/Z Stereo**: Checks that stereo labels are preserved during round-trips even if RDKit sanitization might drop them.<br>• **Fallback Serialization**: Ensuring that custom `MolBlock` generation works when RDKit fails on invalid atoms.<br>• **Direct Mode**: Verifies the 'direct' calculation path for simple conversions.<br>• **Coordinate Mapping**: Verifies that custom pixel-to-Angstrom scaling (`ANGSTROM_PER_PIXEL`) is correctly applied before RDKit processing. | **Not RDKit**: RDKit doesn't have halt flags, `conversion_mode`, or pixel scaling. This tests **MoleditPy's** orchestration. |
| **`test_hydrogen.py`** | **Chemical Validity** | • **H Addition**: Checks intelligent hydrogen placement based on valency.<br>• **H Removal**: Verifies only explicit hydrogens are removed while preserving stereochemistry. | **Not RDKit**: MoleditPy has custom logic to decide *which* H to remove based on user selection/mode. |
| **`test_stereochemistry.py`** | **Stereo Mapping** | • **Chiral Flags**: Ensures chiral centers are correctly identified and mapped back to the UI.<br>• **Bond Stereo**: Verifies E/Z bond attributes are correctly parsed and stored. | **Not RDKit**: Tests the *mapping* between RDKit's internal flags and MoleditPy's UI flags (e.g. `stereo=3` vs `STEREOZ`). |
| **`test_geometry.py`** | **Geometric Transformations** | • **Mirroring**: Verifies atom coordinate inversion along specified axes.<br>• **Inversion**: Checks point inversion logic for creating enantiomers. | **Not RDKit**: RDKit has conformer generation, but MoleditPy handles the immediate geometric inversion of *2D scene coordinates*. |
| **`test_properties.py`** | **Analysis Calculations** | • **MW/LogP**: Verifies accurate calculation of molecular descriptors.<br>• **Error Handling**: Checks adequate fallback when property calculation fails. | **Not RDKit**: Tests error handling and data formatting specific to the Analysis Window. |

### GUI Tests (`tests/gui/`)

These tests verify the integration of the UI with the underlying logic. They run in a headless environment but simulate real user actions.

| Test Function (`test_main_app.py`) | Verified Application Behavior |
| :--- | :--- |
| **`test_draw_molecule`** | **Canvas Interaction**: Verifies that clicking creates atoms and dragging creates bonds using the scene's event logic. |
| **`test_undo_redo`** | **State Management**: Checks that the command history correctly tracks and restores molecular state. |
| **`test_copy_paste`** | **Clipboard Logic**: Verifies correct serialization to/from generic MIME types and internal formats. |
| **`test_2d_to_3d_conversion`** | **Workflow Integration**: Verifies the handoff from 2D canvas -> RDKit embedding -> 3D View state update. |
| **`test_toggle_3d_atom_info`** | **View State**: Asserts that toggling view options correctly updates internal display state and generates label actors (verifies the *effect*). |
| **`test_save_project_as`** | **Persistence**: Verifies the full project save workflow, ensuring JSON structure integrity. |

### Unified Runner (`run_all_tests.py`)

For convenience, a unified runner script is provided to execute **both** Unit and GUI tests in a single pass, ensuring the environment is correctly set up.

```bash
# Run all tests (Unit + GUI)
python tests/run_all_tests.py

# Run all tests in HEADLESS mode (no windows will pop up)
python tests/run_all_tests.py --headless
```

This script:
1.  Sets `PYTHONPATH` to prioritize the local `src` directory.
2.  Runs `tests/unit` first.
3.  Runs `tests/gui` second.
