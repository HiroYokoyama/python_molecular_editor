# MoleditPy Unit Test Suite

This directory contains standalone unit tests for the core scientific logic and application modules of MoleditPy.

## Test Descriptions

### Unit Test Coverage & Logic Justification

| Test File | Component | Verification Scope |
| :--- | :--- | :--- |
| **`test_app_logic.py`** | **Core Application Orchestration** | Verifies the thread safety and state management of the application. Key scenarios include:<br>• **Async Worker Control**: Ensuring calculation threads correctly respond to halt signals.<br>• **Stereochemistry Persistence**: Confirming E/Z flags survive round-trip serialization even when backend sanitization might alter them.<br>• **Robust Fallbacks**: Testing custom `MolBlock` generation robustness when standard libraries fail on invalid valencies.<br>• **Coordinate Systems**: Validating the precise application of `ANGSTROM_PER_PIXEL` scaling during 2D-to-3D conversions. |
| **`test_hydrogen.py`** | **Chemical Intelligence** | Verifies the logic for implicit vs. explicit hydrogen management.<br>• **Smart Addition**: checks that hydrogens are added according to free valency rules.<br>• **Selective Removal**: Ensures only explicit hydrogens are removed while strictly preserving stereochemical configuration. |
| **`test_stereochemistry.py`** | **Stereo Flags & UI Mapping** | Verifies the translation layer between the internal RDKit representation and the UI's visual indicators.<br>• **Chiral Mapping**: Ensures atomic chirality is correctly parsed and visualized.<br>• **Bond Stereo**: Validates that double bond geometry (E/Z) is correctly extracted from the connectivity table and stored in the application model. |
| **`test_geometry.py`** | **2D Geometrical Operations** | Verifies purely geometrical manipulations of the 2D scene coordinates.<br>• **Mirroring/Inversion**: Checks algorithms for coordinate inversion (enantiomer generation) and axis mirroring, ensuring they operate correctly on the 2D canvas independent of 3D conformation. |
| **`test_properties.py`** | **Chemical Properties Analysis** | Verifies the calculation and presentation of molecular descriptors.<br>• **Descriptor Accuracy**: Validates calculations like MW and LogP against known standards.<br>• **Error Stability**: Ensures the analysis window handles calculation failures gracefully without crashing the UI. |
| **`test_parsers.py`** | **File I/O & Parsing** | Verifies the robustness of data ingestion modules.<br>• **XYZ/MOL Loading**: Tests custom parsers designed to handle non-standard or slightly malformed files that strict standard libraries might reject. |
| **`test_io.py`** | **Project Serialization** | Verifies the integrity of the application's native save format.<br>• **.pmeprj Integrity**: Tests the generic JSON structure, ensuring all custom attributes (canvas positions, view states) round-trip correctly. |

## Mocking Strategy
These tests use a `MockMainWindow` (defined in `conftest.py`) which simulates the host environment with standardized `CalculationWorker`, `MolecularData`, and logic implementations from `main_window_app_state.py`.
