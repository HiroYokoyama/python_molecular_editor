# MoleditPy Unit Test Suite

This directory contains standalone unit tests for the core scientific logic and application modules of MoleditPy.

## Test Descriptions

### Unit Test Coverage & Logic Justification

| Test File | Verified Application Logic | App Logic vs Library (Why this test matters) |
| :--- | :--- | :--- |
| **`test_app_logic.py`** | **Async Workflow & Safety** | **Not RDKit**: Verifies MoleditPy's specific worker orchestration, halt flags, and custom fallback serialization that runs when external libraries fail. |
| **`test_hydrogen.py`** | **Chemical Validity** | **Not RDKit**: Verifies MoleditPy's decision-making on *which* hydrogens to add/remove based on selection state. |
| **`test_stereochemistry.py`** | **Stereo Mapping** | **Not RDKit**: Tests the *translation layer* between RDKit's internal flags and MoleditPy's UI bond properties (e.g. `stereo=3` vs `STEREOZ`). |
| **`test_geometry.py`** | **Geometric Transformations** | **Not RDKit**: Verifies immediate 2D point inversions on the canvas, independent of 3D conformer generation. |
| **`test_properties.py`** | **Analysis Calculations** | **Not RDKit**: Verifies error handling wrappers and data formatting for the analysis window. |
| **`test_parsers.py`** | **Data Ingestion** | **Not RDKit**: Verifies custom XYZ/MOL parsing logic that handles non-standard or corrupt files gracefully. |
| **`test_io.py`** | **Serialization** | **Not RDKit**: Verifies the proprietary `.pmeprj` JSON structure and round-trip integrity. |

## Mocking Strategy
These tests use a `MockMainWindow` (defined in `conftest.py`) which simulates the host environment with standardized `CalculationWorker`, `MolecularData`, and logic implementations from `main_window_app_state.py`.
