# Integration Tests

This directory contains integration tests for MoleditPy, focusing on component interactions and real scientific calculations.

## Purpose

While unit tests verify individual components in isolation, integration tests ensure that multiple components work together correctly. Our primary focus here is the **CalculationWorker** and its ability to perform 2D→3D conversions and optimizations using RDKit.

## Key Features

- **Real Calculations**: These tests do *not* mock RDKit's 3D embedding or optimization logic. They execute the full `CalculationWorker.run_calculation` pipeline.
- **Meaningful Geometry Validation**: Instead of hardcoding expected coordinates, these tests validate the *physics* of the result:
  - **Bond Lengths**: Measured against RDKit-computed reference values (MMFF94s).
  - **Bond Angles**: Verified against ideal or reference optimized geometries.
  - **Dihedral Angles**: Checked for convergence to stable staggered conformations (e.g., in n-Butane).
- **Topology-Based Atom Identification**: Since re-parsing MOL blocks can lose atom metadata, these tests use topological identification (connectivity and atomic counts) or `_original_atom_id` to reliably track atoms through the 3D conversion process.

## Running Tests

From the project root:

```bash
# Run only integration tests
python tests/run_all_tests.py --integration

# Or use pytest directly
python -m pytest tests/integration -v
```

## Test Descriptions

| Test File | Component | Verification Scope |
| :--- | :--- | :--- |
| **`test_calculation_worker.py`** | **CalculationWorker Pipeline** | Validates real 3D geometry against physics-based reference values.<br>• **Bond Lengths**: Measured against RDKit MMFF94s reference values.<br>• **Bond Angles**: Verified against ideal or reference optimized geometries.<br>• **Dihedral Angles**: Checked for convergence to stable staggered conformations (e.g., n-Butane). |
| **`test_chiral_labels_integration.py`** | **Chiral Label Display** | Integration tests for the chiral label lifecycle in the 3D viewer.<br>• **Show/Hide Toggle**: Verifies that "Show Chiral Labels" correctly creates and removes R/S label actors.<br>• **Mirror Inversion**: Confirms that R/S assignments are correctly inverted after a mirror transformation. |

## Requirements

- **PyQt6**: Required for the `CalculationWorker` signals and event loop.
- **RDKit**: Required for the 3D structure generation and geometry validation.
- **pytest-qt**: Used to wait for worker signals (`qtbot.waitSignal`).
