# MoleditPy Test Suite

This directory contains the testing infrastructure for MoleditPy.

## Test Architecture

The MoleditPy test suite is divided into two primary categories to ensure both scientific accuracy and user interface robustness.

### 1. Unit Tests (`tests/unit/`)
**Scope:** Core scientific logic, data structure integrity, and backend algorithms.
These tests run in isolation and verify that individual components—such as the molecular data model, hydrogen management, and file parsers—function correctly without requiring a running GUI instance.
*   [**View Detailed Unit Test Documentation**](unit/README.md)

### 2. Integration Tests (`tests/integration/`)
**Scope:** Real 2D→3D conversion, physics-based geometry validation, and component interaction.
These tests perform real calculations using the `CalculationWorker` without mocking, validating molecular coordinates (bonds, angles, dihedrals) against RDKit-computed physical references.
*   [**View Detailed Integration Test Documentation**](integration/README.md)

### 3. GUI Tests (`tests/gui/`)
**Scope:** Application integration, user workflow verification, and event handling.
These tests simulate real user interactions (clicks, key presses, drag-and-drop) in a headless environment to verify that the UI correctly triggers the underlying scientific logic and updates the visual state.
*   [**View Detailed GUI Test Documentation**](gui/README.md)

### Unified Runner (`run_all_tests.py`)

For convenience, a unified runner script is provided to execute **Unit**, **Integration**, and **GUI** tests in a single pass, ensuring the environment is correctly set up.

```bash
# Run all tests (Unit + Integration + GUI)
python tests/run_all_tests.py

# Run all tests in HEADLESS mode (no windows will pop up)
python tests/run_all_tests.py --headless

# Run only a specific suite
python tests/run_all_tests.py --unit
python tests/run_all_tests.py --integration
python tests/run_all_tests.py --gui --headless
```

This script:
1.  Sets `PYTHONPATH` to prioritize the local `src` directory.
2.  Runs `tests/unit` first.
3.  Runs `tests/integration` second.
4.  Runs `tests/gui` third (only if previous suites pass).
