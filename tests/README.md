# MoleditPy Test Suite

This directory contains the testing infrastructure for MoleditPy.

## Quick Start: Running Tests

The recommended way to run tests is using the **Unified Runner**, which handles environment setup and executes suites in the correct order (Unit -> Integration -> GUI).

```bash
# Run all tests (Unit + Integration + GUI)
python tests/run_all_tests.py

# Run all tests in HEADLESS mode (recommended for CI or background runs/wsl)
python tests/run_all_tests.py --headless

# Run specific suites
python tests/run_all_tests.py --unit
python tests/run_all_tests.py --integration
python tests/run_all_tests.py --gui --headless
```

## Test Architecture

The MoleditPy test suite is organized into three layers:

### 1. Unit Tests (`tests/unit/`)
**Scope:** Core scientific logic, data structures, and parsers.
*   Fast, isolated tests mocking all external dependencies.
*   [**Read More**](unit/README.md)

### 2. Integration Tests (`tests/integration/`)
**Scope:** Real scientific calculations and component interactions.
*   Validates `CalculationWorker` pipelines and physics-based geometry (bond lengths/angles) using RDKit.
*   [**Read More**](integration/README.md)

### 3. GUI Tests (`tests/gui/`)
**Scope:** Full application workflows and user interactions.
*   Simulates clicks, key presses, and complex workflows (e.g., drawing, undo/redo, saving) in the actual `MainWindow`.
*   [**Read More**](gui/README.md)

## Test Reports & Catalogs

To simplify human review of the test suite, the following generated reports are available:

*   [**Assertion Catalog**](assertion_catalog.md): A comprehensive list of all 91 test assertions, mapping test names and descriptions to their core verification logic.
*   [**Coverage Report**](coverage_report.md): A Markdown-formatted summary of the combined coverage (Unit + Integration + GUI) across the entire codebase.

### Running Coverage Locally

You can regenerate the coverage report at any time using the provided utility:

```powershell
# From the project root
python tests/print_cov.py
```
