# MoleditPy Test Suite

This directory contains the testing infrastructure for MoleditPy.

## Quick Start: Prerequisites

Before running the tests, ensure you have the required `pytest` plugins installed:

```bash
pip install pytest pytest-qt pytest-cov pytest-timeout
```

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

# Pass custom arguments to pytest (e.g., run tests matching a pattern)
python tests/run_all_tests.py --unit -- -k test_edit

# Skip coverage collection (useful for CI or speed)
python tests/run_all_tests.py --no-cov

# Skip the reporting phase (Markdown/Catalog generation)
python tests/run_all_tests.py --no-report

# Reporting and maintenance
python tests/run_all_tests.py --report-only    # Regenerate all reports without running tests
python tests/run_all_tests.py --catalog-only   # Update only the assertion catalog
```

## Mocking Architecture

To maintain high performance and prevent crashes in headless environments, we use a multi-layered mocking strategy:

### 1. Global "Nuclear Mock" (`tests/conftest.py`)
Applied at the very start of the `pytest` lifecycle. It pre-fills `sys.modules` with mock objects for heavy C++ extensions like **VTK** and **PyVista**. This prevents the real libraries from ever loading in test environments, eliminating GPU-related segmentation faults in CI.

### 2. Standardized App Teardown
All suites use a consistent `app` fixture with platform-aware teardown. It performs aggressive cleanup on Windows/Desktop (to prevent `0x80010108` errors) while staying minimal in CI to avoid instability.

### 3. Local Fixture Mocks
Specific directories (like `tests/gui/`) provide further mocks for UI elements, such as replacing the 3D viewport with a dummy widget and auto-answering dialogs.

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

*   [**Assertion Catalog**](assertion_catalog.md): A comprehensive list of test assertions, mapping test names and descriptions to their core verification logic.
*   [**Coverage Report**](coverage_report.md): A Markdown-formatted summary of the combined coverage (Unit + Integration + GUI) across the entire codebase.

### Running Coverage Locally

You can regenerate the coverage report at any time using the provided utility:

```powershell
# From the project root
python tests/print_cov.py
```
