# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MoleditPy is a Python-based molecular editor (v3.6.0) built with PyQt6, RDKit, and PyVista. It provides 2D editing + 3D visualization of molecules with a plugin system.

## Development Setup

Install for development from the repo root:
```bash
pip install -e moleditpy/
```

Install test dependencies:
```bash
pip install pytest pytest-qt pytest-cov pytest-timeout
```

## Running Tests

Run the full test suite (mirrors CI):
```bash
MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python tests/run_all_tests.py --no-cov --no-report --unit --integration
```

Run a single test file:
```bash
MOLEDITPY_HEADLESS=1 python tests/run_all_tests.py --unit --no-cov --no-report -- tests/unit/test_molecular_data.py
```

Run a single test by name:
```bash
MOLEDITPY_HEADLESS=1 python tests/run_all_tests.py --unit --no-cov --no-report -- -k test_function_name
```

Run with coverage:
```bash
MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python tests/run_all_tests.py --unit --integration
```

Run only unit tests or only integration tests:
```bash
python tests/run_all_tests.py --unit
python tests/run_all_tests.py --integration
```

Speed up dev testing (fail fast / exit on first failure, run each suite at most once, no retries):
```bash
python tests/run_all_tests.py --unit --no-cov --no-report -x
```

Run full suite + pylint in one pass (updates `tests/pylint-score.txt` automatically):
```bash
MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python tests/run_all_tests.py --pylint
```

## Linting

```bash
pylint moleditpy/src/moleditpy/
```

Target score: > 9.0/10. PEP 8 compliance required. Type hints required for all functions and methods.

The `--pylint` flag on `run_all_tests.py` runs pylint automatically after all tests pass and writes the score to `tests/pylint-score.txt`.

## Scope

All code changes go inside `moleditpy/` only. Never edit `moleditpy-linux/` — it is synced from `moleditpy/` separately.

## Code Quality Standards

- Pylint > 9.0/10, PEP 8
- Type hints required for all functions and methods
- **Full Application coverage around 80%** (tracked in `tests/coverage_report.md`)
- All new public methods in `plugins/` must have tests
- Error handling: use `logging.exception()` or `logging.error()` in `except` blocks — never `print`, never bare `pass`
- No unnecessary comments — only add comments where logic is genuinely non-obvious

## Architecture

The package source lives at `moleditpy/src/moleditpy/`. The entry point is `__main__.py` → `main.py`.

### Composition-Over-Inheritance Pattern

`MainWindow` (`ui/main_window.py`) is a thin composition hub. It instantiates specialized manager objects and delegates all work to them. Each manager holds `self.host` pointing back to `MainWindow`.

Key managers:
- `StateManager` (`app_state.py`) — centralized undo/redo via serialized state snapshots
- `IOManager` / `io_logic.py` — file I/O (MOL, SDF, SMILES, etc.)
- `View3DManager` / `view_3d_logic.py` — PyVista 3D visualization
- `DialogManager` / `dialog_logic.py` — dialog lifecycle
- `CalculationWorker` (`calculation_worker.py`) — heavy computation on a `QThread`

### Layer Separation

| Layer | Modules |
|-------|---------|
| Core data model (pure Python, no Qt) | `core/molecular_data.py`, `core/mol_geometry.py` |
| Logic (Qt-aware but not UI) | `*_logic.py`, `app_state.py`, `calculation_worker.py` |
| UI components | `molecule_scene.py`, `atom_item.py`, `bond_item.py`, dialogs |
| Plugin system | `plugins/` |
| Utilities | `utils/constants.py`, `utils/default_settings.py`, `utils/system_utils.py` |

### Plugin System

Plugins are discovered and loaded by `plugins/plugin_manager.py`. The public API for plugin authors is `plugins/plugin_interface.py` (`PluginContext`). Safe mode (`--safe` flag) skips plugin loading.

### VTK/PyVista Mocking in Tests

`tests/conftest.py` aggressively mocks VTK and PyVista at import time so unit/integration tests run headless without GPU. Do not import VTK/PyVista in core logic — keep 3D rendering isolated to `view_3d_logic.py`, `custom_interactor_style.py`, and `custom_qt_interactor.py`.

## Environment Variables

| Variable | Purpose |
|----------|---------|
| `MOLEDITPY_HEADLESS=1` | Disables GUI features that require a display |
| `QT_QPA_PLATFORM=offscreen` | Qt headless rendering for CI/test |

## Running the Application

```bash
moleditpy                  # launch GUI
moleditpy molecule.mol     # open file on startup
moleditpy --safe           # skip plugin loading
```

## Key Documentation

- `docs/ARCHITECTURE.md` — detailed design decisions and manager patterns
- `docs/CORE_DATA_STRUCTURES.md` — molecular data model internals
- `docs/PLUGIN_DEVELOPMENT_MANUAL_V4.md` — plugin API for authors (current)
- `docs/old/PLUGIN_DEVELOPMENT_MANUAL_V3.md` — plugin API for authors (v3, archived)
- `DESIGN_PRINCIPLES.md` — project philosophy
- `CONTRIBUTING.md` — full contribution workflow and PR checklist
