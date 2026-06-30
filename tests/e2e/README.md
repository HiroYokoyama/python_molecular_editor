# E2E Tests

This directory contains end-to-end tests that confirm real application behaviour
across all supported platforms (Windows, Linux, macOS).

## Purpose

E2E tests exercise the full conversion pipeline — from drawing a molecule in 2D
to producing a valid optimised 3D structure — using **zero mocks**.  
They run completely headlessly (no display, no Qt window) and are fast enough
to include in every CI run on every OS.

Unlike the integration tests in `tests/integration/`, which focus on individual
component contracts (worker signal protocols, geometry validation), e2e tests
focus on **real package behaviour under the actual OS conditions** the user
experiences:

| Concern | Integration tests | E2E tests |
|---|---|---|
| RDKit geometry accuracy | ✓ (bond lengths, dihedrals) | — |
| Correct package imported per OS | — | ✓ |
| All conversion modes smoke-tested | — | ✓ |
| Open Babel path (Linux) | — | ✓ (auto-skip elsewhere) |
| Headless, no VTK/PyVista mocking | — | ✓ |

## Package Selection

On **Linux**, the `moleditpy-linux` package is used (mirrors the distributed
installer and CI environment).  On **Windows / macOS**, the standard `moleditpy`
package is used.  Selection is automatic based on `sys.platform`.

```
Linux   → moleditpy-linux/src/moleditpy_linux/
Windows → moleditpy/src/moleditpy/
macOS   → moleditpy/src/moleditpy/
```

## Running

```bash
# E2E suite only
python tests/run_all_tests.py --e2e

# Or directly with pytest
python -m pytest tests/e2e/ -v

# As part of the full suite
python tests/run_all_tests.py --unit --integration --e2e
```

No environment variables are required — the conftest sets
`QT_QPA_PLATFORM=offscreen` automatically.

## Test Descriptions

### `test_ethane_conversion.py`

Ethane (2 C, 1 bond) is the minimal molecule that exercises the full
2D → MOL-block → 3D embedding → optimisation pipeline.

| Test | What it checks |
|---|---|
| `test_ethane_2d_structure` | `MolecularData` builds 2 atoms and 1 bond correctly |
| `test_ethane_molblock_roundtrip` | MOL-block produced by `MolecularData.to_mol_block()` is parseable by RDKit |
| `test_ethane_conversion_rdkit` | `CalculationWorker` with `rdkit` mode returns a valid 3D mol |
| `test_ethane_conversion_direct` | `direct` mode returns a valid 3D mol |
| `test_ethane_conversion_fallback` | `fallback` mode (RDKit → obabel chain) returns a valid 3D mol |
| `test_ethane_conversion_obabel` | `obabel` mode returns a valid 3D mol — **skipped** when Open Babel is absent |
| `test_linux_package_is_active` | On Linux, `moleditpy_linux` is the active import — **skipped** on other OSes |
| `test_main_package_is_active` | On Windows/macOS, `moleditpy` is the active import — **skipped** on Linux |

**Geometry assertion:** every conversion test checks that the 3D mol has exactly
2 carbon atoms and a C–C bond length in the chemically reasonable range
`[1.3, 1.7] Å`.

## Requirements

- **RDKit** — 3D embedding and geometry validation
- **PyQt6** — `QApplication` instance required by `CalculationWorker` signals
- **pytest-qt** — `qtbot` fixture (used to satisfy the `app` fixture dependency)
- **Open Babel** *(optional)* — required only for `test_ethane_conversion_obabel`;
  the test is automatically skipped when absent
