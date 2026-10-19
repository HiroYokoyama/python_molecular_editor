# E2E tests

These tests exercise the unified `moleditpy` package on Windows, Linux, and
macOS, including RDKit conversion, molecule geometry and application workflows.
They use real RDKit chemistry and real PySide6 widgets where a window is needed.
VTK/PyVista are mocked by `conftest.py`, so this tier runs headlessly without
an OpenGL context. It does not test an Open Babel backend.

For real native rendering and the background calculation thread, use the
separate [full-GUI tier](../full_gui/README.md).

## Package selection

The default uses `moleditpy/src/moleditpy/` on every platform. Set
`MOLEDITPY_USE_INSTALLED=1` to validate a pip-installed package instead of the
local source. There is no separate Linux import namespace.

## Running

```sh
python tests/run_all_tests.py --e2e
python -m pytest tests/e2e/ -v
```

Use `QT_QPA_PLATFORM=offscreen` for headless execution and
`PYTEST_QT_API=pyside6` when multiple Qt bindings are installed.
