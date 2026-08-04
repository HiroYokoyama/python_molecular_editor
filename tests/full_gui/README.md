# Full-GUI tests

The fourth tier (`tests/gui`) runs the real `MainWindow` with VTK and PyVista
**mocked**. This tier mocks none of it:

| | `tests/gui`, `tests/e2e` | `tests/full_gui` |
|---|---|---|
| Qt platform | `offscreen` | native (`windows` / `cocoa` / `xcb` under xvfb) |
| 3D widget | dummy `QWidget` | real `CustomQtInteractor` (`pyvistaqt` + VTK) |
| Window | constructed only | `show()`-n and waited on with `waitExposed` |
| 2D→3D conversion | monkeypatched to run synchronously | real `CalculationWorker` on its own `QThread` |
| Plugins | dummy manager | real manager, started in safe mode |

What it covers: the application launches, benzene is drawn through the actual
toolbar/mouse path, and the 2D→3D conversion produces chemically correct
geometry that reaches the real renderer.

## Running it

```bash
# Linux (no monitor): xvfb supplies a real X display, mesa supplies software GL
xvfb-run -a --server-args="-screen 0 1600x1200x24" \
  python -m pytest tests/full_gui -c tests/full_gui/pytest.ini -v

# Windows / macOS, on a machine with a desktop
python -m pytest tests/full_gui -c tests/full_gui/pytest.ini -v

# through the standard runner (never part of a no-flag "run everything")
python tests/run_all_tests.py --full-gui --no-cov --no-report
```

Do **not** export `QT_QPA_PLATFORM=offscreen` for this tier. `conftest.py`
unsets it anyway: with the offscreen plugin the embedded
`QVTKRenderWindowInteractor` has no native window handle, and
`ui_manager._setup_3d_picker`'s `interactor.Initialize()` then crashes the
process with an access violation instead of raising.

## No render context?

A session-scoped probe launches a throwaway subprocess that puts a VTK sphere
in a Qt window. If that subprocess dies — no display, no usable OpenGL — the
whole tier **skips** with the reason attached, because a missing GPU is an
environment fact, not a defect in the application.

Set `MOLEDITPY_FULL_GUI_REQUIRED=1` to turn that skip into a failure. CI's
Linux job sets it, so a broken runner image is reported instead of quietly
disabling every test in this directory.

## Environment variables

| Variable | Default | Effect |
|---|---|---|
| `MOLEDITPY_FULL_GUI_REQUIRED` | unset | `1` = no render context is a failure, not a skip |
| `MOLEDITPY_LAUNCH_TIMEOUT` | `180` | seconds the subprocess launch test waits before killing the app and failing |

`MOLEDITPY_LAUNCH_TIMEOUT` exists because the failure mode worth catching is a
launch that *never finishes* — reported on macOS with the PyQt6 range pinned
for Python 3.9. `test_entry_point_boots_in_a_subprocess` prints staged
markers (`STAGE_QT_IMPORTED`, `STAGE_APP_IMPORTED`,
`STAGE_WINDOW_CONSTRUCTED`, `STAGE_SHOW_CALLED`), so a timeout reports the
last stage reached rather than just "it hung".
