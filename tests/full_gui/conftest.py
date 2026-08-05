# -*- coding: utf-8 -*-
"""Conftest for the full-GUI tier.

Unlike ``tests/gui`` and ``tests/e2e``, nothing is mocked here: the real
PyVista/VTK ``CustomQtInteractor`` is embedded, the real ``CalculationWorker``
QThread performs 2D->3D conversion, and the window is actually ``show()``-n.
The only concessions to CI are a sandboxed home directory, safe mode (no
plugins), and neutralised modal dialogs.

If VTK cannot create a rendering context (no GPU *and* no software GL), every
test in this tier skips rather than fails -- see :func:`_probe_render_context`.
"""

import os
import sys
import importlib
import importlib.util

import pytest

# ---------------------------------------------------------------------------
# Package selection: moleditpy_linux on Linux, moleditpy everywhere else.
# Mirrors tests/e2e/conftest.py so the same suite covers both variants.
# ---------------------------------------------------------------------------
_IS_LINUX = sys.platform.startswith("linux")

_LINUX_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy-linux", "src")
)
_MAIN_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)

if os.environ.get("MOLEDITPY_USE_INSTALLED") == "1":
    _PKG = (
        "moleditpy_linux"
        if importlib.util.find_spec("moleditpy_linux") is not None
        else "moleditpy"
    )
elif _IS_LINUX and os.path.isdir(_LINUX_SRC):
    if _LINUX_SRC not in sys.path:
        sys.path.insert(0, _LINUX_SRC)
    _PKG = "moleditpy_linux"
else:
    if _MAIN_SRC not in sys.path:
        sys.path.insert(0, _MAIN_SRC)
    _PKG = "moleditpy"

PKG = _PKG

# A *real* window is the point of this tier. QT_QPA_PLATFORM=offscreen gives the
# embedded QVTKRenderWindowInteractor no native window handle, and
# `interactor.Initialize()` then dereferences it -- an access violation, not an
# exception. So the platform plugin is forced back to the native default
# (windows/cocoa/xcb) even when the surrounding CI job exported "offscreen";
# the Linux job supplies the display via xvfb-run.
os.environ.pop("QT_QPA_PLATFORM", None)
os.environ.pop("PYVISTA_OFF_SCREEN", None)


_render_probe_error: "str | None" = None
_render_probe_done = False

_PROBE_CODE = """
import os, sys
os.environ.pop("QT_QPA_PLATFORM", None)
from PyQt6.QtWidgets import QApplication
app = QApplication([sys.argv[0]])
from pyvistaqt import QtInteractor
w = QtInteractor()
w.resize(320, 240)
w.show()
w.add_mesh(__import__("pyvista").Sphere())
w.render()
print("PROBE_OK")
"""


def _probe_render_context() -> "str | None":
    """Return None if a real VTK-in-Qt window works here, else why it does not.

    A runner with no display (or no software GL) aborts the process inside
    VTK rather than raising, so the probe runs once, up front, in a subprocess
    where a crash is just a non-zero exit code.
    """
    global _render_probe_done, _render_probe_error
    if _render_probe_done:
        return _render_probe_error

    _render_probe_done = True
    import subprocess

    try:
        res = subprocess.run(
            [sys.executable, "-c", _PROBE_CODE],
            capture_output=True,
            text=True,
            timeout=300,
            env={**os.environ},
        )
    except (OSError, subprocess.SubprocessError) as exc:
        _render_probe_error = f"VTK/Qt render probe could not run: {exc}"
        return _render_probe_error

    if res.returncode != 0 or "PROBE_OK" not in res.stdout:
        tail = (res.stderr or res.stdout or "").strip().splitlines()[-3:]
        _render_probe_error = "no usable VTK-in-Qt render context: " + " | ".join(tail)
    return _render_probe_error


@pytest.fixture(scope="session", autouse=True)
def require_render_context():
    """Skip the whole tier when there is no usable OpenGL context.

    Set ``MOLEDITPY_FULL_GUI_REQUIRED=1`` to turn that skip into a failure.
    The Linux CI job sets it: there, xvfb plus software GL is guaranteed, so a
    skip would mean the environment broke and the tier silently stopped
    testing anything -- exactly the outcome this tier exists to prevent.
    """
    reason = _probe_render_context()
    if not reason:
        return
    if os.environ.get("MOLEDITPY_FULL_GUI_REQUIRED") == "1":
        pytest.fail(
            "MOLEDITPY_FULL_GUI_REQUIRED=1 but the full-GUI tier cannot run: " + reason,
            pytrace=False,
        )
    pytest.skip(reason)


@pytest.fixture(scope="session")
def app(require_render_context):
    from PyQt6.QtWidgets import QApplication

    # Tearing a QVTKRenderWindowInteractor down emits a wall of wglMakeCurrent /
    # glX errors that drown the pytest summary. They are teardown-only and do
    # not affect results.
    try:
        import vtk

        vtk.vtkObject.GlobalWarningDisplayOff()
    except (ImportError, AttributeError):
        pass

    q_app = QApplication.instance() or QApplication([sys.argv[0]])
    yield q_app


@pytest.fixture
def full_window(app, qtbot, monkeypatch, tmp_path):
    """A really-launched MainWindow: real VTK widget, real worker thread, shown."""
    from PyQt6.QtWidgets import QMessageBox, QFileDialog

    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setenv("USERPROFILE", str(fake_home))

    # Any modal that slips through would block the run forever under offscreen.
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes),
    )
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: ("", ""))
    )
    monkeypatch.setattr(
        QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: ("", ""))
    )

    MainWindow = importlib.import_module(f"{PKG}.ui.main_window").MainWindow

    # safe_mode=True: the user's real plugin folder must not influence CI.
    win = MainWindow(initial_file=None, safe_mode=True)
    qtbot.addWidget(win)
    win.resize(1280, 800)
    win.show()
    qtbot.waitExposed(win, timeout=30000)

    try:
        yield win
    finally:
        _teardown_window(win, app)


def _teardown_window(win, app) -> None:
    """Dismantle the window in an order VTK survives.

    Closing the QWidget alone is not enough. Two things outlive it and then
    touch a render window whose GL context is gone -- a segfault, not an
    exception, and it lands in whichever test runs next:

    * ``UIManager._style_watchdog``, a 2 s QTimer that reads
      ``plotter.interactor`` (harmless in the real app, where the window
      outlives the process, fatal across a test boundary);
    * pyvistaqt's own render path -- ``QtInteractor.render`` emits a signal
      handled later, and on macOS that handler runs on another thread.

    So: stop the timer, finalise the plotter, then drain the event queue.
    """
    plotter = getattr(win.view_3d_manager, "plotter", None)

    watchdog = getattr(win.ui_manager, "_style_watchdog", None)
    if watchdog is not None:
        try:
            watchdog.stop()
        except RuntimeError:
            pass

    win.state_manager.has_unsaved_changes = False
    win.closeEvent = lambda e: e.accept()
    win.close()

    # Drain first: the app still has queued callbacks (view_isometric,
    # deferred redraws) that need a live plotter. Finalising before they run
    # gives them a camera of None.
    for _ in range(10):
        app.processEvents()

    if plotter is not None:
        # Silence the render path before finalising. pyvistaqt's `render` emits
        # a signal whose handler is decorated `@threaded` on macOS, so closing
        # can schedule a render on another thread against a window that is
        # already going away -- segfault, inside plotter.close() itself.
        try:
            plotter.render_signal.disconnect()
        except (RuntimeError, TypeError, AttributeError):
            pass
        try:
            plotter.render = lambda *a, **k: None
        except (RuntimeError, AttributeError):
            pass
        try:
            plotter.close()
        except (RuntimeError, AttributeError, TypeError):
            pass

    for _ in range(10):
        app.processEvents()
