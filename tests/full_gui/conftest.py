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
            # Short on purpose: a context that works answers in seconds, and a
            # context that does not should say so quickly rather than stall the
            # whole job.
            timeout=120,
            env={**os.environ},
        )
    except (OSError, subprocess.SubprocessError) as exc:
        _render_probe_error = f"VTK/Qt render probe could not run: {exc}"
        return _render_probe_error

    if res.returncode != 0 or "PROBE_OK" not in res.stdout:
        tail = (res.stderr or res.stdout or "").strip().splitlines()[-12:]
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
    _make_rendering_synchronous(win)
    qtbot.addWidget(win)
    win.resize(1280, 800)
    win.show()
    qtbot.waitExposed(win, timeout=30000)

    try:
        yield win
    finally:
        _teardown_window(win, app)


_live_plotters: list = []


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_teardown(item, nextitem):
    """Disarm every plotter *before* pytest-qt pumps the event queue.

    pytest-qt's own ``pytest_runtest_teardown`` calls ``_process_events()``,
    and it runs ahead of any fixture finaliser. A ``render_signal`` emission
    queued during the test is therefore delivered there, into a plotter this
    fixture has not been given the chance to close yet -- and on macOS that
    delivery segfaults. Disconnecting the signal in the fixture cannot help:
    an already-queued metacall is not undone by disconnecting.

    ``_render`` is the slot on the receiving end, so stubbing it here makes the
    queued delivery a no-op no matter when it arrives.
    """
    for plotter in _live_plotters:
        try:
            plotter._render = lambda *a, **k: None
            plotter.render = lambda *a, **k: None
        except (RuntimeError, AttributeError):
            continue


def _make_rendering_synchronous(win) -> None:
    """Bypass pyvistaqt's deferred render signal for the lifetime of the test.

    ``QtInteractor.render`` normally emits ``render_signal`` and lets the
    handler run later -- and on macOS that handler is decorated ``@threaded``,
    so the render happens on another thread at an unpredictable time. A render
    left in flight when the window goes away is a segfault, and it does not
    even land in our own teardown: pytest-qt pumps the event queue in its
    ``pytest_runtest_teardown`` hook, ahead of any fixture finaliser, so the
    crash appears inside pytest-qt with no way to guard it from here.

    ``_render`` is pyvistaqt's own synchronous path (it is what the signal
    handler calls). Binding it directly keeps rendering real -- the tests still
    drive genuine VTK draws -- while removing the deferral that makes teardown
    unsafe.
    """
    plotter = getattr(win.view_3d_manager, "plotter", None)
    if plotter is not None and hasattr(plotter, "_render"):
        plotter.render = plotter._render
        _live_plotters.append(plotter)


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

    Order matters, and it is not the obvious one:

    1. Stop the watchdog and silence the render path *before* anything else.
       Rendering is what segfaults, and a render can already be sitting in the
       queue -- draining first ran it and crashed inside ``processEvents``.
       Stubbing ``render`` does not break the app's other teardown callbacks;
       ``view_isometric`` and friends still work, they just stop drawing.
    2. Close the window and drain, so those queued callbacks run against a
       plotter that is still alive (finalising first gives them a camera of
       ``None`` and an AttributeError).
    3. Finalise the plotter, then drain again.
    """
    plotter = getattr(win.view_3d_manager, "plotter", None)

    watchdog = getattr(win.ui_manager, "_style_watchdog", None)
    if watchdog is not None:
        try:
            watchdog.stop()
        except RuntimeError:
            pass

    if plotter is not None:
        try:
            plotter.render_signal.disconnect()
        except (RuntimeError, TypeError, AttributeError):
            pass
        try:
            plotter.render = lambda *a, **k: None
        except (RuntimeError, AttributeError):
            pass

    # UIManager schedules `QTimer.singleShot(100, view_3d_manager.fit_to_view)`
    # after a conversion. If the test ends inside that 100 ms the timer fires
    # later, against a manager whose plotter has since been finalised, and
    # `fit_to_view` renders into a dead window. It surfaces in the *next*
    # test's teardown, inside pytest-qt's own event pump, so make the callbacks
    # themselves inert rather than trying to outrun the timer.
    for name in ("fit_to_view", "draw_molecule_3d", "update_3d_view"):
        if hasattr(win.view_3d_manager, name):
            try:
                setattr(win.view_3d_manager, name, lambda *a, **k: None)
            except (AttributeError, RuntimeError):
                pass

    win.state_manager.has_unsaved_changes = False
    win.closeEvent = lambda e: e.accept()
    win.close()

    for _ in range(10):
        app.processEvents()

    if plotter is not None:
        try:
            plotter.close()
        except (RuntimeError, AttributeError, TypeError):
            pass
        if plotter in _live_plotters:
            _live_plotters.remove(plotter)

    for _ in range(10):
        app.processEvents()
