# -*- coding: utf-8 -*-
"""Conftest for the full-GUI tier: real window, real VTK, real worker thread.

See README.md for what this tier covers and how to run it.
"""

import os
import sys
import importlib
import importlib.util

import pytest

# Package selection mirrors tests/e2e/conftest.py.
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

# offscreen leaves the interactor without a window handle, and
# interactor.Initialize() then crashes the process. CI supplies xvfb instead.
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

    Runs in a subprocess: no display aborts VTK outright rather than raising.
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
    """Skip the tier without a usable OpenGL context.

    MOLEDITPY_FULL_GUI_REQUIRED=1 makes that a failure instead; the Linux job
    sets it, where xvfb guarantees a context and a skip means CI broke.
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

    # Teardown emits a wall of wglMakeCurrent/glX errors that drown the summary.
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

    # Any modal that slips through would block the run forever.
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

    # safe_mode keeps the developer's own plugin folder out of the run.
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
    """Disarm every plotter before pytest-qt pumps the queue in its own hook.

    pyvistaqt decorates render with @threaded on macOS, so renders become
    queued signals.  pytest-qt's teardown hook pumps the event queue ahead of
    any fixture finaliser, delivering a render into a plotter that is still
    open.  Python attribute reassignment (_render = lambda) does NOT redirect
    queued metacalls because Qt binds the C++ slot pointer at connect time.

    The only safe defence: finalize the VTK render window *before* pytest-qt
    pumps.  Finalize() detaches the OpenGL context at the C++ level, making
    any subsequent render delivery a harmless no-op inside VTK itself.

    We intentionally do NOT call processEvents() here — doing so would
    deliver the queued metacalls ourselves, which crashes just the same.
    """
    for plotter in list(_live_plotters):
        try:
            plotter._render = lambda *a, **k: None
            plotter.render = lambda *a, **k: None
        except (RuntimeError, AttributeError):
            pass
        try:
            plotter.render_signal.disconnect()
        except (RuntimeError, TypeError, AttributeError):
            pass
        try:
            rw = plotter.render_window
            if rw is not None:
                rw.Finalize()
        except (RuntimeError, AttributeError):
            pass


def _make_rendering_synchronous(win) -> None:
    """Bind pyvistaqt's synchronous `_render` over its deferred `render`.

    `render` is @threaded on macOS, so a draw can land after the window is
    gone. Rendering stays real, it just stops being deferred.
    """
    plotter = getattr(win.view_3d_manager, "plotter", None)
    if plotter is not None and hasattr(plotter, "_render"):
        plotter.render = plotter._render
        _live_plotters.append(plotter)


def _teardown_window(win, app) -> None:
    """Dismantle the window in an order VTK survives.

    Silence rendering first (a queued render crashes inside processEvents),
    then close and drain while the plotter is still alive for the app's own
    callbacks, and only then finalise it.
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

    # UIManager arms singleShot(100, fit_to_view) after a conversion; a test
    # ending inside that window leaves it to fire against a dead plotter.
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

    if plotter is not None and sys.platform != "darwin":
        try:
            plotter.close()
        except (RuntimeError, AttributeError, TypeError):
            pass
        if plotter in _live_plotters:
            _live_plotters.remove(plotter)

    for _ in range(10):
        app.processEvents()


@pytest.hookimpl(trylast=True)
def pytest_sessionfinish(session, exitstatus):
    """Bypass Qt/VTK shutdown on macOS to avoid libqmacstyle Abort trap 6.

    On macOS, Qt unloads style plugins during QApplication destruction and
    then tries to access them, raising SIGABRT (exit code 134).  The tests
    have already run and the exit status is final, so skipping normal cleanup
    via os._exit() is safe.  The guard is macOS-only; other platforms use
    normal teardown.
    """
    if sys.platform == "darwin":
        import sys as _sys
        _sys.stdout.flush()
        _sys.stderr.flush()
        os._exit(int(exitstatus))
