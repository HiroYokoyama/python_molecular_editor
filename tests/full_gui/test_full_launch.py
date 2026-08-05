# -*- coding: utf-8 -*-
"""Full application launch: real window, real VTK widget, real plugin-free boot."""

import os
import subprocess
import sys
import time

import pytest

from conftest import PKG

pytestmark = pytest.mark.full_gui

# A launch that never completes is the failure mode this file exists for.
LAUNCH_TIMEOUT_S = float(os.environ.get("MOLEDITPY_LAUNCH_TIMEOUT", "180"))


def test_window_is_shown_and_exposed(full_window):
    """The window really becomes visible (not just constructed)."""
    assert full_window.isVisible()
    assert full_window.windowTitle()
    handle = full_window.windowHandle()
    assert handle is not None, "shown window has no native handle"
    assert handle.isExposed(), "window was shown but never mapped by the compositor"


def test_core_managers_exist(full_window):
    """Every manager MainWindow composes is live after a real launch."""
    for name in (
        "state_manager",
        "init_manager",
        "ui_manager",
        "io_manager",
        "compute_manager",
        "view_3d_manager",
        "edit_actions_manager",
        "dialog_manager",
    ):
        assert getattr(full_window, name) is not None, f"{name} missing"


def test_real_vtk_interactor_is_embedded(full_window):
    """The 3D pane is the genuine PyVista widget, not a stand-in."""
    from pyvistaqt import QtInteractor

    plotter = full_window.view_3d_manager.plotter
    assert isinstance(plotter, QtInteractor), f"got {type(plotter).__name__}"
    assert plotter.renderer is not None
    # A real render round-trip must not raise.
    plotter.render()


def test_2d_view_and_scene_are_live(full_window):
    scene = full_window.init_manager.scene
    view = full_window.init_manager.view_2d
    assert scene is not None
    assert view is not None
    assert scene.views(), "scene has no attached view"


def test_safe_mode_loads_no_plugins(full_window):
    """Safe mode leaves no plugin manager, so ~/.moleditpy/plugins cannot leak in."""
    assert full_window.plugin_manager is None, (
        f"safe mode built a plugin manager: {full_window.plugin_manager!r}"
    )


def _child_source(src_root: str) -> str:
    """Source for the cold-interpreter launch, driven through the real `main()`.

    Nothing here reimplements the boot: `moleditpy.main.main()` runs verbatim,
    so setup_logging, the excepthook, argparse, qInstallMessageHandler and the
    win32 AppUserModelID call are all on the path a user's `moleditpy` takes.
    MainWindow is subclassed only to emit stage markers and to arm the probe,
    which happens from `show()` — after main() has built the QApplication.
    """
    return f"""
import os, sys
os.environ.pop("QT_QPA_PLATFORM", None)
if os.path.isdir({src_root!r}):
    sys.path.insert(0, {src_root!r})


def mark(name):
    print(name, flush=True)


import importlib
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer

mark("STAGE_QT_IMPORTED")
main_mod = importlib.import_module("{PKG}.main")
mark("STAGE_APP_IMPORTED")

EXPOSE_DEADLINE_MS = 30000
POLL_MS = 100
_elapsed = [0]
_entered = [False]


def _verdict(win):
    handle = win.windowHandle()
    exposed = handle is not None and handle.isExposed()
    plotter = getattr(win.view_3d_manager, "plotter", None)
    rendered = False
    if exposed and plotter is not None:
        # A real round trip through the GL pipeline, not just "the object exists".
        plotter.render()
        rendered = True
    return win.isVisible(), exposed, rendered


def _poll(win):
    if not _entered[0]:
        _entered[0] = True
        mark("STAGE_EVENT_LOOP_ENTERED")
    visible, exposed, rendered = _verdict(win)
    if visible and exposed and rendered:
        # Printed from inside the live loop: after app.exec() the verdict would
        # sit behind Qt/VTK shutdown, which can crash and hide a good launch.
        print("LAUNCH_OK visible=True exposed=True render=True", flush=True)
        QApplication.instance().quit()
        return
    _elapsed[0] += POLL_MS
    if _elapsed[0] >= EXPOSE_DEADLINE_MS:
        print(
            "LAUNCH_INCOMPLETE visible=%s exposed=%s render=%s"
            % (visible, exposed, rendered),
            flush=True,
        )
        QApplication.instance().quit()
        return
    QTimer.singleShot(POLL_MS, lambda: _poll(win))


class _ProbeWindow(main_mod.MainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        mark("STAGE_WINDOW_CONSTRUCTED")

    def show(self):
        super().show()
        mark("STAGE_SHOW_CALLED")
        QTimer.singleShot(POLL_MS, lambda: _poll(self))


main_mod.MainWindow = _ProbeWindow

sys.argv = ["moleditpy", "--safe"]
try:
    main_mod.main()
except SystemExit:
    pass
print("EXIT_REACHED", flush=True)
os._exit(0)
"""


@pytest.mark.timeout(LAUNCH_TIMEOUT_S + 60)
def test_entry_point_boots_in_a_subprocess(tmp_path):
    """Run the shipped entry point in a cold interpreter and watch it come up.

    A launch that never finishes is the failure this exists for, so the child
    is killed and reported rather than left to hang the job.
    """
    src_root = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "moleditpy-linux" if PKG == "moleditpy_linux" else "moleditpy",
            "src",
        )
    )
    env = dict(os.environ)
    env["HOME"] = str(tmp_path)
    env["USERPROFILE"] = str(tmp_path)
    env.pop("QT_QPA_PLATFORM", None)
    env.pop("PYVISTA_OFF_SCREEN", None)

    started = time.monotonic()
    proc = subprocess.Popen(
        [sys.executable, "-c", _child_source(src_root)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    try:
        stdout, stderr = proc.communicate(timeout=LAUNCH_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        proc.kill()
        try:
            # kill() reaps only the direct child; a surviving grandchild holding
            # the pipe would otherwise block here with no bound.
            stdout, stderr = proc.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            stdout, stderr = "", "<pipes never closed after kill>"
        stages = [ln for ln in stdout.splitlines() if ln.startswith("STAGE_")]
        pytest.fail(
            f"the GUI never finished launching within {LAUNCH_TIMEOUT_S}s "
            f"(the process was killed, not left hanging).\n"
            f"last stage reached: {stages[-1] if stages else '<none>'}\n"
            f"stages: {stages}\n"
            f"--- stderr ---\n{stderr[-4000:]}"
        )

    elapsed = time.monotonic() - started
    assert "LAUNCH_OK visible=True exposed=True render=True" in stdout, (
        f"launch failed after {elapsed:.1f}s (rc={proc.returncode})\n"
        f"--- stdout ---\n{stdout}\n--- stderr ---\n{stderr[-4000:]}"
    )
