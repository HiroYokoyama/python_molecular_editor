# -*- coding: utf-8 -*-
"""Full application launch: real window, real VTK widget, real plugin-free boot."""

import os
import subprocess
import sys
import time

import pytest

from conftest import PKG

pytestmark = pytest.mark.full_gui

# A launch that never completes is the failure mode we care about most: on
# macOS the PyQt6 version pinned for Python 3.9 has been reported to hang
# instead of showing a window. The budget is generous enough for a cold CI
# runner (RDKit/VTK imports + font cache) and is overridable per platform.
LAUNCH_TIMEOUT_S = float(os.environ.get("MOLEDITPY_LAUNCH_TIMEOUT", "180"))


def test_window_is_shown_and_exposed(full_window):
    """The window really becomes visible (not just constructed)."""
    assert full_window.isVisible()
    assert full_window.windowTitle()


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
    """Safe mode leaves no plugin manager at all, so CI is deterministic.

    `MainInitManager.__init__` sets `plugin_manager = None` under safe mode; if
    that ever becomes a real manager, the developer's own ~/.moleditpy/plugins
    would start influencing this tier.
    """
    assert full_window.plugin_manager is None, (
        f"safe mode built a plugin manager: {full_window.plugin_manager!r}"
    )


def test_entry_point_boots_in_a_subprocess(tmp_path):
    """`MainWindow` + `show()` via the installed package, in a clean process.

    This is the only test that exercises a cold interpreter, so it catches
    import-order and top-level-side-effect regressions the in-process fixture
    cannot see.
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
    code = f"""
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
MainWindow = importlib.import_module("{PKG}.ui.main_window").MainWindow
mark("STAGE_APP_IMPORTED")

app = QApplication([sys.argv[0]])
win = MainWindow(initial_file=None, safe_mode=True)
mark("STAGE_WINDOW_CONSTRUCTED")
win.show()
mark("STAGE_SHOW_CALLED")

def _check():
    # Reported from inside the live event loop. Printing after app.exec()
    # returns would put the verdict behind Qt/VTK shutdown, which can crash on
    # some platforms and would then hide a launch that actually succeeded.
    print(
        "LAUNCH_OK",
        win.isVisible(),
        win.view_3d_manager.plotter is not None,
        flush=True,
    )
    app.quit()


QTimer.singleShot(3000, _check)
app.exec()
print("EXIT_REACHED", flush=True)
os._exit(0)
"""
    env = dict(os.environ)
    env["HOME"] = str(tmp_path)
    env["USERPROFILE"] = str(tmp_path)
    env.pop("QT_QPA_PLATFORM", None)
    env.pop("PYVISTA_OFF_SCREEN", None)

    started = time.monotonic()
    proc = subprocess.Popen(
        [sys.executable, "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    try:
        stdout, stderr = proc.communicate(timeout=LAUNCH_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, stderr = proc.communicate()
        stages = [ln for ln in stdout.splitlines() if ln.startswith("STAGE_")]
        pytest.fail(
            f"the GUI never finished launching within {LAUNCH_TIMEOUT_S}s "
            f"(the process was killed, not left hanging).\n"
            f"last stage reached: {stages[-1] if stages else '<none>'}\n"
            f"stages: {stages}\n"
            f"--- stderr ---\n{stderr[-4000:]}"
        )

    elapsed = time.monotonic() - started
    assert "LAUNCH_OK True True" in stdout, (
        f"launch failed after {elapsed:.1f}s (rc={proc.returncode})\n"
        f"--- stdout ---\n{stdout}\n--- stderr ---\n{stderr[-4000:]}"
    )
