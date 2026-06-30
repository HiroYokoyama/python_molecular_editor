# -*- coding: utf-8 -*-
"""Conftest for e2e tests.

VTK/PyVista are mocked at import time so the real MainWindow can be
instantiated headlessly.  The `window` fixture selects the right package
(moleditpy_linux on Linux, moleditpy everywhere else) automatically.
"""

import os
import sys
import types
import pytest
from unittest import mock as _mock
from PyQt6.QtWidgets import QApplication

# ---------------------------------------------------------------------------
# Package selection (mirrors test_ethane_conversion.py)
# ---------------------------------------------------------------------------
_IS_LINUX = sys.platform.startswith("linux")

_LINUX_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy-linux", "src")
)
_MAIN_SRC = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)

if _IS_LINUX and os.path.isdir(_LINUX_SRC):
    if _LINUX_SRC not in sys.path:
        sys.path.insert(0, _LINUX_SRC)
    _PKG = "moleditpy_linux"
else:
    if _MAIN_SRC not in sys.path:
        sys.path.insert(0, _MAIN_SRC)
    _PKG = "moleditpy"

# ---------------------------------------------------------------------------
# VTK / PyVista mocking — must happen before any app module is imported
# ---------------------------------------------------------------------------

def _force_mock(name: str) -> types.ModuleType:
    mod = types.ModuleType(name)
    mod.__path__ = []
    sys.modules[name] = mod
    return mod


vtk = _force_mock("vtk")


class _DummyCellPicker:
    def SetTolerance(self, v): pass
    def Pick(self, x, y, z, r): pass
    def GetActor(self): return None


class _DummyAxesWidget:
    def SetOrientationMarker(self, m): pass
    def SetInteractor(self, i): pass
    def SetViewport(self, *a): pass
    def On(self): pass
    def Off(self): pass
    def SetInteractive(self, v): pass


vtk.vtkCellPicker = _DummyCellPicker
vtk.vtkAxesActor = _mock.MagicMock
vtk.vtkOrientationMarkerWidget = _DummyAxesWidget
vtk.vtkInteractorStyleTrackballCamera = _mock.MagicMock

vmod = _force_mock("vtkmodules")
_force_mock("vtkmodules.all")
_force_mock("vtkmodules.vtkRenderingCore")
_force_mock("vtkmodules.vtkCommonCore")
vinter = _force_mock("vtkmodules.vtkInteractionStyle")


class _DummyInteractorStyle:
    def __init__(self, *a, **k): pass
    def AddObserver(self, *a, **k): pass
    def GetInteractor(self): return _mock.MagicMock()
    def OnLeftButtonDown(self, *a): pass
    def OnLeftButtonUp(self, *a): pass
    def OnRightButtonDown(self, *a): pass
    def OnRightButtonUp(self, *a): pass
    def OnMouseMove(self, *a): pass


vinter.vtkInteractorStyleTrackballCamera = _DummyInteractorStyle

pyv = _force_mock("pyvista")
pvqt = _force_mock("pyvistaqt")


class _DummyPlotterBase:
    def __init__(self, *a, **k):
        self.camera_position = [(0, 0, 5), (0, 0, 0), (0, 1, 0)]
    def add_mesh(self, *a, **k): return _mock.MagicMock()
    def add_point_labels(self, *a, **k): return []
    def clear(self): pass
    def render(self): pass
    def reset_camera(self): pass
    def add_light(self, *a, **k): return _mock.MagicMock()
    def remove_actor(self, *a, **k): return True
    @property
    def camera(self):
        c = _mock.MagicMock()
        c.copy = _mock.MagicMock(return_value={})
        return c


class _DummyPolyData(_mock.MagicMock):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.point_data = {}
        self.cell_data = {}
    def __getitem__(self, k): return _mock.MagicMock()
    def __setitem__(self, k, v): pass
    def glyph(self, *a, **k): return _mock.MagicMock()
    def tube(self, *a, **k): return _mock.MagicMock()


pyv.Plotter = _DummyPlotterBase
pyv.PolyData = _DummyPolyData
pyv.Light = _mock.MagicMock
pyv.Sphere = _mock.MagicMock
pyv.Spline = _mock.MagicMock
pyv.Cylinder = _mock.MagicMock
pyv.Box = _mock.MagicMock
pvqt.BackgroundPlotter = _DummyPlotterBase

# QtInteractor must be a QWidget subclass so it can be embedded in the splitter
from PyQt6.QtWidgets import QWidget as _QWidget


class _DummyQtInteractor(_QWidget):
    def __init__(self, parent=None, *a, **k):
        super().__init__(parent)
        self.renderer = _mock.MagicMock()
        self.add_mesh = _mock.MagicMock(return_value=_mock.MagicMock())
        self.add_point_labels = _mock.MagicMock(return_value=[])
        self.clear = _mock.MagicMock()
        self.render = _mock.MagicMock()
        self.reset_camera = _mock.MagicMock()
        self.set_background = _mock.MagicMock()
        self.interactor = _mock.MagicMock()
        self.picker = _mock.MagicMock()
        self.camera = _mock.MagicMock()
        self.camera.copy = _mock.MagicMock(return_value={})
        self.camera_position = [(0, 0, 5), (0, 0, 0), (0, 1, 0)]
        self.add_light = _mock.MagicMock(return_value=_mock.MagicMock())
        self.remove_actor = _mock.MagicMock()


pvqt.QtInteractor = _DummyQtInteractor

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    q_app = QApplication.instance() or QApplication(sys.argv)
    yield q_app
    q_app.quit()


@pytest.fixture
def window(app, qtbot, monkeypatch):
    """Real MainWindow, headless, using the platform-appropriate package."""
    from PyQt6.QtWidgets import QWidget, QMessageBox, QFileDialog

    monkeypatch.setattr(
        f"{_PKG}.ui.main_window_init.MainInitManager.init_worker_thread",
        lambda self: None,
        raising=False,
    )

    class _DummyPM:
        def __init__(self, main_window=None):
            self.plugins = []
            self.menu_actions = []
            self.toolbar_actions = []
            self.drop_handlers = []
            self.export_actions = []
            self.optimization_methods = {}
            self.file_openers = {}
            self.analysis_tools = []
            self.save_handlers = {}
            self.load_handlers = {}
            self.custom_3d_styles = {}
            self.document_reset_handlers = []
        def discover_plugins(self, parent=None): return []
        def update_plugin_menu(self, menu): pass
        def set_main_window(self, mw): pass
        def run_plugin(self, module, mw): pass
        def rebuild_plugin_menus(self): pass

    monkeypatch.setattr(
        f"{_PKG}.ui.main_window_init.PluginManager", _DummyPM, raising=False
    )

    class _DummyPlotter(QWidget):
        def __init__(self, parent=None, *a, **k):
            super().__init__(parent)
            self.renderer = _mock.MagicMock()
            self.add_mesh = _mock.MagicMock(return_value="actor")
            self.add_text = _mock.MagicMock()
            self.add_point_labels = _mock.MagicMock(return_value=["labels"])
            self.clear = _mock.MagicMock()
            self.reset_camera = _mock.MagicMock()
            self.render = _mock.MagicMock()
            self.set_background = _mock.MagicMock()
            self.setAcceptDrops = _mock.MagicMock()
            self.interactor = _mock.MagicMock()
            self.picker = _mock.MagicMock()
            self.camera_position = [(0, 0, 5), (0, 0, 0), (0, 1, 0)]
            self.camera = _mock.MagicMock()
            self.camera.copy = _mock.MagicMock(return_value={})
            self.add_light = _mock.MagicMock(return_value="light")
            self.remove_actor = _mock.MagicMock()

    monkeypatch.setattr(
        f"{_PKG}.ui.main_window_init.CustomQtInteractor", _DummyPlotter, raising=False
    )

    monkeypatch.setattr(
        f"{_PKG}.ui.main_window_init.MainInitManager.apply_initial_settings",
        lambda *a, **k: None, raising=False,
    )
    monkeypatch.setattr(
        f"{_PKG}.ui.view_3d_logic.View3DManager.apply_3d_settings",
        lambda *a, **k: None, raising=False,
    )

    # Real RDKit conversion, synchronous (no QThread)
    def _sync_trigger(self_compute):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from PyQt6.QtCore import QTimer

        mol = self_compute.host.state_manager.data.to_rdkit_mol(use_2d_stereo=False)
        if mol and mol.GetNumAtoms() > 0:
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv2()
            params.randomSeed = 42
            AllChem.EmbedMolecule(mol, params)
            AllChem.MMFFOptimizeMolecule(mol)
            AllChem.AssignAtomChiralTagsFromStructure(mol)
            QTimer.singleShot(0, lambda: self_compute.on_calculation_finished(mol))
        else:
            QTimer.singleShot(0, lambda: self_compute.on_calculation_finished(None))

    monkeypatch.setattr(
        f"{_PKG}.ui.compute_logic.ComputeManager.trigger_conversion",
        _sync_trigger, raising=False,
    )

    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    monkeypatch.setattr(
        QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes
    )
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))

    MainWindow = __import__(
        f"{_PKG}.ui.main_window", fromlist=["MainWindow"]
    ).MainWindow
    win = MainWindow()
    qtbot.addWidget(win)
    monkeypatch.setattr(win, "isVisible", lambda *a, **k: True, raising=False)

    try:
        yield win
    finally:
        try:
            win.closeEvent = lambda e: e.accept()
            for w in QApplication.topLevelWidgets():
                try:
                    w.closeEvent = lambda e: e.accept()
                    w.close()
                except Exception:
                    pass
            for _ in range(5):
                app.processEvents()
        except Exception:
            pass
