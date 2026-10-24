import pytest
import os
import sys
import types
from unittest import mock as _mock

from rendering_stubs import DummyAxesWidget, DummyPolyDataMock

# 1. Path setup
src_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
src_moleditpy = os.path.join(src_root, "moleditpy")

if os.path.isdir(src_root) and src_root not in sys.path:
    sys.path.insert(0, src_root)
if os.path.isdir(src_moleditpy) and src_moleditpy not in sys.path:
    sys.path.insert(0, src_moleditpy)


# 2. Aggressive Mocking at top-level to prevent any COM interaction
def force_mock_module(name):
    mod = types.ModuleType(name)
    sys.modules[name] = mod
    return mod


# Mock VTK
vtk = force_mock_module("vtk")


class DummyCellPicker:
    def __init__(self, *a, **k):
        self._actor = None

    def SetTolerance(self, v):
        self._tol = v

    def Pick(self, x, y, z, renderer):
        return

    def GetActor(self):
        return self._actor


vtk.vtkCellPicker = DummyCellPicker
vtk.vtkAxesActor = _mock.MagicMock


vtk.vtkOrientationMarkerWidget = DummyAxesWidget
vtk.vtkInteractorStyleTrackballCamera = _mock.MagicMock

# Mock vtkmodules
vmod = force_mock_module("vtkmodules")
vmod.__path__ = []
force_mock_module("vtkmodules.all")
force_mock_module("vtkmodules.vtkRenderingCore")
force_mock_module("vtkmodules.vtkCommonCore")
vinter = force_mock_module("vtkmodules.vtkInteractionStyle")


class DummyInteractorStyleBase:
    def __init__(self, *a, **k):
        pass

    def AddObserver(self, *_a, **_k):
        return

    def GetInteractor(self):
        return _mock.MagicMock()

    def OnLeftButtonDown(self, *a, **k):
        pass

    def OnLeftButtonUp(self, *a, **k):
        pass

    def OnRightButtonDown(self, *a, **k):
        pass

    def OnRightButtonUp(self, *a, **k):
        pass

    def OnMouseMove(self, *a, **k):
        pass


vinter.vtkInteractorStyleTrackballCamera = DummyInteractorStyleBase

# Mock PyVista and pyvistaqt
pyv = force_mock_module("pyvista")


class DummyPlotterMinimal:
    def __init__(self, *a, **k):
        self.camera_position = [(0, 0, 5), (0, 0, 0), (0, 1, 0)]

    def add_mesh(self, *a, **k):
        return _mock.MagicMock()

    def add_point_labels(self, *a, **k):
        return []

    def clear(self):
        return

    def render(self):
        return

    def reset_camera(self):
        return

    def add_light(self, *a, **k):
        return _mock.MagicMock()

    def remove_actor(self, *a, **k):
        return True

    @property
    def camera(self):
        c = _mock.MagicMock()
        c.copy = _mock.MagicMock(return_value={})
        return c


pyv.Plotter = DummyPlotterMinimal


pyv.PolyData = DummyPolyDataMock
pyv.Light = _mock.MagicMock
pyv.Sphere = _mock.MagicMock
pyv.Spline = _mock.MagicMock
# Add other common pyvista shapes used in the app
pyv.Cylinder = _mock.MagicMock
pyv.Box = _mock.MagicMock

pvqt = force_mock_module("pyvistaqt")
pvqt.BackgroundPlotter = DummyPlotterMinimal
# Handle QWidget inheritance safely for mocking
_parent_class = object
if "PySide6.QtWidgets" in sys.modules:
    try:
        from PySide6.QtWidgets import QWidget

        _parent_class = QWidget
    except ImportError:
        pass


class DummyQtInteractorMock(_parent_class):
    def __init__(self, parent=None, *a, **k):
        if _parent_class is not object:
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
        self.add_light = _mock.MagicMock(return_value=_mock.MagicMock())
        self.remove_actor = _mock.MagicMock()


pvqt.QtInteractor = DummyQtInteractorMock

from PySide6.QtCore import QCoreApplication, QEvent


@pytest.fixture
def window(app, qtbot, monkeypatch, tmp_path):
    """
    Create a MainWindow instance with necessary mocks for integration tests.
    """
    from PySide6.QtWidgets import QWidget, QMessageBox

    # Pre-patch init_worker_thread to avoid QThread creation
    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.MainInitManager.init_worker_thread",
        lambda self: None,
        raising=False,
    )

    from moleditpy.ui.main_window import MainWindow

    # 1. Mock PluginManager
    class DummyPluginManager:
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

        def discover_plugins(self, parent=None):
            return []

        def invoke_document_reset_handlers(self):
            pass

        def update_plugin_menu(self, menu):
            pass

        def set_main_window(self, mw):
            pass

        def run_plugin(self, module, main_window):
            pass

    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.PluginManager",
        DummyPluginManager,
        raising=False,
    )

    # 2. Dummy Plotter
    class DummyPlotter(QWidget):
        def __init__(self, parent=None, *a, **k):
            super().__init__(parent)
            self.renderer = _mock.MagicMock()
            self.add_mesh = _mock.MagicMock(return_value="dummy")
            self.add_text = _mock.MagicMock()
            self.add_point_labels = _mock.MagicMock(return_value=["point_labels"])
            self.clear = _mock.MagicMock()
            self.reset_camera = _mock.MagicMock()
            self.render = _mock.MagicMock()
            self.set_background = _mock.MagicMock()
            self.setAcceptDrops = _mock.MagicMock()
            self.interactor = _mock.MagicMock()
            self.picker = _mock.MagicMock()
            self.camera = _mock.MagicMock()
            self.camera.copy = _mock.MagicMock(return_value={})
            self.camera_position = [(0, 0, 5), (0, 0, 0), (0, 1, 0)]
            self.add_light = _mock.MagicMock(return_value="light_actor")
            self.remove_actor = _mock.MagicMock()

    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.CustomQtInteractor",
        DummyPlotter,
        raising=False,
    )

    # 3. Disable complex initializations to avoid crashes
    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.MainInitManager.apply_initial_settings",
        lambda *a, **k: None,
        raising=False,
    )
    monkeypatch.setattr(
        "moleditpy.ui.view_3d_logic.View3DManager.apply_3d_settings",
        lambda *a, **k: None,
    )

    # 4. Synchronous trigger_conversion
    def sync_trigger_conversion(self_compute):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from PySide6.QtCore import QTimer

        mol = self_compute.host.state_manager.data.to_rdkit_mol(use_2d_stereo=False)
        if mol and mol.GetNumAtoms() > 0:
            # Honor the drawn wedges as the real worker does (it parses them
            # from the MOL block). Without this the embedding picks an
            # enantiomer at random and the post-conversion chirality check
            # rightly reports the result as wrong. Before AddHs, which puts the
            # new atoms at the origin.
            Chem.AssignChiralTypesFromBondDirs(mol)
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDG()
            params.randomSeed = 42
            AllChem.EmbedMolecule(mol, params)
            AllChem.MMFFOptimizeMolecule(mol)
            AllChem.AssignAtomChiralTagsFromStructure(mol)
            # Use QTimer to avoid COM/re-entrancy issues
            QTimer.singleShot(0, lambda: self_compute.on_calculation_finished(mol))
        else:
            QTimer.singleShot(0, lambda: self_compute.on_calculation_finished(None))

    monkeypatch.setattr(
        "moleditpy.ui.compute_logic.ComputeManager.trigger_conversion",
        sync_trigger_conversion,
    )

    # 5. Mock QMessageBox and QFileDialog to avoid blocking tests
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    monkeypatch.setattr(
        QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes
    )

    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: "")

    # 6. Instantiate
    win = MainWindow()
    # A test that touches settings must never rewrite the real ~/.moleditpy
    win.init_manager.settings_dir = str(tmp_path / ".moleditpy")
    win.init_manager.settings_file = str(tmp_path / ".moleditpy" / "settings.json")
    qtbot.addWidget(win)
    # Don't win.show() if headless to be safer, or win.show() if needed by qtbot
    if os.environ.get("MOLEDITPY_HEADLESS", "0") != "1":
        win.show()

    monkeypatch.setattr(win, "isVisible", lambda *a, **k: True, raising=False)

    yield win

    # Only this fixture's window belongs here. Deleting unrelated top-level
    # dialogs can race their timers and invalidate another fixture's wrappers.
    win.removeEventFilter(win.ui_manager)
    win.closeEvent = lambda event: event.accept()
    win.close()
    win.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
