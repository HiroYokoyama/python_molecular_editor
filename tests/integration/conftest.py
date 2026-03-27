import pytest
import os
import sys
import types
from unittest import mock as _mock

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


class DummyAxesWidget:
    def SetOrientationMarker(self, marker):
        pass

    def SetInteractor(self, interactor):
        pass

    def SetViewport(self, *a):
        pass

    def On(self):
        pass

    def Off(self):
        pass

    def SetInteractive(self, v):
        pass


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


vinter.vtkInteractorStyleTrackballCamera = DummyInteractorStyleBase

# Mock PyVista and pyvistaqt
pyv = force_mock_module("pyvista")


class DummyPlotterMinimal:
    def __init__(self, *a, **k):
        pass

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


class DummyPolyDataMock(_mock.MagicMock):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.point_data = {}
        self.cell_data = {}

    def __getitem__(self, key):
        return _mock.MagicMock()

    def __setitem__(self, key, val):
        pass

    def glyph(self, *a, **k):
        return _mock.MagicMock()

    def tube(self, *a, **k):
        return _mock.MagicMock()


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
if "PyQt6.QtWidgets" in sys.modules:
    try:
        from PyQt6.QtWidgets import QWidget

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

from PyQt6.QtWidgets import QApplication


# 3. Minimal cleanup to avoid COM crashes
@pytest.fixture(scope="session")
def app():
    """QApplication session-wide instance with platform-aware teardown."""
    is_headless = os.environ.get("MOLEDITPY_HEADLESS", "0") == "1"
    is_offscreen = os.environ.get("QT_QPA_PLATFORM") == "offscreen"

    if is_headless:
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)

    yield q_app

    # Platform-aware teardown: prevent 0x80010108/0x8001010d on Windows but avoid CI segfaults
    if not (is_headless or is_offscreen):
        try:
            q_app.closeAllWindows()
            for _ in range(10):
                q_app.processEvents()

            # Force garbage collection to release Python-wrapped Qt objects
            import gc

            gc.collect()
            q_app.processEvents()

            try:
                import colorama

                colorama.deinit()
            except (ImportError, Exception):
                pass
                
            # Wait for any lingering background threads to finish.
            import threading
            import time

            main_thread = threading.main_thread()
            deadline = time.monotonic() + 2.0  # 2 second max wait
            for t in threading.enumerate():
                if t is main_thread or not t.is_alive():
                    continue
                remaining = deadline - time.monotonic()
                if remaining > 0:
                    t.join(timeout=remaining)

            # Final event processing after threads are done
            for _ in range(5):
                q_app.processEvents()
        except Exception:
            pass

    # Suppress Windows fatal exception handler for known teardown crashes.
    if sys.platform == "win32":
        try:
            import faulthandler
            faulthandler.disable()
        except Exception:
            pass

    q_app.quit()


@pytest.fixture
def window(app, qtbot, monkeypatch):
    """
    Create a MainWindow instance with necessary mocks for integration tests.
    """
    from PyQt6.QtWidgets import QWidget, QMessageBox

    # Pre-patch init_worker_thread to avoid QThread creation
    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.MainWindowMainInit.init_worker_thread",
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
            self.add_light = _mock.MagicMock(return_value="light_actor")
            self.remove_actor = _mock.MagicMock()

    monkeypatch.setattr(
        "moleditpy.ui.main_window_init.CustomQtInteractor",
        DummyPlotter,
        raising=False,
    )

    # 3. Disable complex initializations to avoid crashes
    monkeypatch.setattr(
        "moleditpy.ui.main_window.MainWindow.apply_initial_settings",
        lambda *a, **k: None,
    )
    monkeypatch.setattr(
        "moleditpy.ui.view_3d.MainWindowView3d.apply_3d_settings",
        lambda *a, **k: None,
    )

    # 4. Synchronous trigger_conversion
    def sync_trigger_conversion(self_compute):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from PyQt6.QtCore import QTimer

        mol = self_compute.data.to_rdkit_mol(use_2d_stereo=False)
        if mol and mol.GetNumAtoms() > 0:
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
        "moleditpy.core.compute_engine.MainWindowCompute.trigger_conversion",
        sync_trigger_conversion,
    )

    # 5. Mock QMessageBox and QFileDialog to avoid blocking tests
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    monkeypatch.setattr(
        QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes
    )

    from PyQt6.QtWidgets import QFileDialog

    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: "")

    # 6. Instantiate
    win = MainWindow()
    qtbot.addWidget(win)
    # Don't win.show() if headless to be safer, or win.show() if needed by qtbot
    if os.environ.get("MOLEDITPY_HEADLESS", "0") != "1":
        win.show()

    monkeypatch.setattr(win, "isVisible", lambda *a, **k: True, raising=False)

    try:
        yield win
    finally:
        # Clean teardown: stop threads, close window and process events to avoid COM/RPC stale objects
        try:
            # Forcefully prevent "unsaved changes" dialogs
            try:
                win.has_unsaved_changes = False
            except Exception:
                import traceback

                traceback.print_exc()

            # Override closeEvent to avoid any blocking logic
            try:
                win.closeEvent = lambda event: event.accept()
            except Exception:
                import traceback

                traceback.print_exc()

            # Stop any active threads explicitly before closing
            active_threads = list(getattr(win, "_active_calc_threads", []) or [])
            for thr in active_threads:
                try:
                    if hasattr(thr, "isRunning") and thr.isRunning():
                        thr.quit()
                        thr.wait(100)
                except Exception:
                    import traceback

                    traceback.print_exc()

            # Aggressive auto-close: ensure all top-level widgets are closed
            # to satisfy "auto close window" request and prevent COM hangs.
            win.hide()
            try:
                win.close()
            except Exception:
                import traceback

                traceback.print_exc()
            try:
                win.deleteLater()
            except Exception:
                import traceback

                traceback.print_exc()

            from PyQt6.QtWidgets import QApplication

            QApplication.closeAllWindows()

            # Thoroughly process events to clear any pending COM calls or slots
            # Use a small number of iterations with short sleep if needed,
            # but here we just processEvents multiple times.
            for _ in range(10):
                app.processEvents()

        except Exception:
            import traceback

            traceback.print_exc()

        # Attempt to de-initialize colorama if it's wrapping stdout/stderr
        # which can cause 0x80010012/0x80010108 fatal exceptions on Windows teardown
        try:
            import colorama

            colorama.deinit()
        except ImportError:
            pass
        except Exception:
            import traceback

            traceback.print_exc()

        # Force garbage collection to clean up potential circular references or lingering Qt objects
        import gc

        gc.collect()
