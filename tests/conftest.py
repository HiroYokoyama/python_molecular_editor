import sys
import os
import types
from unittest import mock

# 1. Path setup
src_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "moleditpy", "src")
)
src_moleditpy = os.path.join(src_root, "moleditpy")

if os.path.isdir(src_root) and src_root not in sys.path:
    sys.path.insert(0, src_root)
if os.path.isdir(src_moleditpy) and src_moleditpy not in sys.path:
    sys.path.insert(0, src_moleditpy)


# 2. Aggressive Mocking at top-level to prevent ANY library from loading real C++ extensions
modules_to_patch = [
    "vtk",
    "vtkmodules",
    "vtkmodules.all",
    "vtkmodules.vtkRenderingOpenGL2",
    "vtkmodules.vtkRenderingCore",
    "vtkmodules.vtkCommonCore",
    "vtkmodules.vtkInteractionStyle",
    "vtkmodules.vtkWebCore",
    "vtkmodules.vtkCommonMath",
    "vtkmodules.qt.QVTKRenderWindowInteractor",
    "pyvista",
    "pyvistaqt",
]

for name in modules_to_patch:
    m = types.ModuleType(name)
    m.__path__ = []
    sys.modules[name] = m
    if "." in name:
        parent_name, child_name = name.rsplit(".", 1)
        if parent_name in sys.modules:
            setattr(sys.modules[parent_name], child_name, m)


# 3. Provide dummy classes within the mocked submodules to satisfy 'from ... import ...'
class DummyCellPicker:
    def __init__(self, *a, **k):
        self._actor = None

    def SetTolerance(self, v):
        pass

    def Pick(self, x, y, z, renderer):
        pass

    def GetActor(self):
        return None


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


class DummyInteractorStyle:
    def __init__(self, *a, **k):
        pass

    def AddObserver(self, *a, **k):
        pass

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


# Attach to vtk
vtk = sys.modules["vtk"]
vtk.vtkCellPicker = DummyCellPicker
vtk.vtkAxesActor = mock.MagicMock
vtk.vtkOrientationMarkerWidget = DummyAxesWidget
vtk.vtkInteractorStyleTrackballCamera = DummyInteractorStyle

# Attach to specific vtkmodules submodules
if "vtkmodules.vtkRenderingCore" in sys.modules:
    sys.modules["vtkmodules.vtkRenderingCore"].vtkCellPicker = DummyCellPicker
if "vtkmodules.vtkInteractionStyle" in sys.modules:
    sys.modules[
        "vtkmodules.vtkInteractionStyle"
    ].vtkInteractorStyleTrackballCamera = DummyInteractorStyle

# PyVista mocks
pyv = sys.modules["pyvista"]


class DummyPlotterMinimal:
    def __init__(self, *a, **k):
        pass

    def add_mesh(self, *a, **k):
        return mock.MagicMock()

    def add_point_labels(self, *a, **k):
        return []

    def clear(self):
        pass

    def render(self):
        pass

    def reset_camera(self):
        pass

    def add_light(self, *a, **k):
        return mock.MagicMock()

    @property
    def camera(self):
        c = mock.MagicMock()
        c.copy = mock.MagicMock(return_value={})
        return c


pyv.Plotter = DummyPlotterMinimal


class DummyPolyDataMock(mock.MagicMock):
    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.point_data = {}
        self.cell_data = {}

    def __getitem__(self, key):
        return mock.MagicMock()

    def __setitem__(self, key, val):
        pass

    def glyph(self, *a, **k):
        return mock.MagicMock()

    def tube(self, *a, **k):
        return mock.MagicMock()


pyv.PolyData = DummyPolyDataMock

# pyvistaqt mocks
pvqt = sys.modules["pyvistaqt"]
pvqt.BackgroundPlotter = DummyPlotterMinimal


class DummyQtInteractorMock(mock.MagicMock):
    def __init__(self, parent=None, *a, **k):
        super().__init__(*a, **k)
        self.renderer = mock.MagicMock()
        self.add_mesh = mock.MagicMock(return_value=mock.MagicMock())
        self.add_point_labels = mock.MagicMock(return_value=[])
        self.clear = mock.MagicMock()
        self.render = mock.MagicMock()
        self.reset_camera = mock.MagicMock()


pvqt.QtInteractor = DummyQtInteractorMock


def pytest_configure(config):
    """Placeholder hook to ensure pytest sees this as a valid conftest."""
    pass
