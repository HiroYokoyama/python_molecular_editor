from types import SimpleNamespace
import numpy as np
from moleditpy.ui.atom_picking import (
    pick_atom_index_from_screen,
    pick_atom_index_from_screen_vectorized,
    pick_atom_index_from_screen_sequential,
)


class _VtkMatrix:
    def __init__(self, elements):
        self.elements = elements

    def GetElement(self, i, j):
        return self.elements[i][j]


class _Camera:
    def GetCompositeProjectionTransformMatrix(self, aspect_ratio, near, far):
        # We need a matrix that maps x,y directly to display coords via:
        # display_x = (ndc_x + 1) * 0.5 * size_x
        # For size = (200, 200), display_x = 100 * ndc_x + 100
        # To match sequential display: 100.0 + x * 20.0, we want:
        # 100 * ndc_x + 100 = 100 + x * 20 => ndc_x = x * 0.2
        # ndc_x = clip_x / w. For w = 1.0, clip_x = 0.2 * x.
        # So matrix row 0: [0.2, 0, 0, 0]
        # matrix row 1: [0, 0.2, 0, 0]
        # matrix row 3: [0, 0, 0, 1]
        return _VtkMatrix(
            [
                [0.2, 0.0, 0.0, 0.0],
                [0.0, 0.2, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

    def GetParallelProjection(self):
        return True

    def GetParallelScale(self):
        return 1.0

    def GetViewUp(self):
        return (0.0, 1.0, 0.0)


class _Renderer:
    def SetWorldPoint(self, x, y, z, w):
        self._world_point = (x, y, z, w)

    def WorldToDisplay(self):
        pass

    def GetDisplayPoint(self):
        x, y, z, _ = self._world_point
        return (100.0 + x * 20.0, 100.0 + y * 20.0, z)

    def GetActiveCamera(self):
        return _Camera()

    def GetTiledAspectRatio(self):
        return 1.0

    def GetSize(self):
        return (200, 200)


class _Atom:
    def GetSymbol(self):
        return "C"


class _Mol:
    def GetNumAtoms(self):
        return 2

    def GetAtomWithIdx(self, _idx):
        return _Atom()


def _view():
    return SimpleNamespace(
        atom_positions_3d=np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]),
        current_3d_style="ball_and_stick",
        current_mol=_Mol(),
        host=SimpleNamespace(
            init_manager=SimpleNamespace(settings={"ball_stick_atom_scale": 1.0})
        ),
        plotter=SimpleNamespace(renderer=_Renderer()),
    )


def test_pick_atom_index_from_screen_hits_projected_atom_edge():
    assert pick_atom_index_from_screen(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_returns_none_for_background():
    assert pick_atom_index_from_screen(_view(), (200, 200), _Mol()) is None


def test_pick_atom_index_from_screen_vectorized_success():
    assert pick_atom_index_from_screen_vectorized(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_sequential_success():
    assert pick_atom_index_from_screen_sequential(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_fallback():
    view_obj = _view()
    # Deliberately mock GetActiveCamera as None to force fallback to sequential
    view_obj.plotter.renderer.GetActiveCamera = None
    assert pick_atom_index_from_screen(view_obj, (111, 100), _Mol()) == 0
