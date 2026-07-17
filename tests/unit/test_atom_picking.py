"""Unit tests for 3D atom picking from screen coordinates."""

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
    """Click on a projected atom edge returns that atom's index."""
    assert pick_atom_index_from_screen(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_returns_none_for_background():
    """Click on empty background returns None."""
    assert pick_atom_index_from_screen(_view(), (200, 200), _Mol()) is None


def test_pick_atom_index_from_screen_vectorized_success():
    """Vectorized implementation returns the correct atom index on a hit."""
    assert pick_atom_index_from_screen_vectorized(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_sequential_success():
    """Sequential implementation returns the correct atom index on a hit."""
    assert pick_atom_index_from_screen_sequential(_view(), (111, 100), _Mol()) == 0


def test_pick_atom_index_from_screen_fallback():
    """Falls back to sequential picking when the camera is unavailable."""
    view_obj = _view()
    # Deliberately mock GetActiveCamera as None to force fallback to sequential
    view_obj.plotter.renderer.GetActiveCamera = None
    assert pick_atom_index_from_screen(view_obj, (111, 100), _Mol()) == 0


# ---------------------------------------------------------------------------
# Helper internals and guard paths
# ---------------------------------------------------------------------------

from unittest.mock import patch

from moleditpy.ui.atom_picking import (
    _atom_world_radius,
    _projected_radius_px,
    _world_to_display,
)


class _BrokenRenderer:
    def SetWorldPoint(self, *a):
        raise RuntimeError("renderer gone")


def test_world_to_display_returns_none_on_broken_renderer():
    assert _world_to_display(_BrokenRenderer(), (0.0, 0.0, 0.0)) is None


def test_projected_radius_none_when_center_unprojectable():
    assert _projected_radius_px(_BrokenRenderer(), (0, 0, 0), 1.0) is None


class _BrokenMol:
    def GetAtomWithIdx(self, _idx):
        raise RuntimeError("no atoms")


def test_atom_world_radius_defaults_to_carbon_on_broken_mol():
    view_obj = _view()
    radius = _atom_world_radius(view_obj, _BrokenMol(), 0)
    assert radius == _atom_world_radius(view_obj, _Mol(), 0)


def test_atom_world_radius_survives_missing_settings():
    view_obj = _view()
    view_obj.host = None  # settings access raises AttributeError
    radius = _atom_world_radius(view_obj, _Mol(), 0)
    assert radius > 0


def test_atom_world_radius_stick_uses_bond_radius_setting():
    view_obj = _view()
    view_obj.current_3d_style = "stick"
    view_obj.host.init_manager.settings["stick_bond_radius"] = 0.25
    assert _atom_world_radius(view_obj, _Mol(), 0) == 0.25


def test_atom_world_radius_wireframe_is_hairline():
    view_obj = _view()
    view_obj.current_3d_style = "Wireframe"
    assert _atom_world_radius(view_obj, _Mol(), 0) == 0.01


def test_atom_world_radius_cpk_scales_rvdw():
    view_obj = _view()
    view_obj.current_3d_style = "CPK"
    view_obj.host.init_manager.settings["cpk_atom_scale"] = 2.0
    radius = _atom_world_radius(view_obj, _Mol(), 0)
    assert radius > 2.0  # carbon Rvdw (~1.7) doubled


def test_atom_world_radius_style_name_normalized():
    view_obj = _view()
    view_obj.current_3d_style = "Ball and Stick"
    assert _atom_world_radius(view_obj, _Mol(), 0) == _atom_world_radius(
        _view(), _Mol(), 0
    )


# ---------------------------------------------------------------------------
# Sequential guards
# ---------------------------------------------------------------------------


def test_sequential_none_when_no_plotter():
    view_obj = SimpleNamespace(atom_positions_3d=None)  # no plotter attr
    assert pick_atom_index_from_screen_sequential(view_obj, (0, 0)) is None


def test_sequential_none_when_positions_missing():
    view_obj = _view()
    view_obj.atom_positions_3d = None
    assert pick_atom_index_from_screen_sequential(view_obj, (0, 0)) is None


def test_sequential_none_when_positions_not_numeric():
    view_obj = _view()
    view_obj.atom_positions_3d = [["a", "b", "c"]]
    assert pick_atom_index_from_screen_sequential(view_obj, (0, 0)) is None


def test_sequential_none_when_positions_wrong_shape():
    view_obj = _view()
    view_obj.atom_positions_3d = np.array([1.0, 2.0, 3.0])  # 1-D
    assert pick_atom_index_from_screen_sequential(view_obj, (0, 0)) is None


def test_sequential_uses_current_mol_when_mol_omitted():
    view_obj = _view()
    assert pick_atom_index_from_screen_sequential(view_obj, (111, 100)) == 0


def test_sequential_falls_back_to_position_count_on_broken_mol():
    class _NoCountMol(_Mol):
        def GetNumAtoms(self):
            raise RuntimeError("boom")

    view_obj = _view()
    assert (
        pick_atom_index_from_screen_sequential(view_obj, (111, 100), _NoCountMol())
        == 0
    )


# ---------------------------------------------------------------------------
# Vectorized guards and perspective path
# ---------------------------------------------------------------------------


def test_vectorized_none_when_positions_empty():
    view_obj = _view()
    view_obj.atom_positions_3d = np.zeros((0, 3))
    assert pick_atom_index_from_screen_vectorized(view_obj, (0, 0), _Mol()) is None


def test_vectorized_none_when_size_unavailable():
    class _NoSizeRenderer(_Renderer):
        def GetSize(self):
            raise RuntimeError("no window")

    view_obj = _view()
    view_obj.plotter = SimpleNamespace(renderer=_NoSizeRenderer())
    assert pick_atom_index_from_screen_vectorized(view_obj, (111, 100), _Mol()) is None


def test_vectorized_perspective_projection_path():
    class _PerspectiveCamera(_Camera):
        def GetParallelProjection(self):
            return False

        def GetViewAngle(self):
            return 30.0

    class _PerspectiveRenderer(_Renderer):
        def GetActiveCamera(self):
            return _PerspectiveCamera()

    view_obj = _view()
    view_obj.plotter = SimpleNamespace(renderer=_PerspectiveRenderer())
    assert pick_atom_index_from_screen_vectorized(view_obj, (100, 100), _Mol()) == 0


def test_vectorized_pixel_scale_fallback_on_camera_error():
    class _FlakyCamera(_Camera):
        def GetParallelProjection(self):
            raise RuntimeError("camera busy")

    class _FlakyRenderer(_Renderer):
        def GetActiveCamera(self):
            return _FlakyCamera()

    view_obj = _view()
    view_obj.plotter = SimpleNamespace(renderer=_FlakyRenderer())
    # Fallback pixel scale still yields a hit near the projected atom
    assert pick_atom_index_from_screen_vectorized(view_obj, (100, 100), _Mol()) == 0


def test_pick_falls_back_when_vectorized_raises():
    view_obj = _view()
    with patch(
        "moleditpy.ui.atom_picking.pick_atom_index_from_screen_vectorized",
        side_effect=RuntimeError("vtk exploded"),
    ):
        assert pick_atom_index_from_screen(view_obj, (111, 100), _Mol()) == 0
