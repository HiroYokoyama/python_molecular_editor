from types import SimpleNamespace

import numpy as np

from moleditpy.ui.atom_picking import pick_atom_index_from_screen


class _Renderer:
    def SetWorldPoint(self, x, y, z, w):
        self._world_point = (x, y, z, w)

    def WorldToDisplay(self):
        pass

    def GetDisplayPoint(self):
        x, y, z, _ = self._world_point
        return (100.0 + x * 20.0, 100.0 + y * 20.0, z)


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
