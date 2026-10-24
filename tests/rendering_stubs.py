"""Rendering stand-ins shared by the headless test fixtures.

Keep tier-specific picker, plotter, and interactor behavior in each conftest.
This module must stay independent of Qt, VTK, and application imports.
"""

from typing import Any
from unittest import mock


class DummyAxesWidget:
    def SetOrientationMarker(self, marker: object) -> None:
        pass

    def SetInteractor(self, interactor: object) -> None:
        pass

    def SetViewport(self, *a: Any) -> None:
        pass

    def On(self) -> None:
        pass

    def Off(self) -> None:
        pass

    def SetInteractive(self, v: object) -> None:
        pass


class DummyPolyDataMock(mock.MagicMock):
    def __init__(self, *a: Any, **k: Any) -> None:
        super().__init__(*a, **k)
        self.point_data = {}
        self.cell_data = {}

    def __getitem__(self, key: object) -> mock.MagicMock:
        return mock.MagicMock()

    def __setitem__(self, key: object, val: object) -> None:
        pass

    def glyph(self, *a: Any, **k: Any) -> mock.MagicMock:
        return mock.MagicMock()

    def tube(self, *a: Any, **k: Any) -> mock.MagicMock:
        return mock.MagicMock()
