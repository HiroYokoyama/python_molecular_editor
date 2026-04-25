"""
Unit tests for ui/geometry_base_dialog.py (GeometryBaseDialog).

GeometryBaseDialog is abstract; a minimal concrete subclass is used throughout.

Covers:
  - _sync_input_to_slider: valid float syncs slider value
  - _sync_input_to_slider: invalid text is silently ignored
  - _sync_input_to_slider: wrap=True normalises values into (-180, 180]
  - on_slider_pressed: incomplete selection -> no-op
  - on_slider_pressed: complete selection -> sets _slider_dragging, snapshots positions
  - on_slider_released: clears _slider_dragging, calls draw_molecule_3d
  - on_slider_value_changed_click: while dragging -> early return
  - on_slider_value_changed_click: incomplete selection -> early return
  - on_slider_value_changed_click: normal -> updates input box, calls apply_geometry_update
  - on_slider_moved_realtime: updates input box, calls apply_geometry_update
"""

import os
import sys
import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QLineEdit, QSlider

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from moleditpy.ui.geometry_base_dialog import GeometryBaseDialog


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _ethane():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_mw():
    mw = MagicMock()
    mw._picking_consumed = False
    return mw


# ---------------------------------------------------------------------------
# Concrete subclass for testing the abstract base
# ---------------------------------------------------------------------------


class _ConcreteDialog(GeometryBaseDialog):
    """Minimal concrete subclass — no UI, selection flag controlled by tests."""

    def __init__(self, mol, mw):
        super().__init__(mol, mw)
        self._complete = False
        self.applied_values = []

    def _is_selection_complete(self):
        return self._complete

    def apply_geometry_update(self, new_value):
        self.applied_values.append(new_value)


@pytest.fixture
def dlg(qapp):
    mol = _ethane()
    mw = _make_mw()
    d = _ConcreteDialog(mol, mw)
    yield d
    try:
        d.close()
    except Exception:
        pass


def _slider():
    s = QSlider(Qt.Orientation.Horizontal)
    s.setMinimum(0)
    s.setMaximum(1000)
    s.setValue(0)
    return s


def _input():
    return QLineEdit()


# ---------------------------------------------------------------------------
# _sync_input_to_slider
# ---------------------------------------------------------------------------


class TestSyncInputToSlider:
    def test_valid_float_sets_slider(self, dlg):
        s = _slider()
        dlg._sync_input_to_slider("2.50", s, scale=100.0)
        assert s.value() == 250

    def test_invalid_text_is_ignored(self, dlg):
        s = _slider()
        s.setValue(123)
        dlg._sync_input_to_slider("not_a_number", s, scale=1.0)
        assert s.value() == 123  # unchanged

    def test_wrap_true_normalises_into_range(self, dlg):
        s = _slider()
        s.setMinimum(-180)
        s.setMaximum(180)
        # 270 -> (270+180)%360-180 = 450%360-180 = 90-180 = -90
        dlg._sync_input_to_slider("270", s, scale=1.0, wrap=True)
        assert s.value() == -90

    def test_wrap_false_uses_raw_value(self, dlg):
        s = _slider()
        s.setMinimum(-180)
        s.setMaximum(180)
        dlg._sync_input_to_slider("90", s, scale=1.0, wrap=False)
        assert s.value() == 90


# ---------------------------------------------------------------------------
# on_slider_pressed
# ---------------------------------------------------------------------------


class TestOnSliderPressed:
    def test_incomplete_selection_is_noop(self, dlg):
        dlg._complete = False
        dlg.on_slider_pressed()
        assert not dlg._slider_dragging
        assert dlg._snapshot_positions is None

    def test_complete_selection_sets_dragging_and_snapshot(self, dlg):
        dlg._complete = True
        dlg.on_slider_pressed()
        assert dlg._slider_dragging is True
        assert dlg._snapshot_positions is not None
        assert len(dlg._snapshot_positions) == dlg.mol.GetNumAtoms()


# ---------------------------------------------------------------------------
# on_slider_released
# ---------------------------------------------------------------------------


class TestOnSliderReleased:
    def test_clears_dragging_flag(self, dlg):
        dlg._slider_dragging = True
        dlg.on_slider_released()
        assert dlg._slider_dragging is False

    def test_calls_draw_molecule_3d(self, dlg):
        dlg.on_slider_released()
        dlg.main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with(
            dlg.mol
        )


# ---------------------------------------------------------------------------
# on_slider_value_changed_click
# ---------------------------------------------------------------------------


class TestOnSliderValueChangedClick:
    def test_skips_when_dragging(self, dlg):
        dlg._slider_dragging = True
        dlg._complete = True
        inp = _input()
        dlg.on_slider_value_changed_click(200, inp, scale=100.0)
        assert dlg.applied_values == []

    def test_skips_when_selection_incomplete(self, dlg):
        dlg._slider_dragging = False
        dlg._complete = False
        inp = _input()
        dlg.on_slider_value_changed_click(200, inp, scale=100.0)
        assert dlg.applied_values == []

    def test_updates_input_and_calls_apply(self, dlg):
        dlg._slider_dragging = False
        dlg._complete = True
        dlg._snapshot_positions = dlg.mol.GetConformer().GetPositions().copy()
        inp = _input()
        dlg.on_slider_value_changed_click(154, inp, scale=100.0)
        assert inp.text() == "1.540"
        assert dlg.applied_values == [pytest.approx(1.54)]


# ---------------------------------------------------------------------------
# on_slider_moved_realtime
# ---------------------------------------------------------------------------


class TestOnSliderMovedRealtime:
    def test_skips_when_incomplete(self, dlg):
        dlg._complete = False
        inp = _input()
        dlg.on_slider_moved_realtime(100, inp, scale=100.0)
        assert dlg.applied_values == []

    def test_updates_input_and_calls_apply(self, dlg):
        dlg._complete = True
        inp = _input()
        dlg.on_slider_moved_realtime(200, inp, scale=100.0)
        assert inp.text() == "2.000"
        assert dlg.applied_values == [pytest.approx(2.0)]
