"""Unit tests for the drag-to-draw alkyl chain tool (ChainMixin)."""

import math

import pytest
from unittest.mock import MagicMock
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QApplication

from moleditpy.ui.chain_mixin import (
    CHAIN_ANGLE_SNAP_DEG,
    CHAIN_HALF_ANGLE_DEG,
    MAX_CHAIN_LENGTH,
    ChainMixin,
)
from moleditpy.utils.constants import DEFAULT_BOND_LENGTH

AXIS_STEP = DEFAULT_BOND_LENGTH * math.cos(math.radians(CHAIN_HALF_ANGLE_DEG))


class MockChainScene(ChainMixin):
    def __init__(self):
        self.template_preview = MagicMock()
        self.add_molecule_fragment = MagicMock()
        self.current_atom_symbol = "C"
        super().__init__()


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def scene(qapp):
    return MockChainScene()


def bond_lengths(points):
    return [
        math.hypot(b.x() - a.x(), b.y() - a.y()) for a, b in zip(points, points[1:])
    ]


def test_every_bond_keeps_the_standard_length(scene):
    points = scene.calculate_chain_points(QPointF(0, 0), QPointF(5 * AXIS_STEP, 0))
    assert len(points) == 6
    for length in bond_lengths(points):
        assert length == pytest.approx(DEFAULT_BOND_LENGTH)


def test_bond_count_follows_drag_distance(scene):
    for expected in (1, 2, 7):
        points = scene.calculate_chain_points(
            QPointF(0, 0), QPointF(expected * AXIS_STEP, 0)
        )
        assert len(points) - 1 == expected


def test_chain_zigzags_around_the_axis(scene):
    points = scene.calculate_chain_points(QPointF(0, 0), QPointF(4 * AXIS_STEP, 0))
    offsets = [p.y() for p in points]
    assert offsets[0] == pytest.approx(0.0)
    assert offsets[2] == pytest.approx(0.0)
    assert abs(offsets[1]) == pytest.approx(
        DEFAULT_BOND_LENGTH * math.sin(math.radians(CHAIN_HALF_ANGLE_DEG))
    )
    # Vertices advance along the axis by a constant step.
    for i, p in enumerate(points):
        assert p.x() == pytest.approx(i * AXIS_STEP)


def test_zigzag_side_follows_the_cursor(scene):
    below = scene.calculate_chain_points(QPointF(0, 0), QPointF(2 * AXIS_STEP, 3))
    above = scene.calculate_chain_points(QPointF(0, 0), QPointF(2 * AXIS_STEP, -3))
    assert below[1].y() > 0
    assert above[1].y() < 0


def test_axis_is_snapped_to_fixed_angle_steps(scene):
    # 20 deg drag snaps down to the nearest 15 deg step.
    angle = math.radians(20)
    far = QPointF(
        3 * AXIS_STEP * math.cos(angle),
        3 * AXIS_STEP * math.sin(angle),
    )
    points = scene.calculate_chain_points(QPointF(0, 0), far)
    # Even-indexed vertices sit exactly on the snapped axis.
    axis_angle = math.degrees(math.atan2(points[2].y(), points[2].x()))
    assert axis_angle == pytest.approx(CHAIN_ANGLE_SNAP_DEG)


def test_short_or_backward_drag_yields_no_chain(scene):
    assert scene.calculate_chain_points(QPointF(0, 0), QPointF(0, 0)) == []
    assert scene.calculate_chain_points(QPointF(0, 0), QPointF(4, 0)) == []


def test_chain_length_is_capped(scene):
    points = scene.calculate_chain_points(QPointF(0, 0), QPointF(5000 * AXIS_STEP, 0))
    assert len(points) - 1 == MAX_CHAIN_LENGTH


def test_preview_shows_repeat_count_near_the_cursor(scene):
    scene.begin_chain(QPointF(0, 0))
    cursor = QPointF(3 * AXIS_STEP, 0)
    scene.update_chain_preview(cursor)

    points, label, label_pos = scene.template_preview.set_chain_geometry.call_args[0]
    assert label == "n = 3"
    assert len(points) == 4
    assert label_pos.x() > cursor.x()
    scene.template_preview.show.assert_called_once()


def test_preview_hides_when_the_drag_is_too_short(scene):
    scene.begin_chain(QPointF(0, 0))
    scene.update_chain_preview(QPointF(2, 0))

    scene.template_preview.set_chain_geometry.assert_not_called()
    scene.template_preview.hide.assert_called_once()


def test_preview_ignores_moves_outside_a_drag(scene):
    scene.update_chain_preview(QPointF(3 * AXIS_STEP, 0))
    scene.template_preview.set_chain_geometry.assert_not_called()


def test_begin_chain_anchors_on_the_start_atom(scene):
    atom = MagicMock()
    atom.pos.return_value = QPointF(10, 20)
    scene.begin_chain(QPointF(13, 24), atom)

    assert scene.chain_active is True
    assert scene.chain_anchor == QPointF(10, 20)
    assert scene.chain_start_atom is atom


def test_commit_builds_a_single_bonded_path(scene):
    scene.begin_chain(QPointF(0, 0))
    assert scene.commit_chain(QPointF(3 * AXIS_STEP, 0)) is True

    args, kwargs = scene.add_molecule_fragment.call_args
    points, bonds_info = args
    assert len(points) == 4
    assert bonds_info == [(0, 1, 1), (1, 2, 1), (2, 3, 1)]
    assert kwargs["existing_items"] == []
    assert kwargs["symbol"] == "C"


def test_commit_passes_the_start_atom_for_fusing(scene):
    atom = MagicMock()
    atom.pos.return_value = QPointF(0, 0)
    scene.begin_chain(QPointF(0, 0), atom)
    scene.commit_chain(QPointF(2 * AXIS_STEP, 0))

    assert scene.add_molecule_fragment.call_args[1]["existing_items"] == [atom]


def test_commit_is_a_no_op_without_a_usable_drag(scene):
    assert scene.commit_chain(QPointF(50, 0)) is False

    scene.begin_chain(QPointF(0, 0))
    assert scene.commit_chain(QPointF(2, 0)) is False
    scene.add_molecule_fragment.assert_not_called()


def test_clear_chain_preview_resets_the_drag_state(scene):
    scene.begin_chain(QPointF(0, 0))
    scene.update_chain_preview(QPointF(2 * AXIS_STEP, 0))
    scene.clear_chain_preview()

    assert scene.chain_active is False
    assert scene.chain_anchor is None
    assert scene.chain_start_atom is None
    assert scene.chain_points == []
    scene.template_preview.hide.assert_called()
