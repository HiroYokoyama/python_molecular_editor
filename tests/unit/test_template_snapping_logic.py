"""Unit tests for template snapping alignment in TemplateMixin and KeyboardMixin."""

from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QKeyEvent
from moleditpy.ui.molecular_scene_handler import TemplateMixin, KeyboardMixin


class MockTemplateScene(TemplateMixin):
    def __init__(self):
        self.mode = "template_benzene"
        self.template_preview = MagicMock()
        self.views = MagicMock(return_value=[MagicMock()])
        self.settings = {}
        self.find_atom_near_args = []
        self.items_returned = []

    def get_setting(self, key, default=None):
        return self.settings.get(key, default)

    def find_atom_near(self, pos, tol=14.0):
        self.find_atom_near_args.append((pos, tol))
        return None

    def items(self, *args):
        return self.items_returned

    def _calculate_polygon_from_edge(
        self, p0, p1, n, cursor_pos=None, use_existing_length=False
    ):
        return [QPointF(0, 0)] * 6


class MockKeyboardScene(KeyboardMixin):
    def __init__(self):
        self.mode = "select"
        self.window = MagicMock()
        self.window.ui_manager.is_2d_editable = True
        self.views = MagicMock(return_value=[MagicMock()])
        # Mock view mapping
        view = self.views()[0]
        view.mapToScene.return_value = QPointF(100, 100)
        view.mapFromGlobal.return_value = QPointF(0, 0)
        view.transform.return_value = MagicMock()

        self.settings = {}
        self.find_atom_near_args = []
        self.update_all_items = MagicMock()
        self._calculate_polygon_from_edge = MagicMock(return_value=[QPointF(0, 0)] * 6)
        self.add_molecule_fragment = MagicMock()
        self.selectedItems = MagicMock(return_value=[])

    def get_setting(self, key, default=None):
        return self.settings.get(key, default)

    def find_atom_near(self, pos, tol=14.0):
        self.find_atom_near_args.append((pos, tol))
        return None

    def itemAt(self, pos, transform):
        return None


def test_update_template_preview_uses_template_snapping_distance(app):
    """update_template_preview calls find_atom_near with template_snapping_distance_2d."""
    scene = MockTemplateScene()
    scene.settings["template_snapping_distance_2d"] = 25.0
    scene.settings["template_fusing_distance_2d"] = 5.0

    scene.update_template_preview(QPointF(50, 50))

    # Assert that find_atom_near was called with tol=25.0
    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 25.0


def test_update_user_template_preview_uses_template_snapping_distance(app):
    """update_user_template_preview calls find_atom_near with template_snapping_distance_2d."""
    scene = MockTemplateScene()
    scene.settings["template_snapping_distance_2d"] = 30.0
    scene.settings["template_fusing_distance_2d"] = 5.0
    scene.user_template_data = {"atoms": [{"id": 0, "x": 0.0, "y": 0.0}], "bonds": []}

    scene.update_user_template_preview(QPointF(50, 50))

    # Assert that find_atom_near was called with tol=30.0
    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 30.0


def test_keyboard_key_4_uses_template_snapping_distance(app):
    """Key 4 (template ring) uses template_snapping_distance_2d for atom snap."""
    scene = MockKeyboardScene()
    scene.settings["template_snapping_distance_2d"] = 22.0
    scene.settings["template_fusing_distance_2d"] = 5.0

    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_4
    event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0, 0)
        scene.keyPressEvent(event)

    # Assert find_atom_near was called with 22.0 (template snap)
    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 22.0


def test_keyboard_key_4_uses_template_snapping_when_fusing_disabled(app):
    """Key 4 still uses template_snapping_distance_2d even when fusing is disabled."""
    scene = MockKeyboardScene()
    scene.settings["template_fusing_enabled_2d"] = False
    scene.settings["template_snapping_distance_2d"] = 22.0
    scene.settings["template_fusing_distance_2d"] = 5.0

    event = MagicMock(spec=QKeyEvent)
    event.key.return_value = Qt.Key.Key_4
    event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0, 0)
        scene.keyPressEvent(event)

    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 22.0


def test_keyboard_other_keys_use_bond_snapping_distance(app):
    """Non-template key presses use bond_snapping_distance_2d for atom snap."""
    scene = MockKeyboardScene()
    scene.settings["template_snapping_distance_2d"] = 22.0
    scene.settings["bond_snapping_distance_2d"] = 8.0

    event = MagicMock(spec=QKeyEvent)
    # Qt.Key.Key_C or similar to trigger atomic action
    event.key.return_value = Qt.Key.Key_C
    event.modifiers.return_value = Qt.KeyboardModifier.NoModifier

    with patch("moleditpy.ui.molecular_scene_handler.QCursor.pos") as mock_cursor:
        mock_cursor.return_value = QPointF(0, 0)
        scene.keyPressEvent(event)

    # Assert find_atom_near was called with 8.0 (bond snapping)
    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 8.0


# ---------------------------------------------------------------------------
# Ring template bond-length tests
# ---------------------------------------------------------------------------

import math
from moleditpy.utils.constants import DEFAULT_BOND_LENGTH


def _edge_lengths(points):
    """Return a list of distances between consecutive points (cyclic)."""
    n = len(points)
    return [
        math.hypot(
            points[(i + 1) % n].x() - points[i].x(),
            points[(i + 1) % n].y() - points[i].y(),
        )
        for i in range(n)
    ]


class _FreeTemplateMixin(TemplateMixin):
    """Minimal TemplateMixin subclass that captures placed points."""

    def __init__(self, n: int):
        self.mode = f"template_{n}"
        self.template_preview = MagicMock()
        self.views = MagicMock(return_value=[MagicMock()])
        self.settings = {}
        self.captured_points = None

    def get_setting(self, key, default=None):
        return self.settings.get(key, default)

    def find_atom_near(self, pos, tol=14.0):
        return None  # no snap target → free placement

    def items(self, *args):
        return []  # no snap target → free placement


def _free_ring_points(n: int) -> list:
    """Call update_template_preview with no snap target and return points."""
    scene = _FreeTemplateMixin(n)
    scene.update_template_preview(QPointF(0, 0))
    return scene.template_preview.set_geometry.call_args[0][0]


def _assert_all_edges_equal_bond_length(n: int, tol: float = 1e-6) -> None:
    points = _free_ring_points(n)
    assert len(points) == n, f"{n}-ring: expected {n} points, got {len(points)}"
    for length in _edge_lengths(points):
        assert abs(length - DEFAULT_BOND_LENGTH) < tol, (
            f"{n}-ring edge {length:.4f} != DEFAULT_BOND_LENGTH {DEFAULT_BOND_LENGTH}"
        )


def test_free_ring_3_bond_length(app):
    """3-ring (cyclopropane) free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(3)


def test_free_ring_4_bond_length(app):
    """4-ring free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(4)


def test_free_ring_5_bond_length(app):
    """5-ring free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(5)


def test_free_ring_6_bond_length(app):
    """6-ring (benzene) free placement: all edges == DEFAULT_BOND_LENGTH (regression guard)."""
    _assert_all_edges_equal_bond_length(6)


def test_free_ring_7_bond_length(app):
    """7-ring free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(7)


def test_free_ring_8_bond_length(app):
    """8-ring free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(8)


def test_free_ring_9_bond_length(app):
    """9-ring free placement: all edges == DEFAULT_BOND_LENGTH."""
    _assert_all_edges_equal_bond_length(9)


def test_free_ring_circumradius_formula():
    """Circumradius R = L / (2·sin(π/n)) gives correct edge for each ring size."""
    L = DEFAULT_BOND_LENGTH
    for n in range(3, 10):
        R = L / (2 * math.sin(math.pi / n))
        angle_step = 2 * math.pi / n
        # Two adjacent vertices on the unit circle separated by angle_step
        edge = 2 * R * math.sin(angle_step / 2)
        assert abs(edge - L) < 1e-9, f"n={n}: computed edge {edge} != {L}"
