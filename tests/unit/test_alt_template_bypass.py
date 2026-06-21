"""Unit tests for alt-key template bypass in TemplateMixin."""
import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtWidgets import QApplication
from moleditpy.ui.molecular_scene_handler import TemplateMixin


class MockTemplateScene(TemplateMixin):
    def __init__(self):
        self.mode = "template_benzene"
        self.template_preview = MagicMock()
        self.views = MagicMock(return_value=[MagicMock()])
        self.settings = {}
        self.find_atom_near_args = []
        self.items_returned = []
        self.data = MagicMock()
        self.data.atoms = {}
        self.atom_items = {}
        self.bond_items = {}

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
        return [QPointF(10, 20)] * n


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def test_update_template_preview_allows_snapping_but_bypasses_fusing_when_alt_pressed(
    qapp,
):
    scene = MockTemplateScene()
    scene.settings["template_snapping_distance_2d"] = 25.0
    scene.settings["template_fusing_enabled_2d"] = True
    scene.settings["template_fusing_distance_2d"] = 15.0

    # Mock a snapped atom item
    from moleditpy.ui.atom_item import AtomItem

    snapped_atom = MagicMock(spec=AtomItem)
    snapped_atom.pos.return_value = QPointF(0.0, 0.0)
    scene.items_returned = [snapped_atom]

    # Setup database with another atom close to expected preview points, but unclicked
    mock_atom = MagicMock()
    mock_atom.pos.return_value = QPointF(10.5, 20.5)
    scene.atom_items[101] = mock_atom

    # Simulate Alt key pressed
    with patch(
        "PyQt6.QtWidgets.QApplication.keyboardModifiers",
        return_value=Qt.KeyboardModifier.AltModifier,
    ):
        scene.update_template_preview(QPointF(5, 5))

    # Assert that find_atom_near WAS called because cursor snapping is allowed
    assert len(scene.find_atom_near_args) == 1

    # Assert that the preview did NOT snap other vertices to the existing atom at (10.5, 20.5)
    args, kwargs = scene.template_preview.set_geometry.call_args
    passed_points = args[0]
    # Point 1 would have snapped to (10.5, 20.5) if fusing were active. Because Alt is pressed, it stays at (10.0, 20.0)
    assert passed_points[1].x() == 10.0
    assert passed_points[1].y() == 20.0


def test_add_molecule_fragment_bypasses_fusing_when_alt_pressed(qapp):
    scene = MockTemplateScene()
    scene.settings["template_fusing_enabled_2d"] = True
    scene.settings["template_fusing_distance_2d"] = 15.0

    # Mock an atom in the scene close to one of the points
    mock_atom = MagicMock()
    mock_atom.pos.return_value = QPointF(10, 20)
    scene.atom_items[1] = mock_atom

    # Create dummy points for benzene placement
    points = [QPointF(10, 20)] + [QPointF(100, 100)] * 5
    bonds_info = [(0, 1, 1)]

    # Mock create_atom and create_bond; also set up atom_items so add_molecule_fragment can look up the new atom
    new_atom_mock = MagicMock()
    scene.create_atom = MagicMock(return_value=99)
    scene.atom_items[99] = new_atom_mock
    scene.create_bond = MagicMock()

    # Simulate Alt key pressed
    with patch(
        "PyQt6.QtWidgets.QApplication.keyboardModifiers",
        return_value=Qt.KeyboardModifier.AltModifier,
    ):
        scene.add_molecule_fragment(points, bonds_info)

    # If fusing were active, it would map points[0] to mock_atom (index 1), so create_atom wouldn't be called for it.
    # Because Alt is pressed, fusing is bypassed, so create_atom should be called for all 6 points.
    assert scene.create_atom.call_count == 6
