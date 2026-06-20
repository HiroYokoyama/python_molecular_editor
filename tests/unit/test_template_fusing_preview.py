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
        self.template_context = {}

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
        # Return points centered roughly around (0,0)
        return [
            QPointF(0.0, 0.0),
            QPointF(10.0, 0.0),
            QPointF(15.0, 10.0),
            QPointF(10.0, 20.0),
            QPointF(0.0, 20.0),
            QPointF(-5.0, 10.0),
        ]


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def test_update_template_preview_snaps_points_for_fusing(qapp):
    scene = MockTemplateScene()
    scene.settings["template_snapping_distance_2d"] = 15.0
    scene.settings["template_fusing_enabled_2d"] = True
    scene.settings["template_fusing_distance_2d"] = 10.0

    # Mock a snapped atom item
    from moleditpy.ui.atom_item import AtomItem

    snapped_atom = MagicMock(spec=AtomItem)
    snapped_atom.pos.return_value = QPointF(0.0, 0.0)
    scene.items_returned = [snapped_atom]

    # Place an existing atom in the scene at (10.5, 0.5), which is close to point 1 (10.0, 0.0)
    mock_atom = MagicMock()
    mock_atom.pos.return_value = QPointF(10.5, 0.5)
    scene.atom_items[101] = mock_atom

    # Call update_template_preview (simulating mouse move)
    with patch(
        "PyQt6.QtWidgets.QApplication.keyboardModifiers",
        return_value=Qt.KeyboardModifier.NoModifier,
    ):
        scene.update_template_preview(QPointF(5, 5))

    # Get the points passed to set_geometry
    args, kwargs = scene.template_preview.set_geometry.call_args
    passed_points = args[0]

    # Verify that point 0 remained (0.0, 0.0) as it is mapped to snapped_atom
    assert passed_points[0].x() == 0.0
    assert passed_points[0].y() == 0.0

    # Verify that point 1 snapped to the mock_atom position (10.5, 0.5) to reflect fusing in preview
    assert passed_points[1].x() == 10.5
    assert passed_points[1].y() == 0.5
