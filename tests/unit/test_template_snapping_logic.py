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
    scene = MockTemplateScene()
    scene.settings["template_snapping_distance_2d"] = 25.0
    scene.settings["template_fusing_distance_2d"] = 5.0

    scene.update_template_preview(QPointF(50, 50))

    # Assert that find_atom_near was called with tol=25.0
    assert len(scene.find_atom_near_args) == 1
    pos, tol = scene.find_atom_near_args[0]
    assert tol == 25.0


def test_update_user_template_preview_uses_template_snapping_distance(app):
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
