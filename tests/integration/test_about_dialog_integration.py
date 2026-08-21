"""The About image's right-click shortcut has to reach the 3D view."""

from unittest.mock import MagicMock

from PyQt6.QtCore import Qt

from moleditpy.ui.about_dialog import AboutDialog


def test_right_click_loads_and_converts(window, qtbot):
    """Right-clicking the image loads the molecule and builds it in 3D."""
    dialog = AboutDialog(window, window)
    event = MagicMock()
    event.button.return_value = Qt.MouseButton.RightButton

    dialog.image_mouse_press_event(event)

    assert len(window.state_manager.data.atoms) == 12
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None, timeout=15000
    )
