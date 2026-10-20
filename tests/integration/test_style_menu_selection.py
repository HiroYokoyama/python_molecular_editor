"""Exercise native Qt action exclusivity, including Python 3.9 / Qt 6.9."""

import gc

from PySide6.QtCore import QPoint, Qt
from PySide6.QtWidgets import QToolBar


def _style_actions(window):
    return [
        action
        for action in window.init_manager.style_button.menu().actions()
        if action.isCheckable()
    ]


def test_style_menu_clicks_keep_exactly_one_selection(window):
    gc.collect()
    actions = _style_actions(window)
    assert all(action.actionGroup() is not None for action in actions)
    for action in actions + list(reversed(actions)):
        action.trigger()
        assert [candidate for candidate in actions if candidate.isChecked()] == [action]


def test_programmatic_style_changes_update_menu(window):
    actions = _style_actions(window)
    for style, action in zip(["ball_and_stick", "cpk", "wireframe", "stick"], actions):
        window.view_3d_manager.set_3d_style(style)
        assert [candidate for candidate in actions if candidate.isChecked()] == [action]


def test_style_menu_mouse_selection(window, qtbot):
    menu = window.init_manager.style_button.menu()
    actions = _style_actions(window)
    for action in actions:
        menu.popup(QPoint(50, 50))
        qtbot.mouseClick(
            menu, Qt.MouseButton.LeftButton, pos=menu.actionGeometry(action).center()
        )
        assert [candidate for candidate in actions if candidate.isChecked()] == [action]


def test_plugin_style_joins_exclusive_selection(window):
    window.plugin_manager.custom_3d_styles["Test style"] = {"callback": lambda *_: None}
    window.init_manager.plugin_menu_manager.update_style_menu_with_plugins()
    gc.collect()
    actions = _style_actions(window)
    assert actions[-1].text() == "Test style"
    for action in actions:
        action.trigger()
        assert [candidate for candidate in actions if candidate.isChecked()] == [action]


def test_startup_plugin_style_joins_exclusive_selection(window):
    window.plugin_manager.custom_3d_styles["Startup style"] = {
        "callback": lambda *_: None
    }
    toolbar = QToolBar(window)
    window.init_manager._add_3d_edit_actions(toolbar)
    actions = _style_actions(window)
    assert actions[-1].text() == "Startup style"
    actions[-1].trigger()
    assert window.view_3d_manager.current_3d_style == "Startup style"
    assert [action for action in actions if action.isChecked()] == [actions[-1]]
    window.view_3d_manager.set_3d_style("Startup style")
    assert [action for action in actions if action.isChecked()] == [actions[-1]]
