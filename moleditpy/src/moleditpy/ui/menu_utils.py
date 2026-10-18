"""Native menu ownership helpers (LGPL-3.0-or-later)."""

from typing import Union

from PySide6.QtWidgets import QMenu, QMenuBar


def add_menu(parent: Union[QMenu, QMenuBar], title: str) -> QMenu:
    """Create a parent-owned menu whose wrapper survives garbage collection."""
    # Qt 6.9's addMenu(str) overload can invalidate the submenu's Python wrapper.
    # Constructing the widget with its parent avoids that ownership ambiguity.
    menu = QMenu(title, parent)
    parent.addMenu(menu)
    return menu
