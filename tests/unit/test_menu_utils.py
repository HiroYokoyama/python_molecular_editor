"""Submenus must remain usable after temporary Python references disappear."""

import gc

import pytest
from PySide6.QtWidgets import QMenu, QMenuBar
from shiboken6 import isValid

from moleditpy.ui.menu_utils import add_menu


@pytest.mark.parametrize("parent_type", [QMenu, QMenuBar])
def test_parent_owns_submenu_after_python_reference_is_released(app, parent_type):
    parent = parent_type()
    submenu = add_menu(parent, "Tools")
    submenu.addAction("Run")
    del submenu
    gc.collect()

    # Repeated access also releases temporary wrappers, as menu traversal does.
    assert parent.actions()[0].menu().title() == "Tools"
    gc.collect()
    submenu = parent.actions()[0].menu()
    assert isValid(submenu)
    assert submenu.parent() is parent
    assert [action.text() for action in submenu.actions()] == ["Run"]
    parent.deleteLater()
