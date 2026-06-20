"""
Unit tests for ui/periodic_table_dialog.py.
"""

import pytest
from PyQt6.QtWidgets import QWidget
from moleditpy.ui.periodic_table_dialog import PeriodicTableDialog


class MockParent(QWidget):
    def __init__(self):
        super().__init__()
        self.settings = {}


class ErrorMockParent(QWidget):
    def __init__(self):
        super().__init__()

    @property
    def settings(self):
        raise AttributeError("Mock error")


def test_periodic_table_dialog_no_parent(app):
    """Test dialog creation without parent settings."""
    dialog = PeriodicTableDialog()
    assert dialog.windowTitle() == "Select an Element"
    assert dialog.layout() is not None


def test_periodic_table_dialog_with_parent_overrides(app):
    """Test dialog creation with parent settings cpk_colors overrides."""
    parent = MockParent()
    # Mock settings override for H and C
    parent.settings = {
        "cpk_colors": {
            "H": "#111111",
            "C": "#222222",
        }
    }

    dialog = PeriodicTableDialog(parent=parent)
    assert dialog.windowTitle() == "Select an Element"


def test_periodic_table_dialog_with_parent_error_handling(app):
    """Test dialog creation when parent attribute access raises exception."""
    parent = ErrorMockParent()

    dialog = PeriodicTableDialog(parent=parent)
    assert dialog.windowTitle() == "Select an Element"


def test_periodic_table_dialog_element_clicked(app, qtbot):
    """Test that clicking a button emits element_selected and accepts the dialog."""
    dialog = PeriodicTableDialog()

    # We will spy on the element_selected signal
    with qtbot.waitSignal(dialog.element_selected) as blocker:
        # Find a button to click, e.g., the "H" button
        h_button = None
        for i in range(dialog.layout().count()):
            widget = dialog.layout().itemAt(i).widget()
            if widget and widget.text() == "H":
                h_button = widget
                break

        assert h_button is not None
        # Simulate click
        qtbot.mouseClick(
            h_button, pytest.importorskip("PyQt6.QtCore").Qt.MouseButton.LeftButton
        )

    assert blocker.args == ["H"]
