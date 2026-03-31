import os
import sys
import pytest
from unittest.mock import MagicMock
from PyQt6.QtWidgets import QDialog

# Ensure local moleditpy is discoverable
workspace_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(workspace_src_path) and workspace_src_path not in sys.path:
    sys.path.insert(0, workspace_src_path)

from moleditpy.ui.edit_actions_logic import Rotate2DDialog

def test_rotate_2d_dialog_init(window, qtbot):
    """Test Rotate2DDialog GUI initialization in the GUI test environment."""
    # Use integer for initial_angle as QSpinBox expects int
    dialog = Rotate2DDialog(parent=window, initial_angle=45)
    qtbot.add_widget(dialog)
    dialog.show()
    
    assert dialog.windowTitle() == "Rotate 2D"
    assert dialog.angle_spin.value() == 45
    assert dialog.slider.value() == 45
    
    # Test sync: spinbox -> slider
    dialog.angle_spin.setValue(90)
    assert dialog.slider.value() == 90
    
    # Test sync: slider -> spinbox
    dialog.slider.setValue(-30)
    assert dialog.angle_spin.value() == -30
    
    # get_angle returns int (cast to float by type hint)
    assert dialog.get_angle() == -30.0
    dialog.close()
