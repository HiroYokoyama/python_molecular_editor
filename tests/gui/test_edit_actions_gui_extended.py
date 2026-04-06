
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
