"""Unit tests for SettingsDialog multi-tab layout and settings flow."""

from unittest.mock import patch, MagicMock
from moleditpy.ui.settings_dialog import SettingsDialog
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


def _make_parent():
    """Build a minimal mock of the parent MainWindow for SettingsDialog."""
    parent = MagicMock()
    parent.init_manager.settings = dict(DEFAULT_SETTINGS)
    parent.init_manager.settings_dirty = False
    parent.view_3d_manager.current_mol = None
    parent.statusBar.return_value = MagicMock()
    return parent


def test_init_creates_seven_tabs(app):
    """SettingsDialog initialises with exactly seven tabs."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    assert dialog.tab_widget.count() == 7


def test_init_tab_labels(app):
    """Tab titles include the expected section names."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    tab_titles = [
        dialog.tab_widget.tabText(i) for i in range(dialog.tab_widget.count())
    ]
    assert "2D Settings" in tab_titles
    assert "3D Scene" in tab_titles
    assert "Ball & Stick" in tab_titles
    assert "Other" in tab_titles


def test_init_window_title(app):
    """SettingsDialog window title is 'Settings'."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    assert dialog.windowTitle() == "Settings"


def test_update_ui_from_settings_sets_all_tabs(app):
    """update_ui_from_settings propagates settings values to each tab widget."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    custom = dict(DEFAULT_SETTINGS)
    custom["background_color_2d"] = "#aabbcc"
    dialog.update_ui_from_settings(custom)
    assert dialog.tab_2d.current_bg_color_2d == "#aabbcc"


def test_get_settings_aggregates_all_tabs(app):
    """get_settings aggregates keys from every settings tab."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    result = dialog.get_settings()
    # Should contain keys from every tab
    assert "background_color_2d" in result  # 2D tab
    assert "background_color" in result  # 3D scene tab
    assert "ball_stick_atom_scale" in result  # Ball & Stick tab
    assert "cpk_atom_scale" in result  # CPK tab
    assert "wireframe_bond_radius" in result  # Wireframe tab
    assert "stick_bond_radius" in result  # Stick tab
    assert "skip_chemistry_checks" in result  # Other tab


def test_reset_current_tab_calls_reset_on_active_tab(app):
    """reset_current_tab calls reset_to_defaults on the currently visible tab."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.tab_widget.setCurrentIndex(0)

    with (
        patch.object(dialog.tab_2d, "reset_to_defaults") as mock_reset,
        patch("moleditpy.ui.settings_dialog.QMessageBox.information"),
    ):
        dialog.reset_current_tab()
        mock_reset.assert_called_once()


def test_reset_current_tab_shows_info_message(app):
    """reset_current_tab shows an information message naming the reset tab."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.tab_widget.setCurrentIndex(0)

    with (
        patch.object(dialog.tab_2d, "reset_to_defaults"),
        patch("moleditpy.ui.settings_dialog.QMessageBox.information") as mock_info,
    ):
        dialog.reset_current_tab()
        mock_info.assert_called_once()
        args = mock_info.call_args[0]
        assert "2D Settings" in args[2]


def test_reset_all_settings_yes_resets_and_applies(app):
    """reset_all_settings resets UI and applies when user confirms Yes."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = _make_parent()

    with (
        patch(
            "moleditpy.ui.settings_dialog.QMessageBox.question",
            return_value=MagicMock(),
        ) as mock_q,
        patch.object(dialog, "update_ui_from_settings") as mock_update,
        patch.object(dialog, "apply_settings") as mock_apply,
    ):
        from PyQt6.QtWidgets import QMessageBox

        mock_q.return_value = QMessageBox.StandardButton.Yes
        dialog.reset_all_settings()
        mock_update.assert_called_once_with(dialog.default_settings)
        mock_apply.assert_called_once()


def test_reset_all_settings_no_does_nothing(app):
    """reset_all_settings is a no-op when user answers No."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = _make_parent()

    with (
        patch("moleditpy.ui.settings_dialog.QMessageBox.question") as mock_q,
        patch.object(dialog, "update_ui_from_settings") as mock_update,
    ):
        from PyQt6.QtWidgets import QMessageBox

        mock_q.return_value = QMessageBox.StandardButton.No
        dialog.reset_all_settings()
        mock_update.assert_not_called()


def test_apply_settings_no_parent_returns_early(app):
    """apply_settings is a no-op when parent_window is None."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = None
    # Should not raise
    dialog.apply_settings()


def test_apply_settings_updates_parent_settings(app):
    """apply_settings pushes current tab values into the parent's settings dict."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = _make_parent()
    dialog.tab_2d.current_bg_color_2d = "#123456"
    dialog.apply_settings()
    assert (
        dialog.parent_window.init_manager.settings["background_color_2d"] == "#123456"
    )


def test_apply_settings_fires_all_side_effects(app):
    """apply_settings calls save_settings, apply_3d_settings, update_cpk_colors, and status bar."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = _make_parent()
    dialog.apply_settings()
    p = dialog.parent_window
    p.init_manager.save_settings.assert_called()
    p.view_3d_manager.apply_3d_settings.assert_called()
    p.init_manager.update_cpk_colors_from_settings.assert_called()
    p.statusBar.return_value.showMessage.assert_called_with(
        "Settings applied successfully"
    )


def test_apply_settings_redraws_molecule_when_present(app):
    """apply_settings calls draw_molecule_3d when a molecule is currently loaded."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    parent = _make_parent()
    mol = MagicMock()
    parent.view_3d_manager.current_mol = mol
    dialog.parent_window = parent
    dialog.apply_settings()
    parent.view_3d_manager.draw_molecule_3d.assert_called_with(mol)


def test_apply_settings_skips_redraw_when_no_molecule(app):
    """apply_settings skips draw_molecule_3d when no molecule is loaded."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    parent = _make_parent()
    parent.view_3d_manager.current_mol = None
    dialog.parent_window = parent
    dialog.apply_settings()
    parent.view_3d_manager.draw_molecule_3d.assert_not_called()


def test_apply_settings_updates_2d_scene_background(app):
    """apply_settings calls setBackgroundBrush on the 2D scene."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    parent = _make_parent()
    scene = MagicMock()
    scene.items.return_value = []
    parent.init_manager.scene = scene
    dialog.parent_window = parent
    dialog.apply_settings()
    scene.setBackgroundBrush.assert_called()


def test_accept_applies_then_closes(app):
    """accept() calls apply_settings before delegating to QDialog.accept."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.parent_window = _make_parent()

    with (
        patch.object(dialog, "apply_settings") as mock_apply,
        patch("moleditpy.ui.settings_dialog.QDialog.accept") as mock_super_accept,
    ):
        dialog.accept()
        mock_apply.assert_called_once()
        mock_super_accept.assert_called_once()


def test_get_settings_roundtrip_defaults(app):
    """get_settings after loading DEFAULT_SETTINGS round-trips the default values."""
    dialog = SettingsDialog(DEFAULT_SETTINGS, parent=None)
    dialog.update_ui_from_settings(DEFAULT_SETTINGS)
    result = dialog.get_settings()
    assert (
        abs(result["ball_stick_atom_scale"] - DEFAULT_SETTINGS["ball_stick_atom_scale"])
        < 0.01
    )
    assert result["background_color"] == DEFAULT_SETTINGS["background_color"]
    assert result["skip_chemistry_checks"] == DEFAULT_SETTINGS["skip_chemistry_checks"]
