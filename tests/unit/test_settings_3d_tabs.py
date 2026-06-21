"""Unit tests for 3D scene and model settings tabs."""

from unittest.mock import patch, MagicMock
from moleditpy.ui.settings_tabs.settings_3d_tabs import (
    Settings3DSceneTab,
    SettingsModelTab,
)
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


# --- Settings3DSceneTab ---


def test_scene_tab_init(app):
    """Settings3DSceneTab initialises with the default background color."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    assert tab.current_bg_color == DEFAULT_SETTINGS["background_color"]


def test_scene_tab_update_ui(app):
    """update_ui applies all provided settings to the scene tab widgets."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["background_color"] = "#112233"
    settings["show_3d_axes"] = False
    settings["lighting_enabled"] = False
    settings["light_intensity"] = 0.5
    settings["specular"] = 0.1
    settings["specular_power"] = 40
    settings["projection_mode"] = "Orthographic"
    tab.update_ui(settings)
    assert tab.current_bg_color == "#112233"
    assert tab.axes_checkbox.isChecked() is False
    assert tab.light_checkbox.isChecked() is False
    assert tab.intensity_slider.value() == 50
    assert tab.specular_slider.value() == 10
    assert tab.spec_power_slider.value() == 40
    assert tab.projection_combo.currentText() == "Orthographic"


def test_scene_tab_get_settings_keys(app):
    """get_settings returns a dict with all expected scene setting keys."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    result = tab.get_settings()
    expected_keys = {
        "background_color",
        "show_3d_axes",
        "lighting_enabled",
        "light_intensity",
        "specular",
        "specular_power",
        "projection_mode",
    }
    assert expected_keys == set(result.keys())


def test_scene_tab_roundtrip(app):
    """Settings round-trip: update_ui then get_settings recovers the original values."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert result["background_color"] == DEFAULT_SETTINGS["background_color"]
    assert result["show_3d_axes"] == DEFAULT_SETTINGS["show_3d_axes"]
    assert result["projection_mode"] == DEFAULT_SETTINGS["projection_mode"]
    assert abs(result["light_intensity"] - DEFAULT_SETTINGS["light_intensity"]) < 0.01
    assert abs(result["specular"] - DEFAULT_SETTINGS["specular"]) < 0.01
    assert result["specular_power"] == DEFAULT_SETTINGS["specular_power"]


@patch("moleditpy.ui.settings_tabs.settings_3d_tabs.QColorDialog.getColor")
def test_scene_tab_pick_color_valid(mock_get_color, app):
    """_select_color updates current_bg_color when a valid color is chosen."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#abcdef"
    mock_get_color.return_value = mock_color
    tab._select_color()
    assert tab.current_bg_color == "#abcdef"


@patch("moleditpy.ui.settings_tabs.settings_3d_tabs.QColorDialog.getColor")
def test_scene_tab_pick_color_invalid_no_change(mock_get_color, app):
    """_select_color leaves current_bg_color unchanged when an invalid color is returned."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    original = tab.current_bg_color
    mock_color = MagicMock()
    mock_color.isValid.return_value = False
    mock_get_color.return_value = mock_color
    tab._select_color()
    assert tab.current_bg_color == original


def test_scene_tab_reset_to_defaults(app):
    """reset_to_defaults restores projection_mode to the default setting."""
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    custom = dict(DEFAULT_SETTINGS)
    custom["projection_mode"] = "Orthographic"
    tab.update_ui(custom)
    assert tab.projection_combo.currentText() == "Orthographic"
    tab.reset_to_defaults()
    assert tab.projection_combo.currentText() == DEFAULT_SETTINGS["projection_mode"]


# --- SettingsModelTab (ball_stick) ---


def test_model_tab_ball_stick_has_atom_scale(app):
    """ball_stick SettingsModelTab creates an atom_scale_slider widget."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    assert hasattr(tab, "atom_scale_slider")


def test_model_tab_ball_stick_has_bond_radius(app):
    """ball_stick SettingsModelTab creates a bond_radius_slider widget."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    assert hasattr(tab, "bond_radius_slider")


def test_model_tab_ball_stick_has_color_options(app):
    """ball_stick SettingsModelTab creates bond_color_button and use_cpk_checkbox."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    assert hasattr(tab, "bond_color_button")
    assert hasattr(tab, "use_cpk_checkbox")


def test_model_tab_ball_stick_update_ui(app):
    """update_ui sets all ball-and-stick sliders and checkbox from settings."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["ball_stick_atom_scale"] = 1.5
    settings["ball_stick_bond_radius"] = 0.2
    settings["ball_stick_resolution"] = 20
    settings["ball_stick_use_cpk_bond_color"] = True
    tab.update_ui(settings)
    assert tab.atom_scale_slider.value() == 150
    assert tab.bond_radius_slider.value() == 20
    assert tab.res_slider.value() == 20
    assert tab.use_cpk_checkbox.isChecked() is True


def test_model_tab_ball_stick_get_settings_keys(app):
    """get_settings returns all expected ball-and-stick setting keys."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    result = tab.get_settings()
    expected_keys = {
        "ball_stick_atom_scale",
        "ball_stick_bond_radius",
        "ball_stick_double_bond_offset_factor",
        "ball_stick_triple_bond_offset_factor",
        "ball_stick_double_bond_radius_factor",
        "ball_stick_triple_bond_radius_factor",
        "ball_stick_resolution",
        "ball_stick_bond_color",
        "ball_stick_use_cpk_bond_color",
    }
    assert expected_keys == set(result.keys())


def test_model_tab_ball_stick_roundtrip(app):
    """ball_stick settings round-trip recovers atom scale, bond radius, and resolution."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert (
        abs(result["ball_stick_atom_scale"] - DEFAULT_SETTINGS["ball_stick_atom_scale"])
        < 0.01
    )
    assert (
        abs(
            result["ball_stick_bond_radius"]
            - DEFAULT_SETTINGS["ball_stick_bond_radius"]
        )
        < 0.01
    )
    assert result["ball_stick_resolution"] == DEFAULT_SETTINGS["ball_stick_resolution"]


# --- SettingsModelTab (cpk) ---


def test_model_tab_cpk_has_atom_scale_no_bond_radius(app):
    """cpk SettingsModelTab has atom_scale_slider but not bond_radius_slider."""
    tab = SettingsModelTab("cpk", "info", DEFAULT_SETTINGS)
    assert hasattr(tab, "atom_scale_slider")
    assert not hasattr(tab, "bond_radius_slider")


def test_model_tab_cpk_get_settings_keys(app):
    """cpk get_settings returns cpk_atom_scale and cpk_resolution, not ball-stick keys."""
    tab = SettingsModelTab("cpk", "info", DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert "cpk_atom_scale" in result
    assert "cpk_resolution" in result
    assert "ball_stick_bond_color" not in result


def test_model_tab_cpk_roundtrip(app):
    """cpk settings round-trip recovers atom scale and resolution."""
    tab = SettingsModelTab("cpk", "info", DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert abs(result["cpk_atom_scale"] - DEFAULT_SETTINGS["cpk_atom_scale"]) < 0.01
    assert result["cpk_resolution"] == DEFAULT_SETTINGS["cpk_resolution"]


# --- SettingsModelTab (wireframe) ---


def test_model_tab_wireframe_no_atom_scale(app):
    """wireframe SettingsModelTab has bond_radius_slider but not atom_scale_slider."""
    tab = SettingsModelTab("wireframe", "info", DEFAULT_SETTINGS)
    assert not hasattr(tab, "atom_scale_slider")
    assert hasattr(tab, "bond_radius_slider")


def test_model_tab_wireframe_get_settings_keys(app):
    """wireframe get_settings returns wireframe_bond_radius, resolution, and offset keys."""
    tab = SettingsModelTab("wireframe", "info", DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert "wireframe_bond_radius" in result
    assert "wireframe_resolution" in result
    assert "wireframe_double_bond_offset_factor" in result


def test_model_tab_wireframe_roundtrip(app):
    """wireframe settings round-trip recovers bond radius and resolution."""
    tab = SettingsModelTab("wireframe", "info", DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert (
        abs(result["wireframe_bond_radius"] - DEFAULT_SETTINGS["wireframe_bond_radius"])
        < 0.01
    )
    assert result["wireframe_resolution"] == DEFAULT_SETTINGS["wireframe_resolution"]


# --- SettingsModelTab (stick) ---


def test_model_tab_stick_get_settings_keys(app):
    """stick get_settings returns stick_bond_radius, resolution, and offset keys."""
    tab = SettingsModelTab("stick", "info", DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert "stick_bond_radius" in result
    assert "stick_resolution" in result
    assert "stick_double_bond_offset_factor" in result


def test_model_tab_stick_roundtrip(app):
    """stick settings round-trip recovers bond radius and resolution."""
    tab = SettingsModelTab("stick", "info", DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert (
        abs(result["stick_bond_radius"] - DEFAULT_SETTINGS["stick_bond_radius"]) < 0.01
    )
    assert result["stick_resolution"] == DEFAULT_SETTINGS["stick_resolution"]


@patch("moleditpy.ui.settings_tabs.settings_3d_tabs.QColorDialog.getColor")
def test_model_tab_ball_stick_pick_bond_color(mock_get_color, app):
    """_pick_bond_color updates current_bond_color when a valid color is chosen."""
    tab = SettingsModelTab("ball_stick", "info", DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#ff00ff"
    mock_get_color.return_value = mock_color
    tab._pick_bond_color()
    assert tab.current_bond_color == "#ff00ff"
