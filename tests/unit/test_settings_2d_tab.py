from unittest.mock import patch, MagicMock
from moleditpy.ui.settings_tabs.settings_2d_tab import Settings2DTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


def test_init_uses_default_colors(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    assert tab.current_bg_color_2d == DEFAULT_SETTINGS["background_color_2d"]
    assert tab.current_bond_color_2d == DEFAULT_SETTINGS["bond_color_2d"]


def test_update_ui_sets_colors(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["background_color_2d"] = "#112233"
    settings["bond_color_2d"] = "#aabbcc"
    tab.update_ui(settings)
    assert tab.current_bg_color_2d == "#112233"
    assert tab.current_bond_color_2d == "#aabbcc"


def test_update_ui_sets_sliders(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["bond_width_2d"] = 3.0
    settings["bond_spacing_double_2d"] = 4.0
    settings["bond_spacing_triple_2d"] = 5.0
    settings["bond_wedge_width_2d"] = 8.0
    settings["bond_dash_count_2d"] = 12
    settings["atom_font_size_2d"] = 24
    tab.update_ui(settings)
    assert tab.bond_width_2d_slider.value() == 30
    assert tab.bond_spacing_double_2d_slider.value() == 40
    assert tab.bond_spacing_triple_2d_slider.value() == 50
    assert tab.bond_wedge_width_2d_slider.value() == 80
    assert tab.bond_dash_count_2d_slider.value() == 12
    assert tab.atom_font_size_2d_slider.value() == 24


def test_update_ui_sets_cap_style(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["bond_cap_style_2d"] = "Flat"
    tab.update_ui(settings)
    assert tab.bond_cap_style_2d_combo.currentText() == "Flat"


def test_update_ui_sets_use_bond_color_checkbox(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["atom_use_bond_color_2d"] = True
    tab.update_ui(settings)
    assert tab.atom_use_bond_color_2d_checkbox.isChecked() is True


def test_get_settings_roundtrip(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert result["background_color_2d"] == DEFAULT_SETTINGS["background_color_2d"]
    assert result["bond_color_2d"] == DEFAULT_SETTINGS["bond_color_2d"]
    assert abs(result["bond_width_2d"] - DEFAULT_SETTINGS["bond_width_2d"]) < 0.05
    assert result["bond_cap_style_2d"] == DEFAULT_SETTINGS["bond_cap_style_2d"]
    assert result["bond_dash_count_2d"] == DEFAULT_SETTINGS["bond_dash_count_2d"]
    assert result["atom_font_size_2d"] == DEFAULT_SETTINGS["atom_font_size_2d"]
    assert result["atom_use_bond_color_2d"] == DEFAULT_SETTINGS["atom_use_bond_color_2d"]


def test_get_settings_returns_all_keys(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    result = tab.get_settings()
    expected_keys = {
        "background_color_2d", "bond_color_2d", "bond_width_2d",
        "bond_spacing_double_2d", "bond_spacing_triple_2d", "bond_cap_style_2d",
        "bond_wedge_width_2d", "bond_dash_count_2d", "atom_font_family_2d",
        "atom_font_size_2d", "atom_use_bond_color_2d",
    }
    assert expected_keys == set(result.keys())


@patch("moleditpy.ui.settings_tabs.settings_2d_tab.QColorDialog.getColor")
def test_pick_bg_color_updates_on_valid(mock_get_color, app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#ff0000"
    mock_get_color.return_value = mock_color
    tab._pick_bg_color_2d()
    assert tab.current_bg_color_2d == "#ff0000"


@patch("moleditpy.ui.settings_tabs.settings_2d_tab.QColorDialog.getColor")
def test_pick_bg_color_no_change_on_invalid(mock_get_color, app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    original = tab.current_bg_color_2d
    mock_color = MagicMock()
    mock_color.isValid.return_value = False
    mock_get_color.return_value = mock_color
    tab._pick_bg_color_2d()
    assert tab.current_bg_color_2d == original


@patch("moleditpy.ui.settings_tabs.settings_2d_tab.QColorDialog.getColor")
def test_pick_bond_color_updates_on_valid(mock_get_color, app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#00ff00"
    mock_get_color.return_value = mock_color
    tab._pick_bond_color_2d()
    assert tab.current_bond_color_2d == "#00ff00"


def test_reset_to_defaults(app):
    tab = Settings2DTab(DEFAULT_SETTINGS)
    custom = dict(DEFAULT_SETTINGS)
    custom["background_color_2d"] = "#123456"
    custom["bond_cap_style_2d"] = "Square"
    tab.update_ui(custom)
    assert tab.current_bg_color_2d == "#123456"

    tab.reset_to_defaults()
    assert tab.current_bg_color_2d == DEFAULT_SETTINGS["background_color_2d"]
    assert tab.bond_cap_style_2d_combo.currentText() == DEFAULT_SETTINGS["bond_cap_style_2d"]
