from moleditpy.ui.settings_tabs.settings_other_tab import SettingsOtherTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


def test_init_creates_all_controls(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    assert hasattr(tab, "skip_chem_checks_checkbox")
    assert hasattr(tab, "always_ask_charge_checkbox")
    assert hasattr(tab, "kekule_3d_checkbox")
    assert hasattr(tab, "aromatic_circle_checkbox")
    assert hasattr(tab, "aromatic_torus_thickness_slider")


def test_update_ui_sets_checkboxes(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["skip_chemistry_checks"] = True
    settings["always_ask_charge"] = True
    settings["display_kekule_3d"] = False
    settings["display_aromatic_circles_3d"] = False
    tab.update_ui(settings)
    assert tab.skip_chem_checks_checkbox.isChecked() is True
    assert tab.always_ask_charge_checkbox.isChecked() is True
    assert tab.kekule_3d_checkbox.isChecked() is False
    assert tab.aromatic_circle_checkbox.isChecked() is False


def test_update_ui_sets_torus_thickness(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["aromatic_torus_thickness_factor"] = 1.2
    tab.update_ui(settings)
    assert tab.aromatic_torus_thickness_slider.value() == 120


def test_get_settings_returns_correct_keys(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    result = tab.get_settings()
    expected_keys = {
        "skip_chemistry_checks",
        "always_ask_charge",
        "display_kekule_3d",
        "display_aromatic_circles_3d",
        "aromatic_torus_thickness_factor",
    }
    assert expected_keys == set(result.keys())


def test_get_settings_roundtrip_defaults(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert result["skip_chemistry_checks"] == DEFAULT_SETTINGS["skip_chemistry_checks"]
    assert result["display_kekule_3d"] == DEFAULT_SETTINGS["display_kekule_3d"]
    assert (
        result["display_aromatic_circles_3d"]
        == DEFAULT_SETTINGS["display_aromatic_circles_3d"]
    )
    assert (
        abs(
            result["aromatic_torus_thickness_factor"]
            - DEFAULT_SETTINGS["aromatic_torus_thickness_factor"]
        )
        < 0.01
    )


def test_kekule_toggled_disables_aromatic_circle(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.kekule_3d_checkbox.setChecked(True)
    assert tab.aromatic_circle_checkbox.isEnabled() is False


def test_kekule_untoggled_enables_aromatic_circle(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.kekule_3d_checkbox.setChecked(True)
    tab.kekule_3d_checkbox.setChecked(False)
    assert tab.aromatic_circle_checkbox.isEnabled() is True


def test_aromatic_toggled_disables_kekule(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.aromatic_circle_checkbox.setChecked(True)
    assert tab.kekule_3d_checkbox.isEnabled() is False


def test_aromatic_untoggled_enables_kekule(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.aromatic_circle_checkbox.setChecked(True)
    tab.aromatic_circle_checkbox.setChecked(False)
    assert tab.kekule_3d_checkbox.isEnabled() is True


def test_update_ui_kekule_true_disables_aromatic(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["display_kekule_3d"] = True
    settings["display_aromatic_circles_3d"] = False
    tab.update_ui(settings)
    assert tab.kekule_3d_checkbox.isChecked() is True
    assert tab.aromatic_circle_checkbox.isEnabled() is False


def test_update_ui_aromatic_true_disables_kekule(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["display_kekule_3d"] = False
    settings["display_aromatic_circles_3d"] = True
    tab.update_ui(settings)
    assert tab.aromatic_circle_checkbox.isChecked() is True
    assert tab.kekule_3d_checkbox.isEnabled() is False


def test_reset_to_defaults(app):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    custom = dict(DEFAULT_SETTINGS)
    custom["skip_chemistry_checks"] = True
    custom["aromatic_torus_thickness_factor"] = 2.0
    tab.update_ui(custom)
    tab.reset_to_defaults()
    result = tab.get_settings()
    assert result["skip_chemistry_checks"] == DEFAULT_SETTINGS["skip_chemistry_checks"]
    assert (
        abs(
            result["aromatic_torus_thickness_factor"]
            - DEFAULT_SETTINGS["aromatic_torus_thickness_factor"]
        )
        < 0.01
    )
