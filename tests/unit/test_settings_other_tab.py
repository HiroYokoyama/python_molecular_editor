"""Unit tests for SettingsOtherTab."""

from moleditpy.ui.settings_tabs.settings_other_tab import SettingsOtherTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


def test_init_creates_all_controls(app):
    """SettingsOtherTab creates all expected control widgets on init."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    assert hasattr(tab, "skip_chem_checks_checkbox")
    assert hasattr(tab, "always_ask_charge_checkbox")
    assert hasattr(tab, "kekule_3d_checkbox")
    assert hasattr(tab, "aromatic_circle_checkbox")
    assert hasattr(tab, "aromatic_torus_thickness_slider")


def test_update_ui_sets_checkboxes(app):
    """update_ui sets checkbox states from provided settings dict."""
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
    """update_ui maps aromatic_torus_thickness_factor to slider integer value."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["aromatic_torus_thickness_factor"] = 1.2
    tab.update_ui(settings)
    assert tab.aromatic_torus_thickness_slider.value() == 120


def test_get_settings_returns_correct_keys(app):
    """get_settings returns a dict with exactly the expected setting keys."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    result = tab.get_settings()
    expected_keys = {
        "skip_chemistry_checks",
        "always_ask_charge",
        "display_kekule_3d",
        "display_aromatic_circles_3d",
        "aromatic_torus_thickness_factor",
        "log_to_file",
        "log_level_debug",
    }
    assert expected_keys == set(result.keys())


def test_get_settings_roundtrip_defaults(app):
    """get_settings after update_ui(DEFAULT_SETTINGS) round-trips the default values."""
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
    """Checking Kekulé disables the aromatic circle checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.kekule_3d_checkbox.setChecked(True)
    assert tab.aromatic_circle_checkbox.isEnabled() is False


def test_kekule_untoggled_enables_aromatic_circle(app):
    """Unchecking Kekulé re-enables the aromatic circle checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.kekule_3d_checkbox.setChecked(True)
    tab.kekule_3d_checkbox.setChecked(False)
    assert tab.aromatic_circle_checkbox.isEnabled() is True


def test_aromatic_toggled_disables_kekule(app):
    """Checking aromatic circles disables the Kekulé checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.aromatic_circle_checkbox.setChecked(True)
    assert tab.kekule_3d_checkbox.isEnabled() is False


def test_aromatic_untoggled_enables_kekule(app):
    """Unchecking aromatic circles re-enables the Kekulé checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.aromatic_circle_checkbox.setChecked(True)
    tab.aromatic_circle_checkbox.setChecked(False)
    assert tab.kekule_3d_checkbox.isEnabled() is True


def test_update_ui_kekule_true_disables_aromatic(app):
    """update_ui with display_kekule_3d=True disables the aromatic circle checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["display_kekule_3d"] = True
    settings["display_aromatic_circles_3d"] = False
    tab.update_ui(settings)
    assert tab.kekule_3d_checkbox.isChecked() is True
    assert tab.aromatic_circle_checkbox.isEnabled() is False


def test_update_ui_aromatic_true_disables_kekule(app):
    """update_ui with display_aromatic_circles_3d=True disables the Kekulé checkbox."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["display_kekule_3d"] = False
    settings["display_aromatic_circles_3d"] = True
    tab.update_ui(settings)
    assert tab.aromatic_circle_checkbox.isChecked() is True
    assert tab.kekule_3d_checkbox.isEnabled() is False


def test_reset_to_defaults(app):
    """reset_to_defaults restores settings to default values."""
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


def test_sync_slider_from_spinbox(app):
    """Setting the spinbox value syncs the torus thickness slider accordingly."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.aromatic_torus_thickness_label.setValue(1.50)
    assert tab.aromatic_torus_thickness_slider.value() == 150


def test_log_to_file_default_unchecked(app):
    """log_to_file checkbox is unchecked by default."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    assert tab.log_to_file_checkbox.isChecked() is False
    assert tab.log_level_debug_checkbox.isChecked() is False


def test_log_settings_roundtrip(app):
    """Enabling log_to_file and log_level_debug round-trips through get_settings."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["log_to_file"] = True
    settings["log_level_debug"] = True
    tab.update_ui(settings)
    result = tab.get_settings()
    assert result["log_to_file"] is True
    assert result["log_level_debug"] is True


def test_log_settings_off_roundtrip(app):
    """Disabling both log settings returns False in get_settings."""
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    settings = dict(DEFAULT_SETTINGS)
    settings["log_to_file"] = False
    settings["log_level_debug"] = False
    tab.update_ui(settings)
    result = tab.get_settings()
    assert result["log_to_file"] is False
    assert result["log_level_debug"] is False
