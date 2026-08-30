"""Unit tests for Settings2DCleanupTab."""

from moleditpy.ui.settings_tabs.settings_2d_cleanup_tab import Settings2DCleanupTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


CHECKBOX_KEYS = (
    "prefer_coordgen_2d",
    "cleanup_canonical_orientation_2d",
    "cleanup_use_ring_templates_2d",
    "cleanup_straighten_bonds_2d",
    "cleanup_avoid_clashes_2d",
)


def _checkbox_by_key(tab):
    return {
        "prefer_coordgen_2d": tab.prefer_coordgen_2d_checkbox,
        "cleanup_canonical_orientation_2d": tab.cleanup_canonical_orientation_2d_checkbox,
        "cleanup_use_ring_templates_2d": tab.cleanup_use_ring_templates_2d_checkbox,
        "cleanup_straighten_bonds_2d": tab.cleanup_straighten_bonds_2d_checkbox,
        "cleanup_avoid_clashes_2d": tab.cleanup_avoid_clashes_2d_checkbox,
    }


def test_init_uses_default_settings(app):
    """update_ui(DEFAULT_SETTINGS) reflects DEFAULT_SETTINGS in each checkbox."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    checkbox_by_key = _checkbox_by_key(tab)
    for key in CHECKBOX_KEYS:
        assert checkbox_by_key[key].isChecked() == DEFAULT_SETTINGS[key]


def test_get_settings_returns_all_keys(app):
    """get_settings() returns exactly the five cleanup-option keys."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    result = tab.get_settings()
    assert set(result.keys()) == set(CHECKBOX_KEYS)


def test_checkboxes_round_trip(app):
    """update_ui()/get_settings() correctly reflect all five checkboxes."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)
    tab.update_ui(DEFAULT_SETTINGS)
    checkbox_by_key = _checkbox_by_key(tab)

    for key in CHECKBOX_KEYS:
        assert checkbox_by_key[key].isChecked() == DEFAULT_SETTINGS[key]

    custom = dict(DEFAULT_SETTINGS)
    for key in CHECKBOX_KEYS:
        custom[key] = not DEFAULT_SETTINGS[key]
    tab.update_ui(custom)

    result = tab.get_settings()
    for key in CHECKBOX_KEYS:
        assert checkbox_by_key[key].isChecked() == custom[key]
        assert result[key] == custom[key]

    tab.reset_to_defaults()
    for key in CHECKBOX_KEYS:
        assert checkbox_by_key[key].isChecked() == DEFAULT_SETTINGS[key]


def test_prefer_coordgen_disables_inapplicable_options(app):
    """Ring Templates and Avoid Clashes are disabled while CoordGen is preferred.

    Both options are silently ignored by RDKit's CoordGen algorithm, so they
    should be grayed out rather than implying they have an effect.
    """
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)

    tab.prefer_coordgen_2d_checkbox.setChecked(True)
    assert tab.cleanup_use_ring_templates_2d_checkbox.isEnabled() is False
    assert tab.cleanup_avoid_clashes_2d_checkbox.isEnabled() is False
    # Unaffected options stay enabled.
    assert tab.cleanup_canonical_orientation_2d_checkbox.isEnabled() is True
    assert tab.cleanup_straighten_bonds_2d_checkbox.isEnabled() is True

    tab.prefer_coordgen_2d_checkbox.setChecked(False)
    assert tab.cleanup_use_ring_templates_2d_checkbox.isEnabled() is True
    assert tab.cleanup_avoid_clashes_2d_checkbox.isEnabled() is True


def test_update_ui_sets_dependent_option_enabled_state(app):
    """update_ui() applies the disabled state even without a checkbox toggle."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)

    settings = dict(DEFAULT_SETTINGS)
    settings["prefer_coordgen_2d"] = True
    tab.update_ui(settings)
    assert tab.cleanup_use_ring_templates_2d_checkbox.isEnabled() is False
    assert tab.cleanup_avoid_clashes_2d_checkbox.isEnabled() is False

    settings["prefer_coordgen_2d"] = False
    tab.update_ui(settings)
    assert tab.cleanup_use_ring_templates_2d_checkbox.isEnabled() is True
    assert tab.cleanup_avoid_clashes_2d_checkbox.isEnabled() is True
