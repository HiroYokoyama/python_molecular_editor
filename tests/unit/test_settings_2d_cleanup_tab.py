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


def _coordgen_ignored_checkboxes(tab):
    """The options RDKit's CoordGen silently ignores (verified against RDKit)."""
    return (
        tab.cleanup_canonical_orientation_2d_checkbox,
        tab.cleanup_use_ring_templates_2d_checkbox,
        tab.cleanup_avoid_clashes_2d_checkbox,
    )


def test_prefer_coordgen_disables_inapplicable_options(app):
    """Options CoordGen ignores are disabled while CoordGen is preferred.

    canonOrient, useRingTemplates and the clash-avoidance sampling params are
    all silently ignored by RDKit's CoordGen algorithm, so they should be
    grayed out rather than implying they have an effect. StraightenDepiction
    runs after the depiction is built, so it still applies.
    """
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)

    tab.prefer_coordgen_2d_checkbox.setChecked(True)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert checkbox.isEnabled() is False
    assert tab.cleanup_straighten_bonds_2d_checkbox.isEnabled() is True

    tab.prefer_coordgen_2d_checkbox.setChecked(False)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert checkbox.isEnabled() is True
    assert tab.cleanup_straighten_bonds_2d_checkbox.isEnabled() is True


def test_prefer_coordgen_annotates_inapplicable_tooltips(app):
    """Disabled options explain why they don't apply, and revert when re-enabled."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)
    suffix = "not applicable when Prefer CoordGen is checked"

    tab.prefer_coordgen_2d_checkbox.setChecked(True)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert suffix in checkbox.toolTip()

    tab.prefer_coordgen_2d_checkbox.setChecked(False)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert suffix not in checkbox.toolTip()
        assert checkbox.toolTip()  # base tooltip survived the round trip


def test_update_ui_sets_dependent_option_enabled_state(app):
    """update_ui() applies the disabled state even without a checkbox toggle."""
    tab = Settings2DCleanupTab(DEFAULT_SETTINGS)

    settings = dict(DEFAULT_SETTINGS)
    settings["prefer_coordgen_2d"] = True
    tab.update_ui(settings)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert checkbox.isEnabled() is False

    settings["prefer_coordgen_2d"] = False
    tab.update_ui(settings)
    for checkbox in _coordgen_ignored_checkboxes(tab):
        assert checkbox.isEnabled() is True
