"""Opening Settings and pressing OK must not change any value.

Slider-backed settings are stored as slider / 100. Loading them back with
int(x * 100) truncated values such as 0.29 (28.999... -> 28), so every
open-and-OK of the dialog silently lowered them by one step.
"""

import pytest

from moleditpy.ui.settings_tabs.settings_3d_tabs import (
    Settings3DSceneTab,
    SettingsModelTab,
)
from moleditpy.ui.settings_tabs.settings_other_tab import SettingsOtherTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS

# Values whose x * 100 falls just below the integer in floating point.
DRIFTING = [0.29, 0.57, 0.58]


def _roundtrip(tab, key, value):
    settings = dict(DEFAULT_SETTINGS)
    settings[key] = value
    tab.update_ui(settings)
    return tab.get_settings()[key]


@pytest.mark.parametrize("value", DRIFTING)
@pytest.mark.parametrize(
    "key", ["light_intensity", "specular", "mouse_rotation_sensitivity"]
)
def test_scene_tab_keeps_value(app, key, value):
    tab = Settings3DSceneTab(DEFAULT_SETTINGS)
    assert _roundtrip(tab, key, value) == pytest.approx(value)


@pytest.mark.parametrize(
    "prefix, suffix, value",
    [
        ("ball_stick", "atom_scale", 0.57),
        ("ball_stick", "bond_radius", 0.29),
        ("stick", "double_bond_offset_factor", 1.13),
        ("stick", "triple_bond_radius_factor", 0.58),
    ],
)
def test_model_tab_keeps_value(app, prefix, suffix, value):
    tab = SettingsModelTab(prefix, "", DEFAULT_SETTINGS)
    assert _roundtrip(tab, f"{prefix}_{suffix}", value) == pytest.approx(value)


@pytest.mark.parametrize("value", DRIFTING)
def test_other_tab_keeps_torus_thickness(app, value):
    tab = SettingsOtherTab(DEFAULT_SETTINGS)
    assert _roundtrip(tab, "aromatic_torus_thickness_factor", value) == (
        pytest.approx(value)
    )
