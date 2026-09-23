"""Unit tests for utils.label_style: 3D label style settings."""

from unittest.mock import MagicMock

import pytest

from moleditpy.utils.default_settings import DEFAULT_SETTINGS
from moleditpy.utils.label_style import (
    DEFAULT_LABEL_SETTINGS,
    LABEL_KINDS,
    LABEL_SECTIONS,
    color_key,
    label_kwargs,
)

# The color each label had before it could be styled.
PREVIOUS_COLOR = {
    "index": "#003366",
    "original_id": "#009000",
    "xyz_index": "#8B0000",
    "symbol": "#000000",
    "coords": "#000000",
    "chiral": "#0000FF",
    "ez": "#006400",
    "selection": "#FFFF00",
    "measurement": "#FF0000",
    "constraint": "#00FFFF",
}


@pytest.mark.parametrize("kind", [k for k, *_ in LABEL_KINDS])
def test_defaults(kind):
    """Defaults keep each label's previous color; all share one 18 pt size."""
    kwargs = label_kwargs(DEFAULT_SETTINGS, kind)
    assert kwargs["text_color"] == PREVIOUS_COLOR[kind]
    assert kwargs["font_size"] == 18
    assert kwargs["shape_color"] == "#808080"
    assert kwargs["shape_opacity"] == 0.5
    assert kwargs["font_family"] == "arial"
    assert kwargs["bold"] is True
    assert kwargs["italic"] is False


def test_sections_cover_every_kind_once():
    """Every label kind appears in exactly one settings-tab section."""
    listed = [kind for _, kinds in LABEL_SECTIONS for kind in kinds]
    assert sorted(listed) == sorted(k for k, *_ in LABEL_KINDS)
    assert len(listed) == len(set(listed))


def test_every_default_is_in_app_defaults():
    """All label keys are part of DEFAULT_SETTINGS, so Reset restores them."""
    for key, value in DEFAULT_LABEL_SETTINGS.items():
        assert DEFAULT_SETTINGS[key] == value
    assert DEFAULT_SETTINGS["check_chirality_after_conversion"] is True


def test_user_settings_are_applied():
    """Stored values reach the label arguments."""
    settings = {
        color_key("measurement"): "#abcdef",
        "label_background_color_3d": "#010203",
        "label_background_opacity_3d": 0.0,
        "label_font_size_3d": 24,
        "label_font_family_3d": "times",
        "label_font_bold_3d": False,
        "label_font_italic_3d": True,
    }
    assert label_kwargs(settings, "measurement") == {
        "text_color": "#abcdef",
        "shape_color": "#010203",
        "shape_opacity": 0.0,
        "font_size": 24,
        "font_family": "times",
        "bold": False,
        "italic": True,
    }


@pytest.mark.parametrize(
    "settings",
    [
        None,
        MagicMock(),
        {color_key("index"): "", "label_background_color_3d": "red"},
        {color_key("index"): 12, "label_font_family_3d": "comic sans"},
        {"label_background_opacity_3d": "half", "label_font_bold_3d": "yes"},
        {"label_font_size_3d": True},
    ],
    ids=["none", "not-a-mapping", "bad-colors", "wrong-types", "strings", "bool-size"],
)
def test_unusable_values_fall_back_to_defaults(settings):
    """Missing, malformed or non-mapping settings never break a label."""
    assert label_kwargs(settings, "index") == label_kwargs(DEFAULT_SETTINGS, "index")


def test_numbers_are_clamped():
    """Out-of-range opacity and font size are clamped to their ranges."""
    small = label_kwargs(
        {"label_background_opacity_3d": 7, "label_font_size_3d": 1}, "selection"
    )
    assert small["shape_opacity"] == 1.0
    assert small["font_size"] == 6
    big = label_kwargs({"label_font_size_3d": 500}, "chiral")
    assert big["font_size"] == 72


def test_one_font_size_for_every_kind():
    """The font size setting applies to every label kind alike."""
    settings = {"label_font_size_3d": 30}
    assert {label_kwargs(settings, k)["font_size"] for k, *_ in LABEL_KINDS} == {30}
