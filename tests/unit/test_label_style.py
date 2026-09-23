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
    size_key,
)

# How each label looked before it could be styled: (text color, font size).
PREVIOUS_LOOK = {
    "index": ("#003366", 18),
    "original_id": ("#009000", 18),
    "xyz_index": ("#8B0000", 18),
    "symbol": ("#000000", 18),
    "coords": ("#000000", 18),
    "chiral": ("#0000FF", 20),
    "ez": ("#006400", 18),
    "selection": ("#FFFF00", 12),
    "measurement": ("#FF0000", 16),
    "constraint": ("#00FFFF", 12),
}


@pytest.mark.parametrize("kind", [k for k, *_ in LABEL_KINDS])
def test_defaults_reproduce_previous_look(kind):
    """With default settings every label looks as it did when hard-coded."""
    kwargs = label_kwargs(DEFAULT_SETTINGS, kind)
    color, size = PREVIOUS_LOOK[kind]
    assert kwargs["text_color"] == color
    assert kwargs["font_size"] == size
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
        size_key("measurement"): 24,
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
        {size_key("index"): True},
    ],
    ids=["none", "not-a-mapping", "bad-colors", "wrong-types", "strings", "bool-size"],
)
def test_unusable_values_fall_back_to_defaults(settings):
    """Missing, malformed or non-mapping settings never break a label."""
    assert label_kwargs(settings, "index") == label_kwargs(DEFAULT_SETTINGS, "index")


def test_numbers_are_clamped():
    """Out-of-range opacity and font size are clamped to their ranges."""
    small = label_kwargs(
        {"label_background_opacity_3d": 7, size_key("selection"): 1}, "selection"
    )
    assert small["shape_opacity"] == 1.0
    assert small["font_size"] == 6
    big = label_kwargs({size_key("chiral"): 500}, "chiral")
    assert big["font_size"] == 72


def test_font_sizes_are_per_kind():
    """Changing one label's size leaves the others alone."""
    settings = {size_key("measurement"): 30}
    assert label_kwargs(settings, "measurement")["font_size"] == 30
    assert label_kwargs(settings, "selection")["font_size"] == 12
