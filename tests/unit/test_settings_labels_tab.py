"""Unit tests for the 3D Labels settings tab."""

from unittest.mock import patch

from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QLabel

from moleditpy.ui.settings_tabs.settings_labels_tab import SettingsLabelsTab
from moleditpy.utils.default_settings import DEFAULT_SETTINGS
from moleditpy.utils.label_style import (
    DEFAULT_LABEL_SETTINGS,
    LABEL_KINDS,
    LABEL_SECTIONS,
    color_key,
)

_GET_COLOR = "moleditpy.ui.settings_tabs.settings_labels_tab.QColorDialog.getColor"


def test_defaults_round_trip(app):
    """A fresh tab reports exactly the default label settings."""
    out = SettingsLabelsTab(DEFAULT_SETTINGS).get_settings()
    assert out == dict(
        DEFAULT_LABEL_SETTINGS,
        check_stereo_after_conversion=True,
        check_stereo_after_optimization=False,
    )


def test_one_color_row_per_label_kind(app):
    """Every label kind has its own color swatch, plus the background."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    for kind, *_ in LABEL_KINDS:
        assert color_key(kind) in tab.color_buttons
    assert "label_background_color_3d" in tab.color_buttons


def test_sections_are_headed(app):
    """Label colors carry one subsection per group; every section has a heading."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    headings = {lbl.text() for lbl in tab.findChildren(QLabel)}
    assert "<b>Label Colors</b>" in headings
    for title, _ in LABEL_SECTIONS:
        assert f"<i>{title}</i>" in headings
    assert "<b>Label Appearance</b>" in headings
    assert "<b>Stereo Check (R/S, E/Z)</b>" in headings


def test_update_ui_then_get_settings(app):
    """Custom values survive update_ui -> get_settings."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    custom = dict(DEFAULT_SETTINGS)
    custom.update(
        {
            color_key("chiral"): "#112233",
            "label_background_color_3d": "#445566",
            "label_background_opacity_3d": 0.0,
            "label_font_size_3d": 30,
            "label_font_family_3d": "courier",
            "label_font_bold_3d": False,
            "label_font_italic_3d": True,
            "check_stereo_after_conversion": False,
        }
    )
    tab.update_ui(custom)
    out = tab.get_settings()
    for key in (
        color_key("chiral"),
        "label_background_color_3d",
        "label_background_opacity_3d",
        "label_font_size_3d",
        "label_font_family_3d",
        "label_font_bold_3d",
        "label_font_italic_3d",
        "check_stereo_after_conversion",
    ):
        assert out[key] == custom[key], key
    assert "#112233" in tab.color_buttons[color_key("chiral")].styleSheet()
    assert tab.font_bold_btn.isChecked() is False
    assert tab.font_italic_btn.isChecked() is True


def test_malformed_settings_show_defaults(app):
    """A hand-edited settings file with bad values shows the defaults."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    tab.update_ui({"label_font_family_3d": "wingdings", "label_font_size_3d": "x"})
    out = tab.get_settings()
    assert out["label_font_family_3d"] == "arial"
    assert out["label_font_size_3d"] == 18


def test_pick_color(app):
    """Picking a color from the dialog stores it for that label kind."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    with patch(_GET_COLOR, return_value=QColor("#abcdef")):
        tab._pick_color(color_key("selection"))
    assert tab.get_settings()[color_key("selection")] == "#abcdef"


def test_cancelled_color_pick_changes_nothing(app):
    """Cancelling the color dialog keeps the previous color."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    with patch(_GET_COLOR, return_value=QColor()):
        tab._pick_color(color_key("selection"))
    assert tab.get_settings()[color_key("selection")] == "#FFFF00"


def test_reset_to_defaults(app):
    """Reset Current Tab restores every label setting."""
    tab = SettingsLabelsTab(DEFAULT_SETTINGS)
    tab.update_ui({**DEFAULT_SETTINGS, "label_font_size_3d": 40})
    tab.reset_to_defaults()
    assert tab.get_settings()["label_font_size_3d"] == 18


def test_dialog_has_labels_tab_after_scene(app):
    """The Settings dialog shows 3D Labels right after 3D Scene."""
    from moleditpy.ui.settings_dialog import SettingsDialog

    dlg = SettingsDialog(dict(DEFAULT_SETTINGS))
    try:
        names = [dlg.tab_widget.tabText(i) for i in range(dlg.tab_widget.count())]
        assert names.index("3D Labels") == names.index("3D Scene") + 1
        assert dlg.get_settings()["label_font_family_3d"] == "arial"
    finally:
        dlg.close()
