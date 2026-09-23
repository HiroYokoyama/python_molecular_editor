#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from collections.abc import Mapping
from typing import Any, Dict, Optional

from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QWidget,
)

from ...utils.label_style import (
    BACKGROUND_COLOR_KEY,
    BACKGROUND_OPACITY_KEY,
    FONT_BOLD_KEY,
    FONT_FAMILIES,
    FONT_FAMILY_KEY,
    FONT_ITALIC_KEY,
    FONT_SIZE_RANGE,
    LABEL_KINDS,
    LABEL_SECTIONS,
    color_key,
    label_kwargs,
    size_key,
)
from .settings_tab_base import SettingsTabBase


class SettingsLabelsTab(SettingsTabBase):
    """Settings tab for 3D label colors, background and size."""

    def __init__(
        self, default_settings: Mapping[str, Any], parent: Optional[QWidget] = None
    ) -> None:
        """Initialize the 3D labels settings tab."""
        super().__init__(default_settings, parent)
        self.colors: Dict[str, str] = {}
        self.color_buttons: Dict[str, QPushButton] = {}
        self.size_spins: Dict[str, QSpinBox] = {}
        self.chirality_check_checkbox: Any = None
        self._setup_ui()
        self.update_ui(default_settings)

    def _setup_ui(self) -> None:
        """Construct the color rows, appearance sliders and chirality option."""
        form_layout = self._create_form_layout()

        form_layout.addRow(QLabel("<b>Label Colors and Sizes</b>"))
        names = {kind: name for kind, name, _, _ in LABEL_KINDS}
        for title, kinds in LABEL_SECTIONS:
            form_layout.addRow(QLabel(f"<i>{title}</i>"))
            for kind in kinds:
                form_layout.addRow(f"{names[kind]}:", self._make_label_row(kind))
        form_layout.addRow(self._create_separator())

        form_layout.addRow(QLabel("<b>Label Appearance</b>"))
        form_layout.addRow(
            "Background Color:", self._make_color_button(BACKGROUND_COLOR_KEY)
        )

        self.opacity_slider, self.opacity_label = self._create_slider(0, 100, 100.0)
        self.opacity_slider.setToolTip("0 makes the label background transparent")
        form_layout.addRow(
            "Background Opacity:",
            self._wrap_layout(self.opacity_slider, self.opacity_label),
        )

        self.font_family_combo = QComboBox()
        for value, name in FONT_FAMILIES:
            self.font_family_combo.addItem(name, value)
        form_layout.addRow("Label Font Family:", self.font_family_combo)

        # B / I toggles, as for the 2D atom label font (VTK has no underline).
        self.font_bold_btn = self._make_style_button("B", bold=True)
        self.font_italic_btn = self._make_style_button("I", italic=True)
        style_row = QHBoxLayout()
        style_row.setContentsMargins(0, 0, 0, 0)
        style_row.setSpacing(4)
        style_row.addWidget(self.font_bold_btn)
        style_row.addWidget(self.font_italic_btn)
        style_row.addStretch()
        style_widget = QWidget()
        style_widget.setLayout(style_row)
        form_layout.addRow("Label Font Style:", style_widget)

        form_layout.addRow(self._create_separator())
        form_layout.addRow(QLabel("<b>Chirality Check</b>"))
        self.chirality_check_checkbox = QCheckBox()
        self.chirality_check_checkbox.setToolTip(
            "After Convert 2D to 3D, compare every stereocenter drawn with a\n"
            "wedge or hash against the 3D result, and warn if any differs."
        )
        form_layout.addRow(
            "Check Chirality After 3D Conversion:", self.chirality_check_checkbox
        )

    def _make_label_row(self, kind: str) -> QWidget:
        """Color swatch and font size box for one label kind."""
        size = QSpinBox()
        size.setRange(*FONT_SIZE_RANGE)
        size.setSuffix(" pt")
        size.setToolTip("Font size")
        self.size_spins[size_key(kind)] = size

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._make_color_button(color_key(kind)))
        row.addWidget(size)
        row.addStretch()
        widget = QWidget()
        widget.setLayout(row)
        return widget

    @staticmethod
    def _make_style_button(
        text: str, bold: bool = False, italic: bool = False
    ) -> QPushButton:
        """A checkable 28x24 font-style toggle drawn in its own style."""
        button = QPushButton(text)
        button.setCheckable(True)
        button.setFixedSize(28, 24)
        font = QFont()
        font.setBold(bold)
        font.setItalic(italic)
        button.setFont(font)
        return button

    def _make_color_button(self, key: str) -> QPushButton:
        """Create a swatch button that picks the color stored under *key*."""
        button = QPushButton()
        button.setFixedSize(60, 24)
        button.setToolTip("Click to select a color")
        button.clicked.connect(lambda: self._pick_color(key))
        self.color_buttons[key] = button
        return button

    def _pick_color(self, key: str) -> None:
        """Open a color dialog for *key* and store the choice."""
        color = QColorDialog.getColor(QColor(self.colors[key]), self)
        if color.isValid():
            self._set_color(key, color.name())

    def _set_color(self, key: str, value: str) -> None:
        """Store a color and show it on its swatch."""
        self.colors[key] = value
        self.color_buttons[key].setStyleSheet(
            f"background-color: {value}; border: 1px solid #888;"
        )

    def update_ui(self, settings_dict: Mapping[str, Any]) -> None:
        """Update every control from the settings dictionary."""
        for kind, _, _, _ in LABEL_KINDS:
            style = label_kwargs(settings_dict, kind)
            self._set_color(color_key(kind), style["text_color"])
            self.size_spins[size_key(kind)].setValue(style["font_size"])
        # Background, opacity and font are shared, so any kind reads them.
        shared = label_kwargs(settings_dict, LABEL_KINDS[0][0])
        self._set_color(BACKGROUND_COLOR_KEY, shared["shape_color"])
        self.opacity_slider.setValue(int(round(shared["shape_opacity"] * 100)))
        index = self.font_family_combo.findData(shared["font_family"])
        self.font_family_combo.setCurrentIndex(max(index, 0))
        self.font_bold_btn.setChecked(shared["bold"])
        self.font_italic_btn.setChecked(shared["italic"])
        self.chirality_check_checkbox.setChecked(
            settings_dict.get("check_chirality_after_conversion", True)
        )

    def get_settings(self) -> dict[str, Any]:
        """Collect the label settings into a dictionary."""
        return {
            **self.colors,
            **{key: spin.value() for key, spin in self.size_spins.items()},
            BACKGROUND_OPACITY_KEY: self.opacity_slider.value() / 100.0,
            FONT_FAMILY_KEY: self.font_family_combo.currentData(),
            FONT_BOLD_KEY: self.font_bold_btn.isChecked(),
            FONT_ITALIC_KEY: self.font_italic_btn.isChecked(),
            "check_chirality_after_conversion": self.chirality_check_checkbox.isChecked(),
        }
