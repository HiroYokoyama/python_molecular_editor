#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QFont
from PyQt6.QtWidgets import (
    QPushButton,
    QSlider,
    QComboBox,
    QFontComboBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QWidget,
    QColorDialog,
)
from .settings_tab_base import SettingsTabBase


class Settings2DTab(SettingsTabBase):
    def __init__(self, default_settings, parent=None):
        super().__init__(default_settings, parent)
        self.current_bg_color_2d = default_settings["background_color_2d"]
        self.current_bond_color_2d = default_settings["bond_color_2d"]
        self._setup_ui()

    def _setup_ui(self):
        form_layout = QFormLayout(self)

        # --- View Appearance ---
        form_layout.addRow(QLabel("<b>View Appearance</b>"))

        self.bg_color_2d_button = QPushButton()
        self.bg_color_2d_button.setFixedSize(60, 24)
        self.bg_color_2d_button.clicked.connect(self._pick_bg_color_2d)
        form_layout.addRow("Background Color:", self.bg_color_2d_button)

        form_layout.addRow(self._create_separator())

        # --- Bond Settings ---
        form_layout.addRow(QLabel("<b>Bond Settings</b>"))

        self.bond_color_2d_button = QPushButton()
        self.bond_color_2d_button.setFixedSize(60, 24)
        self.bond_color_2d_button.clicked.connect(self._pick_bond_color_2d)
        form_layout.addRow("Bond Color:", self.bond_color_2d_button)

        # Bond Width
        self.bond_width_2d_slider, self.bond_width_2d_label = (
            self._create_slider_with_label(10, 200, 10.0)
        )
        form_layout.addRow(
            "Bond Width:",
            self._wrap_layout(self.bond_width_2d_slider, self.bond_width_2d_label),
        )

        # Double Bond Spacing
        self.bond_spacing_double_2d_slider, self.bond_spacing_double_2d_label = (
            self._create_slider_with_label(10, 200, 10.0)
        )
        form_layout.addRow(
            "Double Bond Spacing:",
            self._wrap_layout(
                self.bond_spacing_double_2d_slider, self.bond_spacing_double_2d_label
            ),
        )

        # Triple Bond Spacing
        self.bond_spacing_triple_2d_slider, self.bond_spacing_triple_2d_label = (
            self._create_slider_with_label(10, 200, 10.0)
        )
        form_layout.addRow(
            "Triple Bond Spacing:",
            self._wrap_layout(
                self.bond_spacing_triple_2d_slider, self.bond_spacing_triple_2d_label
            ),
        )

        self.bond_cap_style_2d_combo = QComboBox()
        self.bond_cap_style_2d_combo.addItems(["Round", "Flat", "Square"])
        form_layout.addRow("Bond Cap Style:", self.bond_cap_style_2d_combo)

        # Wedge Bond Width
        self.bond_wedge_width_2d_slider, self.bond_wedge_width_2d_label = (
            self._create_slider_with_label(10, 300, 10.0)
        )
        form_layout.addRow(
            "Wedge Bond Width:",
            self._wrap_layout(
                self.bond_wedge_width_2d_slider, self.bond_wedge_width_2d_label
            ),
        )

        # Dash Count
        self.bond_dash_count_2d_slider, self.bond_dash_count_2d_label = (
            self._create_slider_with_label(3, 20, 1.0, is_int=True)
        )
        form_layout.addRow(
            "Dash Count:",
            self._wrap_layout(
                self.bond_dash_count_2d_slider, self.bond_dash_count_2d_label
            ),
        )

        form_layout.addRow(self._create_separator())

        # --- Atom Settings ---
        form_layout.addRow(QLabel("<b>Atom Settings</b>"))

        self.atom_font_family_2d_combo = QFontComboBox()
        self.atom_font_family_2d_combo.setEditable(False)
        self.atom_font_family_2d_combo.setFontFilters(
            QFontComboBox.FontFilter.ScalableFonts
        )
        form_layout.addRow("Atom Label Font Family:", self.atom_font_family_2d_combo)

        self.atom_font_size_2d_slider, self.atom_font_size_2d_label = (
            self._create_slider_with_label(8, 72, 1.0, is_int=True)
        )
        form_layout.addRow(
            "Atom Label Font Size:",
            self._wrap_layout(
                self.atom_font_size_2d_slider, self.atom_font_size_2d_label
            ),
        )

        from PyQt6.QtWidgets import QCheckBox

        self.atom_use_bond_color_2d_checkbox = QCheckBox()
        self.atom_use_bond_color_2d_checkbox.setToolTip(
            "If checked, atoms will use the unified Bond Color instead of element-specific colors (CPK)."
        )
        form_layout.addRow(
            "Use Bond Color for Atoms:", self.atom_use_bond_color_2d_checkbox
        )

    def _create_separator(self):
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line

    def _create_slider_with_label(self, min_val, max_val, scale=1.0, is_int=False):
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_val, max_val)
        label = QLabel()
        if is_int:
            slider.valueChanged.connect(lambda v: label.setText(str(v)))
        else:
            slider.valueChanged.connect(lambda v: label.setText(f"{v / scale:.1f}"))
        return slider, label

    def _wrap_layout(self, slider, label):
        layout = QHBoxLayout()
        layout.addWidget(slider)
        layout.addWidget(label)
        container = QWidget()
        container.setLayout(layout)
        return container

    def _pick_bg_color_2d(self):
        color = QColorDialog.getColor(
            QColor(self.current_bg_color_2d), self, "Select 2D Background Color"
        )
        if color.isValid():
            self.current_bg_color_2d = color.name()
            self._update_color_buttons()

    def _pick_bond_color_2d(self):
        color = QColorDialog.getColor(
            QColor(self.current_bond_color_2d), self, "Select 2D Bond Color"
        )
        if color.isValid():
            self.current_bond_color_2d = color.name()
            self._update_color_buttons()

    def _update_color_buttons(self):
        self.bg_color_2d_button.setStyleSheet(
            f"background-color: {self.current_bg_color_2d}; border: 1px solid #888;"
        )
        self.bond_color_2d_button.setStyleSheet(
            f"background-color: {self.current_bond_color_2d}; border: 1px solid #888;"
        )

    def update_ui(self, settings_dict):
        self.current_bg_color_2d = settings_dict.get(
            "background_color_2d", self.default_settings["background_color_2d"]
        )
        self.current_bond_color_2d = settings_dict.get(
            "bond_color_2d", self.default_settings["bond_color_2d"]
        )
        self._update_color_buttons()

        self.bond_width_2d_slider.setValue(
            int(settings_dict.get("bond_width_2d", 2.0) * 10)
        )
        self.bond_spacing_double_2d_slider.setValue(
            int(settings_dict.get("bond_spacing_double_2d", 3.5) * 10)
        )
        self.bond_spacing_triple_2d_slider.setValue(
            int(settings_dict.get("bond_spacing_triple_2d", 3.5) * 10)
        )

        cap_style = settings_dict.get("bond_cap_style_2d", "Round")
        idx = self.bond_cap_style_2d_combo.findText(cap_style)
        if idx >= 0:
            self.bond_cap_style_2d_combo.setCurrentIndex(idx)

        self.bond_wedge_width_2d_slider.setValue(
            int(settings_dict.get("bond_wedge_width_2d", 6.0) * 10)
        )
        self.bond_dash_count_2d_slider.setValue(
            settings_dict.get("bond_dash_count_2d", 8)
        )

        font_family = settings_dict.get("atom_font_family_2d", "Arial")
        self.atom_font_family_2d_combo.setCurrentFont(QFont(font_family))

        self.atom_font_size_2d_slider.setValue(
            settings_dict.get("atom_font_size_2d", 20)
        )
        self.atom_use_bond_color_2d_checkbox.setChecked(
            settings_dict.get("atom_use_bond_color_2d", False)
        )

    def get_settings(self):
        return {
            "background_color_2d": self.current_bg_color_2d,
            "bond_color_2d": self.current_bond_color_2d,
            "bond_width_2d": self.bond_width_2d_slider.value() / 10.0,
            "bond_spacing_double_2d": self.bond_spacing_double_2d_slider.value() / 10.0,
            "bond_spacing_triple_2d": self.bond_spacing_triple_2d_slider.value() / 10.0,
            "bond_cap_style_2d": self.bond_cap_style_2d_combo.currentText(),
            "bond_wedge_width_2d": self.bond_wedge_width_2d_slider.value() / 10.0,
            "bond_dash_count_2d": self.bond_dash_count_2d_slider.value(),
            "atom_font_family_2d": self.atom_font_family_2d_combo.currentFont().family(),
            "atom_font_size_2d": self.atom_font_size_2d_slider.value(),
            "atom_use_bond_color_2d": self.atom_use_bond_color_2d_checkbox.isChecked(),
        }
