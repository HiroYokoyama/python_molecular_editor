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
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QPushButton,
    QSlider,
    QComboBox,
    QCheckBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QWidget,
    QColorDialog,
    QFrame,
)
from .settings_tab_base import SettingsTabBase


class Settings3DSceneTab(SettingsTabBase):
    def __init__(self, default_settings, parent=None):
        super().__init__(default_settings, parent)
        self.current_bg_color = default_settings["background_color"]
        self._setup_ui()

    def _setup_ui(self):
        form_layout = QFormLayout(self)

        self.bg_button = QPushButton()
        self.bg_button.setToolTip("Click to select a color")
        self.bg_button.clicked.connect(self._select_color)
        form_layout.addRow("Background Color:", self.bg_button)

        self.axes_checkbox = QCheckBox()
        form_layout.addRow("Show 3D Axes:", self.axes_checkbox)

        self.light_checkbox = QCheckBox()
        form_layout.addRow("Enable Lighting:", self.light_checkbox)

        self.intensity_slider = QSlider(Qt.Orientation.Horizontal)
        self.intensity_slider.setRange(0, 200)
        self.intensity_label = QLabel("1.0")
        self.intensity_slider.valueChanged.connect(
            lambda v: self.intensity_label.setText(f"{v / 100:.2f}")
        )

        il = QHBoxLayout()
        il.addWidget(self.intensity_slider)
        il.addWidget(self.intensity_label)
        form_layout.addRow("Light Intensity:", il)

        self.specular_slider = QSlider(Qt.Orientation.Horizontal)
        self.specular_slider.setRange(0, 100)
        self.specular_label = QLabel("0.20")
        self.specular_slider.valueChanged.connect(
            lambda v: self.specular_label.setText(f"{v / 100:.2f}")
        )

        sl = QHBoxLayout()
        sl.addWidget(self.specular_slider)
        sl.addWidget(self.specular_label)
        form_layout.addRow("Shininess (Specular):", sl)

        self.spec_power_slider = QSlider(Qt.Orientation.Horizontal)
        self.spec_power_slider.setRange(0, 100)
        self.spec_power_label = QLabel("20")
        self.spec_power_slider.valueChanged.connect(
            lambda v: self.spec_power_label.setText(str(v))
        )

        spl = QHBoxLayout()
        spl.addWidget(self.spec_power_slider)
        spl.addWidget(self.spec_power_label)
        form_layout.addRow("Shininess Power:", spl)

        self.projection_combo = QComboBox()
        self.projection_combo.addItems(["Perspective", "Orthographic"])
        form_layout.addRow("Projection Mode:", self.projection_combo)

    def _select_color(self):
        color = QColorDialog.getColor(QColor(self.current_bg_color), self)
        if color.isValid():
            self.current_bg_color = color.name()
            self._update_color_button()

    def _update_color_button(self):
        self.bg_button.setStyleSheet(
            f"background-color: {self.current_bg_color}; border: 1px solid #888;"
        )

    def update_ui(self, settings_dict):
        self.current_bg_color = settings_dict.get(
            "background_color", self.default_settings["background_color"]
        )
        self._update_color_button()
        self.axes_checkbox.setChecked(settings_dict.get("show_3d_axes", True))
        self.light_checkbox.setChecked(settings_dict.get("lighting_enabled", True))

        int_val = settings_dict.get("light_intensity", 1.0)
        self.intensity_slider.setValue(int(int_val * 100))

        spec_val = settings_dict.get("specular", 0.2)
        self.specular_slider.setValue(int(spec_val * 100))

        pow_val = settings_dict.get("specular_power", 20)
        self.spec_power_slider.setValue(int(pow_val))

        proj = settings_dict.get("projection_mode", "Perspective")
        idx = self.projection_combo.findText(proj)
        if idx >= 0:
            self.projection_combo.setCurrentIndex(idx)

    def get_settings(self):
        return {
            "background_color": self.current_bg_color,
            "show_3d_axes": self.axes_checkbox.isChecked(),
            "lighting_enabled": self.light_checkbox.isChecked(),
            "light_intensity": self.intensity_slider.value() / 100.0,
            "specular": self.specular_slider.value() / 100.0,
            "specular_power": self.spec_power_slider.value(),
            "projection_mode": self.projection_combo.currentText(),
        }


class SettingsModelTab(SettingsTabBase):
    def __init__(self, model_prefix, info_text, default_settings, parent=None):
        self.prefix = model_prefix
        self.info_text = info_text
        super().__init__(default_settings, parent)
        self._setup_ui()

    def _setup_ui(self):
        form_layout = QFormLayout(self)

        info_label = QLabel(self.info_text)
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        form_layout.addRow(info_label)

        # Common Atom Scale (for BS and CPK)
        if self.prefix in ["ball_stick", "cpk"]:
            self.atom_scale_slider, self.atom_scale_label = self._create_slider(
                50, 200, 100
            )
            form_layout.addRow(
                "Atom Size Scale:",
                self._wrap(self.atom_scale_slider, self.atom_scale_label),
            )

        # Common Bond Radius (for BS, WF, Stick)
        if self.prefix in ["ball_stick", "wireframe", "stick"]:
            ranges = {"ball_stick": (1, 50), "wireframe": (1, 10), "stick": (5, 50)}
            r_min, r_max = ranges[self.prefix]
            self.bond_radius_slider, self.bond_radius_label = self._create_slider(
                r_min, r_max, 100
            )
            form_layout.addRow(
                "Bond Radius:",
                self._wrap(self.bond_radius_slider, self.bond_radius_label),
            )

        # Model-specific multi-bond offsets
        if self.prefix in ["ball_stick", "wireframe", "stick"]:
            form_layout.addRow(self._create_separator())
            self.db_offset_slider, self.db_offset_label = self._create_slider(
                50, 400, 100
            )
            form_layout.addRow(
                "Double Bond Offset:",
                self._wrap(self.db_offset_slider, self.db_offset_label),
            )

            self.tr_offset_slider, self.tr_offset_label = self._create_slider(
                50, 400, 100
            )
            form_layout.addRow(
                "Triple Bond Offset:",
                self._wrap(self.tr_offset_slider, self.tr_offset_label),
            )

            self.db_radius_slider, self.db_radius_label = self._create_slider(
                20, 100, 100
            )
            form_layout.addRow(
                "Double Bond Thickness:",
                self._wrap(self.db_radius_slider, self.db_radius_label),
            )

            self.tr_radius_slider, self.tr_radius_label = self._create_slider(
                20, 100, 100
            )
            form_layout.addRow(
                "Triple Bond Thickness:",
                self._wrap(self.tr_radius_slider, self.tr_radius_label),
            )

        # Common Resolution
        form_layout.addRow(self._create_separator())
        res_ranges = {
            "ball_stick": (6, 32),
            "cpk": (8, 64),
            "wireframe": (4, 16),
            "stick": (6, 32),
        }
        r_min, r_max = res_ranges[self.prefix]
        self.res_slider, self.res_label = self._create_slider(
            r_min, r_max, 1, is_int=True
        )
        form_layout.addRow(
            "Resolution (Quality):", self._wrap(self.res_slider, self.res_label)
        )

        # Ball & Stick specific color options
        if self.prefix == "ball_stick":
            form_layout.addRow(self._create_separator())
            self.bond_color_button = QPushButton()
            self.bond_color_button.setFixedSize(36, 24)
            self.bond_color_button.clicked.connect(self._pick_bond_color)
            form_layout.addRow("Bond Color:", self.bond_color_button)

            self.use_cpk_checkbox = QCheckBox("Use CPK colors for bonds")
            form_layout.addRow(self.use_cpk_checkbox)

    def _create_slider(self, min_val, max_val, scale, is_int=False):
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_val, max_val)
        label = QLabel()
        if is_int:
            slider.valueChanged.connect(lambda v: label.setText(str(v)))
        else:
            slider.valueChanged.connect(lambda v: label.setText(f"{v / scale:.2f}"))
        return slider, label

    def _wrap(self, slider, label):
        h = QHBoxLayout()
        h.addWidget(slider)
        h.addWidget(label)
        w = QWidget()
        w.setLayout(h)
        return w

    def _create_separator(self):
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line

    def _pick_bond_color(self):
        cur = getattr(self, "current_bond_color", "#7F7F7F")
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.current_bond_color = color.name()
            self.bond_color_button.setStyleSheet(
                f"background-color: {self.current_bond_color}; border: 1px solid #888;"
            )

    def update_ui(self, settings_dict):
        p = self.prefix
        if p in ["ball_stick", "cpk"]:
            val = settings_dict.get(f"{p}_atom_scale", 1.0)
            self.atom_scale_slider.setValue(int(val * 100))
        if p in ["ball_stick", "wireframe", "stick"]:
            prefix = "wf" if p == "wireframe" else p
            val = settings_dict.get(f"{prefix}_bond_radius", 0.1)
            self.bond_radius_slider.setValue(int(val * 100))

            self.db_offset_slider.setValue(
                int(settings_dict.get(f"{p}_double_bond_offset_factor", 2.0) * 100)
            )
            self.tr_offset_slider.setValue(
                int(settings_dict.get(f"{p}_triple_bond_offset_factor", 2.0) * 100)
            )
            self.db_radius_slider.setValue(
                int(settings_dict.get(f"{p}_double_bond_radius_factor", 0.8) * 100)
            )
            self.tr_radius_slider.setValue(
                int(settings_dict.get(f"{p}_triple_bond_radius_factor", 0.75) * 100)
            )

        self.res_slider.setValue(int(settings_dict.get(f"{p}_resolution", 16)))

        if p == "ball_stick":
            self.current_bond_color = settings_dict.get(
                "ball_stick_bond_color", "#7F7F7F"
            )
            self.bond_color_button.setStyleSheet(
                f"background-color: {self.current_bond_color}; border: 1px solid #888;"
            )
            self.use_cpk_checkbox.setChecked(
                settings_dict.get("ball_stick_use_cpk_bond_color", False)
            )

    def get_settings(self):
        s = {}
        p = self.prefix
        if p in ["ball_stick", "cpk"]:
            s[f"{p}_atom_scale"] = self.atom_scale_slider.value() / 100.0
        if p in ["ball_stick", "wireframe", "stick"]:
            prefix = "wireframe" if p == "wireframe" else p
            s[f"{prefix}_bond_radius"] = self.bond_radius_slider.value() / 100.0
            s[f"{p}_double_bond_offset_factor"] = self.db_offset_slider.value() / 100.0
            s[f"{p}_triple_bond_offset_factor"] = self.tr_offset_slider.value() / 100.0
            s[f"{p}_double_bond_radius_factor"] = self.db_radius_slider.value() / 100.0
            s[f"{p}_triple_bond_radius_factor"] = self.tr_radius_slider.value() / 100.0

        s[f"{p}_resolution"] = self.res_slider.value()

        if p == "ball_stick":
            s["ball_stick_bond_color"] = getattr(self, "current_bond_color", "#7F7F7F")
            s["ball_stick_use_cpk_bond_color"] = self.use_cpk_checkbox.isChecked()
        return s
