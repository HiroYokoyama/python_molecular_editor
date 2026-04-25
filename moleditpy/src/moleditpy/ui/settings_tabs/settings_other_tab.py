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
from typing import Any, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QSlider,
    QWidget,
)
from .settings_tab_base import SettingsTabBase


class SettingsOtherTab(SettingsTabBase):
    def __init__(
        self, default_settings: Mapping[str, Any], parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(default_settings, parent)
        self._setup_ui()

    def _setup_ui(self) -> None:
        form_layout = QFormLayout(self)

        self.skip_chem_checks_checkbox = QCheckBox()
        self.skip_chem_checks_checkbox.setToolTip(
            "When enabled, XYZ file import will try to ignore chemical/sanitization errors."
        )
        form_layout.addRow(
            "Skip chemistry checks on import XYZ file:", self.skip_chem_checks_checkbox
        )

        self.always_ask_charge_checkbox = QCheckBox()
        self.always_ask_charge_checkbox.setToolTip(
            "Prompt for overall molecular charge when importing XYZ files."
        )
        form_layout.addRow(
            "Always ask molecular charge on import XYZ file:",
            self.always_ask_charge_checkbox,
        )

        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        self.kekule_3d_checkbox = QCheckBox()
        self.kekule_3d_checkbox.setToolTip(
            "Enable alternating single/double bonds in 3D view for aromatic rings."
        )
        self.kekule_3d_checkbox.toggled.connect(self._on_kekule_toggled)
        form_layout.addRow("Display Kekulé bonds in 3D:", self.kekule_3d_checkbox)

        self.aromatic_circle_checkbox = QCheckBox()
        self.aromatic_circle_checkbox.setToolTip(
            "Display a circle inside aromatic rings in 3D view."
        )
        self.aromatic_circle_checkbox.toggled.connect(self._on_aromatic_toggled)
        form_layout.addRow(
            "Display aromatic rings as circles in 3D:", self.aromatic_circle_checkbox
        )

        self.aromatic_torus_thickness_slider = QSlider(Qt.Orientation.Horizontal)
        self.aromatic_torus_thickness_slider.setRange(10, 300)
        self.aromatic_torus_thickness_label = QLabel("0.6")
        self.aromatic_torus_thickness_slider.valueChanged.connect(
            lambda v: self.aromatic_torus_thickness_label.setText(f"{v / 100:.1f}")
        )

        atl = QHBoxLayout()
        atl.addWidget(self.aromatic_torus_thickness_slider)
        atl.addWidget(self.aromatic_torus_thickness_label)
        form_layout.addRow("Aromatic torus thickness (× bond radius):", atl)

    def _on_kekule_toggled(self, checked: bool) -> None:
        self.aromatic_circle_checkbox.setEnabled(not checked)

    def _on_aromatic_toggled(self, checked: bool) -> None:
        self.kekule_3d_checkbox.setEnabled(not checked)

    def update_ui(self, settings_dict: Mapping[str, Any]) -> None:
        self.skip_chem_checks_checkbox.setChecked(
            settings_dict.get("skip_chemistry_checks", False)
        )
        self.always_ask_charge_checkbox.setChecked(
            settings_dict.get("always_ask_charge", False)
        )

        display_kekule = settings_dict.get("display_kekule_3d", False)
        display_aromatic = settings_dict.get("display_aromatic_circles_3d", False)

        self.kekule_3d_checkbox.setChecked(display_kekule)
        self.aromatic_circle_checkbox.setChecked(display_aromatic)

        # Set initial enabled state
        self.aromatic_circle_checkbox.setEnabled(not display_kekule)
        self.kekule_3d_checkbox.setEnabled(not display_aromatic)

        thick = settings_dict.get("aromatic_torus_thickness_factor", 0.6)
        self.aromatic_torus_thickness_slider.setValue(int(thick * 100))

    def get_settings(self) -> dict[str, Any]:
        return {
            "skip_chemistry_checks": self.skip_chem_checks_checkbox.isChecked(),
            "always_ask_charge": self.always_ask_charge_checkbox.isChecked(),
            "display_kekule_3d": self.kekule_3d_checkbox.isChecked(),
            "display_aromatic_circles_3d": self.aromatic_circle_checkbox.isChecked(),
            "aromatic_torus_thickness_factor": self.aromatic_torus_thickness_slider.value()
            / 100.0,
        }
