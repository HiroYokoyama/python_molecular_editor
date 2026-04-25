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
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QSlider, QWidget


class SettingsTabBase(QWidget):
    """Base class for all settings tabs."""

    def __init__(
        self, default_settings: Mapping[str, Any], parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.default_settings = default_settings

    def update_ui(self, settings_dict: Mapping[str, Any]) -> None:
        """Update UI based on settings dictionary. Must be implemented by subclass."""
        raise NotImplementedError

    def get_settings(self) -> dict[str, Any]:
        """Get settings values from the current UI. Must be implemented by subclass."""
        raise NotImplementedError

    def reset_to_defaults(self) -> None:
        """Reset only the settings of the current tab to defaults."""
        self.update_ui(self.default_settings)

    def _create_separator(self) -> QFrame:
        """Create a horizontal separator line."""
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line

    def _create_slider(
        self, min_val: int, max_val: int, scale: float = 1.0, is_int: bool = False
    ) -> tuple[QSlider, QLabel]:
        """Create a slider with a linked label showing the value."""
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_val, max_val)
        label = QLabel()
        if is_int:
            slider.valueChanged.connect(lambda v: label.setText(str(v)))
            label.setText(str(slider.value()))
        else:
            slider.valueChanged.connect(lambda v: label.setText(f"{v / scale:.2f}"))
            label.setText(f"{slider.value() / scale:.2f}")
        return slider, label

    def _wrap_layout(self, slider: QWidget, label: QWidget) -> QWidget:
        """Wrap a slider and its label in a horizontal layout/widget."""
        layout = QHBoxLayout()
        layout.addWidget(slider)
        layout.addWidget(label)
        container = QWidget()
        container.setLayout(layout)
        return container
