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
from typing import Any, Optional, Union

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QSlider,
    QWidget,
    QSpinBox,
    QDoubleSpinBox,
)


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
    ) -> tuple[QSlider, Union[QSpinBox, QDoubleSpinBox]]:
        """Create a slider with a linked spinbox showing and setting the value."""
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(min_val, max_val)

        if is_int:
            spin = QSpinBox()
            spin.setRange(min_val, max_val)
            spin.setValue(slider.value())

            def sync_spin(val):
                spin.blockSignals(True)
                spin.setValue(val)
                spin.blockSignals(False)

            def sync_slider(val):
                slider.blockSignals(True)
                slider.setValue(val)
                slider.blockSignals(False)

            slider.valueChanged.connect(sync_spin)
            spin.valueChanged.connect(sync_slider)
        else:
            spin = QDoubleSpinBox()
            spin.setRange(min_val / scale, max_val / scale)
            spin.setSingleStep(1.0 / scale)
            spin.setDecimals(2)
            spin.setValue(slider.value() / scale)

            def sync_spin(val):
                spin.blockSignals(True)
                spin.setValue(val / scale)
                spin.blockSignals(False)

            def sync_slider(val):
                slider.blockSignals(True)
                slider.setValue(int(round(val * scale)))
                slider.blockSignals(False)

            slider.valueChanged.connect(sync_spin)
            spin.valueChanged.connect(sync_slider)

        return slider, spin

    def _wrap_layout(self, slider: QWidget, label: QWidget) -> QHBoxLayout:
        """Wrap a slider and its label in a horizontal layout."""
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(slider)
        layout.addWidget(label)
        return layout
