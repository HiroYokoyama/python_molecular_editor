#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Any, Optional

from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QDialog, QGridLayout, QPushButton

from ..utils.constants import CPK_COLORS, PERIODIC_TABLE_LAYOUT

from PyQt6.QtCore import pyqtSignal


def element_button_style(color: QColor) -> str:
    """Stylesheet for an element button: *color* fill, readable label.

    Label is white on dark fills, black on light ones (perceived
    brightness, R*299 + G*587 + B*114).
    """
    brightness = (color.red() * 299 + color.green() * 587 + color.blue() * 114) / 1000
    text_color = "white" if brightness < 128 else "black"
    return (
        f"background-color: {color.name()}; color: {text_color}; "
        "border: 1px solid #555; font-weight: bold;"
    )


class PeriodicTableDialog(QDialog):
    """Dialog with a grid of element buttons for selecting an atom symbol."""

    element_selected = pyqtSignal(str)

    def __init__(self, parent: Optional[Any] = None) -> None:
        """Initialize periodic table element picker dialog."""
        super().__init__(parent)
        self.setWindowTitle("Select an Element")
        layout = QGridLayout(self)
        self.setLayout(layout)

        for symbol, row, col in PERIODIC_TABLE_LAYOUT:
            b = QPushButton(symbol)
            b.setFixedSize(40, 40)

            # Prefer saved user override (from parent.init_manager.settings), otherwise use CPK_COLORS
            try:
                settings = {}
                if parent:
                    if hasattr(parent, "init_manager") and parent.init_manager:
                        settings = parent.init_manager.settings
                    elif hasattr(parent, "settings"):
                        settings = parent.settings
                overrides = settings.get("cpk_colors", {})
                override = overrides.get(symbol)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                override = None
            q_color = (
                QColor(override)
                if override
                else CPK_COLORS.get(symbol, CPK_COLORS["DEFAULT"])
            )

            b.setStyleSheet(element_button_style(q_color))

            b.clicked.connect(self.on_button_clicked)
            layout.addWidget(b, row, col)

    def on_button_clicked(self) -> None:
        """Emit element_selected with the button's symbol and accept the dialog."""
        b = self.sender()
        self.element_selected.emit(b.text())  # type: ignore[union-attr]
        self.accept()
