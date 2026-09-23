#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
from typing import Any, Dict, Optional

from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QApplication,
    QColorDialog,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from ..utils.constants import (
    CPK_COLORS,
    DEFAULT_CPK_COLORS,
    PERIODIC_TABLE_LAYOUT,
)
from .periodic_table_dialog import element_button_style


class ColorSettingsDialog(QDialog):
    """Dialog to customize CPK element colors.

    - Click an element to pick a new color for the element (CPK colors).
    - Reset All button to restore defaults for everything.
    """

    def __init__(self, current_settings: Any, parent: Optional[QWidget] = None) -> None:
        """Initialize CPK element color customization dialog."""
        super().__init__(parent)
        self.setWindowTitle("CPK Colors")
        self.parent_window: Any = parent
        self.current_settings = current_settings or {}

        self.changed_cpk: Dict[str, str] = {}  # symbol -> hex
        self._reset_all_flag = False

        layout = QVBoxLayout(self)

        # Periodic table grid
        grid = QGridLayout()
        self.element_buttons = {}

        for symbol, row, col in PERIODIC_TABLE_LAYOUT:
            b = QPushButton(symbol)
            b.setFixedSize(40, 40)
            # Choose override color (if present) else default CPK color
            override = self.current_settings.get("cpk_colors", {}).get(symbol)
            if override:
                q_color = QColor(override)
            else:
                q_color = CPK_COLORS.get(symbol, CPK_COLORS["DEFAULT"])

            b.setStyleSheet(element_button_style(q_color))
            b.clicked.connect(self.on_element_clicked)
            grid.addWidget(b, row, col)
            self.element_buttons[symbol] = b

        layout.addLayout(grid)

        # Ball & Stick bond color (3D) picker
        self.changed_bs_color: Optional[str] = None
        bs_h = QHBoxLayout()
        bs_label = QLabel("Ball & Stick bond color:")
        self.bs_button = QPushButton()
        self.bs_button.setFixedSize(36, 24)

        # Safe initialization from settings
        settings = self.current_settings or {}
        cur_bs = settings.get("ball_stick_bond_color")

        if not cur_bs and self.parent_window:
            parent_settings = getattr(self.parent_window, "settings", {})
            cur_bs = parent_settings.get("ball_stick_bond_color", "#7F7F7F")

        cur_bs = cur_bs or "#7F7F7F"

        self.bs_button.setStyleSheet(
            f"background-color: {cur_bs}; border: 1px solid #888;"
        )
        self.bs_button.setToolTip(cur_bs)
        self.bs_button.clicked.connect(self.pick_bs_bond_color)

        bs_h.addWidget(bs_label)
        bs_h.addWidget(self.bs_button)
        bs_h.addStretch(1)
        layout.addLayout(bs_h)

        # Action buttons
        h = QHBoxLayout()
        reset_button = QPushButton("Reset All")
        reset_button.clicked.connect(self.reset_all)
        h.addWidget(reset_button)
        h.addStretch(1)
        apply_button = QPushButton("Apply")
        apply_button.clicked.connect(self.apply_changes)
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)

        h.addWidget(apply_button)
        h.addWidget(ok_button)
        h.addWidget(cancel_button)
        layout.addLayout(h)

    def on_element_clicked(self) -> None:
        """Open a color picker for the clicked element button and update its swatch."""
        btn = self.sender()
        symbol = btn.text()  # type: ignore[union-attr]
        cur = self.current_settings.get("cpk_colors", {}).get(symbol)
        if not cur:
            cur = CPK_COLORS.get(symbol, CPK_COLORS["DEFAULT"]).name()
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.changed_cpk[symbol] = color.name()
            btn.setStyleSheet(element_button_style(color))  # type: ignore[union-attr]

    def reset_all(self) -> None:
        """Revert all CPK and ball-and-stick colors to their defaults."""
        self.changed_cpk = {}
        self._reset_all_flag = True

        # Restore CPK button displays to defaults
        for s, btn in self.element_buttons.items():
            q_color = DEFAULT_CPK_COLORS.get(s, DEFAULT_CPK_COLORS["DEFAULT"])
            btn.setStyleSheet(element_button_style(q_color))

        # Restore B&S bond color
        hexv = "#7F7F7F"
        if self.parent_window:
            default_settings = getattr(self.parent_window, "default_settings", {})
            hexv = default_settings.get("ball_stick_bond_color", "#7F7F7F")

        self.changed_bs_color = hexv
        self.bs_button.setStyleSheet(
            f"background-color: {hexv}; border: 1px solid #888;"
        )
        self.bs_button.setToolTip(hexv)

    def apply_changes(self) -> None:
        """Write changed color settings to the application settings and redraw."""
        if not self.parent_window or not hasattr(self.parent_window, "init_manager"):
            return

        settings = self.parent_window.init_manager.settings

        if self._reset_all_flag:
            if "cpk_colors" in settings:
                try:
                    del settings["cpk_colors"]
                except KeyError:
                    # Suppress if cpk_colors key is already missing or removed during reset.
                    # Safe defensive fallback catching KeyError
                    logging.debug("Suppressed non-critical error", exc_info=True)
        if self.changed_cpk:
            cdict = settings.get("cpk_colors", {}).copy()
            cdict.update(self.changed_cpk)
            settings["cpk_colors"] = cdict
            # init_manager owns the flag; save_settings() skips a clean one.
            self.parent_window.init_manager.settings_dirty = True

            # Persist to disk immediately
            self.parent_window.init_manager.save_settings()

        self.parent_window.init_manager.update_cpk_colors_from_settings()

        # Redraw 3D scene
        vm = self.parent_window.view_3d_manager
        vm.apply_3d_settings(redraw=False)
        mol = getattr(vm, "current_mol", None)
        if mol:
            vm.draw_molecule_3d(mol)

        # Update 2D scene
        scene = getattr(self.parent_window.init_manager, "scene", None)
        if scene:
            for it in scene.items():
                if hasattr(it, "update_style"):
                    try:
                        it.update_style()
                    except (AttributeError, RuntimeError, TypeError):
                        logging.debug("Suppressed non-critical error", exc_info=True)
                else:
                    try:
                        it.update()
                    except (AttributeError, RuntimeError, TypeError):
                        logging.debug("Suppressed non-critical error", exc_info=True)

        # Update button styles
        for s, btn in self.element_buttons.items():
            overrides = settings.get("cpk_colors", {})
            q_color = QColor(
                overrides.get(s, CPK_COLORS.get(s, CPK_COLORS["DEFAULT"]).name())
            )
            btn.setStyleSheet(element_button_style(q_color))

        # Refresh SettingsDialog
        from .settings_dialog import SettingsDialog

        for w in QApplication.topLevelWidgets():
            if isinstance(w, SettingsDialog):
                w.update_ui_from_settings(settings)

        # Persist B&S color
        if getattr(self, "changed_bs_color", None):
            settings["ball_stick_bond_color"] = self.changed_bs_color
            self.parent_window.init_manager.settings_dirty = True
            self.parent_window.view_3d_manager.apply_3d_settings()
            mol = getattr(self.parent_window.view_3d_manager, "current_mol", None)
            if mol:
                self.parent_window.view_3d_manager.draw_molecule_3d(mol)
        elif self._reset_all_flag:
            settings["ball_stick_bond_color"] = "#7F7F7F"
            self.parent_window.init_manager.settings_dirty = True
            self.parent_window.init_manager.update_cpk_colors_from_settings()
            self.parent_window.view_3d_manager.apply_3d_settings()
            mol = getattr(self.parent_window.view_3d_manager, "current_mol", None)
            if mol:
                self.parent_window.view_3d_manager.draw_molecule_3d(mol)

    def accept(self) -> None:
        """Apply changes and close the dialog."""
        self.apply_changes()
        super().accept()

    def pick_bs_bond_color(self) -> None:
        """Open a color picker for the ball-and-stick bond color."""
        settings = self.current_settings or {}
        cur = getattr(self, "changed_bs_color", None) or settings.get(
            "ball_stick_bond_color"
        )
        if not cur and self.parent_window:
            parent_settings = getattr(self.parent_window, "settings", {})
            cur = parent_settings.get("ball_stick_bond_color", "#7F7F7F")

        cur = cur or "#7F7F7F"
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            hexv = color.name()
            self.changed_bs_color = hexv
            self.bs_button.setStyleSheet(
                f"background-color: {hexv}; border: 1px solid #888;"
            )
            self.bs_button.setToolTip(hexv)
