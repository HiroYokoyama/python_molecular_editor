#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging  # [REPORT ERROR MISSING ATTRIBUTE]
from typing import Any, Optional

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

try:
    from ..utils.constants import CPK_COLORS, DEFAULT_CPK_COLORS
except (AttributeError, RuntimeError, TypeError):
    from moleditpy_linux.utils.constants import CPK_COLORS, DEFAULT_CPK_COLORS


class ColorSettingsDialog(QDialog):
    """Dialog to customize CPK element colors.

    - Click an element to pick a new color for the element (CPK colors).
    - Reset All button to restore defaults for everything.
    """

    def __init__(self, current_settings: Any, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("CPK Colors")
        self.parent_window = parent
        self.current_settings = current_settings or {}

        self.changed_cpk = {}  # symbol -> hex
        self._reset_all_flag = False

        layout = QVBoxLayout(self)

        # Periodic table grid
        grid = QGridLayout()
        self.element_buttons = {}
        elements = [
            ("H", 1, 1),
            ("He", 1, 18),
            ("Li", 2, 1),
            ("Be", 2, 2),
            ("B", 2, 13),
            ("C", 2, 14),
            ("N", 2, 15),
            ("O", 2, 16),
            ("F", 2, 17),
            ("Ne", 2, 18),
            ("Na", 3, 1),
            ("Mg", 3, 2),
            ("Al", 3, 13),
            ("Si", 3, 14),
            ("P", 3, 15),
            ("S", 3, 16),
            ("Cl", 3, 17),
            ("Ar", 3, 18),
            ("K", 4, 1),
            ("Ca", 4, 2),
            ("Sc", 4, 3),
            ("Ti", 4, 4),
            ("V", 4, 5),
            ("Cr", 4, 6),
            ("Mn", 4, 7),
            ("Fe", 4, 8),
            ("Co", 4, 9),
            ("Ni", 4, 10),
            ("Cu", 4, 11),
            ("Zn", 4, 12),
            ("Ga", 4, 13),
            ("Ge", 4, 14),
            ("As", 4, 15),
            ("Se", 4, 16),
            ("Br", 4, 17),
            ("Kr", 4, 18),
            ("Rb", 5, 1),
            ("Sr", 5, 2),
            ("Y", 5, 3),
            ("Zr", 5, 4),
            ("Nb", 5, 5),
            ("Mo", 5, 6),
            ("Tc", 5, 7),
            ("Ru", 5, 8),
            ("Rh", 5, 9),
            ("Pd", 5, 10),
            ("Ag", 5, 11),
            ("Cd", 5, 12),
            ("In", 5, 13),
            ("Sn", 5, 14),
            ("Sb", 5, 15),
            ("Te", 5, 16),
            ("I", 5, 17),
            ("Xe", 5, 18),
            ("Cs", 6, 1),
            ("Ba", 6, 2),
            ("Hf", 6, 4),
            ("Ta", 6, 5),
            ("W", 6, 6),
            ("Re", 6, 7),
            ("Os", 6, 8),
            ("Ir", 6, 9),
            ("Pt", 6, 10),
            ("Au", 6, 11),
            ("Hg", 6, 12),
            ("Tl", 6, 13),
            ("Pb", 6, 14),
            ("Bi", 6, 15),
            ("Po", 6, 16),
            ("At", 6, 17),
            ("Rn", 6, 18),
            ("Fr", 7, 1),
            ("Ra", 7, 2),
            ("Rf", 7, 4),
            ("Db", 7, 5),
            ("Sg", 7, 6),
            ("Bh", 7, 7),
            ("Hs", 7, 8),
            ("Mt", 7, 9),
            ("Ds", 7, 10),
            ("Rg", 7, 11),
            ("Cn", 7, 12),
            ("Nh", 7, 13),
            ("Fl", 7, 14),
            ("Mc", 7, 15),
            ("Lv", 7, 16),
            ("Ts", 7, 17),
            ("Og", 7, 18),
            ("La", 8, 3),
            ("Ce", 8, 4),
            ("Pr", 8, 5),
            ("Nd", 8, 6),
            ("Pm", 8, 7),
            ("Sm", 8, 8),
            ("Eu", 8, 9),
            ("Gd", 8, 10),
            ("Tb", 8, 11),
            ("Dy", 8, 12),
            ("Ho", 8, 13),
            ("Er", 8, 14),
            ("Tm", 8, 15),
            ("Yb", 8, 16),
            ("Lu", 8, 17),
            ("Ac", 9, 3),
            ("Th", 9, 4),
            ("Pa", 9, 5),
            ("U", 9, 6),
            ("Np", 9, 7),
            ("Pu", 9, 8),
            ("Am", 9, 9),
            ("Cm", 9, 10),
            ("Bk", 9, 11),
            ("Cf", 9, 12),
            ("Es", 9, 13),
            ("Fm", 9, 14),
            ("Md", 9, 15),
            ("No", 9, 16),
            ("Lr", 9, 17),
        ]

        for symbol, row, col in elements:
            b = QPushButton(symbol)
            b.setFixedSize(40, 40)
            # Choose override color (if present) else default CPK color
            override = self.current_settings.get("cpk_colors", {}).get(symbol)
            if override:
                q_color = QColor(override)
            else:
                q_color = CPK_COLORS.get(symbol, CPK_COLORS["DEFAULT"])

            brightness = (
                q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114
            ) / 1000
            text_color = "white" if brightness < 128 else "black"
            b.setStyleSheet(
                f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;"
            )
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
        btn = self.sender()
        symbol = btn.text()
        cur = self.current_settings.get("cpk_colors", {}).get(symbol)
        if not cur:
            cur = CPK_COLORS.get(symbol, CPK_COLORS["DEFAULT"]).name()
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.changed_cpk[symbol] = color.name()
            brightness = (
                color.red() * 299 + color.green() * 587 + color.blue() * 114
            ) / 1000
            text_color = "white" if brightness < 128 else "black"
            btn.setStyleSheet(
                f"background-color: {color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;"
            )

    def reset_all(self) -> None:
        self.changed_cpk = {}
        self._reset_all_flag = True

        # Restore CPK button displays to defaults
        for s, btn in self.element_buttons.items():
            q_color = DEFAULT_CPK_COLORS.get(s, DEFAULT_CPK_COLORS["DEFAULT"])
            brightness = (
                q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114
            ) / 1000
            text_color = "white" if brightness < 128 else "black"
            btn.setStyleSheet(
                f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;"
            )

        # Restore B&S bond color
        hexv = "#7F7F7F"
        if self.parent_window:
            default_settings = getattr(self.parent_window, "default_settings", {})
            hexv = default_settings.get("ball_stick_bond_color", "#7F7F7F")

        self.changed_bs_color = hexv
        if hasattr(self, "bs_button"):
            self.bs_button.setStyleSheet(
                f"background-color: {hexv}; border: 1px solid #888;"
            )
            self.bs_button.setToolTip(hexv)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error("REPORT ERROR: Missing attribute 'bs_button' on self")

    def apply_changes(self) -> None:
        if not self.parent_window or not hasattr(self.parent_window, "init_manager"):
            return

        settings = self.parent_window.init_manager.settings

        if self._reset_all_flag:
            if "cpk_colors" in settings:
                try:
                    del settings["cpk_colors"]
                except KeyError:
                    # Suppress if cpk_colors key is already missing or removed during reset.
                    pass
        if self.changed_cpk:
            cdict = settings.get("cpk_colors", {}).copy()
            cdict.update(self.changed_cpk)
            settings["cpk_colors"] = cdict
            self.parent_window.settings_dirty = True

            # Persist to disk immediately
            if hasattr(self.parent_window, "init_manager"):
                self.parent_window.init_manager.save_settings()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'init_manager' on self.parent_window"
                )

        if hasattr(self.parent_window.init_manager, "update_cpk_colors_from_settings"):
            self.parent_window.init_manager.update_cpk_colors_from_settings()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'update_cpk_colors_from_settings' on object"
            )

        # Redraw 3D scene
        if hasattr(self.parent_window, "view_3d_manager"):
            vm = self.parent_window.view_3d_manager
            if hasattr(vm, "apply_3d_settings"):
                vm.apply_3d_settings(redraw=False)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'apply_3d_settings' on vm"
                )

            mol = getattr(vm, "current_mol", None)
            if mol and hasattr(vm, "draw_molecule_3d"):
                vm.draw_molecule_3d(mol)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'view_3d_manager' on self.parent_window"
            )

        # Update 2D scene
        scene = getattr(self.parent_window.init_manager, "scene", None)
        if scene:
            for it in scene.items():
                if hasattr(it, "update_style"):
                    try:
                        it.update_style()
                    except (AttributeError, RuntimeError, TypeError):
                        pass
                elif hasattr(it, "update"):
                    try:
                        it.update()
                    except (AttributeError, RuntimeError, TypeError):
                        pass
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error("REPORT ERROR: Missing attribute 'update' on it")

        # Update button styles
        for s, btn in self.element_buttons.items():
            overrides = settings.get("cpk_colors", {})
            q_color = QColor(
                overrides.get(s, CPK_COLORS.get(s, CPK_COLORS["DEFAULT"]).name())
            )
            brightness = (
                q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114
            ) / 1000
            text_color = "white" if brightness < 128 else "black"
            btn.setStyleSheet(
                f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;"
            )

        # Refresh SettingsDialog
        SettingsDialog: Optional[type] = None
        try:
            from .settings_dialog import SettingsDialog  # type: ignore[assignment]
        except (ImportError, ValueError):
            try:
                from moleditpy_linux.ui.settings_dialog import SettingsDialog  # type: ignore[assignment]
            except ImportError:
                SettingsDialog = None

        if SettingsDialog:
            for w in QApplication.topLevelWidgets():
                if isinstance(w, SettingsDialog):
                    if hasattr(w, "update_ui_from_settings"):
                        w.update_ui_from_settings(settings)
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error(
                            "REPORT ERROR: Missing attribute 'update_ui_from_settings' on w"
                        )

        # Persist B&S color
        if getattr(self, "changed_bs_color", None):
            settings["ball_stick_bond_color"] = self.changed_bs_color
            self.parent_window.init_manager.settings_dirty = True
            if hasattr(self.parent_window.view_3d_manager, "apply_3d_settings"):
                self.parent_window.view_3d_manager.apply_3d_settings()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'apply_3d_settings' on object"
                )

            mol = getattr(self.parent_window.view_3d_manager, "current_mol", None)
            if mol and hasattr(self.parent_window.view_3d_manager, "draw_molecule_3d"):
                self.parent_window.view_3d_manager.draw_molecule_3d(mol)
        elif self._reset_all_flag:
            settings["ball_stick_bond_color"] = "#7F7F7F"
            self.parent_window.init_manager.settings_dirty = True
            if hasattr(
                self.parent_window.init_manager, "update_cpk_colors_from_settings"
            ):
                self.parent_window.init_manager.update_cpk_colors_from_settings()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'update_cpk_colors_from_settings' on object"
                )

            if hasattr(self.parent_window.view_3d_manager, "apply_3d_settings"):
                self.parent_window.view_3d_manager.apply_3d_settings()
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'apply_3d_settings' on object"
                )

            mol = getattr(self.parent_window.view_3d_manager, "current_mol", None)
            if mol and hasattr(self.parent_window.view_3d_manager, "draw_molecule_3d"):
                self.parent_window.view_3d_manager.draw_molecule_3d(mol)

    def accept(self) -> None:
        self.apply_changes()
        super().accept()

    def pick_bs_bond_color(self) -> None:
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
            if hasattr(self, "bs_button"):
                self.bs_button.setStyleSheet(
                    f"background-color: {hexv}; border: 1px solid #888;"
                )
                self.bs_button.setToolTip(hexv)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error("REPORT ERROR: Missing attribute 'bs_button' on self")
