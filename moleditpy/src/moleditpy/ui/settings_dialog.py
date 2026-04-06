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

from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QMessageBox,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
)

from .settings_tabs.settings_2d_tab import Settings2DTab
from .settings_tabs.settings_3d_tabs import Settings3DSceneTab, SettingsModelTab
from .settings_tabs.settings_other_tab import SettingsOtherTab

try:
    from ..utils.default_settings import DEFAULT_SETTINGS
except ImportError:
    from moleditpy.utils.default_settings import DEFAULT_SETTINGS


class SettingsDialog(QDialog):
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumSize(650, 750)
        self.parent_window = parent

        self.default_settings = DEFAULT_SETTINGS.copy()

        self._setup_ui(current_settings)

    def _setup_ui(self, current_settings):
        layout = QVBoxLayout(self)
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)

        # Instantiate Tabs
        self.tab_2d = Settings2DTab(self.default_settings, self)
        self.tab_scene = Settings3DSceneTab(self.default_settings, self)
        self.tab_bs = SettingsModelTab(
            "ball_stick",
            "Ball & Stick model shows atoms as spheres and bonds as cylinders.",
            self.default_settings,
            self,
        )
        self.tab_cpk = SettingsModelTab(
            "cpk",
            "CPK model shows atoms as space-filling spheres using van der Waals radii.",
            self.default_settings,
            self,
        )
        self.tab_wf = SettingsModelTab(
            "wireframe",
            "Wireframe model shows molecular structure with thin lines only.",
            self.default_settings,
            self,
        )
        self.tab_stick = SettingsModelTab(
            "stick",
            "Stick model shows bonds as thick cylinders with atoms as small spheres.",
            self.default_settings,
            self,
        )
        self.tab_other = SettingsOtherTab(self.default_settings, self)

        # Add Tabs to Widget
        self.tab_widget.addTab(self.tab_2d, "2D Settings")
        self.tab_widget.addTab(self.tab_scene, "3D Scene")
        self.tab_widget.addTab(self.tab_bs, "Ball & Stick")
        self.tab_widget.addTab(self.tab_cpk, "CPK (Space-filling)")
        self.tab_widget.addTab(self.tab_wf, "Wireframe")
        self.tab_widget.addTab(self.tab_stick, "Stick")
        self.tab_widget.addTab(self.tab_other, "Other")

        # Command Buttons
        buttons_layout = QHBoxLayout()

        reset_tab_btn = QPushButton("Reset Current Tab")
        reset_tab_btn.clicked.connect(self.reset_current_tab)
        buttons_layout.addWidget(reset_tab_btn)

        reset_all_btn = QPushButton("Reset All")
        reset_all_btn.clicked.connect(self.reset_all_settings)
        buttons_layout.addWidget(reset_all_btn)

        buttons_layout.addStretch(1)

        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(self.apply_settings)
        buttons_layout.addWidget(apply_btn)

        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        buttons_layout.addWidget(ok_btn)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        buttons_layout.addWidget(cancel_btn)

        layout.addLayout(buttons_layout)

        # Initialize UI from settings
        self.update_ui_from_settings(current_settings)

    def update_ui_from_settings(self, settings):
        for i in range(self.tab_widget.count()):
            self.tab_widget.widget(i).update_ui(settings)

    def get_settings(self):
        settings = {}
        for i in range(self.tab_widget.count()):
            settings.update(self.tab_widget.widget(i).get_settings())
        return settings

    def reset_current_tab(self):
        self.tab_widget.currentWidget().reset_to_defaults()
        QMessageBox.information(
            self,
            "Reset Complete",
            f"Settings for '{self.tab_widget.tabText(self.tab_widget.currentIndex())}' tab have been reset to defaults.",
        )

    def reset_all_settings(self):
        reply = QMessageBox.question(
            self,
            "Reset All Settings",
            "Are you sure you want to reset all settings to defaults?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.update_ui_from_settings(self.default_settings)
            self.apply_settings()

    def apply_settings(self):
        if not self.parent_window:
            return

        settings = self.get_settings()
        self.parent_window.init_manager.settings.update(settings)

        if hasattr(self.parent_window.init_manager, "settings_dirty"):
            self.parent_window.init_manager.settings_dirty = True
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error("REPORT ERROR: Missing attribute 'settings_dirty' on object")

        # Persist to disk immediately
        if hasattr(self.parent_window, "init_manager"):
            self.parent_window.init_manager.save_settings()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'init_manager' on self.parent_window"
            )

        if hasattr(self.parent_window.view_3d_manager, "apply_3d_settings"):
            self.parent_window.view_3d_manager.apply_3d_settings()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'apply_3d_settings' on object"
            )

        if hasattr(self.parent_window.init_manager, "update_cpk_colors_from_settings"):
            self.parent_window.init_manager.update_cpk_colors_from_settings()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'update_cpk_colors_from_settings' on object"
            )

        # Redraw molecule
        current_mol = getattr(self.parent_window.view_3d_manager, "current_mol", None)
        if current_mol and hasattr(
            self.parent_window.view_3d_manager, "draw_molecule_3d"
        ):
            self.parent_window.view_3d_manager.draw_molecule_3d(current_mol)

        # Apply 2D view settings
        scene = getattr(self.parent_window.init_manager, "scene", None)
        if scene:
            bg_col_2d = settings.get("background_color_2d", "#FFFFFF")
            scene.setBackgroundBrush(QColor(bg_col_2d))
            for item in scene.items():
                if hasattr(item, "update_style"):
                    item.update_style()
                elif hasattr(item, "update"):
                    item.update()
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error("REPORT ERROR: Missing attribute 'update' on item")

        if hasattr(self.parent_window, "statusBar") and self.parent_window.statusBar():
            self.parent_window.statusBar().showMessage("Settings applied successfully")

    def accept(self):
        self.apply_settings()
        super().accept()
