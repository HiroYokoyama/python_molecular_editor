#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software
Refactored SettingsDialog (Phase 2)
"""

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
    from ..utils.constants import CPK_COLORS
except ImportError:
    pass


class SettingsDialog(QDialog):
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumSize(650, 750)
        self.parent_window = parent

        self.default_settings = {
            "background_color": "#919191",
            "projection_mode": "Perspective",
            "lighting_enabled": True,
            "specular": 0.20,
            "specular_power": 20,
            "light_intensity": 1.0,
            "show_3d_axes": True,
            "ball_stick_atom_scale": 1.0,
            "ball_stick_bond_radius": 0.1,
            "ball_stick_resolution": 16,
            "cpk_atom_scale": 1.0,
            "cpk_resolution": 32,
            "wireframe_bond_radius": 0.01,
            "wireframe_resolution": 6,
            "stick_bond_radius": 0.15,
            "stick_resolution": 16,
            "ball_stick_double_bond_offset_factor": 2.0,
            "ball_stick_triple_bond_offset_factor": 2.0,
            "ball_stick_double_bond_radius_factor": 0.8,
            "ball_stick_triple_bond_radius_factor": 0.75,
            "wireframe_double_bond_offset_factor": 3.0,
            "wireframe_triple_bond_offset_factor": 3.0,
            "wireframe_double_bond_radius_factor": 0.8,
            "wireframe_triple_bond_radius_factor": 0.75,
            "stick_double_bond_offset_factor": 1.5,
            "stick_triple_bond_offset_factor": 1.0,
            "stick_double_bond_radius_factor": 0.6,
            "stick_triple_bond_radius_factor": 0.4,
            "aromatic_torus_thickness_factor": 0.6,
            "display_aromatic_circles_3d": False,
            "skip_chemistry_checks": False,
            "always_ask_charge": False,
            "3d_conversion_mode": "fallback",
            "optimization_method": "MMFF_RDKIT",
            "ball_stick_bond_color": "#7F7F7F",
            "cpk_colors": {},
            "display_kekule_3d": False,
            "ball_stick_use_cpk_bond_color": False,
            "bond_width_2d": 2.0,
            "bond_spacing_double_2d": 3.5,
            "bond_spacing_triple_2d": 3.5,
            "atom_font_size_2d": 20,
            "background_color_2d": "#FFFFFF",
            "bond_color_2d": "#222222",
            "atom_use_bond_color_2d": False,
            "bond_cap_style_2d": "Round",
            "bond_wedge_width_2d": 6.0,
            "bond_dash_count_2d": 8,
            "atom_font_family_2d": "Arial",
        }

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
        self.parent_window.settings.update(settings)

        if hasattr(self.parent_window, "settings_dirty"):
            self.parent_window.settings_dirty = True

        if hasattr(self.parent_window.view_3d_manager, "apply_3d_settings"):
            self.parent_window.view_3d_manager.apply_3d_settings()

        if hasattr(self.parent_window.init_manager, "update_cpk_colors_from_settings"):
            self.parent_window.init_manager.update_cpk_colors_from_settings()

        # Redraw molecule
        current_mol = getattr(self.parent_window, "current_mol", None)
        if current_mol and hasattr(self.parent_window.view_3d_manager, "draw_molecule_3d"):
            self.parent_window.view_3d_manager.draw_molecule_3d(current_mol)

        # Apply 2D view settings
        scene = getattr(self.parent_window, "scene", None)
        if scene:
            bg_col_2d = settings.get("background_color_2d", "#FFFFFF")
            scene.setBackgroundBrush(QColor(bg_col_2d))
            for item in scene.items():
                if hasattr(item, "update_style"):
                    item.update_style()
                elif hasattr(item, "update"):
                    item.update()

        if hasattr(self.parent_window, "statusBar") and self.parent_window.statusBar():
            self.parent_window.statusBar().showMessage("Settings applied successfully")

    def accept(self):
        self.apply_settings()
        super().accept()
