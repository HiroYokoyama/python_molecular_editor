#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QTabWidget, QWidget, QFormLayout, QPushButton, QHBoxLayout,
    QCheckBox, QComboBox, QLabel, QColorDialog, QSlider, QFrame, QMessageBox
)
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
try:
    from .constants import CPK_COLORS
except Exception:
    from modules.constants import CPK_COLORS


class SettingsDialog(QDialog):
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumSize(700, 800)
        
        # 親ウィンドウの参照を保存（Apply機能のため）
        self.parent_window = parent
        
        # デフォルト設定をクラス内で定義
        # Multi-bond settings are model-specific now (ball_stick, cpk, wireframe, stick)
        self.default_settings = {
            'background_color': '#919191',
            'projection_mode': 'Perspective',
            'lighting_enabled': True,
            'specular': 0.20,
            'specular_power': 20,
            'light_intensity': 1.0,
            'show_3d_axes': True,
            # Ball and Stick model parameters
            'ball_stick_atom_scale': 1.0,
            'ball_stick_bond_radius': 0.1,
            'ball_stick_resolution': 16,
            # CPK (Space-filling) model parameters
            'cpk_atom_scale': 1.0,
            'cpk_resolution': 32,
            # Wireframe model parameters
            'wireframe_bond_radius': 0.01,
            'wireframe_resolution': 6,
            # Stick model parameters
            'stick_bond_radius': 0.15,
            'stick_resolution': 16,
            # Multiple bond offset parameters (per-model)
            'ball_stick_double_bond_offset_factor': 2.0,
            'ball_stick_triple_bond_offset_factor': 2.0,
            'ball_stick_double_bond_radius_factor': 0.8,
            'ball_stick_triple_bond_radius_factor': 0.75,
            'wireframe_double_bond_offset_factor': 3.0,
            'wireframe_triple_bond_offset_factor': 3.0,
            'wireframe_double_bond_radius_factor': 0.8,
            'wireframe_triple_bond_radius_factor': 0.75,
            'stick_double_bond_offset_factor': 1.5,
            'stick_triple_bond_offset_factor': 1.0,
            'stick_double_bond_radius_factor': 0.6,
            'stick_triple_bond_radius_factor': 0.4,
            'aromatic_torus_thickness_factor': 0.6,
            # Whether to draw an aromatic circle inside rings in 3D
            'display_aromatic_circles_3d': False,
            # If True, attempts to be permissive when RDKit raises chemical/sanitization errors
            # during file import (useful for viewing malformed XYZ/MOL files). When enabled,
            # element symbol recognition will be coerced where possible and Chem.SanitizeMol
            # failures will be ignored so the 3D viewer can still display the structure.
            'skip_chemistry_checks': False,
            # When True, always prompt the user for molecular charge on XYZ import
            # instead of silently trying charge=0 first. Default True to disable
            # the silent 'charge=0' test.
            'always_ask_charge': False,
            # 3D conversion/optimization defaults
            '3d_conversion_mode': 'fallback',
            'optimization_method': 'MMFF_RDKIT',
            'ball_stick_bond_color': '#7F7F7F',
            'cpk_colors': {},
            # If True, RDKit will attempt to kekulize aromatic systems for 3D display
            # (shows alternating single/double bonds rather than aromatic circles)
            'display_kekule_3d': False,
            'ball_stick_use_cpk_bond_color': False,
            
            # --- 2D Settings Defaults ---
            'bond_width_2d': 2.0,
            'bond_spacing_double_2d': 3.5,
            'bond_spacing_triple_2d': 3.5,
            'atom_font_size_2d': 20,
            'background_color_2d': '#FFFFFF',
            'bond_color_2d': '#222222', # Almost black
            'atom_use_bond_color_2d': False,
            'bond_cap_style_2d': 'Round',
        }
        
        # --- 選択された色を管理する専用のインスタンス変数 ---
        self.current_bg_color = None

        # --- UI要素の作成 ---
        layout = QVBoxLayout(self)
        
        # タブウィジェットを作成
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)

        # 2D Settings Tab (First)
        self.create_2d_settings_tab()

        # Scene設定タブ
        self.create_scene_tab()
        
        # Ball and Stick設定タブ
        self.create_ball_stick_tab()
        
        # CPK設定タブ
        self.create_cpk_tab()
        
        # Wireframe設定タブ
        self.create_wireframe_tab()
        
        # Stick設定タブ
        self.create_stick_tab()

        # Other設定タブ
        self.create_other_tab()

        # 渡された設定でUIと内部変数を初期化
        self.update_ui_from_settings(current_settings)
        
        # Initialize aromatic circle checkbox and torus thickness from settings
        self.aromatic_circle_checkbox.setChecked(current_settings.get('display_aromatic_circles_3d', self.default_settings.get('display_aromatic_circles_3d', False)))
        # Thickness factor is stored as a multiplier (e.g., 1.0), slider uses integer 0-300 representing 0.1x-3.0x
        thickness_factor = current_settings.get('aromatic_torus_thickness_factor', self.default_settings.get('aromatic_torus_thickness_factor', 1.0))
        self.aromatic_torus_thickness_slider.setValue(int(thickness_factor * 100))
        self.aromatic_torus_thickness_label.setText(f"{thickness_factor:.1f}")

        # Initialize internal 2D colors
        self.current_bg_color_2d = current_settings.get('background_color_2d', self.default_settings['background_color_2d'])
        self.current_bond_color_2d = current_settings.get('bond_color_2d', self.default_settings['bond_color_2d'])
        
        # Apply initial 2D button styles
        self.update_2d_color_buttons()

        # --- ボタンの配置 ---
        buttons = QHBoxLayout()
        
        # タブごとのリセットボタン
        reset_tab_button = QPushButton("Reset Current Tab")
        reset_tab_button.clicked.connect(self.reset_current_tab)
        reset_tab_button.setToolTip("Reset settings for the currently selected tab only")
        buttons.addWidget(reset_tab_button)
        
        # 全体リセットボタン
        reset_all_button = QPushButton("Reset All")
        reset_all_button.clicked.connect(self.reset_all_settings)
        reset_all_button.setToolTip("Reset all settings to defaults")
        buttons.addWidget(reset_all_button)
        
        buttons.addStretch(1)
        
        # Applyボタンを追加
        apply_button = QPushButton("Apply")
        apply_button.clicked.connect(self.apply_settings)
        apply_button.setToolTip("Apply settings without closing dialog")
        buttons.addWidget(apply_button)
        
        ok_button = QPushButton("OK")
        cancel_button = QPushButton("Cancel")
        ok_button.clicked.connect(self.accept)
        cancel_button.clicked.connect(self.reject)

        buttons.addWidget(ok_button)
        buttons.addWidget(cancel_button)
        layout.addLayout(buttons)

    def create_2d_settings_tab(self):
        """2D Settings Tab"""
        widget = QWidget()
        form_layout = QFormLayout(widget)

        # --- View Settings ---
        form_layout.addRow(QLabel("<b>View Appearance</b>"))
        
        # Background Color
        self.bg_color_2d_button = QPushButton()
        self.bg_color_2d_button.setFixedSize(60, 24)
        self.bg_color_2d_button.clicked.connect(self.pick_bg_color_2d)
        form_layout.addRow("Background Color:", self.bg_color_2d_button)

        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # --- Bond Settings ---
        form_layout.addRow(QLabel("<b>Bond Settings</b>"))

        # Bond Color
        self.bond_color_2d_button = QPushButton()
        self.bond_color_2d_button.setFixedSize(60, 24)
        self.bond_color_2d_button.clicked.connect(self.pick_bond_color_2d)
        form_layout.addRow("Bond Color:", self.bond_color_2d_button)

        # Bond Width
        self.bond_width_2d_slider = QSlider(Qt.Orientation.Horizontal)
        self.bond_width_2d_slider.setRange(10, 200) # 1.0 - 20.0
        self.bond_width_2d_label = QLabel("2.0")
        self.bond_width_2d_slider.valueChanged.connect(lambda v: self.bond_width_2d_label.setText(f"{v/10:.1f}"))
        bw_layout = QHBoxLayout()
        bw_layout.addWidget(self.bond_width_2d_slider)
        bw_layout.addWidget(self.bond_width_2d_label)
        form_layout.addRow("Bond Width:", bw_layout)

        # Double Bond Spacing
        self.bond_spacing_double_2d_slider = QSlider(Qt.Orientation.Horizontal)
        self.bond_spacing_double_2d_slider.setRange(10, 200) # 1.0 - 20.0
        self.bond_spacing_double_2d_label = QLabel("3.5")
        self.bond_spacing_double_2d_slider.valueChanged.connect(lambda v: self.bond_spacing_double_2d_label.setText(f"{v/10:.1f}"))
        bsd_layout = QHBoxLayout()
        bsd_layout.addWidget(self.bond_spacing_double_2d_slider)
        bsd_layout.addWidget(self.bond_spacing_double_2d_label)
        form_layout.addRow("Double Bond Spacing:", bsd_layout)

        # Triple Bond Spacing
        self.bond_spacing_triple_2d_slider = QSlider(Qt.Orientation.Horizontal)
        self.bond_spacing_triple_2d_slider.setRange(10, 200) # 1.0 - 20.0
        self.bond_spacing_triple_2d_label = QLabel("3.5")
        self.bond_spacing_triple_2d_slider.valueChanged.connect(lambda v: self.bond_spacing_triple_2d_label.setText(f"{v/10:.1f}"))
        bst_layout = QHBoxLayout()
        bst_layout.addWidget(self.bond_spacing_triple_2d_slider)
        bst_layout.addWidget(self.bond_spacing_triple_2d_label)
        form_layout.addRow("Triple Bond Spacing:", bst_layout)
        
        # Bond Cap Style
        self.bond_cap_style_2d_combo = QComboBox()
        self.bond_cap_style_2d_combo.addItems(['Round', 'Flat', 'Square'])
        form_layout.addRow("Bond Cap Style:", self.bond_cap_style_2d_combo)


        line2 = QFrame()
        line2.setFrameShape(QFrame.Shape.HLine)
        line2.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line2)

        # --- Atom Settings ---
        form_layout.addRow(QLabel("<b>Atom Settings</b>"))

        # Font Size
        self.atom_font_size_2d_slider = QSlider(Qt.Orientation.Horizontal)
        self.atom_font_size_2d_slider.setRange(8, 72)
        self.atom_font_size_2d_label = QLabel("20")
        self.atom_font_size_2d_slider.valueChanged.connect(lambda v: self.atom_font_size_2d_label.setText(str(v)))
        fs_layout = QHBoxLayout()
        fs_layout.addWidget(self.atom_font_size_2d_slider)
        fs_layout.addWidget(self.atom_font_size_2d_label)
        form_layout.addRow("Atom Label Font Size:", fs_layout)
        
        # Use Bond Color Checkbox
        self.atom_use_bond_color_2d_checkbox = QCheckBox()
        self.atom_use_bond_color_2d_checkbox.setToolTip("If checked, atoms will use the unified Bond Color instead of element-specific colors (CPK).")
        form_layout.addRow("Use Bond Color for Atoms:", self.atom_use_bond_color_2d_checkbox)

        self.tab_widget.addTab(widget, "2D Settings")
        

        
        # Initialize internal variable for combo updates
        # (Usually done in update_ui_from_settings, but we add direct access here for simple binding if needed)
        # We rely on update_ui_from_settings to set current selection.


    def pick_bg_color_2d(self):
        color = QColorDialog.getColor(QColor(self.current_bg_color_2d), self, "Select 2D Background Color")
        if color.isValid():
            self.current_bg_color_2d = color.name()
            self.update_2d_color_buttons()

    def pick_bond_color_2d(self):
        color = QColorDialog.getColor(QColor(self.current_bond_color_2d), self, "Select 2D Bond Color")
        if color.isValid():
            self.current_bond_color_2d = color.name()
            self.update_2d_color_buttons()

    def update_2d_color_buttons(self):
        try:
            self.bg_color_2d_button.setStyleSheet(f"background-color: {self.current_bg_color_2d}; border: 1px solid #888;")
            self.bond_color_2d_button.setStyleSheet(f"background-color: {self.current_bond_color_2d}; border: 1px solid #888;")
        except Exception:
            pass
    
    def create_scene_tab(self):
        """基本設定タブを作成"""
        scene_widget = QWidget()
        form_layout = QFormLayout(scene_widget)

        # 1. 背景色
        self.bg_button = QPushButton()
        self.bg_button.setToolTip("Click to select a color")
        self.bg_button.clicked.connect(self.select_color)
        form_layout.addRow("Background Color:", self.bg_button)

        # 1a. 軸の表示/非表示
        self.axes_checkbox = QCheckBox()
        form_layout.addRow("Show 3D Axes:", self.axes_checkbox)

        # 2. ライトの有効/無効
        self.light_checkbox = QCheckBox()
        form_layout.addRow("Enable Lighting:", self.light_checkbox)

        # 光の強さスライダーを追加
        self.intensity_slider = QSlider(Qt.Orientation.Horizontal)
        self.intensity_slider.setRange(0, 200) # 0.0 ~ 2.0 の範囲
        self.intensity_label = QLabel("1.0")
        self.intensity_slider.valueChanged.connect(lambda v: self.intensity_label.setText(f"{v/100:.2f}"))
        intensity_layout = QHBoxLayout()
        intensity_layout.addWidget(self.intensity_slider)
        intensity_layout.addWidget(self.intensity_label)
        form_layout.addRow("Light Intensity:", intensity_layout)

        # 3. 光沢 (Specular)
        self.specular_slider = QSlider(Qt.Orientation.Horizontal)
        self.specular_slider.setRange(0, 100)
        self.specular_label = QLabel("0.20")
        self.specular_slider.valueChanged.connect(lambda v: self.specular_label.setText(f"{v/100:.2f}"))
        specular_layout = QHBoxLayout()
        specular_layout.addWidget(self.specular_slider)
        specular_layout.addWidget(self.specular_label)
        form_layout.addRow("Shininess (Specular):", specular_layout)
        
        # 4. 光沢の強さ (Specular Power)
        self.spec_power_slider = QSlider(Qt.Orientation.Horizontal)
        self.spec_power_slider.setRange(0, 100)
        self.spec_power_label = QLabel("20")
        self.spec_power_slider.valueChanged.connect(lambda v: self.spec_power_label.setText(str(v)))
        spec_power_layout = QHBoxLayout()
        spec_power_layout.addWidget(self.spec_power_slider)
        spec_power_layout.addWidget(self.spec_power_label)
        form_layout.addRow("Shininess Power:", spec_power_layout)
        
        # Projection mode (Perspective / Orthographic)
        self.projection_combo = QComboBox()
        self.projection_combo.addItem("Perspective")
        self.projection_combo.addItem("Orthographic")
        self.projection_combo.setToolTip("Choose camera projection mode: Perspective (default) or Orthographic")
        form_layout.addRow("Projection Mode:", self.projection_combo)
        
        self.tab_widget.addTab(scene_widget, "3D Scene")
    
    def create_other_tab(self):
        """other設定タブを作成"""
        self.other_widget = QWidget()
        self.other_form_layout = QFormLayout(self.other_widget)

        # 化学チェックスキップオプション（otherタブに移動）
        self.skip_chem_checks_checkbox = QCheckBox()
        self.skip_chem_checks_checkbox.setToolTip("When enabled, XYZ file import will try to ignore chemical/sanitization errors and allow viewing malformed files.")
        # Immediately persist change to settings when user toggles the checkbox
        try:
            self.skip_chem_checks_checkbox.stateChanged.connect(lambda s: self._on_skip_chem_checks_changed(s))
        except Exception:
            pass

        # Add the checkbox to the other tab's form
        try:
            self.other_form_layout.addRow("Skip chemistry checks on import XYZ file:", self.skip_chem_checks_checkbox)
        except Exception:
            pass

        # 3D Kekule display option (under Other) will be added below the
        # 'Always ask molecular charge on import' option so ordering is clear
        # in the UI.
        self.kekule_3d_checkbox = QCheckBox()
        self.kekule_3d_checkbox.setToolTip("When enabled, aromatic bonds will be kekulized in the 3D view (show alternating single/double bonds).")
        # Don't persist kekule state immediately; Apply/OK should commit setting.
        # Always ask charge on XYZ import (skip silent charge=0 test)
        self.always_ask_charge_checkbox = QCheckBox()
        self.always_ask_charge_checkbox.setToolTip("Prompt for overall molecular charge when importing XYZ files instead of silently trying charge=0 first.")
        try:
            self.other_form_layout.addRow("Always ask molecular charge on import XYZ file:", self.always_ask_charge_checkbox)
        except Exception:
            pass

        # Add separator after Kekule bonds option
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)
        self.other_form_layout.addRow(separator)

        # Place the Kekulé option after the always-ask-charge option
        try:
            self.other_form_layout.addRow("Display Kekulé bonds in 3D:", self.kekule_3d_checkbox)
        except Exception:
            pass
        

        # Aromatic ring circle display option
        self.aromatic_circle_checkbox = QCheckBox()
        self.aromatic_circle_checkbox.setToolTip("When enabled, aromatic rings will be displayed with a circle inside the ring in 3D view.")
        try:
            self.other_form_layout.addRow("Display aromatic rings as circles in 3D:", self.aromatic_circle_checkbox)
        except Exception:
            pass
        
        # Aromatic torus thickness factor
        self.aromatic_torus_thickness_slider = QSlider(Qt.Orientation.Horizontal)
        self.aromatic_torus_thickness_slider.setRange(10, 300)  # 0.1x to 3.0x
        self.aromatic_torus_thickness_slider.setValue(60)  # Default 0.6x
        self.aromatic_torus_thickness_label = QLabel("0.6")
        self.aromatic_torus_thickness_slider.valueChanged.connect(
            lambda v: self.aromatic_torus_thickness_label.setText(f"{v/100:.1f}")
        )
        thickness_layout = QHBoxLayout()
        thickness_layout.addWidget(self.aromatic_torus_thickness_slider)
        thickness_layout.addWidget(self.aromatic_torus_thickness_label)
        try:
            self.other_form_layout.addRow("Aromatic torus thickness (× bond radius):", thickness_layout)
        except Exception:
            pass

        # Add Other tab to the tab widget
        self.tab_widget.addTab(self.other_widget, "Other")

    def refresh_ui(self):
        """Refresh periodic table / BS button visuals using current settings.

        Called when settings change externally (e.g., Reset All in main settings) so
        the dialog reflects the current stored overrides.
        """
        try:
            # Update element button colors from parent.settings cpks
            overrides = self.parent_window.settings.get('cpk_colors', {}) if self.parent_window and hasattr(self.parent_window, 'settings') else {}
            for s, btn in self.element_buttons.items():
                try:
                    override = overrides.get(s)
                    q_color = QColor(override) if override else CPK_COLORS.get(s, CPK_COLORS['DEFAULT'])
                    brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
                    text_color = 'white' if brightness < 128 else 'black'
                    btn.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
                except Exception:
                    pass
            # Update BS color button from parent settings
            try:
                if hasattr(self, 'bs_button') and self.parent_window and hasattr(self.parent_window, 'settings'):
                    bs_hex = self.parent_window.settings.get('ball_stick_bond_color', self.parent_window.default_settings.get('ball_stick_bond_color', '#7F7F7F'))
                    self.bs_button.setStyleSheet(f"background-color: {bs_hex}; border: 1px solid #888;")
                    self.bs_button.setToolTip(bs_hex)
            except Exception:
                pass
        except Exception:
            pass
        # Avoid circular import: import ColorSettingsDialog only inside method using it

        # NOTE: Multi-bond offset/thickness settings moved to per-model tabs to allow
        # independent configuration for Ball&Stick/CPK/Wireframe/Stick.
                
        # 'Other' tab is created in create_other_tab; nothing to do here.
    
    def create_ball_stick_tab(self):
        """Ball and Stick設定タブを作成"""
        ball_stick_widget = QWidget()
        form_layout = QFormLayout(ball_stick_widget)

        info_label = QLabel("Ball & Stick model shows atoms as spheres and bonds as cylinders.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow(info_label)
        
        # 原子サイズスケール
        self.bs_atom_scale_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_atom_scale_slider.setRange(10, 200)  # 0.1 ~ 2.0
        self.bs_atom_scale_label = QLabel("1.00")
        self.bs_atom_scale_slider.valueChanged.connect(lambda v: self.bs_atom_scale_label.setText(f"{v/100:.2f}"))
        atom_scale_layout = QHBoxLayout()
        atom_scale_layout.addWidget(self.bs_atom_scale_slider)
        atom_scale_layout.addWidget(self.bs_atom_scale_label)
        form_layout.addRow("Atom Size Scale:", atom_scale_layout)
        
        # ボンド半径
        self.bs_bond_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_bond_radius_slider.setRange(1, 50)  # 0.01 ~ 0.5
        self.bs_bond_radius_label = QLabel("0.10")
        self.bs_bond_radius_slider.valueChanged.connect(lambda v: self.bs_bond_radius_label.setText(f"{v/100:.2f}"))
        bond_radius_layout = QHBoxLayout()
        bond_radius_layout.addWidget(self.bs_bond_radius_slider)
        bond_radius_layout.addWidget(self.bs_bond_radius_label)
        form_layout.addRow("Bond Radius:", bond_radius_layout)

        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # --- Per-model multi-bond controls (Ball & Stick) ---
        # 二重/三重結合のオフセット倍率（Ball & Stick）
        self.bs_double_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_double_offset_slider.setRange(100, 400)
        self.bs_double_offset_label = QLabel("2.00")
        self.bs_double_offset_slider.valueChanged.connect(lambda v: self.bs_double_offset_label.setText(f"{v/100:.2f}"))
        bs_double_offset_layout = QHBoxLayout()
        bs_double_offset_layout.addWidget(self.bs_double_offset_slider)
        bs_double_offset_layout.addWidget(self.bs_double_offset_label)
        form_layout.addRow("Double Bond Offset (Ball & Stick):", bs_double_offset_layout)

        self.bs_triple_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_triple_offset_slider.setRange(100, 400)
        self.bs_triple_offset_label = QLabel("2.00")
        self.bs_triple_offset_slider.valueChanged.connect(lambda v: self.bs_triple_offset_label.setText(f"{v/100:.2f}"))
        bs_triple_offset_layout = QHBoxLayout()
        bs_triple_offset_layout.addWidget(self.bs_triple_offset_slider)
        bs_triple_offset_layout.addWidget(self.bs_triple_offset_label)
        form_layout.addRow("Triple Bond Offset (Ball & Stick):", bs_triple_offset_layout)

        # 半径倍率
        self.bs_double_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_double_radius_slider.setRange(50, 100)
        self.bs_double_radius_label = QLabel("0.80")
        self.bs_double_radius_slider.valueChanged.connect(lambda v: self.bs_double_radius_label.setText(f"{v/100:.2f}"))
        bs_double_radius_layout = QHBoxLayout()
        bs_double_radius_layout.addWidget(self.bs_double_radius_slider)
        bs_double_radius_layout.addWidget(self.bs_double_radius_label)
        form_layout.addRow("Double Bond Thickness (Ball & Stick):", bs_double_radius_layout)

        self.bs_triple_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_triple_radius_slider.setRange(50, 100)
        self.bs_triple_radius_label = QLabel("0.70")
        self.bs_triple_radius_slider.valueChanged.connect(lambda v: self.bs_triple_radius_label.setText(f"{v/100:.2f}"))
        bs_triple_radius_layout = QHBoxLayout()
        bs_triple_radius_layout.addWidget(self.bs_triple_radius_slider)
        bs_triple_radius_layout.addWidget(self.bs_triple_radius_label)
        form_layout.addRow("Triple Bond Thickness (Ball & Stick):", bs_triple_radius_layout)

        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # 解像度
        self.bs_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_resolution_slider.setRange(6, 32)
        self.bs_resolution_label = QLabel("16")
        self.bs_resolution_slider.valueChanged.connect(lambda v: self.bs_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.bs_resolution_slider)
        resolution_layout.addWidget(self.bs_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)

        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # --- Ball & Stick bond color ---
        self.bs_bond_color_button = QPushButton()
        self.bs_bond_color_button.setFixedSize(36, 24)
        self.bs_bond_color_button.clicked.connect(self.pick_bs_bond_color)
        self.bs_bond_color_button.setToolTip("Choose the uniform bond color for Ball & Stick model (3D)")
        form_layout.addRow("Ball & Stick bond color:", self.bs_bond_color_button)

        # Use CPK colors for bonds option
        self.bs_use_cpk_bond_checkbox = QCheckBox()
        self.bs_use_cpk_bond_checkbox.setToolTip("If checked, bonds will be colored using the atom colors (split bonds). If unchecked, a uniform color is used.")
        form_layout.addRow("Use CPK colors for bonds:", self.bs_use_cpk_bond_checkbox)

        self.tab_widget.addTab(ball_stick_widget, "Ball & Stick")
    
    def create_cpk_tab(self):
        """CPK設定タブを作成"""
        cpk_widget = QWidget()
        form_layout = QFormLayout(cpk_widget)

        info_label = QLabel("CPK model shows atoms as space-filling spheres using van der Waals radii.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow(info_label)
        
        # 原子サイズスケール
        self.cpk_atom_scale_slider = QSlider(Qt.Orientation.Horizontal)
        self.cpk_atom_scale_slider.setRange(50, 200)  # 0.5 ~ 2.0
        self.cpk_atom_scale_label = QLabel("1.00")
        self.cpk_atom_scale_slider.valueChanged.connect(lambda v: self.cpk_atom_scale_label.setText(f"{v/100:.2f}"))
        atom_scale_layout = QHBoxLayout()
        atom_scale_layout.addWidget(self.cpk_atom_scale_slider)
        atom_scale_layout.addWidget(self.cpk_atom_scale_label)
        form_layout.addRow("Atom Size Scale:", atom_scale_layout)
        
        # 解像度
        self.cpk_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.cpk_resolution_slider.setRange(8, 64)
        self.cpk_resolution_label = QLabel("32")
        self.cpk_resolution_slider.valueChanged.connect(lambda v: self.cpk_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.cpk_resolution_slider)
        resolution_layout.addWidget(self.cpk_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)

        self.tab_widget.addTab(cpk_widget, "CPK (Space-filling)")
    
    def create_wireframe_tab(self):
        """Wireframe設定タブを作成"""
        wireframe_widget = QWidget()
        form_layout = QFormLayout(wireframe_widget)

        info_label = QLabel("Wireframe model shows molecular structure with thin lines only.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow(info_label)
        
        # ボンド半径
        self.wf_bond_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_bond_radius_slider.setRange(1, 10)  # 0.01 ~ 0.1
        self.wf_bond_radius_label = QLabel("0.01")
        self.wf_bond_radius_slider.valueChanged.connect(lambda v: self.wf_bond_radius_label.setText(f"{v/100:.2f}"))
        bond_radius_layout = QHBoxLayout()
        bond_radius_layout.addWidget(self.wf_bond_radius_slider)
        bond_radius_layout.addWidget(self.wf_bond_radius_label)
        form_layout.addRow("Bond Radius:", bond_radius_layout)


        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)
            
        # --- Per-model multi-bond controls (Wireframe) ---
        self.wf_double_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_double_offset_slider.setRange(100, 400)
        self.wf_double_offset_label = QLabel("2.00")
        self.wf_double_offset_slider.valueChanged.connect(lambda v: self.wf_double_offset_label.setText(f"{v/100:.2f}"))
        wf_double_offset_layout = QHBoxLayout()
        wf_double_offset_layout.addWidget(self.wf_double_offset_slider)
        wf_double_offset_layout.addWidget(self.wf_double_offset_label)
        form_layout.addRow("Double Bond Offset (Wireframe):", wf_double_offset_layout)
    
        self.wf_triple_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_triple_offset_slider.setRange(100, 400)
        self.wf_triple_offset_label = QLabel("2.00")
        self.wf_triple_offset_slider.valueChanged.connect(lambda v: self.wf_triple_offset_label.setText(f"{v/100:.2f}"))
        wf_triple_offset_layout = QHBoxLayout()
        wf_triple_offset_layout.addWidget(self.wf_triple_offset_slider)
        wf_triple_offset_layout.addWidget(self.wf_triple_offset_label)
        form_layout.addRow("Triple Bond Offset (Wireframe):", wf_triple_offset_layout)
    
        self.wf_double_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_double_radius_slider.setRange(50, 100)
        self.wf_double_radius_label = QLabel("0.80")
        self.wf_double_radius_slider.valueChanged.connect(lambda v: self.wf_double_radius_label.setText(f"{v/100:.2f}"))
        wf_double_radius_layout = QHBoxLayout()
        wf_double_radius_layout.addWidget(self.wf_double_radius_slider)
        wf_double_radius_layout.addWidget(self.wf_double_radius_label)
        form_layout.addRow("Double Bond Thickness (Wireframe):", wf_double_radius_layout)
    
        self.wf_triple_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_triple_radius_slider.setRange(50, 100)
        self.wf_triple_radius_label = QLabel("0.70")
        self.wf_triple_radius_slider.valueChanged.connect(lambda v: self.wf_triple_radius_label.setText(f"{v/100:.2f}"))
        wf_triple_radius_layout = QHBoxLayout()
        wf_triple_radius_layout.addWidget(self.wf_triple_radius_slider)
        wf_triple_radius_layout.addWidget(self.wf_triple_radius_label)
        form_layout.addRow("Triple Bond Thickness (Wireframe):", wf_triple_radius_layout)

        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # 解像度
        self.wf_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_resolution_slider.setRange(4, 16)
        self.wf_resolution_label = QLabel("6")
        self.wf_resolution_slider.valueChanged.connect(lambda v: self.wf_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.wf_resolution_slider)
        resolution_layout.addWidget(self.wf_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)
    
        self.tab_widget.addTab(wireframe_widget, "Wireframe")
    
    def create_stick_tab(self):
        """Stick設定タブを作成"""
        stick_widget = QWidget()
        form_layout = QFormLayout(stick_widget)

        info_label = QLabel("Stick model shows bonds as thick cylinders with atoms as small spheres.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow(info_label)
        
        # ボンド半径（原子半径も同じ値を使用）
        self.stick_bond_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_bond_radius_slider.setRange(5, 50)  # 0.05 ~ 0.5
        self.stick_bond_radius_label = QLabel("0.15")
        self.stick_bond_radius_slider.valueChanged.connect(lambda v: self.stick_bond_radius_label.setText(f"{v/100:.2f}"))
        bond_radius_layout = QHBoxLayout()
        bond_radius_layout.addWidget(self.stick_bond_radius_slider)
        bond_radius_layout.addWidget(self.stick_bond_radius_label)
        form_layout.addRow("Bond Radius:", bond_radius_layout)
        
        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)

        # --- Per-model multi-bond controls (Stick) ---
        self.stick_double_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_double_offset_slider.setRange(50, 400)
        self.stick_double_offset_label = QLabel("2.00")
        self.stick_double_offset_slider.valueChanged.connect(lambda v: self.stick_double_offset_label.setText(f"{v/100:.2f}"))
        stick_double_offset_layout = QHBoxLayout()
        stick_double_offset_layout.addWidget(self.stick_double_offset_slider)
        stick_double_offset_layout.addWidget(self.stick_double_offset_label)
        form_layout.addRow("Double Bond Offset (Stick):", stick_double_offset_layout)

        self.stick_triple_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_triple_offset_slider.setRange(50, 400)
        self.stick_triple_offset_label = QLabel("2.00")
        self.stick_triple_offset_slider.valueChanged.connect(lambda v: self.stick_triple_offset_label.setText(f"{v/100:.2f}"))
        stick_triple_offset_layout = QHBoxLayout()
        stick_triple_offset_layout.addWidget(self.stick_triple_offset_slider)
        stick_triple_offset_layout.addWidget(self.stick_triple_offset_label)
        form_layout.addRow("Triple Bond Offset (Stick):", stick_triple_offset_layout)

        self.stick_double_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_double_radius_slider.setRange(20, 100)
        self.stick_double_radius_label = QLabel("0.80")
        self.stick_double_radius_slider.valueChanged.connect(lambda v: self.stick_double_radius_label.setText(f"{v/100:.2f}"))
        stick_double_radius_layout = QHBoxLayout()
        stick_double_radius_layout.addWidget(self.stick_double_radius_slider)
        stick_double_radius_layout.addWidget(self.stick_double_radius_label)
        form_layout.addRow("Double Bond Thickness (Stick):", stick_double_radius_layout)

        self.stick_triple_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_triple_radius_slider.setRange(20, 100)
        self.stick_triple_radius_label = QLabel("0.70")
        self.stick_triple_radius_slider.valueChanged.connect(lambda v: self.stick_triple_radius_label.setText(f"{v/100:.2f}"))
        stick_triple_radius_layout = QHBoxLayout()
        stick_triple_radius_layout.addWidget(self.stick_triple_radius_slider)
        stick_triple_radius_layout.addWidget(self.stick_triple_radius_label)
        form_layout.addRow("Triple Bond Thickness (Stick):", stick_triple_radius_layout)

        # --- 区切り線（水平ライン） ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        form_layout.addRow(line)
        
        # 解像度
        self.stick_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_resolution_slider.setRange(6, 32)
        self.stick_resolution_label = QLabel("16")
        self.stick_resolution_slider.valueChanged.connect(lambda v: self.stick_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.stick_resolution_slider)
        resolution_layout.addWidget(self.stick_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)




        self.tab_widget.addTab(stick_widget, "Stick")

    def reset_current_tab(self):
        """現在選択されているタブの設定のみをデフォルトに戻す"""
        current_tab_index = self.tab_widget.currentIndex()
        tab_name = self.tab_widget.tabText(current_tab_index)
        
        # 各タブの設定項目を定義
        # Each tab settings
        # Note: tab labels must match those added to the QTabWidget ("2D Settings", "Scene", "Ball & Stick", ...
        # "CPK (Space-filling)", "Wireframe", "Stick", "Other"). Use the per-model
        # multi-bond keys present in self.default_settings.
        tab_settings = {
            "2D Settings": {
                'bond_width_2d': self.default_settings['bond_width_2d'],
                'bond_spacing_double_2d': self.default_settings['bond_spacing_double_2d'],
                'bond_spacing_triple_2d': self.default_settings['bond_spacing_triple_2d'],
                'atom_font_size_2d': self.default_settings['atom_font_size_2d'],
                'background_color_2d': self.default_settings['background_color_2d'],
                'bond_color_2d': self.default_settings['bond_color_2d'],
                'atom_use_bond_color_2d': self.default_settings['atom_use_bond_color_2d'],
                'bond_cap_style_2d': self.default_settings['bond_cap_style_2d']
            },
            "3D Scene": {
                'background_color': self.default_settings['background_color'],
                'projection_mode': self.default_settings['projection_mode'],
                'show_3d_axes': self.default_settings['show_3d_axes'],
                'lighting_enabled': self.default_settings['lighting_enabled'],
                'light_intensity': self.default_settings['light_intensity'],
                'specular': self.default_settings['specular'],
                'specular_power': self.default_settings['specular_power']
            },
            "Other": {
                # other options
                'skip_chemistry_checks': self.default_settings.get('skip_chemistry_checks', False),
                'display_kekule_3d': self.default_settings.get('display_kekule_3d', False),
                'always_ask_charge': self.default_settings.get('always_ask_charge', False),
                'display_aromatic_circles_3d': self.default_settings.get('display_aromatic_circles_3d', False),
                'aromatic_torus_thickness_factor': self.default_settings.get('aromatic_torus_thickness_factor', 0.6),
            },
            "Ball & Stick": {
                'ball_stick_atom_scale': self.default_settings['ball_stick_atom_scale'],
                'ball_stick_bond_radius': self.default_settings['ball_stick_bond_radius'],
                'ball_stick_resolution': self.default_settings['ball_stick_resolution'],
                'ball_stick_double_bond_offset_factor': self.default_settings.get('ball_stick_double_bond_offset_factor', 2.0),
                'ball_stick_triple_bond_offset_factor': self.default_settings.get('ball_stick_triple_bond_offset_factor', 2.0),
                'ball_stick_double_bond_radius_factor': self.default_settings.get('ball_stick_double_bond_radius_factor', 0.8),
                'ball_stick_triple_bond_radius_factor': self.default_settings.get('ball_stick_triple_bond_radius_factor', 0.75),
                'ball_stick_use_cpk_bond_color': self.default_settings['ball_stick_use_cpk_bond_color'],
                'ball_stick_bond_color': self.default_settings.get('ball_stick_bond_color', '#7F7F7F')
            },
            "CPK (Space-filling)": {
                'cpk_atom_scale': self.default_settings['cpk_atom_scale'],
                'cpk_resolution': self.default_settings['cpk_resolution'],
                'cpk_colors': {}
            },
            "Wireframe": {
                'wireframe_bond_radius': self.default_settings['wireframe_bond_radius'],
                'wireframe_resolution': self.default_settings['wireframe_resolution'],
                'wireframe_double_bond_offset_factor': self.default_settings.get('wireframe_double_bond_offset_factor', 3.0),
                'wireframe_triple_bond_offset_factor': self.default_settings.get('wireframe_triple_bond_offset_factor', 3.0),
                'wireframe_double_bond_radius_factor': self.default_settings.get('wireframe_double_bond_radius_factor', 0.8),
                'wireframe_triple_bond_radius_factor': self.default_settings.get('wireframe_triple_bond_radius_factor', 0.75)
            },
            "Stick": {
                'stick_bond_radius': self.default_settings['stick_bond_radius'],
                'stick_resolution': self.default_settings['stick_resolution'],
                'stick_double_bond_offset_factor': self.default_settings.get('stick_double_bond_offset_factor', 1.5),
                'stick_triple_bond_offset_factor': self.default_settings.get('stick_triple_bond_offset_factor', 1.0),
                'stick_double_bond_radius_factor': self.default_settings.get('stick_double_bond_radius_factor', 0.6),
                'stick_triple_bond_radius_factor': self.default_settings.get('stick_triple_bond_radius_factor', 0.4)
            }
        }
        
        # 選択されたタブの設定のみを適用
        if tab_name in tab_settings:
            tab_defaults = tab_settings[tab_name]
            
            # 現在の設定を取得
            current_settings = self.get_current_ui_settings()
            
            # 選択されたタブの項目のみをデフォルト値で更新
            updated_settings = current_settings.copy()
            updated_settings.update(tab_defaults)
            
            # UIを更新
            self.update_ui_from_settings(updated_settings)
            
            # If 2D settings were reset, update internal color variables too
            if tab_name == "2D Settings":
                self.current_bg_color_2d = updated_settings.get('background_color_2d', self.default_settings['background_color_2d'])
                self.current_bond_color_2d = updated_settings.get('bond_color_2d', self.default_settings['bond_color_2d'])
                self.update_2d_color_buttons()

            # CPK tab: do not change parent/settings immediately; let Apply/OK persist any changes
            
            # ユーザーへのフィードバック
            QMessageBox.information(self, "Reset Complete", f"Settings for '{tab_name}' tab have been reset to defaults.")
        else:
            QMessageBox.warning(self, "Error", f"Unknown tab: {tab_name}")
    
    def reset_all_settings(self):
        """すべての設定をデフォルトに戻す"""
        reply = QMessageBox.question(
            self, 
            "Reset All Settings", 
            "Are you sure you want to reset all settings to defaults?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            # Update the dialog UI
            self.update_ui_from_settings(self.default_settings)

            # Also persist defaults to the application-level settings if parent is available
            try:
                if self.parent_window and hasattr(self.parent_window, 'settings'):
                    # Update parent settings and save
                    self.parent_window.settings.update(self.default_settings)
                    # defer writing to disk; mark dirty so closeEvent will persist
                    try:
                        self.parent_window.settings_dirty = True
                    except Exception:
                        pass
                    # Also ensure color settings return to defaults and UI reflects them
                    try:
                        # Remove any CPK overrides to restore defaults
                        if 'cpk_colors' in self.parent_window.settings:
                            # Reset to defaults (empty override dict)
                            self.parent_window.settings['cpk_colors'] = {}
                        
                        # Reset 2D Settings (colors need manual refresh in parent scene/view usually handled by apply)
                        
                        # Reset 3D Ball & Stick uniform bond color to default
                        self.parent_window.settings['ball_stick_bond_color'] = self.default_settings.get('ball_stick_bond_color', '#7F7F7F')
                        # Update global CPK colors and reapply 3D settings immediately
                        try:
                            self.parent_window.update_cpk_colors_from_settings()
                        except Exception:
                            pass
                        try:
                            self.parent_window.apply_3d_settings()
                        except Exception:
                            pass
                        # Re-draw current 3D molecule if any
                        try:
                            if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                                self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
                        except Exception:
                            pass
                        # Update 2D scene items to reflect color reset
                        try:
                            if hasattr(self.parent_window, 'scene'):
                                for it in self.parent_window.scene.items():
                                    try:
                                        if hasattr(it, 'update_style'):
                                            it.update_style()
                                    except Exception:
                                        pass
                        except Exception:
                            pass
                        # Mark settings dirty so they'll be saved on exit
                        try:
                            self.parent_window.settings_dirty = True
                        except Exception:
                            pass
                    except Exception:
                        pass

                    # Refresh parent's optimization and conversion menu/action states
                    try:
                        # Optimization method
                        if hasattr(self.parent_window, 'optimization_method'):
                            self.parent_window.optimization_method = self.parent_window.settings.get('optimization_method', 'MMFF_RDKIT')
                        if hasattr(self.parent_window, 'opt3d_actions'):
                            for k, act in self.parent_window.opt3d_actions.items():
                                try:
                                    act.setChecked(k.upper() == (self.parent_window.optimization_method or '').upper())
                                except Exception:
                                    pass

                        # Conversion mode
                        conv_mode = self.parent_window.settings.get('3d_conversion_mode', 'fallback')
                        if hasattr(self.parent_window, 'conv_actions'):
                            for k, act in self.parent_window.conv_actions.items():
                                try:
                                    act.setChecked(k == conv_mode)
                                except Exception:
                                    pass
                    except Exception:
                        pass
            except Exception:
                pass

            QMessageBox.information(self, "Reset Complete", "All settings have been reset to defaults.")

    def get_current_ui_settings(self):
        """現在のUIから設定値を取得"""
        return {
            'background_color': self.current_bg_color,
            'show_3d_axes': self.axes_checkbox.isChecked(),
            'lighting_enabled': self.light_checkbox.isChecked(),
            'light_intensity': self.intensity_slider.value() / 100.0,
            'specular': self.specular_slider.value() / 100.0,
            'specular_power': self.spec_power_slider.value(),
            # Ball and Stick settings
            'ball_stick_atom_scale': self.bs_atom_scale_slider.value() / 100.0,
            'ball_stick_bond_radius': self.bs_bond_radius_slider.value() / 100.0,
            'ball_stick_resolution': self.bs_resolution_slider.value(),
            # CPK settings
            'cpk_atom_scale': self.cpk_atom_scale_slider.value() / 100.0,
            'cpk_resolution': self.cpk_resolution_slider.value(),
            # Wireframe settings
            'wireframe_bond_radius': self.wf_bond_radius_slider.value() / 100.0,
            'wireframe_resolution': self.wf_resolution_slider.value(),
            # Stick settings
            'stick_bond_radius': self.stick_bond_radius_slider.value() / 100.0,
            'stick_resolution': self.stick_resolution_slider.value(),
            # Multi-bond settings (per-model)
            'ball_stick_double_bond_offset_factor': self.bs_double_offset_slider.value() / 100.0,
            'ball_stick_triple_bond_offset_factor': self.bs_triple_offset_slider.value() / 100.0,
            'ball_stick_double_bond_radius_factor': self.bs_double_radius_slider.value() / 100.0,
            'ball_stick_triple_bond_radius_factor': self.bs_triple_radius_slider.value() / 100.0,
            'wireframe_double_bond_offset_factor': self.wf_double_offset_slider.value() / 100.0,
            'wireframe_triple_bond_offset_factor': self.wf_triple_offset_slider.value() / 100.0,
            'wireframe_double_bond_radius_factor': self.wf_double_radius_slider.value() / 100.0,
            'wireframe_triple_bond_radius_factor': self.wf_triple_radius_slider.value() / 100.0,
            'stick_double_bond_offset_factor': self.stick_double_offset_slider.value() / 100.0,
            'stick_triple_bond_offset_factor': self.stick_triple_offset_slider.value() / 100.0,
            'stick_double_bond_radius_factor': self.stick_double_radius_slider.value() / 100.0,
            'stick_triple_bond_radius_factor': self.stick_triple_radius_slider.value() / 100.0,
            # Projection mode
            'projection_mode': self.projection_combo.currentText(),
            # Kekule / aromatic / torus settings
            'display_kekule_3d': self.kekule_3d_checkbox.isChecked(),
            'display_aromatic_circles_3d': self.aromatic_circle_checkbox.isChecked(),
            'aromatic_torus_thickness_factor': self.aromatic_torus_thickness_slider.value() / 100.0,
            'skip_chemistry_checks': self.skip_chem_checks_checkbox.isChecked(),
            'always_ask_charge': self.always_ask_charge_checkbox.isChecked(),
            # Ball & Stick bond color and use-cpk option
            'ball_stick_bond_color': getattr(self, 'bs_bond_color', self.default_settings.get('ball_stick_bond_color')),
            'ball_stick_use_cpk_bond_color': self.bs_use_cpk_bond_checkbox.isChecked(),
        }
    
    def reset_to_defaults(self):
        """UIをデフォルト設定に戻す（後方互換性のため残存）"""
        self.reset_all_settings()

    def update_ui_from_settings(self, settings_dict):
        # 基本設定
        self.current_bg_color = settings_dict.get('background_color', self.default_settings['background_color'])
        self.update_color_button(self.current_bg_color)
        self.axes_checkbox.setChecked(settings_dict.get('show_3d_axes', self.default_settings['show_3d_axes']))
        self.light_checkbox.setChecked(settings_dict.get('lighting_enabled', self.default_settings['lighting_enabled']))
        
        # スライダーの値を設定
        intensity_val = int(settings_dict.get('light_intensity', self.default_settings['light_intensity']) * 100)
        self.intensity_slider.setValue(intensity_val)
        self.intensity_label.setText(f"{intensity_val/100:.2f}")
        
        specular_val = int(settings_dict.get('specular', self.default_settings['specular']) * 100)
        self.specular_slider.setValue(specular_val)
        self.specular_label.setText(f"{specular_val/100:.2f}")
        
        self.spec_power_slider.setValue(settings_dict.get('specular_power', self.default_settings['specular_power']))
        self.spec_power_label.setText(str(settings_dict.get('specular_power', self.default_settings['specular_power'])))
        
        # Ball and Stick設定
        bs_atom_scale = int(settings_dict.get('ball_stick_atom_scale', self.default_settings['ball_stick_atom_scale']) * 100)
        self.bs_atom_scale_slider.setValue(bs_atom_scale)
        self.bs_atom_scale_label.setText(f"{bs_atom_scale/100:.2f}")
        
        bs_bond_radius = int(settings_dict.get('ball_stick_bond_radius', self.default_settings['ball_stick_bond_radius']) * 100)
        self.bs_bond_radius_slider.setValue(bs_bond_radius)
        self.bs_bond_radius_label.setText(f"{bs_bond_radius/100:.2f}")
        
        self.bs_resolution_slider.setValue(settings_dict.get('ball_stick_resolution', self.default_settings['ball_stick_resolution']))
        self.bs_resolution_label.setText(str(settings_dict.get('ball_stick_resolution', self.default_settings['ball_stick_resolution'])))
        # Ball & Stick bond color (uniform gray color for ball-and-stick)
        bs_bond_color = settings_dict.get('ball_stick_bond_color', self.default_settings.get('ball_stick_bond_color', '#7F7F7F'))
        try:
            self.bs_bond_color = QColor(bs_bond_color).name()
        except Exception:
            self.bs_bond_color = self.default_settings.get('ball_stick_bond_color', '#7F7F7F')
        # Ensure color button exists and update its appearance
        try:
            if hasattr(self, 'bs_bond_color_button'):
                self.bs_bond_color_button.setStyleSheet(f"background-color: {self.bs_bond_color}; border: 1px solid #888;")
                try:
                    self.bs_bond_color_button.setToolTip(self.bs_bond_color)
                except Exception:
                    pass
        except Exception:
            pass
        
        # Ball & Stick CPK bond color option
        self.bs_use_cpk_bond_checkbox.setChecked(settings_dict.get('ball_stick_use_cpk_bond_color', self.default_settings['ball_stick_use_cpk_bond_color']))
        
        # CPK設定
        cpk_atom_scale = int(settings_dict.get('cpk_atom_scale', self.default_settings['cpk_atom_scale']) * 100)
        self.cpk_atom_scale_slider.setValue(cpk_atom_scale)
        self.cpk_atom_scale_label.setText(f"{cpk_atom_scale/100:.2f}")
        
        self.cpk_resolution_slider.setValue(settings_dict.get('cpk_resolution', self.default_settings['cpk_resolution']))
        self.cpk_resolution_label.setText(str(settings_dict.get('cpk_resolution', self.default_settings['cpk_resolution'])))
        
        # Wireframe設定
        wf_bond_radius = int(settings_dict.get('wireframe_bond_radius', self.default_settings['wireframe_bond_radius']) * 100)
        self.wf_bond_radius_slider.setValue(wf_bond_radius)
        self.wf_bond_radius_label.setText(f"{wf_bond_radius/100:.2f}")
        
        self.wf_resolution_slider.setValue(settings_dict.get('wireframe_resolution', self.default_settings['wireframe_resolution']))
        self.wf_resolution_label.setText(str(settings_dict.get('wireframe_resolution', self.default_settings['wireframe_resolution'])))
        
        # Stick設定
        stick_bond_radius = int(settings_dict.get('stick_bond_radius', self.default_settings['stick_bond_radius']) * 100)
        self.stick_bond_radius_slider.setValue(stick_bond_radius)
        self.stick_bond_radius_label.setText(f"{stick_bond_radius/100:.2f}")
        
        self.stick_resolution_slider.setValue(settings_dict.get('stick_resolution', self.default_settings['stick_resolution']))
        self.stick_resolution_label.setText(str(settings_dict.get('stick_resolution', self.default_settings['stick_resolution'])))
        
        # 多重結合設定（モデル毎） — 後方互換のため既存のグローバルキーがあればフォールバック
        # Ball & Stick
        bs_double_offset = int(settings_dict.get('ball_stick_double_bond_offset_factor', settings_dict.get('double_bond_offset_factor', self.default_settings.get('ball_stick_double_bond_offset_factor', 2.0))) * 100)
        self.bs_double_offset_slider.setValue(bs_double_offset)
        self.bs_double_offset_label.setText(f"{bs_double_offset/100:.2f}")
    
        bs_triple_offset = int(settings_dict.get('ball_stick_triple_bond_offset_factor', settings_dict.get('triple_bond_offset_factor', self.default_settings.get('ball_stick_triple_bond_offset_factor', 2.0))) * 100)
        self.bs_triple_offset_slider.setValue(bs_triple_offset)
        self.bs_triple_offset_label.setText(f"{bs_triple_offset/100:.2f}")
    
        bs_double_radius = int(settings_dict.get('ball_stick_double_bond_radius_factor', settings_dict.get('double_bond_radius_factor', self.default_settings.get('ball_stick_double_bond_radius_factor', 1.0))) * 100)
        self.bs_double_radius_slider.setValue(bs_double_radius)
        self.bs_double_radius_label.setText(f"{bs_double_radius/100:.2f}")
    
        bs_triple_radius = int(settings_dict.get('ball_stick_triple_bond_radius_factor', settings_dict.get('triple_bond_radius_factor', self.default_settings.get('ball_stick_triple_bond_radius_factor', 0.75))) * 100)
        self.bs_triple_radius_slider.setValue(bs_triple_radius)
        self.bs_triple_radius_label.setText(f"{bs_triple_radius/100:.2f}")
    
        # Wireframe
        wf_double_offset = int(settings_dict.get('wireframe_double_bond_offset_factor', settings_dict.get('double_bond_offset_factor', self.default_settings.get('wireframe_double_bond_offset_factor', 3.0))) * 100)
        self.wf_double_offset_slider.setValue(wf_double_offset)
        self.wf_double_offset_label.setText(f"{wf_double_offset/100:.2f}")
    
        wf_triple_offset = int(settings_dict.get('wireframe_triple_bond_offset_factor', settings_dict.get('triple_bond_offset_factor', self.default_settings.get('wireframe_triple_bond_offset_factor', 3.0))) * 100)
        self.wf_triple_offset_slider.setValue(wf_triple_offset)
        self.wf_triple_offset_label.setText(f"{wf_triple_offset/100:.2f}")
    
        wf_double_radius = int(settings_dict.get('wireframe_double_bond_radius_factor', settings_dict.get('double_bond_radius_factor', self.default_settings.get('wireframe_double_bond_radius_factor', 1.0))) * 100)
        self.wf_double_radius_slider.setValue(wf_double_radius)
        self.wf_double_radius_label.setText(f"{wf_double_radius/100:.2f}")
    
        wf_triple_radius = int(settings_dict.get('wireframe_triple_bond_radius_factor', settings_dict.get('triple_bond_radius_factor', self.default_settings.get('wireframe_triple_bond_radius_factor', 0.75))) * 100)
        self.wf_triple_radius_slider.setValue(wf_triple_radius)
        self.wf_triple_radius_label.setText(f"{wf_triple_radius/100:.2f}")
    
        # Stick
        stick_double_offset = int(settings_dict.get('stick_double_bond_offset_factor', settings_dict.get('double_bond_offset_factor', self.default_settings.get('stick_double_bond_offset_factor', 1.5))) * 100)
        self.stick_double_offset_slider.setValue(stick_double_offset)
        self.stick_double_offset_label.setText(f"{stick_double_offset/100:.2f}")
    
        stick_triple_offset = int(settings_dict.get('stick_triple_bond_offset_factor', settings_dict.get('triple_bond_offset_factor', self.default_settings.get('stick_triple_bond_offset_factor', 1.0))) * 100)
        self.stick_triple_offset_slider.setValue(stick_triple_offset)
        self.stick_triple_offset_label.setText(f"{stick_triple_offset/100:.2f}")
    
        stick_double_radius = int(settings_dict.get('stick_double_bond_radius_factor', settings_dict.get('double_bond_radius_factor', self.default_settings.get('stick_double_bond_radius_factor', 0.60))) * 100)
        self.stick_double_radius_slider.setValue(stick_double_radius)
        self.stick_double_radius_label.setText(f"{stick_double_radius/100:.2f}")
    
        stick_triple_radius = int(settings_dict.get('stick_triple_bond_radius_factor', settings_dict.get('triple_bond_radius_factor', self.default_settings.get('stick_triple_bond_radius_factor', 0.40))) * 100)
        self.stick_triple_radius_slider.setValue(stick_triple_radius)
        self.stick_triple_radius_label.setText(f"{stick_triple_radius/100:.2f}")
        
        # Projection mode
        proj_mode = settings_dict.get('projection_mode', self.default_settings.get('projection_mode', 'Perspective'))
        idx = self.projection_combo.findText(proj_mode)
        self.projection_combo.setCurrentIndex(idx if idx != -1 else 0)
        # skip chemistry checks
        self.skip_chem_checks_checkbox.setChecked(settings_dict.get('skip_chemistry_checks', self.default_settings.get('skip_chemistry_checks', False)))
        # kekule setting
        self.kekule_3d_checkbox.setChecked(settings_dict.get('display_kekule_3d', self.default_settings.get('display_kekule_3d', False)))
        # always ask for charge on XYZ imports
        self.always_ask_charge_checkbox.setChecked(settings_dict.get('always_ask_charge', self.default_settings.get('always_ask_charge', False)))
        # Aromatic ring circle display and torus thickness factor
        self.aromatic_circle_checkbox.setChecked(settings_dict.get('display_aromatic_circles_3d', self.default_settings.get('display_aromatic_circles_3d', False)))
        thickness_factor = float(settings_dict.get('aromatic_torus_thickness_factor', self.default_settings.get('aromatic_torus_thickness_factor', 0.6)))
        try:
            self.aromatic_torus_thickness_slider.setValue(int(thickness_factor * 100))
            self.aromatic_torus_thickness_label.setText(f"{thickness_factor:.1f}")
        except Exception:
            pass
      
    def select_color(self):
        """カラーピッカーを開き、選択された色を内部変数とUIに反映させる"""
        # 内部変数から現在の色を取得してカラーピッカーを初期化
        color = QColorDialog.getColor(QColor(self.current_bg_color), self)
        if color.isValid():
            # 内部変数を更新
            self.current_bg_color = color.name()
            # UIの見た目を更新
            self.update_color_button(self.current_bg_color)

    def update_color_button(self, color_hex):
        """ボタンの背景色と境界線を設定する"""
        self.bg_button.setStyleSheet(f"background-color: {color_hex}; border: 1px solid #888;")

    def get_settings(self):
        return {
            'background_color': self.current_bg_color,
            'projection_mode': self.projection_combo.currentText(),
            'show_3d_axes': self.axes_checkbox.isChecked(),
            'lighting_enabled': self.light_checkbox.isChecked(),
            'light_intensity': self.intensity_slider.value() / 100.0,
            'specular': self.specular_slider.value() / 100.0,
            'specular_power': self.spec_power_slider.value(),
            # Ball and Stick settings
            'ball_stick_atom_scale': self.bs_atom_scale_slider.value() / 100.0,
            'ball_stick_bond_radius': self.bs_bond_radius_slider.value() / 100.0,
            'ball_stick_resolution': self.bs_resolution_slider.value(),
            # CPK settings
            'cpk_atom_scale': self.cpk_atom_scale_slider.value() / 100.0,
            'cpk_resolution': self.cpk_resolution_slider.value(),
            # Wireframe settings
            'wireframe_bond_radius': self.wf_bond_radius_slider.value() / 100.0,
            'wireframe_resolution': self.wf_resolution_slider.value(),
            # Stick settings
            'stick_bond_radius': self.stick_bond_radius_slider.value() / 100.0,
            'stick_resolution': self.stick_resolution_slider.value(),
            # Multiple bond offset settings (per-model)
            'ball_stick_double_bond_offset_factor': self.bs_double_offset_slider.value() / 100.0,
            'ball_stick_triple_bond_offset_factor': self.bs_triple_offset_slider.value() / 100.0,
            'ball_stick_double_bond_radius_factor': self.bs_double_radius_slider.value() / 100.0,
            'ball_stick_triple_bond_radius_factor': self.bs_triple_radius_slider.value() / 100.0,
            'wireframe_double_bond_offset_factor': self.wf_double_offset_slider.value() / 100.0,
            'wireframe_triple_bond_offset_factor': self.wf_triple_offset_slider.value() / 100.0,
            'wireframe_double_bond_radius_factor': self.wf_double_radius_slider.value() / 100.0,
            'wireframe_triple_bond_radius_factor': self.wf_triple_radius_slider.value() / 100.0,
            'stick_double_bond_offset_factor': self.stick_double_offset_slider.value() / 100.0,
            'stick_triple_bond_offset_factor': self.stick_triple_offset_slider.value() / 100.0,
            'stick_double_bond_radius_factor': self.stick_double_radius_slider.value() / 100.0,
            'stick_triple_bond_radius_factor': self.stick_triple_radius_slider.value() / 100.0,
            'display_kekule_3d': self.kekule_3d_checkbox.isChecked(),
            'display_aromatic_circles_3d': self.aromatic_circle_checkbox.isChecked(),
            'aromatic_torus_thickness_factor': self.aromatic_torus_thickness_slider.value() / 100.0,
            'skip_chemistry_checks': self.skip_chem_checks_checkbox.isChecked(),
            'always_ask_charge': self.always_ask_charge_checkbox.isChecked(),
            'ball_stick_bond_color': getattr(self, 'bs_bond_color', self.default_settings.get('ball_stick_bond_color', '#7F7F7F')),
            'ball_stick_use_cpk_bond_color': self.bs_use_cpk_bond_checkbox.isChecked(),
            # 2D Settings
            'bond_width_2d': self.bond_width_2d_slider.value() / 10.0,
            'bond_spacing_double_2d': self.bond_spacing_double_2d_slider.value() / 10.0,
            'bond_spacing_triple_2d': self.bond_spacing_triple_2d_slider.value() / 10.0,
            'atom_font_size_2d': self.atom_font_size_2d_slider.value(),
            'background_color_2d': self.current_bg_color_2d,
            'bond_color_2d': self.current_bond_color_2d,
            'atom_use_bond_color_2d': self.atom_use_bond_color_2d_checkbox.isChecked(),
            'bond_cap_style_2d': self.bond_cap_style_2d_combo.currentText(),
        }

    def pick_bs_bond_color(self):
        """Open QColorDialog to pick Ball & Stick bond color (3D)."""
        cur = getattr(self, 'bs_bond_color', self.default_settings.get('ball_stick_bond_color', '#7F7F7F'))
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.bs_bond_color = color.name()
            try:
                self.bs_bond_color_button.setStyleSheet(f"background-color: {self.bs_bond_color}; border: 1px solid #888;")
                self.bs_bond_color_button.setToolTip(self.bs_bond_color)
            except Exception:
                pass

    def apply_settings(self):
        """設定を適用（ダイアログは開いたまま）"""
        # 親ウィンドウの設定を更新
        if self.parent_window:
            settings = self.get_settings()
            self.parent_window.settings.update(settings)
            # Mark settings dirty; persist on exit to avoid frequent disk writes
            try:
                self.parent_window.settings_dirty = True
            except Exception:
                pass
            # 3Dビューの設定を適用
            self.parent_window.apply_3d_settings()
            # Update CPK colors from settings if present (no-op otherwise)
            try:
                self.parent_window.update_cpk_colors_from_settings()
            except Exception:
                pass
            # Refresh any open CPK color dialogs so they update their UI
            try:
                for w in QApplication.topLevelWidgets():
                    try:
                        # import locally to avoid circular import
                        try:
                            from .color_settings_dialog import ColorSettingsDialog
                        except Exception:
                            from modules.color_settings_dialog import ColorSettingsDialog
                        if isinstance(w, ColorSettingsDialog):
                            try:
                                w.refresh_ui()
                            except Exception:
                                pass
                    except Exception:
                        pass
            except Exception:
                pass
            # 現在の分子を再描画（設定変更を反映）
            if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
            
            # 2Dビューの設定適用 (背景色、結合スタイルなど)
            # update_style() will read the new settings from parent_window.settings
            try:
                if hasattr(self.parent_window, 'scene') and self.parent_window.scene:
                    # Update Background
                    bg_col_2d = self.parent_window.settings.get('background_color_2d', '#FFFFFF')
                    self.parent_window.scene.setBackgroundBrush(QColor(bg_col_2d))
                    
                    # Update Items
                    for item in self.parent_window.scene.items():
                        if hasattr(item, 'update_style'):
                            item.update_style()
                        elif hasattr(item, 'update'):
                            item.update()
                            
                    if hasattr(self.parent_window.view_2d, 'viewport'):
                        self.parent_window.view_2d.viewport().update()
            except Exception:
                pass
            
            # ステータスバーに適用完了を表示
            self.parent_window.statusBar().showMessage("Settings applied successfully")

    def _on_skip_chem_checks_changed(self, state):
        """Handle user toggling of skip chemistry checks: persist and update UI.

        state: Qt.Checked (2) or Qt.Unchecked (0)
        """
        try:
            enabled = bool(state)
            self.settings['skip_chemistry_checks'] = enabled
            # mark dirty instead of immediate save
            try:
                self.settings_dirty = True
            except Exception:
                pass
            # If skip is enabled, allow Optimize button; otherwise, respect chem_check flags

        except Exception:
            pass

    # Note: Kekule display is applied only when user clicks Apply/OK.

    def accept(self):
        """ダイアログの設定を適用してから閉じる"""
        # apply_settingsを呼び出して設定を適用
        self.apply_settings()
        super().accept()

    def update_ui_from_settings(self, settings_dict):
        """設定辞書に基づいてUIを更新"""
        # 1. Scene settings
        self.current_bg_color = settings_dict.get('background_color', self.default_settings['background_color'])
        self.update_color_button(self.current_bg_color)
        
        self.axes_checkbox.setChecked(settings_dict.get('show_3d_axes', self.default_settings['show_3d_axes']))
        self.light_checkbox.setChecked(settings_dict.get('lighting_enabled', self.default_settings['lighting_enabled']))
        
        intensity = settings_dict.get('light_intensity', self.default_settings['light_intensity'])
        self.intensity_slider.setValue(int(intensity * 100))
        self.intensity_label.setText(f"{intensity:.2f}")
        
        specular = settings_dict.get('specular', self.default_settings['specular'])
        self.specular_slider.setValue(int(specular * 100))
        self.specular_label.setText(f"{specular:.2f}")
        
        spec_power = settings_dict.get('specular_power', self.default_settings['specular_power'])
        self.spec_power_slider.setValue(int(spec_power))
        self.spec_power_label.setText(str(spec_power))
        
        proj_mode = settings_dict.get('projection_mode', self.default_settings['projection_mode'])
        idx = self.projection_combo.findText(proj_mode)
        if idx >= 0:
            self.projection_combo.setCurrentIndex(idx)
            
        # 2. Ball & Stick settings
        bs_atom_scale = settings_dict.get('ball_stick_atom_scale', self.default_settings['ball_stick_atom_scale'])
        self.bs_atom_scale_slider.setValue(int(bs_atom_scale * 100))
        self.bs_atom_scale_label.setText(f"{bs_atom_scale:.2f}")
        
        bs_bond_radius = settings_dict.get('ball_stick_bond_radius', self.default_settings['ball_stick_bond_radius'])
        self.bs_bond_radius_slider.setValue(int(bs_bond_radius * 100))
        self.bs_bond_radius_label.setText(f"{bs_bond_radius:.2f}")

        # Multi-bond offsets (Ball & Stick)
        bs_db_offset = settings_dict.get('ball_stick_double_bond_offset_factor', self.default_settings.get('ball_stick_double_bond_offset_factor', 2.0))
        self.bs_double_offset_slider.setValue(int(bs_db_offset * 100))
        self.bs_double_offset_label.setText(f"{bs_db_offset:.2f}")

        bs_tr_offset = settings_dict.get('ball_stick_triple_bond_offset_factor', self.default_settings.get('ball_stick_triple_bond_offset_factor', 2.0))
        self.bs_triple_offset_slider.setValue(int(bs_tr_offset * 100))
        self.bs_triple_offset_label.setText(f"{bs_tr_offset:.2f}")
        
        bs_db_rad = settings_dict.get('ball_stick_double_bond_radius_factor', self.default_settings.get('ball_stick_double_bond_radius_factor', 0.8))
        self.bs_double_radius_slider.setValue(int(bs_db_rad * 100))
        self.bs_double_radius_label.setText(f"{bs_db_rad:.2f}")

        bs_tr_rad = settings_dict.get('ball_stick_triple_bond_radius_factor', self.default_settings.get('ball_stick_triple_bond_radius_factor', 0.75))
        self.bs_triple_radius_slider.setValue(int(bs_tr_rad * 100))
        self.bs_triple_radius_label.setText(f"{bs_tr_rad:.2f}")

        bs_res = settings_dict.get('ball_stick_resolution', self.default_settings['ball_stick_resolution'])
        self.bs_resolution_slider.setValue(int(bs_res))
        self.bs_resolution_label.setText(str(bs_res))

        # Ball & Stick bond color
        self.bs_bond_color = settings_dict.get('ball_stick_bond_color', self.default_settings.get('ball_stick_bond_color', '#7F7F7F'))
        try:
            self.bs_bond_color_button.setStyleSheet(f"background-color: {self.bs_bond_color}; border: 1px solid #888;")
            self.bs_bond_color_button.setToolTip(self.bs_bond_color)
        except Exception:
            pass

        self.bs_use_cpk_bond_checkbox.setChecked(settings_dict.get('ball_stick_use_cpk_bond_color', self.default_settings.get('ball_stick_use_cpk_bond_color', False)))
            
        # 3. CPK settings
        cpk_atom_scale = settings_dict.get('cpk_atom_scale', self.default_settings['cpk_atom_scale'])
        self.cpk_atom_scale_slider.setValue(int(cpk_atom_scale * 100))
        self.cpk_atom_scale_label.setText(f"{cpk_atom_scale:.2f}")
        
        cpk_res = settings_dict.get('cpk_resolution', self.default_settings['cpk_resolution'])
        self.cpk_resolution_slider.setValue(int(cpk_res))
        self.cpk_resolution_label.setText(str(cpk_res))
        
        # 4. Wireframe settings
        wf_bond_radius = settings_dict.get('wireframe_bond_radius', self.default_settings['wireframe_bond_radius'])
        self.wf_bond_radius_slider.setValue(int(wf_bond_radius * 100))
        self.wf_bond_radius_label.setText(f"{wf_bond_radius:.2f}")
        
        wf_res = settings_dict.get('wireframe_resolution', self.default_settings['wireframe_resolution'])
        self.wf_resolution_slider.setValue(int(wf_res))
        self.wf_resolution_label.setText(str(wf_res))

        # Multi-bond offsets (Wireframe)
        wf_db_offset = settings_dict.get('wireframe_double_bond_offset_factor', self.default_settings.get('wireframe_double_bond_offset_factor', 3.0))
        self.wf_double_offset_slider.setValue(int(wf_db_offset * 100))
        self.wf_double_offset_label.setText(f"{wf_db_offset:.2f}")

        wf_tr_offset = settings_dict.get('wireframe_triple_bond_offset_factor', self.default_settings.get('wireframe_triple_bond_offset_factor', 3.0))
        self.wf_triple_offset_slider.setValue(int(wf_tr_offset * 100))
        self.wf_triple_offset_label.setText(f"{wf_tr_offset:.2f}")
        
        wf_db_rad = settings_dict.get('wireframe_double_bond_radius_factor', self.default_settings.get('wireframe_double_bond_radius_factor', 0.8))
        self.wf_double_radius_slider.setValue(int(wf_db_rad * 100))
        self.wf_double_radius_label.setText(f"{wf_db_rad:.2f}")

        wf_tr_rad = settings_dict.get('wireframe_triple_bond_radius_factor', self.default_settings.get('wireframe_triple_bond_radius_factor', 0.75))
        self.wf_triple_radius_slider.setValue(int(wf_tr_rad * 100))
        self.wf_triple_radius_label.setText(f"{wf_tr_rad:.2f}")
            
        # 5. Stick settings
        stick_bond_radius = settings_dict.get('stick_bond_radius', self.default_settings['stick_bond_radius'])
        self.stick_bond_radius_slider.setValue(int(stick_bond_radius * 100))
        self.stick_bond_radius_label.setText(f"{stick_bond_radius:.2f}")
        
        stick_res = settings_dict.get('stick_resolution', self.default_settings['stick_resolution'])
        self.stick_resolution_slider.setValue(int(stick_res))
        self.stick_resolution_label.setText(str(stick_res))

        # Multi-bond offsets (Stick)
        stick_db_offset = settings_dict.get('stick_double_bond_offset_factor', self.default_settings.get('stick_double_bond_offset_factor', 1.5))
        self.stick_double_offset_slider.setValue(int(stick_db_offset * 100))
        self.stick_double_offset_label.setText(f"{stick_db_offset:.2f}")

        stick_tr_offset = settings_dict.get('stick_triple_bond_offset_factor', self.default_settings.get('stick_triple_bond_offset_factor', 1.0))
        self.stick_triple_offset_slider.setValue(int(stick_tr_offset * 100))
        self.stick_triple_offset_label.setText(f"{stick_tr_offset:.2f}")
        
        stick_db_rad = settings_dict.get('stick_double_bond_radius_factor', self.default_settings.get('stick_double_bond_radius_factor', 0.6))
        self.stick_double_radius_slider.setValue(int(stick_db_rad * 100))
        self.stick_double_radius_label.setText(f"{stick_db_rad:.2f}")

        stick_tr_rad = settings_dict.get('stick_triple_bond_radius_factor', self.default_settings.get('stick_triple_bond_radius_factor', 0.4))
        self.stick_triple_radius_slider.setValue(int(stick_tr_rad * 100))
        self.stick_triple_radius_label.setText(f"{stick_tr_rad:.2f}")

        # 6. Other settings
        self.skip_chem_checks_checkbox.setChecked(settings_dict.get('skip_chemistry_checks', self.default_settings.get('skip_chemistry_checks', False)))
        self.kekule_3d_checkbox.setChecked(settings_dict.get('display_kekule_3d', self.default_settings.get('display_kekule_3d', False)))
        # always ask for charge on XYZ imports
        self.always_ask_charge_checkbox.setChecked(settings_dict.get('always_ask_charge', self.default_settings.get('always_ask_charge', False)))
        # Aromatic ring circle display and torus thickness factor
        self.aromatic_circle_checkbox.setChecked(settings_dict.get('display_aromatic_circles_3d', self.default_settings.get('display_aromatic_circles_3d', False)))
        thickness_factor = float(settings_dict.get('aromatic_torus_thickness_factor', self.default_settings.get('aromatic_torus_thickness_factor', 0.6)))
        try:
            self.aromatic_torus_thickness_slider.setValue(int(thickness_factor * 100))
            self.aromatic_torus_thickness_label.setText(f"{thickness_factor:.1f}")
        except Exception:
            pass

        # 7. 2D Settings
        bw_2d = settings_dict.get('bond_width_2d', self.default_settings['bond_width_2d'])
        self.bond_width_2d_slider.setValue(int(bw_2d * 10))
        self.bond_width_2d_label.setText(f"{bw_2d:.1f}")

        bsd_2d = settings_dict.get('bond_spacing_double_2d', self.default_settings['bond_spacing_double_2d'])
        self.bond_spacing_double_2d_slider.setValue(int(bsd_2d * 10))
        self.bond_spacing_double_2d_label.setText(f"{bsd_2d:.1f}")

        bst_2d = settings_dict.get('bond_spacing_triple_2d', self.default_settings['bond_spacing_triple_2d'])
        self.bond_spacing_triple_2d_slider.setValue(int(bst_2d * 10))
        self.bond_spacing_triple_2d_label.setText(f"{bst_2d:.1f}")

        fs_2d = settings_dict.get('atom_font_size_2d', self.default_settings['atom_font_size_2d'])
        self.atom_font_size_2d_slider.setValue(int(fs_2d))
        self.atom_font_size_2d_label.setText(str(fs_2d))

        self.atom_use_bond_color_2d_checkbox.setChecked(settings_dict.get('atom_use_bond_color_2d', self.default_settings.get('atom_use_bond_color_2d', False)))

        self.current_bg_color_2d = settings_dict.get('background_color_2d', self.default_settings['background_color_2d'])

        # Load Bond Cap Style
        cap_style = settings_dict.get('bond_cap_style_2d', self.default_settings.get('bond_cap_style_2d', 'Round'))
        index = self.bond_cap_style_2d_combo.findText(cap_style)
        if index >= 0:
            self.bond_cap_style_2d_combo.setCurrentIndex(index)
        self.current_bond_color_2d = settings_dict.get('bond_color_2d', self.default_settings['bond_color_2d'])
        self.update_2d_color_buttons()
