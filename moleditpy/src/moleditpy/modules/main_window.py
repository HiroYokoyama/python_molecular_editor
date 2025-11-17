#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: Apache-2.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI 10.5281/zenodo.17268532
"""

import sys
import numpy as np
import pickle
import copy
import math
import io
import os
import ctypes
import itertools
import json 
import vtk
import base64
import contextlib
import re
import traceback

from collections import deque

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QSplitter, QGraphicsView, QGraphicsScene, QGraphicsItem,
    QToolBar, QStatusBar, QGraphicsTextItem, QGraphicsLineItem, QDialog, QGridLayout,
    QFileDialog, QSizePolicy, QLabel, QLineEdit, QToolButton, QMenu, QMessageBox, 
    QInputDialog, QDialogButtonBox, QColorDialog, QCheckBox, QSlider, QFormLayout, 
    QRadioButton, QComboBox, QListWidget, QListWidgetItem, QButtonGroup, QTabWidget, 
    QScrollArea, QFrame, QTableWidget, QTableWidgetItem, QAbstractItemView
)

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QPainterPath, QPainterPathStroker, QFontMetrics, QFontMetricsF, QKeySequence, 
    QTransform, QCursor, QPixmap, QIcon, QShortcut, QDesktopServices, QImage
)


from PyQt6.QtCore import (
    Qt, QPointF, QRectF, QLineF, QObject, QThread, pyqtSignal, pyqtSlot, QEvent, 
    QMimeData, QByteArray, QUrl, QTimer, QDateTime
)

import pyvista as pv

# Open Babel Python binding (optional; required for fallback)
try:
    from openbabel import pybel
    OBABEL_AVAILABLE = True
except Exception:
    # pybel (Open Babel Python bindings) is optional. If not present, disable OBabel features.
    pybel = None
    OBABEL_AVAILABLE = False
    print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")

# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import *
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
    from .template_preview_view import TemplatePreviewView
    from .user_template_dialog import UserTemplateDialog
    from .about_dialog import AboutDialog
    from .translation_dialog import TranslationDialog
    from .mirror_dialog import MirrorDialog
    from .move_group_dialog import MoveGroupDialog
    from .align_plane_dialog import AlignPlaneDialog
    from .planarize_dialog import PlanarizeDialog
    from .alignment_dialog import AlignmentDialog
    from .molecular_data import MolecularData
    from .atom_item import AtomItem
    from .bond_item import BondItem
    from .template_preview_item import TemplatePreviewItem
    from .molecule_scene import MoleculeScene
    from .zoomable_view import ZoomableView
    from .calculation_worker import CalculationWorker
    from .periodic_table_dialog import PeriodicTableDialog
    from .color_settings_dialog import ColorSettingsDialog
    from .analysis_window import AnalysisWindow
    from .settings_dialog import SettingsDialog
    from .custom_qt_interactor import CustomQtInteractor
    from .custom_interactor_style import CustomInteractorStyle
    from .bond_length_dialog import BondLengthDialog
    from .angle_dialog import AngleDialog
    from .dihedral_dialog import DihedralDialog
    from .constrained_optimization_dialog import ConstrainedOptimizationDialog
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import *
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin
    from modules.template_preview_view import TemplatePreviewView
    from modules.user_template_dialog import UserTemplateDialog
    from modules.about_dialog import AboutDialog
    from modules.translation_dialog import TranslationDialog
    from modules.mirror_dialog import MirrorDialog
    from modules.move_group_dialog import MoveGroupDialog
    from modules.align_plane_dialog import AlignPlaneDialog
    from modules.planarize_dialog import PlanarizeDialog
    from modules.alignment_dialog import AlignmentDialog
    from modules.molecular_data import MolecularData
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem
    from modules.template_preview_item import TemplatePreviewItem
    from modules.molecule_scene import MoleculeScene
    from modules.zoomable_view import ZoomableView
    from modules.calculation_worker import CalculationWorker
    from modules.periodic_table_dialog import PeriodicTableDialog
    from modules.color_settings_dialog import ColorSettingsDialog
    from modules.analysis_window import AnalysisWindow
    from modules.settings_dialog import SettingsDialog
    from modules.custom_qt_interactor import CustomQtInteractor
    from modules.custom_interactor_style import CustomInteractorStyle
    from modules.bond_length_dialog import BondLengthDialog
    from modules.angle_dialog import AngleDialog
    from modules.dihedral_dialog import DihedralDialog
    from modules.constrained_optimization_dialog import ConstrainedOptimizationDialog

class MainWindow(QMainWindow):

    # start_calculation carries the MOL block and an options object (second arg)
    start_calculation = pyqtSignal(str, object)
    def __init__(self, initial_file=None):
        super().__init__()
        self.setAcceptDrops(True)
        self.settings_dir = os.path.join(os.path.expanduser('~'), '.moleditpy')
        self.settings_file = os.path.join(self.settings_dir, 'settings.json')
        self.settings = {}
        self.load_settings()
        self.initial_settings = self.settings.copy()
        self.setWindowTitle("MoleditPy Ver. " + VERSION); self.setGeometry(100, 100, 1400, 800)
        self.data = MolecularData(); self.current_mol = None
        self.current_3d_style = 'ball_and_stick'
        self.show_chiral_labels = False
        self.atom_info_display_mode = None  # 'id', 'coords', 'symbol', or None
        self.current_atom_info_labels = None  # 現在の原子情報ラベル
        self.is_3d_edit_mode = False
        self.dragged_atom_info = None
        self.atom_actor = None 
        self.is_2d_editable = True
        self.is_xyz_derived = False  # XYZ由来の分子かどうかのフラグ
        # Chemical check flags: whether a chemical/sanitization check was attempted and whether it failed
        self.chem_check_tried = False
        self.chem_check_failed = False
        # 3D最適化のデフォルト手法
        self.optimization_method = self.settings.get('optimization_method', 'MMFF_RDKIT')
        self.axes_actor = None
        self.axes_widget = None
        self._template_dialog = None  # テンプレートダイアログの参照
        self.undo_stack = []
        self.redo_stack = []
        self.constraints_3d = []
        self.mode_actions = {} 
        
        # 保存状態を追跡する変数
        self.has_unsaved_changes = False
        # 設定ファイルのディスク書き込みを遅延するフラグ
        # True に設定された場合、設定はメモリ上で更新され、アプリ終了時にまとめて保存されます。
        self.settings_dirty = True
        self.current_file_path = None  # 現在開いているファイルのパス
        self.initialization_complete = False  # 初期化完了フラグ
        # Token to invalidate pending implicit-hydrogen UI updates
        self._ih_update_counter = 0
        
        # 測定機能用の変数
        self.measurement_mode = False
        self.selected_atoms_for_measurement = []
        self.measurement_labels = []  # (atom_idx, label_text) のタプルのリスト
        self.measurement_text_actor = None
        self.measurement_label_items_2d = []  # 2Dビューの測定ラベルアイテム
        self.atom_id_to_rdkit_idx_map = {}  # 2D原子IDから3D RDKit原子インデックスへのマッピング
        
        # 3D原子選択用の変数
        self.selected_atoms_3d = set()
        self.atom_selection_mode = False
        self.selected_atom_actors = []
        
        # 3D編集用の原子選択状態
        self.selected_atoms_3d = set()  # 3Dビューで選択された原子のインデックス
        
        # 3D編集ダイアログの参照を保持
        self.active_3d_dialogs = []
        
        self.init_ui()
        self.init_worker_thread()
        self._setup_3d_picker() 

        # --- RDKit初回実行コストの事前読み込み（ウォームアップ）---
        try:
            # Create a molecule with a variety of common atoms to ensure
            # the valence/H-count machinery is fully initialized.
            warmup_smiles = "OC(N)C(S)P"
            warmup_mol = Chem.MolFromSmiles(warmup_smiles)
            if warmup_mol:
                for atom in warmup_mol.GetAtoms():
                    atom.GetNumImplicitHs()
        except Exception as e:
            print(f"RDKit warm-up failed: {e}")

        self.reset_undo_stack()
        self.scene.selectionChanged.connect(self.update_edit_menu_actions)
        QApplication.clipboard().dataChanged.connect(self.update_edit_menu_actions)

        self.update_edit_menu_actions()

        if initial_file:
            self.load_command_line_file(initial_file)
        
        QTimer.singleShot(0, self.apply_initial_settings)
        # カメラ初期化フラグ（初回描画時のみリセットを許可する）
        self._camera_initialized = False
        
        # 初期メニューテキストと状態を設定
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()
        
        # 初期化完了を設定
        self.initialization_complete = True
        self.update_window_title()  # 初期化完了後にタイトルを更新
        # Ensure initial keyboard/mouse focus is placed on the 2D view
        # when opening a file or starting the application. This avoids
        # accidental focus landing on toolbar/buttons (e.g. Optimize 2D).
        try:
            QTimer.singleShot(0, self.view_2d.setFocus)
        except Exception:
            pass

    def init_ui(self):
        # 1. 現在のスクリプトがあるディレクトリのパスを取得
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 2. 'assets'フォルダ内のアイコンファイルへのフルパスを構築
        icon_path = os.path.join(script_dir, 'assets', 'icon.png')
        
        # 3. ファイルパスから直接QIconオブジェクトを作成
        if os.path.exists(icon_path): # ファイルが存在するか確認
            app_icon = QIcon(icon_path)
            
            # 4. ウィンドウにアイコンを設定
            self.setWindowIcon(app_icon)
        else:
            print(f"警告: アイコンファイルが見つかりません: {icon_path}")

        self.init_menu_bar()

        self.splitter = QSplitter(Qt.Orientation.Horizontal)
        # スプリッターハンドルを太くして視認性を向上
        self.splitter.setHandleWidth(8)
        # スプリッターハンドルのスタイルを改善
        self.splitter.setStyleSheet("""
            QSplitter::handle {
                background-color: #ccc;
                border: 1px solid #999;
                border-radius: 4px;
                margin: 2px;
            }
            QSplitter::handle:hover {
                background-color: #aaa;
            }
            QSplitter::handle:pressed {
                background-color: #888;
            }
        """)
        self.setCentralWidget(self.splitter)

        left_pane=QWidget()
        left_pane.setAcceptDrops(True)
        left_layout=QVBoxLayout(left_pane)

        self.scene=MoleculeScene(self.data,self)
        self.scene.setSceneRect(-4000,-4000,4000,4000)
        self.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.view_2d=ZoomableView(self.scene, self)
        self.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view_2d.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        left_layout.addWidget(self.view_2d, 1)

        self.view_2d.scale(0.75, 0.75)

        # --- 左パネルのボタンレイアウト ---
        left_buttons_layout = QHBoxLayout()
        self.cleanup_button = QPushButton("Optimize 2D")
        self.cleanup_button.clicked.connect(self.clean_up_2d_structure)
        left_buttons_layout.addWidget(self.cleanup_button)

        self.convert_button = QPushButton("Convert 2D to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        # Allow right-click to open a temporary conversion-mode menu
        try:
            self.convert_button.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            self.convert_button.customContextMenuRequested.connect(self.show_convert_menu)
        except Exception:
            pass
        left_buttons_layout.addWidget(self.convert_button)
        
        left_layout.addLayout(left_buttons_layout)
        self.splitter.addWidget(left_pane)

        # --- 右パネルとボタンレイアウト ---
        right_pane = QWidget()
        # 1. 右パネル全体は「垂直」レイアウトにする
        right_layout = QVBoxLayout(right_pane)
        self.plotter = CustomQtInteractor(right_pane, main_window=self, lighting='none')
        self.plotter.setAcceptDrops(False)
        self.plotter.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        # 2. 垂直レイアウトに3Dビューを追加
        right_layout.addWidget(self.plotter, 1)
        #self.plotter.installEventFilter(self)

        # 3. ボタンをまとめるための「水平」レイアウトを作成
        right_buttons_layout = QHBoxLayout()

        # 3D最適化ボタン
        self.optimize_3d_button = QPushButton("Optimize 3D")
        self.optimize_3d_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.optimize_3d_button.clicked.connect(self.optimize_3d_structure)
        self.optimize_3d_button.setEnabled(False)
        # 初期状態は_enable_3d_features(False)で統一的に設定
        # Allow right-click to open a temporary optimization-method menu
        try:
            self.optimize_3d_button.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            self.optimize_3d_button.customContextMenuRequested.connect(self.show_optimize_menu)
        except Exception:
            pass
        right_buttons_layout.addWidget(self.optimize_3d_button)

        # エクスポートボタン (メニュー付き)
        self.export_button = QToolButton()
        self.export_button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.export_button.setText("Export 3D")
        self.export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.export_button.setEnabled(False) # 初期状態は無効

        export_menu = QMenu(self)
        export_mol_action = QAction("Export as MOL...", self)
        export_mol_action.triggered.connect(self.save_3d_as_mol)
        export_menu.addAction(export_mol_action)

        export_xyz_action = QAction("Export as XYZ...", self)
        export_xyz_action.triggered.connect(self.save_as_xyz)
        export_menu.addAction(export_xyz_action)

        export_png_action = QAction("Export as PNG...", self)
        export_png_action.triggered.connect(self.export_3d_png)
        export_menu.addAction(export_png_action)

        self.export_button.setMenu(export_menu)
        right_buttons_layout.addWidget(self.export_button)

        # 4. 水平のボタンレイアウトを、全体の垂直レイアウトに追加
        right_layout.addLayout(right_buttons_layout)
        self.splitter.addWidget(right_pane)
        
        # スプリッターのサイズ変更をモニターして、フィードバックを提供
        self.splitter.splitterMoved.connect(self.on_splitter_moved)
        
        self.splitter.setSizes([600, 600])
        
        # スプリッターハンドルにツールチップを設定
        QTimer.singleShot(100, self.setup_splitter_tooltip)

        # ステータスバーを左右に分離するための設定
        self.status_bar = self.statusBar()
        self.formula_label = QLabel("")  # 右側に表示するラベルを作成
        # 右端に余白を追加して見栄えを調整
        self.formula_label.setStyleSheet("padding-right: 8px;")
        # ラベルを右側に常時表示ウィジェットとして追加
        self.status_bar.addPermanentWidget(self.formula_label)

        #self.view_2d.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

        # Top/main toolbar (keep 3D Edit controls on the right end of this toolbar)
        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        # Keep a reference to the main toolbar for later updates
        self.toolbar = toolbar

        # Templates toolbar: place it directly below the main toolbar (second row at the top)
        # Use addToolBarBreak to ensure this toolbar appears on the next row under the main toolbar.
        # Some older PyQt/PySide versions may not have addToolBarBreak; fall back silently in that case.
        try:
            # Insert a toolbar break in the Top toolbar area to force the next toolbar onto a new row
            self.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        except Exception:
            # If addToolBarBreak isn't available, continue without raising; placement may still work depending on the platform.
            pass

        toolbar_bottom = QToolBar("Templates Toolbar")
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, toolbar_bottom)
        self.toolbar_bottom = toolbar_bottom
        self.tool_group = QActionGroup(self)
        self.tool_group.setExclusive(True)

        actions_data = [
            ("Select", 'select', 'Space'), ("C", 'atom_C', 'c'), ("H", 'atom_H', 'h'), ("B", 'atom_B', 'b'),
            ("N", 'atom_N', 'n'), ("O", 'atom_O', 'o'), ("S", 'atom_S', 's'), ("Si", 'atom_Si', 'Shift+S'), ("P", 'atom_P', 'p'), 
            ("F", 'atom_F', 'f'), ("Cl", 'atom_Cl', 'Shift+C'), ("Br", 'atom_Br', 'Shift+B'), ("I", 'atom_I', 'i'), 
            ("Other...", 'atom_other', '')
        ]

        for text, mode, shortcut_text in actions_data:
            if text == "C": toolbar.addSeparator()
            
            action = QAction(text, self, checkable=(mode != 'atom_other'))
            if shortcut_text: action.setToolTip(f"{text} ({shortcut_text})")

            if mode == 'atom_other':
                action.triggered.connect(self.open_periodic_table_dialog)
                self.other_atom_action = action
            else:
                action.triggered.connect(lambda c, m=mode: self.set_mode(m))
                self.mode_actions[mode] = action

            toolbar.addAction(action)
            if mode != 'atom_other': self.tool_group.addAction(action)
            
            if text == "Select":
                select_action = action
        
        toolbar.addSeparator()

        # --- アイコン前景色を決めるヘルパー（ダーク/ライトモード対応） ---
        def _icon_foreground_color():
            """Return a QColor for icon foreground (black on light backgrounds, white on dark backgrounds).

            Priority: explicit setting 'icon_foreground' in settings -> infer from configured background color -> infer from application palette.
            """
            try:
                fg_hex = self.settings.get('icon_foreground')
                if fg_hex:
                    c = QColor(fg_hex)
                    if c.isValid():
                        return c
            except Exception:
                pass

            try:
                bg_hex = self.settings.get('background_color')
                if bg_hex:
                    bg = QColor(bg_hex)
                    if bg.isValid():
                        lum = 0.2126 * bg.redF() + 0.7152 * bg.greenF() + 0.0722 * bg.blueF()
                        return QColor('#FFFFFF') if lum < 0.5 else QColor('#000000')
            except Exception:
                pass

            try:
                pal = QApplication.palette()
                # palette.window() returns a QBrush; call color()
                window_bg = pal.window().color()
                lum = 0.2126 * window_bg.redF() + 0.7152 * window_bg.greenF() + 0.0722 * window_bg.blueF()
                return QColor('#FFFFFF') if lum < 0.5 else QColor('#000000')
            except Exception:
                return QColor('#000000')

        # --- 結合ボタンのアイコンを生成するヘルパー関数 ---
        def create_bond_icon(bond_type, size=32):
            fg = _icon_foreground_color()
            pixmap = QPixmap(size, size)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)

            p1 = QPointF(6, size / 2)
            p2 = QPointF(size - 6, size / 2)
            line = QLineF(p1, p2)

            pen = QPen(fg, 2)
            painter.setPen(pen)
            painter.setBrush(QBrush(fg))

            if bond_type == 'single':
                painter.drawLine(line)
            elif bond_type == 'double':
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * 2.5
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))
            elif bond_type == 'triple':
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * 3.0
                painter.drawLine(line)
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))
            elif bond_type == 'wedge':
                vec = line.unitVector()
                normal = vec.normalVector()
                offset = QPointF(normal.dx(), normal.dy()) * 5.0
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)
            elif bond_type == 'dash':
                vec = line.unitVector()
                normal = vec.normalVector()

                num_dashes = NUM_DASHES
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = 10 * t
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.setPen(QPen(fg, 1.5))
                    painter.drawLine(start_pt - offset, start_pt + offset)

            elif bond_type == 'ez_toggle':
                # アイコン下部に二重結合を描画
                p1 = QPointF(6, size * 0.75)
                p2 = QPointF(size - 6, size * 0.75)
                line = QLineF(p1, p2)
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * 2.0
                painter.setPen(QPen(fg, 2))
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))
                # 上部に "Z⇌E" のテキストを描画
                painter.setPen(QPen(fg, 1))
                font = painter.font()
                font.setPointSize(10)
                font.setBold(True)
                painter.setFont(font)
                text_rect = QRectF(0, 0, size, size * 0.6)
                # U+21CC は右向きと左向きのハープーンが重なった記号 (⇌)
                painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, "Z⇌E")

            painter.end()
            return QIcon(pixmap)

        # --- 結合ボタンをツールバーに追加 ---
        bond_actions_data = [
            ("Single Bond", 'bond_1_0', '1', 'single'),
            ("Double Bond", 'bond_2_0', '2', 'double'),
            ("Triple Bond", 'bond_3_0', '3', 'triple'),
            ("Wedge Bond", 'bond_1_1', 'W', 'wedge'),
            ("Dash Bond", 'bond_1_2', 'D', 'dash'),
            ("Toggle E/Z", 'bond_2_5', 'E/Z', 'ez_toggle'),
        ]

        for text, mode, shortcut_text, icon_type in bond_actions_data:
            action = QAction(self)
            action.setIcon(create_bond_icon(icon_type))
            action.setToolTip(f"{text} ({shortcut_text})")
            action.setCheckable(True)
            action.triggered.connect(lambda checked, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            toolbar.addAction(action)
            self.tool_group.addAction(action)
        
        toolbar.addSeparator()

        charge_plus_action = QAction("+ Charge", self, checkable=True)
        charge_plus_action.setToolTip("Increase Atom Charge (+)")
        charge_plus_action.triggered.connect(lambda c, m='charge_plus': self.set_mode(m))
        self.mode_actions['charge_plus'] = charge_plus_action
        toolbar.addAction(charge_plus_action)
        self.tool_group.addAction(charge_plus_action)

        charge_minus_action = QAction("- Charge", self, checkable=True)
        charge_minus_action.setToolTip("Decrease Atom Charge (-)")
        charge_minus_action.triggered.connect(lambda c, m='charge_minus': self.set_mode(m))
        self.mode_actions['charge_minus'] = charge_minus_action
        toolbar.addAction(charge_minus_action)
        self.tool_group.addAction(charge_minus_action)

        radical_action = QAction("Radical", self, checkable=True)
        radical_action.setToolTip("Toggle Radical (0/1/2) (.)")
        radical_action.triggered.connect(lambda c, m='radical': self.set_mode(m))
        self.mode_actions['radical'] = radical_action
        toolbar.addAction(radical_action)
        self.tool_group.addAction(radical_action)

        # We will show template controls in the bottom toolbar to improve layout.
        # Add a small label to the bottom toolbar instead of the main toolbar.
        toolbar_bottom.addWidget(QLabel(" Templates:"))
        
        # --- アイコンを生成するヘルパー関数 ---
        def create_template_icon(n, is_benzene=False):
            size = 32
            fg = _icon_foreground_color()
            pixmap = QPixmap(size, size)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setPen(QPen(fg, 2))

            center = QPointF(size / 2, size / 2)
            radius = size / 2 - 4 # アイコンの余白

            points = []
            angle_step = 2 * math.pi / n
            # ポリゴンが直立するように開始角度を調整
            start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2

            for i in range(n):
                angle = start_angle + i * angle_step
                x = center.x() + radius * math.cos(angle)
                y = center.y() + radius * math.sin(angle)
                points.append(QPointF(x, y))

            painter.drawPolygon(QPolygonF(points))

            if is_benzene:
                painter.drawEllipse(center, radius * 0.6, radius * 0.6)

            if n in [7, 8, 9]:
                font = QFont("Arial", 10, QFont.Weight.Bold)
                painter.setFont(font)
                painter.setPen(QPen(fg, 1))
                painter.drawText(QRectF(0, 0, size, size), Qt.AlignmentFlag.AlignCenter, str(n))

            painter.end()
            return QIcon(pixmap)

        # --- ヘルパー関数を使ってアイコン付きボタンを作成 ---
        templates = [("Benzene", "template_benzene", 6)] + [(f"{i}-Ring", f"template_{i}", i) for i in range(3, 10)]
        for text, mode, n in templates:
            action = QAction(self) # テキストなしでアクションを作成
            action.setCheckable(True)

            is_benzene = (text == "Benzene")
            icon = create_template_icon(n, is_benzene=is_benzene)
            action.setIcon(icon) # アイコンを設定

            if text == "Benzene":
                action.setToolTip(f"{text} Template (4)")
            else:
                action.setToolTip(f"{text} Template")

            action.triggered.connect(lambda c, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            # Add template actions to the bottom toolbar so templates are on the second line
            toolbar_bottom.addAction(action)
            self.tool_group.addAction(action)

        # Add USER button for user templates (placed in bottom toolbar)
        user_template_action = QAction("USER", self)
        user_template_action.setCheckable(True)
        user_template_action.setToolTip("Open User Templates Dialog")
        user_template_action.triggered.connect(self.open_template_dialog_and_activate)
        self.mode_actions['template_user'] = user_template_action
        toolbar_bottom.addAction(user_template_action)
        self.tool_group.addAction(user_template_action)

        # 初期モードを'select'から'atom_C'（炭素原子描画モード）に変更
        self.set_mode('atom_C')
        # 対応するツールバーの'C'ボタンを選択状態にする
        if 'atom_C' in self.mode_actions:
            self.mode_actions['atom_C'].setChecked(True)

        # スペーサーを追加して、次のウィジェットを右端に配置する (keep on top toolbar)
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        toolbar.addWidget(spacer)

        # 測定機能ボタンを追加（"3D Select"に変更）
        self.measurement_action = QAction("3D Select", self, checkable=True)
        self.measurement_action.setToolTip("Enable distance, angle, and dihedral measurement in 3D view")
        # 初期状態でも有効にする
        self.measurement_action.triggered.connect(self.toggle_measurement_mode)
        toolbar.addAction(self.measurement_action)

        self.edit_3d_action = QAction("3D Drag", self, checkable=True)
        self.edit_3d_action.setToolTip("Toggle 3D atom dragging mode (Hold Alt for temporary mode)")
        # 初期状態でも有効にする
        self.edit_3d_action.toggled.connect(self.toggle_3d_edit_mode)
        toolbar.addAction(self.edit_3d_action)

        # 3Dスタイル変更ボタンとメニューを作成

        self.style_button = QToolButton()
        self.style_button.setText("3D Style")
        self.style_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(self.style_button)

        style_menu = QMenu(self)
        self.style_button.setMenu(style_menu)

        style_group = QActionGroup(self)
        style_group.setExclusive(True)

        # Ball & Stick アクション
        bs_action = QAction("Ball & Stick", self, checkable=True)
        bs_action.setChecked(True)
        bs_action.triggered.connect(lambda: self.set_3d_style('ball_and_stick'))
        style_menu.addAction(bs_action)
        style_group.addAction(bs_action)

        # CPK アクション
        cpk_action = QAction("CPK (Space-filling)", self, checkable=True)
        cpk_action.triggered.connect(lambda: self.set_3d_style('cpk'))
        style_menu.addAction(cpk_action)
        style_group.addAction(cpk_action)

        # Wireframe アクション
        wireframe_action = QAction("Wireframe", self, checkable=True)
        wireframe_action.triggered.connect(lambda: self.set_3d_style('wireframe'))
        style_menu.addAction(wireframe_action)
        style_group.addAction(wireframe_action)

        # Stick アクション
        stick_action = QAction("Stick", self, checkable=True)
        stick_action.triggered.connect(lambda: self.set_3d_style('stick'))
        style_menu.addAction(stick_action)
        style_group.addAction(stick_action)

        quit_shortcut = QShortcut(QKeySequence("Ctrl+Q"), self)
        quit_shortcut.activated.connect(self.close)

        self.view_2d.setFocus()

    def init_menu_bar(self):
        menu_bar = self.menuBar()
        
        file_menu = menu_bar.addMenu("&File")
        
        # === プロジェクト操作 ===
        new_action = QAction("&New", self)
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.clear_all)
        file_menu.addAction(new_action)
        
        load_project_action = QAction("&Open Project...", self)
        load_project_action.setShortcut("Ctrl+O")
        load_project_action.triggered.connect(self.open_project_file)
        file_menu.addAction(load_project_action)
        
        save_action = QAction("&Save Project", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_project)
        file_menu.addAction(save_action)
        
        save_as_action = QAction("Save Project &As...", self)
        save_as_action.setShortcut("Ctrl+Shift+S")
        save_as_action.triggered.connect(self.save_project_as)
        file_menu.addAction(save_as_action)
        
        save_template_action = QAction("Save 2D as Template...", self)
        save_template_action.triggered.connect(self.save_2d_as_template)
        file_menu.addAction(save_template_action)
        
        file_menu.addSeparator()
        
        # === インポート ===
        import_menu = file_menu.addMenu("Import")
        
        load_mol_action = QAction("MOL/SDF File...", self)
        load_mol_action.triggered.connect(self.load_mol_file)
        import_menu.addAction(load_mol_action)
        
        import_smiles_action = QAction("SMILES...", self)
        import_smiles_action.triggered.connect(self.import_smiles_dialog)
        import_menu.addAction(import_smiles_action)
        
        import_inchi_action = QAction("InChI...", self)
        import_inchi_action.triggered.connect(self.import_inchi_dialog)
        import_menu.addAction(import_inchi_action)
        
        import_menu.addSeparator()
        
        load_3d_mol_action = QAction("3D MOL/SDF (3D View Only)...", self)
        load_3d_mol_action.triggered.connect(self.load_mol_file_for_3d_viewing)
        import_menu.addAction(load_3d_mol_action)
        
        load_3d_xyz_action = QAction("3D XYZ (3D View Only)...", self)
        load_3d_xyz_action.triggered.connect(self.load_xyz_for_3d_viewing)
        import_menu.addAction(load_3d_xyz_action)
        
        # === エクスポート ===
        export_menu = file_menu.addMenu("Export")
        
        # プロジェクト形式エクスポート
        export_pmeraw_action = QAction("PME Raw Format...", self)
        export_pmeraw_action.triggered.connect(self.save_raw_data)
        export_menu.addAction(export_pmeraw_action)
        
        export_menu.addSeparator()
        
        # 2D エクスポート
        export_2d_menu = export_menu.addMenu("2D Formats")
        save_mol_action = QAction("MOL File...", self)
        save_mol_action.triggered.connect(self.save_as_mol)
        export_2d_menu.addAction(save_mol_action)
        
        export_2d_png_action = QAction("PNG Image...", self)
        export_2d_png_action.triggered.connect(self.export_2d_png)
        export_2d_menu.addAction(export_2d_png_action)
        
        # 3D エクスポート
        export_3d_menu = export_menu.addMenu("3D Formats")
        save_3d_mol_action = QAction("MOL File...", self)
        save_3d_mol_action.triggered.connect(self.save_3d_as_mol)
        export_3d_menu.addAction(save_3d_mol_action)
        
        save_xyz_action = QAction("XYZ File...", self)
        save_xyz_action.triggered.connect(self.save_as_xyz)
        export_3d_menu.addAction(save_xyz_action)
        
        export_3d_png_action = QAction("PNG Image...", self)
        export_3d_png_action.triggered.connect(self.export_3d_png)
        export_3d_menu.addAction(export_3d_png_action)
        
        export_3d_menu.addSeparator()
        
        export_stl_action = QAction("STL File...", self)
        export_stl_action.triggered.connect(self.export_stl)
        export_3d_menu.addAction(export_stl_action)
        
        export_obj_action = QAction("OBJ/MTL (with colors)...", self)
        export_obj_action.triggered.connect(self.export_obj_mtl)
        export_3d_menu.addAction(export_obj_action)
        
        file_menu.addSeparator()
        quit_action = QAction("Quit", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)
        
        edit_menu = menu_bar.addMenu("&Edit")
        self.undo_action = QAction("Undo", self); self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.triggered.connect(self.undo); edit_menu.addAction(self.undo_action)
        
        self.redo_action = QAction("Redo", self); self.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.redo_action.triggered.connect(self.redo); edit_menu.addAction(self.redo_action)
        
        edit_menu.addSeparator()

        self.cut_action = QAction("Cut", self)
        self.cut_action.setShortcut(QKeySequence.StandardKey.Cut)
        self.cut_action.triggered.connect(self.cut_selection)
        edit_menu.addAction(self.cut_action)

        self.copy_action = QAction("Copy", self)
        self.copy_action.setShortcut(QKeySequence.StandardKey.Copy)
        self.copy_action.triggered.connect(self.copy_selection)
        edit_menu.addAction(self.copy_action)
        
        self.paste_action = QAction("Paste", self)
        self.paste_action.setShortcut(QKeySequence.StandardKey.Paste)
        self.paste_action.triggered.connect(self.paste_from_clipboard)
        edit_menu.addAction(self.paste_action)

        edit_menu.addSeparator()

        add_hydrogen_action = QAction("Add Hydrogens", self)
        add_hydrogen_action.setToolTip("Add explicit hydrogens based on RDKit implicit counts")
        add_hydrogen_action.triggered.connect(self.add_hydrogen_atoms)
        edit_menu.addAction(add_hydrogen_action)
    
        remove_hydrogen_action = QAction("Remove Hydrogens", self)
        remove_hydrogen_action.triggered.connect(self.remove_hydrogen_atoms)
        edit_menu.addAction(remove_hydrogen_action)

        edit_menu.addSeparator()

        optimize_2d_action = QAction("Optimize 2D", self)
        optimize_2d_action.setShortcut(QKeySequence("Ctrl+J"))
        optimize_2d_action.triggered.connect(self.clean_up_2d_structure)
        edit_menu.addAction(optimize_2d_action)
        
        convert_3d_action = QAction("Convert 2D to 3D", self)
        convert_3d_action.setShortcut(QKeySequence("Ctrl+K"))
        convert_3d_action.triggered.connect(self.trigger_conversion)
        edit_menu.addAction(convert_3d_action)

        optimize_3d_action = QAction("Optimize 3D", self)
        optimize_3d_action.setShortcut(QKeySequence("Ctrl+L")) 
        optimize_3d_action.triggered.connect(self.optimize_3d_structure)
        edit_menu.addAction(optimize_3d_action)

        # Note: 3D Optimization Settings moved to Settings -> "3D Optimization Settings"
        # to avoid duplicating the same submenu in both Edit and Settings.

        # Note: Open Babel-based optimization menu entries were intentionally
        # removed above. Open Babel (pybel) is still available for conversion
        # fallback elsewhere in the code, so we don't disable menu items here.

        edit_menu.addSeparator()
        
        select_all_action = QAction("Select All", self); select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.select_all); edit_menu.addAction(select_all_action)
        
        clear_all_action = QAction("Clear All", self)
        clear_all_action.setShortcut(QKeySequence("Ctrl+Shift+C"))
        clear_all_action.triggered.connect(self.clear_all); edit_menu.addAction(clear_all_action)

        view_menu = menu_bar.addMenu("&View")

        zoom_in_action = QAction("Zoom In", self)
        zoom_in_action.setShortcut(QKeySequence.StandardKey.ZoomIn) # Ctrl +
        zoom_in_action.triggered.connect(self.zoom_in)
        view_menu.addAction(zoom_in_action)

        zoom_out_action = QAction("Zoom Out", self)
        zoom_out_action.setShortcut(QKeySequence.StandardKey.ZoomOut) # Ctrl -
        zoom_out_action.triggered.connect(self.zoom_out)
        view_menu.addAction(zoom_out_action)

        reset_zoom_action = QAction("Reset Zoom", self)
        reset_zoom_action.setShortcut(QKeySequence("Ctrl+0"))
        reset_zoom_action.triggered.connect(self.reset_zoom)
        view_menu.addAction(reset_zoom_action)
        
        fit_action = QAction("Fit to View", self)
        fit_action.setShortcut(QKeySequence("Ctrl+9"))
        fit_action.triggered.connect(self.fit_to_view)
        view_menu.addAction(fit_action)

        view_menu.addSeparator()

        reset_3d_view_action = QAction("Reset 3D View", self)
        reset_3d_view_action.triggered.connect(lambda: self.plotter.reset_camera() if hasattr(self, 'plotter') else None)
        reset_3d_view_action.setShortcut(QKeySequence("Ctrl+R"))
        view_menu.addAction(reset_3d_view_action)
        
        view_menu.addSeparator()

        # Panel Layout submenu
        layout_menu = view_menu.addMenu("Panel Layout")
        
        equal_panels_action = QAction("Equal Panels (50:50)", self)
        equal_panels_action.setShortcut(QKeySequence("Ctrl+1"))
        equal_panels_action.triggered.connect(lambda: self.set_panel_layout(50, 50))
        layout_menu.addAction(equal_panels_action)
        
        layout_2d_focus_action = QAction("2D Focus (70:30)", self)
        layout_2d_focus_action.setShortcut(QKeySequence("Ctrl+2"))
        layout_2d_focus_action.triggered.connect(lambda: self.set_panel_layout(70, 30))
        layout_menu.addAction(layout_2d_focus_action)
        
        layout_3d_focus_action = QAction("3D Focus (30:70)", self)
        layout_3d_focus_action.setShortcut(QKeySequence("Ctrl+3"))
        layout_3d_focus_action.triggered.connect(lambda: self.set_panel_layout(30, 70))
        layout_menu.addAction(layout_3d_focus_action)
        
        layout_menu.addSeparator()
        
        toggle_2d_panel_action = QAction("Toggle 2D Panel", self)
        toggle_2d_panel_action.setShortcut(QKeySequence("Ctrl+H"))
        toggle_2d_panel_action.triggered.connect(self.toggle_2d_panel)
        layout_menu.addAction(toggle_2d_panel_action)

        view_menu.addSeparator()

        self.toggle_chiral_action = QAction("Show Chiral Labels", self, checkable=True)
        self.toggle_chiral_action.setChecked(self.show_chiral_labels)
        self.toggle_chiral_action.triggered.connect(self.toggle_chiral_labels_display)
        view_menu.addAction(self.toggle_chiral_action)

        view_menu.addSeparator()

        # 3D Atom Info submenu
        atom_info_menu = view_menu.addMenu("3D Atom Info Display")
        
        self.show_atom_id_action = QAction("Show Original ID / Index", self, checkable=True)
        self.show_atom_id_action.triggered.connect(lambda: self.toggle_atom_info_display('id'))
        atom_info_menu.addAction(self.show_atom_id_action)
        
        self.show_rdkit_id_action = QAction("Show RDKit Index", self, checkable=True)
        self.show_rdkit_id_action.triggered.connect(lambda: self.toggle_atom_info_display('rdkit_id'))
        atom_info_menu.addAction(self.show_rdkit_id_action)
        
        self.show_atom_coords_action = QAction("Show Coordinates (X,Y,Z)", self, checkable=True)
        self.show_atom_coords_action.triggered.connect(lambda: self.toggle_atom_info_display('coords'))
        atom_info_menu.addAction(self.show_atom_coords_action)
        
        self.show_atom_symbol_action = QAction("Show Element Symbol", self, checkable=True)
        self.show_atom_symbol_action.triggered.connect(lambda: self.toggle_atom_info_display('symbol'))
        atom_info_menu.addAction(self.show_atom_symbol_action)

        analysis_menu = menu_bar.addMenu("&Analysis")
        self.analysis_action = QAction("Show Analysis...", self)
        self.analysis_action.triggered.connect(self.open_analysis_window)
        self.analysis_action.setEnabled(False)
        analysis_menu.addAction(self.analysis_action)

        # 3D Edit menu
        edit_3d_menu = menu_bar.addMenu("3D &Edit")
        
        # Translation action
        translation_action = QAction("Translation...", self)
        translation_action.triggered.connect(self.open_translation_dialog)
        translation_action.setEnabled(False)
        edit_3d_menu.addAction(translation_action)
        self.translation_action = translation_action
        
        # Move Group action
        move_group_action = QAction("Move Group...", self)
        move_group_action.triggered.connect(self.open_move_group_dialog)
        move_group_action.setEnabled(False)
        edit_3d_menu.addAction(move_group_action)
        self.move_group_action = move_group_action
        
        edit_3d_menu.addSeparator()
        
        # Alignment submenu (統合)
        align_menu = edit_3d_menu.addMenu("Align to")
        align_menu.setEnabled(False)
        self.align_menu = align_menu
        
        # Axis alignment submenu
        axis_align_menu = align_menu.addMenu("Axis")
        
        align_x_action = QAction("X-axis", self)
        align_x_action.triggered.connect(lambda: self.open_alignment_dialog('x'))
        align_x_action.setEnabled(False)
        axis_align_menu.addAction(align_x_action)
        self.align_x_action = align_x_action
        
        align_y_action = QAction("Y-axis", self)
        align_y_action.triggered.connect(lambda: self.open_alignment_dialog('y'))
        align_y_action.setEnabled(False)
        axis_align_menu.addAction(align_y_action)
        self.align_y_action = align_y_action
        
        align_z_action = QAction("Z-axis", self)
        align_z_action.triggered.connect(lambda: self.open_alignment_dialog('z'))
        align_z_action.setEnabled(False)
        axis_align_menu.addAction(align_z_action)
        self.align_z_action = align_z_action
        
        # Plane alignment submenu (旧align)
        plane_align_menu = align_menu.addMenu("Plane")
        
        alignplane_xy_action = QAction("XY-plane", self)
        alignplane_xy_action.triggered.connect(lambda: self.open_align_plane_dialog('xy'))
        alignplane_xy_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xy_action)
        self.alignplane_xy_action = alignplane_xy_action

        alignplane_xz_action = QAction("XZ-plane", self)
        alignplane_xz_action.triggered.connect(lambda: self.open_align_plane_dialog('xz'))
        alignplane_xz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_xz_action)
        self.alignplane_xz_action = alignplane_xz_action

        alignplane_yz_action = QAction("YZ-plane", self)
        alignplane_yz_action.triggered.connect(lambda: self.open_align_plane_dialog('yz'))
        alignplane_yz_action.setEnabled(False)
        plane_align_menu.addAction(alignplane_yz_action)
        self.alignplane_yz_action = alignplane_yz_action

        edit_3d_menu.addSeparator()

        # Mirror action
        mirror_action = QAction("Mirror...", self)
        mirror_action.triggered.connect(self.open_mirror_dialog)
        mirror_action.setEnabled(False)
        edit_3d_menu.addAction(mirror_action)
        self.mirror_action = mirror_action

        edit_3d_menu.addSeparator()
        
        # Planarize selection (best-fit plane)
        planarize_action = QAction("Planarize...", self)
        planarize_action.triggered.connect(lambda: self.open_planarize_dialog(None))
        planarize_action.setEnabled(False)
        edit_3d_menu.addAction(planarize_action)
        self.planarize_action = planarize_action
        
        edit_3d_menu.addSeparator()
        
        # Bond length conversion
        bond_length_action = QAction("Adjust Bond Length...", self)
        bond_length_action.triggered.connect(self.open_bond_length_dialog)
        bond_length_action.setEnabled(False)
        edit_3d_menu.addAction(bond_length_action)
        self.bond_length_action = bond_length_action
        
        # Angle conversion
        angle_action = QAction("Adjust Angle...", self)
        angle_action.triggered.connect(self.open_angle_dialog)
        angle_action.setEnabled(False)
        edit_3d_menu.addAction(angle_action)
        self.angle_action = angle_action
        
        # Dihedral angle conversion
        dihedral_action = QAction("Adjust Dihedral Angle...", self)
        dihedral_action.triggered.connect(self.open_dihedral_dialog)
        dihedral_action.setEnabled(False)
        edit_3d_menu.addAction(dihedral_action)
        self.dihedral_action = dihedral_action

        edit_3d_menu.addSeparator()
        
        # Constrained Optimization action
        constrained_opt_action = QAction("Constrained Optimization...", self)
        constrained_opt_action.triggered.connect(self.open_constrained_optimization_dialog)
        constrained_opt_action.setEnabled(False)  # 3Dモデルロード時に有効化
        edit_3d_menu.addAction(constrained_opt_action)
        self.constrained_opt_action = constrained_opt_action

        settings_menu = menu_bar.addMenu("&Settings")
        # 1) 3D View settings (existing)
        view_settings_action = QAction("3D View Settings...", self)
        view_settings_action.triggered.connect(self.open_settings_dialog)
        settings_menu.addAction(view_settings_action)
        
        # Color settings (CPK/Bond) — keep with other settings
        color_action = QAction("CPK Colors...", self)
        color_action.triggered.connect(lambda: ColorSettingsDialog(self.settings, parent=self).exec_())
        settings_menu.addAction(color_action)
    
        # 2) 3D Conversion settings — submenu with radio/check actions
        conversion_menu = settings_menu.addMenu("3D Conversion")
        conv_group = QActionGroup(self)
        conv_group.setExclusive(True)
        # helper to set conversion mode and persist
        def _set_conv_mode(mode):
            try:
                self.settings['3d_conversion_mode'] = mode
                # defer disk write
                try:
                    self.settings_dirty = True
                except Exception:
                    pass
                self.statusBar().showMessage(f"3D conversion mode set to: {mode}")
            except Exception:
                pass

        conv_options = [
            ("RDKit -> Open Babel (fallback)", 'fallback'),
            ("RDKit only", 'rdkit'),
            ("Open Babel only", 'obabel'),
            ("Direct (use 2D coords + add H)", 'direct')
        ]
        self.conv_actions = {}
        for label, key in conv_options:
            a = QAction(label, self)
            a.setCheckable(True)
            # If Open Babel isn't available, disable the Open Babel-only option
            # and also disable the fallback option since it depends on Open Babel.
            if not OBABEL_AVAILABLE:
                if key == 'obabel' or key == 'fallback':
                    a.setEnabled(False)
            a.triggered.connect(lambda checked, m=key: _set_conv_mode(m))
            conversion_menu.addAction(a)
            conv_group.addAction(a)
            self.conv_actions[key] = a

        # Initialize checked state from settings (fallback default)
        # Determine saved conversion mode. If Open Babel is not available,
        # prefer 'rdkit' as the default rather than 'fallback'. Also ensure
        # the settings reflect the actual enabled choice.
        try:
            default_mode = 'rdkit' if not OBABEL_AVAILABLE else 'fallback'
            saved_conv = self.settings.get('3d_conversion_mode', default_mode)
        except Exception:
            saved_conv = 'rdkit' if not OBABEL_AVAILABLE else 'fallback'

        # If the saved mode is disabled/unavailable, fall back to an enabled option.
        if saved_conv not in self.conv_actions or not self.conv_actions[saved_conv].isEnabled():
            # Prefer 'rdkit' if available, else pick whichever action is enabled
            preferred = 'rdkit' if 'rdkit' in self.conv_actions and self.conv_actions['rdkit'].isEnabled() else None
            if not preferred:
                for k, act in self.conv_actions.items():
                    if act.isEnabled():
                        preferred = k
                        break
            saved_conv = preferred or 'rdkit'

        # Set the checked state and persist the chosen conversion mode
        try:
            if saved_conv in self.conv_actions:
                try:
                    self.conv_actions[saved_conv].setChecked(True)
                except Exception:
                    pass
            self.settings['3d_conversion_mode'] = saved_conv
            try:
                self.settings_dirty = True
            except Exception:
                pass
        except Exception:
            pass

        # 3) 3D Optimization Settings (single location under Settings menu)
        optimization_menu = settings_menu.addMenu("3D Optimization Settings")

        # Only RDKit-backed optimization methods are offered here.
        opt_methods = [
            ("MMFF94s", "MMFF_RDKIT"),
            ("MMFF94", "MMFF94_RDKIT"),
            ("UFF", "UFF_RDKIT"),
        ]

        # Map key -> human-readable label for status messages and later lookups
        try:
            self.opt3d_method_labels = {key.upper(): label for (label, key) in opt_methods}
        except Exception:
            self.opt3d_method_labels = {}

        opt_group = QActionGroup(self)
        opt_group.setExclusive(True)
        opt_actions = {}
        for label, key in opt_methods:
            action = QAction(label, self)
            action.setCheckable(True)
            try:
                action.setActionGroup(opt_group)
            except Exception:
                pass
            action.triggered.connect(lambda checked, m=key: self.set_optimization_method(m))
            optimization_menu.addAction(action)
            opt_group.addAction(action)
            opt_actions[key] = action

        # Persist the actions mapping so other methods can update the checked state
        self.opt3d_actions = opt_actions

        # Determine the initial checked menu item from saved settings (fall back to MMFF_RDKIT)
        try:
            saved_opt = (self.settings.get('optimization_method') or self.optimization_method or 'MMFF_RDKIT').upper()
        except Exception:
            saved_opt = 'MMFF_RDKIT'

        try:
            if saved_opt in self.opt3d_actions and self.opt3d_actions[saved_opt].isEnabled():
                self.opt3d_actions[saved_opt].setChecked(True)
                self.optimization_method = saved_opt
            else:
                if 'MMFF_RDKIT' in self.opt3d_actions:
                    self.opt3d_actions['MMFF_RDKIT'].setChecked(True)
                    self.optimization_method = 'MMFF_RDKIT'
        except Exception:
            pass
    
        # 4) Reset all settings to defaults
        settings_menu.addSeparator()
        reset_settings_action = QAction("Reset All Settings", self)
        reset_settings_action.triggered.connect(self.reset_all_settings_menu)
        settings_menu.addAction(reset_settings_action)

        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(about_action)

        github_action = QAction("GitHub", self)
        github_action.triggered.connect(
            lambda: QDesktopServices.openUrl(QUrl("https://github.com/HiroYokoyama/python_molecular_editor"))
        )
        help_menu.addAction(github_action)

        github_wiki_action = QAction("GitHub Wiki", self)
        github_wiki_action.triggered.connect(
            lambda: QDesktopServices.openUrl(QUrl("https://github.com/HiroYokoyama/python_molecular_editor/wiki"))
        )
        help_menu.addAction(github_wiki_action)

        manual_action = QAction("User Manual", self)
        manual_action.triggered.connect(
            lambda: QDesktopServices.openUrl(QUrl("https://hiroyokoyama.github.io/python_molecular_editor/manual/manual"))
        )
        help_menu.addAction(manual_action)

        # 3D関連機能の初期状態を統一的に設定
        self._enable_3d_features(False)
        
    def init_worker_thread(self):
        # Initialize shared state for calculation runs.
        # NOTE: we no longer create a persistent worker/thread here. Instead,
        # each conversion run will create its own CalculationWorker + QThread
        # so multiple conversions may run in parallel.
        # Shared halt id set used to request early termination of specific worker runs
        self.halt_ids = set()
        # IDs used to correlate start/halt/finish
        self.next_conversion_id = 1
        # Track currently-active conversion worker IDs so Halt can target all
        # running conversions. Use a set because multiple conversions may run
        # concurrently.
        self.active_worker_ids = set()
        # Track active threads for diagnostics/cleanup (weak references ok)
        try:
            self._active_calc_threads = []
        except Exception:
            self._active_calc_threads = []


    def update_status_bar(self, message):
        """ワーカースレッドからのメッセージでステータスバーを更新するスロット"""
        self.statusBar().showMessage(message)

    def set_mode(self, mode_str):
        prev_mode = getattr(self.scene, 'mode', None)
        self.scene.mode = mode_str
        self.view_2d.setMouseTracking(True)
        # テンプレートモードから離れる場合はゴーストを消す
        if prev_mode and prev_mode.startswith('template') and not mode_str.startswith('template'):
            self.scene.clear_template_preview()
        elif not mode_str.startswith('template'):
            self.scene.template_preview.hide()

        # カーソル形状の設定
        if mode_str == 'select':
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)
        elif mode_str.startswith(('atom', 'bond', 'template')):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        elif mode_str.startswith(('charge', 'radical')):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)

        if mode_str.startswith('atom'): 
            self.scene.current_atom_symbol = mode_str.split('_')[1]
            self.statusBar().showMessage(f"Mode: Draw Atom ({self.scene.current_atom_symbol})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True) 
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0
        elif mode_str.startswith('bond'):
            self.scene.current_atom_symbol = 'C'
            parts = mode_str.split('_')
            self.scene.bond_order = int(parts[1])
            self.scene.bond_stereo = int(parts[2]) if len(parts) > 2 else 0
            stereo_text = {0: "", 1: " (Wedge)", 2: " (Dash)"}.get(self.scene.bond_stereo, "")
            self.statusBar().showMessage(f"Mode: Draw Bond (Order: {self.scene.bond_order}{stereo_text})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True)
        elif mode_str.startswith('template'):
            if mode_str.startswith('template_user'):
                # User template mode
                template_name = mode_str.replace('template_user_', '')
                self.statusBar().showMessage(f"Mode: User Template ({template_name})")
            else:
                # Built-in template mode
                self.statusBar().showMessage(f"Mode: {mode_str.split('_')[1].capitalize()} Template")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'charge_plus':
            self.statusBar().showMessage("Mode: Increase Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'charge_minus':
            self.statusBar().showMessage("Mode: Decrease Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'radical':
            self.statusBar().showMessage("Mode: Toggle Radical (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)

        else: # Select mode
            self.statusBar().showMessage("Mode: Select")
            self.view_2d.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0

    def set_mode_and_update_toolbar(self, mode_str):
        self.set_mode(mode_str)
        # QAction→QToolButtonのマッピングを取得
        toolbar = getattr(self, 'toolbar', None)
        action_to_button = {}
        if toolbar:
            for key, action in self.mode_actions.items():
                btn = toolbar.widgetForAction(action)
                if btn:
                    action_to_button[action] = btn

        # すべてのモードボタンの選択解除＆色リセット
        for key, action in self.mode_actions.items():
            action.setChecked(False)
            btn = action_to_button.get(action)
            if btn:
                btn.setStyleSheet("")

        # テンプレート系（User含む）は全て同じスタイル適用
        if mode_str in self.mode_actions:
            action = self.mode_actions[mode_str]
            action.setChecked(True)
            btn = action_to_button.get(action)
            if btn:
                # テンプレート系は青、それ以外はクリア
                if mode_str.startswith('template'):
                    btn.setStyleSheet("background-color: #2196F3; color: white;")
                else:
                    btn.setStyleSheet("")

    def set_3d_style(self, style_name):
        """3D表示スタイルを設定し、ビューを更新する"""
        if self.current_3d_style == style_name:
            return

        # 描画モード変更時に測定モードと3D編集モードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化
        
        # 3D原子選択をクリア
        self.clear_3d_selection()

        self.current_3d_style = style_name
        self.statusBar().showMessage(f"3D style set to: {style_name}")
        
        # 現在表示中の分子があれば、新しいスタイルで再描画する
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)

    def set_optimization_method(self, method_name):
        """Set preferred 3D optimization method and persist to settings.

        Supported values: 'GAFF', 'MMFF'
        """
        # Normalize input and validate
        if not method_name:
            return
        method = str(method_name).strip().upper()
        valid_methods = (
            'MMFF_RDKIT', 'MMFF94_RDKIT', 'UFF_RDKIT',
            'UFF_OBABEL', 'GAFF_OBABEL', 'MMFF94_OBABEL', 'GHEMICAL_OBABEL'
        )
        if method not in valid_methods:
            # Unknown method: ignore but notify
            self.statusBar().showMessage(f"Unknown 3D optimization method: {method_name}")
            return

        # Update internal state (store canonical uppercase key)
        self.optimization_method = method

        # Persist to settings
        try:
                self.settings['optimization_method'] = self.optimization_method
                try:
                    self.settings_dirty = True
                except Exception:
                    pass
        except Exception:
            pass

        # Update menu checked state if actions mapping exists
        try:
            if hasattr(self, 'opt3d_actions') and self.opt3d_actions:
                for k, act in self.opt3d_actions.items():
                    try:
                        # keys in opt3d_actions may be mixed-case; compare uppercased
                        act.setChecked(k.upper() == method)
                    except Exception:
                        pass
        except Exception:
            pass

        # Also show user-friendly label if available
        try:
            label = self.opt3d_method_labels.get(self.optimization_method, self.optimization_method)
        except Exception:
            label = self.optimization_method
        self.statusBar().showMessage(f"3D optimization method set to: {label}")

    def copy_selection(self):
        """選択された原子と結合をクリップボードにコピーする"""
        try:
            selected_atoms = [item for item in self.scene.selectedItems() if isinstance(item, AtomItem)]
            if not selected_atoms:
                return

            # 選択された原子のIDセットを作成
            selected_atom_ids = {atom.atom_id for atom in selected_atoms}
            
            # 選択された原子の幾何学的中心を計算
            center = QPointF(
                sum(atom.pos().x() for atom in selected_atoms) / len(selected_atoms),
                sum(atom.pos().y() for atom in selected_atoms) / len(selected_atoms)
            )
            
            # コピー対象の原子データをリストに格納（位置は中心からの相対座標）
            # 同時に、元のatom_idから新しいインデックス(0, 1, 2...)へのマッピングを作成
            atom_id_to_idx_map = {}
            fragment_atoms = []
            for i, atom in enumerate(selected_atoms):
                atom_id_to_idx_map[atom.atom_id] = i
                fragment_atoms.append({
                    'symbol': atom.symbol,
                    'rel_pos': atom.pos() - center,
                    'charge': atom.charge,
                    'radical': atom.radical,
                })
                
            # 選択された原子同士を結ぶ結合のみをリストに格納
            fragment_bonds = []
            for (id1, id2), bond_data in self.data.bonds.items():
                if id1 in selected_atom_ids and id2 in selected_atom_ids:
                    fragment_bonds.append({
                        'idx1': atom_id_to_idx_map[id1],
                        'idx2': atom_id_to_idx_map[id2],
                        'order': bond_data['order'],
                        'stereo': bond_data.get('stereo', 0),  # E/Z立体化学情報も保存
                    })

            # pickleを使ってデータをバイト配列にシリアライズ
            data_to_pickle = {'atoms': fragment_atoms, 'bonds': fragment_bonds}
            byte_array = QByteArray()
            buffer = io.BytesIO()
            pickle.dump(data_to_pickle, buffer)
            byte_array.append(buffer.getvalue())

            # カスタムMIMEタイプでクリップボードに設定
            mime_data = QMimeData()
            mime_data.setData(CLIPBOARD_MIME_TYPE, byte_array)
            QApplication.clipboard().setMimeData(mime_data)
            self.statusBar().showMessage(f"Copied {len(fragment_atoms)} atoms and {len(fragment_bonds)} bonds.")
            
        except Exception as e:
            print(f"Error during copy operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during copy operation: {e}")

    def cut_selection(self):
        """選択されたアイテムを切り取り（コピーしてから削除）"""
        try:
            selected_items = self.scene.selectedItems()
            if not selected_items:
                return
            
            # 最初にコピー処理を実行
            self.copy_selection()
            
            if self.scene.delete_items(set(selected_items)):
                self.push_undo_state()
                self.statusBar().showMessage("Cut selection.", 2000)
                
        except Exception as e:
            print(f"Error during cut operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during cut operation: {e}")

    def paste_from_clipboard(self):
        """クリップボードから分子フラグメントを貼り付け"""
        try:
            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            if not mime_data.hasFormat(CLIPBOARD_MIME_TYPE):
                return

            byte_array = mime_data.data(CLIPBOARD_MIME_TYPE)
            buffer = io.BytesIO(byte_array)
            try:
                fragment_data = pickle.load(buffer)
            except pickle.UnpicklingError:
                self.statusBar().showMessage("Error: Invalid clipboard data format")
                return
            
            paste_center_pos = self.view_2d.mapToScene(self.view_2d.mapFromGlobal(QCursor.pos()))
            self.scene.clearSelection()

            new_atoms = []
            for atom_data in fragment_data['atoms']:
                pos = paste_center_pos + atom_data['rel_pos']
                new_id = self.scene.create_atom(
                    atom_data['symbol'], pos,
                    charge=atom_data.get('charge', 0),
                    radical=atom_data.get('radical', 0)
                )
                new_item = self.data.atoms[new_id]['item']
                new_atoms.append(new_item)
                new_item.setSelected(True)

            for bond_data in fragment_data['bonds']:
                atom1 = new_atoms[bond_data['idx1']]
                atom2 = new_atoms[bond_data['idx2']]
                self.scene.create_bond(
                    atom1, atom2,
                    bond_order=bond_data.get('order', 1),
                    bond_stereo=bond_data.get('stereo', 0)  # E/Z立体化学情報も復元
                )
            
            self.push_undo_state()
            self.statusBar().showMessage(f"Pasted {len(fragment_data['atoms'])} atoms and {len(fragment_data['bonds'])} bonds.", 2000)
            
        except Exception as e:
            print(f"Error during paste operation: {e}")
            
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during paste operation: {e}")
        self.statusBar().showMessage(f"Pasted {len(new_atoms)} atoms.", 2000)
        self.activate_select_mode()

    def remove_hydrogen_atoms(self):
        """2Dビューで水素原子とその結合を削除する"""
        try:
            # Collect hydrogen atom items robustly (store atom_id -> item)
            hydrogen_map = {}

            # Iterate over a snapshot of atoms to avoid "dictionary changed size"
            for atom_id, atom_data in list(self.data.atoms.items()):
                try:
                    if atom_data.get('symbol') != 'H':
                        continue
                    item = atom_data.get('item')
                    # Only collect live AtomItem wrappers
                    if item is None:
                        continue
                    if sip_isdeleted_safe(item):
                        continue
                    if not isinstance(item, AtomItem):
                        continue
                    # Prefer storing by original atom id to detect actual removals later
                    hydrogen_map[atom_id] = item
                except Exception:
                    # Ignore problematic entries and continue scanning
                    continue

            if not hydrogen_map:
                self.statusBar().showMessage("No hydrogen atoms found to remove.", 2000)
                return

            # To avoid blocking the UI or causing large, monolithic deletions that may
            # trigger internal re-entrancy issues, delete in batches and process UI events
            items = list(hydrogen_map.values())
            total = len(items)
            batch_size = 200  # tuned conservative batch size
            deleted_any = False

            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                batch = set()
                # Filter out items that are already deleted or invalid just before deletion
                for it in items[start:end]:
                    try:
                        if it is None:
                            continue
                        if sip_isdeleted_safe(it):
                            continue
                        if not isinstance(it, AtomItem):
                            continue
                        batch.add(it)
                    except Exception:
                        continue

                if not batch:
                    # Nothing valid to delete in this batch
                    continue

                try:
                    # scene.delete_items is expected to handle bond cleanup; call it per-batch
                    success = False
                    try:
                        success = bool(self.scene.delete_items(batch))
                    except Exception:
                        # If scene.delete_items raises for a batch, attempt a safe per-item fallback
                        success = False

                    if not success:
                        # Fallback: try deleting items one-by-one to isolate problematic items
                        for it in list(batch):
                            try:
                                # Use scene.delete_items for single-item as well
                                ok = bool(self.scene.delete_items({it}))
                                if ok:
                                    deleted_any = True
                            except Exception:
                                # If single deletion also fails, skip that item
                                continue
                    else:
                        deleted_any = True

                except Exception:
                    # Continue with next batch on unexpected errors
                    continue

                # Allow the GUI to process events between batches to remain responsive
                try:
                    QApplication.processEvents()
                except Exception:
                    pass

            # Determine how many hydrogens actually were removed by re-scanning data
            remaining_h = 0
            try:
                for _, atom_data in list(self.data.atoms.items()):
                    try:
                        if atom_data.get('symbol') == 'H':
                            remaining_h += 1
                    except Exception:
                        continue
            except Exception:
                remaining_h = 0

            removed_count = max(0, len(hydrogen_map) - remaining_h)

            if removed_count > 0:
                # Only push a single undo state once for the whole operation
                try:
                    self.push_undo_state()
                except Exception:
                    # Do not allow undo stack problems to crash the app
                    pass
                self.statusBar().showMessage(f"Removed {removed_count} hydrogen atoms.", 2000)
            else:
                # If nothing removed but we attempted, show an informative message
                if deleted_any:
                    # Deleted something but couldn't determine count reliably
                    self.statusBar().showMessage("Removed hydrogen atoms (count unknown).", 2000)
                else:
                    self.statusBar().showMessage("Failed to remove hydrogen atoms or none found.")

        except Exception as e:
            # Capture and log unexpected errors but don't let them crash the UI
            print(f"Error during hydrogen removal: {e}")
            traceback.print_exc()
            try:
                self.statusBar().showMessage(f"Error removing hydrogen atoms: {e}")
            except Exception:
                pass

    def add_hydrogen_atoms(self):
        """RDKitで各原子の暗黙の水素数を調べ、その数だけ明示的な水素原子と単結合を作成する（2Dビュー）。

        実装上の仮定:
        - `self.data.to_rdkit_mol()` は各RDKit原子に `_original_atom_id` プロパティを設定している。
        - 原子の2D座標は `self.data.atoms[orig_id]['item'].pos()` で得られる。
        - 新しい原子は `self.scene.create_atom(symbol, pos, ...)` で追加し、
          結合は `self.scene.create_bond(atom_item, hydrogen_item, bond_order=1)` で作成する。
        """
        try:
            

            mol = self.data.to_rdkit_mol(use_2d_stereo=False)
            if not mol or mol.GetNumAtoms() == 0:
                self.statusBar().showMessage("No molecule available to compute hydrogens.", 2000)
                return

            added_count = 0
            added_items = []

            # すべてのRDKit原子について暗黙水素数を確認
            for idx in range(mol.GetNumAtoms()):
                rd_atom = mol.GetAtomWithIdx(idx)
                try:
                    orig_id = rd_atom.GetIntProp("_original_atom_id")
                except Exception:
                    # 元のエディタ側のIDがない場合はスキップ
                    continue

                if orig_id not in self.data.atoms:
                    continue

                # 暗黙水素数を優先して取得。存在しない場合は総水素数 - 明示水素数を使用
                implicit_h = int(rd_atom.GetNumImplicitHs()) if hasattr(rd_atom, 'GetNumImplicitHs') else 0
                if implicit_h is None or implicit_h < 0:
                    implicit_h = 0
                if implicit_h == 0:
                    # フォールバック
                    try:
                        total_h = int(rd_atom.GetTotalNumHs())
                        explicit_h = int(rd_atom.GetNumExplicitHs()) if hasattr(rd_atom, 'GetNumExplicitHs') else 0
                        implicit_h = max(0, total_h - explicit_h)
                    except Exception:
                        implicit_h = 0

                if implicit_h <= 0:
                    continue

                parent_item = self.data.atoms[orig_id]['item']
                parent_pos = parent_item.pos()

                # 周囲の近接原子の方向を取得して、水素を邪魔しないように角度を決定
                neighbor_angles = []
                try:
                    for (a1, a2), bdata in self.data.bonds.items():
                        # 対象原子に結合している近傍の原子角度を収集する。
                        # ただし既存の水素は配置に影響させない（すでにあるHで埋めない）。
                        try:
                            if a1 == orig_id and a2 in self.data.atoms:
                                neigh = self.data.atoms[a2]
                                if neigh.get('symbol') == 'H':
                                    continue
                                if neigh.get('item') is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get('item')):
                                    continue
                                vec = neigh['item'].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                            elif a2 == orig_id and a1 in self.data.atoms:
                                neigh = self.data.atoms[a1]
                                if neigh.get('symbol') == 'H':
                                    continue
                                if neigh.get('item') is None:
                                    continue
                                if sip_isdeleted_safe(neigh.get('item')):
                                    continue
                                vec = neigh['item'].pos() - parent_pos
                                neighbor_angles.append(math.atan2(vec.y(), vec.x()))
                        except Exception:
                            # 個々の近傍読み取りの問題は無視して続行
                            continue
                except Exception:
                    neighbor_angles = []

                # 画面上の適当な結合長（ピクセル）を使用
                bond_length = 75

                # ヘルパー: 指定インデックスの水素に使うbond_stereoを決定
                def _choose_stereo(i):
                    # 0: plain, 1: wedge, 2: dash, 3: plain, 4+: all plain
                    if i == 0:
                        return 0
                    if i == 1:
                        return 1
                    if i == 2:
                        return 2
                    return 0  #4th+ hydrogens are all plain

                # 角度配置を改善: 既存の結合角度の最大ギャップを見つけ、
                # そこに水素を均等配置する。既存結合が無ければ全周に均等配置。
                target_angles = []
                try:
                    if not neighbor_angles:
                        # 既存結合が無い -> 全円周に均等配置
                        for h_idx in range(implicit_h):
                            angle = (2.0 * math.pi * h_idx) / implicit_h
                            target_angles.append(angle)
                    else:
                        # 正規化してソート
                        angs = [((a + 2.0 * math.pi) if a < 0 else a) for a in neighbor_angles]
                        angs = sorted(angs)
                        # ギャップを計算（循環含む）
                        gaps = []  # list of (gap_size, start_angle, end_angle)
                        for i in range(len(angs)):
                            a1 = angs[i]
                            a2 = angs[(i + 1) % len(angs)]
                            if i == len(angs) - 1:
                                # wrap-around gap
                                gap = (a2 + 2.0 * math.pi) - a1
                                start = a1
                                end = a2 + 2.0 * math.pi
                            else:
                                gap = a2 - a1
                                start = a1
                                end = a2
                            gaps.append((gap, start, end))

                        # 最大ギャップを選ぶ
                        gaps.sort(key=lambda x: x[0], reverse=True)
                        max_gap, gstart, gend = gaps[0]
                        # もし最大ギャップが小さい（つまり周りに均等に原子がある）でも
                        # そのギャップ内に均等配置することで既存結合と重ならないようにする
                        # ギャップ内に implicit_h 個を等間隔で配置（分割数 = implicit_h + 1）
                        for i in range(implicit_h):
                            seg = max_gap / (implicit_h + 1)
                            angle = gstart + (i + 1) * seg
                            # 折り返しを戻して 0..2pi に正規化
                            angle = angle % (2.0 * math.pi)
                            target_angles.append(angle)
                except Exception:
                    # フォールバック: 単純な等間隔配置
                    for h_idx in range(implicit_h):
                        angle = (2.0 * math.pi * h_idx) / implicit_h
                        target_angles.append(angle)

                # 角度から位置を計算して原子と結合を追加
                for h_idx, angle in enumerate(target_angles):
                    dx = bond_length * math.cos(angle)
                    dy = bond_length * math.sin(angle)
                    pos = QPointF(parent_pos.x() + dx, parent_pos.y() + dy)

                    # 新しい水素原子を作成
                    try:
                        new_id = self.scene.create_atom('H', pos)
                        new_item = self.data.atoms[new_id]['item']
                        # bond_stereo を指定（最初は plain=0, 次に wedge/dash）
                        stereo = _choose_stereo(h_idx)
                        self.scene.create_bond(parent_item, new_item, bond_order=1, bond_stereo=stereo)
                        added_items.append(new_item)
                        added_count += 1
                    except Exception as e:
                        # 個々の追加失敗はログに残して続行
                        print(f"Failed to add H for atom {orig_id}: {e}")

            if added_count > 0:
                self.push_undo_state()
                self.statusBar().showMessage(f"Added {added_count} hydrogen atoms.", 2000)
                # 選択を有効化して追加した原子を選択状態にする
                try:
                    self.scene.clearSelection()
                    for it in added_items:
                        it.setSelected(True)
                except Exception:
                    pass
            else:
                self.statusBar().showMessage("No implicit hydrogens found to add.", 2000)

        except Exception as e:
            print(f"Error during hydrogen addition: {e}")
            traceback.print_exc()
            self.statusBar().showMessage(f"Error adding hydrogen atoms: {e}")

    def update_edit_menu_actions(self):
        """選択状態やクリップボードの状態に応じて編集メニューを更新"""
        try:
            has_selection = len(self.scene.selectedItems()) > 0
            self.cut_action.setEnabled(has_selection)
            self.copy_action.setEnabled(has_selection)
            
            clipboard = QApplication.clipboard()
            mime_data = clipboard.mimeData()
            self.paste_action.setEnabled(mime_data is not None and mime_data.hasFormat(CLIPBOARD_MIME_TYPE))
        except RuntimeError:
            pass


    def show_convert_menu(self, pos):
        """右クリックで表示する一時的な3D変換メニュー。
        選択したモードは一時フラグとして保持され、その後の変換で使用されます（永続化しません）。
        """
        try:
            menu = QMenu(self)
            conv_options = [
                ("RDKit -> Open Babel (fallback)", 'fallback'),
                ("RDKit only", 'rdkit'),
                ("Open Babel only", 'obabel'),
                ("Direct (use 2D coords + add H)", 'direct')
            ]
            for label, key in conv_options:
                a = QAction(label, self)
                # If Open Babel is not available, disable actions that depend on it
                if key in ('obabel', 'fallback') and not globals().get('OBABEL_AVAILABLE', False):
                    a.setEnabled(False)
                a.triggered.connect(lambda checked=False, k=key: self._trigger_conversion_with_temp_mode(k))
                menu.addAction(a)

            # Show menu at button position
            menu.exec_(self.convert_button.mapToGlobal(pos))
        except Exception as e:
            print(f"Error showing convert menu: {e}")


    def activate_select_mode(self):
        self.set_mode('select')
        if 'select' in self.mode_actions:
            self.mode_actions['select'].setChecked(True)


    def _trigger_conversion_with_temp_mode(self, mode_key):
        try:
            # store temporary override and invoke conversion
            self._temp_conv_mode = mode_key
            # Call the normal conversion entry point (it will consume the temp)
            QTimer.singleShot(0, self.trigger_conversion)
        except Exception as e:
            print(f"Failed to start conversion with temp mode {mode_key}: {e}")


    def show_optimize_menu(self, pos):
        """右クリックで表示する一時的な3D最適化メニュー。
        選択したメソッドは一時フラグとして保持され、その後の最適化で使用されます（永続化しません）。
        """
        try:
            menu = QMenu(self)
            opt_list = [
                ("MMFF94s", 'MMFF_RDKIT'),
                ("MMFF94", 'MMFF94_RDKIT'),
                ("UFF", 'UFF_RDKIT')
            ]
            for label, key in opt_list:
                a = QAction(label, self)
                # If opt3d_actions exist, reflect their enabled state
                try:
                    if hasattr(self, 'opt3d_actions') and key in self.opt3d_actions:
                        a.setEnabled(self.opt3d_actions[key].isEnabled())
                except Exception:
                    pass
                a.triggered.connect(lambda checked=False, k=key: self._trigger_optimize_with_temp_method(k))
                menu.addAction(a)

            menu.exec_(self.optimize_3d_button.mapToGlobal(pos))
        except Exception as e:
            print(f"Error showing optimize menu: {e}")


    def _trigger_optimize_with_temp_method(self, method_key):
        try:
            # store temporary override and invoke optimization
            self._temp_optimization_method = method_key
            # Run optimize on next event loop turn so UI updates first
            QTimer.singleShot(0, self.optimize_3d_structure)
        except Exception as e:
            print(f"Failed to start optimization with temp method {method_key}: {e}")

    def trigger_conversion(self):
        # Reset last successful optimization method at start of new conversion
        self.last_successful_optimization_method = None
        
        # 3D変換時に既存の3D制約をクリア
        self.constraints_3d = []

        # 2Dエディタに原子が存在しない場合は3Dビューをクリア
        if not self.data.atoms:
            self.plotter.clear()
            self.current_mol = None
            self.analysis_action.setEnabled(False)
            self.statusBar().showMessage("3D view cleared.")
            self.view_2d.setFocus() 
            return

        # 描画モード変更時に測定モードと3D編集モードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化

        mol = self.data.to_rdkit_mol(use_2d_stereo=False)

        # 分子オブジェクトが作成できない場合でも化学的問題をチェック
        if not mol or mol.GetNumAtoms() == 0:
            # RDKitでの変換に失敗した場合は、独自の化学的問題チェックを実行
            self.check_chemistry_problems_fallback()
            return

        # 原子プロパティを保存（ワーカープロセスで失われるため）
        self.original_atom_properties = {}
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            try:
                original_id = atom.GetIntProp("_original_atom_id")
                self.original_atom_properties[i] = original_id
            except KeyError:
                pass

        problems = Chem.DetectChemistryProblems(mol)
        if problems:
            # 化学的問題が見つかった場合は既存のフラグをクリアしてから新しい問題を表示
            self.scene.clear_all_problem_flags()
            self.statusBar().showMessage(f"Error: {len(problems)} chemistry problem(s) found.")
            # 既存の選択状態をクリア
            self.scene.clearSelection() 

            # 問題のある原子に赤枠フラグを立てる
            for prob in problems:
                atom_idx = prob.GetAtomIdx()
                rdkit_atom = mol.GetAtomWithIdx(atom_idx)
                # エディタ側での原子IDの取得と存在確認
                if rdkit_atom.HasProp("_original_atom_id"):
                    original_id = rdkit_atom.GetIntProp("_original_atom_id")
                    if original_id in self.data.atoms and self.data.atoms[original_id]['item']:
                        item = self.data.atoms[original_id]['item']
                        item.has_problem = True 
                        item.update()

            self.view_2d.setFocus()
            return

        # 化学的問題がない場合のみフラグをクリアして3D変換を実行
        self.scene.clear_all_problem_flags()

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus() 
            return

        # 複数分子の処理に対応
        num_frags = len(Chem.GetMolFrags(mol))
        if num_frags > 1:
            self.statusBar().showMessage(f"Converting {num_frags} molecules to 3D with collision detection...")
        else:
            self.statusBar().showMessage("Calculating 3D structure...")
            
        # CRITICAL FIX: Use the 2D editor's MOL block instead of RDKit's to preserve
        # wedge/dash stereo information that is stored in the 2D editor data.
        # RDKit's MolToMolBlock() doesn't preserve this information.
        mol_block = self.data.to_mol_block()
        if not mol_block:
            mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
        
        # Additional E/Z stereo enhancement: add M CFG lines for explicit E/Z bonds
        mol_lines = mol_block.split('\n')
        
        # Find bonds with explicit E/Z labels from our data and map to RDKit bond indices
        ez_bond_info = {}
        for (id1, id2), bond_data in self.data.bonds.items():
            if bond_data.get('stereo') in [3, 4]:  # E/Z labels
                # Find corresponding atoms in RDKit molecule by _original_atom_id property
                rdkit_idx1 = None
                rdkit_idx2 = None
                for atom in mol.GetAtoms():
                    if atom.HasProp("_original_atom_id"):
                        orig_id = atom.GetIntProp("_original_atom_id")
                        if orig_id == id1:
                            rdkit_idx1 = atom.GetIdx()
                        elif orig_id == id2:
                            rdkit_idx2 = atom.GetIdx()
                
                if rdkit_idx1 is not None and rdkit_idx2 is not None:
                    rdkit_bond = mol.GetBondBetweenAtoms(rdkit_idx1, rdkit_idx2)
                    if rdkit_bond and rdkit_bond.GetBondType() == Chem.BondType.DOUBLE:
                        ez_bond_info[rdkit_bond.GetIdx()] = bond_data['stereo']
        
        # Add M  CFG lines for E/Z stereo if needed
        if ez_bond_info:
            insert_idx = len(mol_lines) - 1  # Before M  END
            for bond_idx, stereo_type in ez_bond_info.items():
                cfg_value = 1 if stereo_type == 3 else 2  # 1=Z, 2=E in MOL format
                cfg_line = f"M  CFG  1 {bond_idx + 1:3d}   {cfg_value}"
                mol_lines.insert(insert_idx, cfg_line)
                insert_idx += 1
            mol_block = '\n'.join(mol_lines)
        
        # Assign a unique ID for this conversion run so it can be halted/validated
        try:
            run_id = int(self.next_conversion_id)
        except Exception:
            run_id = 1
        try:
            self.next_conversion_id = run_id + 1
        except Exception:
            self.next_conversion_id = getattr(self, 'next_conversion_id', 1) + 1

        # Record this run as active. Use a set to track all active worker ids
        # so a Halt request can target every running conversion.
        try:
            self.active_worker_ids.add(run_id)
        except Exception:
            # Ensure attribute exists in case of weird states
            self.active_worker_ids = set([run_id])

        # Change the convert button to a Halt button so user can cancel
        try:
            # keep it enabled so the user can click Halt
            self.convert_button.setText("Halt conversion")
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.clicked.connect(self.halt_conversion)
        except Exception:
            pass

        # Keep cleanup disabled while conversion is in progress
        self.cleanup_button.setEnabled(False)
        # Disable 3D features during calculation
        self._enable_3d_features(False)
        self.statusBar().showMessage("Calculating 3D structure...")
        self.plotter.clear() 
        bg_color_hex = self.settings.get('background_color', '#919191')
        bg_qcolor = QColor(bg_color_hex)
        
        if bg_qcolor.isValid():
            luminance = bg_qcolor.toHsl().lightness()
            text_color = 'black' if luminance > 128 else 'white'
        else:
            text_color = 'white'
        
        text_actor = self.plotter.add_text(
            "Calculating...",
            position='lower_right',
            font_size=15,
            color=text_color,
            name='calculating_text'
        )
        # Keep a reference so we can reliably remove the text actor later
        try:
            self._calculating_text_actor = text_actor
        except Exception:
            # Best-effort: if storing fails, ignore — cleanup will still attempt renderer removal
            pass
        text_actor.GetTextProperty().SetOpacity(1)
        self.plotter.render()
        # Emit skip flag so the worker can ignore sanitization errors if user requested
        # Determine conversion_mode from settings (default: 'fallback').
        # If the user invoked conversion via the right-click menu, a temporary
        # override may be set on self._temp_conv_mode and should be used once.
        conv_mode = getattr(self, '_temp_conv_mode', None)
        if conv_mode:
            try:
                del self._temp_conv_mode
            except Exception:
                try:
                    delattr(self, '_temp_conv_mode')
                except Exception:
                    pass
        else:
            conv_mode = self.settings.get('3d_conversion_mode', 'fallback')

        # Allow a temporary optimization method override as well (used when
        # Optimize 3D is invoked via right-click menu). Do not persist here.
        opt_method = getattr(self, '_temp_optimization_method', None) or self.optimization_method
        if hasattr(self, '_temp_optimization_method'):
            try:
                del self._temp_optimization_method
            except Exception:
                try:
                    delattr(self, '_temp_optimization_method')
                except Exception:
                    pass

        options = {'conversion_mode': conv_mode, 'optimization_method': opt_method}
        # Attach the run id so the worker and main thread can correlate
        try:
            # Attach the concrete run id rather than the single waiting id
            options['worker_id'] = run_id
        except Exception:
            pass

        # Create a fresh CalculationWorker + QThread for this run so multiple
        # conversions can execute in parallel. The worker will be cleaned up
        # automatically after it finishes/errors.
        try:
            thread = QThread()
            worker = CalculationWorker()
            # Share the halt_ids set so user can request cancellation
            try:
                worker.halt_ids = self.halt_ids
            except Exception:
                pass

            worker.moveToThread(thread)

            # Forward status signals to main window handlers
            try:
                worker.status_update.connect(self.update_status_bar)
            except Exception:
                pass

            # When the worker finishes, call existing handler and then clean up
            def _on_worker_finished(result, w=worker, t=thread):
                try:
                    # deliver result to existing handler
                    self.on_calculation_finished(result)
                finally:
                    # Clean up signal connections to avoid stale references
                    # worker used its own start_work signal; no shared-signal
                    # disconnect necessary here.
                    # Remove thread from active threads list
                    try:
                        self._active_calc_threads.remove(t)
                    except Exception:
                        pass
                    try:
                        # ask thread to quit; it will finish as worker returns
                        t.quit()
                    except Exception:
                        pass
                    try:
                        # ensure thread object is deleted when finished
                        t.finished.connect(t.deleteLater)
                    except Exception:
                        pass
                    try:
                        # schedule worker deletion
                        w.deleteLater()
                    except Exception:
                        pass

            # When the worker errors (or halts), call existing handler and then clean up
            def _on_worker_error(error_msg, w=worker, t=thread):
                try:
                    # deliver error to existing handler
                    self.on_calculation_error(error_msg)
                finally:
                    # Clean up signal connections to avoid stale references
                    # worker used its own start_work signal; no shared-signal
                    # disconnect necessary here.
                    # Remove thread from active threads list
                    try:
                        self._active_calc_threads.remove(t)
                    except Exception:
                        pass
                    try:
                        # ask thread to quit; it will finish as worker returns
                        t.quit()
                    except Exception:
                        pass
                    try:
                        # ensure thread object is deleted when finished
                        t.finished.connect(t.deleteLater)
                    except Exception:
                        pass
                    try:
                        # schedule worker deletion
                        w.deleteLater()
                    except Exception:
                        pass

            try:
                worker.error.connect(_on_worker_error)
            except Exception:
                pass

            try:
                worker.finished.connect(_on_worker_finished)
            except Exception:
                pass

            # Start the thread
            thread.start()

            # Start the worker calculation via the worker's own start_work signal
            # (queued to the worker thread). Capture variables into lambda defaults
            # to avoid late-binding issues.
            QTimer.singleShot(10, lambda w=worker, m=mol_block, o=options: w.start_work.emit(m, o))

            # Track the thread so it isn't immediately garbage-collected (diagnostics)
            try:
                self._active_calc_threads.append(thread)
            except Exception:
                pass
        except Exception as e:
            # Fall back: if thread/worker creation failed, create a local
            # worker and start it (runs in main thread). This preserves
            # functionality without relying on the shared MainWindow signal.
            try:
                fallback_worker = CalculationWorker()
                QTimer.singleShot(10, lambda w=fallback_worker, m=mol_block, o=options: w.start_work.emit(m, o))
            except Exception:
                # surface the original error via existing UI path
                self.on_calculation_error(str(e))

        # 状態をUndo履歴に保存
        self.push_undo_state()
        self.update_chiral_labels()
        
        self.view_2d.setFocus()

    def halt_conversion(self):
        """User requested to halt the in-progress conversion.

        This will mark the current waiting_worker_id as halted (added to halt_ids),
        clear the waiting_worker_id, and immediately restore the UI (button text
        and handlers). The worker thread will observe halt_ids and should stop.
        """
        try:
            # Halt all currently-active workers by adding their ids to halt_ids
            wids_to_halt = set(getattr(self, 'active_worker_ids', set()))
            if wids_to_halt:
                try:
                    self.halt_ids.update(wids_to_halt)
                except Exception:
                    pass

            # Clear the active set immediately so UI reflects cancellation
            try:
                if hasattr(self, 'active_worker_ids'):
                    self.active_worker_ids.clear()
            except Exception:
                pass

            # Restore UI immediately
            try:
                try:
                    self.convert_button.clicked.disconnect()
                except Exception:
                    pass
                self.convert_button.setText("Convert 2D to 3D")
                self.convert_button.clicked.connect(self.trigger_conversion)
                self.convert_button.setEnabled(True)
            except Exception:
                pass

            try:
                self.cleanup_button.setEnabled(True)
            except Exception:
                pass

            # Remove any calculating text actor if present
            try:
                actor = getattr(self, '_calculating_text_actor', None)
                if actor is not None:
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
            except Exception:
                pass

            # Give immediate feedback
            self.statusBar().showMessage("3D conversion halted. Waiting for the thread to finish")
        except Exception:
            pass

    def check_chemistry_problems_fallback(self):
        """RDKit変換が失敗した場合の化学的問題チェック（独自実装）"""
        try:
            # 既存のフラグをクリア
            self.scene.clear_all_problem_flags()
            
            # 簡易的な化学的問題チェック
            problem_atoms = []
            
            for atom_id, atom_data in self.data.atoms.items():
                atom_item = atom_data.get('item')
                if not atom_item:
                    continue
                
                symbol = atom_data['symbol']
                charge = atom_data.get('charge', 0)
                
                # 結合数を計算
                bond_count = 0
                for (id1, id2), bond_data in self.data.bonds.items():
                    if id1 == atom_id or id2 == atom_id:
                        bond_count += bond_data.get('order', 1)
                
                # 基本的な価数チェック
                is_problematic = False
                if symbol == 'C' and bond_count > 4:
                    is_problematic = True
                elif symbol == 'N' and bond_count > 3 and charge == 0:
                    is_problematic = True
                elif symbol == 'O' and bond_count > 2 and charge == 0:
                    is_problematic = True
                elif symbol == 'H' and bond_count > 1:
                    is_problematic = True
                elif symbol in ['F', 'Cl', 'Br', 'I'] and bond_count > 1 and charge == 0:
                    is_problematic = True
                
                if is_problematic:
                    problem_atoms.append(atom_item)
            
            if problem_atoms:
                # 問題のある原子に赤枠を設定
                for atom_item in problem_atoms:
                    atom_item.has_problem = True
                    atom_item.update()
                
                self.statusBar().showMessage(f"Error: {len(problem_atoms)} chemistry problem(s) found (valence issues).")
            else:
                self.statusBar().showMessage("Error: Invalid chemical structure (RDKit conversion failed).")
            
            self.scene.clearSelection()
            self.view_2d.setFocus()
            
        except Exception as e:
            print(f"Error in fallback chemistry check: {e}")
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus()

    def optimize_3d_structure(self):
        """現在の3D分子構造を力場で最適化する"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to optimize.")
            return

        # If a prior chemical/sanitization check was attempted and failed, do not run optimization
        if getattr(self, 'chem_check_tried', False) and getattr(self, 'chem_check_failed', False):
            self.statusBar().showMessage("3D optimization disabled: molecule failed chemical sanitization.")
            # Ensure the Optimize 3D button is disabled to reflect this
            if hasattr(self, 'optimize_3d_button'):
                try:
                    self.optimize_3d_button.setEnabled(False)
                except Exception:
                    pass
            return

        self.statusBar().showMessage("Optimizing 3D structure...")
        QApplication.processEvents() # UIの更新を確実に行う

        try:
            # Allow a temporary optimization method override (right-click menu)
            method = getattr(self, '_temp_optimization_method', None) or getattr(self, 'optimization_method', 'MMFF_RDKIT')
            # Clear temporary override if present
            if hasattr(self, '_temp_optimization_method'):
                try:
                    del self._temp_optimization_method
                except Exception:
                    try:
                        delattr(self, '_temp_optimization_method')
                    except Exception:
                        pass
            method = method.upper() if method else 'MMFF_RDKIT'
            # 事前チェック：コンフォーマがあるか
            if self.current_mol.GetNumConformers() == 0:
                self.statusBar().showMessage("No conformer found: cannot optimize. Embed molecule first.")
                return
            if method in ('MMFF_RDKIT', 'MMFF94_RDKIT'):
                try:
                    # Choose concrete mmffVariant string
                    mmff_variant = "MMFF94s" if method == 'MMFF_RDKIT' else "MMFF94"
                    res = AllChem.MMFFOptimizeMolecule(self.current_mol, maxIters=4000, mmffVariant=mmff_variant)
                    if res != 0:
                        # 非収束や何らかの問題が起きた可能性 -> ForceField API で詳細に試す
                        try:
                            mmff_props = AllChem.MMFFGetMoleculeProperties(self.current_mol)
                            ff = AllChem.MMFFGetMoleculeForceField(self.current_mol, mmff_props, confId=0)
                            ff_ret = ff.Minimize(maxIts=4000)
                            if ff_ret != 0:
                                self.statusBar().showMessage(f"{mmff_variant} minimize returned non-zero status: {ff_ret}")
                                return
                        except Exception as e:
                            self.statusBar().showMessage(f"{mmff_variant} parameterization/minimize failed: {e}")
                            return
                except Exception as e:
                    self.statusBar().showMessage(f"{mmff_variant} (RDKit) optimization error: {e}")
                    return
            elif method == 'UFF_RDKIT':
                try:
                    res = AllChem.UFFOptimizeMolecule(self.current_mol, maxIters=4000)
                    if res != 0:
                        try:
                            ff = AllChem.UFFGetMoleculeForceField(self.current_mol, confId=0)
                            ff_ret = ff.Minimize(maxIts=4000)
                            if ff_ret != 0:
                                self.statusBar().showMessage(f"UFF minimize returned non-zero status: {ff_ret}")
                                return
                        except Exception as e:
                            self.statusBar().showMessage(f"UFF parameterization/minimize failed: {e}")
                            return
                except Exception as e:
                    self.statusBar().showMessage(f"UFF (RDKit) optimization error: {e}")
                    return
            else:
                self.statusBar().showMessage("Selected optimization method is not available. Use MMFF94 (RDKit) or UFF (RDKit).")
                return
        except Exception as e:
            self.statusBar().showMessage(f"3D optimization error: {e}")
        
        # 最適化後の構造で3Dビューを再描画
        try:
            # Remember which concrete optimizer variant succeeded so it
            # can be saved with the project. Normalize internal flags to
            # a human-friendly label: MMFF94s, MMFF94, or UFF.
            try:
                norm_method = None
                m = method.upper() if method else None
                if m in ('MMFF_RDKIT', 'MMFF94_RDKIT'):
                    # The code above uses mmffVariant="MMFF94s" when
                    # method == 'MMFF_RDKIT' and "MMFF94" otherwise.
                    norm_method = 'MMFF94s' if m == 'MMFF_RDKIT' else 'MMFF94'
                elif m == 'UFF_RDKIT' or m == 'UFF':
                    norm_method = 'UFF'
                else:
                    norm_method = getattr(self, 'optimization_method', None)

                # store for later serialization
                if norm_method:
                    self.last_successful_optimization_method = norm_method
            except Exception:
                pass
            # 3D最適化後は3D座標から立体化学を再計算（2回目以降は3D優先）
            if self.current_mol.GetNumConformers() > 0:
                Chem.AssignAtomChiralTagsFromStructure(self.current_mol, confId=0)
            self.update_chiral_labels() # キラル中心のラベルも更新
        except Exception:
            pass
            
        self.draw_molecule_3d(self.current_mol)
        
        # Show which method was used in the status bar (prefer human-readable label).
        # Prefer the actual method used during this run (last_successful_optimization_method
        # set earlier), then any temporary/local override used for this call (method),
        # and finally the persisted preference (self.optimization_method).
        try:
            used_method = (
                getattr(self, 'last_successful_optimization_method', None)
                or locals().get('method', None)
                or getattr(self, 'optimization_method', None)
            )
            used_label = None
            if used_method:
                # opt3d_method_labels keys are stored upper-case; normalize for lookup
                used_label = (getattr(self, 'opt3d_method_labels', {}) or {}).get(str(used_method).upper(), used_method)
        except Exception:
            used_label = None

        if used_label:
            self.statusBar().showMessage(f"3D structure optimization successful. Method: {used_label}")
        else:
            self.statusBar().showMessage("3D structure optimization successful.")
        self.push_undo_state() # Undo履歴に保存
        self.view_2d.setFocus()

    def on_calculation_finished(self, result):
        # Accept either (worker_id, mol) tuple or legacy single mol arg
        worker_id = None
        mol = None
        try:
            if isinstance(result, tuple) and len(result) == 2:
                worker_id, mol = result
            else:
                mol = result
        except Exception:
            mol = result

        # If this finished result is from a stale/halting run, discard it
        try:
            if worker_id is not None:
                # If this worker_id is not in the active set, it's stale/halting
                if worker_id not in getattr(self, 'active_worker_ids', set()):
                    # Cleanup calculating UI and ignore
                    try:
                        actor = getattr(self, '_calculating_text_actor', None)
                        if actor is not None:
                            if hasattr(self.plotter, 'remove_actor'):
                                try:
                                    self.plotter.remove_actor(actor)
                                except Exception:
                                    pass
                            else:
                                if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                    try:
                                        self.plotter.renderer.RemoveActor(actor)
                                    except Exception:
                                        pass
                            try:
                                delattr(self, '_calculating_text_actor')
                            except Exception:
                                try:
                                    del self._calculating_text_actor
                                except Exception:
                                    pass
                    except Exception:
                        pass
                    # Ensure Convert button is restored
                    try:
                        try:
                            self.convert_button.clicked.disconnect()
                        except Exception:
                            pass
                        self.convert_button.setText("Convert 2D to 3D")
                        self.convert_button.clicked.connect(self.trigger_conversion)
                        self.convert_button.setEnabled(True)
                    except Exception:
                        pass
                    try:
                        self.cleanup_button.setEnabled(True)
                    except Exception:
                        pass
                    self.statusBar().showMessage("Ignored result from stale conversion.")
                    return
        except Exception:
            pass

        # Remove the finished worker id from the active set and any halt set
        try:
            if worker_id is not None:
                try:
                    self.active_worker_ids.discard(worker_id)
                except Exception:
                    pass
            # Also remove id from halt set if present
            if worker_id is not None:
                try:
                    if worker_id in getattr(self, 'halt_ids', set()):
                        try:
                            self.halt_ids.discard(worker_id)
                        except Exception:
                            pass
                except Exception:
                    pass
        except Exception:
            pass

        self.dragged_atom_info = None
        self.current_mol = mol
        self.is_xyz_derived = False  # 2Dから生成した3D構造はXYZ由来ではない
        # Record the optimization method used for this conversion if available.
        try:
            opt_method = None
            try:
                # Worker or molecule may have attached a prop with the used method
                if hasattr(mol, 'HasProp') and mol is not None:
                    try:
                        if mol.HasProp('_pme_optimization_method'):
                            opt_method = mol.GetProp('_pme_optimization_method')
                    except Exception:
                        # not all Mol objects support HasProp/GetProp safely
                        pass
            except Exception:
                pass
            if not opt_method:
                opt_method = getattr(self, 'optimization_method', None)
            # normalize common forms
            if opt_method:
                om = str(opt_method).upper()
                if 'MMFF94S' in om or 'MMFF_RDKIT' in om:
                    self.last_successful_optimization_method = 'MMFF94s'
                elif 'MMFF94' in om:
                    self.last_successful_optimization_method = 'MMFF94'
                elif 'UFF' in om:
                    self.last_successful_optimization_method = 'UFF'
                else:
                    # store raw value otherwise
                    self.last_successful_optimization_method = opt_method
        except Exception:
            # non-fatal
            pass
        
        # 原子プロパティを復元（ワーカープロセスで失われたため）
        if hasattr(self, 'original_atom_properties'):
            for i, original_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(i)
                    atom.SetIntProp("_original_atom_id", original_id)
        
        # 原子IDマッピングを作成
        self.create_atom_id_mapping()
        
        # キラル中心を初回変換時は2Dの立体情報を考慮して設定
        try:
            if mol.GetNumConformers() > 0:
                # 初回変換では、2Dで設定したwedge/dashボンドの立体情報を保持
                # 立体化学の割り当てを行うが、既存の2D立体情報を尊重
                Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
            
            self.update_chiral_labels()
        except Exception:
            # 念のためエラーを握り潰して UI を壊さない
            pass

        self.draw_molecule_3d(mol)
        
        # 複数分子の場合、衝突検出と配置調整を実行
        try:
            frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
            if len(frags) > 1:
                self.statusBar().showMessage(f"Detecting collisions among {len(frags)} molecules...")
                QApplication.processEvents()
                self.adjust_molecule_positions_to_avoid_collisions(mol, frags)
                self.draw_molecule_3d(mol)
                self.update_chiral_labels()
                self.statusBar().showMessage(f"{len(frags)} molecules converted with collision avoidance.")
        except Exception as e:
            print(f"Warning: Collision detection failed: {e}")
            # 衝突検出に失敗してもエラーにはしない

        # Ensure any 'Calculating...' text is removed and the plotter is refreshed
        try:
            actor = getattr(self, '_calculating_text_actor', None)
            if actor is not None:
                try:
                    # Prefer plotter API if available
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            # Some pyvista versions use renderer.RemoveActor
                            if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                try:
                                    self.plotter.renderer.RemoveActor(actor)
                                except Exception:
                                    pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                finally:
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
            # Re-render to ensure the UI updates immediately
            try:
                self.plotter.render()
            except Exception:
                pass
        except Exception:
            pass

        #self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
        # Restore Convert button text/handler in case it was changed to Halt
        try:
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.setText("Convert 2D to 3D")
            self.convert_button.clicked.connect(self.trigger_conversion)
        except Exception:
            pass
        self.push_undo_state()
        self.view_2d.setFocus()
        self.cleanup_button.setEnabled(True)
        
        # 3D関連機能を統一的に有効化
        self._enable_3d_features(True)
            
        self.plotter.reset_camera()
        
        # 3D原子情報ホバー表示を再設定
        self.setup_3d_hover()
        
        # メニューテキストと状態を更新
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

    def create_atom_id_mapping(self):
        """2D原子IDから3D RDKit原子インデックスへのマッピングを作成する（RDKitの原子プロパティ使用）"""
        if not self.current_mol:
            return
            
        self.atom_id_to_rdkit_idx_map = {}
        
        # RDKitの原子プロパティから直接マッピングを作成
        for i in range(self.current_mol.GetNumAtoms()):
            rdkit_atom = self.current_mol.GetAtomWithIdx(i)
            try:
                original_atom_id = rdkit_atom.GetIntProp("_original_atom_id")
                self.atom_id_to_rdkit_idx_map[original_atom_id] = i
            except KeyError:
                # プロパティが設定されていない場合（外部ファイル読み込み時など）
                continue

    @pyqtSlot(object)
    def on_calculation_error(self, result):
        """ワーカースレッドからのエラー（またはHalt）を処理する"""
        worker_id = None
        error_message = ""
        try:
            if isinstance(result, tuple) and len(result) == 2:
                worker_id, error_message = result
            else:
                error_message = str(result)
        except Exception:
            error_message = str(result)

        # If this error is from a stale/previous worker (not in active set), ignore it.
        if worker_id is not None and worker_id not in getattr(self, 'active_worker_ids', set()):
            # Stale/late error from a previously-halted worker; ignore to avoid clobbering newer runs
            print(f"Ignored stale error from worker {worker_id}: {error_message}")
            return

        # Clear temporary plotter content and remove calculating text if present
        try:
            self.plotter.clear()
        except Exception:
            pass

        # Also attempt to explicitly remove the calculating text actor if it was stored
        try:
            actor = getattr(self, '_calculating_text_actor', None)
            if actor is not None:
                try:
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                try:
                                    self.plotter.renderer.RemoveActor(actor)
                                except Exception:
                                    pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                finally:
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
        except Exception:
            pass

        self.dragged_atom_info = None
        # Remove this worker id from active set (error belongs to this worker)
        try:
            if worker_id is not None:
                try:
                    self.active_worker_ids.discard(worker_id)
                except Exception:
                    pass
        except Exception:
            pass

        # If this error was caused by an intentional halt and the main thread
        # already cleared waiting_worker_id earlier for other reasons, suppress the error noise.
        try:
            low = (error_message or '').lower()
            # If a halt message and there are no active workers left, the user
            # already saw the halt message — suppress duplicate noise.
            if 'halt' in low and not getattr(self, 'active_worker_ids', set()):
                return
        except Exception:
            pass

        self.statusBar().showMessage(f"Error: {error_message}")
        
        try:
            self.cleanup_button.setEnabled(True)
        except Exception:
            pass
        try:
            # Restore Convert button text/handler
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.setText("Convert 2D to 3D")
            self.convert_button.clicked.connect(self.trigger_conversion)
            self.convert_button.setEnabled(True)
        except Exception:
            pass

        # On calculation error we should NOT enable 3D-only features.
        # Explicitly disable Optimize and Export so the user can't try to operate
        # on an invalid or missing 3D molecule.
        try:
            if hasattr(self, 'optimize_3d_button'):
                self.optimize_3d_button.setEnabled(False)
        except Exception:
            pass
        try:
            if hasattr(self, 'export_button'):
                self.export_button.setEnabled(False)
        except Exception:
            pass

        # Keep 3D feature buttons disabled to avoid inconsistent UI state
        try:
            self._enable_3d_features(False)
        except Exception:
            pass

        # Keep 3D edit actions disabled (no molecule to edit)
        try:
            self._enable_3d_edit_actions(False)
        except Exception:
            pass
        # Some menu items are explicitly disabled on error
        try:
            if hasattr(self, 'analysis_action'):
                self.analysis_action.setEnabled(False)
        except Exception:
            pass
        try:
            if hasattr(self, 'edit_3d_action'):
                self.edit_3d_action.setEnabled(False)
        except Exception:
            pass

        # Force a UI refresh
        try:
            self.plotter.render()
        except Exception:
            pass

        # Ensure focus returns to 2D editor
        self.view_2d.setFocus()

    def eventFilter(self, obj, event):
        if obj is self.plotter and event.type() == QEvent.Type.MouseButtonPress:
            self.view_2d.setFocus()
        return super().eventFilter(obj, event)

    def get_current_state(self):
        atoms = {atom_id: {'symbol': data['symbol'],
                           'pos': (data['item'].pos().x(), data['item'].pos().y()),
                           'charge': data.get('charge', 0),
                           'radical': data.get('radical', 0)} 
                 for atom_id, data in self.data.atoms.items()}
        bonds = {key: {'order': data['order'], 'stereo': data.get('stereo', 0)} for key, data in self.data.bonds.items()}
        state = {'atoms': atoms, 'bonds': bonds, '_next_atom_id': self.data._next_atom_id}

        state['version'] = VERSION 
        
        if self.current_mol: state['mol_3d'] = self.current_mol.ToBinary()

        state['is_3d_viewer_mode'] = not self.is_2d_editable

        json_safe_constraints = []
        try:
            for const in self.constraints_3d:
                # (Type, (Idx...), Value, Force) -> [Type, [Idx...], Value, Force]
                if len(const) == 4:
                    json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                else:
                    # 後方互換性: 3要素の場合はデフォルトForceを追加
                    json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
        except Exception:
            pass # 失敗したら空リスト
        state['constraints_3d'] = json_safe_constraints
            
        return state

    def set_state_from_data(self, state_data):
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)
        
        loaded_data = copy.deepcopy(state_data)

        # ファイルのバージョンを取得（存在しない場合は '0.0.0' とする）
        file_version_str = loaded_data.get('version', '0.0.0')

        try:
            app_version_parts = tuple(map(int, VERSION.split('.')))
            file_version_parts = tuple(map(int, file_version_str.split('.')))

            # ファイルのバージョンがアプリケーションのバージョンより新しい場合に警告
            if file_version_parts > app_version_parts:
                QMessageBox.warning(
                    self,
                    "Version Mismatch",
                    f"The file you are opening was saved with a newer version of MoleditPy (ver. {file_version_str}).\n\n"
                    f"Your current version is {VERSION}.\n\n"
                    "Some features may not load or work correctly."
                )
        except (ValueError, AttributeError):
            pass

        raw_atoms = loaded_data.get('atoms', {})
        raw_bonds = loaded_data.get('bonds', {})

        # 制約データの復元 (pmeraw)
        try:
            loaded_constraints = loaded_data.get("constraints_3d", [])
            # pmerawもJSON互換形式 [Type, [Idx...], Value, Force] で保存されている想定
            self.constraints_3d = []
            for const in loaded_constraints:
                if isinstance(const, list):
                    if len(const) == 4:
                        # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                        self.constraints_3d.append((const[0], tuple(const[1]), const[2], const[3]))
                    elif len(const) == 3:
                        # 後方互換性: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                        self.constraints_3d.append((const[0], tuple(const[1]), const[2], 1.0e5))
        except Exception:
            self.constraints_3d = [] # 読み込み失敗時はリセット

        for atom_id, data in raw_atoms.items():
            pos = QPointF(data['pos'][0], data['pos'][1])
            charge = data.get('charge', 0)
            radical = data.get('radical', 0)  # <-- ラジカル情報を取得
            # AtomItem生成時にradicalを渡す
            atom_item = AtomItem(atom_id, data['symbol'], pos, charge=charge, radical=radical)
            # self.data.atomsにもradical情報を格納する
            self.data.atoms[atom_id] = {'symbol': data['symbol'], 'pos': pos, 'item': atom_item, 'charge': charge, 'radical': radical}
            self.scene.addItem(atom_item)
        
        self.data._next_atom_id = loaded_data.get('_next_atom_id', max(self.data.atoms.keys()) + 1 if self.data.atoms else 0)

        for key_tuple, data in raw_bonds.items():
            id1, id2 = key_tuple
            if id1 in self.data.atoms and id2 in self.data.atoms:
                atom1_item = self.data.atoms[id1]['item']; atom2_item = self.data.atoms[id2]['item']
                bond_item = BondItem(atom1_item, atom2_item, data.get('order', 1), data.get('stereo', 0))
                self.data.bonds[key_tuple] = {'order': data.get('order', 1), 'stereo': data.get('stereo', 0), 'item': bond_item}
                atom1_item.bonds.append(bond_item); atom2_item.bonds.append(bond_item)
                self.scene.addItem(bond_item)

        for atom_data in self.data.atoms.values():
            if atom_data['item']: atom_data['item'].update_style()
        self.scene.update()

        if 'mol_3d' in loaded_data and loaded_data['mol_3d'] is not None:
            try:
                self.current_mol = Chem.Mol(loaded_data['mol_3d'])
                # デバッグ：3D構造が有効かチェック
                if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                    self.draw_molecule_3d(self.current_mol)
                    self.plotter.reset_camera()
                    # 3D関連機能を統一的に有効化
                    self._enable_3d_features(True)
                    
                    # 3D原子情報ホバー表示を再設定
                    self.setup_3d_hover()
                else:
                    # 無効な3D構造の場合
                    self.current_mol = None
                    self.plotter.clear()
                    # 3D関連機能を統一的に無効化
                    self._enable_3d_features(False)
            except Exception as e:
                self.statusBar().showMessage(f"Could not load 3D model from project: {e}")
                self.current_mol = None
                # 3D関連機能を統一的に無効化
                self._enable_3d_features(False)
        else:
            self.current_mol = None; self.plotter.clear(); self.analysis_action.setEnabled(False)
            self.optimize_3d_button.setEnabled(False)
            # 3D関連機能を統一的に無効化
            self._enable_3d_features(False)

        self.update_implicit_hydrogens()
        self.update_chiral_labels()

        if loaded_data.get('is_3d_viewer_mode', False):
            self._enter_3d_viewer_ui_mode()
            self.statusBar().showMessage("Project loaded in 3D Viewer Mode.")
        else:
            self.restore_ui_for_editing()
            # 3D分子がある場合は、2Dエディタモードでも3D編集機能を有効化
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                self._enable_3d_edit_actions(True)
        
        # undo/redo後に測定ラベルの位置を更新
        self.update_2d_measurement_labels()
        

    def push_undo_state(self):
        current_state_for_comparison = {
            'atoms': {k: (v['symbol'], v['item'].pos().x(), v['item'].pos().y(), v.get('charge', 0), v.get('radical', 0)) for k, v in self.data.atoms.items()},
            'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in self.data.bonds.items()},
            '_next_atom_id': self.data._next_atom_id,
            'mol_3d': self.current_mol.ToBinary() if self.current_mol else None
        }
        
        last_state_for_comparison = None
        if self.undo_stack:
            last_state = self.undo_stack[-1]
            last_atoms = last_state.get('atoms', {})
            last_bonds = last_state.get('bonds', {})
            last_state_for_comparison = {
                'atoms': {k: (v['symbol'], v['pos'][0], v['pos'][1], v.get('charge', 0), v.get('radical', 0)) for k, v in last_atoms.items()},
                'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in last_bonds.items()},
                '_next_atom_id': last_state.get('_next_atom_id'),
                'mol_3d': last_state.get('mol_3d', None)
            }

        if not last_state_for_comparison or current_state_for_comparison != last_state_for_comparison:
            state = self.get_current_state()
            self.undo_stack.append(state)
            self.redo_stack.clear()
            # 初期化完了後のみ変更があったことを記録
            if self.initialization_complete:
                self.has_unsaved_changes = True
                self.update_window_title()
        
        self.update_implicit_hydrogens()
        self.update_realtime_info()
        self.update_undo_redo_actions()

    def update_window_title(self):
        """ウィンドウタイトルを更新（保存状態を反映）"""
        base_title = f"MoleditPy Ver. {VERSION}"
        if self.current_file_path:
            filename = os.path.basename(self.current_file_path)
            title = f"{filename} - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        else:
            # Untitledファイルとして扱う
            title = f"Untitled - {base_title}"
            if self.has_unsaved_changes:
                title = f"*{title}"
        self.setWindowTitle(title)

    def check_unsaved_changes(self):
        """未保存の変更があるかチェックし、警告ダイアログを表示"""
        if not self.has_unsaved_changes:
            return True  # 保存済みまたは変更なし
        
        if not self.data.atoms and self.current_mol is None:
            return True  # 空のドキュメント
        
        reply = QMessageBox.question(
            self,
            "Unsaved Changes",
            "You have unsaved changes. Do you want to save them?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            # 拡張子がPMEPRJでなければ「名前を付けて保存」
            file_path = self.current_file_path
            if not file_path or not file_path.lower().endswith('.pmeprj'):
                self.save_project_as()
            else:
                self.save_project()
            return not self.has_unsaved_changes  # 保存に成功した場合のみTrueを返す
        elif reply == QMessageBox.StandardButton.No:
            return True  # 保存せずに続行
        else:
            return False  # キャンセル

    def reset_undo_stack(self):
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.push_undo_state()

    def undo(self):
        if len(self.undo_stack) > 1:
            self.redo_stack.append(self.undo_stack.pop())
            state = self.undo_stack[-1]
            self.set_state_from_data(state)
            
            # Undo後に3D構造の状態に基づいてメニューを再評価
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D構造がある場合は3D編集機能を有効化
                self._enable_3d_edit_actions(True)
            else:
                # 3D構造がない場合は3D編集機能を無効化
                self._enable_3d_edit_actions(False)
                    
        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus() 

    def redo(self):
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self.set_state_from_data(state)
            
            # Redo後に3D構造の状態に基づいてメニューを再評価
            if self.current_mol and self.current_mol.GetNumAtoms() > 0:
                # 3D構造がある場合は3D編集機能を有効化
                self._enable_3d_edit_actions(True)
            else:
                # 3D構造がない場合は3D編集機能を無効化
                self._enable_3d_edit_actions(False)
                    
        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus() 
        
    def update_undo_redo_actions(self):
        self.undo_action.setEnabled(len(self.undo_stack) > 1)
        self.redo_action.setEnabled(len(self.redo_stack) > 0)

    def update_realtime_info(self):
        """ステータスバーの右側に現在の分子情報を表示する"""
        if not self.data.atoms:
            self.formula_label.setText("")  # 原子がなければ右側のラベルをクリア
            return

        try:
            mol = self.data.to_rdkit_mol()
            if mol:
                # 水素原子を明示的に追加した分子オブジェクトを生成
                mol_with_hs = Chem.AddHs(mol)
                mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                # 水素を含む分子オブジェクトから原子数を取得
                num_atoms = mol_with_hs.GetNumAtoms()
                # 右側のラベルのテキストを更新
                self.formula_label.setText(f"Formula: {mol_formula}   |   Atoms: {num_atoms}")
        except Exception:
            # 計算に失敗してもアプリは継続
            self.formula_label.setText("Invalid structure")

    def select_all(self):
        for item in self.scene.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.setSelected(True)

    def show_about_dialog(self):
        """Show the custom About dialog with Easter egg functionality"""
        dialog = AboutDialog(self, self)
        dialog.exec()

    def clear_all(self):
        # 未保存の変更があるかチェック
        if not self.check_unsaved_changes():
            return  # ユーザーがキャンセルした場合は何もしない

        self.restore_ui_for_editing()

        # データが存在しない場合は何もしない
        if not self.data.atoms and self.current_mol is None:
            return
        
        # 3Dモードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化
        
        # 3D原子選択をクリア
        self.clear_3d_selection()
        
        self.dragged_atom_info = None
            
        # 2Dエディタをクリアする（Undoスタックにはプッシュしない）
        self.clear_2d_editor(push_to_undo=False)
        
        # 3Dモデルをクリアする
        self.current_mol = None
        self.plotter.clear()
        self.constraints_3d = []
        
        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        # Undo/Redoスタックをリセットする
        self.reset_undo_stack()
        
        # ファイル状態をリセット（新規ファイル状態に）
        self.has_unsaved_changes = False
        self.current_file_path = None
        self.update_window_title()
        
        # 2Dビューのズームをリセット
        self.reset_zoom()
        
        # シーンとビューの明示的な更新
        self.scene.update()
        if self.view_2d:
            self.view_2d.viewport().update()

        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        # 3Dプロッターの再描画
        self.plotter.render()
        
        # メニューテキストと状態を更新（分子がクリアされたので通常の表示に戻す）
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()
        
        # アプリケーションのイベントループを強制的に処理し、画面の再描画を確実に行う
        QApplication.processEvents()
        
        self.statusBar().showMessage("Cleared all data.")
        
    def clear_2d_editor(self, push_to_undo=True):
        self.data = MolecularData()
        self.scene.data = self.data
        self.scene.clear()
        self.scene.reinitialize_items()
        self.is_xyz_derived = False  # 2Dエディタをクリアする際にXYZ由来フラグもリセット
        
        # 測定ラベルもクリア
        self.clear_2d_measurement_labels()
        
        # Clear 3D data and disable 3D-related menus
        self.current_mol = None
        self.plotter.clear()
        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        if push_to_undo:
            self.push_undo_state()

    def update_implicit_hydrogens(self):
        """現在の2D構造に基づいて各原子の暗黙の水素数を計算し、AtomItemに反映する"""
        # Quick guards: nothing to do if no atoms or no QApplication
        if not self.data.atoms:
            return

        # If called from non-GUI thread, schedule the heavy RDKit work here but
        # always perform UI mutations on the main thread via QTimer.singleShot.
        try:
            # Bump a local token to identify this request. The closure we
            # schedule below will capture `my_token` and will only apply UI
            # changes if the token still matches the most recent global
            # counter. This avoids applying stale updates after deletions or
            # teardown.
            try:
                self._ih_update_counter += 1
            except Exception:
                self._ih_update_counter = getattr(self, '_ih_update_counter', 0) or 1
            my_token = self._ih_update_counter

            mol = None
            try:
                mol = self.data.to_rdkit_mol()
            except Exception:
                mol = None

            # Build a mapping of original_id -> hydrogen count without touching Qt items
            h_count_map = {}

            if mol is None:
                # Invalid/unsanitizable structure: reset all counts to 0
                for atom_id in list(self.data.atoms.keys()):
                    h_count_map[atom_id] = 0
            else:
                for atom in mol.GetAtoms():
                    try:
                        if not atom.HasProp("_original_atom_id"):
                            continue
                        original_id = atom.GetIntProp("_original_atom_id")

                        # Robust retrieval of H counts: prefer implicit, fallback to total or 0
                        try:
                            h_count = int(atom.GetNumImplicitHs())
                        except Exception:
                            try:
                                h_count = int(atom.GetTotalNumHs())
                            except Exception:
                                h_count = 0

                        h_count_map[int(original_id)] = h_count
                    except Exception:
                        # Skip problematic RDKit atoms
                        continue

            # Compute a per-atom problem map (original_id -> bool) so the
            # UI closure can safely set AtomItem.has_problem on the main thread.
            problem_map = {}
            try:
                if mol is not None:
                    try:
                        problems = Chem.DetectChemistryProblems(mol)
                    except Exception:
                        problems = None

                    if problems:
                        for prob in problems:
                            try:
                                atom_idx = prob.GetAtomIdx()
                                rd_atom = mol.GetAtomWithIdx(atom_idx)
                                if rd_atom and rd_atom.HasProp("_original_atom_id"):
                                    orig = int(rd_atom.GetIntProp("_original_atom_id"))
                                    problem_map[orig] = True
                            except Exception:
                                continue
                else:
                    # Fallback: use a lightweight valence heuristic similar to
                    # check_chemistry_problems_fallback() so we still flag atoms
                    # when RDKit conversion wasn't possible.
                    for atom_id, atom_data in self.data.atoms.items():
                        try:
                            symbol = atom_data.get('symbol')
                            charge = atom_data.get('charge', 0)
                            bond_count = 0
                            for (id1, id2), bond_data in self.data.bonds.items():
                                if id1 == atom_id or id2 == atom_id:
                                    bond_count += bond_data.get('order', 1)

                            is_problematic = False
                            if symbol == 'C' and bond_count > 4:
                                is_problematic = True
                            elif symbol == 'N' and bond_count > 3 and charge == 0:
                                is_problematic = True
                            elif symbol == 'O' and bond_count > 2 and charge == 0:
                                is_problematic = True
                            elif symbol == 'H' and bond_count > 1:
                                is_problematic = True
                            elif symbol in ['F', 'Cl', 'Br', 'I'] and bond_count > 1 and charge == 0:
                                is_problematic = True

                            if is_problematic:
                                problem_map[atom_id] = True
                        except Exception:
                            continue
            except Exception:
                # If any unexpected error occurs while building the map, fall back
                # to an empty map so we don't accidentally crash the UI.
                problem_map = {}

            # Schedule UI updates on the main thread to avoid calling Qt methods from
            # background threads or during teardown (which can crash the C++ layer).
            def _apply_ui_updates():
                # If the global counter changed since this closure was
                # created, bail out — the update is stale.
                try:
                    if my_token != getattr(self, '_ih_update_counter', None):
                        return
                except Exception:
                    # If anything goes wrong checking the token, be conservative
                    # and skip the update to avoid touching possibly-damaged
                    # Qt wrappers.
                    return

                # Work on a shallow copy/snapshot of the data.atoms mapping so
                # that concurrent mutations won't raise KeyError during
                # iteration. We still defensively check each item below.
                try:
                    atoms_snapshot = dict(self.data.atoms)
                except Exception:
                    atoms_snapshot = {}
                # Prefer the module-level SIP helper to avoid repeated imports
                # and centralize exception handling. _sip_isdeleted is set at
                # import time above; fall back to None if unavailable.
                is_deleted_func = _sip_isdeleted if _sip_isdeleted is not None else None

                items_to_update = []
                for atom_id, atom_data in atoms_snapshot.items():
                    try:
                        item = atom_data.get('item')
                        if not item:
                            continue

                        # If sip.isdeleted is available, skip deleted C++ wrappers
                        try:
                            if is_deleted_func and is_deleted_func(item):
                                continue
                        except Exception:
                            # If sip check itself fails, continue with other lightweight guards
                            pass

                        # If the item is no longer in a scene, skip updating it to avoid
                        # touching partially-deleted objects during scene teardown.
                        try:
                            sc = item.scene() if hasattr(item, 'scene') else None
                            if sc is None:
                                continue
                        except Exception:
                            # Accessing scene() might fail for a damaged object; skip it
                            continue

                        # Desired new count (default to 0 if not computed)
                        new_count = h_count_map.get(atom_id, 0)

                        current = getattr(item, 'implicit_h_count', None)
                        current_prob = getattr(item, 'has_problem', False)
                        desired_prob = problem_map.get(atom_id, False)

                        # If neither the implicit-H count nor the problem flag
                        # changed, skip this item.
                        if current == new_count and current_prob == desired_prob:
                            continue

                        # Only prepare a geometry change if the implicit H count
                        # changes (this may affect the item's bounding rect).
                        need_geometry = (current != new_count)
                        try:
                            if need_geometry and hasattr(item, 'prepareGeometryChange'):
                                try:
                                    item.prepareGeometryChange()
                                except Exception:
                                    pass

                            # Apply implicit hydrogen count (guarded)
                            try:
                                item.implicit_h_count = new_count
                            except Exception:
                                # If setting the count fails, continue but still
                                # attempt to set the problem flag below.
                                pass

                            # Apply problem flag (visual red-outline)
                            try:
                                item.has_problem = bool(desired_prob)
                            except Exception:
                                pass

                            # Ensure the item is updated in the scene so paint() runs
                            # when either geometry or problem-flag changed.
                            items_to_update.append(item)
                        except Exception:
                            # Non-fatal: skip problematic items
                            continue

                    except Exception:
                        continue

                # Trigger updates once for unique items; wrap in try/except to avoid crashes
                # Trigger updates once for unique items; dedupe by object id so
                # we don't attempt to hash QGraphicsItem wrappers which may
                # behave oddly when partially deleted.
                seen = set()
                for it in items_to_update:
                    try:
                        if it is None:
                            continue
                        oid = id(it)
                        if oid in seen:
                            continue
                        seen.add(oid)
                        if hasattr(it, 'update'):
                            try:
                                it.update()
                            except Exception:
                                # ignore update errors for robustness
                                pass
                    except Exception:
                        # Ignore any unexpected errors when touching the item
                        continue

            # Always schedule on main thread asynchronously
            try:
                QTimer.singleShot(0, _apply_ui_updates)
            except Exception:
                # Fallback: try to call directly (best-effort)
                try:
                    _apply_ui_updates()
                except Exception:
                    pass

        except Exception:
            # Make sure update failures never crash the application
            pass


    def import_smiles_dialog(self):
        """ユーザーにSMILES文字列の入力を促すダイアログを表示する"""
        smiles, ok = QInputDialog.getText(self, "Import SMILES", "Enter SMILES string:")
        if ok and smiles:
            self.load_from_smiles(smiles)

    def import_inchi_dialog(self):
        """ユーザーにInChI文字列の入力を促すダイアログを表示する"""
        inchi, ok = QInputDialog.getText(self, "Import InChI", "Enter InChI string:")
        if ok and inchi:
            self.load_from_inchi(inchi)

    def load_from_smiles(self, smiles_string):
        """SMILES文字列から分子を読み込み、2Dエディタに表示する"""
        try:
            if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない

            cleaned_smiles = smiles_string.strip()
            
            mol = Chem.MolFromSmiles(cleaned_smiles)
            if mol is None:
                if not cleaned_smiles:
                    raise ValueError("SMILES string was empty.")
                raise ValueError("Invalid SMILES string.")

            AllChem.Compute2DCoords(mol)
            Chem.Kekulize(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)

            conf = mol.GetConformer()
            SCALE_FACTOR = 50.0
            
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = sum(p.x for p in positions) / len(positions) if positions else 0.0
            mol_center_y = sum(p.y for p in positions) / len(positions) if positions else 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id
            

            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                # 単結合の立体
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1 # Wedge
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2 # Dash
                # 二重結合のE/Z
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3 # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4 # E

                if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                    a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                    a1_item = self.data.atoms[a1_id]['item']
                    a2_item = self.data.atoms[a2_id]['item']
                    self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            self.statusBar().showMessage(f"Successfully loaded from SMILES.")
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)
            
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid SMILES: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading from SMILES: {e}")
            
            traceback.print_exc()

    def load_from_inchi(self, inchi_string):
        """InChI文字列から分子を読み込み、2Dエディタに表示する"""
        try:
            if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない
            cleaned_inchi = inchi_string.strip()
            
            mol = Chem.MolFromInchi(cleaned_inchi)
            if mol is None:
                if not cleaned_inchi:
                    raise ValueError("InChI string was empty.")
                raise ValueError("Invalid InChI string.")

            AllChem.Compute2DCoords(mol)
            Chem.Kekulize(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)

            conf = mol.GetConformer()
            SCALE_FACTOR = 50.0
            
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = sum(p.x for p in positions) / len(positions) if positions else 0.0
            mol_center_y = sum(p.y for p in positions) / len(positions) if positions else 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id
            
            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                # 単結合の立体
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1 # Wedge
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2 # Dash
                # 二重結合のE/Z
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3 # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4 # E

                if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                    a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                    a1_item = self.data.atoms[a1_id]['item']
                    a2_item = self.data.atoms[a2_id]['item']
                    self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            self.statusBar().showMessage(f"Successfully loaded from InChI.")
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)
            
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid InChI: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading from InChI: {e}")
            
            traceback.print_exc()
    
    def fix_mol_counts_line(self, line: str) -> str:
        """
        Check and fix the CTAB counts line in a MOL file.
        If the line already contains 'V3000' or 'V2000' it is left unchanged.
        Otherwise the line is treated as V2000 and the proper 39-character
        format (33 chars of counts + ' V2000') is returned.
        """
        # If already V3000 or V2000, leave as-is
        if 'V3000' in line or 'V2000' in line:
            return line

        # Prepare prefix (first 33 characters for the 11 * I3 fields)
        prefix = line.rstrip().ljust(33)[0:33]
        version_str = ' V2000'
        return prefix + version_str

    def fix_mol_block(self, mol_block: str) -> str:
        """
        Given an entire MOL block as a string, ensure the 4th line (CTAB counts
        line) is valid. If the file has fewer than 4 lines, return as-is.
        """
        lines = mol_block.splitlines()
        if len(lines) < 4:
            # Not a valid MOL block — return unchanged
            return mol_block

        counts_line = lines[3]
        fixed_counts_line = self.fix_mol_counts_line(counts_line)
        lines[3] = fixed_counts_line
        return "\n".join(lines)

    def load_mol_file(self, file_path=None):
        if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Import MOL File", "", "Chemical Files (*.mol *.sdf);;All Files (*)")
            if not file_path: 
                return

        try:
            self.dragged_atom_info = None
            # If this is a single-record .mol file, read & fix the counts line
            # before parsing. For multi-record .sdf files, keep using SDMolSupplier.
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ''
            if ext == '.mol':
                # Read file text, fix CTAB counts line if needed, then parse
                with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
                if mol is None:
                    raise ValueError("Failed to read molecule from .mol file after fixing counts line.")
            else:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
                if mol is None:
                    raise ValueError("Failed to read molecule from file.")

            Chem.Kekulize(mol)

            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)
            
            # 1. 座標がなければ2D座標を生成する
            if mol.GetNumConformers() == 0: 
                AllChem.Compute2DCoords(mol)
            
            # 2. 座標の有無にかかわらず、常に立体化学を割り当て、2D表示用にくさび結合を設定する
            # これにより、3D座標を持つMOLファイルからでも正しく2Dの立体表現が生成される
            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            conf = mol.GetConformer()

            SCALE_FACTOR = 50.0
            
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())

            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if positions:
                mol_center_x = sum(p.x for p in positions) / len(positions)
                mol_center_y = sum(p.y for p in positions) / len(positions)
            else:
                mol_center_x, mol_center_y = 0.0, 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                
                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y
                
                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()
                
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id
                        
            for bond in mol.GetBonds():
                b_idx,e_idx=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble(); b_dir = bond.GetBondDir()
                stereo = 0
                # Check for single bond Wedge/Dash
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2
                # ADDED: Check for double bond E/Z stereochemistry
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3 # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4 # E

                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item,a2_item=self.data.atoms[a1_id]['item'],self.data.atoms[a2_id]['item']

                self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            self.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.reset_undo_stack()
            # NEWファイル扱い: ファイルパスをクリアし未保存状態はFalse（変更なければ保存警告なし）
            self.current_file_path = file_path
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid MOL file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading file: {e}")
            
            traceback.print_exc()
    
    def load_mol_for_3d_viewing(self):
        # moved to load_mol_file_for_3d_viewing
        pass
        '''
        file_path, _ = QFileDialog.getOpenFileName(self, "Load 3D MOL (View Only)", "", "Chemical Files (*.mol *.sdf);;All Files (*)")
        if not file_path:
            return

        try:
            # For single-record .mol files, read & fix counts line before parsing.
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ''
            if ext == '.mol':
                with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
                if mol is None:
                    raise ValueError("Failed to read .mol molecule after fixing counts line.")
            else:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
                if mol is None:
                    raise ValueError("Failed to read molecule.")

            if mol.GetNumConformers() == 0:
                raise ValueError("MOL file has no 3D coordinates.")

            # 2Dエディタをクリア
            self.clear_2d_editor(push_to_undo=False)
            
            # 3D構造をセットして描画
            self.current_mol = mol
            
            # 3Dファイル読み込み時はマッピングをクリア（2D構造がないため）
            self.atom_id_to_rdkit_idx_map = {}

            # Clear any leftover XYZ-derived flags on this molecule to ensure
            # Optimize 3D and related UI reflects the true source format.
            try:
                self._clear_xyz_flags(self.current_mol)
            except Exception:
                pass
            
            self.draw_molecule_3d(self.current_mol)
            self.plotter.reset_camera()

            # UIを3Dビューアモードに設定
            self._enter_3d_viewer_ui_mode()
            
            # 3D関連機能を統一的に有効化
            self._enable_3d_features(True)
            
            # メニューテキストと状態を更新
            self.update_atom_id_menu_text()
            self.update_atom_id_menu_state()
            
            self.statusBar().showMessage(f"3D Viewer Mode: Loaded {os.path.basename(file_path)}")
            self.reset_undo_stack()
            # NEWファイル扱い: ファイルパスをクリアし未保存状態はFalse（変更なければ保存警告なし）
            self.current_file_path = None
            self.has_unsaved_changes = False
            self.update_window_title()

        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
            self.restore_ui_for_editing()
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid 3D MOL file: {e}")
            self.restore_ui_for_editing()
        except Exception as e:
            self.statusBar().showMessage(f"Error loading 3D file: {e}")
            self.restore_ui_for_editing()
            
            traceback.print_exc()
            '''

    def load_xyz_for_3d_viewing(self, file_path=None):
        """XYZファイルを読み込んで3Dビューアで表示する"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Load 3D XYZ (View Only)", "", "XYZ Files (*.xyz);;All Files (*)")
            if not file_path:
                return

        try:
            mol = self.load_xyz_file(file_path)
            if mol is None:
                raise ValueError("Failed to create molecule from XYZ file.")
            if mol.GetNumConformers() == 0:
                raise ValueError("XYZ file has no 3D coordinates.")

            # 2Dエディタをクリア
            self.clear_2d_editor(push_to_undo=False)
            
            # 3D構造をセットして描画
            # Set the molecule. If bonds were determined (mol has bonds),
            # treat this the same as loading a MOL file: clear the XYZ-derived
            # flag and enable 3D optimization. Only mark as XYZ-derived and
            # disable 3D optimization when the molecule has no bond information.
            self.current_mol = mol

            # XYZファイル読み込み時はマッピングをクリア（2D構造がないため）
            self.atom_id_to_rdkit_idx_map = {}

            # If the loader marked the molecule as produced under skip_chemistry_checks,
            # always treat it as XYZ-derived and disable optimization. Otherwise
            # fall back to the existing behavior based on bond presence.
            skip_flag = False
            try:
                # Prefer RDKit int prop
                skip_flag = bool(self.current_mol.GetIntProp("_xyz_skip_checks"))
            except Exception:
                try:
                    skip_flag = bool(getattr(self.current_mol, '_xyz_skip_checks', False))
                except Exception:
                    skip_flag = False

            if skip_flag:
                self.is_xyz_derived = True
                if hasattr(self, 'optimize_3d_button'):
                    try:
                        self.optimize_3d_button.setEnabled(False)
                    except Exception:
                        pass
            else:
                try:
                    has_bonds = (self.current_mol.GetNumBonds() > 0)
                except Exception:
                    has_bonds = False

                if has_bonds:
                    self.is_xyz_derived = False
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            # Only enable optimize if the molecule is not considered XYZ-derived
                            if not getattr(self, 'is_xyz_derived', False):
                                self.optimize_3d_button.setEnabled(True)
                            else:
                                self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass
                else:
                    self.is_xyz_derived = True
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass

            self.draw_molecule_3d(self.current_mol)
            self.plotter.reset_camera()

            # UIを3Dビューアモードに設定
            self._enter_3d_viewer_ui_mode()

            # 3D関連機能を統一的に有効化
            self._enable_3d_features(True)
            
            # メニューテキストと状態を更新
            self.update_atom_id_menu_text()
            self.update_atom_id_menu_state()
            
            self.statusBar().showMessage(f"3D Viewer Mode: Loaded {os.path.basename(file_path)}")
            self.reset_undo_stack()
            # XYZファイル名をcurrent_file_pathにセットし、未保存状態はFalse
            self.current_file_path = file_path
            self.has_unsaved_changes = False
            self.update_window_title()

        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
            self.restore_ui_for_editing()
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid XYZ file: {e}")
            self.restore_ui_for_editing()
        except Exception as e:
            self.statusBar().showMessage(f"Error loading XYZ file: {e}")
            self.restore_ui_for_editing()
            
            traceback.print_exc()

    def load_xyz_file(self, file_path):
        """XYZファイルを読み込んでRDKitのMolオブジェクトを作成する"""
        
        if not self.check_unsaved_changes():
            return  # ユーザーがキャンセルした場合は何もしない

        try:
            # We will attempt one silent load with default charge=0 (no dialog).
            # If RDKit emits chemistry warnings (for example "Explicit valence ..."),
            # prompt the user once for an overall charge and retry. Only one retry is allowed.


            # Helper: prompt for charge once when needed
            # Returns a tuple: (charge_value_or_0, accepted:bool, skip_chemistry:bool)
            def prompt_for_charge():
                try:
                    # Create a custom dialog so we can provide a "Skip chemistry" button
                    dialog = QDialog(self)
                    dialog.setWindowTitle("Import XYZ Charge")
                    layout = QVBoxLayout(dialog)

                    label = QLabel("Enter total molecular charge:")
                    line_edit = QLineEdit(dialog)
                    line_edit.setText("")

                    # Standard OK/Cancel buttons
                    btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, parent=dialog)

                    # Additional Skip chemistry button
                    skip_btn = QPushButton("Skip chemistry", dialog)

                    # Horizontal layout for buttons
                    hl = QHBoxLayout()
                    hl.addWidget(btn_box)
                    hl.addWidget(skip_btn)

                    layout.addWidget(label)
                    layout.addWidget(line_edit)
                    layout.addLayout(hl)

                    result = {"accepted": False, "skip": False}

                    def on_ok():
                        result["accepted"] = True
                        dialog.accept()

                    def on_cancel():
                        dialog.reject()

                    def on_skip():
                        # Mark skip and accept so caller can proceed with skip behavior
                        result["skip"] = True
                        dialog.accept()

                    try:
                        btn_box.button(QDialogButtonBox.Ok).clicked.connect(on_ok)
                        btn_box.button(QDialogButtonBox.Cancel).clicked.connect(on_cancel)
                    except Exception:
                        # Fallback if button lookup fails
                        btn_box.accepted.connect(on_ok)
                        btn_box.rejected.connect(on_cancel)

                    skip_btn.clicked.connect(on_skip)

                    # Execute dialog modally
                    if dialog.exec_() != QDialog.Accepted:
                        return None, False, False

                    if result["skip"]:
                        # User chose to skip chemistry checks; return skip flag
                        return 0, True, True

                    if not result["accepted"]:
                        return None, False, False

                    charge_text = line_edit.text()
                except Exception:
                    # On any dialog creation error, fall back to simple input dialog
                    try:
                        charge_text, ok = QInputDialog.getText(self, "Import XYZ Charge", "Enter total molecular charge:", text="0")
                    except Exception:
                        return 0, True, False
                    if not ok:
                        return None, False, False
                    try:
                        return int(str(charge_text).strip()), True, False
                    except Exception:
                        try:
                            return int(float(str(charge_text).strip())), True, False
                        except Exception:
                            return 0, True, False

                if charge_text is None:
                    return None, False, False

                try:
                    return int(str(charge_text).strip()), True, False
                except Exception:
                    try:
                        return int(float(str(charge_text).strip())), True, False
                    except Exception:
                        return 0, True, False

            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 空行とコメント行を除去（但し、先頭2行は保持）
            non_empty_lines = []
            for i, line in enumerate(lines):
                stripped = line.strip()
                if i < 2:  # 最初の2行は原子数とコメント行なので保持
                    non_empty_lines.append(stripped)
                elif stripped and not stripped.startswith('#'):  # 空行とコメント行をスキップ
                    non_empty_lines.append(stripped)
            
            if len(non_empty_lines) < 2:
                raise ValueError("XYZ file format error: too few lines")
            
            # 原子数を読み取り
            try:
                num_atoms = int(non_empty_lines[0])
            except ValueError:
                raise ValueError("XYZ file format error: invalid atom count")
            
            if num_atoms <= 0:
                raise ValueError("XYZ file format error: atom count must be positive")
            
            # コメント行（2行目）
            comment = non_empty_lines[1] if len(non_empty_lines) > 1 else ""
            
            # 原子データを読み取り
            atoms_data = []
            data_lines = non_empty_lines[2:]
            
            if len(data_lines) < num_atoms:
                raise ValueError(f"XYZ file format error: expected {num_atoms} atom lines, found {len(data_lines)}")
            
            for i, line in enumerate(data_lines[:num_atoms]):
                parts = line.split()
                if len(parts) < 4:
                    raise ValueError(f"XYZ file format error: invalid atom data at line {i+3}")
                
                symbol = parts[0].strip()
                
                # 元素記号の妥当性をチェック
                try:
                    # RDKitで認識される元素かどうかをチェック
                    test_atom = Chem.Atom(symbol)
                except Exception:
                    # 認識されない場合、最初の文字を大文字にして再試行
                    symbol = symbol.capitalize()
                    try:
                        test_atom = Chem.Atom(symbol)
                    except Exception:
                        # If user requested to skip chemistry checks, coerce unknown symbols to C
                        if self.settings.get('skip_chemistry_checks', False):
                            symbol = 'C'
                        else:
                            raise ValueError(f"Unrecognized element symbol: {parts[0]} at line {i+3}")
                
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError:
                    raise ValueError(f"XYZ file format error: invalid coordinates at line {i+3}")
                
                atoms_data.append((symbol, x, y, z))
            
            if len(atoms_data) == 0:
                raise ValueError("XYZ file format error: no atoms found")
            
            # RDKitのMolオブジェクトを作成
            mol = Chem.RWMol()
            
            # 原子を追加
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                # XYZファイルでの原子のUniqueID（0ベースのインデックス）を保存
                atom.SetIntProp("xyz_unique_id", i)
                mol.AddAtom(atom)
            
            # 3D座標を設定
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)
            # If user requested to skip chemistry checks, bypass RDKit's
            # DetermineBonds/sanitization flow entirely and use only the
            # distance-based bond estimation. Treat the resulting molecule
            # as "XYZ-derived" (disable 3D optimization) and return it.
            try:
                skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
            except Exception:
                skip_checks = False

            if skip_checks:
                used_rd_determine = False
                try:
                    # Use the conservative distance-based heuristic to add bonds
                    self.estimate_bonds_from_distances(mol)
                except Exception:
                    # Non-fatal: continue even if distance-based estimation fails
                    pass

                # Finalize and return a plain Mol object
                try:
                    candidate_mol = mol.GetMol()
                except Exception:
                    try:
                        candidate_mol = Chem.Mol(mol)
                    except Exception:
                        candidate_mol = None

                if candidate_mol is None:
                    raise ValueError("Failed to create valid molecule object when skip_chemistry_checks=True")

                # Attach a default charge property
                try:
                    candidate_mol.SetIntProp("_xyz_charge", 0)
                except Exception:
                    try:
                        candidate_mol._xyz_charge = 0
                    except Exception:
                        pass

                # Mark that this molecule was produced via the skip-chemistry path
                try:
                    candidate_mol.SetIntProp("_xyz_skip_checks", 1)
                except Exception:
                    try:
                        candidate_mol._xyz_skip_checks = True
                    except Exception:
                        pass

                # Set UI flags consistently: mark as XYZ-derived and disable optimize
                try:
                    self.current_mol = candidate_mol
                    self.is_xyz_derived = True
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            self.optimize_3d_button.setEnabled(False)
                        except Exception:
                            pass
                except Exception:
                    pass

                # Store atom data for later analysis and return
                candidate_mol._xyz_atom_data = atoms_data
                return candidate_mol
            # We'll attempt silently first with charge=0 and only prompt the user
            # for a charge when the RDKit processing block fails (raises an
            # exception). If the user provides a charge, retry; allow repeated
            # prompts until the user cancels. This preserves the previous
            # fallback behaviors (skip_chemistry_checks, distance-based bond
            # estimation) and property attachments.
            used_rd_determine = False
            final_mol = None

            # First, try silently with charge=0. If that raises an exception we
            # will enter a loop prompting the user for a charge and retrying as
            # long as the user provides values. If the user cancels, return None.
            def _process_with_charge(charge_val):
                """Inner helper: attempt to build/finalize molecule with given charge.

                Returns the finalized RDKit Mol on success. May raise exceptions
                which will be propagated to the caller.
                """
                nonlocal used_rd_determine
                # Capture RDKit stderr while we run the processing to avoid
                # spamming the console. We won't treat warnings specially here;
                # only exceptions will trigger a prompt/retry. We also want to
                # distinguish failures originating from DetermineBonds so the
                # outer logic can decide whether to prompt the user repeatedly
                # for different charge values.
                buf = io.StringIO()
                determine_failed = False
                with contextlib.redirect_stderr(buf):
                    # Try DetermineBonds if available
                    try:
                        from rdkit.Chem import rdDetermineBonds
                        try:
                            try:
                                mol_candidate = Chem.RWMol(Chem.Mol(mol))
                            except Exception:
                                mol_candidate = Chem.RWMol(mol)

                            # This call may raise. If it does, mark determine_failed
                            # so the caller can prompt for a different charge.
                            rdDetermineBonds.DetermineBonds(mol_candidate, charge=charge_val)
                            mol_to_finalize = mol_candidate
                            used_rd_determine = True
                        except Exception:
                            # DetermineBonds failed for this charge value. We
                            # should allow the caller to prompt for another
                            # charge (or cancel). Mark the flag and re-raise a
                            # dedicated exception to be handled by the outer
                            # loop.
                            determine_failed = True
                            used_rd_determine = False
                            mol_to_finalize = mol
                            # Raise a sentinel exception to indicate DetermineBonds failure
                            raise RuntimeError("DetermineBondsFailed")
                    except RuntimeError:
                        # Propagate our sentinel so outer code can catch it.
                        raise
                    except Exception:
                        # rdDetermineBonds not available or import failed; use
                        # distance-based fallback below.
                        used_rd_determine = False
                        mol_to_finalize = mol

                    if not used_rd_determine:
                        # distance-based fallback
                        self.estimate_bonds_from_distances(mol_to_finalize)

                    # Finalize molecule
                    try:
                        candidate_mol = mol_to_finalize.GetMol()
                    except Exception:
                        candidate_mol = None

                    if candidate_mol is None:
                        # Try salvage path
                        try:
                            candidate_mol = mol.GetMol()
                        except Exception:
                            candidate_mol = None

                    if candidate_mol is None:
                        raise ValueError("Failed to create valid molecule object")

                    # Attach charge property if possible
                    try:
                        try:
                            candidate_mol.SetIntProp("_xyz_charge", int(charge_val))
                        except Exception:
                            try:
                                candidate_mol._xyz_charge = int(charge_val)
                            except Exception:
                                pass
                    except Exception:
                        pass

                    # Preserve whether the user requested skip_chemistry_checks
                    try:
                        if bool(self.settings.get('skip_chemistry_checks', False)):
                            try:
                                candidate_mol.SetIntProp("_xyz_skip_checks", 1)
                            except Exception:
                                try:
                                    candidate_mol._xyz_skip_checks = True
                                except Exception:
                                    pass
                    except Exception:
                        pass

                    # Run chemistry checks which may emit warnings to stderr
                    self._apply_chem_check_and_set_flags(candidate_mol, source_desc='XYZ')

                    # Accept the candidate
                    return candidate_mol

            # Silent first attempt
            try:
                final_mol = _process_with_charge(0)
            except RuntimeError as e_sentinel:
                # DetermineBonds explicitly failed for charge=0. In this
                # situation, repeatedly prompt the user for charges until
                # DetermineBonds succeeds or the user cancels.
                while True:
                    charge_val, ok, skip_flag = prompt_for_charge()
                    if not ok:
                        # user cancelled the prompt -> abort
                        return None
                    if skip_flag:
                        # User selected Skip chemistry: attempt distance-based salvage
                        try:
                            self.estimate_bonds_from_distances(mol)
                        except Exception:
                            pass
                        salvaged = None
                        try:
                            salvaged = mol.GetMol()
                        except Exception:
                            salvaged = None

                        if salvaged is not None:
                            try:
                                salvaged.SetIntProp("_xyz_skip_checks", 1)
                            except Exception:
                                try:
                                    salvaged._xyz_skip_checks = True
                                except Exception:
                                    pass
                            final_mol = salvaged
                            break
                        else:
                            # Could not salvage; abort
                            try:
                                self.statusBar().showMessage("Skip chemistry selected but failed to create salvaged molecule.")
                            except Exception:
                                pass
                            return None

                    try:
                        final_mol = _process_with_charge(charge_val)
                        # success -> break out of prompt loop
                        break
                    except RuntimeError:
                        # DetermineBonds still failing for this charge -> loop again
                        try:
                            self.statusBar().showMessage("DetermineBonds failed for that charge; please try a different total charge or cancel.")
                        except Exception:
                            pass
                        continue
                    except Exception as e_prompt:
                        # Some other failure occurred after DetermineBonds or in
                        # finalization. If skip_chemistry_checks is enabled we
                        # try the salvaged mol once; otherwise prompt again.
                        try:
                            skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
                        except Exception:
                            skip_checks = False

                        salvaged = None
                        try:
                            salvaged = mol.GetMol()
                        except Exception:
                            salvaged = None

                        if skip_checks and salvaged is not None:
                            final_mol = salvaged
                            # mark salvaged molecule as produced under skip_checks
                            try:
                                final_mol.SetIntProp("_xyz_skip_checks", 1)
                            except Exception:
                                try:
                                    final_mol._xyz_skip_checks = True
                                except Exception:
                                    pass
                            break
                        else:
                            try:
                                self.statusBar().showMessage(f"Retry failed: {e_prompt}")
                            except Exception:
                                pass
                            # Continue prompting
                            continue
            except Exception as e_silent:
                # If the silent attempt failed for reasons other than
                # DetermineBonds failing (e.g., finalization errors), fall
                # back to salvaging or prompting depending on settings.
                salvaged = None
                try:
                    salvaged = mol.GetMol()
                except Exception:
                    salvaged = None

                try:
                    skip_checks = bool(self.settings.get('skip_chemistry_checks', False))
                except Exception:
                    skip_checks = False

                if skip_checks and salvaged is not None:
                    final_mol = salvaged
                else:
                    # Repeatedly prompt until the user cancels or processing
                    # succeeds.
                    while True:
                        charge_val, ok, skip_flag = prompt_for_charge()
                        if not ok:
                            # user cancelled the prompt -> abort
                            return None
                        if skip_flag:
                            # User selected Skip chemistry: attempt distance-based salvage
                            try:
                                self.estimate_bonds_from_distances(mol)
                            except Exception:
                                pass
                            salvaged = None
                            try:
                                salvaged = mol.GetMol()
                            except Exception:
                                salvaged = None

                            if salvaged is not None:
                                try:
                                    salvaged.SetIntProp("_xyz_skip_checks", 1)
                                except Exception:
                                    try:
                                        salvaged._xyz_skip_checks = True
                                    except Exception:
                                        pass
                                final_mol = salvaged
                                break
                            else:
                                try:
                                    self.statusBar().showMessage("Skip chemistry selected but failed to create salvaged molecule.")
                                except Exception:
                                    pass
                                return None

                        try:
                            final_mol = _process_with_charge(charge_val)
                            # success -> break out of prompt loop
                            break
                        except RuntimeError:
                            # DetermineBonds failed for this charge -> let the
                            # user try another
                            try:
                                self.statusBar().showMessage("DetermineBonds failed for that charge; please try a different total charge or cancel.")
                            except Exception:
                                pass
                            continue
                        except Exception as e_prompt:
                            try:
                                self.statusBar().showMessage(f"Retry failed: {e_prompt}")
                            except Exception:
                                pass
                            continue

            # If we have a finalized molecule, apply the same UI flags and return
            if final_mol is not None:
                mol = final_mol
                try:
                    self.current_mol = mol

                    self.is_xyz_derived = not used_rd_determine
                    if hasattr(self, 'optimize_3d_button'):
                        try:
                            has_bonds = mol.GetNumBonds() > 0
                            # Respect the XYZ-derived flag: if the molecule is XYZ-derived,
                            # keep Optimize disabled regardless of bond detection.
                            if getattr(self, 'is_xyz_derived', False):
                                self.optimize_3d_button.setEnabled(False)
                            else:
                                self.optimize_3d_button.setEnabled(bool(has_bonds))
                        except Exception:
                            pass
                except Exception:
                    pass

                # Store original atom data for analysis
                mol._xyz_atom_data = atoms_data
                return mol
            
            # 元のXYZ原子データを分子オブジェクトに保存（分析用）
            mol._xyz_atom_data = atoms_data
            
            return mol
            
        except (OSError, IOError) as e:
            raise ValueError(f"File I/O error: {e}")
        except Exception as e:
            if "XYZ file format error" in str(e) or "Unrecognized element" in str(e):
                raise e
            else:
                raise ValueError(f"Error parsing XYZ file: {e}")

    def estimate_bonds_from_distances(self, mol):
        """原子間距離に基づいて結合を推定する"""
        
        # 一般的な共有結合半径（Ångström）- より正確な値
        covalent_radii = {
            'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76,
            'N': 0.75, 'O': 0.73, 'F': 0.71, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41,
            'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
            'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.14, 'Kr': 1.16,
            'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
            'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.33, 'Xe': 1.40
        }
        
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        # 追加された結合をトラッキング
        bonds_added = []
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)
                
                # 原子間距離を計算
                distance = rdMolTransforms.GetBondLength(conf, i, j)
                
                # 期待される結合距離を計算
                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()
                
                radius_i = covalent_radii.get(symbol_i, 1.0)  # デフォルト半径
                radius_j = covalent_radii.get(symbol_j, 1.0)
                
                expected_bond_length = radius_i + radius_j
                
                # 結合タイプによる許容範囲を調整
                # 水素結合は通常の共有結合より短い
                if symbol_i == 'H' or symbol_j == 'H':
                    tolerance_factor = 1.2  # 水素は結合が短くなりがち
                else:
                    tolerance_factor = 1.3  # 他の原子は少し余裕を持たせる
                
                max_bond_length = expected_bond_length * tolerance_factor
                min_bond_length = expected_bond_length * 0.5  # 最小距離も設定
                
                # 距離が期待値の範囲内なら結合を追加
                if min_bond_length <= distance <= max_bond_length:
                    try:
                        mol.AddBond(i, j, Chem.BondType.SINGLE)
                        bonds_added.append((i, j, distance))
                    except:
                        # 既に結合が存在する場合はスキップ
                        pass
        
        # デバッグ情報（オプション）
        # print(f"Added {len(bonds_added)} bonds based on distance analysis")
        
        return len(bonds_added)

    def save_project(self):
        """上書き保存（Ctrl+S）- デフォルトでPMEPRJ形式"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
        # 非ネイティブ形式（.mol, .sdf, .xyz など）は上書き保存せず、必ず「名前を付けて保存」にする
        native_exts = ['.pmeprj', '.pmeraw']
        if self.current_file_path and any(self.current_file_path.lower().endswith(ext) for ext in native_exts):
            # 既存のPMEPRJ/PMERAWファイルの場合は上書き保存
            try:
                if self.current_file_path.lower().endswith('.pmeraw'):
                    # 既存のPMERAWファイルの場合はPMERAW形式で保存
                    save_data = self.get_current_state()
                    with open(self.current_file_path, 'wb') as f: 
                        pickle.dump(save_data, f)
                else:
                    # PMEPRJ形式で保存
                    json_data = self.create_json_data()
                    with open(self.current_file_path, 'w', encoding='utf-8') as f: 
                        json.dump(json_data, f, indent=2, ensure_ascii=False)
                
                # 保存成功時に状態をリセット
                self.has_unsaved_changes = False
                self.update_window_title()
                
                self.statusBar().showMessage(f"Project saved to {self.current_file_path}")
                
            except (OSError, IOError) as e:
                self.statusBar().showMessage(f"File I/O error: {e}")
            except (pickle.PicklingError, TypeError, ValueError) as e:
                self.statusBar().showMessage(f"Data serialization error: {e}")
            except Exception as e: 
                self.statusBar().showMessage(f"Error saving project file: {e}")
                
                traceback.print_exc()
        else:
            # MOL/SDF/XYZなどは上書き保存せず、必ず「名前を付けて保存」にする
            self.save_project_as()

    def save_project_as(self):
        """名前を付けて保存（Ctrl+Shift+S）- デフォルトでPMEPRJ形式"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            # Determine a sensible default filename based on current file (strip extension)
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # Prefer the directory of the currently opened file as default
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save Project As", default_path, 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeprj'): 
                file_path += '.pmeprj'
            
            # JSONデータを保存
            json_data = self.create_json_data()
            with open(file_path, 'w', encoding='utf-8') as f: 
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Replace current file with the newly saved file so subsequent saves go to this path
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving project file: {e}")
            
            traceback.print_exc()

    def save_raw_data(self):
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            save_data = self.get_current_state()
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save Project File", default_path, "Project Files (*.pmeraw);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeraw'): 
                file_path += '.pmeraw'
                
            with open(file_path, 'wb') as f: 
                pickle.dump(save_data, f)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Update current file to the newly saved raw file
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving project file: {e}")
            
            traceback.print_exc()


    def load_raw_data(self, file_path=None):
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)")
            if not file_path: 
                return
        
        try:
            with open(file_path, 'rb') as f: 
                loaded_data = pickle.load(f)
            self.restore_ui_for_editing()
            self.set_state_from_data(loaded_data)
            
            # ファイル読み込み時に状態をリセット
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"Project loaded from {file_path}")
            
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.UnpicklingError as e:
            self.statusBar().showMessage(f"Invalid project file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading project file: {e}")
            
            traceback.print_exc()

    def save_as_json(self):
        """PMEJSONファイル形式で保存 (3D MOL情報含む)"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save as PME Project", default_path, 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeprj'): 
                file_path += '.pmeprj'
            
            # JSONデータを作成
            json_data = self.create_json_data()
            
            # JSON形式で保存（美しい整形付き）
            with open(file_path, 'w', encoding='utf-8') as f: 
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Replace current file with the newly saved PME Project
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"PME Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except (TypeError, ValueError) as e:
            self.statusBar().showMessage(f"JSON serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving PME Project file: {e}")
            
            traceback.print_exc()

    def create_json_data(self):
        """現在の状態をPMEJSON形式のデータに変換"""
        # 基本的なメタデータ
        json_data = {
            "format": "PME Project",
            "version": "1.0",
            "application": "MoleditPy",
            "application_version": VERSION,
            "created": str(QDateTime.currentDateTime().toString(Qt.DateFormat.ISODate)),
            "is_3d_viewer_mode": not self.is_2d_editable
        }
        
        # 2D構造データ
        if self.data.atoms:
            atoms_2d = []
            for atom_id, data in self.data.atoms.items():
                pos = data['item'].pos()
                atom_data = {
                    "id": atom_id,
                    "symbol": data['symbol'],
                    "x": pos.x(),
                    "y": pos.y(),
                    "charge": data.get('charge', 0),
                    "radical": data.get('radical', 0)
                }
                atoms_2d.append(atom_data)
            
            bonds_2d = []
            for (atom1_id, atom2_id), bond_data in self.data.bonds.items():
                bond_info = {
                    "atom1": atom1_id,
                    "atom2": atom2_id,
                    "order": bond_data['order'],
                    "stereo": bond_data.get('stereo', 0)
                }
                bonds_2d.append(bond_info)
            
            json_data["2d_structure"] = {
                "atoms": atoms_2d,
                "bonds": bonds_2d,
                "next_atom_id": self.data._next_atom_id
            }
        
        # 3D分子データ
        if self.current_mol and self.current_mol.GetNumConformers() > 0:
            try:
                # MOLデータをBase64エンコードで保存（バイナリデータの安全な保存）
                mol_binary = self.current_mol.ToBinary()
                mol_base64 = base64.b64encode(mol_binary).decode('ascii')
                
                # 3D座標を抽出
                atoms_3d = []
                if self.current_mol.GetNumConformers() > 0:
                    conf = self.current_mol.GetConformer()
                    for i in range(self.current_mol.GetNumAtoms()):
                        atom = self.current_mol.GetAtomWithIdx(i)
                        pos = conf.GetAtomPosition(i)

                        # Try to preserve original editor atom ID (if present) so it can be
                        # restored when loading PMEPRJ files. RDKit atom properties may
                        # contain _original_atom_id when the molecule was created from
                        # the editor's 2D structure.
                        original_id = None
                        try:
                            if atom.HasProp("_original_atom_id"):
                                original_id = atom.GetIntProp("_original_atom_id")
                        except Exception:
                            original_id = None

                        atom_3d = {
                            "index": i,
                            "symbol": atom.GetSymbol(),
                            "atomic_number": atom.GetAtomicNum(),
                            "x": pos.x,
                            "y": pos.y,
                            "z": pos.z,
                            "formal_charge": atom.GetFormalCharge(),
                            "num_explicit_hs": atom.GetNumExplicitHs(),
                            "num_implicit_hs": atom.GetNumImplicitHs(),
                            # include original editor atom id when available for round-trip
                            "original_id": original_id
                        }
                        atoms_3d.append(atom_3d)
                
                # 結合情報を抽出
                bonds_3d = []
                for bond in self.current_mol.GetBonds():
                    bond_3d = {
                        "atom1": bond.GetBeginAtomIdx(),
                        "atom2": bond.GetEndAtomIdx(),
                        "order": int(bond.GetBondType()),
                        "is_aromatic": bond.GetIsAromatic(),
                        "stereo": int(bond.GetStereo())
                    }
                    bonds_3d.append(bond_3d)
                
                # constraints_3dをJSON互換形式に変換
                json_safe_constraints = []
                try:
                    for const in self.constraints_3d:
                        if len(const) == 4:
                            json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                        else:
                            json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
                except Exception:
                    json_safe_constraints = []
                
                json_data["3d_structure"] = {
                    "mol_binary_base64": mol_base64,
                    "atoms": atoms_3d,
                    "bonds": bonds_3d,
                    "num_conformers": self.current_mol.GetNumConformers(),
                    "constraints_3d": json_safe_constraints
                }
                
                # 分子の基本情報
                json_data["molecular_info"] = {
                    "num_atoms": self.current_mol.GetNumAtoms(),
                    "num_bonds": self.current_mol.GetNumBonds(),
                    "molecular_weight": Descriptors.MolWt(self.current_mol),
                    "formula": rdMolDescriptors.CalcMolFormula(self.current_mol)
                }
                
                # SMILESとInChI（可能であれば）
                try:
                    json_data["identifiers"] = {
                        "smiles": Chem.MolToSmiles(self.current_mol),
                        "canonical_smiles": Chem.MolToSmiles(self.current_mol, canonical=True)
                    }
                    
                    # InChI生成を試行
                    try:
                        inchi = Chem.MolToInchi(self.current_mol)
                        inchi_key = Chem.MolToInchiKey(self.current_mol)
                        json_data["identifiers"]["inchi"] = inchi
                        json_data["identifiers"]["inchi_key"] = inchi_key
                    except:
                        pass  # InChI生成に失敗した場合は無視
                        
                except Exception as e:
                    print(f"Warning: Could not generate molecular identifiers: {e}")
                    
            except Exception as e:
                print(f"Warning: Could not process 3D molecular data: {e}")
        else:
            # 3D情報がない場合の記録
            json_data["3d_structure"] = None
            json_data["note"] = "No 3D structure available. Generate 3D coordinates first."

        # Record the last-successful optimization method (if any)
        # This is a convenience field so saved projects remember which
        # optimizer variant was last used (e.g. "MMFF94s", "MMFF94", "UFF").
        try:
            json_data["last_successful_optimization_method"] = getattr(self, 'last_successful_optimization_method', None)
        except Exception:
            json_data["last_successful_optimization_method"] = None
        
        return json_data

    def load_json_data(self, file_path=None):
        """PME Projectファイル形式を読み込み"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open PME Project File", "", 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path: 
                return
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f: 
                json_data = json.load(f)
            
            # フォーマット検証
            if json_data.get("format") != "PME Project":
                QMessageBox.warning(
                    self, "Invalid Format", 
                    "This file is not a valid PME Project format."
                )
                return
            
            # バージョン確認
            file_version = json_data.get("version", "1.0")
            if file_version != "1.0":
                QMessageBox.information(
                    self, "Version Notice", 
                    f"This file was created with PME Project version {file_version}.\n"
                    "Loading will be attempted but some features may not work correctly."
                )
            
            self.restore_ui_for_editing()
            self.load_from_json_data(json_data)
            # ファイル読み込み時に状態をリセット
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"PME Project loaded from {file_path}")
            
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except json.JSONDecodeError as e:
            self.statusBar().showMessage(f"Invalid JSON format: {e}")
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading PME Project file: {e}")
            
            traceback.print_exc()

    def open_project_file(self, file_path=None):
        """プロジェクトファイルを開く（.pmeprjと.pmerawの両方に対応）"""
        # Check for unsaved changes before opening a new project file.
        # Previously this function opened .pmeprj/.pmeraw without prompting the
        # user to save current unsaved work. Ensure we honor the global
        # unsaved-change check like other loaders (SMILES/MOL/etc.).
        if not self.check_unsaved_changes():
            return
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open Project File", "", 
                "PME Project Files (*.pmeprj);;PME Raw Files (*.pmeraw);;All Files (*)", 
            )
            if not file_path: 
                return
        
        # 拡張子に応じて適切な読み込み関数を呼び出し
        if file_path.lower().endswith('.pmeprj'):
            self.load_json_data(file_path)
        elif file_path.lower().endswith('.pmeraw'):
            self.load_raw_data(file_path)
        else:
            # 拡張子不明の場合はJSONとして試行
            try:
                self.load_json_data(file_path)
            except:
                try:
                    self.load_raw_data(file_path)
                except:
                    self.statusBar().showMessage("Error: Unable to determine file format.")

    def load_from_json_data(self, json_data):
        """JSONデータから状態を復元"""
        self.dragged_atom_info = None
        self.clear_2d_editor(push_to_undo=False)
        self._enable_3d_edit_actions(False)
        self._enable_3d_features(False)

        # 3Dビューアーモードの設定
        is_3d_mode = json_data.get("is_3d_viewer_mode", False)
        # Restore last successful optimization method if present in file
        try:
            self.last_successful_optimization_method = json_data.get("last_successful_optimization_method", None)
        except Exception:
            self.last_successful_optimization_method = None


        # 2D構造データの復元
        if "2d_structure" in json_data:
            structure_2d = json_data["2d_structure"]
            atoms_2d = structure_2d.get("atoms", [])
            bonds_2d = structure_2d.get("bonds", [])

            # 原子の復元
            for atom_data in atoms_2d:
                atom_id = atom_data["id"]
                symbol = atom_data["symbol"]
                pos = QPointF(atom_data["x"], atom_data["y"])
                charge = atom_data.get("charge", 0)
                radical = atom_data.get("radical", 0)

                atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
                self.data.atoms[atom_id] = {
                    'symbol': symbol,
                    'pos': pos,
                    'item': atom_item,
                    'charge': charge,
                    'radical': radical
                }
                self.scene.addItem(atom_item)

            # next_atom_idの復元
            self.data._next_atom_id = structure_2d.get(
                "next_atom_id",
                max([atom["id"] for atom in atoms_2d]) + 1 if atoms_2d else 0
            )

            # 結合の復元
            for bond_data in bonds_2d:
                atom1_id = bond_data["atom1"]
                atom2_id = bond_data["atom2"]

                if atom1_id in self.data.atoms and atom2_id in self.data.atoms:
                    atom1_item = self.data.atoms[atom1_id]['item']
                    atom2_item = self.data.atoms[atom2_id]['item']

                    bond_order = bond_data["order"]
                    stereo = bond_data.get("stereo", 0)

                    bond_item = BondItem(atom1_item, atom2_item, bond_order, stereo=stereo)
                    # 原子の結合リストに追加（重要：炭素原子の可視性判定で使用）
                    atom1_item.bonds.append(bond_item)
                    atom2_item.bonds.append(bond_item)

                    self.data.bonds[(atom1_id, atom2_id)] = {
                        'order': bond_order,
                        'item': bond_item,
                        'stereo': stereo
                    }
                    self.scene.addItem(bond_item)

            # --- ここで全AtomItemのスタイルを更新（炭素原子の可視性を正しく反映） ---
            for atom in self.data.atoms.values():
                atom['item'].update_style()
        # 3D構造データの復元
        if "3d_structure" in json_data:
            structure_3d = json_data["3d_structure"]

            # 制約データの復元 (JSONはタプルをリストとして保存するので、タプルに再変換)
            try:
                loaded_constraints = structure_3d.get("constraints_3d", [])
                self.constraints_3d = []
                for const in loaded_constraints:
                    if isinstance(const, list):
                        if len(const) == 4:
                            # [Type, [Idx...], Value, Force] -> (Type, (Idx...), Value, Force)
                            self.constraints_3d.append((const[0], tuple(const[1]), const[2], const[3]))
                        elif len(const) == 3:
                            # 後方互換性: [Type, [Idx...], Value] -> (Type, (Idx...), Value, 1.0e5)
                            self.constraints_3d.append((const[0], tuple(const[1]), const[2], 1.0e5))
            except Exception:
                self.constraints_3d = [] # 読み込み失敗時はリセット

            try:
                # バイナリデータの復元
                mol_base64 = structure_3d.get("mol_binary_base64")
                if mol_base64:
                    mol_binary = base64.b64decode(mol_base64.encode('ascii'))
                    self.current_mol = Chem.Mol(mol_binary)
                    if self.current_mol:
                        # 3D座標の設定
                        if self.current_mol.GetNumConformers() > 0:
                            conf = self.current_mol.GetConformer()
                            atoms_3d = structure_3d.get("atoms", [])
                            self.atom_positions_3d = np.zeros((len(atoms_3d), 3))
                            for atom_data in atoms_3d:
                                idx = atom_data["index"]
                                if idx < len(self.atom_positions_3d):
                                    self.atom_positions_3d[idx] = [
                                        atom_data["x"], 
                                        atom_data["y"], 
                                        atom_data["z"]
                                    ]
                                # Restore original editor atom id into RDKit atom property
                                try:
                                    original_id = atom_data.get("original_id", None)
                                    if original_id is not None and idx < self.current_mol.GetNumAtoms():
                                        rd_atom = self.current_mol.GetAtomWithIdx(idx)
                                        # set as int prop so other code expecting _original_atom_id works
                                        rd_atom.SetIntProp("_original_atom_id", int(original_id))
                                except Exception:
                                    pass
                            # Build mapping from original 2D atom IDs to RDKit indices so
                            # 3D picks can be synchronized back to 2D AtomItems.
                            try:
                                self.create_atom_id_mapping()
                                # update menu and UI states that depend on original IDs
                                try:
                                    self.update_atom_id_menu_text()
                                    self.update_atom_id_menu_state()
                                except Exception:
                                    pass
                            except Exception:
                                # non-fatal if mapping creation fails
                                pass

                        # 3D分子があれば必ず3D表示
                        self.draw_molecule_3d(self.current_mol)
                        # ViewerモードならUIも切り替え
                        if is_3d_mode:
                            self._enter_3d_viewer_ui_mode()
                        else:
                            self.is_2d_editable = True
                        self.plotter.reset_camera()

                        # 成功的に3D分子が復元されたので、3D関連UIを有効にする
                        try:
                            self._enable_3d_edit_actions(True)
                            self._enable_3d_features(True)
                        except Exception:
                            pass
                            
            except Exception as e:
                print(f"Warning: Could not restore 3D molecular data: {e}")
                self.current_mol = None

    def save_as_mol(self):
        try:
            mol_block = self.data.to_mol_block()
            if not mol_block: 
                self.statusBar().showMessage("Error: No 2D data to save."); 
                return
                
            lines = mol_block.split('\n')
            if len(lines) > 1 and 'RDKit' in lines[1]:
                lines[1] = '  MoleditPy Ver. ' + VERSION + '  2D'
            modified_mol_block = '\n'.join(lines)
            
            # default filename: based on current_file_path, append -2d for 2D mol
            default_name = "untitled-2d"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    name = os.path.splitext(base)[0]
                    default_name = f"{name}-2d"
            except Exception:
                default_name = "untitled-2d"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save 2D MOL File", default_path, "MOL Files (*.mol);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.mol'): 
                file_path += '.mol'
                
            with open(file_path, 'w', encoding='utf-8') as f: 
                f.write(modified_mol_block)
            self.statusBar().showMessage(f"2D data saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except UnicodeEncodeError as e:
            self.statusBar().showMessage(f"Text encoding error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving file: {e}")
            
            traceback.print_exc()
            
    def save_3d_as_mol(self):
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        try:
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    name = os.path.splitext(base)[0]
                    default_name = f"{name}"
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save 3D MOL File", default_path, "MOL Files (*.mol);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.mol'):
                file_path += '.mol'

            mol_to_save = Chem.Mol(self.current_mol)

            if mol_to_save.HasProp("_2D"):
                mol_to_save.ClearProp("_2D")

            mol_block = Chem.MolToMolBlock(mol_to_save, includeStereo=True)
            lines = mol_block.split('\n')
            if len(lines) > 1 and 'RDKit' in lines[1]:
                lines[1] = '  MoleditPy Ver. ' + VERSION + '  3D'
            modified_mol_block = '\n'.join(lines)
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(modified_mol_block)
            self.statusBar().showMessage(f"3D data saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except UnicodeEncodeError as e:
            self.statusBar().showMessage(f"Text encoding error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving 3D MOL file: {e}")
            
            traceback.print_exc()

    def save_as_xyz(self):
        if not self.current_mol: self.statusBar().showMessage("Error: Please generate a 3D structure first."); return
        # default filename based on current file
        default_name = "untitled"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}"
        except Exception:
            default_name = "untitled"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
        except Exception:
            default_path = default_name

        file_path,_=QFileDialog.getSaveFileName(self,"Save 3D XYZ File",default_path,"XYZ Files (*.xyz);;All Files (*)")
        if file_path:
            if not file_path.lower().endswith('.xyz'): file_path += '.xyz'
            try:
                conf=self.current_mol.GetConformer(); num_atoms=self.current_mol.GetNumAtoms()
                xyz_lines=[str(num_atoms)]
                # 電荷と多重度を計算
                try:
                    charge = Chem.GetFormalCharge(self.current_mol)
                except Exception:
                    charge = 0 # 取得失敗時は0
                
                try:
                    # 全原子のラジカル電子の合計を取得
                    num_radicals = Descriptors.NumRadicalElectrons(self.current_mol)
                    # スピン多重度を計算 (M = N + 1, N=ラジカル電子数)
                    multiplicity = num_radicals + 1
                except Exception:
                    multiplicity = 1 # 取得失敗時は 1 (singlet)

                smiles=Chem.MolToSmiles(Chem.RemoveHs(self.current_mol))
                xyz_lines.append(f"chrg = {charge}  mult = {multiplicity} | Generated by MoleditPy Ver. {VERSION}")
                for i in range(num_atoms):
                    pos=conf.GetAtomPosition(i); symbol=self.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
                with open(file_path,'w') as f: f.write("\n".join(xyz_lines) + "\n")
                self.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")

    def export_stl(self):
        """STLファイルとしてエクスポート（色なし）"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから直接データを取得（色情報なし）
            combined_mesh = self.export_from_3d_view_no_color()
            
            if combined_mesh is None or combined_mesh.n_points == 0:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            if not file_path.lower().endswith('.stl'):
                file_path += '.stl'
            
            combined_mesh.save(file_path, binary=True)
            self.statusBar().showMessage(f"STL exported to {file_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting STL: {e}")

    def export_obj_mtl(self):
        """OBJ/MTLファイルとしてエクスポート（表示中のモデルベース、色付き）"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as OBJ/MTL (with colors)", default_dir, "OBJ Files (*.obj);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから表示中のメッシュデータを色情報とともに取得
            meshes_with_colors = self.export_from_3d_view_with_colors()
            
            if not meshes_with_colors:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            # ファイル拡張子を確認・追加
            if not file_path.lower().endswith('.obj'):
                file_path += '.obj'
            
            # OBJ+MTL形式で保存（オブジェクトごとに色分け）
            mtl_path = file_path.replace('.obj', '.mtl')
            
            self.create_multi_material_obj(meshes_with_colors, file_path, mtl_path)
            
            self.statusBar().showMessage(f"OBJ+MTL files with individual colors exported to {file_path} and {mtl_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting OBJ/MTL: {e}")

            return meshes_with_colors
            
        except Exception as e:
            return []

    def create_multi_material_obj(self, meshes_with_colors, obj_path, mtl_path):
        """複数のマテリアルを持つOBJファイルとMTLファイルを作成（改良版）"""
        try:
            
            # MTLファイルを作成
            with open(mtl_path, 'w') as mtl_file:
                mtl_file.write(f"# Material file for {os.path.basename(obj_path)}\n")
                mtl_file.write(f"# Generated with individual object colors\n\n")
                
                for i, mesh_data in enumerate(meshes_with_colors):
                    color = mesh_data['color']
                    material_name = f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    
                    mtl_file.write(f"newmtl {material_name}\n")
                    mtl_file.write("Ka 0.2 0.2 0.2\n")  # Ambient
                    mtl_file.write(f"Kd {color[0]/255.0:.3f} {color[1]/255.0:.3f} {color[2]/255.0:.3f}\n")  # Diffuse
                    mtl_file.write("Ks 0.5 0.5 0.5\n")  # Specular
                    mtl_file.write("Ns 32.0\n")          # Specular exponent
                    mtl_file.write("illum 2\n")          # Illumination model
                    mtl_file.write("\n")
            
            # OBJファイルを作成
            with open(obj_path, 'w') as obj_file:
                obj_file.write(f"# OBJ file with multiple materials\n")
                obj_file.write(f"# Generated with individual object colors\n")
                obj_file.write(f"mtllib {os.path.basename(mtl_path)}\n\n")
                
                vertex_offset = 1  # OBJファイルの頂点インデックスは1から始まる
                
                for i, mesh_data in enumerate(meshes_with_colors):
                    mesh = mesh_data['mesh']
                    material_name = f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    
                    obj_file.write(f"# Object {i}: {mesh_data['name']}\n")
                    obj_file.write(f"# Color: RGB({mesh_data['color'][0]}, {mesh_data['color'][1]}, {mesh_data['color'][2]})\n")
                    obj_file.write(f"o object_{i}\n")
                    obj_file.write(f"usemtl {material_name}\n")
                    
                    # 頂点を書き込み
                    points = mesh.points
                    for point in points:
                        obj_file.write(f"v {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
                    
                    # 面を書き込み
                    for j in range(mesh.n_cells):
                        cell = mesh.get_cell(j)
                        if cell.type == 5:  # VTK_TRIANGLE
                            points_in_cell = cell.point_ids
                            v1 = points_in_cell[0] + vertex_offset
                            v2 = points_in_cell[1] + vertex_offset
                            v3 = points_in_cell[2] + vertex_offset
                            obj_file.write(f"f {v1} {v2} {v3}\n")
                        elif cell.type == 9:  # VTK_QUAD
                            points_in_cell = cell.point_ids
                            v1 = points_in_cell[0] + vertex_offset
                            v2 = points_in_cell[1] + vertex_offset
                            v3 = points_in_cell[2] + vertex_offset
                            v4 = points_in_cell[3] + vertex_offset
                            obj_file.write(f"f {v1} {v2} {v3} {v4}\n")
                    
                    vertex_offset += mesh.n_points
                    obj_file.write("\n")
                    
        except Exception as e:
            raise Exception(f"Failed to create multi-material OBJ: {e}")

    def export_color_stl(self):
        """カラーSTLファイルとしてエクスポート"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as Color STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから直接データを取得
            combined_mesh = self.export_from_3d_view()
            
            if combined_mesh is None or combined_mesh.n_points == 0:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            # STL形式で保存
            if not file_path.lower().endswith('.stl'):
                file_path += '.stl'
            combined_mesh.save(file_path, binary=True)
            self.statusBar().showMessage(f"STL exported to {file_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting STL: {e}")
    
    def export_from_3d_view(self):
        """現在の3Dビューから直接メッシュデータを取得"""
        try:
            
            # PyVistaプロッターから全てのアクターを取得
            combined_mesh = pv.PolyData()
            
            # プロッターのレンダラーからアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得する複数の方法を試行
                    mesh = None
                    
                    # 方法1: mapperのinputから取得
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        if hasattr(actor.mapper, 'input') and actor.mapper.input is not None:
                            mesh = actor.mapper.input
                        elif hasattr(actor.mapper, 'GetInput') and actor.mapper.GetInput() is not None:
                            mesh = actor.mapper.GetInput()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    # 方法3: PyVistaのメッシュデータベースから検索
                    if mesh is None:
                        for mesh_name, mesh_data in self.plotter.mesh.items():
                            if mesh_data is not None and mesh_data.n_points > 0:
                                mesh = mesh_data
                                break
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # 元のメッシュを変更しないようにコピーを作成
                        mesh_copy = mesh.copy()
                        
                        # コピーしたメッシュにカラー情報を追加
                        if hasattr(actor, 'prop') and hasattr(actor.prop, 'color'):
                            color = actor.prop.color
                            # RGB値を0-255の範囲に変換
                            rgb = np.array([int(c * 255) for c in color], dtype=np.uint8)
                            
                            # Blender対応のPLY形式用カラー属性を設定
                            mesh_copy.point_data['diffuse_red'] = np.full(mesh_copy.n_points, rgb[0], dtype=np.uint8)
                            mesh_copy.point_data['diffuse_green'] = np.full(mesh_copy.n_points, rgb[1], dtype=np.uint8) 
                            mesh_copy.point_data['diffuse_blue'] = np.full(mesh_copy.n_points, rgb[2], dtype=np.uint8)
                            
                            # 標準的なPLY形式もサポート
                            mesh_copy.point_data['red'] = np.full(mesh_copy.n_points, rgb[0], dtype=np.uint8)
                            mesh_copy.point_data['green'] = np.full(mesh_copy.n_points, rgb[1], dtype=np.uint8) 
                            mesh_copy.point_data['blue'] = np.full(mesh_copy.n_points, rgb[2], dtype=np.uint8)
                            
                            # 従来の colors 配列も保持（STL用）
                            mesh_colors = np.tile(rgb, (mesh_copy.n_points, 1))
                            mesh_copy.point_data['colors'] = mesh_colors
                        
                        # メッシュを結合
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)
                            
                except Exception:
                    continue
            
            return combined_mesh
            
        except Exception:
            return None

    def export_from_3d_view_no_color(self):
        """現在の3Dビューから直接メッシュデータを取得（色情報なし）"""
        try:
            
            # PyVistaプロッターから全てのアクターを取得
            combined_mesh = pv.PolyData()
            
            # プロッターのレンダラーからアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得する複数の方法を試行
                    mesh = None
                    
                    # 方法1: mapperのinputから取得
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        if hasattr(actor.mapper, 'input') and actor.mapper.input is not None:
                            mesh = actor.mapper.input
                        elif hasattr(actor.mapper, 'GetInput') and actor.mapper.GetInput() is not None:
                            mesh = actor.mapper.GetInput()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    # 方法3: PyVistaのメッシュデータベースから検索
                    if mesh is None:
                        for mesh_name, mesh_data in self.plotter.mesh.items():
                            if mesh_data is not None and mesh_data.n_points > 0:
                                mesh = mesh_data
                                break
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # 元のメッシュを変更しないようにコピーを作成（色情報は追加しない）
                        mesh_copy = mesh.copy()
                        
                        # メッシュを結合
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)
                            
                except Exception:
                    continue
            
            return combined_mesh
            
        except Exception:
            return None

    def export_from_3d_view_with_colors(self):
        """現在の3Dビューから直接メッシュデータを色情報とともに取得"""
        try:
            
            meshes_with_colors = []
            
            # PyVistaプロッターから全てのアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            actor_count = 0
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得
                    mesh = None
                    
                    # 方法1: mapperのinputから取得
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        if hasattr(actor.mapper, 'input') and actor.mapper.input is not None:
                            mesh = actor.mapper.input
                        elif hasattr(actor.mapper, 'GetInput') and actor.mapper.GetInput() is not None:
                            mesh = actor.mapper.GetInput()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # アクターから色情報を取得
                        color = [128, 128, 128]  # デフォルト色（グレー）
                        
                        try:
                            # VTKアクターのプロパティから色を取得
                            if hasattr(actor, 'prop') and actor.prop is not None:
                                vtk_color = actor.prop.GetColor()
                                color = [int(c * 255) for c in vtk_color]
                            elif hasattr(actor, 'GetProperty'):
                                prop = actor.GetProperty()
                                if prop is not None:
                                    vtk_color = prop.GetColor()
                                    color = [int(c * 255) for c in vtk_color]
                        except:
                            # 色取得に失敗した場合はデフォルト色をそのまま使用
                            pass
                        
                        # メッシュのコピーを作成
                        mesh_copy = mesh.copy()

                        # もしメッシュに頂点ごとの色情報が含まれている場合、
                        # それぞれの色ごとにサブメッシュに分割して個別マテリアルを作る。
                        # これにより、glyphs（すべての原子が一つのメッシュにまとめられる場合）でも
                        # 各原子の色を保持してOBJ/MTLへ出力できる。
                        try:
                            colors = None
                            pd = mesh_copy.point_data
                            # 優先的にred/green/blue配列を使用
                            if 'red' in pd and 'green' in pd and 'blue' in pd:
                                r = np.asarray(pd['red']).reshape(-1)
                                g = np.asarray(pd['green']).reshape(-1)
                                b = np.asarray(pd['blue']).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # diffuse_* のキーもサポート
                            elif 'diffuse_red' in pd and 'diffuse_green' in pd and 'diffuse_blue' in pd:
                                r = np.asarray(pd['diffuse_red']).reshape(-1)
                                g = np.asarray(pd['diffuse_green']).reshape(-1)
                                b = np.asarray(pd['diffuse_blue']).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # 単一の colors 配列があればそれを使う
                            elif 'colors' in pd:
                                colors = np.asarray(pd['colors'])

                            if colors is not None and colors.size > 0:
                                # 整数に変換。colors が 0-1 の float の場合は 255 倍して正規化する。
                                colors_arr = np.asarray(colors)
                                # 期待形状に整形
                                if colors_arr.ndim == 1:
                                    # 1次元の場合は単一チャンネルとして扱う
                                    colors_arr = colors_arr.reshape(-1, 1)

                                # float かどうか判定して正規化
                                if np.issubdtype(colors_arr.dtype, np.floating):
                                    # 値の最大が1付近なら0-1レンジとみなして255倍
                                    if colors_arr.max() <= 1.01:
                                        colors_int = np.clip((colors_arr * 255.0).round(), 0, 255).astype(np.int32)
                                    else:
                                        # 既に0-255レンジのfloatならそのまま丸める
                                        colors_int = np.clip(colors_arr.round(), 0, 255).astype(np.int32)
                                else:
                                    colors_int = np.clip(colors_arr, 0, 255).astype(np.int32)
                                # Ensure shape is (n_points, 3)
                                if colors_int.ndim == 1:
                                    # 単一値が入っている場合は同一RGBとして扱う
                                    colors_int = np.vstack([colors_int, colors_int, colors_int]).T

                                # 一意な色ごとにサブメッシュを抽出して追加
                                unique_colors, inverse = np.unique(colors_int, axis=0, return_inverse=True)
                                if unique_colors.shape[0] > 1:
                                    for uc_idx, uc in enumerate(unique_colors):
                                        point_inds = np.where(inverse == uc_idx)[0]
                                        if point_inds.size == 0:
                                            continue
                                        try:
                                            submesh = mesh_copy.extract_points(point_inds, adjacent_cells=True)
                                        except Exception:
                                            # extract_points が利用できない場合はスキップ
                                            continue
                                        if submesh is None or getattr(submesh, 'n_points', 0) == 0:
                                            continue
                                        color_rgb = [int(uc[0]), int(uc[1]), int(uc[2])]
                                        meshes_with_colors.append({
                                            'mesh': submesh,
                                            'color': color_rgb,
                                            'name': f'{actor_name}_color_{uc_idx}',
                                            'type': 'display_actor',
                                            'actor_name': actor_name
                                        })
                                    actor_count += 1
                                    # 分割したので以下の通常追加は行わない
                                    continue
                        except Exception:
                            # 分割処理に失敗した場合はフォールバックで単体メッシュを追加
                            pass

                        meshes_with_colors.append({
                            'mesh': mesh_copy,
                            'color': color,
                            'name': f'actor_{actor_count}_{actor_name}',
                            'type': 'display_actor',
                            'actor_name': actor_name
                        })
                        
                        actor_count += 1
                            
                except Exception as e:
                    print(f"Error processing actor {actor_name}: {e}")
                    continue
            
            return meshes_with_colors
            
        except Exception as e:
            print(f"Error in export_from_3d_view_with_colors: {e}")
            return []

    def export_2d_png(self):
        if not self.data.atoms:
            self.statusBar().showMessage("Nothing to export.")
            return

        # default filename: based on current file, append -2d for 2D exports
        default_name = "untitled-2d"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}-2d"
        except Exception:
            default_name = "untitled-2d"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
        except Exception:
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(self, "Export 2D as PNG", default_path, "PNG Files (*.png)")
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for a white background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.Yes)

        if reply == QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = (reply == QMessageBox.StandardButton.Yes)

        QApplication.processEvents()

        items_to_restore = {}
        original_background = self.scene.backgroundBrush()

        try:
            all_items = list(self.scene.items())
            for item in all_items:
                is_mol_part = isinstance(item, (AtomItem, BondItem))
                if not (is_mol_part and item.isVisible()):
                    items_to_restore[item] = item.isVisible()
                    item.hide()

            molecule_bounds = QRectF()
            for item in self.scene.items():
                if isinstance(item, (AtomItem, BondItem)) and item.isVisible():
                    molecule_bounds = molecule_bounds.united(item.sceneBoundingRect())

            if molecule_bounds.isEmpty() or not molecule_bounds.isValid():
                self.statusBar().showMessage("Error: Could not determine molecule bounds for export.")
                return

            if is_transparent:
                self.scene.setBackgroundBrush(QBrush(Qt.BrushStyle.NoBrush))
            else:
                self.scene.setBackgroundBrush(QBrush(QColor("#FFFFFF")))

            rect_to_render = molecule_bounds.adjusted(-20, -20, 20, 20)

            w = max(1, int(math.ceil(rect_to_render.width())))
            h = max(1, int(math.ceil(rect_to_render.height())))

            if w <= 0 or h <= 0:
                self.statusBar().showMessage("Error: Invalid image size calculated.")
                return

            image = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
            if is_transparent:
                image.fill(Qt.GlobalColor.transparent)
            else:
                image.fill(Qt.GlobalColor.white)

            painter = QPainter()
            ok = painter.begin(image)
            if not ok or not painter.isActive():
                self.statusBar().showMessage("Failed to start QPainter for image rendering.")
                return

            try:
                painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                target_rect = QRectF(0, 0, w, h)
                source_rect = rect_to_render
                self.scene.render(painter, target_rect, source_rect)
            finally:
                painter.end()

            saved = image.save(filePath, "PNG")
            if saved:
                self.statusBar().showMessage(f"2D view exported to {filePath}")
            else:
                self.statusBar().showMessage(f"Failed to save image. Check file path or permissions.")

        except Exception as e:
            self.statusBar().showMessage(f"An unexpected error occurred during 2D export: {e}")

        finally:
            for item, was_visible in items_to_restore.items():
                item.setVisible(was_visible)
            self.scene.setBackgroundBrush(original_background)
            if self.view_2d:
                self.view_2d.viewport().update()

    def export_3d_png(self):
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to export.", 2000)
            return

        # default filename: match XYZ/MOL naming (use base name without suffix)
        default_name = "untitled"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}"
        except Exception:
            default_name = "untitled"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
        except Exception:
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(self, "Export 3D as PNG", default_path, "PNG Files (*.png)")
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for current background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.Yes)

        if reply == QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = (reply == QMessageBox.StandardButton.Yes)

        try:
            self.plotter.screenshot(filePath, transparent_background=is_transparent)
            self.statusBar().showMessage(f"3D view exported to {filePath}", 3000)
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting 3D PNG: {e}")


    def open_periodic_table_dialog(self):
        dialog=PeriodicTableDialog(self); dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action=self.tool_group.checkedAction()
        if checked_action: self.tool_group.setExclusive(False); checked_action.setChecked(False); self.tool_group.setExclusive(True)
        dialog.exec()

    def set_atom_from_periodic_table(self, symbol): 
        self.set_mode(f'atom_{symbol}')

   
    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")
        
        # 最初に既存の化学的問題フラグをクリア
        self.scene.clear_all_problem_flags()
        
        # 2Dエディタに原子が存在しない場合
        if not self.data.atoms:
            self.statusBar().showMessage("Error: No atoms to optimize.")
            return
        
        mol = self.data.to_rdkit_mol()
        if mol is None or mol.GetNumAtoms() == 0:
            # RDKit変換が失敗した場合は化学的問題をチェック
            self.check_chemistry_problems_fallback()
            return

        try:
            # 安定版：原子IDとRDKit座標の確実なマッピング
            view_center = self.view_2d.mapToScene(self.view_2d.viewport().rect().center())
            new_positions_map = {}
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            for rdkit_atom in mol.GetAtoms():
                original_id = rdkit_atom.GetIntProp("_original_atom_id")
                new_positions_map[original_id] = conf.GetAtomPosition(rdkit_atom.GetIdx())

            if not new_positions_map:
                self.statusBar().showMessage("Optimization failed to generate coordinates."); return

            target_atom_items = [self.data.atoms[atom_id]['item'] for atom_id in new_positions_map.keys() if atom_id in self.data.atoms and 'item' in self.data.atoms[atom_id]]
            if not target_atom_items:
                self.statusBar().showMessage("Error: Atom items not found for optimized atoms."); return

            # 元の図形の中心を維持
            #original_center_x = sum(item.pos().x() for item in target_atom_items) / len(target_atom_items)
            #original_center_y = sum(item.pos().y() for item in target_atom_items) / len(target_atom_items)

            positions = list(new_positions_map.values())
            rdkit_cx = sum(p.x for p in positions) / len(positions)
            rdkit_cy = sum(p.y for p in positions) / len(positions)

            SCALE = 50.0

            # 新しい座標を適用
            for atom_id, rdkit_pos in new_positions_map.items():
                if atom_id in self.data.atoms:
                    item = self.data.atoms[atom_id]['item']
                    sx = ((rdkit_pos.x - rdkit_cx) * SCALE) + view_center.x()
                    sy = (-(rdkit_pos.y - rdkit_cy) * SCALE) + view_center.y()
                    new_scene_pos = QPointF(sx, sy)
                    item.setPos(new_scene_pos)
                    self.data.atoms[atom_id]['pos'] = new_scene_pos

            # 最終的な座標に基づき、全ての結合表示を一度に更新
            # Guard against partially-deleted Qt wrappers: skip items that
            # SIP reports as deleted or which are no longer in a scene.
            for bond_data in self.data.bonds.values():
                item = bond_data.get('item') if bond_data else None
                if not item:
                    continue
                try:
                    # If SIP is available, skip wrappers whose C++ object is gone
                    if sip_isdeleted_safe(item):
                        continue
                except Exception:
                    # If the sip check fails, continue with other lightweight guards
                    pass
                try:
                    sc = None
                    try:
                        sc = item.scene() if hasattr(item, 'scene') else None
                    except Exception:
                        sc = None
                    if sc is None:
                        continue
                    try:
                        item.update_position()
                    except Exception:
                        # Best-effort: skip any bond items that raise when updating
                        continue
                except Exception:
                    continue

            # 重なり解消ロジックを実行
            self. resolve_overlapping_groups()
            
            # 測定ラベルの位置を更新
            self.update_2d_measurement_labels()
            
            # シーン全体の再描画を要求
            self.scene.update()

            self.statusBar().showMessage("2D structure optimization successful.")
            self.push_undo_state()

        except Exception as e:
            self.statusBar().showMessage(f"Error during 2D optimization: {e}")
        finally:
            self.view_2d.setFocus()

    def resolve_overlapping_groups(self):
        """
        誤差範囲で完全に重なっている原子のグループを検出し、
        IDが大きい方のフラグメントを左下に平行移動して解消する。
        """

        # --- パラメータ設定 ---
        # 重なっているとみなす距離の閾値。構造に合わせて調整してください。
        OVERLAP_THRESHOLD = 0.5  
        # 左下へ移動させる距離。
        MOVE_DISTANCE = 20

        # self.data.atoms.values() から item を安全に取得
        all_atom_items = [
            data['item'] for data in self.data.atoms.values() 
            if data and 'item' in data
        ]

        if len(all_atom_items) < 2:
            return

        # --- ステップ1: 重なっている原子ペアを全てリストアップ ---
        overlapping_pairs = []
        for item1, item2 in itertools.combinations(all_atom_items, 2):
            # 結合で直接結ばれているペアは重なりと見なさない
            if self.scene.find_bond_between(item1, item2):
                continue

            dist = QLineF(item1.pos(), item2.pos()).length()
            if dist < OVERLAP_THRESHOLD:
                overlapping_pairs.append((item1, item2))

        if not overlapping_pairs:
            self.statusBar().showMessage("No overlapping atoms found.", 2000)
            return

        # --- ステップ2: Union-Findアルゴリズムで重なりグループを構築 ---
        # 各原子がどのグループに属するかを管理する
        parent = {item.atom_id: item.atom_id for item in all_atom_items}

        def find_set(atom_id):
            # atom_idが属するグループの代表（ルート）を見つける
            if parent[atom_id] == atom_id:
                return atom_id
            parent[atom_id] = find_set(parent[atom_id])  # 経路圧縮による最適化
            return parent[atom_id]

        def unite_sets(id1, id2):
            # 2つの原子が属するグループを統合する
            root1 = find_set(id1)
            root2 = find_set(id2)
            if root1 != root2:
                parent[root2] = root1

        for item1, item2 in overlapping_pairs:
            unite_sets(item1.atom_id, item2.atom_id)

        # --- ステップ3: グループごとに移動計画を立てる ---
        # 同じ代表を持つ原子でグループを辞書にまとめる
        groups_by_root = {}
        for item in all_atom_items:
            root_id = find_set(item.atom_id)
            if root_id not in groups_by_root:
                groups_by_root[root_id] = []
            groups_by_root[root_id].append(item.atom_id)

        move_operations = []
        processed_roots = set()

        for root_id, group_atom_ids in groups_by_root.items():
            # 処理済みのグループや、メンバーが1つしかないグループはスキップ
            if root_id in processed_roots or len(group_atom_ids) < 2:
                continue
            processed_roots.add(root_id)

            # 3a: グループを、結合に基づいたフラグメントに分割する (BFSを使用)
            fragments = []
            visited_in_group = set()
            group_atom_ids_set = set(group_atom_ids)

            for atom_id in group_atom_ids:
                if atom_id not in visited_in_group:
                    current_fragment = set()
                    q = deque([atom_id])
                    visited_in_group.add(atom_id)
                    current_fragment.add(atom_id)

                    while q:
                        current_id = q.popleft()
                        # 隣接リスト self.adjacency_list があれば、ここでの探索が高速になります
                        for neighbor_id in self.data.adjacency_list.get(current_id, []):
                            if neighbor_id in group_atom_ids_set and neighbor_id not in visited_in_group:
                                visited_in_group.add(neighbor_id)
                                current_fragment.add(neighbor_id)
                                q.append(neighbor_id)
                    fragments.append(current_fragment)

            if len(fragments) < 2:
                continue  # 複数のフラグメントが重なっていない場合

            # 3b: 移動するフラグメントを決定する
            # このグループの重なりの原因となった代表ペアを一つ探す
            rep_item1, rep_item2 = None, None
            for i1, i2 in overlapping_pairs:
                if find_set(i1.atom_id) == root_id:
                    rep_item1, rep_item2 = i1, i2
                    break

            if not rep_item1: continue

            # 代表ペアがそれぞれどのフラグメントに属するかを見つける
            frag1 = next((f for f in fragments if rep_item1.atom_id in f), None)
            frag2 = next((f for f in fragments if rep_item2.atom_id in f), None)

            # 同一フラグメント内の重なりなどはスキップ
            if not frag1 or not frag2 or frag1 == frag2:
                continue

            # 仕様: IDが大きい方の原子が含まれるフラグメントを動かす
            if rep_item1.atom_id > rep_item2.atom_id:
                ids_to_move = frag1
            else:
                ids_to_move = frag2

            # 3c: 移動計画を作成
            translation_vector = QPointF(-MOVE_DISTANCE, MOVE_DISTANCE)  # 左下方向へのベクトル
            move_operations.append((ids_to_move, translation_vector))

        # --- ステップ4: 計画された移動を一度に実行 ---
        if not move_operations:
            self.statusBar().showMessage("No actionable overlaps found.", 2000)
            return

        for group_ids, vector in move_operations:
            for atom_id in group_ids:
                item = self.data.atoms[atom_id]['item']
                new_pos = item.pos() + vector
                item.setPos(new_pos)
                self.data.atoms[atom_id]['pos'] = new_pos

        # --- ステップ5: 表示と状態を更新 ---
        for bond_data in self.data.bonds.values():
            item = bond_data.get('item') if bond_data else None
            if not item:
                continue
            try:
                if sip_isdeleted_safe(item):
                    continue
            except Exception:
                pass
            try:
                sc = None
                try:
                    sc = item.scene() if hasattr(item, 'scene') else None
                except Exception:
                    sc = None
                if sc is None:
                    continue
                try:
                    item.update_position()
                except Exception:
                    continue
            except Exception:
                continue
        
        # 重なり解消後に測定ラベルの位置を更新
        self.update_2d_measurement_labels()
        
        self.scene.update()
        self.push_undo_state()
        self.statusBar().showMessage("Resolved overlapping groups.", 2000)


    def adjust_molecule_positions_to_avoid_collisions(self, mol, frags):
        """
        複数分子の位置を調整して、衝突を回避する（バウンディングボックス最適化版）
        """
        if len(frags) <= 1:
            return
        
        conf = mol.GetConformer()
        pt = Chem.GetPeriodicTable()
        
        # --- 1. 各フラグメントの情報（原子インデックス、VDW半径）を事前計算 ---
        frag_info = []
        for frag_indices in frags:
            positions = []
            vdw_radii = []
            for idx in frag_indices:
                pos = conf.GetAtomPosition(idx)
                positions.append(np.array([pos.x, pos.y, pos.z]))
                
                atom = mol.GetAtomWithIdx(idx)
                # GetRvdw() はファンデルワールス半径を返す
                vdw_radii.append(pt.GetRvdw(atom.GetAtomicNum()))

            positions_np = np.array(positions)
            vdw_radii_np = np.array(vdw_radii)
            
            # このフラグメントで最大のVDW半径を計算（ボックスのマージンとして使用）
            max_vdw = np.max(vdw_radii_np) if len(vdw_radii_np) > 0 else 0.0

            frag_info.append({
                'indices': frag_indices,
                'centroid': np.mean(positions_np, axis=0),
                'positions_np': positions_np, # Numpy配列として保持
                'vdw_radii_np': vdw_radii_np,  # Numpy配列として保持
                'max_vdw_radius': max_vdw,
                'bbox_min': np.zeros(3), # 後で計算
                'bbox_max': np.zeros(3)  # 後で計算
            })
        
        # --- 2. 衝突判定のパラメータ ---
        collision_scale = 1.2  # VDW半径の120%
        max_iterations = 100
        moved = True
        iteration = 0
        
        while moved and iteration < max_iterations:
            moved = False
            iteration += 1
            
            # --- 3. フラグメントのバウンディングボックスを毎イテレーション更新 ---
            for i in range(len(frag_info)):
                # 現在の座標からボックスを再計算
                current_positions = []
                for idx in frag_info[i]['indices']:
                    pos = conf.GetAtomPosition(idx)
                    current_positions.append([pos.x, pos.y, pos.z])
                
                positions_np = np.array(current_positions)
                frag_info[i]['positions_np'] = positions_np # 座標情報を更新
                
                # VDW半径とスケールを考慮したマージンを計算
                # (最大VDW半径 * スケール) をマージンとして使う
                margin = frag_info[i]['max_vdw_radius'] * collision_scale
                
                frag_info[i]['bbox_min'] = np.min(positions_np, axis=0) - margin
                frag_info[i]['bbox_max'] = np.max(positions_np, axis=0) + margin

            # --- 4. 衝突判定ループ ---
            for i in range(len(frag_info)):
                for j in range(i + 1, len(frag_info)):
                    frag_i = frag_info[i]
                    frag_j = frag_info[j]
                    
                    # === バウンディングボックス判定 ===
                    # 2つのボックスが重なっているかチェック (AABB交差判定)
                    # X, Y, Zの各軸で重なりをチェック
                    overlap_x = (frag_i['bbox_min'][0] <= frag_j['bbox_max'][0] and frag_i['bbox_max'][0] >= frag_j['bbox_min'][0])
                    overlap_y = (frag_i['bbox_min'][1] <= frag_j['bbox_max'][1] and frag_i['bbox_max'][1] >= frag_j['bbox_min'][1])
                    overlap_z = (frag_i['bbox_min'][2] <= frag_j['bbox_max'][2] and frag_i['bbox_max'][2] >= frag_j['bbox_min'][2])
                    
                    # ボックスがX, Y, Zのいずれかの軸で離れている場合、原子間の詳細なチェックをスキップ
                    if not (overlap_x and overlap_y and overlap_z):
                        continue
                    # =================================

                    # ボックスが重なっている場合のみ、高コストな原子間の総当たりチェックを実行
                    total_push_vector = np.zeros(3)
                    collision_count = 0
                    
                    # 事前計算したNumpy配列を使用
                    positions_i = frag_i['positions_np']
                    positions_j = frag_j['positions_np']
                    vdw_i_all = frag_i['vdw_radii_np']
                    vdw_j_all = frag_j['vdw_radii_np']

                    for k, idx_i in enumerate(frag_i['indices']):
                        pos_i = positions_i[k]
                        vdw_i = vdw_i_all[k]
                        
                        for l, idx_j in enumerate(frag_j['indices']):
                            pos_j = positions_j[l]
                            vdw_j = vdw_j_all[l]
                            
                            distance_vec = pos_i - pos_j
                            distance_sq = np.dot(distance_vec, distance_vec) # 平方根を避けて高速化
                            
                            min_distance = (vdw_i + vdw_j) * collision_scale
                            min_distance_sq = min_distance * min_distance
                            
                            if distance_sq < min_distance_sq and distance_sq > 0.0001:
                                distance = np.sqrt(distance_sq)
                                push_direction = distance_vec / distance
                                push_magnitude = (min_distance - distance) / 2 # 押し出し量は半分ずつ
                                total_push_vector += push_direction * push_magnitude
                                collision_count += 1
                    
                    if collision_count > 0:
                        # 平均的な押し出しベクトルを適用
                        avg_push_vector = total_push_vector / collision_count
                        
                        # Conformerの座標を更新
                        for idx in frag_i['indices']:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos + avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())
                        
                        for idx in frag_j['indices']:
                            pos = np.array(conf.GetAtomPosition(idx))
                            new_pos = pos - avg_push_vector
                            conf.SetAtomPosition(idx, new_pos.tolist())
                        
                        moved = True
                        # (この移動により、このイテレーションで使う frag_info の座標キャッシュが古くなりますが、
                        #  次のイテレーションの最初でボックスと共に再計算されるため問題ありません)

    def draw_molecule_3d(self, mol):
        """3D 分子を描画し、軸アクターの参照をクリアする（軸の再制御は apply_3d_settings に任せる）"""
        
        # 測定選択をクリア（分子が変更されたため）
        if hasattr(self, 'measurement_mode'):
            self.clear_measurement_selection()
        
        # 色情報追跡のための辞書を初期化
        if not hasattr(self, '_3d_color_map'):
            self._3d_color_map = {}
        self._3d_color_map.clear()
        
        # 1. カメラ状態とクリア
        camera_state = self.plotter.camera.copy()

        # **残留防止のための強制削除**
        if self.axes_actor is not None:
            try:
                self.plotter.remove_actor(self.axes_actor)
            except Exception:
                pass 
            self.axes_actor = None

        self.plotter.clear()
            
        # 2. 背景色の設定
        self.plotter.set_background(self.settings.get('background_color', '#4f4f4f'))

        # 3. mol が None または原子数ゼロの場合は、背景と軸のみで終了
        if mol is None or mol.GetNumAtoms() == 0:
            self.atom_actor = None
            self.current_mol = None
            self.plotter.render()
            return
            
        # 4. ライティングの設定
        is_lighting_enabled = self.settings.get('lighting_enabled', True)

        if is_lighting_enabled:
            light = pv.Light(
                position=(1, 1, 2),
                light_type='cameralight',
                intensity=self.settings.get('light_intensity', 1.2)
            )
            self.plotter.add_light(light)
            
        # 5. 分子描画ロジック
        conf = mol.GetConformer()

        self.atom_positions_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

        sym = [a.GetSymbol() for a in mol.GetAtoms()]
        col = np.array([CPK_COLORS_PV.get(s, [0.5, 0.5, 0.5]) for s in sym])

        # スタイルに応じて原子の半径を設定（設定から読み込み）
        if self.current_3d_style == 'cpk':
            atom_scale = self.settings.get('cpk_atom_scale', 1.0)
            resolution = self.settings.get('cpk_resolution', 32)
            rad = np.array([pt.GetRvdw(pt.GetAtomicNumber(s)) * atom_scale for s in sym])
        elif self.current_3d_style == 'wireframe':
            # Wireframeでは原子を描画しないので、この設定は実際には使用されない
            resolution = self.settings.get('wireframe_resolution', 6)
            rad = np.array([0.01 for s in sym])  # 極小値（使用されない）
        elif self.current_3d_style == 'stick':
            atom_radius = self.settings.get('stick_atom_radius', 0.15)
            resolution = self.settings.get('stick_resolution', 16)
            rad = np.array([atom_radius for s in sym])
        else:  # ball_and_stick
            atom_scale = self.settings.get('ball_stick_atom_scale', 1.0)
            resolution = self.settings.get('ball_stick_resolution', 16)
            rad = np.array([VDW_RADII.get(s, 0.4) * atom_scale for s in sym])

        self.glyph_source = pv.PolyData(self.atom_positions_3d)
        self.glyph_source['colors'] = col
        self.glyph_source['radii'] = rad

        # メッシュプロパティを共通で定義
        mesh_props = dict(
            smooth_shading=True,
            specular=self.settings.get('specular', 0.2),
            specular_power=self.settings.get('specular_power', 20),
            lighting=is_lighting_enabled,
        )

        # Wireframeスタイルの場合は原子を描画しない
        if self.current_3d_style != 'wireframe':
            glyphs = self.glyph_source.glyph(scale='radii', geom=pv.Sphere(radius=1.0, theta_resolution=resolution, phi_resolution=resolution), orient=False)

            if is_lighting_enabled:
                self.atom_actor = self.plotter.add_mesh(glyphs, scalars='colors', rgb=True, **mesh_props)
            else:
                self.atom_actor = self.plotter.add_mesh(
                    glyphs, scalars='colors', rgb=True, 
                    style='surface', show_edges=True, edge_color='grey',
                    **mesh_props
                )
                self.atom_actor.GetProperty().SetEdgeOpacity(0.3)
            
            # 原子の色情報を記録
            for i, atom_color in enumerate(col):
                atom_rgb = [int(c * 255) for c in atom_color]
                self._3d_color_map[f'atom_{i}'] = atom_rgb


        # ボンドの描画（ball_and_stick、wireframe、stickで描画）
        if self.current_3d_style in ['ball_and_stick', 'wireframe', 'stick']:
            # スタイルに応じてボンドの太さと解像度を設定（設定から読み込み）
            if self.current_3d_style == 'wireframe':
                cyl_radius = self.settings.get('wireframe_bond_radius', 0.01)
                bond_resolution = self.settings.get('wireframe_resolution', 6)
            elif self.current_3d_style == 'stick':
                cyl_radius = self.settings.get('stick_bond_radius', 0.15)
                bond_resolution = self.settings.get('stick_resolution', 16)
            else:  # ball_and_stick
                cyl_radius = self.settings.get('ball_stick_bond_radius', 0.1)
                bond_resolution = self.settings.get('ball_stick_resolution', 16)
            
            bond_counter = 0  # 結合の個別識別用
            
            # Ball and Stick用のシリンダーリストを準備（高速化のため）
            if self.current_3d_style == 'ball_and_stick':
                bond_cylinders = []
                # Compute the configured grey/uniform bond color for Ball & Stick
                try:
                    bs_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                    q = QColor(bs_hex)
                    bs_bond_rgb = [q.red(), q.green(), q.blue()]
                except Exception:
                    bs_bond_rgb = [127, 127, 127]
            
            for bond in mol.GetBonds():
                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()
                sp = np.array(conf.GetAtomPosition(begin_atom_idx))
                ep = np.array(conf.GetAtomPosition(end_atom_idx))
                bt = bond.GetBondType()
                c = (sp + ep) / 2
                d = ep - sp
                h = np.linalg.norm(d)
                if h == 0: continue

                # ボンドの色を原子の色から決定（各半分で異なる色）
                begin_color = col[begin_atom_idx]
                end_color = col[end_atom_idx]
                
                # 結合の色情報を記録
                begin_color_rgb = [int(c * 255) for c in begin_color]
                end_color_rgb = [int(c * 255) for c in end_color]

                # UI応答性維持のためイベント処理
                QApplication.processEvents()
                if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                    if self.current_3d_style == 'ball_and_stick':
                        # Ball and stickは全結合をまとめて処理（高速化）
                        cyl = pv.Cylinder(center=c, direction=d, radius=cyl_radius, height=h, resolution=bond_resolution)
                        bond_cylinders.append(cyl)
                        self._3d_color_map[f'bond_{bond_counter}'] = bs_bond_rgb  # グレー (configurable)
                    else:
                        # その他（stick, wireframe）は中央で色が変わる2つの円柱
                        mid_point = (sp + ep) / 2
                        
                        # 前半（開始原子の色）
                        cyl1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                        actor1 = self.plotter.add_mesh(cyl1, color=begin_color, **mesh_props)
                        self._3d_color_map[f'bond_{bond_counter}_start'] = begin_color_rgb
                        
                        # 後半（終了原子の色）
                        cyl2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                        actor2 = self.plotter.add_mesh(cyl2, color=end_color, **mesh_props)
                        self._3d_color_map[f'bond_{bond_counter}_end'] = end_color_rgb
                else:
                    v1 = d / h
                    # モデルごとの半径ファクターを適用
                    if self.current_3d_style == 'ball_and_stick':
                        double_radius_factor = self.settings.get('ball_stick_double_bond_radius_factor', 0.8)
                        triple_radius_factor = self.settings.get('ball_stick_triple_bond_radius_factor', 0.75)
                    elif self.current_3d_style == 'wireframe':
                        double_radius_factor = self.settings.get('wireframe_double_bond_radius_factor', 1.0)
                        triple_radius_factor = self.settings.get('wireframe_triple_bond_radius_factor', 0.75)
                    elif self.current_3d_style == 'stick':
                        double_radius_factor = self.settings.get('stick_double_bond_radius_factor', 0.60)
                        triple_radius_factor = self.settings.get('stick_triple_bond_radius_factor', 0.40)
                    else:
                        double_radius_factor = 1.0
                        triple_radius_factor = 0.75
                    r = cyl_radius * 0.8  # fallback, will be overridden below
                    # 設定からオフセットファクターを取得（モデルごと）
                    if self.current_3d_style == 'ball_and_stick':
                        double_offset_factor = self.settings.get('ball_stick_double_bond_offset_factor', 2.0)
                        triple_offset_factor = self.settings.get('ball_stick_triple_bond_offset_factor', 2.0)
                    elif self.current_3d_style == 'wireframe':
                        double_offset_factor = self.settings.get('wireframe_double_bond_offset_factor', 3.0)
                        triple_offset_factor = self.settings.get('wireframe_triple_bond_offset_factor', 3.0)
                    elif self.current_3d_style == 'stick':
                        double_offset_factor = self.settings.get('stick_double_bond_offset_factor', 1.5)
                        triple_offset_factor = self.settings.get('stick_triple_bond_offset_factor', 1.0)
                    else:
                        double_offset_factor = 2.0
                        triple_offset_factor = 2.0
                    s = cyl_radius * 2.0  # デフォルト値

                    if bt == Chem.rdchem.BondType.DOUBLE:
                        r = cyl_radius * double_radius_factor
                        # 二重結合の場合、結合している原子の他の結合を考慮してオフセット方向を決定
                        off_dir = self._calculate_double_bond_offset(mol, bond, conf)
                        # 設定から二重結合のオフセットファクターを適用
                        s_double = cyl_radius * double_offset_factor
                        c1, c2 = c + off_dir * (s_double / 2), c - off_dir * (s_double / 2)
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickは全結合をまとめて処理（高速化）
                            cyl1 = pv.Cylinder(center=c1, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c2, direction=d, radius=r, height=h, resolution=bond_resolution)
                            bond_cylinders.extend([cyl1, cyl2])
                            self._3d_color_map[f'bond_{bond_counter}_1'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2'] = bs_bond_rgb
                        else:
                            # その他（stick, wireframe）は中央で色が変わる
                            mid_point = (sp + ep) / 2
                            
                            # 第一の結合線（前半・後半）
                            cyl1_1 = pv.Cylinder(center=(sp + mid_point) / 2 + off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl1_2 = pv.Cylinder(center=(mid_point + ep) / 2 + off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl1_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                            
                            # 第二の結合線（前半・後半）
                            cyl2_1 = pv.Cylinder(center=(sp + mid_point) / 2 - off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl2_2 = pv.Cylinder(center=(mid_point + ep) / 2 - off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl2_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl2_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_2_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_end'] = end_color_rgb
                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        r = cyl_radius * triple_radius_factor
                        # 三重結合
                        v_arb = np.array([0, 0, 1])
                        if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                        off_dir = np.cross(v1, v_arb)
                        off_dir /= np.linalg.norm(off_dir)
                        
                        # 設定から三重結合のオフセットファクターを適用
                        s_triple = cyl_radius * triple_offset_factor
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickは全結合をまとめて処理（高速化）
                            cyl1 = pv.Cylinder(center=c, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c + off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl3 = pv.Cylinder(center=c - off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            bond_cylinders.extend([cyl1, cyl2, cyl3])
                            self._3d_color_map[f'bond_{bond_counter}_1'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_3'] = bs_bond_rgb
                        else:
                            # その他（stick, wireframe）は中央で色が変わる
                            mid_point = (sp + ep) / 2
                            
                            # 中央の結合線（前半・後半）
                            cyl1_1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl1_2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl1_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                            
                            # 上側の結合線（前半・後半）
                            cyl2_1 = pv.Cylinder(center=(sp + mid_point) / 2 + off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl2_2 = pv.Cylinder(center=(mid_point + ep) / 2 + off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl2_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl2_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_2_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_end'] = end_color_rgb
                            
                            # 下側の結合線（前半・後半）
                            cyl3_1 = pv.Cylinder(center=(sp + mid_point) / 2 - off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl3_2 = pv.Cylinder(center=(mid_point + ep) / 2 - off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl3_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl3_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_3_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_3_end'] = end_color_rgb

                bond_counter += 1
            
            # Ball and Stick用：全結合をまとめて一括描画（高速化）
            if self.current_3d_style == 'ball_and_stick' and bond_cylinders:
                # 全シリンダーを結合してMultiBlockを作成
                combined_bonds = pv.MultiBlock(bond_cylinders)
                combined_mesh = combined_bonds.combine()
                
                # 一括でグレーで描画
                # Use the configured Ball & Stick bond color (hex) for the combined bonds
                try:
                    bs_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                    q = QColor(bs_hex)
                    # Use normalized RGB for pyvista (r,g,b) floats in [0,1]
                    bond_color = (q.redF(), q.greenF(), q.blueF())
                    bond_actor = self.plotter.add_mesh(combined_mesh, color=bond_color, **mesh_props)
                except Exception:
                    bond_actor = self.plotter.add_mesh(combined_mesh, color='grey', **mesh_props)
                
                # まとめて色情報を記録
                self._3d_color_map['bonds_combined'] = bs_bond_rgb

        if getattr(self, 'show_chiral_labels', False):
            try:
                # 3D座標からキラル中心を計算
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if chiral_centers:
                    pts, labels = [], []
                    z_off = 0
                    for idx, lbl in chiral_centers:
                        coord = self.atom_positions_3d[idx].copy(); coord[2] += z_off
                        pts.append(coord); labels.append(lbl if lbl is not None else '?')
                    try: self.plotter.remove_actor('chiral_labels')
                    except Exception: pass
                    self.plotter.add_point_labels(np.array(pts), labels, font_size=20, point_size=0, text_color='blue', name='chiral_labels', always_visible=True, tolerance=0.01, show_points=False)
            except Exception as e: self.statusBar().showMessage(f"3D chiral label drawing error: {e}")

        # E/Zラベルも表示
        if getattr(self, 'show_chiral_labels', False):
            try:
                self.show_ez_labels_3d(mol)
            except Exception as e: 
                self.statusBar().showMessage(f"3D E/Z label drawing error: {e}")

        self.plotter.camera = camera_state

        # Ensure the underlying VTK camera's parallel/projection flag matches
        # the saved application setting. draw_molecule_3d restores a PyVista
        # camera object which may not propagate the ParallelProjection flag
        # to the VTK renderer camera; enforce it here to guarantee the
        # projection mode selected in settings actually takes effect.
        try:
            proj_mode = self.settings.get('projection_mode', 'Perspective')
            if hasattr(self.plotter, 'renderer') and hasattr(self.plotter.renderer, 'GetActiveCamera'):
                vcam = self.plotter.renderer.GetActiveCamera()
                if vcam:
                    if proj_mode == 'Orthographic':
                        vcam.SetParallelProjection(True)
                    else:
                        vcam.SetParallelProjection(False)
                    try:
                        # Force a render so the change is visible immediately
                        self.plotter.render()
                    except Exception:
                        pass
        except Exception:
            pass
        
        # AtomIDまたは他の原子情報が表示されている場合は再表示
        if hasattr(self, 'atom_info_display_mode') and self.atom_info_display_mode is not None:
            self.show_all_atom_info()
        
        # メニューテキストと状態を現在の分子の種類に応じて更新
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

    def _calculate_double_bond_offset(self, mol, bond, conf):
        """
        二重結合のオフセット方向を計算する。
        結合している原子の他の結合を考慮して、平面的になるようにする。
        """
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        
        begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
        end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
        
        bond_vec = end_pos - begin_pos
        bond_length = np.linalg.norm(bond_vec)
        if bond_length == 0:
            # フォールバック: Z軸基準
            return np.array([0, 0, 1])
        
        bond_unit = bond_vec / bond_length
        
        # 両端の原子の隣接原子を調べる
        begin_neighbors = []
        end_neighbors = []
        
        for neighbor in begin_atom.GetNeighbors():
            if neighbor.GetIdx() != bond.GetEndAtomIdx():
                neighbor_pos = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                begin_neighbors.append(neighbor_pos)
        
        for neighbor in end_atom.GetNeighbors():
            if neighbor.GetIdx() != bond.GetBeginAtomIdx():
                neighbor_pos = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                end_neighbors.append(neighbor_pos)
        
        # 平面の法線ベクトルを計算
        normal_candidates = []
        
        # 開始原子の隣接原子から平面を推定
        if len(begin_neighbors) >= 1:
            for neighbor_pos in begin_neighbors:
                vec_to_neighbor = neighbor_pos - begin_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # bond_vec と neighbor_vec の外積が平面の法線
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)
        
        # 終了原子の隣接原子から平面を推定
        if len(end_neighbors) >= 1:
            for neighbor_pos in end_neighbors:
                vec_to_neighbor = neighbor_pos - end_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # bond_vec と neighbor_vec の外積が平面の法線
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)
        
        # 複数の法線ベクトルがある場合は平均を取る
        if normal_candidates:
            # 方向を統一するため、最初のベクトルとの内積が正になるように調整
            reference_normal = normal_candidates[0]
            aligned_normals = []
            
            for normal in normal_candidates:
                if np.dot(normal, reference_normal) < 0:
                    normal = -normal
                aligned_normals.append(normal)
            
            avg_normal = np.mean(aligned_normals, axis=0)
            norm_length = np.linalg.norm(avg_normal)
            if norm_length > 1e-6:
                avg_normal /= norm_length
                
                # 法線ベクトルと結合ベクトルに垂直な方向を二重結合のオフセット方向とする
                offset_dir = np.cross(bond_unit, avg_normal)
                offset_length = np.linalg.norm(offset_dir)
                if offset_length > 1e-6:
                    return offset_dir / offset_length
        
        # フォールバック: 結合ベクトルに垂直な任意の方向
        v_arb = np.array([0, 0, 1])
        if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0):
            v_arb = np.array([0, 1, 0])
        
        off_dir = np.cross(bond_unit, v_arb)
        off_dir /= np.linalg.norm(off_dir)
        return off_dir

    def show_ez_labels_3d(self, mol):
        """3DビューでE/Zラベルを表示する（RDKitのステレオ化学判定を使用）"""
        if not mol:
            return
        
        try:
            # 既存のE/Zラベルを削除
            self.plotter.remove_actor('ez_labels')
        except:
            pass
        
        pts, labels = [], []
        
        # 3D座標が存在するかチェック
        if mol.GetNumConformers() == 0:
            return
            
        conf = mol.GetConformer()
        
        # RDKitに3D座標からステレオ化学を計算させる
        try:
            # 3D座標からステレオ化学を再計算
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        except:
            pass
        
        # 二重結合でRDKitが判定したE/Z立体化学を表示
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                    # 結合の中心座標を計算
                    begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
                    end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
                    center_pos = (begin_pos + end_pos) / 2
                    
                    # RDKitの判定結果を使用
                    label = 'E' if stereo == Chem.BondStereo.STEREOE else 'Z'
                    pts.append(center_pos)
                    labels.append(label)
        
        if pts and labels:
            self.plotter.add_point_labels(
                np.array(pts), 
                labels, 
                font_size=18,
                point_size=0,
                text_color='darkgreen',  # 暗い緑色
                name='ez_labels',
                always_visible=True,
                tolerance=0.01,
                show_points=False
            )


    def toggle_chiral_labels_display(self, checked):
        """Viewメニューのアクションに応じてキラルラベル表示を切り替える"""
        self.show_chiral_labels = checked
        
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol) 
        
        if checked:
            self.statusBar().showMessage("Chiral labels: will be (re)computed after Convert→3D.")
        else:
            self.statusBar().showMessage("Chiral labels disabled.")


    def update_chiral_labels(self):
        """分子のキラル中心を計算し、2Dビューの原子アイテムにR/Sラベルを設定/解除する
        ※ 可能なら 3D（self.current_mol）を優先して計算し、なければ 2D から作った RDKit 分子を使う。
        """
        # まず全てのアイテムからラベルをクリア
        for atom_data in self.data.atoms.values():
            if atom_data.get('item'):
                atom_data['item'].chiral_label = None

        if not self.show_chiral_labels:
            self.scene.update()
            return

        # 3D の RDKit Mol（コンフォマーを持つもの）を使う
        mol_for_chirality = None
        if getattr(self, 'current_mol', None) is not None:
            mol_for_chirality = self.current_mol
        else:
            return

        if mol_for_chirality is None or mol_for_chirality.GetNumAtoms() == 0:
            self.scene.update()
            return

        try:
            # --- 重要：3D コンフォマーがあるなら、それを使って原子のキラルタグを割り当てる ---
            if mol_for_chirality.GetNumConformers() > 0:
                # confId=0（最初のコンフォマー）を指定して、原子のキラリティータグを3D座標由来で設定
                try:
                    Chem.AssignAtomChiralTagsFromStructure(mol_for_chirality, confId=0)
                except Exception:
                    # 古い RDKit では関数が無い場合があるので（念のため保護）
                    pass

            # RDKit の通常の stereochemistry 割当（念のため）
            #Chem.AssignStereochemistry(mol_for_chirality, cleanIt=True, force=True, flagPossibleStereoCenters=True)

            # キラル中心の取得（(idx, 'R'/'S'/'?') のリスト）
            chiral_centers = Chem.FindMolChiralCenters(mol_for_chirality, includeUnassigned=True)

            # RDKit atom index -> エディタ側 atom_id へのマッピング
            rdkit_idx_to_my_id = {}
            for atom in mol_for_chirality.GetAtoms():
                if atom.HasProp("_original_atom_id"):
                    rdkit_idx_to_my_id[atom.GetIdx()] = atom.GetIntProp("_original_atom_id")

            # 見つかったキラル中心を対応する AtomItem に設定
            for idx, label in chiral_centers:
                if idx in rdkit_idx_to_my_id:
                    atom_id = rdkit_idx_to_my_id[idx]
                    if atom_id in self.data.atoms and self.data.atoms[atom_id].get('item'):
                        # 'R' / 'S' / '?'
                        self.data.atoms[atom_id]['item'].chiral_label = label

        except Exception as e:
            self.statusBar().showMessage(f"Update chiral labels error: {e}")

        # 最後に 2D シーンを再描画
        self.scene.update()

    def toggle_atom_info_display(self, mode):
        """原子情報表示モードを切り替える"""
        # 現在の表示をクリア
        self.clear_all_atom_info_labels()
        
        # 同じモードが選択された場合はOFFにする
        if self.atom_info_display_mode == mode:
            self.atom_info_display_mode = None
            # 全てのアクションのチェックを外す
            self.show_atom_id_action.setChecked(False)
            self.show_rdkit_id_action.setChecked(False)
            self.show_atom_coords_action.setChecked(False)
            self.show_atom_symbol_action.setChecked(False)
            self.statusBar().showMessage("Atom info display disabled.")
        else:
            # 新しいモードを設定
            self.atom_info_display_mode = mode
            # 該当するアクションのみチェック
            self.show_atom_id_action.setChecked(mode == 'id')
            self.show_rdkit_id_action.setChecked(mode == 'rdkit_id')
            self.show_atom_coords_action.setChecked(mode == 'coords')
            self.show_atom_symbol_action.setChecked(mode == 'symbol')
            
            mode_names = {'id': 'Atom ID', 'rdkit_id': 'RDKit Index', 'coords': 'Coordinates', 'symbol': 'Element Symbol'}
            self.statusBar().showMessage(f"Displaying: {mode_names[mode]}")
            
            # すべての原子に情報を表示
            self.show_all_atom_info()

    def is_xyz_derived_molecule(self):
        """現在の分子がXYZファイル由来かどうかを判定"""
        if not self.current_mol:
            return False
        try:
            # 最初の原子がxyz_unique_idプロパティを持っているかチェック
            if self.current_mol.GetNumAtoms() > 0:
                return self.current_mol.GetAtomWithIdx(0).HasProp("xyz_unique_id")
        except Exception:
            pass
        return False

    def has_original_atom_ids(self):
        """現在の分子がOriginal Atom IDsを持っているかどうかを判定"""
        if not self.current_mol:
            return False
        try:
            # いずれかの原子が_original_atom_idプロパティを持っているかチェック
            for atom_idx in range(self.current_mol.GetNumAtoms()):
                atom = self.current_mol.GetAtomWithIdx(atom_idx)
                if atom.HasProp("_original_atom_id"):
                    return True
        except Exception:
            pass
        return False

    def update_atom_id_menu_text(self):
        """原子IDメニューのテキストを現在の分子の種類に応じて更新"""
        if hasattr(self, 'show_atom_id_action'):
            if self.is_xyz_derived_molecule():
                self.show_atom_id_action.setText("Show XYZ Unique ID")
            else:
                self.show_atom_id_action.setText("Show Original ID / Index")

    def update_atom_id_menu_state(self):
        """原子IDメニューの有効/無効状態を更新"""
        if hasattr(self, 'show_atom_id_action'):
            has_original_ids = self.has_original_atom_ids()
            has_xyz_ids = self.is_xyz_derived_molecule()
            
            # Original IDまたはXYZ IDがある場合のみ有効化
            self.show_atom_id_action.setEnabled(has_original_ids or has_xyz_ids)
            
            # 現在選択されているモードが無効化される場合は解除
            if not (has_original_ids or has_xyz_ids) and self.atom_info_display_mode == 'id':
                self.atom_info_display_mode = None
                self.show_atom_id_action.setChecked(False)
                self.clear_all_atom_info_labels()


    def show_all_atom_info(self):
        """すべての原子に情報を表示"""
        if self.atom_info_display_mode is None or not hasattr(self, 'atom_positions_3d') or self.atom_positions_3d is None:
            return
        
        # 既存のラベルをクリア
        self.clear_all_atom_info_labels()

        # ラベルを表示するためにタイプ別に分けてリストを作る
        rdkit_positions = []
        rdkit_texts = []
        id_positions = []
        id_texts = []
        xyz_positions = []
        xyz_texts = []
        other_positions = []
        other_texts = []

        for atom_idx, pos in enumerate(self.atom_positions_3d):
            # default: skip if no display mode
            if self.atom_info_display_mode is None:
                continue

            if self.atom_info_display_mode == 'id':
                # Original IDがある場合は優先表示、なければXYZのユニークID、最後にRDKitインデックス
                try:
                    if self.current_mol:
                        atom = self.current_mol.GetAtomWithIdx(atom_idx)
                        if atom.HasProp("_original_atom_id"):
                            original_id = atom.GetIntProp("_original_atom_id")
                            # プレフィックスを削除して数値だけ表示
                            id_positions.append(pos)
                            id_texts.append(str(original_id))
                        elif atom.HasProp("xyz_unique_id"):
                            unique_id = atom.GetIntProp("xyz_unique_id")
                            xyz_positions.append(pos)
                            xyz_texts.append(str(unique_id))
                        else:
                            rdkit_positions.append(pos)
                            rdkit_texts.append(str(atom_idx))
                    else:
                        rdkit_positions.append(pos)
                        rdkit_texts.append(str(atom_idx))
                except Exception:
                    rdkit_positions.append(pos)
                    rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == 'rdkit_id':
                rdkit_positions.append(pos)
                rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == 'coords':
                other_positions.append(pos)
                other_texts.append(f"({pos[0]:.2f},{pos[1]:.2f},{pos[2]:.2f})")

            elif self.atom_info_display_mode == 'symbol':
                if self.current_mol:
                    symbol = self.current_mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    other_positions.append(pos)
                    other_texts.append(symbol)
                else:
                    other_positions.append(pos)
                    other_texts.append("?")

            else:
                continue

        # 色の定義（暗めの青/緑/赤）
        rdkit_color = '#003366'   # 暗めの青
        id_color = '#006400'      # 暗めの緑
        xyz_color = '#8B0000'     # 暗めの赤
        other_color = 'black'

        # それぞれのグループごとにラベルを追加し、参照をリストで保持する
        self.current_atom_info_labels = []
        try:
            if rdkit_positions:
                a = self.plotter.add_point_labels(
                    np.array(rdkit_positions), rdkit_texts,
                    point_size=12, font_size=18, text_color=rdkit_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_rdkit'
                )
                self.current_atom_info_labels.append(a)

            if id_positions:
                a = self.plotter.add_point_labels(
                    np.array(id_positions), id_texts,
                    point_size=12, font_size=18, text_color=id_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_id'
                )
                self.current_atom_info_labels.append(a)

            if xyz_positions:
                a = self.plotter.add_point_labels(
                    np.array(xyz_positions), xyz_texts,
                    point_size=12, font_size=18, text_color=xyz_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_xyz'
                )
                self.current_atom_info_labels.append(a)

            if other_positions:
                a = self.plotter.add_point_labels(
                    np.array(other_positions), other_texts,
                    point_size=12, font_size=18, text_color=other_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_other'
                )
                self.current_atom_info_labels.append(a)
        except Exception as e:
            print(f"Error adding atom info labels: {e}")

        # 右上に凡例を表示（既存の凡例は消す）
        try:
            # 古い凡例削除
            if hasattr(self, 'atom_label_legend_names') and self.atom_label_legend_names:
                for nm in self.atom_label_legend_names:
                    try:
                        self.plotter.remove_actor(nm)
                    except:
                        pass
            self.atom_label_legend_names = []

            # 凡例テキストを右上に縦並びで追加（背景なし、太字のみ）
            legend_entries = []
            if rdkit_positions:
                legend_entries.append(('RDKit', rdkit_color, 'legend_rdkit'))
            if id_positions:
                legend_entries.append(('ID', id_color, 'legend_id'))
            if xyz_positions:
                legend_entries.append(('XYZ', xyz_color, 'legend_xyz'))
            # Do not show 'Other' in the legend per UI requirement
            # (other_positions are still labeled in-scene but not listed in the legend)

            # 左下に凡例ラベルを追加（背景なし、太字のみ）
            # Increase spacing to avoid overlapping when short labels like 'RDKit' and 'ID' appear
            spacing = 30
            for i, (label_text, label_color, label_name) in enumerate(legend_entries):
                # 左下基準でy座標を上げる
                # Add a small horizontal offset for very short adjacent labels so they don't visually collide
                y = 0.0 + i * spacing
                x_offset = 0.0
                # If both RDKit and ID are present, nudge the second entry slightly to the right to avoid overlap
                try:
                    if label_text == 'ID' and any(e[0] == 'RDKit' for e in legend_entries):
                        x_offset = 0.06
                except Exception:
                    x_offset = 0.0
                try:
                    actor = self.plotter.add_text(
                        label_text,
                        position=(0.0 + x_offset, y),
                        font_size=12,
                        color=label_color,
                        name=label_name,
                        font='arial'
                    )
                    self.atom_label_legend_names.append(label_name)
                    # 太字のみ設定（背景は設定しない）
                    try:
                        if hasattr(actor, 'GetTextProperty'):
                            tp = actor.GetTextProperty()
                            try:
                                tp.SetBold(True)
                            except Exception:
                                pass
                    except Exception:
                        pass
                except Exception:
                    continue

        except Exception:
            pass

    def clear_all_atom_info_labels(self):
        """すべての原子情報ラベルをクリア"""
        # Remove label actors (may be a single actor, a list, or None)
        try:
            if hasattr(self, 'current_atom_info_labels') and self.current_atom_info_labels:
                if isinstance(self.current_atom_info_labels, (list, tuple)):
                    for a in list(self.current_atom_info_labels):
                        try:
                            self.plotter.remove_actor(a)
                        except:
                            pass
                else:
                    try:
                        self.plotter.remove_actor(self.current_atom_info_labels)
                    except:
                        pass
        except Exception:
            pass
        finally:
            self.current_atom_info_labels = None

        # Remove legend text actors if present
        try:
            if hasattr(self, 'atom_label_legend_names') and self.atom_label_legend_names:
                for nm in list(self.atom_label_legend_names):
                    try:
                        self.plotter.remove_actor(nm)
                    except:
                        pass
        except Exception:
            pass
        finally:
            self.atom_label_legend_names = []

    def setup_3d_hover(self):
        """3Dビューでの表示を設定（常時表示に変更）"""
        if self.atom_info_display_mode is not None:
            self.show_all_atom_info()

    def open_analysis_window(self):
        if self.current_mol:
            dialog = AnalysisWindow(self.current_mol, self, is_xyz_derived=self.is_xyz_derived)
            dialog.exec()
        else:
            self.statusBar().showMessage("Please generate a 3D structure first to show analysis.")

    def closeEvent(self, event):
        # Persist settings on exit only when explicitly modified (deferred save)
        try:
            if getattr(self, 'settings_dirty', False) or self.settings != self.initial_settings:
                self.save_settings()
                self.settings_dirty = False
        except Exception:
            pass
        
        # 未保存の変更がある場合の処理
        if self.has_unsaved_changes:
            reply = QMessageBox.question(
                self, "Unsaved Changes",
                "You have unsaved changes. Do you want to save them?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Yes
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # 保存処理
                self.save_project()
                
                # 保存がキャンセルされた場合は終了もキャンセル
                if self.has_unsaved_changes:
                    event.ignore()
                    return
                    
            elif reply == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
            # No の場合はそのまま終了処理へ
        
        # 開いているすべてのダイアログウィンドウを閉じる
        try:
            for widget in QApplication.topLevelWidgets():
                if widget != self and isinstance(widget, (QDialog, QMainWindow)):
                    try:
                        widget.close()
                    except Exception:
                        pass
        except Exception:
            pass
        
        # 終了処理
        if self.scene and self.scene.template_preview:
            self.scene.template_preview.hide()

        # Clean up any active per-run calculation threads we spawned.
        try:
            for thr in list(getattr(self, '_active_calc_threads', []) or []):
                try:
                    thr.quit()
                except Exception:
                    pass
                try:
                    thr.wait(200)
                except Exception:
                    pass
        except Exception:
            pass
        
        event.accept()

    def zoom_in(self):
        """ ビューを 20% 拡大する """
        self.view_2d.scale(1.2, 1.2)

    def zoom_out(self):
        """ ビューを 20% 縮小する """
        self.view_2d.scale(1/1.2, 1/1.2)
        
    def reset_zoom(self):
        """ ビューの拡大率をデフォルト (75%) にリセットする """
        transform = QTransform()
        transform.scale(0.75, 0.75)
        self.view_2d.setTransform(transform)

    def fit_to_view(self):
        """ シーン上のすべてのアイテムがビューに収まるように調整する """
        if not self.scene.items():
            self.reset_zoom()
            return
            
        # 合計の表示矩形（目に見えるアイテムのみ）を計算
        visible_items_rect = QRectF()
        for item in self.scene.items():
            if item.isVisible() and not isinstance(item, TemplatePreviewItem):
                if visible_items_rect.isEmpty():
                    visible_items_rect = item.sceneBoundingRect()
                else:
                    visible_items_rect = visible_items_rect.united(item.sceneBoundingRect())

        if visible_items_rect.isEmpty():
            self.reset_zoom()
            return

        # 少し余白を持たせる（パディング）
        padding_factor = 1.10  # 10% の余裕
        cx = visible_items_rect.center().x()
        cy = visible_items_rect.center().y()
        w = visible_items_rect.width() * padding_factor
        h = visible_items_rect.height() * padding_factor
        padded = QRectF(cx - w / 2.0, cy - h / 2.0, w, h)

        # フィット時にマウス位置に依存するアンカーが原因でジャンプすることがあるため
        # 一時的にトランスフォームアンカーをビュー中心にしてから fitInView を呼ぶ
        try:
            old_ta = self.view_2d.transformationAnchor()
            old_ra = self.view_2d.resizeAnchor()
        except Exception:
            old_ta = old_ra = None

        try:
            self.view_2d.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
            self.view_2d.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
            self.view_2d.fitInView(padded, Qt.AspectRatioMode.KeepAspectRatio)
        finally:
            # 元のアンカーを復元
            try:
                if old_ta is not None:
                    self.view_2d.setTransformationAnchor(old_ta)
                if old_ra is not None:
                    self.view_2d.setResizeAnchor(old_ra)
            except Exception:
                pass

    def toggle_3d_edit_mode(self, checked):
        """「3D Drag」ボタンの状態に応じて編集モードを切り替える"""
        if checked:
            # 3D Editモードをオンにする時は、Measurementモードを無効化
            if self.measurement_mode:
                self.measurement_action.setChecked(False)
                self.toggle_measurement_mode(False)
        
        self.is_3d_edit_mode = checked
        if checked:
            self.statusBar().showMessage("3D Drag Mode: ON.")
        else:
            self.statusBar().showMessage("3D Drag Mode: OFF.")
        self.view_2d.setFocus()

    def _setup_3d_picker(self):
        self.plotter.picker = vtk.vtkCellPicker()
        self.plotter.picker.SetTolerance(0.025)

        # 新しいカスタムスタイル（原子移動用）のインスタンスを作成
        style = CustomInteractorStyle(self)
        
        # 調査の結果、'style' プロパティへの代入が正しい設定方法と判明
        self.plotter.interactor.SetInteractorStyle(style)
        self.plotter.interactor.Initialize()

    def _apply_chem_check_and_set_flags(self, mol, source_desc=None):
        """Central helper to apply chemical sanitization (or skip it) and set
        chem_check_tried / chem_check_failed flags consistently.

        When sanitization fails, a warning is shown and the Optimize 3D button
        is disabled. If the user setting 'skip_chemistry_checks' is True, no
        sanitization is attempted and both flags remain False.
        """
        try:
            self.chem_check_tried = False
            self.chem_check_failed = False
        except Exception:
            # Ensure attributes exist even if called very early
            self.chem_check_tried = False
            self.chem_check_failed = False

        if self.settings.get('skip_chemistry_checks', False):
            # User asked to skip chemistry checks entirely
            return

        try:
            Chem.SanitizeMol(mol)
            self.chem_check_tried = True
            self.chem_check_failed = False
        except Exception:
            # Mark that we tried sanitization and it failed
            self.chem_check_tried = True
            self.chem_check_failed = True
            try:
                desc = f" ({source_desc})" if source_desc else ''
                self.statusBar().showMessage(f"Molecule sanitization failed{desc}; file may be malformed.")
            except Exception:
                pass
            # Disable 3D optimization UI to prevent running on invalid molecules
            if hasattr(self, 'optimize_3d_button'):
                try:
                    self.optimize_3d_button.setEnabled(False)
                except Exception:
                    pass
        
    def _clear_xyz_flags(self, mol=None):
        """Clear XYZ-derived markers from a molecule (or current_mol) and
        reset UI flags accordingly.

        This is a best-effort cleanup to remove properties like
        _xyz_skip_checks and _xyz_atom_data that may have been attached when
        an XYZ file was previously loaded. After clearing molecule-level
        markers, the UI flag self.is_xyz_derived is set to False and the
        Optimize 3D button is re-evaluated (enabled unless chem_check_failed
        is True).
        """
        target = mol if mol is not None else getattr(self, 'current_mol', None)
        try:
            if target is not None:
                # Remove RDKit property if present
                try:
                    if hasattr(target, 'HasProp') and target.HasProp('_xyz_skip_checks'):
                        try:
                            target.ClearProp('_xyz_skip_checks')
                        except Exception:
                            try:
                                target.SetIntProp('_xyz_skip_checks', 0)
                            except Exception:
                                pass
                except Exception:
                    pass

                # Remove attribute-style markers if present
                try:
                    if hasattr(target, '_xyz_skip_checks'):
                        try:
                            delattr(target, '_xyz_skip_checks')
                        except Exception:
                            try:
                                del target._xyz_skip_checks
                            except Exception:
                                try:
                                    target._xyz_skip_checks = False
                                except Exception:
                                    pass
                except Exception:
                    pass

                try:
                    if hasattr(target, '_xyz_atom_data'):
                        try:
                            delattr(target, '_xyz_atom_data')
                        except Exception:
                            try:
                                del target._xyz_atom_data
                            except Exception:
                                try:
                                    target._xyz_atom_data = None
                                except Exception:
                                    pass
                except Exception:
                    pass

        except Exception:
            # best-effort only
            pass

        # Reset UI flags
        try:
            self.is_xyz_derived = False
        except Exception:
            pass

        # Enable Optimize 3D unless sanitization failed
        try:
            if hasattr(self, 'optimize_3d_button'):
                if getattr(self, 'chem_check_failed', False):
                    try:
                        self.optimize_3d_button.setEnabled(False)
                    except Exception:
                        pass
                else:
                    try:
                        self.optimize_3d_button.setEnabled(True)
                    except Exception:
                        pass
        except Exception:
            pass
            
    def load_mol_file_for_3d_viewing(self, file_path=None):
        """MOL/SDFファイルを3Dビューアーで開く"""
        if not self.check_unsaved_changes():
                return  # ユーザーがキャンセルした場合は何もしない
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open MOL/SDF File", "", 
                "MOL/SDF Files (*.mol *.sdf);;All Files (*)"
            )
            if not file_path:
                return
        
        try:
            
            # Determine extension early and handle .mol specially by reading the
            # raw block and running it through fix_mol_block before parsing.
            _, ext = os.path.splitext(file_path)
            ext = ext.lower() if ext else ''

            if ext == '.sdf':
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)

            elif ext == '.mol':
                # Read the file contents and attempt to fix malformed counts lines
                with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)

                # If parsing the fixed block fails, fall back to RDKit's file reader
                # as a last resort (keeps behavior conservative).
                if mol is None:
                    try:
                        mol = Chem.MolFromMolFile(file_path, removeHs=False)
                    except Exception:
                        mol = None

                if mol is None:
                    self.statusBar().showMessage(f"Failed to load molecule from {file_path}")
                    return

            else:
                # Default: let RDKit try to read the file (most common case)
                if file_path.lower().endswith('.sdf'):
                    suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                    mol = next(suppl, None)
                else:
                    mol = Chem.MolFromMolFile(file_path, removeHs=False)

                if mol is None:
                    self.statusBar().showMessage(f"Failed to load molecule from {file_path}")
                    return

            # 3D座標がない場合は2Dから3D変換（最適化なし）
            if mol.GetNumConformers() == 0:
                self.statusBar().showMessage("No 3D coordinates found. Converting to 3D...")
                try:
                    try:
                        AllChem.EmbedMolecule(mol)
                        # 最適化は実行しない
                        # 3D変換直後にUndoスタックに積む
                        self.current_mol = mol
                        self.push_undo_state()
                    except Exception as e_embed:
                        # If skipping chemistry checks, allow molecule to be displayed without 3D embedding
                        if self.settings.get('skip_chemistry_checks', False):
                            self.statusBar().showMessage("Warning: failed to generate 3D coordinates but skip_chemistry_checks is enabled; continuing.")
                            # Keep mol as-is (may lack conformer); downstream code checks for conformers
                        else:
                            raise
                except:
                    self.statusBar().showMessage("Failed to generate 3D coordinates")
                    return
            
            # Clear XYZ markers on the newly loaded MOL/SDF so Optimize 3D is
            # correctly enabled when appropriate.
            try:
                self._clear_xyz_flags(mol)
            except Exception:
                pass

            # 3Dビューアーに表示
            # Centralized chemical/sanitization handling
            # Ensure the skip_chemistry_checks setting is respected and flags are set
            self._apply_chem_check_and_set_flags(mol, source_desc='MOL/SDF')

            self.current_mol = mol
            self.draw_molecule_3d(mol)
            
            # カメラをリセット
            self.plotter.reset_camera()
            
            # UIを3Dビューアーモードに設定
            self._enter_3d_viewer_ui_mode()
            
            # メニューテキストと状態を更新
            self.update_atom_id_menu_text()
            self.update_atom_id_menu_state()
            
            self.statusBar().showMessage(f"Loaded {file_path} in 3D viewer")
            
            self.reset_undo_stack()
            self.has_unsaved_changes = False  # ファイル読込直後は未変更扱い
            self.current_file_path = file_path
            self.update_window_title()
            

        except Exception as e:
            self.statusBar().showMessage(f"Error loading MOL/SDF file: {e}")
    
    def load_command_line_file(self, file_path):
        """コマンドライン引数で指定されたファイルを開く"""
        if not file_path or not os.path.exists(file_path):
            return
        
        file_ext = file_path.lower().split('.')[-1]
        
        if file_ext in ['mol', 'sdf']:
            self.load_mol_file_for_3d_viewing(file_path)
        elif file_ext == 'xyz':
            self.load_xyz_for_3d_viewing(file_path)
        elif file_ext in ['pmeraw', 'pmeprj']:
            self.open_project_file(file_path=file_path)
        else:
            self.statusBar().showMessage(f"Unsupported file type: {file_ext}")
        
    def dragEnterEvent(self, event):
        """ウィンドウ全体で .pmeraw、.pmeprj、.mol、.sdf、.xyz ファイルのドラッグを受け入れる"""
        # Accept if any dragged local file has a supported extension
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            for url in urls:
                try:
                    if url.isLocalFile():
                        file_path = url.toLocalFile()
                        if file_path.lower().endswith(('.pmeraw', '.pmeprj', '.mol', '.sdf', '.xyz')):
                            event.acceptProposedAction()
                            return
                except Exception:
                    continue
        event.ignore()

    def dropEvent(self, event):
        """ファイルがウィンドウ上でドロップされたときに呼び出される"""
        urls = event.mimeData().urls()
        # Find the first local file from the dropped URLs
        file_path = None
        if urls:
            for url in urls:
                try:
                    if url.isLocalFile():
                        file_path = url.toLocalFile()
                        break
                except Exception:
                    continue

        if file_path:
            # ドロップ位置を取得
            drop_pos = event.position().toPoint()
            # 拡張子に応じて適切な読み込みメソッドを呼び出す
            if file_path.lower().endswith((".pmeraw", ".pmeprj")):
                self.open_project_file(file_path=file_path)
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            elif file_path.lower().endswith((".mol", ".sdf")):
                plotter_widget = self.splitter.widget(1)  # 3Dビューアーウィジェット
                plotter_rect = plotter_widget.geometry()
                if plotter_rect.contains(drop_pos):
                    self.load_mol_file_for_3d_viewing(file_path=file_path)
                else:
                    if hasattr(self, "load_mol_file"):
                        self.load_mol_file(file_path=file_path)
                    else:
                        self.statusBar().showMessage("MOL file import not implemented for 2D editor.")
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            elif file_path.lower().endswith(".xyz"):
                self.load_xyz_for_3d_viewing(file_path=file_path)
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            else:
                self.statusBar().showMessage(f"Unsupported file type: {file_path}")
                event.ignore()
        else:
            event.ignore()

    def _enable_3d_edit_actions(self, enabled=True):
        """3D編集機能のアクションを統一的に有効/無効化する"""
        actions = [
            'translation_action',
            'move_group_action',
            'alignplane_xy_action',
            'alignplane_xz_action',
            'alignplane_yz_action',
            'align_x_action',
            'align_y_action', 
            'align_z_action',
            'bond_length_action',
            'angle_action',
            'dihedral_action',
            'mirror_action',
            'planarize_action',
            'constrained_opt_action'
        ]
        
        # メニューとサブメニューも有効/無効化
        menus = [
            'align_menu'
        ]
        
        for action_name in actions:
            if hasattr(self, action_name):
                getattr(self, action_name).setEnabled(enabled)
        
        for menu_name in menus:
            if hasattr(self, menu_name):
                getattr(self, menu_name).setEnabled(enabled)

    def _enable_3d_features(self, enabled=True):
        """3D関連機能を統一的に有効/無効化する"""
        # 基本的な3D機能（3D SelectとEditは除外して常に有効にする）
        basic_3d_actions = [
            'optimize_3d_button',
            'export_button', 
            'analysis_action'
        ]
        
        for action_name in basic_3d_actions:
            if hasattr(self, action_name):
                # If enabling globally but chemical sanitization failed earlier, keep Optimize 3D disabled
                # Keep Optimize disabled when any of these conditions are true:
                # - we're globally disabling 3D features (enabled==False)
                # - the current molecule was created via the "skip chemistry checks" XYZ path
                # - a prior chemistry check was attempted and failed
                if action_name == 'optimize_3d_button':
                    try:
                        # If we're disabling all 3D features, ensure Optimize is disabled
                        if not enabled:
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # If the current molecule was marked as XYZ-derived (skip path), keep Optimize disabled
                        if getattr(self, 'is_xyz_derived', False):
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # If a chemistry check was tried and failed, keep Optimize disabled
                        if getattr(self, 'chem_check_tried', False) and getattr(self, 'chem_check_failed', False):
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # Otherwise enable/disable according to the requested global flag
                        getattr(self, action_name).setEnabled(bool(enabled))
                    except Exception:
                        pass
                else:
                    try:
                        getattr(self, action_name).setEnabled(enabled)
                    except Exception:
                        pass
        
        # 3D Selectボタンは常に有効にする
        if hasattr(self, 'measurement_action'):
            self.measurement_action.setEnabled(True)
        
        # 3D Dragボタンも常に有効にする
        if hasattr(self, 'edit_3d_action'):
            self.edit_3d_action.setEnabled(True)
        
        # 3D編集機能も含める
        if enabled:
            self._enable_3d_edit_actions(True)
        else:
            self._enable_3d_edit_actions(False)

    def _enter_3d_viewer_ui_mode(self):
        """3DビューアモードのUI状態に設定する"""
        self.is_2d_editable = False
        self.cleanup_button.setEnabled(False)
        self.convert_button.setEnabled(False)
        for action in self.tool_group.actions():
            action.setEnabled(False)
        if hasattr(self, 'other_atom_action'):
            self.other_atom_action.setEnabled(False)
        
        self.minimize_2d_panel()

        # 3D関連機能を統一的に有効化
        self._enable_3d_features(True)

    def restore_ui_for_editing(self):
        """Enables all 2D editing UI elements."""
        self.is_2d_editable = True
        self.restore_2d_panel()
        self.cleanup_button.setEnabled(True)
        self.convert_button.setEnabled(True)

        for action in self.tool_group.actions():
            action.setEnabled(True)
        
        if hasattr(self, 'other_atom_action'):
            self.other_atom_action.setEnabled(True)
            
        # 2Dモードに戻る時は3D編集機能を統一的に無効化
        self._enable_3d_edit_actions(False)

    def minimize_2d_panel(self):
        """2Dパネルを最小化（非表示に）する"""
        sizes = self.splitter.sizes()
        # すでに最小化されていなければ実行
        if sizes[0] > 0:
            total_width = sum(sizes)
            self.splitter.setSizes([0, total_width])

    def restore_2d_panel(self):
        """最小化された2Dパネルを元のサイズに戻す"""
        sizes = self.splitter.sizes()
        
        # sizesリストが空でないことを確認してからアクセスする
        if sizes and sizes[0] == 0:
            self.splitter.setSizes([600, 600])

    def set_panel_layout(self, left_percent, right_percent):
        """パネルレイアウトを指定した比率に設定する"""
        if left_percent + right_percent != 100:
            return
        
        total_width = self.splitter.width()
        if total_width <= 0:
            total_width = 1200  # デフォルト幅
        
        left_width = int(total_width * left_percent / 100)
        right_width = int(total_width * right_percent / 100)
        
        self.splitter.setSizes([left_width, right_width])
        
        # ユーザーにフィードバック表示
        self.statusBar().showMessage(
            f"Panel layout set to {left_percent}% : {right_percent}%", 
            2000
        )

    def toggle_2d_panel(self):
        """2Dパネルの表示/非表示を切り替える"""
        sizes = self.splitter.sizes()
        if not sizes:
            return
            
        if sizes[0] == 0:
            # 2Dパネルが非表示の場合は表示
            self.restore_2d_panel()
            self.statusBar().showMessage("2D panel restored", 1500)
        else:
            # 2Dパネルが表示されている場合は非表示
            self.minimize_2d_panel()
            self.statusBar().showMessage("2D panel minimized", 1500)

    def on_splitter_moved(self, pos, index):
        """スプリッターが移動された時のフィードバック表示"""
        sizes = self.splitter.sizes()
        if len(sizes) >= 2:
            total = sum(sizes)
            if total > 0:
                left_percent = round(sizes[0] * 100 / total)
                right_percent = round(sizes[1] * 100 / total)
                
                # 現在の比率をツールチップで表示
                if hasattr(self.splitter, 'handle'):
                    handle = self.splitter.handle(1)
                    if handle:
                        handle.setToolTip(f"2D: {left_percent}% | 3D: {right_percent}%")

    def open_template_dialog(self):
        """テンプレートダイアログを開く"""
        dialog = UserTemplateDialog(self, self)
        dialog.exec()
    
    def open_template_dialog_and_activate(self):
        """テンプレートダイアログを開き、テンプレートがメイン画面で使用できるようにする"""
        # 既存のダイアログがあるかチェック
        if hasattr(self, '_template_dialog') and self._template_dialog and not self._template_dialog.isHidden():
            # 既存のダイアログを前面に表示
            self._template_dialog.raise_()
            self._template_dialog.activateWindow()
            return
        
        # 新しいダイアログを作成
        self._template_dialog = UserTemplateDialog(self, self)
        self._template_dialog.show()  # モードレスで表示
        
        # ダイアログが閉じられた後、テンプレートが選択されていればアクティブ化
        def on_dialog_finished():
            if hasattr(self._template_dialog, 'selected_template') and self._template_dialog.selected_template:
                template_name = self._template_dialog.selected_template.get('name', 'user_template')
                mode_name = f"template_user_{template_name}"
                
                # Store template data for the scene to use
                self.scene.user_template_data = self._template_dialog.selected_template
                self.set_mode(mode_name)
                
                # Update status
                self.statusBar().showMessage(f"Template mode: {template_name}")
        
        self._template_dialog.finished.connect(on_dialog_finished)
    
    def save_2d_as_template(self):
        """現在の2D構造をテンプレートとして保存"""
        if not self.data.atoms:
            QMessageBox.warning(self, "Warning", "No structure to save as template.")
            return
        
        # Get template name
        name, ok = QInputDialog.getText(self, "Save Template", "Enter template name:")
        if not ok or not name.strip():
            return
        
        name = name.strip()
        
        try:
            # Template directory
            template_dir = os.path.join(self.settings_dir, 'user-templates')
            if not os.path.exists(template_dir):
                os.makedirs(template_dir)
            
            # Convert current structure to template format
            atoms_data = []
            bonds_data = []
            
            # Convert atoms
            for atom_id, atom_info in self.data.atoms.items():
                pos = atom_info['pos']
                atoms_data.append({
                    'id': atom_id,
                    'symbol': atom_info['symbol'],
                    'x': pos.x(),
                    'y': pos.y(),
                    'charge': atom_info.get('charge', 0),
                    'radical': atom_info.get('radical', 0)
                })
            
            # Convert bonds
            for (atom1_id, atom2_id), bond_info in self.data.bonds.items():
                bonds_data.append({
                    'atom1': atom1_id,
                    'atom2': atom2_id,
                    'order': bond_info['order'],
                    'stereo': bond_info.get('stereo', 0)
                })
            
            # Create template data
            template_data = {
                'format': "PME Template",
                'version': "1.0",
                'application': "MoleditPy",
                'application_version': VERSION,
                'name': name,
                'created': str(QDateTime.currentDateTime().toString()),
                'atoms': atoms_data,
                'bonds': bonds_data
            }
            
            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(template_dir, filename)
            
            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self, "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return
            
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)
            
            # Mark as saved (no unsaved changes for this operation)
            self.has_unsaved_changes = False
            self.update_window_title()
            
            QMessageBox.information(self, "Success", f"Template '{name}' saved successfully.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save template: {str(e)}")

    def setup_splitter_tooltip(self):
        """スプリッターハンドルの初期ツールチップを設定"""
        handle = self.splitter.handle(1)
        if handle:
            handle.setToolTip("Drag to resize panels | Ctrl+1/2/3 for presets | Ctrl+H to toggle 2D panel")
            # 初期サイズ比率も表示
            self.on_splitter_moved(0, 0)

            
    def apply_initial_settings(self):
        """UIの初期化が完了した後に、保存された設定を3Dビューに適用する"""
        
        try:
            self.update_cpk_colors_from_settings()
        except Exception:
            pass

        if self.plotter and self.plotter.renderer:
            bg_color = self.settings.get('background_color', '#919191')
            self.plotter.set_background(bg_color)
            self.apply_3d_settings()
        
        try:
            if hasattr(self, 'scene') and self.scene:
                for it in list(self.scene.items()):
                    if hasattr(it, 'update_style'):
                        it.update_style()
                self.scene.update()
                for v in list(self.scene.views()):
                    v.viewport().update()
        except Exception:
            pass

    def update_cpk_colors_from_settings(self):
        """Update global CPK_COLORS and CPK_COLORS_PV from saved settings overrides.

        This modifies the in-memory CPK_COLORS mapping (not persisted until settings are saved).
        Only keys present in self.settings['cpk_colors'] are changed; other elements keep the defaults.
        """
        try:
            overrides = self.settings.get('cpk_colors', {}) or {}
            # Start from a clean copy of the defaults
            global CPK_COLORS, CPK_COLORS_PV
            CPK_COLORS = {k: QColor(v) if not isinstance(v, QColor) else QColor(v) for k, v in DEFAULT_CPK_COLORS.items()}
            for k, hexv in overrides.items():
                if isinstance(hexv, str) and hexv:
                    CPK_COLORS[k] = QColor(hexv)
            # Rebuild PV version
            CPK_COLORS_PV = {k: [c.redF(), c.greenF(), c.blueF()] for k, c in CPK_COLORS.items()}
        except Exception as e:
            print(f"Failed to update CPK colors from settings: {e}")


    def apply_3d_settings(self, redraw=True):
        # Projection mode
        proj_mode = self.settings.get('projection_mode', 'Perspective')
        if hasattr(self.plotter, 'renderer') and hasattr(self.plotter.renderer, 'GetActiveCamera'):
            cam = self.plotter.renderer.GetActiveCamera()
            if cam:
                if proj_mode == 'Orthographic':
                    cam.SetParallelProjection(True)
                else:
                    cam.SetParallelProjection(False)
        """3Dビューの視覚設定を適用する"""
        if not hasattr(self, 'plotter'):
            return  
        
        # レンダラーのレイヤー設定を有効化（テキストオーバーレイ用）
        renderer = self.plotter.renderer
        if renderer and hasattr(renderer, 'SetNumberOfLayers'):
            try:
                renderer.SetNumberOfLayers(2)  # レイヤー0:3Dオブジェクト、レイヤー1:2Dオーバーレイ
            except:
                pass  # PyVistaのバージョンによってはサポートされていない場合がある  

        # --- 3D軸ウィジェットの設定 ---
        show_axes = self.settings.get('show_3d_axes', True) 

        # ウィジェットがまだ作成されていない場合は作成する
        if self.axes_widget is None and hasattr(self.plotter, 'interactor'):
            axes = vtk.vtkAxesActor()
            self.axes_widget = vtk.vtkOrientationMarkerWidget()
            self.axes_widget.SetOrientationMarker(axes)
            self.axes_widget.SetInteractor(self.plotter.interactor)
            # 左下隅に設定 (幅・高さ20%)
            self.axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)

        # 設定に応じてウィジェットを有効化/無効化
        if self.axes_widget:
            if show_axes:
                self.axes_widget.On()
                self.axes_widget.SetInteractive(False)  
            else:
                self.axes_widget.Off()  

        if redraw:
            self.draw_molecule_3d(self.current_mol)

        # 設定変更時にカメラ位置をリセットしない（初回のみリセット）
        if not getattr(self, '_camera_initialized', False):
            try:
                self.plotter.reset_camera()
            except Exception:
                pass
            self._camera_initialized = True
        
        # 強制的にプロッターを更新
        try:
            self.plotter.render()
            if hasattr(self.plotter, 'update'):
                self.plotter.update()
        except Exception:
            pass



    def open_settings_dialog(self):
        dialog = SettingsDialog(self.settings, self)
        # accept()メソッドで設定の適用と3Dビューの更新を行うため、ここでは不要
        dialog.exec()


    def reset_all_settings_menu(self):
        # Expose the same functionality as SettingsDialog.reset_all_settings
        dlg = QMessageBox(self)
        dlg.setIcon(QMessageBox.Icon.Warning)
        dlg.setWindowTitle("Reset All Settings")
        dlg.setText("Are you sure you want to reset all settings to defaults?")
        dlg.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        res = dlg.exec()
        if res == QMessageBox.StandardButton.Yes:
            try:
                # Remove settings file and reload defaults
                if os.path.exists(self.settings_file):
                    os.remove(self.settings_file)
                self.load_settings()
                # Do not write to disk immediately; mark dirty so settings will be saved on exit
                try:
                    self.settings_dirty = True
                except Exception:
                    pass
                # If ColorSettingsDialog is open, refresh its UI to reflect the reset
                try:
                    for w in QApplication.topLevelWidgets():
                        try:
                            if isinstance(w, ColorSettingsDialog):
                                try:
                                    w.refresh_ui()
                                except Exception:
                                    pass
                        except Exception:
                            pass
                except Exception:
                    pass
                # Ensure global CPK mapping is rebuilt from defaults and UI is updated
                try:
                    self.update_cpk_colors_from_settings()
                except Exception:
                    pass
                # Refresh UI/menu state for conversion and optimization
                try:
                    # update optimization method
                    self.optimization_method = self.settings.get('optimization_method', 'MMFF_RDKIT')
                    if hasattr(self, 'opt3d_actions') and self.optimization_method:
                        key = (self.optimization_method or '').upper()
                        if key in self.opt3d_actions:
                            # uncheck all then check the saved one
                            for act in self.opt3d_actions.values():
                                act.setChecked(False)
                            try:
                                self.opt3d_actions[key].setChecked(True)
                            except Exception:
                                pass

                    # update conversion mode
                    conv_mode = self.settings.get('3d_conversion_mode', 'fallback')
                    if hasattr(self, 'conv_actions') and conv_mode in self.conv_actions:
                        try:
                            for act in self.conv_actions.values():
                                act.setChecked(False)
                            self.conv_actions[conv_mode].setChecked(True)
                        except Exception:
                            pass

                    # 3Dビューの設定を適用
                    self.apply_3d_settings()
                    # 現在の分子を再描画（設定変更を反映）
                    if hasattr(self, 'current_mol') and self.current_mol:
                        self.draw_molecule_3d(self.current_mol)
                    
                    QMessageBox.information(self, "Reset Complete", "All settings have been reset to defaults.")
                    
                except Exception:
                    pass
                # Update 2D scene styling to reflect default CPK colors
                try:
                    if hasattr(self, 'scene') and self.scene:
                        for it in list(self.scene.items()):
                            try:
                                if hasattr(it, 'update_style'):
                                    it.update_style()
                            except Exception:
                                pass
                        try:
                            # Force a full scene update and viewport repaint for all views
                            self.scene.update()
                            for v in list(self.scene.views()):
                                try:
                                    v.viewport().update()
                                except Exception:
                                    pass
                        except Exception:
                            pass
                except Exception:
                    pass
                # Also refresh any open SettingsDialog instances so their UI matches
                try:
                    for w in QApplication.topLevelWidgets():
                        try:
                            if isinstance(w, SettingsDialog):
                                try:
                                    w.update_ui_from_settings(self.settings)
                                except Exception:
                                    pass
                        except Exception:
                            pass
                except Exception:
                    pass
            except Exception as e:
                QMessageBox.warning(self, "Reset Failed", f"Could not reset settings: {e}")
            

    def load_settings(self):
        default_settings = {
            'background_color': '#919191',
            'projection_mode': 'Perspective',
            'lighting_enabled': True,
            'specular': 0.2,
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
            'stick_atom_radius': 0.15,
            'stick_bond_radius': 0.15,
            'stick_resolution': 16,
            # Multiple bond offset parameters (per-model)
            'ball_stick_double_bond_offset_factor': 2.0,
            'ball_stick_triple_bond_offset_factor': 2.0,
            'ball_stick_double_bond_radius_factor': 0.8,
            'ball_stick_triple_bond_radius_factor': 0.75,
            'wireframe_double_bond_offset_factor': 3.0,
            'wireframe_triple_bond_offset_factor': 3.0,
            'wireframe_double_bond_radius_factor': 1.0,
            'wireframe_triple_bond_radius_factor': 0.75,
            'stick_double_bond_offset_factor': 1.5,
            'stick_triple_bond_offset_factor': 1.0,
            'stick_double_bond_radius_factor': 0.60,
            'stick_triple_bond_radius_factor': 0.40,
            # Ensure conversion/optimization defaults are present
            # If True, attempts to be permissive when RDKit raises chemical/sanitization errors
            # during file import (useful for viewing malformed XYZ/MOL files).
            'skip_chemistry_checks': False,
            '3d_conversion_mode': 'fallback',
            'optimization_method': 'MMFF_RDKIT',
            # Color overrides
            'ball_stick_bond_color': '#7F7F7F',
            'cpk_colors': {},  # symbol->hex overrides
        }

        try:
            if os.path.exists(self.settings_file):
                with open(self.settings_file, 'r') as f:
                    loaded_settings = json.load(f)
                # Ensure any missing default keys are inserted and persisted.
                changed = False
                for key, value in default_settings.items():
                    if key not in loaded_settings:
                        loaded_settings[key] = value
                        changed = True

                self.settings = loaded_settings

                # Migration: if older global multi-bond keys exist, copy them to per-model keys
                legacy_keys = ['double_bond_offset_factor', 'triple_bond_offset_factor', 'double_bond_radius_factor', 'triple_bond_radius_factor']
                migrated = False
                # If legacy keys exist, propagate to per-model keys when per-model keys missing
                if any(k in self.settings for k in legacy_keys):
                    # For each per-model key, if missing, set from legacy fallback
                    def copy_if_missing(new_key, legacy_key, default_val):
                        nonlocal migrated
                        if new_key not in self.settings:
                            if legacy_key in self.settings:
                                self.settings[new_key] = self.settings[legacy_key]
                                migrated = True
                            else:
                                self.settings[new_key] = default_val
                                migrated = True

                    per_model_map = [
                        ('ball_stick_double_bond_offset_factor', 'double_bond_offset_factor', 2.0),
                        ('ball_stick_triple_bond_offset_factor', 'triple_bond_offset_factor', 2.0),
                        ('ball_stick_double_bond_radius_factor', 'double_bond_radius_factor', 0.8),
                        ('ball_stick_triple_bond_radius_factor', 'triple_bond_radius_factor', 0.75),
                        ('wireframe_double_bond_offset_factor', 'double_bond_offset_factor', 3.0),
                        ('wireframe_triple_bond_offset_factor', 'triple_bond_offset_factor', 3.0),
                        ('wireframe_double_bond_radius_factor', 'double_bond_radius_factor', 1.0),
                        ('wireframe_triple_bond_radius_factor', 'triple_bond_radius_factor', 0.75),
                        ('stick_double_bond_offset_factor', 'double_bond_offset_factor', 1.5),
                        ('stick_triple_bond_offset_factor', 'triple_bond_offset_factor', 1.0),
                        ('stick_double_bond_radius_factor', 'double_bond_radius_factor', 0.60),
                        ('stick_triple_bond_radius_factor', 'triple_bond_radius_factor', 0.40),
                    ]
                    for new_k, legacy_k, default_v in per_model_map:
                        copy_if_missing(new_k, legacy_k, default_v)

                    # Optionally remove legacy keys to avoid confusion (keep them for now but mark dirty)
                    if migrated:
                        changed = True

                # If we added any defaults (e.g. skip_chemistry_checks) or migrated keys, write them back so
                # the configuration file reflects the effective defaults without requiring
                # the user to edit the file manually.
                if changed:
                    # Don't write immediately; mark dirty and let closeEvent persist
                    try:
                        self.settings_dirty = True
                    except Exception:
                        pass
            
            else:
                # No settings file - use defaults. Mark dirty so defaults will be written on exit.
                self.settings = default_settings
                try:
                    self.settings_dirty = True
                except Exception:
                    pass
        
        except Exception:
            self.settings = default_settings

    def save_settings(self):
        try:
            if not os.path.exists(self.settings_dir):
                os.makedirs(self.settings_dir)
            with open(self.settings_file, 'w') as f:
                json.dump(self.settings, f, indent=4)
        except Exception as e:
            print(f"Error saving settings: {e}")

    def toggle_measurement_mode(self, checked):
        """測定モードのオン/オフを切り替える"""
        if checked:
            # 測定モードをオンにする時は、3D Dragモードを無効化
            if self.is_3d_edit_mode:
                self.edit_3d_action.setChecked(False)
                self.toggle_3d_edit_mode(False)
            
            # アクティブな3D編集ダイアログを閉じる
            self.close_all_3d_edit_dialogs()
        
        self.measurement_mode = checked
        
        if not checked:
            self.clear_measurement_selection()
        
        # ボタンのテキストとステータスメッセージを更新
        if checked:
            self.statusBar().showMessage("Measurement mode enabled. Click atoms to measure distances/angles/dihedrals.")
        else:
            self.statusBar().showMessage("Measurement mode disabled.")
    
    def close_all_3d_edit_dialogs(self):
        """すべてのアクティブな3D編集ダイアログを閉じる"""
        dialogs_to_close = self.active_3d_dialogs.copy()
        for dialog in dialogs_to_close:
            try:
                dialog.close()
            except:
                pass
        self.active_3d_dialogs.clear()

    def handle_measurement_atom_selection(self, atom_idx):
        """測定用の原子選択を処理する"""
        # 既に選択されている原子の場合は除外
        if atom_idx in self.selected_atoms_for_measurement:
            return
        
        self.selected_atoms_for_measurement.append(atom_idx)
        
        '''
        # 4つ以上選択された場合はクリア
        if len(self.selected_atoms_for_measurement) > 4:
            self.clear_measurement_selection()
            self.selected_atoms_for_measurement.append(atom_idx)
        '''
        
        # 原子にラベルを追加
        self.add_measurement_label(atom_idx, len(self.selected_atoms_for_measurement))
        
        # 測定値を計算して表示
        self.calculate_and_display_measurements()

    def add_measurement_label(self, atom_idx, label_number):
        """原子に数字ラベルを追加する"""
        if not self.current_mol or atom_idx >= self.current_mol.GetNumAtoms():
            return
        
        # 測定ラベルリストを更新
        self.measurement_labels.append((atom_idx, str(label_number)))
        
        # 3Dビューの測定ラベルを再描画
        self.update_measurement_labels_display()
        
        # 2Dビューの測定ラベルも更新
        self.update_2d_measurement_labels()

    def update_measurement_labels_display(self):
        """測定ラベルを3D表示に描画する（原子中心配置）"""
        try:
            # 既存の測定ラベルを削除
            self.plotter.remove_actor('measurement_labels')
        except:
            pass
        
        if not self.measurement_labels or not self.current_mol:
            return
        
        # ラベル位置とテキストを準備
        pts, labels = [], []
        for atom_idx, label_text in self.measurement_labels:
            if atom_idx < len(self.atom_positions_3d):
                coord = self.atom_positions_3d[atom_idx].copy()
                # オフセットを削除して原子中心に配置
                pts.append(coord)
                labels.append(label_text)
        
        if pts and labels:
            # PyVistaのpoint_labelsを使用（赤色固定）
            self.plotter.add_point_labels(
                np.array(pts), 
                labels, 
                font_size=16,
                point_size=0,
                text_color='red',  # 測定時は常に赤色
                name='measurement_labels',
                always_visible=True,
                tolerance=0.01,
                show_points=False
            )

    def clear_measurement_selection(self):
        """測定選択をクリアする"""
        self.selected_atoms_for_measurement.clear()
        
        # 3Dビューのラベルを削除
        self.measurement_labels.clear()
        try:
            self.plotter.remove_actor('measurement_labels')
        except:
            pass
        
        # 2Dビューの測定ラベルも削除
        self.clear_2d_measurement_labels()
        
        # 測定結果のテキストを削除
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
                self.measurement_text_actor = None
            except:
                pass
        
        self.plotter.render()

    def update_2d_measurement_labels(self):
        """2Dビューで測定ラベルを更新表示する"""
        # 既存の2D測定ラベルを削除
        self.clear_2d_measurement_labels()
        
        # 現在の分子から原子-AtomItemマッピングを作成
        if not self.current_mol or not hasattr(self, 'data') or not self.data.atoms:
            return
            
        # RDKit原子インデックスから2D AtomItemへのマッピングを作成
        atom_idx_to_item = {}
        
        # シーンからAtomItemを取得してマッピング
        if hasattr(self, 'scene'):
            for item in self.scene.items():
                if hasattr(item, 'atom_id') and hasattr(item, 'symbol'):  # AtomItemかチェック
                    # 原子IDから対応するRDKit原子インデックスを見つける
                    rdkit_idx = self.find_rdkit_atom_index(item)
                    if rdkit_idx is not None:
                        atom_idx_to_item[rdkit_idx] = item
        
        # 測定ラベルを2Dビューに追加
        if not hasattr(self, 'measurement_label_items_2d'):
            self.measurement_label_items_2d = []
            
        for atom_idx, label_text in self.measurement_labels:
            if atom_idx in atom_idx_to_item:
                atom_item = atom_idx_to_item[atom_idx]
                self.add_2d_measurement_label(atom_item, label_text)

    def add_2d_measurement_label(self, atom_item, label_text):
        """特定のAtomItemに測定ラベルを追加する"""
        # ラベルアイテムを作成
        label_item = QGraphicsTextItem(label_text)
        label_item.setDefaultTextColor(QColor(255, 0, 0))  # 赤色
        label_item.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        
        # Z値を設定して最前面に表示（原子ラベルより上）
        label_item.setZValue(2000)  # より高い値で確実に最前面に配置
        
        # 原子の右上により近く配置
        atom_pos = atom_item.pos()
        atom_rect = atom_item.boundingRect()
        label_pos = QPointF(
            atom_pos.x() + atom_rect.width() / 4 + 2,
            atom_pos.y() - atom_rect.height() / 4 - 8
        )
        label_item.setPos(label_pos)
        
        # シーンに追加
        self.scene.addItem(label_item)
        self.measurement_label_items_2d.append(label_item)

    def clear_2d_measurement_labels(self):
        """2Dビューの測定ラベルを全て削除する"""
        if hasattr(self, 'measurement_label_items_2d'):
            for label_item in self.measurement_label_items_2d:
                try:
                    # Avoid touching partially-deleted wrappers
                    if sip_isdeleted_safe(label_item):
                        continue
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except Exception:
                        # Scene access or removal failed; skip
                        continue
                except Exception:
                    # If sip check itself fails, fall back to best-effort removal
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except Exception:
                        continue
            self.measurement_label_items_2d.clear()

    def find_rdkit_atom_index(self, atom_item):
        """AtomItemから対応するRDKit原子インデックスを見つける"""
        if not self.current_mol or not atom_item:
            return None
        
        # マッピング辞書を使用（最も確実）
        if hasattr(self, 'atom_id_to_rdkit_idx_map') and atom_item.atom_id in self.atom_id_to_rdkit_idx_map:
            return self.atom_id_to_rdkit_idx_map[atom_item.atom_id]
        
        # マッピングが存在しない場合はNone（外部ファイル読み込み時など）
        return None

    def calculate_and_display_measurements(self):
        """選択された原子に基づいて測定値を計算し表示する"""
        num_selected = len(self.selected_atoms_for_measurement)
        if num_selected < 2:
            return
        
        measurement_text = []
        
        if num_selected >= 2:
            # 距離の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            distance = self.calculate_distance(atom1_idx, atom2_idx)
            measurement_text.append(f"Distance 1-2: {distance:.3f} Å")
        
        if num_selected >= 3:
            # 角度の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1] 
            atom3_idx = self.selected_atoms_for_measurement[2]
            angle = self.calculate_angle(atom1_idx, atom2_idx, atom3_idx)
            measurement_text.append(f"Angle 1-2-3: {angle:.2f}°")
        
        if num_selected >= 4:
            # 二面角の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            atom3_idx = self.selected_atoms_for_measurement[2]
            atom4_idx = self.selected_atoms_for_measurement[3]
            dihedral = self.calculate_dihedral(atom1_idx, atom2_idx, atom3_idx, atom4_idx)
            measurement_text.append(f"Dihedral 1-2-3-4: {dihedral:.2f}°")
        
        # 測定結果を3D画面の右上に表示
        self.display_measurement_text(measurement_text)

    def calculate_distance(self, atom1_idx, atom2_idx):
        """2原子間の距離を計算する"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])
        return np.linalg.norm(pos2 - pos1)

    def calculate_angle(self, atom1_idx, atom2_idx, atom3_idx):
        """3原子の角度を計算する（中央が頂点）"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])  # 頂点
        pos3 = np.array(self.atom_positions_3d[atom3_idx])
        
        # ベクトルを計算
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        # 角度を計算（ラジアンから度に変換）
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        # 数値誤差による範囲外の値をクリップ
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)

    def calculate_dihedral(self, atom1_idx, atom2_idx, atom3_idx, atom4_idx):
        """4原子の二面角を計算する（正しい公式を使用）"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])
        pos3 = np.array(self.atom_positions_3d[atom3_idx])
        pos4 = np.array(self.atom_positions_3d[atom4_idx])
        
        # Vectors between consecutive atoms
        v1 = pos2 - pos1  # 1->2
        v2 = pos3 - pos2  # 2->3 (central bond)
        v3 = pos4 - pos3  # 3->4
        
        # Normalize the central bond vector
        v2_norm = v2 / np.linalg.norm(v2)
        
        # Calculate plane normal vectors
        n1 = np.cross(v1, v2)  # Normal to plane 1-2-3
        n2 = np.cross(v2, v3)  # Normal to plane 2-3-4
        
        # Normalize the normal vectors
        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)
        
        if n1_norm == 0 or n2_norm == 0:
            return 0.0  # Atoms are collinear
        
        n1 = n1 / n1_norm
        n2 = n2 / n2_norm
        
        # Calculate the cosine of the dihedral angle
        cos_angle = np.dot(n1, n2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        # Calculate the sine for proper sign determination
        sin_angle = np.dot(np.cross(n1, n2), v2_norm)
        
        # Calculate the dihedral angle with correct sign
        angle_rad = np.arctan2(sin_angle, cos_angle)
        return np.degrees(angle_rad)

    def display_measurement_text(self, measurement_lines):
        """測定結果のテキストを3D画面の左上に表示する（小さな等幅フォント）"""
        # 既存のテキストを削除
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
            except:
                pass
        
        if not measurement_lines:
            self.measurement_text_actor = None
            return
        
        # テキストを結合
        text = '\n'.join(measurement_lines)
        
        # 背景色から適切なテキスト色を決定
        try:
            bg_color_hex = self.settings.get('background_color', '#919191')
            bg_qcolor = QColor(bg_color_hex)
            if bg_qcolor.isValid():
                luminance = bg_qcolor.toHsl().lightness()
                text_color = 'black' if luminance > 128 else 'white'
            else:
                text_color = 'white'
        except:
            text_color = 'white'
        
        # 左上に表示（小さな等幅フォント）
        self.measurement_text_actor = self.plotter.add_text(
            text,
            position='upper_left',
            font_size=10,  # より小さく
            color=text_color,  # 背景に合わせた色
            font='courier',  # 等幅フォント
            name='measurement_display'
        )
        
        self.plotter.render()

    # --- 3D Drag functionality ---
    
    def toggle_atom_selection_3d(self, atom_idx):
        """3Dビューで原子の選択状態をトグルする"""
        if atom_idx in self.selected_atoms_3d:
            self.selected_atoms_3d.remove(atom_idx)
        else:
            self.selected_atoms_3d.add(atom_idx)
        
        # 選択状態のビジュアルフィードバックを更新
        self.update_3d_selection_display()
    
    def clear_3d_selection(self):
        """3Dビューでの原子選択をクリア"""
        self.selected_atoms_3d.clear()
        self.update_3d_selection_display()
    
    def update_3d_selection_display(self):
        """3Dビューでの選択原子のハイライト表示を更新"""
        try:
            # 既存の選択ハイライトを削除
            self.plotter.remove_actor('selection_highlight')
        except:
            pass
        
        if not self.selected_atoms_3d or not self.current_mol:
            self.plotter.render()
            return
        
        # 選択された原子のインデックスリストを作成
        selected_indices = list(self.selected_atoms_3d)
        
        # 選択された原子の位置を取得
        selected_positions = self.atom_positions_3d[selected_indices]
        
        # 原子の半径を少し大きくしてハイライト表示
        selected_radii = np.array([VDW_RADII.get(
            self.current_mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3 
            for i in selected_indices])
        
        # ハイライト用のデータセットを作成
        highlight_source = pv.PolyData(selected_positions)
        highlight_source['radii'] = selected_radii
        
        # 黄色の半透明球でハイライト
        highlight_glyphs = highlight_source.glyph(
            scale='radii', 
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16), 
            orient=False
        )
        
        self.plotter.add_mesh(
            highlight_glyphs, 
            color='yellow', 
            opacity=0.3, 
            name='selection_highlight'
        )
        
        self.plotter.render()
    
    '''
    def planarize_selection(self, plane):
        """選択された原子群を指定された平面にPlanarizeする"""
        if not self.selected_atoms_3d or not self.current_mol:
            self.statusBar().showMessage("No atoms selected for align.")
            return

        if len(self.selected_atoms_3d) < 3:
            self.statusBar().showMessage("Please select at least 3 atoms for align.")
            return

        try:
            # 選択された原子の位置を取得
            selected_indices = list(self.selected_atoms_3d)
            selected_positions = self.atom_positions_3d[selected_indices].copy()

            # 重心を計算
            centroid = np.mean(selected_positions, axis=0)

            # 重心を原点に移動
            centered_positions = selected_positions - centroid

            # SVDで法線を取得
            u, s, vh = np.linalg.svd(centered_positions, full_matrices=False)
            normal = vh[-1]
            # 法線を正規化
            norm = np.linalg.norm(normal)
            if norm == 0:
                self.statusBar().showMessage("Cannot determine fit plane (degenerate positions).")
                return
            normal = normal / norm

            # 各点を重心を通る平面へ直交射影する
            projections = centered_positions - np.outer(np.dot(centered_positions, normal), normal)
            new_positions = projections + centroid
            plane_name = "best-fit"
        except Exception as e:
            self.statusBar().showMessage(f"Error computing fit plane: {e}")
            return
            
            # 分子の座標を更新
            conf = self.current_mol.GetConformer()
            for i, new_pos in zip(selected_indices, new_positions):
                conf.SetAtomPosition(i, new_pos.tolist())
                self.atom_positions_3d[i] = new_pos
            
            # 3Dビューを更新
            self.draw_molecule_3d(self.current_mol)
            
            # 選択状態を維持
            temp_selection = self.selected_atoms_3d.copy()
            self.selected_atoms_3d = temp_selection
            self.update_3d_selection_display()

            self.statusBar().showMessage(f"Planarized {len(selected_indices)} atoms to {plane_name} plane.")
            self.push_undo_state()
            
        except Exception as e:
            self.statusBar().showMessage(f"Error during planarization: {e}")
        '''
    
    def open_translation_dialog(self):
        """平行移動ダイアログを開く"""
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = TranslationDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Translation applied."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
    def open_move_group_dialog(self):
        """Move Groupダイアログを開く"""
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = MoveGroupDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Group transformation applied."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    
    def open_align_plane_dialog(self, plane):
        """alignダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AlignPlaneDialog(self.current_mol, self, plane, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage(f"Atoms alignd to {plane.upper()} plane."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
        
    def open_planarize_dialog(self, plane=None):
        """選択原子群を最適平面へ投影するダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = PlanarizeDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Selection planarized to best-fit plane."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    
    def open_alignment_dialog(self, axis):
        """アライメントダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AlignmentDialog(self.current_mol, self, axis, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage(f"Atoms aligned to {axis.upper()}-axis."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
    def open_bond_length_dialog(self):
        """結合長変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = BondLengthDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Bond length adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
    def open_angle_dialog(self):
        """角度変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AngleDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Angle adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
    def open_dihedral_dialog(self):
        """二面角変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = DihedralDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Dihedral angle adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
    def open_mirror_dialog(self):
        """ミラー機能ダイアログを開く"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = MirrorDialog(self.current_mol, self)
        dialog.exec()  # モーダルダイアログとして表示

    def open_constrained_optimization_dialog(self):
        """制約付き最適化ダイアログを開く"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = ConstrainedOptimizationDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # モードレス表示
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    
    def remove_dialog_from_list(self, dialog):
        """ダイアログをアクティブリストから削除"""
        if dialog in self.active_3d_dialogs:
            self.active_3d_dialogs.remove(dialog)
