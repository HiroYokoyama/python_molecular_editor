#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_window_main_init.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowMainInit
"""


import math
import os
import json 


# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QSplitter, QToolBar, QSizePolicy, QLabel, QToolButton, QMenu, QMessageBox
)

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QKeySequence, 
    QPixmap, QIcon, QShortcut, QDesktopServices
)


from PyQt6.QtCore import (
    Qt, QPointF, QRectF, QLineF, QUrl, QTimer
)
import platform
import subprocess
try:
    import winreg
except Exception:
    winreg = None


def detect_system_dark_mode():
    """Return True if the OS prefers dark app theme, False if light, or None if unknown.

    This is a best-effort, cross-platform check supporting Windows (registry),
    macOS (defaults read), and GNOME/GTK-based Linux (gsettings). Return
    None if no reliable information is available.
    """
    # Delegate detailed OS detection to `detect_system_theme` and map
    # 'dark' -> True, 'light' -> False. This avoids duplicating the
    # registry and subprocess calls in two places.
    theme = detect_system_theme()
    if theme == 'dark':
        return True
    if theme == 'light':
        return False
    return None

def detect_system_theme():
    """OSの優先テーマ設定を 'dark', 'light', または None として返す。

    This is a best-effort, cross-platform check.
    """
    try:
        # Windows: AppsUseLightTheme (0 = dark, 1 = light)
        if platform.system() == 'Windows' and winreg is not None:
            try:
                with winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                                    r'Software\Microsoft\Windows\CurrentVersion\Themes\Personalize') as k:
                    val, _ = winreg.QueryValueEx(k, 'AppsUseLightTheme')
                    return 'dark' if int(val) == 0 else 'light'
            except Exception:
                pass

        # macOS: 'defaults read -g AppleInterfaceStyle'
        if platform.system() == 'Darwin':
            return 'light'
            '''
            try:
                # 'defaults read ...' が成功すればダークモード
                p = subprocess.run(
                    ['defaults', 'read', '-g', 'AppleInterfaceStyle'],
                    capture_output=True, text=True, check=True, encoding='utf-8'
                )
                if p.stdout.strip().lower() == 'dark':
                    return 'dark'
                
            except subprocess.CalledProcessError:
                # コマンド失敗 (キーが存在しない) = ライトモード
                return 'light'
            except Exception:
                # その他のエラー
                pass
            '''

        # Linux / GNOME: try color-scheme gsetting; fallback to gtk-theme detection
        if platform.system() == 'Linux':
            try:
                p = subprocess.run(['gsettings', 'get', 'org.gnome.desktop.interface', 'color-scheme'], capture_output=True, text=True)
                if p.returncode == 0:
                    out = p.stdout.strip().strip("'\n ")
                    if 'dark' in out.lower():
                        return 'dark'
                    if 'light' in out.lower():
                        return 'light'
            except Exception:
                pass

            try:
                p = subprocess.run(['gsettings', 'get', 'org.gnome.desktop.interface', 'gtk-theme'], capture_output=True, text=True)
                if p.returncode == 0 and '-dark' in p.stdout.lower():
                    return 'dark'
            except Exception:
                pass
    except Exception:
        pass
    return None


# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None
    
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
    from .constants import NUM_DASHES, VERSION
    from .molecular_data import MolecularData
    from .molecule_scene import MoleculeScene
    from .zoomable_view import ZoomableView
    from .color_settings_dialog import ColorSettingsDialog
    from .settings_dialog import SettingsDialog
    from .custom_qt_interactor import CustomQtInteractor
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import NUM_DASHES, VERSION
    from modules.molecular_data import MolecularData
    from modules.molecule_scene import MoleculeScene
    from modules.zoomable_view import ZoomableView
    from modules.color_settings_dialog import ColorSettingsDialog
    from modules.settings_dialog import SettingsDialog
    from modules.custom_qt_interactor import CustomQtInteractor

    

# --- クラス定義 ---
class MainWindowMainInit(object):
    """ main_window.py から分離された機能クラス """

    # __init__ は main_window.py からコピーされます


    def __init__(self, initial_file=None):
        # This helper is not used as a mixin in this project; initialization
        # happens on the `MainWindow` instance. Avoid calling super() here
        # because we initialize the `QMainWindow` base class in
        # `MainWindow.__init__` directly.
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
        
        # 3D編集用の原子選択状態 (3Dビューで選択された原子のインデックス)
        self.selected_atoms_3d = set()
        
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
        # Use module-level detector `detect_system_dark_mode()` so tests and other
        # modules can reuse the logic.


        def _icon_foreground_color():
            """Return a QColor for icon foreground.

            NOTE: choose icon foreground to contrast against the background
            (i.e., white on dark backgrounds, black on light backgrounds). This
            matches common conventions. Priority: explicit setting in
            'icon_foreground' -> OS theme preference -> configured 3D
            theme preference -> configured 3D background -> application palette.
            """
            try:
                fg_hex = self.settings.get('icon_foreground')
                if fg_hex:
                    c = QColor(fg_hex)
                    if c.isValid():
                        return c
            except Exception:
                pass

            # 1) Prefer the system/OS dark-mode preference if available.
            try:
                os_pref = detect_system_dark_mode()
                # Standard mapping: dark -> white, light -> black
                if os_pref is not None:
                    return QColor('#FFFFFF') if os_pref else QColor('#000000')
            except Exception:
                pass

            try:
                # Keep background_color as a fallback: if system preference isn't
                # available we'll use the configured 3D view background from settings.
                bg_hex = self.settings.get('background_color')
                if bg_hex:
                    bg = QColor(bg_hex)
                    if bg.isValid():
                        lum = 0.2126 * bg.redF() + 0.7152 * bg.greenF() + 0.0722 * bg.blueF()
                        # Return white on dark (lum<0.5), black on light
                        return QColor('#FFFFFF') if lum < 0.5 else QColor('#000000')
            except Exception:
                pass

            try:
                pal = QApplication.palette()
                # palette.window() returns a QBrush; call color()
                window_bg = pal.window().color()
                lum = 0.2126 * window_bg.redF() + 0.7152 * window_bg.greenF() + 0.0722 * window_bg.blueF()
                # Palette-based mapping: white on dark palette background
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
            # Whether to kekulize aromatic systems for 3D display
            'display_kekule_3d': False,
            'always_ask_charge': False,
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

