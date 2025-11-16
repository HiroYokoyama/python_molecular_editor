#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: Apache-2.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI 10.5281/zenodo.17268532
"""

#Version
VERSION = '1.15.0'

print("-----------------------------------------------------")
print("MoleditPy — A Python-based molecular editing software")
print("-----------------------------------------------------\n")

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


from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdGeometry
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdMolTransforms
from rdkit.DistanceGeometry import DoTriangleSmoothing


# Open Babel Python binding (optional; required for fallback)
try:
    from openbabel import pybel
    OBABEL_AVAILABLE = True
except Exception:
    # pybel (Open Babel Python bindings) is optional. If not present, disable OBabel features.
    pybel = None
    OBABEL_AVAILABLE = False
    print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")

# PyVista
import pyvista as pv
from pyvistaqt import QtInteractor

# --- Constants ---
ATOM_RADIUS = 18
BOND_OFFSET = 3.5
DEFAULT_BOND_LENGTH = 75 # テンプレートで使用する標準結合長
CLIPBOARD_MIME_TYPE = "application/x-moleditpy-fragment"

# Physical bond length (approximate) used to convert scene pixels to angstroms.
# DEFAULT_BOND_LENGTH is the length in pixels used in the editor UI for a typical bond.
# Many molecular file formats expect coordinates in angstroms; use ~1.5 Å as a typical single-bond length.
DEFAULT_BOND_LENGTH_ANGSTROM = 1.5
# Multiply pixel coordinates by this to get angstroms: ANGSTROM_PER_PIXEL = 1.5Å / DEFAULT_BOND_LENGTH(px)
ANGSTROM_PER_PIXEL = DEFAULT_BOND_LENGTH_ANGSTROM / DEFAULT_BOND_LENGTH

# UI / drawing / behavior constants (centralized for maintainability)
FONT_FAMILY = "Arial"
FONT_SIZE_LARGE = 20
FONT_SIZE_SMALL = 12
FONT_WEIGHT_BOLD = QFont.Weight.Bold

# Hit / visual sizes (in pixels at scale=1)
DESIRED_ATOM_PIXEL_RADIUS = 15.0
DESIRED_BOND_PIXEL_WIDTH = 18.0

# Bond/EZ label
EZ_LABEL_TEXT_OUTLINE = 2.5
EZ_LABEL_MARGIN = 16
EZ_LABEL_BOX_SIZE = 28

# Interaction thresholds
SNAP_DISTANCE = 14.0
SUM_TOLERANCE = 5.0

# Misc drawing
NUM_DASHES = 8
HOVER_PEN_WIDTH = 8

CPK_COLORS = {
    'H': QColor('#FFFFFF'), 'C': QColor('#222222'), 'N': QColor('#3377FF'), 'O': QColor('#FF3333'), 'F': QColor('#99E6E6'),
    'Cl': QColor('#33FF33'), 'Br': QColor('#A52A2A'), 'I': QColor('#9400D3'), 'S': QColor('#FFC000'), 'P': QColor('#FF8000'),
    'Si': QColor('#DAA520'), 'B': QColor('#FA8072'), 'He': QColor('#D9FFFF'), 'Ne': QColor('#B3E3F5'), 'Ar': QColor('#80D1E3'),
    'Kr': QColor('#5CACC8'), 'Xe': QColor('#429EB0'), 'Rn': QColor('#298FA2'), 'Li': QColor('#CC80FF'), 'Na': QColor('#AB5CF2'),
    'K': QColor('#8F44D7'), 'Rb': QColor('#702EBC'), 'Cs': QColor('#561B9E'), 'Fr': QColor('#421384'), 'Be': QColor('#C2FF00'),
    'Mg': QColor('#8AFF00'), 'Ca': QColor('#3DFF00'), 'Sr': QColor('#00FF00'), 'Ba': QColor('#00E600'), 'Ra': QColor('#00B800'),
    'Sc': QColor('#E6E6E6'), 'Ti': QColor('#BFC2C7'), 'V': QColor('#A6A6AB'), 'Cr': QColor('#8A99C7'), 'Mn': QColor('#9C7AC7'),
    'Fe': QColor('#E06633'), 'Co': QColor('#F090A0'), 'Ni': QColor('#50D050'), 'Cu': QColor('#C88033'), 'Zn': QColor('#7D80B0'),
    'Ga': QColor('#C28F8F'), 'Ge': QColor('#668F8F'), 'As': QColor('#BD80E3'), 'Se': QColor('#FFA100'), 'Tc': QColor('#3B9E9E'),
    'Ru': QColor('#248F8F'), 'Rh': QColor('#0A7D8F'), 'Pd': QColor('#006985'), 'Ag': QColor('#C0C0C0'), 'Cd': QColor('#FFD700'),
    'In': QColor('#A67573'), 'Sn': QColor('#668080'), 'Sb': QColor('#9E63B5'), 'Te': QColor('#D47A00'), 'La': QColor('#70D4FF'),
    'Ce': QColor('#FFFFC7'), 'Pr': QColor('#D9FFC7'), 'Nd': QColor('#C7FFC7'), 'Pm': QColor('#A3FFC7'), 'Sm': QColor('#8FFFC7'),
    'Eu': QColor('#61FFC7'), 'Gd': QColor('#45FFC7'), 'Tb': QColor('#30FFC7'), 'Dy': QColor('#1FFFC7'), 'Ho': QColor('#00FF9C'),
    'Er': QColor('#00E675'), 'Tm': QColor('#00D452'), 'Yb': QColor('#00BF38'), 'Lu': QColor('#00AB24'), 'Hf': QColor('#4DC2FF'),
    'Ta': QColor('#4DA6FF'), 'W': QColor('#2194D6'), 'Re': QColor('#267DAB'), 'Os': QColor('#266696'), 'Ir': QColor('#175487'),
    'Pt': QColor('#D0D0E0'), 'Au': QColor('#FFD123'), 'Hg': QColor('#B8B8D0'), 'Tl': QColor('#A6544D'), 'Pb': QColor('#575961'),
    'Bi': QColor('#9E4FB5'), 'Po': QColor('#AB5C00'), 'At': QColor('#754F45'), 'Ac': QColor('#70ABFA'), 'Th': QColor('#00BAFF'),
    'Pa': QColor('#00A1FF'), 'U': QColor('#008FFF'), 'Np': QColor('#0080FF'), 'Pu': QColor('#006BFF'), 'Am': QColor('#545CF2'),
    'Cm': QColor('#785CE3'), 'Bk': QColor('#8A4FE3'), 'Cf': QColor('#A136D4'), 'Es': QColor('#B31FD4'), 'Fm': QColor('#B31FBA'),
    'Md': QColor('#B30DA6'), 'No': QColor('#BD0D87'), 'Lr': QColor('#C70066'), 'Al': QColor('#B3A68F'), 'Y': QColor('#99FFFF'), 
    'Zr': QColor('#7EE7E7'), 'Nb': QColor('#68CFCE'), 'Mo': QColor('#52B7B7'), 'DEFAULT': QColor('#FF1493') # Pink fallback
}
CPK_COLORS_PV = {
    k: [c.redF(), c.greenF(), c.blueF()] for k, c in CPK_COLORS.items()
}

# Keep a copy of the original default map so we can restore it when user resets
DEFAULT_CPK_COLORS = {k: QColor(v) if not isinstance(v, QColor) else v for k, v in CPK_COLORS.items()}

pt = Chem.GetPeriodicTable()
VDW_RADII = {pt.GetElementSymbol(i): pt.GetRvdw(i) * 0.3 for i in range(1, 119)}


def sip_isdeleted_safe(obj) -> bool:
    """Return True if sip reports the given wrapper object as deleted.

    This function is conservative: if SIP isn't available or any error
    occurs while checking, it returns False (i.e. not deleted) so that the
    caller can continue other lightweight guards (like checking scene()).
    """
    try:
        if _sip_isdeleted is None:
            return False
        return bool(_sip_isdeleted(obj))
    except Exception:
        return False

class Dialog3DPickingMixin:
    """3D原子選択のための共通機能を提供するMixin"""
    
    def __init__(self):
        """Mixinの初期化"""
        self.picking_enabled = False
    
    def eventFilter(self, obj, event):
        """3Dビューでのマウスクリックをキャプチャする（元の3D editロジックを正確に再現）"""
        if (obj == self.main_window.plotter.interactor and 
            event.type() == QEvent.Type.MouseButtonPress and 
            event.button() == Qt.MouseButton.LeftButton):
            
            try:
                # VTKイベント座標を取得（元のロジックと同じ）
                interactor = self.main_window.plotter.interactor
                click_pos = interactor.GetEventPosition()
                picker = self.main_window.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer)

                if picker.GetActor() is self.main_window.atom_actor:
                    picked_position = np.array(picker.GetPickPosition())
                    distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                    closest_atom_idx = np.argmin(distances)

                    # 範囲チェックを追加
                    if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                        # クリック閾値チェック（元のロジックと同じ）
                        atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                        if atom:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            click_threshold = vdw_radius * 1.5

                            if distances[closest_atom_idx] < click_threshold:
                                # We handled the pick (atom clicked) -> consume the event so
                                # other UI elements (including the VTK interactor observers)
                                # don't also process it. Set a flag on the main window so
                                # the VTK-based handlers can ignore the same logical click
                                # when it arrives via the VTK event pipeline.
                                try:
                                    self.main_window._picking_consumed = True
                                except Exception:
                                    pass
                                self.on_atom_picked(int(closest_atom_idx))
                                return True
                
                # 原子以外をクリックした場合は選択をクリア（Measurementモードと同じロジック）
                if hasattr(self, 'clear_selection'):
                    self.clear_selection()
                # We did not actually pick an atom; do NOT consume the event here so
                # the interactor and CustomInteractorStyle can handle camera rotation
                # and other behaviors. Returning False (or calling the base
                # implementation) allows normal processing to continue.
                return False
                    
            except Exception as e:
                print(f"Error in eventFilter: {e}")
                # On exception, don't swallow the event either — let the normal
                # event pipeline continue so the UI remains responsive.
                return False

        return super().eventFilter(obj, event)
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
        # Ensure the main window flag exists
        try:
            self.main_window._picking_consumed = False
        except Exception:
            pass
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, 'picking_enabled') and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        try:
            # Clear any leftover flag when picking is disabled
            if hasattr(self.main_window, '_picking_consumed'):
                self.main_window._picking_consumed = False
        except Exception:
            pass
    
    def try_alternative_picking(self, x, y):
        """代替のピッキング方法（使用しない）"""
        pass

class TemplatePreviewView(QGraphicsView):
    """テンプレートプレビュー用のカスタムビュークラス"""
    
    def __init__(self, scene):
        super().__init__(scene)
        self.original_scene_rect = None
        self.template_data = None  # Store template data for dynamic redrawing
        self.parent_dialog = None  # Reference to parent dialog for redraw access
    
    def set_template_data(self, template_data, parent_dialog):
        """テンプレートデータと親ダイアログの参照を設定"""
        self.template_data = template_data
        self.parent_dialog = parent_dialog
    
    def resizeEvent(self, event):
        """リサイズイベントを処理してプレビューを再フィット"""
        super().resizeEvent(event)
        if self.original_scene_rect and not self.original_scene_rect.isEmpty():
            # Delay the fitInView call to ensure proper widget sizing
            QTimer.singleShot(10, self.refit_view)
    
    def refit_view(self):
        """ビューを再フィット"""
        try:
            if self.original_scene_rect and not self.original_scene_rect.isEmpty():
                self.fitInView(self.original_scene_rect, Qt.AspectRatioMode.KeepAspectRatio)
        except Exception as e:
            print(f"Warning: Failed to refit template preview: {e}")
    
    def showEvent(self, event):
        """表示イベントを処理"""
        super().showEvent(event)
        # Ensure proper fitting when widget becomes visible
        if self.original_scene_rect:
            QTimer.singleShot(50, self.refit_view)
    
    def redraw_with_current_size(self):
        """現在のサイズに合わせてテンプレートを再描画"""
        if self.template_data and self.parent_dialog:
            try:
                # Clear current scene
                self.scene().clear()
                
                # Redraw with current view size for proper fit-based scaling
                view_size = (self.width(), self.height())
                self.parent_dialog.draw_template_preview(self.scene(), self.template_data, view_size)
                
                # Refit the view
                bounding_rect = self.scene().itemsBoundingRect()
                if not bounding_rect.isEmpty() and bounding_rect.width() > 0 and bounding_rect.height() > 0:
                    content_size = max(bounding_rect.width(), bounding_rect.height())
                    padding = max(20, content_size * 0.2)
                    padded_rect = bounding_rect.adjusted(-padding, -padding, padding, padding)
                    self.scene().setSceneRect(padded_rect)
                    self.original_scene_rect = padded_rect
                    QTimer.singleShot(10, lambda: self.fitInView(padded_rect, Qt.AspectRatioMode.KeepAspectRatio))
            except Exception as e:
                print(f"Warning: Failed to redraw template preview: {e}")

class UserTemplateDialog(QDialog):
    """ユーザーテンプレート管理ダイアログ"""
    
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.user_templates = []
        self.selected_template = None
        self.init_ui()
        self.load_user_templates()
    
    def init_ui(self):
        self.setWindowTitle("User Templates")
        self.setModal(False)  # モードレスに変更
        self.resize(800, 600)
        
        # ウィンドウを右上に配置
        if self.parent():
            parent_geometry = self.parent().geometry()
            x = parent_geometry.right() - self.width() - 20
            y = parent_geometry.top() + 50
            self.move(x, y)
        
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Create and manage your custom molecular templates. Click a template to use it in the editor.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Template grid
        self.template_widget = QWidget()
        self.template_layout = QGridLayout(self.template_widget)
        self.template_layout.setSpacing(10)
        
        scroll_area = QScrollArea()
        scroll_area.setWidget(self.template_widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumHeight(400)
        layout.addWidget(scroll_area)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.save_current_button = QPushButton("Save Current 2D as Template")
        self.save_current_button.clicked.connect(self.save_current_as_template)
        button_layout.addWidget(self.save_current_button)
        
        button_layout.addStretch()
        
        self.delete_button = QPushButton("Delete Selected")
        self.delete_button.clicked.connect(self.delete_selected_template)
        self.delete_button.setEnabled(False)
        button_layout.addWidget(self.delete_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
    
    def resizeEvent(self, event):
        """ダイアログリサイズ時にテンプレートプレビューを再フィット"""
        super().resizeEvent(event)
        # Delay the refit to ensure proper widget sizing
        QTimer.singleShot(100, self.refit_all_previews)
    
    def refit_all_previews(self):
        """すべてのテンプレートプレビューを再フィット"""
        try:
            for i in range(self.template_layout.count()):
                item = self.template_layout.itemAt(i)
                if item and item.widget():
                    widget = item.widget()
                    # Find the TemplatePreviewView within this widget
                    for child in widget.findChildren(TemplatePreviewView):
                        if hasattr(child, 'redraw_with_current_size'):
                            # Use redraw for better scaling adaptation
                            child.redraw_with_current_size()
                        elif hasattr(child, 'refit_view'):
                            child.refit_view()
        except Exception as e:
            print(f"Warning: Failed to refit template previews: {e}")
    
    def showEvent(self, event):
        """ダイアログ表示時にプレビューを適切にフィット"""
        super().showEvent(event)
        # Ensure all previews are properly fitted when dialog becomes visible
        QTimer.singleShot(300, self.refit_all_previews)
    
    def get_template_directory(self):
        """テンプレートディレクトリのパスを取得"""
        template_dir = os.path.join(self.main_window.settings_dir, 'user-templates')
        if not os.path.exists(template_dir):
            os.makedirs(template_dir)
        return template_dir
    
    def load_user_templates(self):
        """ユーザーテンプレートを読み込み"""
        template_dir = self.get_template_directory()
        self.user_templates.clear()
        
        try:
            for filename in os.listdir(template_dir):
                if filename.endswith('.pmetmplt'):
                    filepath = os.path.join(template_dir, filename)
                    template_data = self.load_template_file(filepath)
                    if template_data:
                        template_data['filename'] = filename
                        template_data['filepath'] = filepath
                        self.user_templates.append(template_data)
        except Exception as e:
            print(f"Error loading user templates: {e}")
        
        self.update_template_grid()
    
    def load_template_file(self, filepath):
        """テンプレートファイルを読み込み"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading template file {filepath}: {e}")
            return None
    
    def save_template_file(self, filepath, template_data):
        """テンプレートファイルを保存"""
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)
            return True
        except Exception as e:
            print(f"Error saving template file {filepath}: {e}")
            return False
    
    def update_template_grid(self):
        """テンプレートグリッドを更新"""
        # Clear existing widgets
        for i in reversed(range(self.template_layout.count())):
            self.template_layout.itemAt(i).widget().setParent(None)
        
        # Add template previews (left-to-right, top-to-bottom ordering)
        cols = 4
        for i, template in enumerate(self.user_templates):
            row = i // cols
            col = i % cols  # Left-to-right ordering
            
            preview_widget = self.create_template_preview(template)
            self.template_layout.addWidget(preview_widget, row, col)
        
        # Ensure all previews are properly fitted after grid update
        QTimer.singleShot(200, self.refit_all_previews)
    
    def create_template_preview(self, template_data):
        """テンプレートプレビューウィジェットを作成"""
        widget = QWidget()
        widget.setFixedSize(180, 200)
        widget.setStyleSheet("""
            QWidget {
                border: 2px solid #ccc;
                border-radius: 8px;
                background-color: white;
            }
            QWidget:hover {
                border-color: #007acc;
                background-color: #f0f8ff;
            }
        """)
        
        layout = QVBoxLayout(widget)
        
        # Preview graphics - use custom view class for better resize handling
        preview_scene = QGraphicsScene()
        preview_view = TemplatePreviewView(preview_scene)
        preview_view.setFixedSize(160, 140)
        preview_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        preview_view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        preview_view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        
        # Set template data for dynamic redrawing
        preview_view.set_template_data(template_data, self)
        
        # Draw template structure with view size for proper scaling
        view_size = (preview_view.width(), preview_view.height())
        self.draw_template_preview(preview_scene, template_data, view_size)
        
        # Improved fitting approach with better error handling
        bounding_rect = preview_scene.itemsBoundingRect()
        if not bounding_rect.isEmpty() and bounding_rect.width() > 0 and bounding_rect.height() > 0:
            # Calculate appropriate padding based on content size
            content_size = max(bounding_rect.width(), bounding_rect.height())
            padding = max(20, content_size * 0.2)  # At least 20 units or 20% of content
            
            padded_rect = bounding_rect.adjusted(-padding, -padding, padding, padding)
            preview_scene.setSceneRect(padded_rect)
            
            # Store original scene rect for proper fitting on resize
            preview_view.original_scene_rect = padded_rect
            
            # Use QTimer to ensure fitInView happens after widget is fully initialized
            QTimer.singleShot(0, lambda: self.fit_preview_view_safely(preview_view, padded_rect))
        else:
            # Default view for empty or invalid content
            default_rect = QRectF(-50, -50, 100, 100)
            preview_scene.setSceneRect(default_rect)
            preview_view.original_scene_rect = default_rect
            QTimer.singleShot(0, lambda: self.fit_preview_view_safely(preview_view, default_rect))
        
        layout.addWidget(preview_view)
        
        # Template name
        name = template_data.get('name', 'Unnamed Template')
        name_label = QLabel(name)
        name_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        name_label.setWordWrap(True)
        layout.addWidget(name_label)
        
        # Mouse events
        widget.mousePressEvent = lambda event: self.select_template(template_data, widget)
        widget.mouseDoubleClickEvent = lambda event: self.use_template(template_data)
        
        return widget
    
    def fit_preview_view_safely(self, view, rect):
        """プレビュービューを安全にフィット"""
        try:
            if view and not rect.isEmpty():
                view.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)
        except Exception as e:
            print(f"Warning: Failed to fit preview view: {e}")
    
    def draw_template_preview(self, scene, template_data, view_size=None):
        """テンプレートプレビューを描画 - fitInView縮小率に基づく動的スケーリング"""
        atoms = template_data.get('atoms', [])
        bonds = template_data.get('bonds', [])
        
        if not atoms:
            # Add placeholder text when no atoms
            text = scene.addText("No structure", QFont('Arial', 12))
            text.setDefaultTextColor(QColor('gray'))
            return
        
        # Calculate molecular dimensions
        positions = [QPointF(atom['x'], atom['y']) for atom in atoms]
        min_x = min(pos.x() for pos in positions)
        max_x = max(pos.x() for pos in positions)
        min_y = min(pos.y() for pos in positions)
        max_y = max(pos.y() for pos in positions)
        
        mol_width = max_x - min_x
        mol_height = max_y - min_y
        mol_size = max(mol_width, mol_height)
        
        # Calculate fit scale factor (how much fitInView will shrink the content)
        if view_size is None:
            view_size = (160, 140)  # Default preview view size
        
        view_width, view_height = view_size
        
        if mol_size > 0 and mol_width > 0 and mol_height > 0:
            # Calculate the padding that will be added
            padding = max(20, mol_size * 0.2)
            padded_width = mol_width + 2 * padding
            padded_height = mol_height + 2 * padding
            
            # Calculate how much fitInView will scale down the content
            # fitInView fits the padded rectangle into the view while maintaining aspect ratio
            fit_scale_x = view_width / padded_width
            fit_scale_y = view_height / padded_height
            fit_scale = min(fit_scale_x, fit_scale_y)  # fitInView uses the smaller scale
            
            # Compensate for the fit scaling to maintain visual thickness
            # When fit_scale is small (content heavily shrunk), we need thicker lines/fonts
            if fit_scale > 0:
                scale_factor = max(0.4, min(4.0, 1.0 / fit_scale))
            else:
                scale_factor = 4.0
            
            # Debug info (can be removed in production)
            # print(f"Mol size: {mol_size:.1f}, Fit scale: {fit_scale:.3f}, Scale factor: {scale_factor:.2f}")
        else:
            scale_factor = 1.0
        
        # Base sizes that look good at 1:1 scale
        base_bond_width = 1.8
        base_font_size = 11
        base_ellipse_width = 18
        base_ellipse_height = 14
        base_double_bond_offset = 3.5
        base_triple_bond_offset = 2.5
        
        # Apply inverse fit scaling to maintain visual consistency
        bond_width = max(1.0, min(8.0, base_bond_width * scale_factor))
        font_size = max(8, min(24, int(base_font_size * scale_factor)))
        ellipse_width = max(10, min(40, base_ellipse_width * scale_factor))
        ellipse_height = max(8, min(30, base_ellipse_height * scale_factor))
        double_bond_offset = max(2.0, min(10.0, base_double_bond_offset * scale_factor))
        triple_bond_offset = max(1.5, min(8.0, base_triple_bond_offset * scale_factor))
        
        # Create atom ID to index mapping for bond drawing
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get('id', i)  # Use id if available, otherwise use index
            atom_id_to_index[atom_id] = i
        
        # Draw bonds first using original coordinates with dynamic sizing
        for bond in bonds:
            atom1_id, atom2_id = bond['atom1'], bond['atom2']
            
            # Get atom indices from IDs
            atom1_idx = atom_id_to_index.get(atom1_id)
            atom2_idx = atom_id_to_index.get(atom2_id)
            
            if atom1_idx is not None and atom2_idx is not None and atom1_idx < len(atoms) and atom2_idx < len(atoms):
                pos1 = QPointF(atoms[atom1_idx]['x'], atoms[atom1_idx]['y'])
                pos2 = QPointF(atoms[atom2_idx]['x'], atoms[atom2_idx]['y'])
                
                # Draw bonds with proper order - dynamic thickness
                bond_order = bond.get('order', 1)
                pen = QPen(QColor('black'), bond_width)
                
                if bond_order == 2:
                    # Double bond - draw two parallel lines
                    line = QLineF(pos1, pos2)
                    if line.length() > 0:
                        normal = line.normalVector()
                        normal.setLength(double_bond_offset)
                        
                        line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                        line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                        
                        scene.addLine(line1, pen)
                        scene.addLine(line2, pen)
                    else:
                        scene.addLine(line, pen)
                elif bond_order == 3:
                    # Triple bond - draw three parallel lines
                    line = QLineF(pos1, pos2)
                    if line.length() > 0:
                        normal = line.normalVector()
                        normal.setLength(triple_bond_offset)
                        
                        # Center line
                        scene.addLine(line, pen)
                        # Side lines
                        line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                        line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                        
                        scene.addLine(line1, pen)
                        scene.addLine(line2, pen)
                    else:
                        scene.addLine(line, pen)
                else:
                    # Single bond
                    scene.addLine(QLineF(pos1, pos2), pen)
        
        # Draw only non-carbon atom labels with dynamic sizing
        for i, atom in enumerate(atoms):
            try:
                pos = QPointF(atom['x'], atom['y'])
                symbol = atom.get('symbol', 'C')
                
                # Draw atoms - white ellipse background to hide bonds, then CPK colored text
                if symbol != 'C':
                    # All non-carbon atoms including hydrogen: white background ellipse + CPK colored text
                    color = CPK_COLORS.get(symbol, CPK_COLORS.get('DEFAULT', QColor('#FF1493')))
                    
                    # Add white background ellipse to hide bonds - dynamic size
                    pen = QPen(Qt.GlobalColor.white, 0)  # No border
                    brush = QBrush(Qt.GlobalColor.white)
                    ellipse_x = pos.x() - ellipse_width/2
                    ellipse_y = pos.y() - ellipse_height/2
                    ellipse = scene.addEllipse(ellipse_x, ellipse_y, ellipse_width, ellipse_height, pen, brush)
                    
                    # Add CPK colored text label on top - dynamic font size
                    font = QFont("Arial", font_size, QFont.Weight.Bold)
                    text = scene.addText(symbol, font)
                    text.setDefaultTextColor(color)  # CPK colored text
                    text_rect = text.boundingRect()
                    text.setPos(pos.x() - text_rect.width()/2, pos.y() - text_rect.height()/2)
                    
            except Exception as e:
                continue
    
    def select_template(self, template_data, widget):
        """テンプレートを選択してテンプレートモードに切り替え"""
        # Clear previous selection styling
        for i in range(self.template_layout.count()):
            item = self.template_layout.itemAt(i)
            if item and item.widget():
                item.widget().setStyleSheet("""
                    QWidget {
                        border: 2px solid #ccc;
                        border-radius: 8px;
                        background-color: white;
                    }
                    QWidget:hover {
                        border-color: #007acc;
                        background-color: #f0f8ff;
                    }
                """)
        
        # Highlight selected widget - only border, no background change
        widget.setStyleSheet("""
            QWidget {
                border: 3px solid #007acc;
                border-radius: 8px;
                background-color: white;
            }
        """)
        
        self.selected_template = template_data
        self.delete_button.setEnabled(True)

        # Automatically switch to template mode when template is selected
        template_name = template_data.get('name', 'user_template')
        mode_name = f"template_user_{template_name}"

        # Store template data for the scene to use
        try:
            self.main_window.scene.user_template_data = template_data
        except Exception:
            # Best-effort: ignore if scene or attribute missing
            pass

        # Force the main window into the template mode.
        # Clear or uncheck any existing mode actions if present to avoid staying in another mode.
        try:
            # Uncheck all mode actions first (if a dict of QAction exists)
            if hasattr(self.main_window, 'mode_actions') and isinstance(self.main_window.mode_actions, dict):
                for act in self.main_window.mode_actions.values():
                    try:
                        act.setChecked(False)
                    except Exception:
                        continue

            # If main_window has a set_mode method, call it. Otherwise, try to set a mode attribute.
            if hasattr(self.main_window, 'set_mode') and callable(self.main_window.set_mode):
                self.main_window.set_mode(mode_name)
            else:
                # Fallback: set an attribute and try to update UI
                setattr(self.main_window, 'mode', mode_name)

            # Update UI
            try:
                self.main_window.statusBar().showMessage(f"Template mode: {template_name}")
            except Exception:
                # ignore status bar failures
                pass

            # If there is a matching QAction in mode_actions, check it
            try:
                if hasattr(self.main_window, 'mode_actions') and f"template_user_{template_name}" in self.main_window.mode_actions:
                    self.main_window.mode_actions[f"template_user_{template_name}"].setChecked(True)
            except Exception:
                pass
        except Exception as e:
            print(f"Warning: Failed to switch main window to template mode: {e}")
    
    def use_template(self, template_data):
        """テンプレートを使用（エディタに適用）"""
        try:
            # Switch to template mode
            template_name = template_data.get('name', 'user_template')
            mode_name = f"template_user_{template_name}"
            
            # Store template data for the scene to use
            try:
                self.main_window.scene.user_template_data = template_data
            except Exception:
                pass

            # Force the main window into the template mode (same approach as select_template)
            try:
                if hasattr(self.main_window, 'mode_actions') and isinstance(self.main_window.mode_actions, dict):
                    for act in self.main_window.mode_actions.values():
                        try:
                            act.setChecked(False)
                        except Exception:
                            continue

                if hasattr(self.main_window, 'set_mode') and callable(self.main_window.set_mode):
                    self.main_window.set_mode(mode_name)
                else:
                    setattr(self.main_window, 'mode', mode_name)

                try:
                    self.main_window.statusBar().showMessage(f"Template mode: {template_name}")
                except Exception:
                    pass

                # Mark selected and keep dialog open
                self.selected_template = template_data
            except Exception as e:
                print(f"Warning: Failed to switch main window to template mode: {e}")

            # Don't close dialog - keep it open for easy template switching
            # self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply template: {str(e)}")
    
    def save_current_as_template(self):
        """現在の2D構造をテンプレートとして保存"""
        if not self.main_window.data.atoms:
            QMessageBox.warning(self, "Warning", "No structure to save as template.")
            return
        
        # Get template name
        name, ok = QInputDialog.getText(self, "Save Template", "Enter template name:")
        if not ok or not name.strip():
            return
        
        name = name.strip()
        
        try:
            # Convert current structure to template format
            template_data = self.convert_structure_to_template(name)
            
            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(self.get_template_directory(), filename)
            
            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self, "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return
            
            if self.save_template_file(filepath, template_data):
                # Mark main window as saved
                self.main_window.has_unsaved_changes = False
                self.main_window.update_window_title()
                
                QMessageBox.information(self, "Success", f"Template '{name}' saved successfully.")
                self.load_user_templates()  # Refresh the display
            else:
                QMessageBox.critical(self, "Error", "Failed to save template.")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save template: {str(e)}")
    
    def convert_structure_to_template(self, name):
        """現在の構造をテンプレート形式に変換"""
        atoms_data = []
        bonds_data = []
        
        # Convert atoms
        for atom_id, atom_info in self.main_window.data.atoms.items():
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
        for (atom1_id, atom2_id), bond_info in self.main_window.data.bonds.items():
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

        return template_data
    
    def delete_selected_template(self):
        """選択されたテンプレートを削除"""
        if not self.selected_template:
            return
        
        name = self.selected_template.get('name', 'Unknown')
        reply = QMessageBox.question(
            self, "Delete Template",
            f"Are you sure you want to delete template '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            try:
                filepath = self.selected_template['filepath']
                os.remove(filepath)
                QMessageBox.information(self, "Success", f"Template '{name}' deleted successfully.")
                self.load_user_templates()  # Refresh the display
                self.selected_template = None
                self.delete_button.setEnabled(False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to delete template: {str(e)}")

class AboutDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("About MoleditPy")
        self.setFixedSize(250, 300)
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Create a clickable image label
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        # Load the original icon image
        icon_path = os.path.join(os.path.dirname(__file__), 'assets', 'icon.png')
        if os.path.exists(icon_path):
            original_pixmap = QPixmap(icon_path)
            # Scale to 2x size (160x160)
            pixmap = original_pixmap.scaled(160, 160, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
        else:
            # Fallback: create a simple placeholder if icon.png not found
            pixmap = QPixmap(160, 160)
            pixmap.fill(Qt.GlobalColor.lightGray)
            painter = QPainter(pixmap)
            painter.setPen(QPen(Qt.GlobalColor.black, 2))
            painter.drawText(pixmap.rect(), Qt.AlignmentFlag.AlignCenter, "MoleditPy")
            painter.end()
        
        self.image_label.setPixmap(pixmap)
        try:
            self.image_label.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
        except Exception:
            pass

        self.image_label.mousePressEvent = self.image_mouse_press_event
        
        layout.addWidget(self.image_label)
        
        # Add text information
        info_text = f"MoleditPy Ver. {VERSION}\nAuthor: Hiromichi Yokoyama\nLicense: Apache-2.0"
        info_label = QLabel(info_text)
        info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(info_label)
        
        # Add OK button
        ok_button = QPushButton("OK")
        ok_button.setFixedSize(80, 30)  # 小さいサイズに固定
        ok_button.clicked.connect(self.accept)
        
        # Center the button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(ok_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)
    
    def image_clicked(self, event):
        """Easter egg: Clear all and load bipyrimidine from SMILES"""
        # Clear the current scene
        self.main_window.clear_all()

        bipyrimidine_smiles = "C1=CN=C(N=C1)C2=NC=CC=N2"
        self.main_window.load_from_smiles(bipyrimidine_smiles)

        # Close the dialog
        self.accept()

    def image_mouse_press_event(self, event):
        """Handle mouse press on the image: trigger easter egg only for right-click."""
        try:
            if event.button() == Qt.MouseButton.RightButton:
                self.image_clicked(event)
            else:
                event.ignore()
        except Exception:
            try:
                event.ignore()
            except Exception:
                pass

class TranslationDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()  # 複数原子選択用
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Translation")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click atoms in the 3D view to select them. The centroid of selected atoms will be moved to the target coordinates, translating the entire molecule.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Coordinate inputs
        coord_layout = QGridLayout()
        coord_layout.addWidget(QLabel("Target X:"), 0, 0)
        self.x_input = QLineEdit("0.0")
        coord_layout.addWidget(self.x_input, 0, 1)
        
        coord_layout.addWidget(QLabel("Target Y:"), 1, 0)
        self.y_input = QLineEdit("0.0")
        coord_layout.addWidget(self.y_input, 1, 1)
        
        coord_layout.addWidget(QLabel("Target Z:"), 2, 0)
        self.z_input = QLineEdit("0.0")
        coord_layout.addWidget(self.z_input, 2, 1)
        
        layout.addLayout(coord_layout)

        # Translation target toggle: Entire molecule (default) or Selected atoms only
        self.translate_selected_only_checkbox = QCheckBox("Translate selected atoms only")
        self.translate_selected_only_checkbox.setToolTip(
            "When checked, only the atoms you selected will be moved so their centroid matches the target.\n"
            "When unchecked (default), the entire molecule will be translated so the selected atoms' centroid moves to the target."
        )
        self.translate_selected_only_checkbox.setChecked(False)  # default: entire molecule
        layout.addWidget(self.translate_selected_only_checkbox)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
    
        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for translation")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Translation")
        self.apply_button.clicked.connect(self.apply_translation)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_translation()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def update_display(self):
        """表示を更新"""
        if not self.selected_atoms:
            self.selection_label.setText("No atoms selected")
            self.apply_button.setEnabled(False)
        else:
            # 分子の有効性チェック
            if not self.mol or self.mol.GetNumConformers() == 0:
                self.selection_label.setText("Error: No valid molecule or conformer")
                self.apply_button.setEnabled(False)
                return
            
            try:
                conf = self.mol.GetConformer()
                # 選択原子の重心を計算
                centroid = self.calculate_centroid()
                
                # 選択原子の情報を表示
                atom_info = []
                for atom_idx in sorted(self.selected_atoms):
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    atom_info.append(f"{symbol}({atom_idx})")
                
                self.selection_label.setText(
                    f"Selected atoms: {', '.join(atom_info)}\n"
                    f"Centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})"
                )
                self.apply_button.setEnabled(True)
            except Exception as e:
                self.selection_label.setText(f"Error accessing atom data: {str(e)}")
                self.apply_button.setEnabled(False)

        # Update the coordinate input fields when selection changes
        # If there are selected atoms, fill the inputs with the computed centroid (or single atom pos).
        # If no atoms selected, clear or reset to 0.0.
        try:
            if self.selected_atoms:
                # Use the centroid we just computed if available; otherwise compute now.
                try:
                    coords = centroid
                except NameError:
                    coords = self.calculate_centroid()

                # Format with reasonable precision
                self.x_input.setText(f"{coords[0]:.4f}")
                self.y_input.setText(f"{coords[1]:.4f}")
                self.z_input.setText(f"{coords[2]:.4f}")
            else:
                # No selection: reset fields to default
                self.x_input.setText("0.0")
                self.y_input.setText("0.0")
                self.z_input.setText("0.0")
        except Exception:
            # Be tolerant: do not crash the UI if inputs cannot be updated
            pass
    
    def calculate_centroid(self):
        """選択原子の重心を計算"""
        if not self.selected_atoms:
            return np.array([0.0, 0.0, 0.0])
        
        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.selected_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        
        return np.mean(positions, axis=0)
    
    def apply_translation(self):
        """平行移動を適用"""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select at least one atom.")
            return

        # 分子の有効性チェック
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Warning", "No valid molecule or conformer available.")
            return

        try:
            target_x = float(self.x_input.text())
            target_y = float(self.y_input.text())
            target_z = float(self.z_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid coordinates.")
            return

        try:
            # 選択原子の重心を計算
            current_centroid = self.calculate_centroid()
            target_pos = np.array([target_x, target_y, target_z])

            # 移動ベクトルを計算
            translation_vector = target_pos - current_centroid

            conf = self.mol.GetConformer()

            if self.translate_selected_only_checkbox.isChecked():
                # Move only the selected atoms: shift selected atoms by translation_vector
                for i in range(self.mol.GetNumAtoms()):
                    if i in self.selected_atoms:
                        atom_pos = np.array(conf.GetAtomPosition(i))
                        new_pos = atom_pos + translation_vector
                        conf.SetAtomPosition(i, new_pos.tolist())
                        # Update 3d positions for this atom only
                        try:
                            self.main_window.atom_positions_3d[i] = new_pos
                        except Exception:
                            pass
                    else:
                        # leave other atoms unchanged
                        continue
            else:
                # Default: translate entire molecule so centroid moves to target
                for i in range(self.mol.GetNumAtoms()):
                    atom_pos = np.array(conf.GetAtomPosition(i))
                    new_pos = atom_pos + translation_vector
                    conf.SetAtomPosition(i, new_pos.tolist())
                    self.main_window.atom_positions_3d[i] = new_pos


            # 3D表示を更新
            self.main_window.draw_molecule_3d(self.mol)

            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Apply後に選択解除
            self.clear_selection()

            # Undo状態を保存
            self.main_window.push_undo_state()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply translation: {str(e)}")
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
    
    def select_all_atoms(self):
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, 'mol') and self.mol is not None:
                try:
                    n = self.mol.GetNumAtoms()
                    # create a set of indices [0..n-1]
                    self.selected_atoms = set(range(n))
                except Exception:
                    # fallback to main_window data map
                    self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()
            else:
                # fallback to main_window data map
                self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()

        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []

        if self.selected_atoms:
            positions = []
            labels = []

            for i, atom_idx in enumerate(sorted(self.selected_atoms)):
                pos = self.main_window.atom_positions_3d[atom_idx]
                positions.append(pos)
                labels.append(f"S{i+1}")

            # 重心位置も表示
            if len(self.selected_atoms) > 1:
                centroid = self.calculate_centroid()
                positions.append(centroid)
                labels.append("CEN")

            # ラベルを追加
            if positions:
                label_actor = self.main_window.plotter.add_point_labels(
                    positions, labels,
                    point_size=20,
                    font_size=12,
                    text_color='cyan',
                    always_visible=True
                )
                # add_point_labelsがリストを返す場合も考慮
                if isinstance(label_actor, list):
                    self.selection_labels.extend(label_actor)
                else:
                    self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except Exception:
                    pass
            self.selection_labels = []
        # ラベル消去後に再描画を強制
        try:
            self.main_window.plotter.render()
        except Exception:
            pass
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().accept()


class MirrorDialog(QDialog):
    """分子の鏡像を作成するダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.main_window = main_window
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Mirror Molecule")
        self.setMinimumSize(300, 200)
        
        layout = QVBoxLayout(self)
        
        # 説明テキスト
        info_label = QLabel("Select the mirror plane to create molecular mirror image:")
        layout.addWidget(info_label)
        
        # ミラー平面選択のラジオボタン
        self.plane_group = QButtonGroup(self)
        
        self.xy_radio = QRadioButton("XY plane (Z = 0)")
        self.xz_radio = QRadioButton("XZ plane (Y = 0)")
        self.yz_radio = QRadioButton("YZ plane (X = 0)")
        
        self.xy_radio.setChecked(True)  # デフォルト選択
        
        self.plane_group.addButton(self.xy_radio, 0)
        self.plane_group.addButton(self.xz_radio, 1)
        self.plane_group.addButton(self.yz_radio, 2)
        
        layout.addWidget(self.xy_radio)
        layout.addWidget(self.xz_radio)
        layout.addWidget(self.yz_radio)
        
        layout.addSpacing(20)
        
        # ボタン
        button_layout = QHBoxLayout()
        
        apply_button = QPushButton("Apply Mirror")
        apply_button.clicked.connect(self.apply_mirror)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)

        button_layout.addWidget(apply_button)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
    
    def apply_mirror(self):
        """選択された平面に対してミラー変換を適用"""
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No 3D coordinates available.")
            return
        
        # 選択された平面を取得
        plane_id = self.plane_group.checkedId()
        
        try:
            conf = self.mol.GetConformer()
            
            # 各原子の座標を変換
            for atom_idx in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(atom_idx)
                
                if plane_id == 0:  # XY平面（Z軸に対してミラー）
                    new_pos = [pos.x, pos.y, -pos.z]
                elif plane_id == 1:  # XZ平面（Y軸に対してミラー）
                    new_pos = [pos.x, -pos.y, pos.z]
                elif plane_id == 2:  # YZ平面（X軸に対してミラー）
                    new_pos = [-pos.x, pos.y, pos.z]
                
                # 新しい座標を設定
                from rdkit.Geometry import Point3D
                conf.SetAtomPosition(atom_idx, Point3D(new_pos[0], new_pos[1], new_pos[2]))
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # ミラー変換後にキラルタグを強制的に再計算
            try:
                if self.mol.GetNumConformers() > 0:
                    # 既存のキラルタグをクリア
                    for atom in self.mol.GetAtoms():
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                    # 3D座標から新しいキラルタグを計算
                    Chem.AssignAtomChiralTagsFromStructure(self.mol, confId=0)
            except Exception as e:
                print(f"Error updating chiral tags: {e}")
            
            # キラルラベルを更新（鏡像変換でキラリティが変わる可能性があるため）
            self.main_window.update_chiral_labels()
            
            self.main_window.push_undo_state()
            
            plane_names = ["XY", "XZ", "YZ"]
            self.main_window.statusBar().showMessage(f"Molecule mirrored across {plane_names[plane_id]} plane.")
        

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply mirror transformation: {str(e)}")

class MoveGroupDialog(Dialog3DPickingMixin, QDialog):
    """結合している分子グループを選択して並行移動・回転するダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()
        self.group_atoms = set()  # 選択原子に結合している全原子
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Move Group")
        self.setModal(False)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        self.resize(300,400)  # ウィンドウサイズを設定
        layout = QVBoxLayout(self)
        
        # ドラッグ状態管理
        self.is_dragging_group = False
        self.drag_start_pos = None
        self.mouse_moved_during_drag = False  # ドラッグ中にマウスが動いたかを追跡
        
        # Instructions
        instruction_label = QLabel("Click an atom in the 3D view to select its connected molecule group.\n"
                                   "Left-drag: Move the group\n"
                                   "Right-drag: Rotate the group around its center")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected group display
        self.selection_label = QLabel("No group selected")
        layout.addWidget(self.selection_label)
        
        # Translation controls
        trans_group = QLabel("Translation (Å):")
        trans_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(trans_group)
        
        trans_layout = QGridLayout()
        self.x_trans_input = QLineEdit("0.0")
        self.y_trans_input = QLineEdit("0.0")
        self.z_trans_input = QLineEdit("0.0")
        
        # Enterキーでapply_translationを実行
        self.x_trans_input.returnPressed.connect(self.apply_translation)
        self.y_trans_input.returnPressed.connect(self.apply_translation)
        self.z_trans_input.returnPressed.connect(self.apply_translation)
        
        trans_layout.addWidget(QLabel("X:"), 0, 0)
        trans_layout.addWidget(self.x_trans_input, 0, 1)
        trans_layout.addWidget(QLabel("Y:"), 1, 0)
        trans_layout.addWidget(self.y_trans_input, 1, 1)
        trans_layout.addWidget(QLabel("Z:"), 2, 0)
        trans_layout.addWidget(self.z_trans_input, 2, 1)
        
        trans_button_layout = QHBoxLayout()
        reset_trans_button = QPushButton("Reset")
        reset_trans_button.clicked.connect(self.reset_translation_inputs)
        trans_button_layout.addWidget(reset_trans_button)
        
        apply_trans_button = QPushButton("Apply Translation")
        apply_trans_button.clicked.connect(self.apply_translation)
        trans_button_layout.addWidget(apply_trans_button)
        
        trans_layout.addLayout(trans_button_layout, 3, 0, 1, 2)
        
        layout.addLayout(trans_layout)
        
        layout.addSpacing(10)
        
        # Rotation controls
        rot_group = QLabel("Rotation (degrees):")
        rot_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(rot_group)
        
        rot_layout = QGridLayout()
        self.x_rot_input = QLineEdit("0.0")
        self.y_rot_input = QLineEdit("0.0")
        self.z_rot_input = QLineEdit("0.0")
        
        # Enterキーでapply_rotationを実行
        self.x_rot_input.returnPressed.connect(self.apply_rotation)
        self.y_rot_input.returnPressed.connect(self.apply_rotation)
        self.z_rot_input.returnPressed.connect(self.apply_rotation)
        
        rot_layout.addWidget(QLabel("Around X:"), 0, 0)
        rot_layout.addWidget(self.x_rot_input, 0, 1)
        rot_layout.addWidget(QLabel("Around Y:"), 1, 0)
        rot_layout.addWidget(self.y_rot_input, 1, 1)
        rot_layout.addWidget(QLabel("Around Z:"), 2, 0)
        rot_layout.addWidget(self.z_rot_input, 2, 1)
        
        rot_button_layout = QHBoxLayout()
        reset_rot_button = QPushButton("Reset")
        reset_rot_button.clicked.connect(self.reset_rotation_inputs)
        rot_button_layout.addWidget(reset_rot_button)
        
        apply_rot_button = QPushButton("Apply Rotation")
        apply_rot_button.clicked.connect(self.apply_rotation)
        rot_button_layout.addWidget(apply_rot_button)
        
        rot_layout.addLayout(rot_button_layout, 3, 0, 1, 2)
        
        layout.addLayout(rot_layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Enable picking to handle atom selection
        self.enable_picking()
    
    def eventFilter(self, obj, event):
        """3Dビューでのマウスイベント処理 - グループが選択されている場合はCustomInteractorStyleに任せる"""
        if obj == self.main_window.plotter.interactor:
            # ダブルクリック/トリプルクリックで状態が混乱するのを防ぐ
            if event.type() == QEvent.Type.MouseButtonDblClick:
                # ダブルクリックは無視し、状態をリセット
                self.is_dragging_group = False
                self.drag_start_pos = None
                self.mouse_moved_during_drag = False
                self.potential_drag = False
                if hasattr(self, 'clicked_atom_for_toggle'):
                    delattr(self, 'clicked_atom_for_toggle')
                return False
            
            if event.type() == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
                # 前回の状態をクリーンアップ（トリプルクリック対策）
                self.is_dragging_group = False
                self.potential_drag = False
                if hasattr(self, 'clicked_atom_for_toggle'):
                    delattr(self, 'clicked_atom_for_toggle')
                # グループが既に選択されている場合は、CustomInteractorStyleに処理を任せる
                if self.group_atoms:
                    return False
                
                # マウスプレス時の処理
                # マウスプレス時の処理
                try:
                    interactor = self.main_window.plotter.interactor
                    click_pos = interactor.GetEventPosition()
                    
                    # まずピッキングしてどの原子がクリックされたか確認
                    picker = self.main_window.plotter.picker
                    picker.Pick(click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer)
                    
                    clicked_atom_idx = None
                    if picker.GetActor() is self.main_window.atom_actor:
                        picked_position = np.array(picker.GetPickPosition())
                        distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                        closest_atom_idx = np.argmin(distances)
                        
                        # 閾値チェック
                        if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                            atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                            if atom:
                                atomic_num = atom.GetAtomicNum()
                                pt = Chem.GetPeriodicTable()
                                vdw_radius = pt.GetRvdw(atomic_num)
                                click_threshold = vdw_radius * 1.5
                                
                                if distances[closest_atom_idx] < click_threshold:
                                    clicked_atom_idx = int(closest_atom_idx)
                    
                    
                    # クリックされた原子の処理
                    if clicked_atom_idx is not None:
                        if self.group_atoms and clicked_atom_idx in self.group_atoms:
                            # 既存のグループ内の原子 - ドラッグ準備（まだドラッグとは確定しない）
                            self.is_dragging_group = False  # まだドラッグ中ではない
                            self.drag_start_pos = click_pos
                            self.drag_atom_idx = clicked_atom_idx
                            self.mouse_moved_during_drag = False
                            self.potential_drag = True  # ドラッグの可能性がある
                            self.clicked_atom_for_toggle = clicked_atom_idx  # トグル用に保存
                            # イベントを消費せず、カメラ操作を許可（閾値超えたらドラッグ開始）
                            return False
                        else:
                            # グループ外の原子 - 新しいグループを選択
                            # 親クラス（Mixin）のon_atom_pickedを手動で呼ぶ
                            self.on_atom_picked(clicked_atom_idx)
                            return True
                    else:
                        # 原子以外をクリック
                        # グループがあっても通常のカメラ操作を許可
                        return False
                    
                except Exception as e:
                    print(f"Error in mouse press: {e}")
                    return False
            
            elif event.type() == QEvent.Type.MouseMove:
                # マウス移動時の処理
                if getattr(self, 'potential_drag', False) and self.drag_start_pos and not self.is_dragging_group:
                    # potential_drag状態：閾値チェック
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]
                        
                        # 閾値を超えたらドラッグ開始
                        drag_threshold = 5  # ピクセル
                        if abs(dx) > drag_threshold or abs(dy) > drag_threshold:
                            # ドラッグ開始を確定
                            self.is_dragging_group = True
                            self.potential_drag = False
                            try:
                                self.main_window.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                            except Exception:
                                pass
                    except Exception:
                        pass
                    
                    # 閾値以下の場合はカメラ操作を許可
                    if not self.is_dragging_group:
                        return False
                
                if self.is_dragging_group and self.drag_start_pos:
                    # ドラッグモード中 - 移動距離を記録するのみ（リアルタイム更新なし）
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]
                        
                        if abs(dx) > 2 or abs(dy) > 2:
                            self.mouse_moved_during_drag = True
                    except Exception:
                        pass
                    
                    # ドラッグ中はイベントを消費してカメラ回転を防ぐ
                    return True
                
                # ホバー処理（ドラッグ中でない場合）
                if self.group_atoms:
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        picker = self.main_window.plotter.picker
                        picker.Pick(current_pos[0], current_pos[1], 0, self.main_window.plotter.renderer)
                        
                        if picker.GetActor() is self.main_window.atom_actor:
                            picked_position = np.array(picker.GetPickPosition())
                            distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                            closest_atom_idx = np.argmin(distances)
                            
                            if closest_atom_idx in self.group_atoms:
                                self.main_window.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                            else:
                                self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                        else:
                            self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                    except Exception:
                        pass
                
                # ドラッグ中でない場合はカメラ回転を許可
                return False
            
            elif event.type() == QEvent.Type.MouseButtonRelease and event.button() == Qt.MouseButton.LeftButton:
                # マウスリリース時の処理
                if getattr(self, 'potential_drag', False) or (self.is_dragging_group and self.drag_start_pos):
                    try:
                        if self.is_dragging_group and self.mouse_moved_during_drag:
                            # ドラッグが実行された - CustomInteractorStyleに任せる（何もしない）
                            pass
                        else:
                            # マウスが閾値以下の移動 = 単なるクリック
                            # グループ内の原子をクリックした場合は選択/解除をトグル
                            if hasattr(self, 'clicked_atom_for_toggle'):
                                clicked_atom = self.clicked_atom_for_toggle
                                delattr(self, 'clicked_atom_for_toggle')
                                # ドラッグ状態をリセットしてからトグル処理
                                self.is_dragging_group = False
                                self.drag_start_pos = None
                                self.mouse_moved_during_drag = False
                                self.potential_drag = False
                                if hasattr(self, 'last_drag_positions'):
                                    delattr(self, 'last_drag_positions')
                                # トグル処理を実行
                                self.on_atom_picked(clicked_atom)
                                try:
                                    self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                                except Exception:
                                    pass
                                return True
                        
                    except Exception:
                        pass
                    finally:
                        # ドラッグ状態をリセット
                        self.is_dragging_group = False
                        self.drag_start_pos = None
                        self.mouse_moved_during_drag = False
                        self.potential_drag = False
                        # 保存していた位置情報をクリア
                        if hasattr(self, 'last_drag_positions'):
                            delattr(self, 'last_drag_positions')
                        try:
                            self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                        except Exception:
                            pass
                    
                    return True  # イベントを消費
                
                # ドラッグ中でない場合は通常のリリース処理
                return False
        
        # その他のイベントは親クラスに渡す
        return False
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときに、その原子が属する連結成分全体を選択（複数グループ対応）"""
        # ドラッグ中は選択を変更しない（ただしリリース時のトグルは許可）
        if getattr(self, 'is_dragging_group', False):
            return
        
        # BFS/DFSで連結成分を探索
        visited = set()
        queue = [atom_idx]
        visited.add(atom_idx)
        
        while queue:
            current_idx = queue.pop(0)
            for bond_idx in range(self.mol.GetNumBonds()):
                bond = self.mol.GetBondWithIdx(bond_idx)
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                
                if begin_idx == current_idx and end_idx not in visited:
                    visited.add(end_idx)
                    queue.append(end_idx)
                elif end_idx == current_idx and begin_idx not in visited:
                    visited.add(begin_idx)
                    queue.append(begin_idx)
        
        # 新しいグループとして追加または解除
        if visited.issubset(self.group_atoms):
            # すでに選択されている - 解除
            self.group_atoms -= visited
        else:
            # 新しいグループを追加
            self.group_atoms |= visited
        
        self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()
    
    def update_display(self):
        if not self.group_atoms:
            self.selection_label.setText("No group selected")
        else:
            atom_info = []
            for atom_idx in sorted(self.group_atoms):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_info.append(f"{symbol}({atom_idx})")
            
            self.selection_label.setText(f"Selected group: {len(self.group_atoms)} atoms - {', '.join(atom_info[:5])}{' ...' if len(atom_info) > 5 else ''}")
    
    def show_atom_labels(self):
        """選択されたグループの原子をハイライト表示（Ctrlクリックと同じスタイル）"""
        self.clear_atom_labels()
        
        if not self.group_atoms:
            return
        
        # 選択された原子のインデックスリストを作成
        selected_indices = list(self.group_atoms)
        
        # 選択された原子の位置を取得
        selected_positions = self.main_window.atom_positions_3d[selected_indices]
        
        # 原子の半径を少し大きくしてハイライト表示
        selected_radii = np.array([VDW_RADII.get(
            self.mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3 
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
        
        # ハイライトアクターを追加して保存（ピッキング不可に設定）
        self.highlight_actor = self.main_window.plotter.add_mesh(
            highlight_glyphs, 
            color='yellow', 
            opacity=0.3, 
            name='move_group_highlight',
            pickable=False  # ピッキングを無効化
        )
        
        self.main_window.plotter.render()
    
    def clear_atom_labels(self):
        """原子ハイライトをクリア"""
        try:
            self.main_window.plotter.remove_actor('move_group_highlight')
        except Exception:
            pass
        
        if hasattr(self, 'highlight_actor'):
            try:
                self.main_window.plotter.remove_actor(self.highlight_actor)
            except Exception:
                pass
            self.highlight_actor = None
        
        try:
            self.main_window.plotter.render()
        except Exception:
            pass
    
    def reset_translation_inputs(self):
        """Translation入力フィールドをリセット"""
        self.x_trans_input.setText("0.0")
        self.y_trans_input.setText("0.0")
        self.z_trans_input.setText("0.0")
    
    def apply_translation(self):
        """選択したグループを並行移動"""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return
        
        try:
            dx = float(self.x_trans_input.text())
            dy = float(self.y_trans_input.text())
            dz = float(self.z_trans_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid translation values.")
            return
        
        translation_vector = np.array([dx, dy, dz])
        
        conf = self.mol.GetConformer()
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            new_pos = atom_pos + translation_vector
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # ラベルを再描画
        self.main_window.push_undo_state()
    
    def reset_rotation_inputs(self):
        """Rotation入力フィールドをリセット"""
        self.x_rot_input.setText("0.0")
        self.y_rot_input.setText("0.0")
        self.z_rot_input.setText("0.0")
    
    def apply_rotation(self):
        """選択したグループを回転"""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return
        
        try:
            rx = float(self.x_rot_input.text())
            ry = float(self.y_rot_input.text())
            rz = float(self.z_rot_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid rotation values.")
            return
        
        # 度をラジアンに変換
        rx_rad = np.radians(rx)
        ry_rad = np.radians(ry)
        rz_rad = np.radians(rz)
        
        # グループの重心を計算
        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.group_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        centroid = np.mean(positions, axis=0)
        
        # 回転行列を作成
        # X軸周り
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx_rad), -np.sin(rx_rad)],
            [0, np.sin(rx_rad), np.cos(rx_rad)]
        ])
        # Y軸周り
        Ry = np.array([
            [np.cos(ry_rad), 0, np.sin(ry_rad)],
            [0, 1, 0],
            [-np.sin(ry_rad), 0, np.cos(ry_rad)]
        ])
        # Z軸周り
        Rz = np.array([
            [np.cos(rz_rad), -np.sin(rz_rad), 0],
            [np.sin(rz_rad), np.cos(rz_rad), 0],
            [0, 0, 1]
        ])
        
        # 合成回転行列 (Z * Y * X)
        R = Rz @ Ry @ Rx
        
        # 各原子を回転
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            # 重心を原点に移動
            centered_pos = atom_pos - centroid
            # 回転
            rotated_pos = R @ centered_pos
            # 重心を元に戻す
            new_pos = rotated_pos + centroid
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # ラベルを再描画
        self.main_window.push_undo_state()
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.group_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        # ドラッグ関連のフラグもリセット
        self.is_dragging_group = False
        self.drag_start_pos = None
        if hasattr(self, 'last_drag_positions'):
            delattr(self, 'last_drag_positions')
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().reject()


class AlignPlaneDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, plane, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.plane = plane
        self.selected_atoms = set()
        
        # 事前選択された原子を追加
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)
        
        self.init_ui()
        
        # 事前選択された原子にラベルを追加
        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()
    
    def init_ui(self):
        plane_names = {'xy': 'XY', 'xz': 'XZ', 'yz': 'YZ'}
        self.setWindowTitle(f"Align to {plane_names[self.plane]} Plane")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(f"Click atoms in the 3D view to select them for align to the {plane_names[self.plane]} plane. At least 3 atoms are required.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for alignment")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply align")
        self.apply_button.clicked.connect(self.apply_PlaneAlign)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, 'picking_enabled') and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_PlaneAlign()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
    
    def select_all_atoms(self):
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, 'mol') and self.mol is not None:
                try:
                    n = self.mol.GetNumAtoms()
                    # create a set of indices [0..n-1]
                    self.selected_atoms = set(range(n))
                except Exception:
                    # fallback to main_window data map
                    self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()
            else:
                # fallback to main_window data map
                self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")
    
    def update_display(self):
        """表示を更新"""
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select for align (minimum 3 required)")
            self.apply_button.setEnabled(False)
        else:
            atom_list = sorted(self.selected_atoms)
            atom_display = []
            for i, atom_idx in enumerate(atom_list):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_display.append(f"#{i+1}: {symbol}({atom_idx})")
            
            self.selection_label.setText(f"Selected {count} atoms: {', '.join(atom_display)}")
            self.apply_button.setEnabled(count >= 3)
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
            
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            
            for i, atom_idx in enumerate(sorted_atoms):
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"#{i+1}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color='blue',
                    always_visible=True
                )
                self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def apply_PlaneAlign(self):
        """alignを適用（回転ベース）"""
        if len(self.selected_atoms) < 3:
            QMessageBox.warning(self, "Warning", "Please select at least 3 atoms for align.")
            return
        try:

            # 選択された原子の位置を取得
            selected_indices = list(self.selected_atoms)
            selected_positions = self.main_window.atom_positions_3d[selected_indices].copy()

            # 重心を計算
            centroid = np.mean(selected_positions, axis=0)

            # 重心を原点に移動
            centered_positions = selected_positions - centroid

            # 主成分分析で最適な平面を見つける
            # 選択された原子の座標の共分散行列を計算
            cov_matrix = np.cov(centered_positions.T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

            # 固有値が最も小さい固有ベクトルが平面の法線方向
            normal_vector = eigenvectors[:, 0]  # 最小固有値に対応する固有ベクトル

            # 目標の平面の法線ベクトルを定義
            if self.plane == 'xy':
                target_normal = np.array([0, 0, 1])  # Z軸方向
            elif self.plane == 'xz':
                target_normal = np.array([0, 1, 0])  # Y軸方向
            elif self.plane == 'yz':
                target_normal = np.array([1, 0, 0])  # X軸方向

            # 法線ベクトルの向きを調整（内積が正になるように）
            if np.dot(normal_vector, target_normal) < 0:
                normal_vector = -normal_vector

            # 回転軸と回転角度を計算
            rotation_axis = np.cross(normal_vector, target_normal)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-10:  # 回転が必要な場合
                rotation_axis = rotation_axis / rotation_axis_norm
                cos_angle = np.dot(normal_vector, target_normal)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                rotation_angle = np.arccos(cos_angle)

                # Rodrigues回転公式を使用して全分子を回転
                def rodrigues_rotation(v, axis, angle):
                    cos_a = np.cos(angle)
                    sin_a = np.sin(angle)
                    return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)

                # 分子全体を回転させる
                conf = self.mol.GetConformer()
                for i in range(self.mol.GetNumAtoms()):
                    current_pos = np.array(conf.GetAtomPosition(i))
                    # 重心基準で回転
                    centered_pos = current_pos - centroid
                    rotated_pos = rodrigues_rotation(centered_pos, rotation_axis, rotation_angle)
                    new_pos = rotated_pos + centroid
                    conf.SetAtomPosition(i, new_pos.tolist())
                    self.main_window.atom_positions_3d[i] = new_pos

            # 3D表示を更新
            self.main_window.draw_molecule_3d(self.mol)

            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Undo状態を保存
            self.main_window.push_undo_state()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply align: {str(e)}")
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()


class PlanarizeDialog(Dialog3DPickingMixin, QDialog):

    """選択原子群を最適フィット平面へ投影して planarize するダイアログ
    AlignPlane を参考にした選択UIを持ち、Apply ボタンで選択原子を平面へ直交射影する。
    """
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()

        if preselected_atoms:
            # 事前選択された原子を追加
            self.selected_atoms.update(preselected_atoms)

        self.init_ui()

        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    def init_ui(self):
        self.setWindowTitle("Planarize")
        self.setModal(False)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        layout = QVBoxLayout(self)

        instruction_label = QLabel("Click atoms in the 3D view to select them for planarization (minimum 3 required).")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)


        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
    
        # Select All Atoms ボタンを追加
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for planarization")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)

        self.apply_button = QPushButton("Apply planarize")
        self.apply_button.clicked.connect(self.apply_planarize)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
    
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        button_layout.addStretch()
    
        layout.addLayout(button_layout)

        # enable picking
        self.picker_connection = None
        self.enable_picking()

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def clear_selection(self):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def update_display(self):
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select for planarize (minimum 3 required)")
            self.apply_button.setEnabled(False)
        else:
            atom_list = sorted(self.selected_atoms)
            atom_display = []
            for i, atom_idx in enumerate(atom_list):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_display.append(f"#{i+1}: {symbol}({atom_idx})")
            self.selection_label.setText(f"Selected {count} atoms: {', '.join(atom_display)}")
            self.apply_button.setEnabled(count >= 3)

    def select_all_atoms(self):
        """Select all atoms in the current molecule (or fallback) and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, 'mol') and self.mol is not None:
                try:
                    n = self.mol.GetNumAtoms()
                    # create a set of indices [0..n-1]
                    self.selected_atoms = set(range(n))
                except Exception:
                    # fallback to main_window data map
                    self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()
            else:
                # fallback to main_window data map
                self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")

    def show_atom_labels(self):
        self.clear_atom_labels()
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        if self.selected_atoms:
            for i, atom_idx in enumerate(sorted(self.selected_atoms)):
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"#{i+1}"
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text],
                    point_size=20,
                    font_size=12,
                    text_color='cyan',
                    always_visible=True
                )
                self.selection_labels.append(label_actor)

    def clear_atom_labels(self):
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []

    def apply_planarize(self):
        if not self.selected_atoms or len(self.selected_atoms) < 3:
            QMessageBox.warning(self, "Warning", "Please select at least 3 atoms for planarize.")
            return

        try:
            selected_indices = list(sorted(self.selected_atoms))
            selected_positions = self.main_window.atom_positions_3d[selected_indices].copy()

            centroid = np.mean(selected_positions, axis=0)
            centered_positions = selected_positions - centroid

            # SVDによる最小二乗平面の法線取得
            u, s, vh = np.linalg.svd(centered_positions, full_matrices=False)
            normal = vh[-1]
            norm = np.linalg.norm(normal)
            if norm == 0:
                QMessageBox.warning(self, "Warning", "Cannot determine fit plane (degenerate positions).")
                return
            normal = normal / norm

            # 各点を重心を通る平面へ直交射影
            projections = centered_positions - np.outer(np.dot(centered_positions, normal), normal)
            new_positions = projections + centroid

            # 分子座標を更新
            conf = self.mol.GetConformer()
            for i, new_pos in zip(selected_indices, new_positions):
                conf.SetAtomPosition(int(i), new_pos.tolist())
                self.main_window.atom_positions_3d[int(i)] = new_pos

            # 3Dビュー更新
            self.main_window.draw_molecule_3d(self.mol)
            self.main_window.update_chiral_labels()
            self.main_window.push_undo_state()

            QMessageBox.information(self, "Success", f"Planarized {len(selected_indices)} atoms to best-fit plane.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to planarize: {e}")

    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()

    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()



class AlignmentDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, axis, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.axis = axis
        self.selected_atoms = set()
        
        # 事前選択された原子を追加（最大2個まで）
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms[:2])
        
        self.init_ui()
        
        # 事前選択された原子にラベルを追加
        if self.selected_atoms:
            for i, atom_idx in enumerate(sorted(self.selected_atoms), 1):
                self.add_selection_label(atom_idx, f"Atom {i}")
            self.update_display()
    
    def init_ui(self):
        axis_names = {'x': 'X-axis', 'y': 'Y-axis', 'z': 'Z-axis'}
        self.setWindowTitle(f"Align to {axis_names[self.axis]}")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(f"Click atoms in the 3D view to select them for alignment to the {axis_names[self.axis]}. Exactly 2 atoms are required. The first atom will be moved to the origin, and the second atom will be positioned on the {axis_names[self.axis]}.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Alignment")
        self.apply_button.clicked.connect(self.apply_alignment)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        # Dialog3DPickingMixinの機能を使用
        super().enable_picking()
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        # Dialog3DPickingMixinの機能を使用
        super().disable_picking()
    
    def on_atom_picked(self, atom_idx):
        """原子がクリックされた時の処理"""
        if self.main_window.current_mol is None:
            return
            
        if atom_idx in self.selected_atoms:
            # 既に選択されている場合は選択解除
            self.selected_atoms.remove(atom_idx)
            self.remove_atom_label(atom_idx)
        else:
            # 2つまでしか選択できない
            if len(self.selected_atoms) < 2:
                self.selected_atoms.add(atom_idx)
                # ラベルの順番を示す
                label_text = f"Atom {len(self.selected_atoms)}"
                self.add_selection_label(atom_idx, label_text)
        
        self.update_display()
    
    def update_display(self):
        """選択状態の表示を更新"""
        if len(self.selected_atoms) == 0:
            self.selection_label.setText("Click atoms to select for alignment (exactly 2 required)")
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 1:
            selected_list = list(self.selected_atoms)
            atom = self.mol.GetAtomWithIdx(selected_list[0])
            self.selection_label.setText(f"Selected 1 atom: {atom.GetSymbol()}{selected_list[0]+1}")
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 2:
            selected_list = sorted(list(self.selected_atoms))
            atom1 = self.mol.GetAtomWithIdx(selected_list[0])
            atom2 = self.mol.GetAtomWithIdx(selected_list[1])
            self.selection_label.setText(f"Selected 2 atoms: {atom1.GetSymbol()}{selected_list[0]+1}, {atom2.GetSymbol()}{selected_list[1]+1}")
            self.apply_button.setEnabled(True)
    
    def clear_selection(self):
        """選択をクリア"""
        self.clear_selection_labels()
        self.selected_atoms.clear()
        self.update_display()
    
    def add_selection_label(self, atom_idx, label_text):
        """選択された原子にラベルを追加"""
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        # 原子の位置を取得
        pos = self.main_window.atom_positions_3d[atom_idx]
        
        # ラベルを追加
        label_actor = self.main_window.plotter.add_point_labels(
            [pos], [label_text], 
            point_size=20, 
            font_size=12,
            text_color='yellow',
            always_visible=True
        )
        self.selection_labels.append(label_actor)
    
    def remove_atom_label(self, atom_idx):
        """特定の原子のラベルを削除"""
        # 簡単化のため、全ラベルをクリアして再描画
        self.clear_selection_labels()
        for i, idx in enumerate(sorted(self.selected_atoms), 1):
            if idx != atom_idx:
                self.add_selection_label(idx, f"Atom {i}")
    
    def clear_selection_labels(self):
        """選択ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def apply_alignment(self):
        """アライメントを適用"""
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Warning", "Please select exactly 2 atoms for alignment.")
            return
        try:

            selected_list = sorted(list(self.selected_atoms))
            atom1_idx, atom2_idx = selected_list[0], selected_list[1]

            conf = self.mol.GetConformer()

            # 原子の現在位置を取得
            pos1 = np.array(conf.GetAtomPosition(atom1_idx))
            pos2 = np.array(conf.GetAtomPosition(atom2_idx))

            # 最初に全分子を移動して、atom1を原点に配置
            translation = -pos1
            for i in range(self.mol.GetNumAtoms()):
                current_pos = np.array(conf.GetAtomPosition(i))
                new_pos = current_pos + translation
                conf.SetAtomPosition(i, new_pos.tolist())

            # atom2の新しい位置を取得（移動後）
            pos2_translated = pos2 + translation

            # atom2を選択した軸上に配置するための回転を計算
            axis_vectors = {
                'x': np.array([1.0, 0.0, 0.0]),
                'y': np.array([0.0, 1.0, 0.0]),
                'z': np.array([0.0, 0.0, 1.0])
            }
            target_axis = axis_vectors[self.axis]
            
            # atom2から原点への方向ベクトル
            current_vector = pos2_translated
            current_length = np.linalg.norm(current_vector)
            
            if current_length > 1e-10:  # ゼロベクトルでない場合
                current_vector_normalized = current_vector / current_length
                
                # 回転軸と角度を計算
                rotation_axis = np.cross(current_vector_normalized, target_axis)
                rotation_axis_length = np.linalg.norm(rotation_axis)
                
                if rotation_axis_length > 1e-10:  # 回転が必要
                    rotation_axis = rotation_axis / rotation_axis_length
                    cos_angle = np.dot(current_vector_normalized, target_axis)
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    rotation_angle = np.arccos(cos_angle)
                    
                    # ロドリゲスの回転公式を使用
                    def rodrigues_rotation(v, k, theta):
                        cos_theta = np.cos(theta)
                        sin_theta = np.sin(theta)
                        return (v * cos_theta + 
                               np.cross(k, v) * sin_theta + 
                               k * np.dot(k, v) * (1 - cos_theta))
                    
                    # 全ての原子に回転を適用
                    for i in range(self.mol.GetNumAtoms()):
                        current_pos = np.array(conf.GetAtomPosition(i))
                        rotated_pos = rodrigues_rotation(current_pos, rotation_axis, rotation_angle)
                        conf.SetAtomPosition(i, rotated_pos.tolist())
            
            # 3D座標を更新
            self.main_window.atom_positions_3d = np.array([
                list(conf.GetAtomPosition(i)) for i in range(self.mol.GetNumAtoms())
            ])
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Undo状態を保存
            self.main_window.push_undo_state()
            
            QMessageBox.information(self, "Success", f"Alignment to {self.axis.upper()}-axis completed.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply alignment: {str(e)}")
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().accept()


def main():
    # --- Windows タスクバーアイコンのための追加処理 ---
    if sys.platform == 'win32':
        myappid = 'hyoko.moleditpy.1.0' # アプリケーション固有のID（任意）
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QApplication(sys.argv)
    file_path = sys.argv[1] if len(sys.argv) > 1 else None
    window = MainWindow(initial_file=file_path)
    window.show()
    sys.exit(app.exec())


# --- Data Model ---
class MolecularData:
    def __init__(self):
        self.atoms = {}
        self.bonds = {}
        self._next_atom_id = 0
        self.adjacency_list = {} 

    def add_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self._next_atom_id
        self.atoms[atom_id] = {'symbol': symbol, 'pos': pos, 'item': None, 'charge': charge, 'radical': radical}
        self.adjacency_list[atom_id] = [] 
        self._next_atom_id += 1
        return atom_id

    def add_bond(self, id1, id2, order=1, stereo=0):
        # 立体結合の場合、IDの順序は方向性を意味するため、ソートしない。
        # 非立体結合の場合は、キーを正規化するためにソートする。
        if stereo == 0:
            if id1 > id2: id1, id2 = id2, id1

        bond_data = {'order': order, 'stereo': stereo, 'item': None}
        
        # 逆方向のキーも考慮して、新規結合かどうかをチェック
        is_new_bond = (id1, id2) not in self.bonds and (id2, id1) not in self.bonds
        if is_new_bond:
            if id1 in self.adjacency_list and id2 in self.adjacency_list:
                self.adjacency_list[id1].append(id2)
                self.adjacency_list[id2].append(id1)

        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)].update(bond_data)
            return (id1, id2), 'updated'
        else:
            self.bonds[(id1, id2)] = bond_data
            return (id1, id2), 'created'

    def remove_atom(self, atom_id):
        if atom_id in self.atoms:
            try:
                # Safely get neighbors before deleting the atom's own entry
                neighbors = self.adjacency_list.get(atom_id, [])
                for neighbor_id in neighbors:
                    if neighbor_id in self.adjacency_list and atom_id in self.adjacency_list[neighbor_id]:
                        self.adjacency_list[neighbor_id].remove(atom_id)

                # Now, safely delete the atom's own entry from the adjacency list
                if atom_id in self.adjacency_list:
                    del self.adjacency_list[atom_id]

                del self.atoms[atom_id]
                
                # Remove bonds involving this atom
                bonds_to_remove = [key for key in self.bonds if atom_id in key]
                for key in bonds_to_remove:
                    del self.bonds[key]
                    
            except Exception as e:
                print(f"Error removing atom {atom_id}: {e}")
                
                traceback.print_exc()

    def remove_bond(self, id1, id2):
        try:
            # 方向性のある立体結合(順方向/逆方向)と、正規化された非立体結合のキーを探す
            key_to_remove = None
            if (id1, id2) in self.bonds:
                key_to_remove = (id1, id2)
            elif (id2, id1) in self.bonds:
                key_to_remove = (id2, id1)

            if key_to_remove:
                if id1 in self.adjacency_list and id2 in self.adjacency_list[id1]:
                    self.adjacency_list[id1].remove(id2)
                if id2 in self.adjacency_list and id1 in self.adjacency_list[id2]:
                    self.adjacency_list[id2].remove(id1)
                del self.bonds[key_to_remove]
                
        except Exception as e:
            print(f"Error removing bond {id1}-{id2}: {e}")
            
            traceback.print_exc()


    def to_rdkit_mol(self, use_2d_stereo=True):
        """
        use_2d_stereo: Trueなら2D座標からE/Zを推定（従来通り）。FalseならE/Zラベル優先、ラベルがない場合のみ2D座標推定。
        3D変換時はuse_2d_stereo=Falseで呼び出すこと。
        """
        if not self.atoms:
            return None
        mol = Chem.RWMol()

        # --- Step 1: atoms ---
        atom_id_to_idx_map = {}
        for atom_id, data in self.atoms.items():
            atom = Chem.Atom(data['symbol'])
            atom.SetFormalCharge(data.get('charge', 0))
            atom.SetNumRadicalElectrons(data.get('radical', 0))
            atom.SetIntProp("_original_atom_id", atom_id)
            idx = mol.AddAtom(atom)
            atom_id_to_idx_map[atom_id] = idx

        # --- Step 2: bonds & stereo info保存（ラベル情報はここで保持） ---
        bond_stereo_info = {}  # bond_idx -> {'type': int, 'atom_ids': (id1,id2), 'bond_data': bond_data}
        for (id1, id2), bond_data in self.bonds.items():
            if id1 not in atom_id_to_idx_map or id2 not in atom_id_to_idx_map:
                continue
            idx1, idx2 = atom_id_to_idx_map[id1], atom_id_to_idx_map[id2]

            order_val = float(bond_data['order'])
            order = {1.0: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
                     2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE}.get(order_val, Chem.BondType.SINGLE)

            bond_idx = mol.AddBond(idx1, idx2, order) - 1

            # stereoラベルがあれば、bond_idxに対して詳細を保持（あとで使う）
            if 'stereo' in bond_data and bond_data['stereo'] in [1, 2, 3, 4]:
                bond_stereo_info[bond_idx] = {
                    'type': int(bond_data['stereo']),
                    'atom_ids': (id1, id2),
                    'bond_data': bond_data
                }

        # --- Step 3: sanitize ---
        final_mol = mol.GetMol()
        try:
            Chem.SanitizeMol(final_mol)
        except Exception as e:
            return None

        # --- Step 4: add 2D conformer ---
        # Convert from scene pixels to angstroms when creating RDKit conformer.
        conf = Chem.Conformer(final_mol.GetNumAtoms())
        conf.Set3D(False)
        for atom_id, data in self.atoms.items():
            if atom_id in atom_id_to_idx_map:
                idx = atom_id_to_idx_map[atom_id]
                pos = data.get('pos')
                if pos:
                    ax = pos.x() * ANGSTROM_PER_PIXEL
                    ay = -pos.y() * ANGSTROM_PER_PIXEL  # Y座標を反転（画面座標系→化学座標系）
                    conf.SetAtomPosition(idx, (ax, ay, 0.0))
        final_mol.AddConformer(conf)

        # --- Step 5: E/Zラベル優先の立体設定 ---
        # まず、E/Zラベルがあるbondを記録
        ez_labeled_bonds = set()
        for bond_idx, info in bond_stereo_info.items():
            if info['type'] in [3, 4]:
                ez_labeled_bonds.add(bond_idx)

        # 2D座標からE/Zを推定するのは、use_2d_stereo=True かつE/Zラベルがないbondのみ
        if use_2d_stereo:
            Chem.SetDoubleBondNeighborDirections(final_mol, final_mol.GetConformer(0))
        else:
            # 3D変換時: E/Zラベルがある場合は座標ベースの推定を完全に無効化
            if ez_labeled_bonds:
                # E/Zラベルがある場合は、すべての結合のBondDirをクリアして座標ベースの推定を無効化
                for b in final_mol.GetBonds():
                    b.SetBondDir(Chem.BondDir.NONE)
            else:
                # E/Zラベルがない場合のみ座標ベースの推定を実行
                Chem.SetDoubleBondNeighborDirections(final_mol, final_mol.GetConformer(0))

        # ヘルパー: 重原子優先で近傍を選ぶ
        def pick_preferred_neighbor(atom, exclude_idx):
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                if nbr.GetAtomicNum() > 1:
                    return nbr.GetIdx()
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() != exclude_idx:
                    return nbr.GetIdx()
            return None

        # --- Step 6: ラベルベースで上書き（E/Z を最優先） ---
        for bond_idx, info in bond_stereo_info.items():
            stereo_type = info['type']
            bond = final_mol.GetBondWithIdx(bond_idx)

            # 単結合の wedge/dash ラベル（1/2）がある場合
            if stereo_type in [1, 2]:
                if stereo_type == 1:
                    bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif stereo_type == 2:
                    bond.SetBondDir(Chem.BondDir.BEGINDASH)
                continue

            # 二重結合の E/Z ラベル（3/4）
            if stereo_type in [3, 4]:
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue

                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()

                bond_data = info.get('bond_data', {}) or {}
                stereo_atoms_specified = bond_data.get('stereo_atoms')

                if stereo_atoms_specified:
                    try:
                        a1_id, a2_id = stereo_atoms_specified
                        neigh1_idx = atom_id_to_idx_map.get(a1_id)
                        neigh2_idx = atom_id_to_idx_map.get(a2_id)
                    except Exception:
                        neigh1_idx = None
                        neigh2_idx = None
                else:
                    neigh1_idx = pick_preferred_neighbor(final_mol.GetAtomWithIdx(begin_atom_idx), end_atom_idx)
                    neigh2_idx = pick_preferred_neighbor(final_mol.GetAtomWithIdx(end_atom_idx), begin_atom_idx)

                if neigh1_idx is None or neigh2_idx is None:
                    continue

                bond.SetStereoAtoms(neigh1_idx, neigh2_idx)
                if stereo_type == 3:
                    bond.SetStereo(Chem.BondStereo.STEREOZ)
                elif stereo_type == 4:
                    bond.SetStereo(Chem.BondStereo.STEREOE)

                # 座標ベースでつけられた隣接単結合の BondDir（wedge/dash）がラベルと矛盾する可能性があるので消す
                b1 = final_mol.GetBondBetweenAtoms(begin_atom_idx, neigh1_idx)
                b2 = final_mol.GetBondBetweenAtoms(end_atom_idx, neigh2_idx)
                if b1 is not None:
                    b1.SetBondDir(Chem.BondDir.NONE)
                if b2 is not None:
                    b2.SetBondDir(Chem.BondDir.NONE)

        # Step 7: 最終化（キャッシュ更新 + 立体割当の再実行）
        final_mol.UpdatePropertyCache(strict=False)
        
        # 3D変換時（use_2d_stereo=False）でE/Zラベルがある場合は、force=Trueで強制適用
        if not use_2d_stereo and ez_labeled_bonds:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=True)
        else:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=False)
        return final_mol

    def to_mol_block(self):
        try:
            mol = self.to_rdkit_mol()
            if mol:
                return Chem.MolToMolBlock(mol, includeStereo=True)
        except Exception:
            pass
        if not self.atoms: return None
        atom_map = {old_id: new_id for new_id, old_id in enumerate(self.atoms.keys())}
        num_atoms, num_bonds = len(self.atoms), len(self.bonds)
        mol_block = "\n  MoleditPy\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
        for old_id, atom in self.atoms.items():
            # Convert scene pixel coordinates to angstroms when emitting MOL block
            x_px = atom['item'].pos().x()
            y_px = -atom['item'].pos().y()
            x, y = x_px * ANGSTROM_PER_PIXEL, y_px * ANGSTROM_PER_PIXEL
            z, symbol = 0.0, atom['symbol']
            charge = atom.get('charge', 0)

            chg_code = 0
            if charge == 3: chg_code = 1
            elif charge == 2: chg_code = 2
            elif charge == 1: chg_code = 3
            elif charge == -1: chg_code = 5
            elif charge == -2: chg_code = 6
            elif charge == -3: chg_code = 7

            mol_block += f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0{chg_code:3d}  0  0  0  0  0  0  0\n"

        for (id1, id2), bond in self.bonds.items():
            idx1, idx2, order = atom_map[id1] + 1, atom_map[id2] + 1, bond['order']
            stereo_code = 0
            bond_stereo = bond.get('stereo', 0)
            if bond_stereo == 1:
                stereo_code = 1
            elif bond_stereo == 2:
                stereo_code = 6

            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}{stereo_code:3d}  0  0  0\n"
            
        mol_block += "M  END\n"
        return mol_block


class AtomItem(QGraphicsItem):
    def __init__(self, atom_id, symbol, pos, charge=0, radical=0):
        super().__init__()
        self.atom_id, self.symbol, self.charge, self.radical, self.bonds, self.chiral_label = atom_id, symbol, charge, radical, [], None
        self.setPos(pos)
        self.implicit_h_count = 0 
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(1); self.font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD); self.update_style()
        self.setAcceptHoverEvents(True)
        self.hovered = False
        self.has_problem = False 

    def boundingRect(self):
        # --- paint()メソッドと完全に同じロジックでテキストの位置とサイズを計算 ---
        font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD)
        fm = QFontMetricsF(font)

        hydrogen_part = ""
        if self.implicit_h_count > 0:
            is_skeletal_carbon = (self.symbol == 'C' and 
                                      self.charge == 0 and 
                                      self.radical == 0 and 
                                      len(self.bonds) > 0)
            if not is_skeletal_carbon:
                hydrogen_part = "H"
                if self.implicit_h_count > 1:
                    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                    hydrogen_part += str(self.implicit_h_count).translate(subscript_map)

        flip_text = False
        if hydrogen_part and self.bonds:
            my_pos_x = self.pos().x()
            total_dx = 0.0
            # Defensive: some bonds may have missing atom references (None) or C++ wrappers
            # that have been deleted. Iterate and accumulate only valid partner positions.
            for b in self.bonds:
                # partner is the atom at the other end of the bond
                partner = b.atom2 if b.atom1 is self else b.atom1
                try:
                    if partner is None:
                        continue
                    # If SIP reports the wrapper as deleted, skip it
                    if sip_isdeleted_safe(partner):
                        continue
                    partner_pos = partner.pos()
                    if partner_pos is None:
                        continue
                    total_dx += partner_pos.x() - my_pos_x
                except Exception:
                    # Skip any bond that raises while inspecting; keep UI tolerant
                    continue

            if total_dx > 0:
                flip_text = True
        
        if flip_text:
            display_text = hydrogen_part + self.symbol
        else:
            display_text = self.symbol + hydrogen_part

        text_rect = fm.boundingRect(display_text)
        text_rect.adjust(-2, -2, 2, 2)
        if hydrogen_part:
            symbol_rect = fm.boundingRect(self.symbol)
            if flip_text:
                offset_x = symbol_rect.width() // 2
                text_rect.moveTo(offset_x - text_rect.width(), -text_rect.height() / 2)
            else:
                offset_x = -symbol_rect.width() // 2
                text_rect.moveTo(offset_x, -text_rect.height() / 2)
        else:
            text_rect.moveCenter(QPointF(0, 0))

        # 1. paint()で描画される背景の矩形(bg_rect)を計算する
        bg_rect = text_rect.adjusted(-5, -8, 5, 8)
        
        # 2. このbg_rectを基準として全体の描画領域を構築する
        full_visual_rect = QRectF(bg_rect)

        # 電荷記号の領域を計算に含める
        if self.charge != 0:
            # Chemical convention: single charge as "+"/"-", multiple as "2+"/"2-"
            if self.charge == 1:
                charge_str = "+"
            elif self.charge == -1:
                charge_str = "-"
            else:
                sign = '+' if self.charge > 0 else '-'
                charge_str = f"{abs(self.charge)}{sign}"
            charge_font = QFont("Arial", 12, QFont.Weight.Bold)
            charge_fm = QFontMetricsF(charge_font)
            charge_rect = charge_fm.boundingRect(charge_str)

            if flip_text:
                charge_pos = QPointF(text_rect.left() - charge_rect.width() - 2, text_rect.top())
            else:
                charge_pos = QPointF(text_rect.right() + 2, text_rect.top())
            charge_rect.moveTopLeft(charge_pos)
            full_visual_rect = full_visual_rect.united(charge_rect)

        # ラジカル記号の領域を計算に含める
        if self.radical > 0:
            radical_area = QRectF(text_rect.center().x() - 8, text_rect.top() - 8, 16, 8)
            full_visual_rect = full_visual_rect.united(radical_area)

        # 3. 選択ハイライト等のための最終的なマージンを追加する
        return full_visual_rect.adjusted(-3, -3, 3, 3)

    def shape(self):
        scene = self.scene()
        if not scene or not scene.views():
            path = QPainterPath()
            hit_r = max(4.0, ATOM_RADIUS - 6.0) * 2
            path.addEllipse(QRectF(-hit_r, -hit_r, hit_r * 2.0, hit_r * 2.0))
            return path

        view = scene.views()[0]
        scale = view.transform().m11()

        scene_radius = DESIRED_ATOM_PIXEL_RADIUS / scale

        path = QPainterPath()
        path.addEllipse(QPointF(0, 0), scene_radius, scene_radius)
        return path

    def paint(self, painter, option, widget):
        color = CPK_COLORS.get(self.symbol, CPK_COLORS['DEFAULT'])
        if self.is_visible:
            # 1. 描画の準備
            painter.setFont(self.font)
            fm = painter.fontMetrics()

            # --- 水素部分のテキストを作成 ---
            hydrogen_part = ""
            if self.implicit_h_count > 0:
                is_skeletal_carbon = (self.symbol == 'C' and 
                                      self.charge == 0 and 
                                      self.radical == 0 and 
                                      len(self.bonds) > 0)
                if not is_skeletal_carbon:
                    hydrogen_part = "H"
                    if self.implicit_h_count > 1:
                        subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                        hydrogen_part += str(self.implicit_h_count).translate(subscript_map)

            # --- テキストを反転させるか決定 ---
            flip_text = False
            # 水素ラベルがあり、結合が1本以上ある場合のみ反転を考慮
            if hydrogen_part and self.bonds:

                    # 相対的なX座標で、結合が左右どちらに偏っているか判定
                    my_pos_x = self.pos().x()
                    total_dx = 0.0
                    # Defensive: some bonds may have missing atom references (None) or
                    # wrappers that were deleted by SIP. Only accumulate valid partner positions.
                    for bond in self.bonds:
                        try:
                            other_atom = bond.atom1 if bond.atom2 is self else bond.atom2
                            if other_atom is None:
                                continue
                            # If SIP reports the wrapper as deleted, skip it
                            try:
                                if sip_isdeleted_safe(other_atom):
                                    continue
                            except Exception:
                                # If sip check fails, continue defensively
                                pass

                            other_pos = None
                            try:
                                other_pos = other_atom.pos()
                            except Exception:
                                # Accessing .pos() may raise if the C++ object was destroyed
                                other_pos = None

                            if other_pos is None:
                                continue

                            total_dx += (other_pos.x() - my_pos_x)
                        except Exception:
                            # Skip any problematic bond/partner rather than crashing the paint
                            continue

                    # 結合が主に右側にある場合はテキストを反転させる
                    if total_dx > 0:
                        flip_text = True

            # --- 表示テキストとアライメントを最終決定 ---
            if flip_text:
                display_text = hydrogen_part + self.symbol
                alignment_flag = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            else:
                display_text = self.symbol + hydrogen_part
                alignment_flag = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter

            text_rect = fm.boundingRect(display_text)
            text_rect.adjust(-2, -2, 2, 2)
            symbol_rect = fm.boundingRect(self.symbol) # 主元素のみの幅を計算

            # --- テキストの描画位置を決定 ---
            # 水素ラベルがない場合 (従来通り中央揃え)
            if not hydrogen_part:
                alignment_flag = Qt.AlignmentFlag.AlignCenter
                text_rect.moveCenter(QPointF(0, 0).toPoint())
            # 水素ラベルがあり、反転する場合 (右揃え)
            elif flip_text:
                # 主元素の中心が原子の中心に来るように、矩形の右端を調整
                offset_x = symbol_rect.width() // 2
                text_rect.moveTo(offset_x - text_rect.width(), -text_rect.height() // 2)
            # 水素ラベルがあり、反転しない場合 (左揃え)
            else:
                # 主元素の中心が原子の中心に来るように、矩形の左端を調整
                offset_x = -symbol_rect.width() // 2
                text_rect.moveTo(offset_x, -text_rect.height() // 2)

            # 2. 原子記号の背景を白で塗りつぶす
            if self.scene():
                bg_brush = self.scene().backgroundBrush()
                bg_rect = text_rect.adjusted(-5, -8, 5, 8)
                painter.setBrush(bg_brush)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawEllipse(bg_rect)
            
            # 3. 原子記号自体を描画
            if self.symbol == 'H':
                painter.setPen(QPen(Qt.GlobalColor.black))
            else:
                painter.setPen(QPen(color))
            painter.drawText(text_rect, int(alignment_flag), display_text)
            
            # --- 電荷とラジカルの描画  ---
            if self.charge != 0:
                # Chemical convention: single charge as "+"/"-", multiple as "2+"/"2-"
                if self.charge == 1:
                    charge_str = "+"
                elif self.charge == -1:
                    charge_str = "-"
                else:
                    sign = '+' if self.charge > 0 else '-'
                    charge_str = f"{abs(self.charge)}{sign}"
                charge_font = QFont("Arial", 12, QFont.Weight.Bold)
                painter.setFont(charge_font)
                charge_rect = painter.fontMetrics().boundingRect(charge_str)
                # 電荷の位置も反転に対応
                if flip_text:
                    charge_pos = QPointF(text_rect.left() - charge_rect.width() -2, text_rect.top() + charge_rect.height() - 2)
                else:
                    charge_pos = QPointF(text_rect.right() + 2, text_rect.top() + charge_rect.height() - 2)
                painter.setPen(Qt.GlobalColor.black)
                painter.drawText(charge_pos, charge_str)
            
            if self.radical > 0:
                painter.setBrush(QBrush(Qt.GlobalColor.black))
                painter.setPen(Qt.PenStyle.NoPen)
                radical_pos_y = text_rect.top() - 5
                if self.radical == 1:
                    painter.drawEllipse(QPointF(text_rect.center().x(), radical_pos_y), 3, 3)
                elif self.radical == 2:
                    painter.drawEllipse(QPointF(text_rect.center().x() - 5, radical_pos_y), 3, 3)
                    painter.drawEllipse(QPointF(text_rect.center().x() + 5, radical_pos_y), 3, 3)


        # --- 選択時のハイライトなど ---
        if self.has_problem:
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(255, 0, 0, 200), 4))
            painter.drawRect(self.boundingRect())
        elif self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.boundingRect())
        if (not self.isSelected()) and getattr(self, 'hovered', False):
            pen = QPen(QColor(144, 238, 144, 200), 3)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(pen)
            painter.drawRect(self.boundingRect())

    def update_style(self):
        self.is_visible = not (self.symbol == 'C' and len(self.bonds) > 0 and self.charge == 0 and self.radical == 0)
        self.update()


    # 約203行目 AtomItem クラス内

    def itemChange(self, change, value):
        res = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            if self.flags() & QGraphicsItem.GraphicsItemFlag.ItemIsMovable:
                # Prevent cascading updates during batch operations
                if not getattr(self, '_updating_position', False):
                    for bond in self.bonds: 
                        if bond.scene():  # Only update if bond is still in scene
                            bond.update_position()
            
        return res

    def hoverEnterEvent(self, event):
        # シーンのモードにかかわらず、ホバー時にハイライトを有効にする
        self.hovered = True
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        if self.hovered:
            self.hovered = False
            self.update()
        super().hoverLeaveEvent(event)

class BondItem(QGraphicsItem):

    def get_ez_label_rect(self):
        """E/Zラベルの描画範囲（シーン座標）を返す。ラベルが無い場合はNone。"""
        if self.order != 2 or self.stereo not in [3, 4]:
            return None
        line = self.get_line_in_local_coords()
        center = line.center()
        label_width = EZ_LABEL_BOX_SIZE
        label_height = EZ_LABEL_BOX_SIZE
        label_rect = QRectF(center.x() - label_width/2, center.y() - label_height/2, label_width, label_height)
        # シーン座標に変換
        return self.mapToScene(label_rect).boundingRect()
    def set_stereo(self, new_stereo):
        try:
            # ラベルを消す場合は、消す前のboundingRectをscene().invalidateで強制的に無効化
            if new_stereo == 0 and self.stereo in [3, 4] and self.scene():
                rect = self.mapToScene(self.boundingRect()).boundingRect()
                self.scene().invalidate(rect, QGraphicsScene.SceneLayer.BackgroundLayer | QGraphicsScene.SceneLayer.ForegroundLayer)
            
            self.prepareGeometryChange()
            self.stereo = new_stereo
            self.update()
            
            if self.scene() and self.scene().views():
                try:
                    self.scene().views()[0].viewport().update()
                except (IndexError, RuntimeError):
                    # Handle case where views are being destroyed
                    pass
                    
        except Exception as e:
            print(f"Error in BondItem.set_stereo: {e}")
            # Continue without crashing
            self.stereo = new_stereo

    def set_order(self, new_order):
        self.prepareGeometryChange()
        self.order = new_order
        self.update()
        if self.scene() and self.scene().views():
            self.scene().views()[0].viewport().update()
    def __init__(self, atom1_item, atom2_item, order=1, stereo=0):
        super().__init__()
        # Validate input parameters
        if atom1_item is None or atom2_item is None:
            raise ValueError("BondItem requires non-None atom items")
        self.atom1, self.atom2, self.order, self.stereo = atom1_item, atom2_item, order, stereo
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        self.update_position()
        self.setAcceptHoverEvents(True)
        self.hovered = False


    def get_line_in_local_coords(self):
        if self.atom1 is None or self.atom2 is None:
            return QLineF(0, 0, 0, 0)
        try:
            p2 = self.mapFromItem(self.atom2, 0, 0)
            return QLineF(QPointF(0, 0), p2)
        except (RuntimeError, TypeError):
            # Handle case where atoms are deleted from scene
            return QLineF(0, 0, 0, 0)

    def boundingRect(self):
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            line = QLineF(0, 0, 0, 0)
        bond_offset = globals().get('BOND_OFFSET', 2)
        extra = (getattr(self, 'order', 1) - 1) * bond_offset + 20
        rect = QRectF(line.p1(), line.p2()).normalized().adjusted(-extra, -extra, extra, extra)

        # E/Zラベルの描画範囲も考慮して拡張（QFontMetricsFで正確に）
        if self.order == 2 and self.stereo in [3, 4]:
            font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD)
            font.setItalic(True)
            text = "Z" if self.stereo == 3 else "E"
            fm = QFontMetricsF(font)
            text_rect = fm.boundingRect(text)
            outline = EZ_LABEL_TEXT_OUTLINE  # 輪郭の太さ分
            margin = EZ_LABEL_MARGIN   # 追加余白
            center = line.center()
            label_rect = QRectF(center.x() - text_rect.width()/2 - outline - margin,
                                center.y() - text_rect.height()/2 - outline - margin,
                                text_rect.width() + 2*outline + 2*margin,
                                text_rect.height() + 2*outline + 2*margin)
            rect = rect.united(label_rect)
        return rect

    def shape(self):
        path = QPainterPath()
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            return path 
        if line.length() == 0:
            return path

        scene = self.scene()
        if not scene or not scene.views():
            return super().shape()

        view = scene.views()[0]
        scale = view.transform().m11()

        scene_width = DESIRED_BOND_PIXEL_WIDTH / scale

        stroker = QPainterPathStroker()
        stroker.setWidth(scene_width)
        stroker.setCapStyle(Qt.PenCapStyle.RoundCap)  
        stroker.setJoinStyle(Qt.PenJoinStyle.RoundJoin) 

        center_line_path = QPainterPath(line.p1())
        center_line_path.lineTo(line.p2())
        
        return stroker.createStroke(center_line_path)

    def paint(self, painter, option, widget):
        if self.atom1 is None or self.atom2 is None:
            return
        line = self.get_line_in_local_coords()
        if line.length() == 0: return

        # --- 1. 選択状態に応じてペンとブラシを準備 ---
        if self.isSelected():
            selection_color = QColor("blue")
            painter.setPen(QPen(selection_color, 3))
            painter.setBrush(QBrush(selection_color))
        else:
            # Allow bond color override from app settings (2D color)
            try:
                sc = self.scene()
                if sc is not None and hasattr(sc, 'window') and sc.window is not None:
                    bond_hex = sc.window.settings.get('bond_color', '#222222')
                    bond_color = QColor(bond_hex)
                    painter.setPen(QPen(bond_color, 2))
                else:
                    painter.setPen(self.pen)
            except Exception:
                painter.setPen(self.pen)
            painter.setBrush(QBrush(Qt.GlobalColor.black))

        # --- 立体化学 (Wedge/Dash) の描画 ---
        if self.order == 1 and self.stereo in [1, 2]:
            vec = line.unitVector()
            normal = vec.normalVector()
            p1 = line.p1() + vec.p2() * 5
            p2 = line.p2() - vec.p2() * 5

            if self.stereo == 1: # Wedge (くさび形)
                offset = QPointF(normal.dx(), normal.dy()) * 6.0
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)
            
            elif self.stereo == 2: # Dash (破線)
                painter.save()
                if not self.isSelected():
                    pen = painter.pen()
                    pen.setWidthF(2.5) 
                    painter.setPen(pen)
                
                num_dashes = 8
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = 12.0 * t
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.drawLine(start_pt - offset, start_pt + offset)
                painter.restore()
        
        # --- 通常の結合 (単/二重/三重) の描画 ---
        else:
            if self.order == 1:
                painter.drawLine(line)
            else:
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * BOND_OFFSET

                if self.order == 2:
                    # -------------------- ここから差し替え --------------------)
                    line1 = line.translated(offset)
                    line2 = line.translated(-offset)
                    painter.drawLine(line1)
                    painter.drawLine(line2)

                    # E/Z ラベルの描画処理
                    if self.stereo in [3, 4]:
                        painter.save() # 現在の描画設定を保存

                        # --- ラベルの設定 ---
                        font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD)
                        font.setItalic(True)
                        text_color = QColor("gray")
                        # 輪郭の色を背景色と同じにする（scene()がNoneのときは安全なフォールバックを使う）
                        outline_color = None
                        try:
                            sc = self.scene()
                            if sc is not None:
                                outline_color = sc.backgroundBrush().color()
                        except Exception:
                            outline_color = None
                        if outline_color is None:
                            # デフォルトでは白背景を想定して黒系の輪郭が見やすい
                            outline_color = QColor(255, 255, 255)

                        # --- 描画パスの作成 ---
                        text = "Z" if self.stereo == 3 else "E"
                        path = QPainterPath()
                        
                        # テキストが正確に中央に来るように位置を計算
                        fm = QFontMetricsF(font)
                        text_rect = fm.boundingRect(text)
                        text_rect.moveCenter(line.center())
                        path.addText(text_rect.topLeft(), font, text)

                        # --- 輪郭の描画 ---
                        stroker = QPainterPathStroker()
                        stroker.setWidth(EZ_LABEL_TEXT_OUTLINE) # 輪郭の太さ
                        outline_path = stroker.createStroke(path)
                        
                        painter.setBrush(outline_color)
                        painter.setPen(Qt.PenStyle.NoPen)
                        painter.drawPath(outline_path)

                        # --- 文字本体の描画 ---
                        painter.setBrush(text_color)
                        painter.setPen(text_color)
                        painter.drawPath(path)

                        painter.restore() # 描画設定を元に戻す

                elif self.order == 3:
                    painter.drawLine(line)
                    painter.drawLine(line.translated(offset))
                    painter.drawLine(line.translated(-offset))

        # --- 2. ホバー時のエフェクトを上から重ねて描画 ---
        if (not self.isSelected()) and getattr(self, 'hovered', False):
            try:
                # ホバー時のハイライトを太めの半透明な線で描画
                hover_pen = QPen(QColor(144, 238, 144, 180), HOVER_PEN_WIDTH) # LightGreen, 半透明
                hover_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                painter.setPen(hover_pen)
                painter.drawLine(line) 
            except Exception:
                pass



    def update_position(self):
        try:
            self.prepareGeometryChange()
            if self.atom1:
                self.setPos(self.atom1.pos())
            self.update()
        except Exception as e:
            print(f"Error updating bond position: {e}")
            # Continue without crashing


    def hoverEnterEvent(self, event):
        scene = self.scene()
        mode = getattr(scene, 'mode', '')
        self.hovered = True
        self.update()
        if self.scene():
            self.scene().set_hovered_item(self)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        if self.hovered:
            self.hovered = False
            self.update()
        if self.scene():
            self.scene().set_hovered_item(None)
        super().hoverLeaveEvent(event)


class TemplatePreviewItem(QGraphicsItem):
    def __init__(self):
        super().__init__()
        self.setZValue(2)
        self.pen = QPen(QColor(80, 80, 80, 180), 2)
        self.polygon = QPolygonF()
        self.is_aromatic = False
        self.user_template_points = []
        self.user_template_bonds = []
        self.user_template_atoms = []
        self.is_user_template = False

    def set_geometry(self, points, is_aromatic=False):
        self.prepareGeometryChange()
        self.polygon = QPolygonF(points)
        self.is_aromatic = is_aromatic
        self.is_user_template = False
        self.update()
    
    def set_user_template_geometry(self, points, bonds_info, atoms_data):
        self.prepareGeometryChange()
        self.user_template_points = points
        self.user_template_bonds = bonds_info
        self.user_template_atoms = atoms_data
        self.is_user_template = True
        self.is_aromatic = False
        self.polygon = QPolygonF()
        self.update()

    def boundingRect(self):
        if self.is_user_template and self.user_template_points:
            # Calculate bounding rect for user template
            min_x = min(p.x() for p in self.user_template_points)
            max_x = max(p.x() for p in self.user_template_points)
            min_y = min(p.y() for p in self.user_template_points)
            max_y = max(p.y() for p in self.user_template_points)
            return QRectF(min_x - 20, min_y - 20, max_x - min_x + 40, max_y - min_y + 40)
        return self.polygon.boundingRect().adjusted(-5, -5, 5, 5)

    def paint(self, painter, option, widget):
        if self.is_user_template:
            self.paint_user_template(painter)
        else:
            self.paint_regular_template(painter)
    
    def paint_regular_template(self, painter):
        painter.setPen(self.pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        if not self.polygon.isEmpty():
            painter.drawPolygon(self.polygon)
            if self.is_aromatic:
                center = self.polygon.boundingRect().center()
                radius = QLineF(center, self.polygon.first()).length() * 0.6
                painter.drawEllipse(center, radius, radius)
    
    def paint_user_template(self, painter):
        if not self.user_template_points:
            return
        
        # Draw bonds first with better visibility
        # Bond preview color may also follow the configured bond color
        try:
            bond_hex = self.scene().window.settings.get('bond_color', '#222222')
            bond_pen = QPen(QColor(bond_hex), 2.5)
        except Exception:
            bond_pen = QPen(QColor(100, 100, 100, 200), 2.5)
        painter.setPen(bond_pen)
        
        for bond_info in self.user_template_bonds:
            if len(bond_info) >= 3:
                atom1_idx, atom2_idx, order = bond_info[:3]
            else:
                atom1_idx, atom2_idx = bond_info[:2]
                order = 1
                
            if atom1_idx < len(self.user_template_points) and atom2_idx < len(self.user_template_points):
                pos1 = self.user_template_points[atom1_idx]
                pos2 = self.user_template_points[atom2_idx]
                
                if order == 2:
                    # Double bond - draw two parallel lines
                    line = QLineF(pos1, pos2)
                    normal = line.normalVector()
                    normal.setLength(4)
                    
                    line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                    line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                    
                    painter.drawLine(line1)
                    painter.drawLine(line2)
                elif order == 3:
                    # Triple bond - draw three parallel lines
                    line = QLineF(pos1, pos2)
                    normal = line.normalVector()
                    normal.setLength(6)
                    
                    painter.drawLine(line)
                    line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                    line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                    
                    painter.drawLine(line1)
                    painter.drawLine(line2)
                else:
                    # Single bond
                    painter.drawLine(QLineF(pos1, pos2))
        
        # Draw atoms - white ellipse background to hide bonds, then CPK colored text
        for i, pos in enumerate(self.user_template_points):
            if i < len(self.user_template_atoms):
                atom_data = self.user_template_atoms[i]
                symbol = atom_data.get('symbol', 'C')
                
                # Draw all non-carbon atoms including hydrogen with white background ellipse + CPK colored text
                if symbol != 'C':
                    # Get CPK color for text
                    color = CPK_COLORS.get(symbol, CPK_COLORS.get('DEFAULT', QColor('#FF1493')))
                    
                    # Draw white background ellipse to hide bonds
                    painter.setPen(QPen(Qt.GlobalColor.white, 0))  # No border
                    painter.setBrush(QBrush(Qt.GlobalColor.white))
                    painter.drawEllipse(int(pos.x() - 12), int(pos.y() - 8), 24, 16)
                    
                    # Draw CPK colored text on top
                    painter.setPen(QPen(color))
                    font = QFont("Arial", 12, QFont.Weight.Bold)  # Larger font
                    painter.setFont(font)
                    metrics = painter.fontMetrics()
                    text_rect = metrics.boundingRect(symbol)
                    text_pos = QPointF(pos.x() - text_rect.width()/2, pos.y() + text_rect.height()/3)
                    painter.drawText(text_pos, symbol)

class MoleculeScene(QGraphicsScene):
    def clear_template_preview(self):
        """テンプレートプレビュー用のゴースト線を全て消す"""
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(item, '_is_template_preview', False):
                try:
                    # If SIP reports the wrapper as deleted, skip it. Otherwise
                    # ensure it is still in a scene before attempting removal.
                    if sip_isdeleted_safe(item):
                        continue
                    sc = None
                    try:
                        sc = item.scene() if hasattr(item, 'scene') else None
                    except Exception:
                        sc = None
                    if sc is None:
                        continue
                    try:
                        self.removeItem(item)
                    except Exception:
                        # Best-effort: ignore removal errors to avoid crashes during teardown
                        pass
                except Exception:
                    # Non-fatal: continue with other items
                    continue
        self.template_context = {}
        if hasattr(self, 'template_preview'):
            self.template_preview.hide()

    def __init__(self, data, window):
        super().__init__()
        self.data, self.window = data, window
        self.mode, self.current_atom_symbol = 'select', 'C'
        self.bond_order, self.bond_stereo = 1, 0
        self.start_atom, self.temp_line, self.start_pos = None, None, None; self.press_pos = None
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.hovered_item = None
        
        self.key_to_symbol_map = {
            Qt.Key.Key_C: 'C', Qt.Key.Key_N: 'N', Qt.Key.Key_O: 'O', Qt.Key.Key_S: 'S',
            Qt.Key.Key_F: 'F', Qt.Key.Key_B: 'B', Qt.Key.Key_I: 'I', Qt.Key.Key_H: 'H',
            Qt.Key.Key_P: 'P',
        }
        self.key_to_symbol_map_shift = { Qt.Key.Key_C: 'Cl', Qt.Key.Key_B: 'Br', Qt.Key.Key_S: 'Si',}

        self.key_to_bond_mode_map = {
            Qt.Key.Key_1: 'bond_1_0',
            Qt.Key.Key_2: 'bond_2_0',
            Qt.Key.Key_3: 'bond_3_0',
            Qt.Key.Key_W: 'bond_1_1',
            Qt.Key.Key_D: 'bond_1_2',
        }
        self.reinitialize_items()

    def reinitialize_items(self):
        self.template_preview = TemplatePreviewItem(); self.addItem(self.template_preview)
        self.template_preview.hide(); self.template_preview_points = []; self.template_context = {}
        # Hold strong references to deleted wrappers for the lifetime of the scene
        # to avoid SIP/C++ finalization causing segfaults when Python still
        # briefly touches those objects elsewhere in the app. Items collected
        # here are hidden and never accessed again by normal code paths.
        self._deleted_items = []
        # Ensure we purge any held deleted-wrapper references when the
        # application is shutting down. Connecting here is safe even if
        # multiple scenes exist; the slot is defensive and idempotent.
        try:
            app = QApplication.instance()
            if app is not None:
                try:
                    app.aboutToQuit.connect(self.purge_deleted_items)
                except Exception:
                    # If connecting fails for any reason, continue without
                    # the connection — at worst holders will be freed by
                    # process teardown.
                    pass
        except Exception:
            pass

    def clear_all_problem_flags(self):
        """全ての AtomItem の has_problem フラグをリセットし、再描画する"""
        needs_update = False
        for atom_data in self.data.atoms.values():
            item = atom_data.get('item')
            # hasattr は安全性のためのチェック
            if item and hasattr(item, 'has_problem') and item.has_problem: 
                item.has_problem = False
                item.update()
                needs_update = True
        return needs_update

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        
        # 削除されたオブジェクトを安全にチェックして初期位置を記録
        self.initial_positions_in_event = {}
        for item in self.items():
            if isinstance(item, AtomItem):
                try:
                    self.initial_positions_in_event[item] = item.pos()
                except RuntimeError:
                    # オブジェクトが削除されている場合はスキップ
                    continue

        if not self.window.is_2d_editable:
            return

        if event.button() == Qt.MouseButton.RightButton:
            item = self.itemAt(event.scenePos(), self.views()[0].transform())
            if not isinstance(item, (AtomItem, BondItem)):
                return # 対象外のものをクリックした場合は何もしない
            data_changed = False
            # If the user has a rectangular multi-selection and the clicked item
            # is part of that selection, delete all selected items (atoms/bonds).
            try:
                selected_items = [it for it in self.selectedItems() if isinstance(it, (AtomItem, BondItem))]
            except Exception:
                selected_items = []

            if len(selected_items) > 1 and item in selected_items and not self.mode.startswith(('template', 'charge', 'radical')):
                # Delete the entire rectangular selection
                data_changed = self.delete_items(set(selected_items))
                if data_changed:
                    self.window.push_undo_state()
                self.press_pos = None
                event.accept()
                return
            # --- E/Zモード専用処理 ---
            if self.mode == 'bond_2_5':
                if isinstance(item, BondItem):
                    try:
                        # E/Zラベルを消す（ノーマルに戻す）
                        if item.stereo in [3, 4]:
                            item.set_stereo(0)
                            # データモデルも更新
                            for (id1, id2), bdata in self.data.bonds.items():
                                if bdata.get('item') is item:
                                    bdata['stereo'] = 0
                                    break
                            self.window.push_undo_state()
                            data_changed = False  # ここでundo済みなので以降で積まない
                    except Exception as e:
                        print(f"Error clearing E/Z label: {e}")
                        
                        traceback.print_exc()
                        if hasattr(self.window, 'statusBar'):
                            self.window.statusBar().showMessage(f"Error clearing E/Z label: {e}", 5000)
                # AtomItemは何もしない
            # --- 通常の処理 ---
            elif isinstance(item, AtomItem):
                # ラジカルモードの場合、ラジカルを0にする
                if self.mode == 'radical' and item.radical != 0:
                    item.prepareGeometryChange()
                    item.radical = 0
                    self.data.atoms[item.atom_id]['radical'] = 0
                    item.update_style()
                    data_changed = True
                # 電荷モードの場合、電荷を0にする
                elif self.mode in ['charge_plus', 'charge_minus'] and item.charge != 0:
                    item.prepareGeometryChange()
                    item.charge = 0
                    self.data.atoms[item.atom_id]['charge'] = 0
                    item.update_style()
                    data_changed = True
                # 上記以外のモード（テンプレート、電荷、ラジカルを除く）では原子を削除
                elif not self.mode.startswith(('template', 'charge', 'radical')):
                    data_changed = self.delete_items({item})
            elif isinstance(item, BondItem):
                # テンプレート、電荷、ラジカルモード以外で結合を削除
                if not self.mode.startswith(('template', 'charge', 'radical')):
                    data_changed = self.delete_items({item})

            if data_changed:
                self.window.push_undo_state()
            self.press_pos = None
            event.accept()
            return # 右クリック処理を完了し、左クリックの処理へ進ませない

        if self.mode.startswith('template'):
            self.clearSelection() # テンプレートモードでは選択処理を一切行わず、クリック位置の記録のみ行う
            return

        # Z,Eモードの時は選択処理を行わないようにする
        if self.mode in ['bond_2_5']:
            self.clearSelection()
            event.accept()
            return

        if getattr(self, "mode", "") != "select":
            self.clearSelection()
            event.accept()

        item = self.itemAt(self.press_pos, self.views()[0].transform())

        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != 'select':
                self.clearSelection()
                self.temp_line = QGraphicsLineItem(QLineF(self.start_atom.pos(), self.press_pos))
                self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
                self.addItem(self.temp_line)
            else:
                super().mousePressEvent(event)
        elif item is None and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos)); self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine)); self.addItem(self.temp_line)
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.window.is_2d_editable:
            return 

        if self.mode.startswith('template'):
            self.update_template_preview(event.scenePos())
        
        if not self.mouse_moved_since_press and self.press_pos:
            if (event.scenePos() - self.press_pos).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True
        
        if self.temp_line and not self.mode.startswith('template'):
            start_point = self.start_atom.pos() if self.start_atom else self.start_pos
            if not start_point:
                super().mouseMoveEvent(event)
                return

            current_pos = event.scenePos()
            end_point = current_pos

            target_atom = None
            for item in self.items(current_pos):
                if isinstance(item, AtomItem):
                    target_atom = item
                    break
            
            is_valid_snap_target = (
                target_atom is not None and
                (self.start_atom is None or target_atom is not self.start_atom)
            )

            if is_valid_snap_target:
                end_point = target_atom.pos()
            
            self.temp_line.setLine(QLineF(start_point, end_point))
        else: 
            # テンプレートモードであっても、ホバーイベントはここで伝播する
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if not self.window.is_2d_editable:
            return 

        end_pos = event.scenePos()
        is_click = self.press_pos and (end_pos - self.press_pos).manhattanLength() < QApplication.startDragDistance()

        if self.temp_line:
            try:
                if not sip_isdeleted_safe(self.temp_line):
                    try:
                        if getattr(self.temp_line, 'scene', None) and self.temp_line.scene():
                            self.removeItem(self.temp_line)
                    except Exception:
                        pass
            except Exception:
                try:
                    self.removeItem(self.temp_line)
                except Exception:
                    pass
            finally:
                self.temp_line = None

        if self.mode.startswith('template') and is_click:
            if self.template_context and self.template_context.get('points'):
                context = self.template_context
                # Check if this is a user template
                if self.mode.startswith('template_user'):
                    self.add_user_template_fragment(context)
                else:
                    self.add_molecule_fragment(context['points'], context['bonds_info'], existing_items=context.get('items', []))
                self.data_changed_in_event = True
                # イベント処理をここで完了させ、下のアイテムが選択されるのを防ぐ
                self.start_atom=None; self.start_pos = None; self.press_pos = None
                if self.data_changed_in_event: self.window.push_undo_state()
                return

        released_item = self.itemAt(end_pos, self.views()[0].transform())

        # 1. 特殊モード（ラジカル/電荷）の処理
        if (self.mode == 'radical') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            atom.prepareGeometryChange()
            # ラジカルの状態をトグル (0 -> 1 -> 2 -> 0)
            atom.radical = (atom.radical + 1) % 3 
            self.data.atoms[atom.atom_id]['radical'] = atom.radical
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom=None; self.start_pos = None; self.press_pos = None
            if self.data_changed_in_event: self.window.push_undo_state()
            return
        elif (self.mode == 'charge_plus' or self.mode == 'charge_minus') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            atom.prepareGeometryChange()
            delta = 1 if self.mode == 'charge_plus' else -1
            atom.charge += delta
            self.data.atoms[atom.atom_id]['charge'] = atom.charge
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom=None; self.start_pos = None; self.press_pos = None
            if self.data_changed_in_event: self.window.push_undo_state()
            return

        elif self.mode.startswith('bond') and is_click and isinstance(released_item, BondItem):
            b = released_item 
            if self.mode == 'bond_2_5':
                try:
                    if b.order == 2:
                        current_stereo = b.stereo
                        if current_stereo not in [3, 4]:
                            new_stereo = 3  # None -> Z
                        elif current_stereo == 3:
                            new_stereo = 4  # Z -> E
                        else:  # current_stereo == 4
                            new_stereo = 0  # E -> None
                        self.update_bond_stereo(b, new_stereo)
                        self.window.push_undo_state()  # ここでUndo stackに積む
                except Exception as e:
                    print(f"Error in E/Z stereo toggle: {e}")
                    
                    traceback.print_exc()
                    if hasattr(self.window, 'statusBar'):
                        self.window.statusBar().showMessage(f"Error changing E/Z stereochemistry: {e}", 5000)
                return # この後の処理は行わない
            elif self.bond_stereo != 0 and b.order == self.bond_order and b.stereo == self.bond_stereo:
                # 方向性を反転させる
                old_id1, old_id2 = b.atom1.atom_id, b.atom2.atom_id
                # 1. 古い方向の結合をデータから削除
                self.data.remove_bond(old_id1, old_id2)
                # 2. 逆方向で結合をデータに再追加
                new_key, _ = self.data.add_bond(old_id2, old_id1, self.bond_order, self.bond_stereo)
                # 3. BondItemの原子参照を入れ替え、新しいデータと関連付ける
                b.atom1, b.atom2 = b.atom2, b.atom1
                self.data.bonds[new_key]['item'] = b
                # 4. 見た目を更新
                b.update_position()
            else:
                # 既存の結合を一度削除
                self.data.remove_bond(b.atom1.atom_id, b.atom2.atom_id)
                # BondItemが記憶している方向(b.atom1 -> b.atom2)で、新しい結合様式を再作成
                # これにより、修正済みのadd_bondが呼ばれ、正しい方向で保存される
                new_key, _ = self.data.add_bond(b.atom1.atom_id, b.atom2.atom_id, self.bond_order, self.bond_stereo)
                # BondItemの見た目とデータ参照を更新
                b.prepareGeometryChange()
                b.order = self.bond_order
                b.stereo = self.bond_stereo
                self.data.bonds[new_key]['item'] = b
                b.update()
            self.clearSelection()
            self.data_changed_in_event = True
        # 3. 新規原子・結合の作成処理 (atom_* モード および すべての bond_* モードで許可)
        elif self.start_atom and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            line = QLineF(self.start_atom.pos(), end_pos); end_item = self.itemAt(end_pos, self.views()[0].transform())
            # 使用する結合様式を決定
            # atomモードの場合は bond_order/stereo を None にして create_bond にデフォルト値(1, 0)を適用
            # bond_* モードの場合は現在の設定 (self.bond_order/stereo) を使用
            order_to_use = self.bond_order if self.mode.startswith('bond') else None
            stereo_to_use = self.bond_stereo if self.mode.startswith('bond') else None
            if is_click:
                # 短いクリック: 既存原子のシンボル更新 (atomモードのみ)
                if self.mode.startswith('atom') and self.start_atom.symbol != self.current_atom_symbol:
                    self.start_atom.symbol=self.current_atom_symbol; self.data.atoms[self.start_atom.atom_id]['symbol']=self.current_atom_symbol; self.start_atom.update_style()
                    self.data_changed_in_event = True
            else:
                # ドラッグ: 新規結合または既存原子への結合
                if isinstance(end_item, AtomItem) and self.start_atom!=end_item: 
                    self.create_bond(self.start_atom, end_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                else:
                    new_id = self.create_atom(self.current_atom_symbol, end_pos); new_item = self.data.atoms[new_id]['item']
                    self.create_bond(self.start_atom, new_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                self.data_changed_in_event = True
        # 4. 空白領域からの新規作成処理 (atom_* モード および すべての bond_* モードで許可)
        elif self.start_pos and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            line = QLineF(self.start_pos, end_pos)
            # 使用する結合様式を決定
            order_to_use = self.bond_order if self.mode.startswith('bond') else None
            stereo_to_use = self.bond_stereo if self.mode.startswith('bond') else None
            if line.length() < 10:
                self.create_atom(self.current_atom_symbol, end_pos); self.data_changed_in_event = True
            else:
                end_item = self.itemAt(end_pos, self.views()[0].transform())
                if isinstance(end_item, AtomItem):
                    start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                    start_item = self.data.atoms[start_id]['item']
                    self.create_bond(start_item, end_item, bond_order=order_to_use, bond_stereo=stereo_to_use)
                else:
                    start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                    end_id = self.create_atom(self.current_atom_symbol, end_pos)
                    self.create_bond(
                        self.data.atoms[start_id]['item'], 
                        self.data.atoms[end_id]['item'], 
                        bond_order=order_to_use, 
                        bond_stereo=stereo_to_use
                    )
                self.data_changed_in_event = True 
        # 5. それ以外の処理 (Selectモードなど)
        else: super().mouseReleaseEvent(event)

        # 削除されたオブジェクトを安全にチェック
        moved_atoms = []
        for item, old_pos in self.initial_positions_in_event.items():
            try:
                # オブジェクトが有効で、シーンに存在し、位置が変更されているかチェック
                if item.scene() and item.pos() != old_pos:
                    moved_atoms.append(item)
            except RuntimeError:
                # オブジェクトが削除されている場合はスキップ
                continue
        if moved_atoms:
            self.data_changed_in_event = True
            bonds_to_update = set()
            for atom in moved_atoms:
                try:
                    self.data.atoms[atom.atom_id]['pos'] = atom.pos()
                    bonds_to_update.update(atom.bonds)
                except RuntimeError:
                    # オブジェクトが削除されている場合はスキップ
                    continue
            for bond in bonds_to_update: bond.update_position()
            # 原子移動後に測定ラベルの位置を更新
            self.window.update_2d_measurement_labels()
            if self.views(): self.views()[0].viewport().update()
        self.start_atom=None; self.start_pos = None; self.press_pos = None; self.temp_line = None
        self.template_context = {}
        # Clear user template data when switching modes
        if hasattr(self, 'user_template_data'):
            self.user_template_data = None
        if self.data_changed_in_event: self.window.push_undo_state()

    def mouseDoubleClickEvent(self, event):
        """ダブルクリックイベントを処理する"""
        item = self.itemAt(event.scenePos(), self.views()[0].transform())

        if self.mode in ['charge_plus', 'charge_minus', 'radical'] and isinstance(item, AtomItem):
            if self.mode == 'radical':
                item.prepareGeometryChange()
                item.radical = (item.radical + 1) % 3
                self.data.atoms[item.atom_id]['radical'] = item.radical
                item.update_style()
            else:
                item.prepareGeometryChange()
                delta = 1 if self.mode == 'charge_plus' else -1
                item.charge += delta
                self.data.atoms[item.atom_id]['charge'] = item.charge
                item.update_style()

            self.window.push_undo_state()

            event.accept()
            return
        
        # Select-mode: double-click should select the clicked atom/bond and
        # only the atoms/bonds connected to it (the connected component).
        if self.mode == 'select' and isinstance(item, (AtomItem, BondItem)):
            try:
                start_atoms = set()
                if isinstance(item, AtomItem):
                    start_atoms.add(item)
                else:
                    # BondItem: start from both ends if available
                    a1 = getattr(item, 'atom1', None)
                    a2 = getattr(item, 'atom2', None)
                    if a1 is not None:
                        start_atoms.add(a1)
                    if a2 is not None:
                        start_atoms.add(a2)

                # BFS/DFS over atoms via bond references (defensive checks)
                atoms_to_visit = list(start_atoms)
                connected_atoms = set()
                connected_bonds = set()

                while atoms_to_visit:
                    a = atoms_to_visit.pop()
                    if a is None:
                        continue
                    if a in connected_atoms:
                        continue
                    connected_atoms.add(a)
                    # iterate bonds attached to atom
                    for b in getattr(a, 'bonds', []) or []:
                        if b is None:
                            continue
                        connected_bonds.add(b)
                        # find the other atom at the bond
                        other = None
                        try:
                            if getattr(b, 'atom1', None) is a:
                                other = getattr(b, 'atom2', None)
                            else:
                                other = getattr(b, 'atom1', None)
                        except Exception:
                            other = None
                        if other is not None and other not in connected_atoms:
                            atoms_to_visit.append(other)

                # Apply selection: clear previous and select only these
                try:
                    self.clearSelection()
                except Exception:
                    pass

                for a in connected_atoms:
                    try:
                        a.setSelected(True)
                    except Exception:
                        try:
                            # fallback: set selected attribute if exists
                            setattr(a, 'selected', True)
                        except Exception:
                            pass
                for b in connected_bonds:
                    try:
                        b.setSelected(True)
                    except Exception:
                        try:
                            setattr(b, 'selected', True)
                        except Exception:
                            pass

                event.accept()
                return
            except Exception:
                # On any unexpected error, fall back to default handling
                pass

        elif self.mode in ['bond_2_5']:
                event.accept()
                return

        super().mouseDoubleClickEvent(event)

    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self.data.add_atom(symbol, pos, charge=charge, radical=radical)
        atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
        self.data.atoms[atom_id]['item'] = atom_item; self.addItem(atom_item); return atom_id


    def create_bond(self, start_atom, end_atom, bond_order=None, bond_stereo=None):
        try:
            if start_atom is None or end_atom is None:
                print("Error: Cannot create bond with None atoms")
                return
                
            exist_b = self.find_bond_between(start_atom, end_atom)
            if exist_b:
                return

            # 引数で次数が指定されていればそれを使用し、なければ現在のモードの値を使用する
            order_to_use = self.bond_order if bond_order is None else bond_order
            stereo_to_use = self.bond_stereo if bond_stereo is None else bond_stereo

            key, status = self.data.add_bond(start_atom.atom_id, end_atom.atom_id, order_to_use, stereo_to_use)
            if status == 'created':
                bond_item = BondItem(start_atom, end_atom, order_to_use, stereo_to_use)
                self.data.bonds[key]['item'] = bond_item
                if hasattr(start_atom, 'bonds'):
                    start_atom.bonds.append(bond_item)
                if hasattr(end_atom, 'bonds'):
                    end_atom.bonds.append(bond_item)
                self.addItem(bond_item)
            
            if hasattr(start_atom, 'update_style'):
                start_atom.update_style()
            if hasattr(end_atom, 'update_style'):
                end_atom.update_style()
                
        except Exception as e:
            print(f"Error creating bond: {e}")
            
            traceback.print_exc()

    def add_molecule_fragment(self, points, bonds_info, existing_items=None, symbol='C'):
        """
        add_molecule_fragment の最終確定版。
        - 既存の結合次数を変更しないポリシーを徹底（最重要）。
        - ベンゼン環テンプレートは、フューズされる既存結合の次数に基づき、
          「新規に作られる二重結合が2本になるように」回転を決定するロジックを適用（条件分岐あり）。
        """
    
        num_points = len(points)
        atom_items = [None] * num_points

        is_benzene_template = (num_points == 6 and any(o == 2 for _, _, o in bonds_info))

    
        def coords(p):
            if hasattr(p, 'x') and hasattr(p, 'y'):
                return (p.x(), p.y())
            try:
                return (p[0], p[1])
            except Exception:
                raise ValueError("point has no x/y")
    
        def dist_pts(a, b):
            ax, ay = coords(a); bx, by = coords(b)
            return math.hypot(ax - bx, ay - by)
    
        # --- 1) 既にクリックされた existing_items をテンプレート頂点にマップ ---
        existing_items = existing_items or []
        used_indices = set()
        ref_lengths = [dist_pts(points[i], points[j]) for i, j, _ in bonds_info if i < num_points and j < num_points]
        avg_len = (sum(ref_lengths) / len(ref_lengths)) if ref_lengths else 20.0
        map_threshold = max(0.5 * avg_len, 8.0)
    
        for ex_item in existing_items:
            try:
                ex_pos = ex_item.pos()
                best_idx, best_d = -1, float('inf')
                for i, p in enumerate(points):
                    if i in used_indices: continue
                    d = dist_pts(p, ex_pos)
                    if best_d is None or d < best_d:
                        best_d, best_idx = d, i
                if best_idx != -1 and best_d <= max(map_threshold, 1.5 * avg_len):
                    atom_items[best_idx] = ex_item
                    used_indices.add(best_idx)
            except Exception:
                pass
    
        # --- 2) シーン内既存原子を self.data.atoms から列挙してマップ ---
        mapped_atoms = {it for it in atom_items if it is not None}
        for i, p in enumerate(points):
            if atom_items[i] is not None: continue
            
            nearby = None
            best_d = float('inf')
            
            for atom_data in self.data.atoms.values():
                a_item = atom_data.get('item')
                if not a_item or a_item in mapped_atoms: continue
                try:
                    d = dist_pts(p, a_item.pos())
                except Exception:
                    continue
                if d < best_d:
                    best_d, nearby = d, a_item

            if nearby and best_d <= map_threshold:
                atom_items[i] = nearby
                mapped_atoms.add(nearby)
    
        # --- 3) 足りない頂点は新規作成　---
        for i, p in enumerate(points):
            if atom_items[i] is None:
                atom_id = self.create_atom(symbol, p)
                atom_items[i] = self.data.atoms[atom_id]['item']
    
        # --- 4) テンプレートのボンド配列を決定（ベンゼン回転合わせの処理） ---
        template_bonds_to_use = list(bonds_info)
        is_6ring = (num_points == 6 and len(bonds_info) == 6)
        template_has_double = any(o == 2 for (_, _, o) in bonds_info)
    
        if is_6ring and template_has_double:
            existing_orders = {} # key: bonds_infoのインデックス, value: 既存の結合次数
            for k, (i_idx, j_idx, _) in enumerate(bonds_info):
                if i_idx < len(atom_items) and j_idx < len(atom_items):
                    a, b = atom_items[i_idx], atom_items[j_idx]
                    if a is None or b is None: continue
                    eb = self.find_bond_between(a, b)
                    if eb:
                        existing_orders[k] = getattr(eb, 'order', 1) 

            if existing_orders:
                orig_orders = [o for (_, _, o) in bonds_info]
                best_rot = 0
                max_score = -999 # スコアは「適合度」を意味する

                # --- フューズされた辺の数による条件分岐 ---
                if len(existing_orders) >= 2:
                    # 2辺以上フューズ: 単純に既存の辺の次数とテンプレートの辺の次数が一致するものを最優先する
                    # (この場合、新しい環を交互配置にするのは難しく、単に既存の構造を壊さないことを優先)
                    for rot in range(num_points):
                        current_score = sum(100 for k, exist_order in existing_orders.items() 
                                            if orig_orders[(k + rot) % num_points] == exist_order)
                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot

                elif len(existing_orders) == 1:
                    # 1辺フューズ: 既存の辺を維持しつつ、その両隣で「反転一致」を達成し、新しい環を交互配置にする
                    
                    # フューズされた辺のインデックスと次数を取得
                    k_fuse = next(iter(existing_orders.keys()))
                    exist_order = existing_orders[k_fuse]
                    
                    # 目標: フューズされた辺の両隣（k-1とk+1）に来るテンプレートの次数が、既存の辺の次数と逆であること
                    # k_adj_1 -> (k_fuse - 1) % 6
                    # k_adj_2 -> (k_fuse + 1) % 6
                    
                    for rot in range(num_points):
                        current_score = 0
                        rotated_template_order = orig_orders[(k_fuse + rot) % num_points]

                        # 1. まず、フューズされた辺自体が次数を反転させられる位置にあるかチェック（必須ではないが、回転を絞る）
                        if (exist_order == 1 and rotated_template_order == 2) or \
                           (exist_order == 2 and rotated_template_order == 1):
                            current_score += 100 # 大幅ボーナス: 理想的な回転

                        # 2. 次に、両隣の辺の次数をチェック（交互配置維持の主目的）
                        # 既存辺の両隣は、新規に作成されるため、テンプレートの次数でボンドが作成されます。
                        # ここで、テンプレートの次数が既存辺の次数と逆になる回転を選ぶ必要があります。
                        
                        # テンプレートの辺は、回転後のk_fuseの両隣（m_adj1, m_adj2）
                        m_adj1 = (k_fuse - 1 + rot) % num_points 
                        m_adj2 = (k_fuse + 1 + rot) % num_points
                        
                        neighbor_order_1 = orig_orders[m_adj1]
                        neighbor_order_2 = orig_orders[m_adj2]

                        # 既存が単結合(1)の場合、両隣は二重結合(2)であってほしい
                        if exist_order == 1:
                            if neighbor_order_1 == 2: current_score += 50
                            if neighbor_order_2 == 2: current_score += 50
                        
                        # 既存が二重結合(2)の場合、両隣は単結合(1)であってほしい
                        elif exist_order == 2:
                            if neighbor_order_1 == 1: current_score += 50
                            if neighbor_order_2 == 1: current_score += 50
                            
                        # 3. タイブレーク: その他の既存結合（フューズ辺ではない）との次数一致度も加味
                        for k, e_order in existing_orders.items():
                             if k != k_fuse:
                                r_t_order = orig_orders[(k + rot) % num_points]
                                if r_t_order == e_order: current_score += 10 # 既存構造維持のボーナス
                        
                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot
                
                # 最終的な回転を反映
                new_tb = []
                for m in range(num_points):
                    i_idx, j_idx, _ = bonds_info[m]
                    new_order = orig_orders[(m + best_rot) % num_points]
                    new_tb.append((i_idx, j_idx, new_order))
                template_bonds_to_use = new_tb
    
        # --- 5) ボンド作成／更新---
        for id1_idx, id2_idx, order in template_bonds_to_use:
            if id1_idx < len(atom_items) and id2_idx < len(atom_items):
                a_item, b_item = atom_items[id1_idx], atom_items[id2_idx]
                if not a_item or not b_item or a_item is b_item: continue

                id1, id2 = a_item.atom_id, b_item.atom_id
                if id1 > id2: id1, id2 = id2, id1

                exist_b = self.find_bond_between(a_item, b_item)

                if exist_b:
                    # デフォルトでは既存の結合を維持する
                    should_overwrite = False

                    # 条件1: ベンゼン環テンプレートであること
                    # 条件2: 接続先が単結合であること
                    if is_benzene_template and exist_b.order == 1:

                        # 条件3: 接続先の単結合が共役系の一部ではないこと
                        # (つまり、両端の原子が他に二重結合を持たないこと)
                        atom1 = exist_b.atom1
                        atom2 = exist_b.atom2

                        # atom1が他に二重結合を持つかチェック
                        atom1_has_other_double_bond = any(b.order == 2 for b in atom1.bonds if b is not exist_b)

                        # atom2が他に二重結合を持つかチェック
                        atom2_has_other_double_bond = any(b.order == 2 for b in atom2.bonds if b is not exist_b)

                        # 両方の原子が他に二重結合を持たない「孤立した単結合」の場合のみ上書きフラグを立てる
                        if not atom1_has_other_double_bond and not atom2_has_other_double_bond:
                            should_overwrite = True

                    if should_overwrite:
                        # 上書き条件が全て満たされた場合にのみ、結合次数を更新
                        exist_b.order = order
                        exist_b.stereo = 0
                        self.data.bonds[(id1, id2)]['order'] = order
                        self.data.bonds[(id1, id2)]['stereo'] = 0
                        exist_b.update()
                    else:
                        # 上書き条件を満たさない場合は、既存の結合を維持する
                        continue
                else:
                    # 新規ボンド作成
                    self.create_bond(a_item, b_item, bond_order=order, bond_stereo=0)
        
        # --- 6) 表示更新　---
        for at in atom_items:
            try:
                if at: at.update_style() 
            except Exception:
                pass
    
        return atom_items


    def update_template_preview(self, pos):
        mode_parts = self.mode.split('_')
        
        # Check if this is a user template
        if len(mode_parts) >= 3 and mode_parts[1] == 'user':
            self.update_user_template_preview(pos)
            return
        
        is_aromatic = False
        if mode_parts[1] == 'benzene':
            n = 6
            is_aromatic = True
        else:
            try: n = int(mode_parts[1])
            except ValueError: return

        items_under = self.items(pos)  # top-most first
        item = None
        for it in items_under:
            if isinstance(it, (AtomItem, BondItem)):
                item = it
                break

        points, bonds_info = [], []
        l = DEFAULT_BOND_LENGTH
        self.template_context = {}


        if isinstance(item, AtomItem):
            p0 = item.pos()
            continuous_angle = math.atan2(pos.y() - p0.y(), pos.x() - p0.x())
            snap_angle_rad = math.radians(15)
            snapped_angle = round(continuous_angle / snap_angle_rad) * snap_angle_rad
            p1 = p0 + QPointF(l * math.cos(snapped_angle), l * math.sin(snapped_angle))
            points = self._calculate_polygon_from_edge(p0, p1, n)
            self.template_context['items'] = [item]

        elif isinstance(item, BondItem):
            # 結合にスナップ
            p0, p1 = item.atom1.pos(), item.atom2.pos()
            points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=pos, use_existing_length=True)
            self.template_context['items'] = [item.atom1, item.atom2]

        else:
            angle_step = 2 * math.pi / n
            start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2
            points = [
                pos + QPointF(l * math.cos(start_angle + i * angle_step), l * math.sin(start_angle + i * angle_step))
                for i in range(n)
            ]

        if points:
            if is_aromatic:
                bonds_info = [(i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)]
            else:
                bonds_info = [(i, (i + 1) % n, 1) for i in range(n)]

            self.template_context['points'] = points
            self.template_context['bonds_info'] = bonds_info

            self.template_preview.set_geometry(points, is_aromatic)

            self.template_preview.show()
            if self.views():
                self.views()[0].viewport().update()
        else:
            self.template_preview.hide()
            if self.views():
                self.views()[0].viewport().update()

    def _calculate_polygon_from_edge(self, p0, p1, n, cursor_pos=None, use_existing_length=False):
        if n < 3: return []
        v_edge = p1 - p0
        edge_length = math.sqrt(v_edge.x()**2 + v_edge.y()**2)
        if edge_length == 0: return []
        
        target_length = edge_length if use_existing_length else DEFAULT_BOND_LENGTH
        
        v_edge = (v_edge / edge_length) * target_length
        
        if not use_existing_length:
             p1 = p0 + v_edge

        points = [p0, p1]
        
        interior_angle = (n - 2) * math.pi / n
        rotation_angle = math.pi - interior_angle
        
        if cursor_pos:
            # Note: v_edgeは正規化済みだが、方向は同じなので判定には問題ない
            v_cursor = cursor_pos - p0
            cross_product_z = (p1 - p0).x() * v_cursor.y() - (p1 - p0).y() * v_cursor.x()
            if cross_product_z < 0:
                rotation_angle = -rotation_angle

        cos_a, sin_a = math.cos(rotation_angle), math.sin(rotation_angle)
        
        current_p, current_v = p1, v_edge
        for _ in range(n - 2):
            new_vx = current_v.x() * cos_a - current_v.y() * sin_a
            new_vy = current_v.x() * sin_a + current_v.y() * cos_a
            current_v = QPointF(new_vx, new_vy)
            current_p = current_p + current_v
            points.append(current_p)
        return points

    def delete_items(self, items_to_delete):
        """指定されたアイテムセット（原子・結合）を安全な順序で削除する修正版"""
        # Hardened deletion: perform data-model removals first, then scene removals,
        # and always defensively check attributes to avoid accessing partially-deleted objects.
        if not items_to_delete:
            return False

        # First sanitize the incoming collection: only keep live, expected QGraphics wrappers
        try:
            sanitized = set()
            for it in items_to_delete:
                try:
                    if it is None:
                        continue
                    # Skip SIP-deleted wrappers early to avoid native crashes
                    if sip_isdeleted_safe(it):
                        continue
                    # Only accept AtomItem/BondItem or other QGraphicsItem subclasses
                    if isinstance(it, (AtomItem, BondItem, QGraphicsItem)):
                        sanitized.add(it)
                except Exception:
                    # If isinstance or sip check raises, skip this entry
                    continue
            items_to_delete = sanitized
        except Exception:
            # If sanitization fails, fall back to original input and proceed defensively
            pass

        try:
            atoms_to_delete = {item for item in items_to_delete if isinstance(item, AtomItem)}
            bonds_to_delete = {item for item in items_to_delete if isinstance(item, BondItem)}

            # Include bonds attached to atoms being deleted
            for atom in list(atoms_to_delete):
                try:
                    if hasattr(atom, 'bonds') and atom.bonds:
                        for b in list(atom.bonds):
                            bonds_to_delete.add(b)
                except Exception:
                    # If accessing bonds raises (item partially deleted), skip
                    continue

            # Determine atoms that will remain but whose bond lists must be updated
            atoms_to_update = set()
            for bond in list(bonds_to_delete):
                try:
                    a1 = getattr(bond, 'atom1', None)
                    a2 = getattr(bond, 'atom2', None)
                    if a1 and a1 not in atoms_to_delete:
                        atoms_to_update.add(a1)
                    if a2 and a2 not in atoms_to_delete:
                        atoms_to_update.add(a2)
                except Exception:
                    continue

            # 1) Update surviving atoms' bond lists to remove references to bonds_to_delete
            #    (Important: remove BondItem references so atoms properly reflect
            #     that they have no remaining bonds and update visibility accordingly.)
            for atom in list(atoms_to_update):
                try:
                    if sip_isdeleted_safe(atom):
                        continue
                    # Defensive: if the atom has a bonds list, filter out bonds being deleted
                    if hasattr(atom, 'bonds') and atom.bonds:
                        try:
                            # Replace in-place to preserve any other references.
                            # Avoid touching SIP-deleted bond wrappers: build a set
                            # of live bonds-to-delete and also prune any SIP-deleted
                            # entries that may exist in atom.bonds.
                            live_btd = {b for b in bonds_to_delete if not sip_isdeleted_safe(b)}

                            # First, remove any SIP-deleted bond wrappers from atom.bonds
                            atom.bonds[:] = [b for b in atom.bonds if not sip_isdeleted_safe(b)]

                            # Then remove bonds which are in the live_btd set
                            if live_btd:
                                atom.bonds[:] = [b for b in atom.bonds if b not in live_btd]
                        except Exception:
                            # Fall back to iterative removal if list comprehension fails
                            try:
                                live_btd = [b for b in list(bonds_to_delete) if not sip_isdeleted_safe(b)]
                                for b in live_btd:
                                    if b in atom.bonds:
                                        atom.bonds.remove(b)
                            except Exception:
                                pass

                    # After pruning bond references, update visual style so carbons without
                    # bonds become visible again.
                    if hasattr(atom, 'update_style'):
                        atom.update_style()
                except Exception:
                    continue

            # 2) Remove bonds/atoms from the data model first (so other code reading the model
            #    doesn't encounter stale entries while we are removing graphics)
            for bond in list(bonds_to_delete):
                try:
                    a1 = getattr(bond, 'atom1', None)
                    a2 = getattr(bond, 'atom2', None)
                    if a1 and a2 and hasattr(self, 'data'):
                        try:
                            self.data.remove_bond(a1.atom_id, a2.atom_id)
                        except Exception:
                            # try reverse order if remove_bond expects ordered tuple
                            try:
                                self.data.remove_bond(a2.atom_id, a1.atom_id)
                            except Exception:
                                pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    if hasattr(atom, 'atom_id') and hasattr(self, 'data'):
                        try:
                            self.data.remove_atom(atom.atom_id)
                        except Exception:
                            pass
                except Exception:
                    continue

            # Invalidate any pending implicit-hydrogen UI updates because the
            # underlying data model changed. This prevents a scheduled
            # update_implicit_hydrogens closure from touching atoms/bonds that
            # were just removed. Do a single increment rather than one per-atom.
            try:
                self._ih_update_counter += 1
            except Exception:
                try:
                    self._ih_update_counter = 0
                except Exception:
                    pass

            # 3) Remove graphic items from the scene (bonds first)
            # To avoid calling into methods on wrappers that may refer to
            # already-deleted C++ objects (which can cause a native crash when
            # SIP is not available), take a snapshot of the current scene's
            # items and use membership tests instead of calling item.scene().
            try:
                current_scene_items = set(self.items())
            except Exception:
                # If for any reason items() fails, fall back to an empty set
                current_scene_items = set()

            for bond in list(bonds_to_delete):
                try:
                    # If the SIP wrapper is already deleted, skip it.
                    if sip_isdeleted_safe(bond):
                        continue
                    # Only attempt to remove the bond if it is present in the
                    # scene snapshot. This avoids calling bond.scene() which
                    # may invoke C++ on a deleted object.
                    if bond in current_scene_items:
                        try:
                            self.removeItem(bond)
                        except Exception:
                            pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    # Skip if wrapper is reported deleted by SIP
                    if sip_isdeleted_safe(atom):
                        continue
                    if atom in current_scene_items:
                        try:
                            self.removeItem(atom)
                        except Exception:
                            pass
                except Exception:
                    continue

            # 4) Instead of aggressively nullling object attributes (which can
            #    lead to C++/SIP finalization races and segfaults), keep a
            #    strong reference to the deleted wrappers for the lifetime of
            #    the scene. This prevents their underlying SIP wrappers from
            #    being finalized while other code may still touch them.
            try:
                if not hasattr(self, '_deleted_items') or self._deleted_items is None:
                    self._deleted_items = []
            except Exception:
                self._deleted_items = []

            for bond in list(bonds_to_delete):
                try:
                    # Hide the graphics item if possible and stash it
                    if not sip_isdeleted_safe(bond):
                        try:
                            bond.hide()
                        except Exception:
                            pass
                        try:
                            self._deleted_items.append(bond)
                        except Exception:
                            # Swallow any error while stashing
                            pass
                except Exception:
                    continue

            for atom in list(atoms_to_delete):
                try:
                    if not sip_isdeleted_safe(atom):
                        try:
                            atom.hide()
                        except Exception:
                            pass
                        try:
                            self._deleted_items.append(atom)
                        except Exception:
                            pass
                except Exception:
                    continue

            # 5) Final visual updates for surviving atoms
            for atom in list(atoms_to_update):
                try:
                    if hasattr(atom, 'update_style'):
                        atom.update_style()
                except Exception:
                    continue

            return True

        except Exception as e:
            # Keep the application alive on unexpected errors
            print(f"Error during delete_items operation: {e}")
            
            traceback.print_exc()
            return False
    def purge_deleted_items(self):
        """Purge and release any held deleted-wrapper references.

        This is intended to be invoked on application shutdown to allow
        the process to release references to SIP/C++ wrappers that were
        kept around to avoid finalization races during normal runtime.
        The method is defensive: it tolerates partially-deleted wrappers
        and any SIP unavailability.
        """
        try:
            if not hasattr(self, '_deleted_items') or not self._deleted_items:
                return

            # Iterate a copy since we will clear the list.
            for obj in list(self._deleted_items):
                try:
                    # If the wrapper is still alive, attempt to hide it so
                    # the graphics subsystem isn't holding on to resources.
                    if not sip_isdeleted_safe(obj):
                        try:
                            obj.hide()
                        except Exception:
                            pass

                    # Try to clear container attributes that may hold refs
                    # to other scene objects (bonds, etc.) to help GC.
                    try:
                        if hasattr(obj, 'bonds') and getattr(obj, 'bonds') is not None:
                            try:
                                obj.bonds.clear()
                            except Exception:
                                # Try assignment fallback
                                try:
                                    obj.bonds = []
                                except Exception:
                                    pass
                    except Exception:
                        pass

                except Exception:
                    # Continue purging remaining items even if one fails.
                    continue

            # Finally, drop our references.
            try:
                self._deleted_items.clear()
            except Exception:
                try:
                    self._deleted_items = []
                except Exception:
                    pass

        except Exception as e:
            # Never raise during shutdown
            try:
                print(f"Error purging deleted items: {e}")
            except Exception:
                pass
    
    def add_user_template_fragment(self, context):
        """ユーザーテンプレートフラグメントを配置"""
        points = context.get('points', [])
        bonds_info = context.get('bonds_info', [])
        atoms_data = context.get('atoms_data', [])
        attachment_atom = context.get('attachment_atom')
        
        if not points or not atoms_data:
            return
        
        # Create atoms
        atom_id_map = {}  # template id -> scene atom id
        
        for i, (pos, atom_data) in enumerate(zip(points, atoms_data)):
            # Skip first atom if attaching to existing atom
            if i == 0 and attachment_atom:
                atom_id_map[atom_data['id']] = attachment_atom.atom_id
                continue
            
            symbol = atom_data.get('symbol', 'C')
            charge = atom_data.get('charge', 0)
            radical = atom_data.get('radical', 0)
            
            atom_id = self.data.add_atom(symbol, pos, charge, radical)
            atom_id_map[atom_data['id']] = atom_id
            
            # Create visual atom item
            atom_item = AtomItem(atom_id, symbol, pos, charge, radical)
            self.data.atoms[atom_id]['item'] = atom_item
            self.addItem(atom_item)
        
        # Create bonds (bonds_infoは必ずidベースで扱う)
        # まずindex→id変換テーブルを作る
        index_to_id = [atom_data.get('id', i) for i, atom_data in enumerate(atoms_data)]
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                # bonds_infoの0,1番目がindexならidに変換
                atom1_idx = bond_info[0]
                atom2_idx = bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                stereo = bond_info[3] if len(bond_info) > 3 else 0

                # index→id変換（すでにidならそのまま）
                if isinstance(atom1_idx, int) and atom1_idx < len(index_to_id):
                    template_atom1_id = index_to_id[atom1_idx]
                else:
                    template_atom1_id = atom1_idx
                if isinstance(atom2_idx, int) and atom2_idx < len(index_to_id):
                    template_atom2_id = index_to_id[atom2_idx]
                else:
                    template_atom2_id = atom2_idx

                atom1_id = atom_id_map.get(template_atom1_id)
                atom2_id = atom_id_map.get(template_atom2_id)

                if atom1_id is not None and atom2_id is not None:
                    # Skip if bond already exists
                    existing_bond = None
                    if (atom1_id, atom2_id) in self.data.bonds:
                        existing_bond = (atom1_id, atom2_id)
                    elif (atom2_id, atom1_id) in self.data.bonds:
                        existing_bond = (atom2_id, atom1_id)

                    if not existing_bond:
                        bond_key, _ = self.data.add_bond(atom1_id, atom2_id, order, stereo)
                        # Create visual bond item
                        atom1_item = self.data.atoms[atom1_id]['item']
                        atom2_item = self.data.atoms[atom2_id]['item']
                        if atom1_item and atom2_item:
                            bond_item = BondItem(atom1_item, atom2_item, order, stereo)
                            self.data.bonds[bond_key]['item'] = bond_item
                            self.addItem(bond_item)
                            atom1_item.bonds.append(bond_item)
                            atom2_item.bonds.append(bond_item)
        
        # Update atom visuals
        for atom_id in atom_id_map.values():
            if atom_id in self.data.atoms and self.data.atoms[atom_id]['item']:
                self.data.atoms[atom_id]['item'].update_style()
    
    def update_user_template_preview(self, pos):
        """ユーザーテンプレートのプレビューを更新"""
        # Robust user template preview: do not access self.data.atoms for preview-only atoms
        if not hasattr(self, 'user_template_data') or not self.user_template_data:
            return

        template_data = self.user_template_data
        atoms = template_data.get('atoms', [])
        bonds = template_data.get('bonds', [])

        if not atoms:
            return

        # Find attachment point (first atom or clicked item)
        items_under = self.items(pos)
        attachment_atom = None
        for item in items_under:
            if isinstance(item, AtomItem):
                attachment_atom = item
                break

        # Calculate template positions
        points = []
        # Find template bounds for centering
        if atoms:
            min_x = min(atom['x'] for atom in atoms)
            max_x = max(atom['x'] for atom in atoms)
            min_y = min(atom['y'] for atom in atoms)
            max_y = max(atom['y'] for atom in atoms)
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
        # Position template
        if attachment_atom:
            # Attach to existing atom
            attach_pos = attachment_atom.pos()
            offset_x = attach_pos.x() - atoms[0]['x']
            offset_y = attach_pos.y() - atoms[0]['y']
        else:
            # Center at cursor position
            offset_x = pos.x() - center_x
            offset_y = pos.y() - center_y
        # Calculate atom positions
        for atom in atoms:
            new_pos = QPointF(atom['x'] + offset_x, atom['y'] + offset_y)
            points.append(new_pos)
        # Create atom ID to index mapping (for preview only)
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get('id', i)
            atom_id_to_index[atom_id] = i
        # bonds_info をテンプレートの bonds から生成
        bonds_info = []
        for bond in bonds:
            atom1_idx = atom_id_to_index.get(bond['atom1'])
            atom2_idx = atom_id_to_index.get(bond['atom2'])
            if atom1_idx is not None and atom2_idx is not None:
                order = bond.get('order', 1)
                stereo = bond.get('stereo', 0)
                bonds_info.append((atom1_idx, atom2_idx, order, stereo))
        # プレビュー用: points, bonds_info から線を描画
        # 設置用 context を保存
        self.template_context = {
            'points': points,
            'bonds_info': bonds_info,
            'atoms_data': atoms,
            'attachment_atom': attachment_atom,
        }
        # 既存のプレビューアイテムを一旦クリア
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(item, '_is_template_preview', False):
                self.removeItem(item)

        # Draw preview lines only using calculated points (do not access self.data.atoms)
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                i, j = bond_info[0], bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                # stereo = bond_info[3] if len(bond_info) > 3 else 0
                if i < len(points) and j < len(points):
                    line = QGraphicsLineItem(QLineF(points[i], points[j]))
                    pen = QPen(Qt.black, 2 if order == 2 else 1)
                    line.setPen(pen)
                    line._is_template_preview = True  # フラグで区別
                    self.addItem(line)
        # Never access self.data.atoms here for preview-only atoms

    def leaveEvent(self, event):
        self.template_preview.hide(); super().leaveEvent(event)

    def set_hovered_item(self, item):
        """BondItemから呼ばれ、ホバー中のアイテムを記録する"""
        self.hovered_item = item

    def keyPressEvent(self, event):
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())
        key = event.key()
        modifiers = event.modifiers()
        
        if not self.window.is_2d_editable:
            return    


        if key == Qt.Key.Key_4:
            # --- 動作1: カーソルが原子/結合上にある場合 (ワンショットでテンプレート配置) ---
            if isinstance(item_at_cursor, (AtomItem, BondItem)):
                
                # ベンゼンテンプレートのパラメータを設定
                n, is_aromatic = 6, True
                points, bonds_info, existing_items = [], [], []
                
                # update_template_preview と同様のロジックで配置情報を計算
                if isinstance(item_at_cursor, AtomItem):
                    p0 = item_at_cursor.pos()
                    l = DEFAULT_BOND_LENGTH
                    direction = QLineF(p0, cursor_pos).unitVector()
                    p1 = p0 + direction.p2() * l if direction.length() > 0 else p0 + QPointF(l, 0)
                    points = self._calculate_polygon_from_edge(p0, p1, n)
                    existing_items = [item_at_cursor]

                elif isinstance(item_at_cursor, BondItem):
                    p0, p1 = item_at_cursor.atom1.pos(), item_at_cursor.atom2.pos()
                    points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=cursor_pos, use_existing_length=True)
                    existing_items = [item_at_cursor.atom1, item_at_cursor.atom2]
                
                if points:
                    bonds_info = [(i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)]
                    
                    # 計算した情報を使って、その場にフラグメントを追加
                    self.add_molecule_fragment(points, bonds_info, existing_items=existing_items)
                    self.window.push_undo_state()

            # --- 動作2: カーソルが空白領域にある場合 (モード切替) ---
            else:
                self.window.set_mode_and_update_toolbar('template_benzene')

            event.accept()
            return

        # --- 0a. ラジカルの変更 (.) ---
        if key == Qt.Key.Key_Period:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                for atom in target_atoms:
                    # ラジカルの状態をトグル (0 -> 1 -> 2 -> 0)
                    atom.prepareGeometryChange()
                    atom.radical = (atom.radical + 1) % 3
                    self.data.atoms[atom.atom_id]['radical'] = atom.radical
                    atom.update_style()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 0b. 電荷の変更 (+/-キー) ---
        if key == Qt.Key.Key_Plus or key == Qt.Key.Key_Minus:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                delta = 1 if key == Qt.Key.Key_Plus else -1
                for atom in target_atoms:
                    atom.prepareGeometryChange()
                    atom.charge += delta
                    self.data.atoms[atom.atom_id]['charge'] = atom.charge
                    atom.update_style()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 1. Atomに対する操作 (元素記号の変更) ---
        if isinstance(item_at_cursor, AtomItem):
            new_symbol = None
            if modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_symbol_map:
                new_symbol = self.key_to_symbol_map[key]
            elif modifiers == Qt.KeyboardModifier.ShiftModifier and key in self.key_to_symbol_map_shift:
                new_symbol = self.key_to_symbol_map_shift[key]

            if new_symbol and item_at_cursor.symbol != new_symbol:
                item_at_cursor.prepareGeometryChange()
                
                item_at_cursor.symbol = new_symbol
                self.data.atoms[item_at_cursor.atom_id]['symbol'] = new_symbol
                item_at_cursor.update_style()


                atoms_to_update = {item_at_cursor}
                for bond in item_at_cursor.bonds:
                    bond.update()
                    other_atom = bond.atom1 if bond.atom2 is item_at_cursor else bond.atom2
                    atoms_to_update.add(other_atom)

                for atom in atoms_to_update:
                    atom.update_style()

                self.window.push_undo_state()
                event.accept()
                return

        # --- 2. Bondに対する操作 (次数・立体化学の変更) ---
        target_bonds = []
        if isinstance(item_at_cursor, BondItem):
            target_bonds = [item_at_cursor]
        else:
            target_bonds = [it for it in self.selectedItems() if isinstance(it, BondItem)]

        if target_bonds:
            any_bond_changed = False
            for bond in target_bonds:
                # 1. 結合の向きを考慮して、データ辞書内の現在のキーを正しく特定する
                id1, id2 = bond.atom1.atom_id, bond.atom2.atom_id
                current_key = None
                if (id1, id2) in self.data.bonds:
                    current_key = (id1, id2)
                elif (id2, id1) in self.data.bonds:
                    current_key = (id2, id1)
                
                if not current_key: continue

                # 2. 変更前の状態を保存
                old_order, old_stereo = bond.order, bond.stereo

                # 3. キー入力に応じてBondItemのプロパティを変更
                if key == Qt.Key.Key_W:
                    if bond.stereo == 1:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1; bond.stereo = 1

                elif key == Qt.Key.Key_D:
                    if bond.stereo == 2:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1; bond.stereo = 2

                elif key == Qt.Key.Key_1 and (bond.order != 1 or bond.stereo != 0):
                    bond.order = 1; bond.stereo = 0
                elif key == Qt.Key.Key_2 and (bond.order != 2 or bond.stereo != 0):
                    bond.order = 2; bond.stereo = 0; needs_update = True
                elif key == Qt.Key.Key_3 and bond.order != 3:
                    bond.order = 3; bond.stereo = 0; needs_update = True

                # 4. 実際に変更があった場合のみデータモデルを更新
                if old_order != bond.order or old_stereo != bond.stereo:
                    any_bond_changed = True
                    
                    # 5. 古いキーでデータを辞書から一度削除
                    bond_data = self.data.bonds.pop(current_key)
                    bond_data['order'] = bond.order
                    bond_data['stereo'] = bond.stereo

                    # 6. 変更後の種類に応じて新しいキーを決定し、再登録する
                    new_key_id1, new_key_id2 = bond.atom1.atom_id, bond.atom2.atom_id
                    if bond.stereo == 0:
                        if new_key_id1 > new_key_id2:
                            new_key_id1, new_key_id2 = new_key_id2, new_key_id1
                    
                    new_key = (new_key_id1, new_key_id2)
                    self.data.bonds[new_key] = bond_data
                    
                    bond.update()

            if any_bond_changed:
                self.window.push_undo_state()
            
            if key in [Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3, Qt.Key.Key_W, Qt.Key.Key_D]:
                event.accept()
                return

        if isinstance(self.hovered_item, BondItem) and self.hovered_item.order == 2:
            if event.key() == Qt.Key.Key_Z:
                self.update_bond_stereo(self.hovered_item, 3)  # Z-isomer
                self.window.push_undo_state()
                event.accept()
                return
            elif event.key() == Qt.Key.Key_E:
                self.update_bond_stereo(self.hovered_item, 4)  # E-isomer
                self.window.push_undo_state()
                event.accept()
                return
                    
        # --- 3. Atomに対する操作 (原子の追加 - マージされた機能) ---
        if key == Qt.Key.Key_1:
            start_atom = None
            if isinstance(item_at_cursor, AtomItem):
                start_atom = item_at_cursor
            else:
                selected_atoms = [item for item in self.selectedItems() if isinstance(item, AtomItem)]
                if len(selected_atoms) == 1:
                    start_atom = selected_atoms[0]

            if start_atom:
                start_pos = start_atom.pos()
                l = DEFAULT_BOND_LENGTH
                new_pos_offset = QPointF(0, -l) # デフォルトのオフセット (上)

                # 接続している原子のリストを取得 (H原子以外)
                neighbor_positions = []
                for bond in start_atom.bonds:
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    if other_atom.symbol != 'H': # 水素原子を無視 (四面体構造の考慮のため)
                        neighbor_positions.append(other_atom.pos())

                num_non_H_neighbors = len(neighbor_positions)
                
                if num_non_H_neighbors == 0:
                    # 結合ゼロ: デフォルト方向
                    new_pos_offset = QPointF(0, -l)
                
                elif num_non_H_neighbors == 1:
                    # 結合1本: 既存結合と約120度（または60度）の角度
                    bond = start_atom.bonds[0]
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    existing_bond_vector = start_pos - other_atom.pos()
                    
                    # 既存の結合から時計回り60度回転 (ベンゼン環のような構造にしやすい)
                    angle_rad = math.radians(60) 
                    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
                    vx, vy = existing_bond_vector.x(), existing_bond_vector.y()
                    new_vx, new_vy = vx * cos_a - vy * sin_a, vx * sin_a + vy * cos_a
                    rotated_vector = QPointF(new_vx, new_vy)
                    line = QLineF(QPointF(0, 0), rotated_vector)
                    line.setLength(l)
                    new_pos_offset = line.p2()

                elif num_non_H_neighbors == 3:

                    bond_vectors_sum = QPointF(0, 0)
                    for pos in neighbor_positions:
                        # start_pos から neighbor_pos へのベクトル
                        vec = pos - start_pos 
                        # 単位ベクトルに変換
                        line_to_other = QLineF(QPointF(0,0), vec)
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2()
                    
                    # SUM_TOLERANCE is now a module-level constant
                    if bond_vectors_sum.manhattanLength() > SUM_TOLERANCE:
                        new_direction_line = QLineF(QPointF(0,0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        new_pos_offset = QPointF(l * 0.7071, -l * 0.7071) 


                else: # 2本または4本以上の場合 (一般的な骨格の継続、または過結合)
                    bond_vectors_sum = QPointF(0, 0)
                    for bond in start_atom.bonds:
                        other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                        line_to_other = QLineF(start_pos, other_atom.pos())
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2() - line_to_other.p1()
                    
                    if bond_vectors_sum.manhattanLength() > 0.01:
                        new_direction_line = QLineF(QPointF(0,0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        # 総和がゼロの場合は、デフォルト（上）
                        new_pos_offset = QPointF(0, -l)


                # SNAP_DISTANCE is a module-level constant
                target_pos = start_pos + new_pos_offset
                
                # 近くに原子を探す
                near_atom = self.find_atom_near(target_pos, tol=SNAP_DISTANCE)
                
                if near_atom and near_atom is not start_atom:
                    # 近くに既存原子があれば結合
                    self.create_bond(start_atom, near_atom)
                else:
                    # 新規原子を作成し結合
                    new_atom_id = self.create_atom('C', target_pos)
                    new_atom_item = self.data.atoms[new_atom_id]['item']
                    self.create_bond(start_atom, new_atom_item)

                self.clearSelection()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 4. 全体に対する操作 (削除、モード切替など) ---
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            if self.temp_line:
                try:
                    if not sip_isdeleted_safe(self.temp_line):
                        try:
                            if getattr(self.temp_line, 'scene', None) and self.temp_line.scene():
                                self.removeItem(self.temp_line)
                        except Exception:
                            pass
                except Exception:
                    try:
                        self.removeItem(self.temp_line)
                    except Exception:
                        pass
                self.temp_line = None; self.start_atom = None; self.start_pos = None
                self.initial_positions_in_event = {}
                event.accept()
                return

            items_to_process = set(self.selectedItems()) 
            # カーソル下のアイテムも削除対象に加える
            if item_at_cursor and isinstance(item_at_cursor, (AtomItem, BondItem)):
                items_to_process.add(item_at_cursor)

            if self.delete_items(items_to_process):
                self.window.push_undo_state()
                self.window.statusBar().showMessage("Deleted selected items.")

            # もしデータモデル内の原子が全て無くなっていたら、シーンをクリアして初期状態に戻す
            if not self.data.atoms:
                # 1. シーン上の全グラフィックアイテムを削除する
                self.clear() 

                # 2. テンプレートプレビューなど、初期状態で必要なアイテムを再生成する
                self.reinitialize_items()
                
                # 3. 結合描画中などの一時的な状態も完全にリセットする
                self.temp_line = None
                self.start_atom = None
                self.start_pos = None
                self.initial_positions_in_event = {}
                
                # このイベントはここで処理完了とする
                event.accept()
                return
    
            # 描画の強制更新
            if self.views():
                self.views()[0].viewport().update() 
                QApplication.processEvents()
    
                event.accept()
                return
        

        if key == Qt.Key.Key_Space:
            if self.mode != 'select':
                self.window.activate_select_mode()
            else:
                self.window.select_all()
            event.accept()
            return

        # グローバルな描画モード切替
        mode_to_set = None

        # 1. 原子描画モードへの切り替え
        symbol_for_mode_change = None
        if modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_symbol_map:
            symbol_for_mode_change = self.key_to_symbol_map[key]
        elif modifiers == Qt.KeyboardModifier.ShiftModifier and key in self.key_to_symbol_map_shift:
            symbol_for_mode_change = self.key_to_symbol_map_shift[key]
        
        if symbol_for_mode_change:
            mode_to_set = f'atom_{symbol_for_mode_change}'

        # 2. 結合描画モードへの切り替え
        elif modifiers == Qt.KeyboardModifier.NoModifier and key in self.key_to_bond_mode_map:
            mode_to_set = self.key_to_bond_mode_map[key]

        # モードが決定されていれば、モード変更を実行
        if mode_to_set:
            if hasattr(self.window, 'set_mode_and_update_toolbar'):
                 self.window.set_mode_and_update_toolbar(mode_to_set)
                 event.accept()
                 return
        
        # --- どの操作にも当てはまらない場合 ---
        super().keyPressEvent(event)
        
    def find_atom_near(self, pos, tol=14.0):
        # Create a small search rectangle around the position
        search_rect = QRectF(pos.x() - tol, pos.y() - tol, 2 * tol, 2 * tol)
        nearby_items = self.items(search_rect)

        for it in nearby_items:
            if isinstance(it, AtomItem):
                # Check the precise distance only for candidate items
                if QLineF(it.pos(), pos).length() <= tol:
                    return it
        return None

    def find_bond_between(self, atom1, atom2):
        for b in atom1.bonds:
            if (b.atom1 is atom1 and b.atom2 is atom2) or \
               (b.atom1 is atom2 and b.atom2 is atom1):
                return b
        return None

    def update_bond_stereo(self, bond_item, new_stereo):
        """結合の立体化学を更新する共通メソッド"""
        try:
            if bond_item is None:
                print("Error: bond_item is None in update_bond_stereo")
                return
                
            if bond_item.order != 2 or bond_item.stereo == new_stereo:
                return

            if not hasattr(bond_item, 'atom1') or not hasattr(bond_item, 'atom2'):
                print("Error: bond_item missing atom references")
                return
                
            if bond_item.atom1 is None or bond_item.atom2 is None:
                print("Error: bond_item has None atom references")
                return
                
            if not hasattr(bond_item.atom1, 'atom_id') or not hasattr(bond_item.atom2, 'atom_id'):
                print("Error: bond atoms missing atom_id")
                return

            id1, id2 = bond_item.atom1.atom_id, bond_item.atom2.atom_id

            # E/Z結合は方向性を持つため、キーは(id1, id2)のまま探す
            key_to_update = (id1, id2)
            if key_to_update not in self.data.bonds:
                # Wedge/Dashなど、逆順で登録されている可能性も考慮
                key_to_update = (id2, id1)
                if key_to_update not in self.data.bonds:
                    # Log error instead of printing to console
                    if hasattr(self.window, 'statusBar'):
                        self.window.statusBar().showMessage(f"Warning: Bond between atoms {id1} and {id2} not found in data model.", 3000)
                    print(f"Error: Bond key not found: {id1}-{id2} or {id2}-{id1}")
                    return
                    
            # Update data model
            self.data.bonds[key_to_update]['stereo'] = new_stereo
            
            # Update visual representation
            bond_item.set_stereo(new_stereo)
            
            self.data_changed_in_event = True
            
        except Exception as e:
            print(f"Error in update_bond_stereo: {e}")
            
            traceback.print_exc()
            if hasattr(self.window, 'statusBar'):
                self.window.statusBar().showMessage(f"Error updating bond stereochemistry: {e}", 5000)

class ZoomableView(QGraphicsView):
    """ マウスホイールでのズームと、中ボタン or Shift+左ドラッグでのパン機能を追加したQGraphicsView """
    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self.main_window = parent
        self.setAcceptDrops(False)

        self._is_panning = False
        self._pan_start_pos = QPointF()
        self._pan_start_scroll_h = 0
        self._pan_start_scroll_v = 0

    def wheelEvent(self, event):
        """ マウスホイールを回した際のイベント """
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            zoom_in_factor = 1.1
            zoom_out_factor = 1 / zoom_in_factor

            transform = self.transform()
            current_scale = transform.m11()
            min_scale, max_scale = 0.05, 20.0

            if event.angleDelta().y() > 0:
                if max_scale > current_scale:
                    self.scale(zoom_in_factor, zoom_in_factor)
            else:
                if min_scale < current_scale:
                    self.scale(zoom_out_factor, zoom_out_factor)
            
            event.accept() 
        else:
            super().wheelEvent(event)

    def mousePressEvent(self, event):
        """ 中ボタン or Shift+左ボタンが押されたらパン（視点移動）モードを開始 """
        is_middle_button = event.button() == Qt.MouseButton.MiddleButton
        is_shift_left_button = (event.button() == Qt.MouseButton.LeftButton and
                                event.modifiers() & Qt.KeyboardModifier.ShiftModifier)

        if is_middle_button or is_shift_left_button:
            self._is_panning = True
            self._pan_start_pos = event.pos() # ビューポート座標で開始点を記録
            # 現在のスクロールバーの位置を記録
            self._pan_start_scroll_h = self.horizontalScrollBar().value()
            self._pan_start_scroll_v = self.verticalScrollBar().value()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        """ パンモード中にマウスを動かしたら、スクロールバーを操作して視点を移動させる """
        if self._is_panning:
            delta = event.pos() - self._pan_start_pos # マウスの移動量を計算
            # 開始時のスクロール位置から移動量を引いた値を新しいスクロール位置に設定
            self.horizontalScrollBar().setValue(self._pan_start_scroll_h - delta.x())
            self.verticalScrollBar().setValue(self._pan_start_scroll_v - delta.y())
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        """ パンに使用したボタンが離されたらパンモードを終了 """
        # パンを開始したボタン（中 or 左）のどちらかが離されたかをチェック
        is_middle_button_release = event.button() == Qt.MouseButton.MiddleButton
        is_left_button_release = event.button() == Qt.MouseButton.LeftButton

        if self._is_panning and (is_middle_button_release or is_left_button_release):
            self._is_panning = False
            # 現在の描画モードに応じたカーソルに戻す
            current_mode = self.scene().mode if self.scene() else 'select'
            if current_mode == 'select':
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif current_mode.startswith(('atom', 'bond', 'template')):
                self.setCursor(Qt.CursorShape.CrossCursor)
            elif current_mode.startswith(('charge', 'radical')):
                self.setCursor(Qt.CursorShape.CrossCursor)
            else:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
        else:
            super().mouseReleaseEvent(event)

class CalculationWorker(QObject):
    status_update = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(object)  # emit (worker_id, msg) tuples for robustness
    # Per-worker start signal to avoid sharing a single MainWindow signal
    # among many worker instances (which causes race conditions and stale
    # workers being started on a single emission).
    start_work = pyqtSignal(str, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Connect the worker's own start signal to its run slot. This
        # guarantees that only this worker will respond when start_work
        # is emitted (prevents cross-talk between workers).
        try:
            self.start_work.connect(self.run_calculation)
        except Exception:
            # Be defensive: if connection fails, continue; the caller may
            # fallback to emitting directly.
            pass

    @pyqtSlot(str, object)
    def run_calculation(self, mol_block, options=None):
        try:
            # The worker may be asked to halt via a shared set `halt_ids` and
            # identifies its own run by options['worker_id'] (int).
            worker_id = None
            try:
                worker_id = options.get('worker_id') if options else None
            except Exception:
                worker_id = None

            # If a caller starts a worker without providing a worker_id, treat
            # it as a "global" worker that can still be halted via a global
            # halt flag. Emit a single status warning so callers know that
            # the worker was started without an identifier.
            _warned_no_worker_id = False
            if worker_id is None:
                try:
                    # best-effort, swallow any errors (signals may not be connected)
                    self.status_update.emit("Warning: worker started without 'worker_id'; will listen for global halt signals.")
                except Exception:
                    pass
                _warned_no_worker_id = True

            def _check_halted():
                try:
                    halt_ids = getattr(self, 'halt_ids', None)
                    # If worker_id is None, allow halting via a global mechanism:
                    #  - an explicit attribute `halt_all` set to True on the worker
                    #  - the shared `halt_ids` set containing None or the sentinel 'ALL'
                    if worker_id is None:
                        if getattr(self, 'halt_all', False):
                            return True
                        if halt_ids is None:
                            return False
                        # Support both None-in-set and string sentinel for compatibility
                        return (None in halt_ids) or ('ALL' in halt_ids)

                    if halt_ids is None:
                        return False
                    return (worker_id in halt_ids)
                except Exception:
                    return False

            # Safe-emission helpers: do nothing if this worker has been halted.
            def _safe_status(msg):
                try:
                    if _check_halted():
                        return
                    self.status_update.emit(msg)
                except Exception:
                    # Swallow any signal-emission errors to avoid crashing the worker
                    pass

            def _safe_finished(payload):
                try:
                    # Attempt to emit the payload; preserve existing fallback behavior
                    try:
                        self.finished.emit(payload)
                    except TypeError:
                        # Some slots/old code may expect a single-molecule arg; try that too
                        try:
                            # If payload was a tuple like (worker_id, mol), try sending the second element
                            if isinstance(payload, (list, tuple)) and len(payload) >= 2:
                                self.finished.emit(payload[1])
                            else:
                                self.finished.emit(payload)
                        except Exception:
                            pass
                except Exception:
                    pass

            def _safe_error(msg):
                try:
                    # Emit a tuple containing the worker_id (may be None) and the message
                    try:
                        self.error.emit((worker_id, msg))
                    except Exception:
                        # Fallback to emitting the raw message if tuple emission fails for any reason
                        try:
                            self.error.emit(msg)
                        except Exception:
                            pass
                except Exception:
                    pass

            # options: dict-like with keys: 'conversion_mode' -> 'fallback'|'rdkit'|'obabel'|'direct'
            if options is None:
                options = {}
            conversion_mode = options.get('conversion_mode', 'fallback')
            # Ensure params exists in all code paths (some RDKit calls below
            # reference `params` and earlier editing introduced a path where
            # it might not be defined). Initialize to None here and assign
            # a proper ETKDG params object later where needed.
            params = None
            if not mol_block:
                raise ValueError("No atoms to convert.")
            
            _safe_status("Creating 3D structure...")

            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None:
                raise ValueError("Failed to create molecule from MOL block.")

            # Check early whether this run has been requested to halt
            if _check_halted():
                raise RuntimeError("Halted")

            # CRITICAL FIX: Extract and restore explicit E/Z labels from MOL block
            # Parse M CFG lines to get explicit stereo labels
            explicit_stereo = {}
            mol_lines = mol_block.split('\n')
            for line in mol_lines:
                if line.startswith('M  CFG'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            bond_idx = int(parts[3]) - 1  # MOL format is 1-indexed
                            cfg_value = int(parts[4])
                            # cfg_value: 1=Z, 2=E in MOL format
                            if cfg_value == 1:
                                explicit_stereo[bond_idx] = Chem.BondStereo.STEREOZ
                            elif cfg_value == 2:
                                explicit_stereo[bond_idx] = Chem.BondStereo.STEREOE
                        except (ValueError, IndexError):
                            continue

            # Force explicit stereo labels regardless of coordinates
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Find suitable stereo atoms
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        
                        # Pick heavy atom neighbors preferentially
                        begin_neighbors = [nbr for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_atom.GetIdx()]
                        end_neighbors = [nbr for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_atom.GetIdx()]
                        
                        if begin_neighbors and end_neighbors:
                            # Prefer heavy atoms
                            begin_heavy = [n for n in begin_neighbors if n.GetAtomicNum() > 1]
                            end_heavy = [n for n in end_neighbors if n.GetAtomicNum() > 1]
                            
                            stereo_atom1 = (begin_heavy[0] if begin_heavy else begin_neighbors[0]).GetIdx()
                            stereo_atom2 = (end_heavy[0] if end_heavy else end_neighbors[0]).GetIdx()
                            
                            bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                            bond.SetStereo(stereo_type)

            # Do NOT call AssignStereochemistry here as it overrides our explicit labels

            mol = Chem.AddHs(mol)

            # Check after adding Hs (may be a long operation)
            if _check_halted():
                raise RuntimeError("Halted")

            # CRITICAL: Re-apply explicit stereo after AddHs which may renumber atoms
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Re-find suitable stereo atoms after hydrogen addition
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        
                        # Pick heavy atom neighbors preferentially
                        begin_neighbors = [nbr for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_atom.GetIdx()]
                        end_neighbors = [nbr for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_atom.GetIdx()]
                        
                        if begin_neighbors and end_neighbors:
                            # Prefer heavy atoms
                            begin_heavy = [n for n in begin_neighbors if n.GetAtomicNum() > 1]
                            end_heavy = [n for n in end_neighbors if n.GetAtomicNum() > 1]
                            
                            stereo_atom1 = (begin_heavy[0] if begin_heavy else begin_neighbors[0]).GetIdx()
                            stereo_atom2 = (end_heavy[0] if end_heavy else end_neighbors[0]).GetIdx()
                            
                            bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                            bond.SetStereo(stereo_type)

            # Direct mode: construct a 3D conformer without embedding by using
            # the 2D coordinates from the MOL block (z=0) and placing added
            # hydrogens close to their parent heavy atoms with a small z offset.
            # This avoids 3D embedding entirely and is useful for quick viewing
            # when stereochemistry/geometry refinement is not desired.
            if conversion_mode == 'direct':
                _safe_status("Direct conversion: using 2D coordinates + adding missing H (no embedding).")
                try:
                    # 1) Parse MOL block *with* existing hydrogens (removeHs=False)
                    #    to get coordinates for *all existing* atoms.
                    parsed_coords = []  # all-atom coordinates (x, y, z)
                    stereo_dirs = []    # list of (begin_idx, end_idx, stereo_flag)
                    
                    base2d_all = None
                    try:
                        # H原子を含めてパース
                        base2d_all = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
                    except Exception:
                        try:
                            base2d_all = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
                        except Exception:
                            base2d_all = None

                    if base2d_all is not None and base2d_all.GetNumConformers() > 0:
                        oconf = base2d_all.GetConformer()
                        for i in range(base2d_all.GetNumAtoms()):
                            p = oconf.GetAtomPosition(i)
                            parsed_coords.append((float(p.x), float(p.y), 0.0))
                    
                    # 2) Parse wedge/dash bond information (using all atoms)
                    try:
                        lines = mol_block.splitlines()
                        counts_idx = None
                        
                        for i, ln in enumerate(lines[:40]):
                            if re.match(r"^\s*\d+\s+\d+", ln):
                                counts_idx = i
                                break
                        
                        if counts_idx is not None:
                            parts = lines[counts_idx].split()
                            try:
                                natoms = int(parts[0])
                                nbonds = int(parts[1])
                            except Exception:
                                natoms = nbonds = 0
                            
                            # 全原子マップ (MOL 1-based index -> 0-based index)
                            atom_map = {i + 1: i for i in range(natoms)}
                            
                            bond_start = counts_idx + 1 + natoms
                            for j in range(min(nbonds, max(0, len(lines) - bond_start))):
                                bond_line = lines[bond_start + j]
                                try:
                                    m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)(?:\s+(-?\d+))?", bond_line)
                                    if m:
                                        try:
                                            atom1_mol = int(m.group(1))  # 1-based MOL index
                                            atom2_mol = int(m.group(2))  # 1-based MOL index
                                        except Exception:
                                            continue
                                        try:
                                            stereo_raw = int(m.group(4)) if m.group(4) is not None else 0
                                        except Exception:
                                            stereo_raw = 0
                                    else:
                                        fields = bond_line.split()
                                        if len(fields) >= 4:
                                            try:
                                                atom1_mol = int(fields[0])  # 1-based MOL index
                                                atom2_mol = int(fields[1])  # 1-based MOL index
                                            except Exception:
                                                continue
                                            try:
                                                stereo_raw = int(fields[3]) if len(fields) > 3 else 0
                                            except Exception:
                                                stereo_raw = 0
                                        else:
                                            continue

                                    # V2000の立体表記を正規化
                                    if stereo_raw == 1:
                                        stereo_flag = 1 # Wedge
                                    elif stereo_raw == 2:
                                        stereo_flag = 6 # Dash (V2000では 6 がDash)
                                    else:
                                        stereo_flag = stereo_raw

                                    # 全原子マップでチェック
                                    if atom1_mol in atom_map and atom2_mol in atom_map:
                                        idx1 = atom_map[atom1_mol]
                                        idx2 = atom_map[atom2_mol]
                                        if stereo_flag in (1, 6): # Wedge (1) or Dash (6)
                                            stereo_dirs.append((idx1, idx2, stereo_flag))
                                except Exception:
                                    continue
                    except Exception:
                        stereo_dirs = []
                
                    # Fallback for parsed_coords (if RDKit parse failed)
                    if not parsed_coords:
                        try:
                            lines = mol_block.splitlines()
                            counts_idx = None
                            for i, ln in enumerate(lines[:40]):
                                if re.match(r"^\s*\d+\s+\d+", ln):
                                    counts_idx = i
                                    break
                            if counts_idx is not None:
                                parts = lines[counts_idx].split()
                                try:
                                    natoms = int(parts[0])
                                except Exception:
                                    natoms = 0
                                atom_start = counts_idx + 1
                                for j in range(min(natoms, max(0, len(lines) - atom_start))):
                                    atom_line = lines[atom_start + j]
                                    try:
                                        x = float(atom_line[0:10].strip()); y = float(atom_line[10:20].strip()); z = float(atom_line[20:30].strip())
                                    except Exception:
                                        fields = atom_line.split()
                                        if len(fields) >= 4:
                                            try:
                                                x = float(fields[0]); y = float(fields[1]); z = float(fields[2])
                                            except Exception:
                                                continue
                                        else:
                                            continue
                                    # H原子もスキップしない
                                    parsed_coords.append((x, y, z))
                        except Exception:
                            parsed_coords = []
                    
                    if not parsed_coords:
                        raise ValueError("Failed to parse coordinates from MOL block for direct conversion.")

                    # 3) `mol` は既に AddHs された状態
                    #    元の原子数 (H含む) を parsed_coords の長さから取得
                    num_existing_atoms = len(parsed_coords)

                    # 4) コンフォーマを作成
                    conf = Chem.Conformer(mol.GetNumAtoms())

                    for i in range(mol.GetNumAtoms()):
                        if i < num_existing_atoms:
                            # 既存原子 (H含む): 2D座標 (z=0) を設定
                            x, y, z_ignored = parsed_coords[i]
                            try:
                                conf.SetAtomPosition(i, rdGeometry.Point3D(float(x), float(y), 0.0))
                            except Exception:
                                pass
                        else:
                            # 新規追加されたH原子: 親原子の近くに配置
                            atom = mol.GetAtomWithIdx(i)
                            if atom.GetAtomicNum() == 1:
                                neighs = [n for n in atom.GetNeighbors() if n.GetIdx() < num_existing_atoms]
                                heavy_pos_found = False
                                for nb in neighs: # 親原子 (重原子または既存H)
                                    try:
                                        nb_idx = nb.GetIdx()
                                        # if nb_idx < num_existing_atoms: # チェックは不要 (neighs で既にフィルタ済み)
                                        nbpos = conf.GetAtomPosition(nb_idx)
                                        # Geometry-based placement:
                                        # Compute an "empty" direction around the parent atom by
                                        # summing existing bond unit vectors and taking the
                                        # opposite. If degenerate, pick a perpendicular or
                                        # fallback vector. Rotate slightly if multiple Hs already
                                        # attached to avoid overlap.
                                        parent_idx = nb_idx
                                        try:
                                            parent_pos = conf.GetAtomPosition(parent_idx)
                                            parent_atom = mol.GetAtomWithIdx(parent_idx)
                                            # collect unit vectors to already-placed neighbors (idx < i)
                                            vecs = []
                                            for nbr in parent_atom.GetNeighbors():
                                                nidx = nbr.GetIdx()
                                                if nidx == i:
                                                    continue
                                                # only consider neighbors whose positions are already set
                                                if nidx < i:
                                                    try:
                                                        p = conf.GetAtomPosition(nidx)
                                                        vx = float(p.x) - float(parent_pos.x)
                                                        vy = float(p.y) - float(parent_pos.y)
                                                        nrm = math.hypot(vx, vy)
                                                        if nrm > 1e-6:
                                                            vecs.append((vx / nrm, vy / nrm))
                                                    except Exception:
                                                        continue

                                            if vecs:
                                                sx = sum(v[0] for v in vecs)
                                                sy = sum(v[1] for v in vecs)
                                                fx = -sx
                                                fy = -sy
                                                fn = math.hypot(fx, fy)
                                                if fn < 1e-6:
                                                    # degenerate: pick a perpendicular to first bond
                                                    fx = -vecs[0][1]
                                                    fy = vecs[0][0]
                                                    fn = math.hypot(fx, fy)
                                                fx /= fn; fy /= fn

                                                # Avoid placing multiple Hs at identical directions
                                                existing_h_count = sum(1 for nbr in parent_atom.GetNeighbors()
                                                                       if nbr.GetIdx() < i and nbr.GetAtomicNum() == 1)
                                                angle = existing_h_count * (math.pi / 6.0)  # 30deg steps
                                                cos_a = math.cos(angle); sin_a = math.sin(angle)
                                                rx = fx * cos_a - fy * sin_a
                                                ry = fx * sin_a + fy * cos_a

                                                bond_length = 1.0
                                                conf.SetAtomPosition(i, rdGeometry.Point3D(
                                                    float(parent_pos.x) + rx * bond_length,
                                                    float(parent_pos.y) + ry * bond_length,
                                                    0.3
                                                ))
                                            else:
                                                # No existing placed neighbors: fallback to small offset
                                                conf.SetAtomPosition(i, rdGeometry.Point3D(
                                                    float(parent_pos.x) + 0.5,
                                                    float(parent_pos.y) + 0.5,
                                                    0.3
                                                ))

                                            heavy_pos_found = True
                                            break
                                        except Exception:
                                            # fall back to trying the next neighbor if any
                                            continue
                                    except Exception:
                                        continue
                                if not heavy_pos_found:
                                    # フォールバック (原点近く)
                                    try:
                                        conf.SetAtomPosition(i, rdGeometry.Point3D(0.0, 0.0, 0.10))
                                    except Exception:
                                        pass
                    
                    # 5) Wedge/Dash の Zオフセットを適用
                    try:
                        stereo_z_offset = 1.5  # wedge -> +1.5, dash -> -1.5
                        for begin_idx, end_idx, stereo_flag in stereo_dirs:
                            try:
                                # インデックスは既存原子内のはず
                                if begin_idx >= num_existing_atoms or end_idx >= num_existing_atoms:
                                    continue
                                    
                                if stereo_flag not in (1, 6):
                                    continue
                                
                                sign = 1.0 if stereo_flag == 1 else -1.0
                                
                                # end_idx (立体表記の終点側原子) にZオフセットを適用
                                pos = conf.GetAtomPosition(end_idx)
                                newz = float(pos.z) + (stereo_z_offset * sign) # 既存のZ=0にオフセットを加算
                                conf.SetAtomPosition(end_idx, rdGeometry.Point3D(float(pos.x), float(pos.y), float(newz)))
                            except Exception:
                                continue
                    except Exception:
                        pass
                    
                    # コンフォーマを入れ替えて終了
                    try:
                        mol.RemoveAllConformers()
                    except Exception:
                        pass
                    mol.AddConformer(conf, assignId=True)
                    
                    if _check_halted():
                        raise RuntimeError("Halted (after optimization)")
                    try:
                        _safe_finished((worker_id, mol))
                    except Exception:
                        _safe_finished(mol)
                    _safe_status("Direct conversion completed.")
                    return
                except Exception as e:
                    _safe_status(f"Direct conversion failed: {e}")

            params = AllChem.ETKDGv2()
            params.randomSeed = 42
            # CRITICAL: Force ETKDG to respect the existing stereochemistry
            params.useExpTorsionAnglePrefs = True
            params.useBasicKnowledge = True
            params.enforceChirality = True  # This is critical for stereo preservation
            
            # Store original stereochemistry before embedding (prioritizing explicit labels)
            original_stereo_info = []
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        stereo_atoms = bond.GetStereoAtoms()
                        original_stereo_info.append((bond.GetIdx(), stereo_type, stereo_atoms))
            
            # Also store any other stereo bonds not in explicit_stereo
            for bond in mol.GetBonds():
                if (bond.GetBondType() == Chem.BondType.DOUBLE and 
                    bond.GetStereo() != Chem.BondStereo.STEREONONE and
                    bond.GetIdx() not in explicit_stereo):
                    stereo_atoms = bond.GetStereoAtoms()
                    original_stereo_info.append((bond.GetIdx(), bond.GetStereo(), stereo_atoms))
            
            # Only report RDKit-specific messages when RDKit embedding will be
            # attempted. For other conversion modes, emit clearer, non-misleading
            # status messages so the UI doesn't show "RDKit" when e.g. direct
            # coordinates or Open Babel will be used.
            if conversion_mode in ('fallback', 'rdkit'):
                _safe_status("RDKit: Embedding 3D coordinates...")
            elif conversion_mode == 'obabel':
                pass
            else:
                # direct mode (or any other explicit non-RDKit mode)
                pass
            if _check_halted():
                raise RuntimeError("Halted")
            
            # Try multiple times with different approaches if needed
            conf_id = -1
            
            # First attempt: Standard ETKDG with stereo enforcement
            try:
                # Only attempt RDKit embedding if mode allows
                if conversion_mode in ('fallback', 'rdkit'):
                    conf_id = AllChem.EmbedMolecule(mol, params)
                else:
                    conf_id = -1
                # Final check before returning success
                if _check_halted():
                    raise RuntimeError("Halted")
            except Exception as e:
                # Standard embedding failed; report and continue to fallback attempts
                _safe_status(f"Standard embedding failed: {e}")

                # Second attempt: Use constraint embedding if available (only when RDKit is allowed)
                if conf_id == -1 and conversion_mode in ('fallback', 'rdkit'):
                    try:
                        # Create distance constraints for double bonds to enforce E/Z geometry
                        bounds_matrix = AllChem.GetMoleculeBoundsMatrix(mol)

                        # Add constraints for E/Z bonds
                        for bond_idx, stereo, stereo_atoms in original_stereo_info:
                            bond = mol.GetBondWithIdx(bond_idx)
                            if len(stereo_atoms) == 2:
                                atom1_idx = bond.GetBeginAtomIdx()
                                atom2_idx = bond.GetEndAtomIdx()
                                neighbor1_idx = stereo_atoms[0]
                                neighbor2_idx = stereo_atoms[1]

                                # For Z (cis): neighbors should be closer
                                # For E (trans): neighbors should be farther
                                if stereo == Chem.BondStereo.STEREOZ:
                                    # Z configuration: set shorter distance constraint
                                    target_dist = 3.0  # Angstroms
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = min(bounds_matrix[neighbor1_idx][neighbor2_idx], target_dist)
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = min(bounds_matrix[neighbor2_idx][neighbor1_idx], target_dist)
                                elif stereo == Chem.BondStereo.STEREOE:
                                    # E configuration: set longer distance constraint  
                                    target_dist = 5.0  # Angstroms
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = max(bounds_matrix[neighbor1_idx][neighbor2_idx], target_dist)
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = max(bounds_matrix[neighbor2_idx][neighbor1_idx], target_dist)

                        DoTriangleSmoothing(bounds_matrix)
                        conf_id = AllChem.EmbedMolecule(mol, bounds_matrix, params)
                        _safe_status("Constraint-based embedding succeeded")
                    except Exception:
                        # Constraint embedding failed: only raise error if mode is 'rdkit', otherwise allow fallback
                        _safe_status("RDKit: Constraint embedding failed")
                        if conversion_mode == 'rdkit':
                            raise RuntimeError("RDKit: Constraint embedding failed")
                        conf_id = -1
                    
            # Fallback: Try basic embedding
            if conf_id == -1:
                try:
                    if conversion_mode in ('fallback', 'rdkit'):
                        basic_params = AllChem.ETKDGv2()
                        basic_params.randomSeed = 42
                        conf_id = AllChem.EmbedMolecule(mol, basic_params)
                    else:
                        conf_id = -1
                except Exception:
                    pass
            '''
            if conf_id == -1:
                        _safe_status("Initial embedding failed, retrying with ignoreSmoothingFailures=True...")
                # Try again with ignoreSmoothingFailures instead of random-seed retries
                params.ignoreSmoothingFailures = True
                # Use a deterministic seed to avoid random-coordinate behavior here
                params.randomSeed = 0
                conf_id = AllChem.EmbedMolecule(mol, params)

            if conf_id == -1:
                self.status_update.emit("Random-seed retry failed, attempting with random coordinates...")
                try:
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True, ignoreSmoothingFailures=True)
                except TypeError:
                    # Some RDKit versions expect useRandomCoords in params
                    params.useRandomCoords = True
                    conf_id = AllChem.EmbedMolecule(mol, params)
            '''

            # Determine requested MMFF variant from options (fall back to MMFF94s)
            opt_method = None
            try:
                opt_method = options.get('optimization_method') if options else None
            except Exception:
                opt_method = None

            if conf_id != -1:
                # Success with RDKit: optimize and finish
                # CRITICAL: Restore original stereochemistry after embedding (explicit labels first)
                for bond_idx, stereo, stereo_atoms in original_stereo_info:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if len(stereo_atoms) == 2:
                        bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                    bond.SetStereo(stereo)
                
                try:
                    mmff_variant = "MMFF94s"
                    if opt_method and str(opt_method).upper() == 'MMFF94_RDKIT':
                        mmff_variant = "MMFF94"
                    if _check_halted():
                        raise RuntimeError("Halted")
                    AllChem.MMFFOptimizeMolecule(mol, mmffVariant=mmff_variant)
                except Exception:
                    # fallback to UFF if MMFF fails
                    try:
                        if _check_halted():
                            raise RuntimeError("Halted")
                        AllChem.UFFOptimizeMolecule(mol)
                    except Exception:
                        pass
                
                # CRITICAL: Restore stereochemistry again after optimization (explicit labels priority)
                for bond_idx, stereo, stereo_atoms in original_stereo_info:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if len(stereo_atoms) == 2:
                        bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                    bond.SetStereo(stereo)
                
                # Do NOT call AssignStereochemistry here as it would override our explicit labels
                # Include worker_id so the main thread can ignore stale results
                # CRITICAL: Check for halt *before* emitting finished signal
                if _check_halted():
                    raise RuntimeError("Halted (after optimization)")
                try:
                    _safe_finished((worker_id, mol))
                except Exception:
                    # Fallback to legacy single-arg emit
                    _safe_finished(mol)
                _safe_status("RDKit 3D conversion succeeded.")
                return

            # If RDKit did not produce a conf and OBabel is allowed, try Open Babel
            if conf_id == -1 and conversion_mode in ('fallback', 'obabel'):
                _safe_status("RDKit embedding failed or disabled. Attempting Open Babel...")
                try:
                    if not OBABEL_AVAILABLE:
                        raise RuntimeError("Open Babel (pybel) is not available in this Python environment.")
                    ob_mol = pybel.readstring("mol", mol_block)
                    try:
                        ob_mol.addh()
                    except Exception:
                        pass
                    ob_mol.make3D()
                    try:
                        _safe_status("Optimizing with Open Babel (MMFF94)...")
                        if _check_halted():
                            raise RuntimeError("Halted")
                        ob_mol.localopt(forcefield='mmff94', steps=500)
                    except Exception:
                        try:
                            _safe_status("MMFF94 failed, falling back to UFF...")
                            if _check_halted():
                                raise RuntimeError("Halted")
                            ob_mol.localopt(forcefield='uff', steps=500)
                        except Exception:
                            _safe_status("UFF optimization also failed.")
                            pass
                    molblock_ob = ob_mol.write("mol")
                    rd_mol = Chem.MolFromMolBlock(molblock_ob, removeHs=False)
                    if rd_mol is None:
                        raise ValueError("Open Babel produced invalid MOL block.")
                    rd_mol = Chem.AddHs(rd_mol)
                    try:
                        mmff_variant = "MMFF94s"
                        if opt_method and str(opt_method).upper() == 'MMFF94_RDKIT':
                            mmff_variant = "MMFF94"
                        if _check_halted():
                            raise RuntimeError("Halted")
                        AllChem.MMFFOptimizeMolecule(rd_mol, mmffVariant=mmff_variant)
                    except Exception:
                        try:
                            if _check_halted():
                                raise RuntimeError("Halted")
                            AllChem.UFFOptimizeMolecule(rd_mol)
                        except Exception:
                            pass
                    _safe_status("Open Babel embedding succeeded. Warning: Conformation accuracy may be limited.")
                    # CRITICAL: Check for halt *before* emitting finished signal
                    if _check_halted():
                        raise RuntimeError("Halted (after optimization)")
                    try:
                        _safe_finished((worker_id, rd_mol))
                    except Exception:
                        _safe_finished(rd_mol)
                    return
                except Exception as ob_err:
                    raise RuntimeError(f"Open Babel 3D conversion failed: {ob_err}")

            if conf_id == -1 and conversion_mode == 'rdkit':
                raise RuntimeError("RDKit 3D conversion failed (rdkit-only mode)")

        except Exception as e:
            _safe_error(str(e))


class PeriodicTableDialog(QDialog):
    element_selected = pyqtSignal(str)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select an Element")
        layout = QGridLayout(self)
        self.setLayout(layout)

        elements = [
            ('H',1,1), ('He',1,18),
            ('Li',2,1), ('Be',2,2), ('B',2,13), ('C',2,14), ('N',2,15), ('O',2,16), ('F',2,17), ('Ne',2,18),
            ('Na',3,1), ('Mg',3,2), ('Al',3,13), ('Si',3,14), ('P',3,15), ('S',3,16), ('Cl',3,17), ('Ar',3,18),
            ('K',4,1), ('Ca',4,2), ('Sc',4,3), ('Ti',4,4), ('V',4,5), ('Cr',4,6), ('Mn',4,7), ('Fe',4,8),
            ('Co',4,9), ('Ni',4,10), ('Cu',4,11), ('Zn',4,12), ('Ga',4,13), ('Ge',4,14), ('As',4,15), ('Se',4,16),
            ('Br',4,17), ('Kr',4,18),
            ('Rb',5,1), ('Sr',5,2), ('Y',5,3), ('Zr',5,4), ('Nb',5,5), ('Mo',5,6), ('Tc',5,7), ('Ru',5,8),
            ('Rh',5,9), ('Pd',5,10), ('Ag',5,11), ('Cd',5,12), ('In',5,13), ('Sn',5,14), ('Sb',5,15), ('Te',5,16),
            ('I',5,17), ('Xe',5,18),
            ('Cs',6,1), ('Ba',6,2), ('Hf',6,4), ('Ta',6,5), ('W',6,6), ('Re',6,7), ('Os',6,8),
            ('Ir',6,9), ('Pt',6,10), ('Au',6,11), ('Hg',6,12), ('Tl',6,13), ('Pb',6,14), ('Bi',6,15), ('Po',6,16),
            ('At',6,17), ('Rn',6,18),
            ('Fr',7,1), ('Ra',7,2), ('Rf',7,4), ('Db',7,5), ('Sg',7,6), ('Bh',7,7), ('Hs',7,8),
            ('Mt',7,9), ('Ds',7,10), ('Rg',7,11), ('Cn',7,12), ('Nh',7,13), ('Fl',7,14), ('Mc',7,15), ('Lv',7,16),
            ('Ts',7,17), ('Og',7,18),
            # Lanthanides (placed on a separate row)
            ('La',8,3), ('Ce',8,4), ('Pr',8,5), ('Nd',8,6), ('Pm',8,7), ('Sm',8,8), ('Eu',8,9), ('Gd',8,10), ('Tb',8,11),
            ('Dy',8,12), ('Ho',8,13), ('Er',8,14), ('Tm',8,15), ('Yb',8,16), ('Lu',8,17),
            # Actinides (separate row)
            ('Ac',9,3), ('Th',9,4), ('Pa',9,5), ('U',9,6), ('Np',9,7), ('Pu',9,8), ('Am',9,9), ('Cm',9,10), ('Bk',9,11),
            ('Cf',9,12), ('Es',9,13), ('Fm',9,14), ('Md',9,15), ('No',9,16), ('Lr',9,17),
        ]
        for symbol, row, col in elements:
            b = QPushButton(symbol)
            b.setFixedSize(40,40)

            # Prefer saved user override (from parent.settings), otherwise use CPK_COLORS
            try:
                overrides = parent.settings.get('cpk_colors', {}) if parent and hasattr(parent, 'settings') else {}
                override = overrides.get(symbol)
            except Exception:
                override = None
            q_color = QColor(override) if override else CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT'])

            # 背景色の輝度を計算して、文字色を黒か白に決定
            # 輝度 = (R*299 + G*587 + B*114) / 1000
            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = "white" if brightness < 128 else "black"

            # ボタンのスタイルシートを設定
            b.setStyleSheet(
                f"background-color: {q_color.name()};"
                f"color: {text_color};"
                "border: 1px solid #555;"
                "font-weight: bold;"
            )

            b.clicked.connect(self.on_button_clicked)
            layout.addWidget(b, row, col)

    def on_button_clicked(self):
        b=self.sender()
        self.element_selected.emit(b.text())
        self.accept()


class ColorSettingsDialog(QDialog):
    """Dialog to customize CPK element colors.

    - Click an element to pick a new color for the element (CPK colors).
    - Reset All button to restore defaults for everything.
    """
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CPK Colors")
        self.parent_window = parent
        self.current_settings = current_settings or {}

        self.changed_cpk = {}  # symbol -> hex
        self._reset_all_flag = False

        layout = QVBoxLayout(self)

        # Color picking for CPK is available in the periodic table and CPK dialog

        # Periodic table grid (buttons like PeriodicTableDialog)
        grid = QGridLayout()
        self.element_buttons = {}
        elements = [
            ('H',1,1), ('He',1,18),
            ('Li',2,1), ('Be',2,2), ('B',2,13), ('C',2,14), ('N',2,15), ('O',2,16), ('F',2,17), ('Ne',2,18),
            ('Na',3,1), ('Mg',3,2), ('Al',3,13), ('Si',3,14), ('P',3,15), ('S',3,16), ('Cl',3,17), ('Ar',3,18),
            ('K',4,1), ('Ca',4,2), ('Sc',4,3), ('Ti',4,4), ('V',4,5), ('Cr',4,6), ('Mn',4,7), ('Fe',4,8),
            ('Co',4,9), ('Ni',4,10), ('Cu',4,11), ('Zn',4,12), ('Ga',4,13), ('Ge',4,14), ('As',4,15), ('Se',4,16),
            ('Br',4,17), ('Kr',4,18),
            ('Rb',5,1), ('Sr',5,2), ('Y',5,3), ('Zr',5,4), ('Nb',5,5), ('Mo',5,6), ('Tc',5,7), ('Ru',5,8),
            ('Rh',5,9), ('Pd',5,10), ('Ag',5,11), ('Cd',5,12), ('In',5,13), ('Sn',5,14), ('Sb',5,15), ('Te',5,16),
            ('I',5,17), ('Xe',5,18),
            ('Cs',6,1), ('Ba',6,2), ('Hf',6,4), ('Ta',6,5), ('W',6,6), ('Re',6,7), ('Os',6,8),
            ('Ir',6,9), ('Pt',6,10), ('Au',6,11), ('Hg',6,12), ('Tl',6,13), ('Pb',6,14), ('Bi',6,15), ('Po',6,16),
            ('At',6,17), ('Rn',6,18),
            ('Fr',7,1), ('Ra',7,2), ('Rf',7,4), ('Db',7,5), ('Sg',7,6), ('Bh',7,7), ('Hs',7,8),
            ('Mt',7,9), ('Ds',7,10), ('Rg',7,11), ('Cn',7,12), ('Nh',7,13), ('Fl',7,14), ('Mc',7,15), ('Lv',7,16),
            ('Ts',7,17), ('Og',7,18),
            ('La',8,3), ('Ce',8,4), ('Pr',8,5), ('Nd',8,6), ('Pm',8,7), ('Sm',8,8), ('Eu',8,9), ('Gd',8,10), ('Tb',8,11),
            ('Dy',8,12), ('Ho',8,13), ('Er',8,14), ('Tm',8,15), ('Yb',8,16), ('Lu',8,17),
            ('Ac',9,3), ('Th',9,4), ('Pa',9,5), ('U',9,6), ('Np',9,7), ('Pu',9,8), ('Am',9,9), ('Cm',9,10), ('Bk',9,11),
            ('Cf',9,12), ('Es',9,13), ('Fm',9,14), ('Md',9,15), ('No',9,16), ('Lr',9,17),
        ]

        for symbol, row, col in elements:
            b = QPushButton(symbol)
            b.setFixedSize(40, 40)
            # Choose override color (if present) else default CPK color
            override = self.current_settings.get('cpk_colors', {}).get(symbol)
            if override:
                q_color = QColor(override)
            else:
                q_color = CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT'])

            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            b.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
            b.clicked.connect(self.on_element_clicked)
            grid.addWidget(b, row, col)
            self.element_buttons[symbol] = b

        layout.addLayout(grid)

        # Ball & Stick bond color (3D) picker - placed near the periodic table for CPK settings
        self.changed_bs_color = None
        try:
            bs_h = QHBoxLayout()
            bs_label = QLabel("Ball & Stick bond color:")
            self.bs_button = QPushButton()
            self.bs_button.setFixedSize(36, 24)
            # initialize from current settings (if provided)
            try:
                cur_bs = self.current_settings.get('ball_stick_bond_color') if self.current_settings else None
            except Exception:
                cur_bs = None
            if not cur_bs and self.parent_window and hasattr(self.parent_window, 'settings'):
                cur_bs = self.parent_window.settings.get('ball_stick_bond_color', '#7F7F7F')
            try:
                self.bs_button.setStyleSheet(f"background-color: {cur_bs}; border: 1px solid #888;")
                self.bs_button.setToolTip(cur_bs)
            except Exception:
                pass
            self.bs_button.clicked.connect(self.pick_bs_bond_color)
            bs_h.addWidget(bs_label)
            bs_h.addWidget(self.bs_button)
            bs_h.addStretch(1)
            layout.addLayout(bs_h)
        except Exception:
            pass

        # Reset button and action buttons
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

        h.addWidget(apply_button); h.addWidget(ok_button); h.addWidget(cancel_button)
        layout.addLayout(h)

        # initialize
        # No 2D bond color control here

    # 2D bond color picker removed — 2D bond color is fixed and not configurable here

    def on_element_clicked(self):
        b = self.sender()
        symbol = b.text()
        # get current color (override if exists else default)
        cur = self.current_settings.get('cpk_colors', {}).get(symbol)
        if not cur:
            cur = CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT']).name()
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.changed_cpk[symbol] = color.name()
            # Update button appearance
            brightness = (color.red() * 299 + color.green() * 587 + color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            b.setStyleSheet(f"background-color: {color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")

    def reset_all(self):
        # Clear overrides
        self.changed_cpk = {}
        self._reset_all_flag = True
        
        # 1. B&S結合色もリセット対象（デフォルト値）に設定
        try:
            self.changed_bs_color = self.parent_window.default_settings.get('ball_stick_bond_color', '#7F7F7F') if hasattr(self.parent_window, 'default_settings') else '#7F7F7F'
        except Exception:
            self.changed_bs_color = '#7F7F7F'
            
        # 2. ダイアログ内のCPKボタンの表示をデフォルトに戻す
        for s, btn in self.element_buttons.items():
            q_color = DEFAULT_CPK_COLORS.get(s, DEFAULT_CPK_COLORS['DEFAULT'])
            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            btn.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
        
        # 3. 3Dプレビューを更新する L.3337〜L.3386 の try...finally ブロックは削除

        # 4. ダイアログ内のB&S結合色ボタンの表示をデフォルトに戻す
        try:
            if hasattr(self, 'bs_button'):
                # self.changed_bs_color に設定したデフォルト値を反映
                hexv = self.changed_bs_color
                self.bs_button.setStyleSheet(f"background-color: {hexv}; border: 1px solid #888;")
                self.bs_button.setToolTip(hexv)
        except Exception:
            pass

    def apply_changes(self):
        # Persist only changed keys
        if self.parent_window:
            if self._reset_all_flag:
                # Remove any cpk overrides
                try:
                    if 'cpk_colors' in self.parent_window.settings:
                        del self.parent_window.settings['cpk_colors']
                except Exception:
                    pass
            if self.changed_cpk:
                # Merge with existing overrides
                cdict = self.parent_window.settings.get('cpk_colors', {}).copy()
                cdict.update(self.changed_cpk)
                self.parent_window.settings['cpk_colors'] = cdict
                self.parent_window.settings_dirty = True
            # After changing settings, update global CPK color map and refresh views
            try:
                self.parent_window.update_cpk_colors_from_settings()
            except Exception:
                pass
            try:
                self.parent_window.apply_3d_settings()
            except Exception:
                pass
            try:
                if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                    self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
            except Exception:
                pass
            # update 2D scene objects
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
            # update periodic table button styles in the dialog to reflect any overrides
            try:
                for s, btn in self.element_buttons.items():
                    try:
                        q_color = QColor(self.parent_window.settings.get('cpk_colors', {}).get(s, CPK_COLORS.get(s, CPK_COLORS['DEFAULT']).name()))
                        brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
                        text_color = 'white' if brightness < 128 else 'black'
                        btn.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
                    except Exception:
                        pass
            except Exception:
                pass
            # Refresh any open SettingsDialog instances so the ball & stick color preview updates
            try:
                for w in QApplication.topLevelWidgets():
                    try:
                        if isinstance(w, SettingsDialog):
                            try:
                                w.update_ui_from_settings(self.parent_window.settings)
                            except Exception:
                                pass
                    except Exception:
                        pass
            except Exception:
                pass
            # Persist changed Ball & Stick color if the user changed it from the CPK dialog
            if getattr(self, 'changed_bs_color', None):
                try:
                    self.parent_window.settings['ball_stick_bond_color'] = self.changed_bs_color
                    try:
                        self.parent_window.settings_dirty = True
                    except Exception:
                        pass
                    # After changing ball-stick color, ensure 3D view updates
                    try:
                        self.parent_window.apply_3d_settings()
                    except Exception:
                        pass
                    try:
                        if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                            self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
                    except Exception:
                        pass
                except Exception:
                    pass
            # Removed 2D bond color control — nothing to persist here
            elif self._reset_all_flag:
                # Reset Ball & Stick 3D bond color to default
                try:
                    # Use a stable default instead of relying on parent_window.default_settings
                    self.parent_window.settings['ball_stick_bond_color'] = '#7F7F7F'
                    try:
                        self.parent_window.settings_dirty = True
                    except Exception:
                        pass
                except Exception:
                    pass
                self.parent_window.update_cpk_colors_from_settings()
                self.parent_window.apply_3d_settings()
                if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                    self.parent_window.draw_molecule_3d(self.parent_window.current_mol)

            # update 2D scene
            try:
                if hasattr(self.parent_window, 'scene'):
                    for it in self.parent_window.scene.items():
                        try:
                            # AtomItem.update_style uses CPK_COLORS map
                            if hasattr(it, 'update_style'):
                                it.update_style()
                        except Exception:
                            pass
            except Exception:
                pass

    def accept(self):
        self.apply_changes()
        super().accept()

    def pick_bs_bond_color(self):
        """Pick Ball & Stick 3D bond color from the CPK dialog and update preview immediately."""
        try:
            cur = getattr(self, 'changed_bs_color', None) or (self.current_settings.get('ball_stick_bond_color') if self.current_settings else None)
        except Exception:
            cur = None
        if not cur and self.parent_window and hasattr(self.parent_window, 'settings'):
            cur = self.parent_window.settings.get('ball_stick_bond_color', '#7F7F7F')
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            hexv = color.name()
            self.changed_bs_color = hexv
            try:
                self.bs_button.setStyleSheet(f"background-color: {hexv}; border: 1px solid #888;")
                self.bs_button.setToolTip(hexv)
            except Exception:
                pass

# --- 最終版 AnalysisWindow クラス ---
class AnalysisWindow(QDialog):
    def __init__(self, mol, parent=None, is_xyz_derived=False):
        super().__init__(parent)
        self.mol = mol
        self.is_xyz_derived = is_xyz_derived  # XYZ由来かどうかのフラグ
        self.setWindowTitle("Molecule Analysis")
        self.setMinimumWidth(400)
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)
        grid_layout = QGridLayout()
        
        # --- 分子特性を計算 ---
        try:
            # RDKitのモジュールをインポート


            if self.is_xyz_derived:
                # XYZ由来の場合：元のXYZファイルの原子情報から直接計算
                # （結合推定の影響を受けない）
                
                # XYZファイルから読み込んだ元の原子情報を取得
                if hasattr(self.mol, '_xyz_atom_data'):
                    xyz_atoms = self.mol._xyz_atom_data
                else:
                    # フォールバック: RDKitオブジェクトから取得
                    xyz_atoms = [(atom.GetSymbol(), 0, 0, 0) for atom in self.mol.GetAtoms()]
                
                # 原子数と元素種を集計
                atom_counts = {}
                total_atoms = len(xyz_atoms)
                num_heavy_atoms = 0
                
                for symbol, x, y, z in xyz_atoms:
                    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                    if symbol != 'H':  # 水素以外
                        num_heavy_atoms += 1
                
                # 化学式を手動で構築（元素順序を考慮）
                element_order = ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
                formula_parts = []
                
                # 定義された順序で元素を追加
                remaining_counts = atom_counts.copy()
                for element in element_order:
                    if element in remaining_counts:
                        count = remaining_counts[element]
                        if count == 1:
                            formula_parts.append(element)
                        else:
                            formula_parts.append(f"{element}{count}")
                        del remaining_counts[element]
                
                # 残りの元素をアルファベット順で追加
                for element in sorted(remaining_counts.keys()):
                    count = remaining_counts[element]
                    if count == 1:
                        formula_parts.append(element)
                    else:
                        formula_parts.append(f"{element}{count}")
                
                mol_formula = ''.join(formula_parts)
                
                # 分子量と精密質量をRDKitから取得
                
                mol_wt = 0.0
                exact_mw = 0.0
                pt = Chem.GetPeriodicTable()
                
                for symbol, count in atom_counts.items():
                    try:
                        # RDKitの周期表から原子量と精密質量を取得
                        atomic_num = pt.GetAtomicNumber(symbol)
                        atomic_weight = pt.GetAtomicWeight(atomic_num)
                        exact_mass = pt.GetMostCommonIsotopeMass(atomic_num)
                        
                        mol_wt += atomic_weight * count
                        exact_mw += exact_mass * count
                    except (ValueError, RuntimeError):
                        # 認識されない元素の場合はスキップ
                        print(f"Warning: Unknown element {symbol}, skipping in mass calculation")
                        continue
                
                # 表示するプロパティを辞書にまとめる（XYZ元データから計算）
                properties = {
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Total Atoms:": str(total_atoms),
                }
                
                # 注意メッセージを追加
                note_label = QLabel("<i>Note: SMILES and structure-dependent properties are not available for XYZ-derived structures due to potential bond estimation inaccuracies.</i>")
                note_label.setWordWrap(True)
                main_layout.addWidget(note_label)
                
            else:
                # 通常の分子（MOLファイルや2Dエディタ由来）の場合：全てのプロパティを計算
                
                # SMILES生成用に、一時的に水素原子を取り除いた分子オブジェクトを作成
                mol_for_smiles = Chem.RemoveHs(self.mol)
                # 水素を取り除いた分子からSMILESを生成（常に簡潔な表記になる）
                smiles = Chem.MolToSmiles(mol_for_smiles, isomericSmiles=True)

                # 各種プロパティを計算
                mol_formula = rdMolDescriptors.CalcMolFormula(self.mol)
                mol_wt = Descriptors.MolWt(self.mol)
                exact_mw = Descriptors.ExactMolWt(self.mol)
                num_heavy_atoms = self.mol.GetNumHeavyAtoms()
                num_rings = rdMolDescriptors.CalcNumRings(self.mol)
                log_p = Descriptors.MolLogP(self.mol)
                tpsa = Descriptors.TPSA(self.mol)
                num_h_donors = rdMolDescriptors.CalcNumHBD(self.mol)
                num_h_acceptors = rdMolDescriptors.CalcNumHBA(self.mol)
                
                # InChIを生成
                try:
                    inchi = Chem.MolToInchi(self.mol)
                except Exception:
                    inchi = "N/A"

                # InChIKeyを生成（RDKitのinchi APIが無い場合に備えてフォールバック）
                try:
                    # Prefer Chem.MolToInchiKey when available
                    inchi_key = None
                    try:
                        inchi_key = Chem.MolToInchiKey(self.mol)
                    except Exception:
                        # Fallback to rdkit.Chem.inchi if present
                        try:
                            inchi_key = rd_inchi.MolToInchiKey(self.mol)
                        except Exception:
                            inchi_key = None

                    if not inchi_key:
                        inchi_key = "N/A"
                except Exception:
                    inchi_key = "N/A"

                # 表示するプロパティを辞書にまとめる
                properties = {
                    "SMILES:": smiles,
                    "InChI:": inchi,
                    "InChIKey:": inchi_key,
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Ring Count:": str(num_rings),
                    "LogP (o/w):": f"{log_p:.3f}",
                    "TPSA (Å²):": f"{tpsa:.2f}",
                    "H-Bond Donors:": str(num_h_donors),
                    "H-Bond Acceptors:": str(num_h_acceptors),
                }
        except Exception as e:
            main_layout.addWidget(QLabel(f"Error calculating properties: {e}"))
            return

        # --- 計算結果をUIに表示 ---
        row = 0
        for label_text, value_text in properties.items():
            label = QLabel(f"<b>{label_text}</b>")
            value = QLineEdit(value_text)
            value.setReadOnly(True)
            
            copy_btn = QPushButton("Copy")
            copy_btn.clicked.connect(lambda _, v=value: self.copy_to_clipboard(v.text()))

            grid_layout.addWidget(label, row, 0)
            grid_layout.addWidget(value, row, 1)
            grid_layout.addWidget(copy_btn, row, 2)
            row += 1
            
        main_layout.addLayout(grid_layout)
        
        # --- OKボタン ---
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        main_layout.addWidget(ok_button, 0, Qt.AlignmentFlag.AlignCenter)
        
        self.setLayout(main_layout)

    def copy_to_clipboard(self, text):
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        if self.parent() and hasattr(self.parent(), 'statusBar'):
            self.parent().statusBar().showMessage(f"Copied '{text}' to clipboard.", 2000)


class SettingsDialog(QDialog):
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("3D View Settings")
        self.setMinimumSize(500, 600)
        
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
            'wireframe_double_bond_radius_factor': 0.8,
            'wireframe_triple_bond_radius_factor': 0.75,
            'stick_double_bond_offset_factor': 1.5,
            'stick_triple_bond_offset_factor': 1.0,
            'stick_double_bond_radius_factor': 0.6,
            'stick_triple_bond_radius_factor': 0.4,
            # If True, attempts to be permissive when RDKit raises chemical/sanitization errors
            # during file import (useful for viewing malformed XYZ/MOL files). When enabled,
            # element symbol recognition will be coerced where possible and Chem.SanitizeMol
            # failures will be ignored so the 3D viewer can still display the structure.
            'skip_chemistry_checks': False,
            # 3D conversion/optimization defaults
            '3d_conversion_mode': 'fallback',
            'optimization_method': 'MMFF_RDKIT',
            'ball_stick_bond_color': '#7F7F7F',
            'cpk_colors': {},
        }
        
        # --- 選択された色を管理する専用のインスタンス変数 ---
        self.current_bg_color = None

        # --- UI要素の作成 ---
        layout = QVBoxLayout(self)
        
        # タブウィジェットを作成
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)

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
        
        self.tab_widget.addTab(scene_widget, "Scene")
    
    def create_other_tab(self):
        """other設定タブを作成"""
        other_widget = QWidget()
        form_layout = QFormLayout(other_widget)

        # 化学チェックスキップオプション（otherタブに移動）
        self.skip_chem_checks_checkbox = QCheckBox()
        self.skip_chem_checks_checkbox.setToolTip("When enabled, file import will try to ignore chemical/sanitization errors and allow viewing malformed files.")
        # Immediately persist change to settings when user toggles the checkbox
        try:
            self.skip_chem_checks_checkbox.stateChanged.connect(lambda s: self._on_skip_chem_checks_changed(s))
        except Exception:
            pass

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

        # NOTE: Multi-bond offset/thickness settings moved to per-model tabs to allow
        # independent configuration for Ball&Stick/CPK/Wireframe/Stick.
                
        # Add the checkbox to the form
        form_layout.addRow("Skip chemistry checks on import xyz file:", self.skip_chem_checks_checkbox)

    
        self.tab_widget.addTab(other_widget, "Other")
    
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

        # --- Ball & Stick bond color ---
        self.bs_bond_color_button = QPushButton()
        self.bs_bond_color_button.setFixedSize(36, 24)
        self.bs_bond_color_button.clicked.connect(self.pick_bs_bond_color)
        self.bs_bond_color_button.setToolTip("Choose the uniform bond color for Ball & Stick model (3D)")
        form_layout.addRow("Ball & Stick bond color:", self.bs_bond_color_button)

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
        
        # 原子半径
        self.stick_atom_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_atom_radius_slider.setRange(5, 50)  # 0.05 ~ 0.5
        self.stick_atom_radius_label = QLabel("0.15")
        self.stick_atom_radius_slider.valueChanged.connect(lambda v: self.stick_atom_radius_label.setText(f"{v/100:.2f}"))
        atom_radius_layout = QHBoxLayout()
        atom_radius_layout.addWidget(self.stick_atom_radius_slider)
        atom_radius_layout.addWidget(self.stick_atom_radius_label)
        form_layout.addRow("Atom Radius:", atom_radius_layout)
        
        # ボンド半径
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
        # Note: tab labels must match those added to the QTabWidget ("Scene", "Ball & Stick",
        # "CPK (Space-filling)", "Wireframe", "Stick", "Other"). Use the per-model
        # multi-bond keys present in self.default_settings.
        tab_settings = {
            "Scene": {
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
            },
            "Ball & Stick": {
                'ball_stick_atom_scale': self.default_settings['ball_stick_atom_scale'],
                'ball_stick_bond_radius': self.default_settings['ball_stick_bond_radius'],
                'ball_stick_resolution': self.default_settings['ball_stick_resolution'],
                'ball_stick_double_bond_offset_factor': self.default_settings.get('ball_stick_double_bond_offset_factor', 2.0),
                'ball_stick_triple_bond_offset_factor': self.default_settings.get('ball_stick_triple_bond_offset_factor', 2.0),
                'ball_stick_double_bond_radius_factor': self.default_settings.get('ball_stick_double_bond_radius_factor', 0.8),
                'ball_stick_triple_bond_radius_factor': self.default_settings.get('ball_stick_triple_bond_radius_factor', 0.75)
            },
            "CPK (Space-filling)": {
                'cpk_atom_scale': self.default_settings['cpk_atom_scale'],
                'cpk_resolution': self.default_settings['cpk_resolution'],
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
                'stick_atom_radius': self.default_settings['stick_atom_radius'],
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
                        # 2D bond color is fixed and not part of the settings; do not modify it here
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
            'stick_atom_radius': self.stick_atom_radius_slider.value() / 100.0,
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
        stick_atom_radius = int(settings_dict.get('stick_atom_radius', self.default_settings['stick_atom_radius']) * 100)
        self.stick_atom_radius_slider.setValue(stick_atom_radius)
        self.stick_atom_radius_label.setText(f"{stick_atom_radius/100:.2f}")
        
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
            'stick_atom_radius': self.stick_atom_radius_slider.value() / 100.0,
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
            'skip_chemistry_checks': self.skip_chem_checks_checkbox.isChecked(),
            # Ball & Stick bond color (3D grey/uniform color)
            'ball_stick_bond_color': getattr(self, 'bs_bond_color', self.default_settings.get('ball_stick_bond_color', '#7F7F7F')),
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

    def accept(self):
        """ダイアログの設定を適用してから閉じる"""
        # apply_settingsを呼び出して設定を適用
        self.apply_settings()
        super().accept()


class CustomQtInteractor(QtInteractor):
    def __init__(self, parent=None, main_window=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.main_window = main_window
        # Track recent clicks so we can detect and swallow triple-clicks
        # Triple-clicks are not a distinct Qt event on all platforms, so we
        # implement a small timing-based counter here and accept the event
        # when 3 rapid clicks are detected to prevent them from reaching
        # the VTK interactor and causing unexpected behaviour in the 3D view.
        self._last_click_time = 0.0
        self._click_count = 0

    def wheelEvent(self, event):
        """
        マウスホイールイベントをオーバーライドする。
        """
        # 最初に親クラスのイベントを呼び、通常のズーム処理を実行させる
        super().wheelEvent(event)
        
        # ズーム処理の完了後、2Dビューにフォーカスを強制的に戻す
        if self.main_window and hasattr(self.main_window, 'view_2d'):
            self.main_window.view_2d.setFocus()

    def mouseReleaseEvent(self, event):
        """
        Qtのマウスリリースイベントをオーバーライドし、
        3Dビューでの全ての操作完了後に2Dビューへフォーカスを戻す。
        """
        super().mouseReleaseEvent(event) # 親クラスのイベントを先に処理
        if self.main_window and hasattr(self.main_window, 'view_2d'):
            self.main_window.view_2d.setFocus()

    def mouseDoubleClickEvent(self, event):
        """Ignore mouse double-clicks on the 3D widget to avoid accidental actions.

        Swallow the double-click event so it doesn't trigger selection, editing,
        or camera jumps. We intentionally do not call the superclass handler.
        """
        try:
            # Accept the event to mark it handled and prevent further processing.
            event.accept()
        except Exception:
            # If event doesn't support accept for some reason, just return.
            return

# --- 3Dインタラクションを管理する専用クラス ---
class CustomInteractorStyle(vtkInteractorStyleTrackballCamera):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        # カスタム状態を管理するフラグを一つに絞ります
        self._is_dragging_atom = False
        # undoスタックのためのフラグ
        self.is_dragging = False
        # 回転操作を検出するためのフラグ
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None

        self.AddObserver("LeftButtonPressEvent", self.on_left_button_down)
        self.AddObserver("RightButtonPressEvent", self.on_right_button_down)
        self.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.AddObserver("LeftButtonReleaseEvent", self.on_left_button_up)
        self.AddObserver("RightButtonReleaseEvent", self.on_right_button_up)

    def on_left_button_down(self, obj, event):
        """
        クリック時の処理を振り分けます。
        原子を掴めた場合のみカスタム動作に入り、それ以外は親クラス（カメラ回転）に任せます。
        """
        mw = self.main_window
        
        # 前回のドラッグ状態をクリア（トリプルクリック/ダブルクリック対策）
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None
        
        # Move Groupダイアログが開いている場合の処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and move_group_dialog.group_atoms:
            # グループが選択されている場合、グループドラッグ処理
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)
                
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5
                        
                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)
            
            # グループ内の原子がクリックされた場合
            if clicked_atom_idx is not None:
                if clicked_atom_idx in move_group_dialog.group_atoms:
                    # 既存グループ内の原子 - ドラッグ準備
                    move_group_dialog._is_dragging_group_vtk = True
                    move_group_dialog._drag_atom_idx = clicked_atom_idx
                    move_group_dialog._drag_start_pos = click_pos
                    move_group_dialog._mouse_moved = False
                    # 初期位置を保存
                    move_group_dialog._initial_positions = {}
                    conf = mw.current_mol.GetConformer()
                    for atom_idx in move_group_dialog.group_atoms:
                        pos = conf.GetAtomPosition(atom_idx)
                        move_group_dialog._initial_positions[atom_idx] = np.array([pos.x, pos.y, pos.z])
                    mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                    return  # カメラ回転を無効化
                else:
                    # グループ外の原子をクリック - BFS/DFSで連結成分を探索
                    visited = set()
                    queue = [clicked_atom_idx]
                    visited.add(clicked_atom_idx)
                    
                    while queue:
                        current_idx = queue.pop(0)
                        for bond_idx in range(mw.current_mol.GetNumBonds()):
                            bond = mw.current_mol.GetBondWithIdx(bond_idx)
                            begin_idx = bond.GetBeginAtomIdx()
                            end_idx = bond.GetEndAtomIdx()
                            
                            if begin_idx == current_idx and end_idx not in visited:
                                visited.add(end_idx)
                                queue.append(end_idx)
                            elif end_idx == current_idx and begin_idx not in visited:
                                visited.add(begin_idx)
                                queue.append(begin_idx)
                    
                    # Ctrlキーが押されている場合のみ複数グループ選択
                    is_ctrl_pressed = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier)
                    
                    if is_ctrl_pressed:
                        # Ctrl + クリック: 追加または解除
                        if visited.issubset(move_group_dialog.group_atoms):
                            # すでに選択されている - 解除
                            move_group_dialog.group_atoms -= visited
                        else:
                            # 新しいグループを追加
                            move_group_dialog.group_atoms |= visited
                    else:
                        # 通常のクリック: 既存の選択を置き換え
                        move_group_dialog.group_atoms = visited.copy()
                    
                    move_group_dialog.selected_atoms.add(clicked_atom_idx)
                    move_group_dialog.show_atom_labels()
                    move_group_dialog.update_display()
                    return
            else:
                # 原子以外をクリック - 全選択を解除
                move_group_dialog.group_atoms.clear()
                move_group_dialog.selected_atoms.clear()
                move_group_dialog.clear_atom_labels()
                move_group_dialog.update_display()
                # カメラ回転を許可
                super(CustomInteractorStyle, self).OnLeftButtonDown()
                return
        
        is_temp_mode = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.AltModifier)
        is_edit_active = mw.is_3d_edit_mode or is_temp_mode
        
        # Ctrl+クリックで原子選択（3D編集用）
        is_ctrl_click = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier)

        # 測定モードが有効な場合の処理
        if mw.measurement_mode and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            self._mouse_press_pos = click_pos  # マウスプレス位置を記録
            self._mouse_moved_during_drag = False  # 移動フラグをリセット
            
            picker = mw.plotter.picker
            
            # 通常のピック処理を実行
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            # 原子がクリックされた場合のみ特別処理
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)

                # 範囲チェックを追加
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # クリック閾値チェック
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            mw.handle_measurement_atom_selection(int(closest_atom_idx))
                            return  # 原子選択処理完了、カメラ回転は無効
            
            # 測定モードで原子以外をクリックした場合は計測選択をクリア
            # ただし、これは通常のカメラ回転も許可する
            self._is_dragging_atom = False
            super().OnLeftButtonDown()
            return
        
        # Ctrl+クリックの原子選択機能は無効化（Move Group機能で代替）
        # if is_ctrl_click and mw.current_mol:
        #     ... (無効化)

        # 3D分子(mw.current_mol)が存在する場合のみ、原子の選択処理を実行
        if is_edit_active and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)

                # 範囲チェックを追加
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # RDKitのMolオブジェクトから原子を安全に取得
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            # 原子を掴むことに成功した場合
                            self._is_dragging_atom = True
                        self.is_dragging = False 
                        mw.dragged_atom_info = {'id': int(closest_atom_idx)}
                        mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                        return  # 親クラスのカメラ回転を呼ばない

        self._is_dragging_atom = False
        super().OnLeftButtonDown()

    def on_right_button_down(self, obj, event):
        """
        右クリック時の処理。Move Groupダイアログが開いている場合はグループ回転を開始。
        """
        mw = self.main_window
        
        # Move Groupダイアログが開いているか確認
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and move_group_dialog.group_atoms:
            # グループが選択されている場合、回転ドラッグを開始
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)
                
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5
                        
                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)
            
            # グループ内の原子がクリックされた場合、回転ドラッグを開始
            if clicked_atom_idx is not None and clicked_atom_idx in move_group_dialog.group_atoms:
                move_group_dialog._is_rotating_group_vtk = True
                move_group_dialog._rotation_start_pos = click_pos
                move_group_dialog._rotation_mouse_moved = False
                move_group_dialog._rotation_atom_idx = clicked_atom_idx  # 掴んだ原子を記録
                
                # 初期位置と重心を保存
                move_group_dialog._initial_positions = {}
                conf = mw.current_mol.GetConformer()
                centroid = np.zeros(3)
                for atom_idx in move_group_dialog.group_atoms:
                    pos = conf.GetAtomPosition(atom_idx)
                    pos_array = np.array([pos.x, pos.y, pos.z])
                    move_group_dialog._initial_positions[atom_idx] = pos_array
                    centroid += pos_array
                centroid /= len(move_group_dialog.group_atoms)
                move_group_dialog._group_centroid = centroid
                
                mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                return  # カメラ回転を無効化
        
        # 通常の右クリック処理
        super().OnRightButtonDown()

    def on_mouse_move(self, obj, event):
        """
        マウス移動時の処理。原子ドラッグ中か、それ以外（カメラ回転＋ホバー）かをハンドリングします。
        """
        mw = self.main_window
        
        # Move Groupダイアログのドラッグ処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and getattr(move_group_dialog, '_is_dragging_group_vtk', False):
            # グループをドラッグ中 - 移動距離を記録するのみ
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()
            
            dx = current_pos[0] - move_group_dialog._drag_start_pos[0]
            dy = current_pos[1] - move_group_dialog._drag_start_pos[1]
            
            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._mouse_moved = True
            
            return  # カメラ回転を無効化
        
        # グループ回転中の処理
        if move_group_dialog and getattr(move_group_dialog, '_is_rotating_group_vtk', False):
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()
            
            dx = current_pos[0] - move_group_dialog._rotation_start_pos[0]
            dy = current_pos[1] - move_group_dialog._rotation_start_pos[1]
            
            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._rotation_mouse_moved = True
            
            return  # カメラ回転を無効化
        
        interactor = self.GetInteractor()

        # マウス移動があったことを記録
        if self._mouse_press_pos is not None:
            current_pos = interactor.GetEventPosition()
            if abs(current_pos[0] - self._mouse_press_pos[0]) > 3 or abs(current_pos[1] - self._mouse_press_pos[1]) > 3:
                self._mouse_moved_during_drag = True

        if self._is_dragging_atom and mw.dragged_atom_info is not None:
            # カスタムの原子ドラッグ処理
            self.is_dragging = True
            atom_id = mw.dragged_atom_info['id']
            # We intentionally do NOT update visible coordinates or the
            # authoritative atom position during mouse-move while dragging.
            # The UX requirement here is that atoms need not visibly move
            # while the mouse is being dragged. Compute and apply the final
            # world-coordinate only once on mouse release (on_left_button_up).
            # Keep minimal state: mark that a drag occurred (is_dragging)
            # and allow the release handler to compute the final position.
            # This avoids duplicate updates and simplifies event ordering.
        else:
            # カメラ回転処理を親クラスに任せます
            super().OnMouseMove()

            # その後、カーソルの表示を更新します
            is_edit_active = mw.is_3d_edit_mode or interactor.GetAltKey()
            if is_edit_active:
                # 編集がアクティブな場合のみ、原子のホバーチェックを行う
                atom_under_cursor = False
                click_pos = interactor.GetEventPosition()
                picker = mw.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
                if picker.GetActor() is mw.atom_actor:
                    atom_under_cursor = True

                if atom_under_cursor:
                    mw.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                else:
                    mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            else:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)

    def on_left_button_up(self, obj, event):
        """
        クリック終了時の処理。状態をリセットします。
        """
        mw = self.main_window
        
        # Move Groupダイアログのドラッグ終了処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        # ダブルクリック/トリプルクリックで状態が混乱するのを防ぐ（Move Group用）
        if move_group_dialog:
            if getattr(move_group_dialog, '_is_dragging_group_vtk', False) and not getattr(move_group_dialog, '_mouse_moved', False):
                # ドラッグしていない状態で複数クリックされた場合は状態をリセット
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, '_initial_positions'):
                    delattr(move_group_dialog, '_initial_positions')
        
        if move_group_dialog and getattr(move_group_dialog, '_is_dragging_group_vtk', False):
            if getattr(move_group_dialog, '_mouse_moved', False):
                # ドラッグが実行された - リリース時に座標を更新
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()
                    
                    # ドラッグ原子の初期位置
                    drag_atom_initial_pos = move_group_dialog._initial_positions[move_group_dialog._drag_atom_idx]
                    
                    # スクリーン座標からワールド座標への変換
                    renderer.SetWorldPoint(drag_atom_initial_pos[0], drag_atom_initial_pos[1], drag_atom_initial_pos[2], 1.0)
                    renderer.WorldToDisplay()
                    display_coords = renderer.GetDisplayPoint()
                    
                    new_display_pos = (current_pos[0], current_pos[1], display_coords[2])
                    renderer.SetDisplayPoint(new_display_pos[0], new_display_pos[1], new_display_pos[2])
                    renderer.DisplayToWorld()
                    new_world_coords = renderer.GetWorldPoint()
                    
                    # 移動ベクトル
                    translation_vector = np.array([
                        new_world_coords[0] - drag_atom_initial_pos[0],
                        new_world_coords[1] - drag_atom_initial_pos[1],
                        new_world_coords[2] - drag_atom_initial_pos[2]
                    ])
                    
                    # グループ全体を移動
                    for atom_idx in move_group_dialog.group_atoms:
                        initial_pos = move_group_dialog._initial_positions[atom_idx]
                        new_pos = initial_pos + translation_vector
                        conf.SetAtomPosition(atom_idx, new_pos.tolist())
                        mw.atom_positions_3d[atom_idx] = new_pos
                    
                    # 3D表示を更新
                    mw.draw_molecule_3d(mw.current_mol)
                    mw.update_chiral_labels()
                    move_group_dialog.show_atom_labels()
                    mw.push_undo_state()
                except Exception as e:
                    print(f"Error finalizing group drag: {e}")
            else:
                # ドラッグがなかった = クリックのみ → トグル処理
                if hasattr(move_group_dialog, '_drag_atom_idx'):
                    clicked_atom = move_group_dialog._drag_atom_idx
                    try:
                        move_group_dialog.on_atom_picked(clicked_atom)
                    except Exception as e:
                        print(f"Error in toggle: {e}")
            
            # 状態をリセット（完全なクリーンアップ）
            move_group_dialog._is_dragging_group_vtk = False
            move_group_dialog._drag_start_pos = None
            move_group_dialog._mouse_moved = False
            if hasattr(move_group_dialog, '_initial_positions'):
                delattr(move_group_dialog, '_initial_positions')
            if hasattr(move_group_dialog, '_drag_atom_idx'):
                delattr(move_group_dialog, '_drag_atom_idx')
            
            # CustomInteractorStyleの状態もクリア
            self._is_dragging_atom = False
            self.is_dragging = False
            self._mouse_moved_during_drag = False
            self._mouse_press_pos = None
            
            try:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except Exception:
                pass
            return

        # 計測モードで、マウスが動いていない場合（つまりクリック）の処理
        if mw.measurement_mode and not self._mouse_moved_during_drag and self._mouse_press_pos is not None:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            # 原子がクリックされていない場合は測定選択をクリア
            if picker.GetActor() is not mw.atom_actor:
                mw.clear_measurement_selection()

        if self._is_dragging_atom:
            # カスタムドラッグの後始末
            if self.is_dragging:
                if mw.current_mol and mw.current_mol.GetNumConformers() > 0:
                    try:
                        # Before applying conformer updates, compute the final
                        # world coordinates for the dragged atom based on the
                        # release pointer position. During the drag we did not
                        # update mw.atom_positions_3d (to keep the visuals
                        # static). Now compute the final position for the
                        # dragged atom and store it into mw.atom_positions_3d
                        # so the conformer update loop below will pick it up.
                        atom_id = None
                        try:
                            atom_id = mw.dragged_atom_info.get('id') if mw.dragged_atom_info else None
                        except Exception:
                            atom_id = None

                        if atom_id is not None:
                            try:
                                interactor = self.GetInteractor()
                                renderer = mw.plotter.renderer
                                current_display_pos = interactor.GetEventPosition()
                                conf = mw.current_mol.GetConformer()
                                # Use the atom's current 3D position to obtain a
                                # display-space depth (z) value, then replace the
                                # x/y with the pointer position to project back to
                                # world coordinates at that depth.
                                pos_3d = conf.GetAtomPosition(atom_id)
                                renderer.SetWorldPoint(pos_3d.x, pos_3d.y, pos_3d.z, 1.0)
                                renderer.WorldToDisplay()
                                display_coords = renderer.GetDisplayPoint()
                                new_display_pos = (current_display_pos[0], current_display_pos[1], display_coords[2])
                                renderer.SetDisplayPoint(new_display_pos[0], new_display_pos[1], new_display_pos[2])
                                renderer.DisplayToWorld()
                                new_world_coords_tuple = renderer.GetWorldPoint()
                                new_world_coords = list(new_world_coords_tuple)[:3]
                                # Ensure the container supports assignment
                                try:
                                    mw.atom_positions_3d[atom_id] = new_world_coords
                                except Exception:
                                    # If atom_positions_3d is immutable or shaped
                                    # differently, attempt a safe conversion.
                                    try:
                                        ap = list(mw.atom_positions_3d)
                                        ap[atom_id] = new_world_coords
                                        mw.atom_positions_3d = ap
                                    except Exception:
                                        pass
                            except Exception:
                                # If final-position computation fails, continue
                                # and apply whatever state is available.
                                pass

                        # Apply the (now updated) positions to the RDKit conformer
                        # exactly once. This ensures the conformer is
                        # authoritative and avoids double-moves.
                        conf = mw.current_mol.GetConformer()
                        for i in range(mw.current_mol.GetNumAtoms()):
                            try:
                                pos = mw.atom_positions_3d[i]
                                conf.SetAtomPosition(i, pos.tolist())
                            except Exception:
                                # Skip individual failures but continue applying
                                # other atom positions.
                                pass
                    except Exception:
                        # If applying positions fails, continue to redraw from
                        # whatever authoritative state is available.
                        pass

                    # Redraw once and push undo state
                    try:
                        mw.draw_molecule_3d(mw.current_mol)
                    except Exception:
                        pass
                    mw.push_undo_state()
            mw.dragged_atom_info = None
            # Refresh overlays and labels that depend on atom_positions_3d. Do
            # not overwrite mw.atom_positions_3d here — it already reflects the
            # positions the user dragged to. Only update dependent displays.
            try:
                mw.update_3d_selection_display()
            except Exception:
                pass
            try:
                mw.update_measurement_labels_display()
            except Exception:
                pass
            try:
                mw.update_2d_measurement_labels()
            except Exception:
                pass
            try:
                mw.show_all_atom_info()
            except Exception:
                pass
            except Exception:
                # Do not allow a failure here to interrupt release flow
                pass
        else:
            # カメラ回転の後始末を親クラスに任せます
            super().OnLeftButtonUp()

        # 状態をリセット（完全なクリーンアップ）
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_press_pos = None
        self._mouse_moved_during_drag = False
        
        # Move Group関連の状態もクリア
        try:
            move_group_dialog = None
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
            
            if move_group_dialog:
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, '_initial_positions'):
                    delattr(move_group_dialog, '_initial_positions')
                if hasattr(move_group_dialog, '_drag_atom_idx'):
                    delattr(move_group_dialog, '_drag_atom_idx')
        except Exception:
            pass
        
        # ピックリセットは測定モードで実際に問題が発生した場合のみ行う
        # （通常のドラッグ回転では行わない）
        
        # ボタンを離した後のカーソル表示を最新の状態に更新
        self.on_mouse_move(obj, event)

        # 2Dビューにフォーカスを戻し、ショートカットキーなどが使えるようにする
        if mw and mw.view_2d:
            mw.view_2d.setFocus()

    def on_right_button_up(self, obj, event):
        """
        右クリック終了時の処理。グループ回転を確定。
        """
        mw = self.main_window
        
        # Move Groupダイアログの回転終了処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and getattr(move_group_dialog, '_is_rotating_group_vtk', False):
            # 回転モードで右クリックリリース - 選択を保持
            if getattr(move_group_dialog, '_rotation_mouse_moved', False):
                # 回転が実行された - リリース時に回転を適用
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()
                    centroid = move_group_dialog._group_centroid
                    
                    # 掴んだ原子の初期位置
                    if not hasattr(move_group_dialog, '_rotation_atom_idx'):
                        # 最初に掴んだ原子のインデックスを保存
                        move_group_dialog._rotation_atom_idx = next(iter(move_group_dialog.group_atoms))
                    
                    grabbed_atom_idx = move_group_dialog._rotation_atom_idx
                    grabbed_initial_pos = move_group_dialog._initial_positions[grabbed_atom_idx]
                    
                    # 開始位置のスクリーン座標を取得
                    renderer.SetWorldPoint(grabbed_initial_pos[0], grabbed_initial_pos[1], grabbed_initial_pos[2], 1.0)
                    renderer.WorldToDisplay()
                    start_display = renderer.GetDisplayPoint()
                    
                    # 現在のマウス位置をワールド座標に変換（同じ深度で）
                    renderer.SetDisplayPoint(current_pos[0], current_pos[1], start_display[2])
                    renderer.DisplayToWorld()
                    target_world = renderer.GetWorldPoint()
                    target_pos = np.array([target_world[0], target_world[1], target_world[2]])
                    
                    # 重心から見た、掴んだ原子の初期ベクトルと目標ベクトル
                    v1 = grabbed_initial_pos - centroid
                    v2 = target_pos - centroid
                    
                    # ベクトルを正規化
                    v1_norm = np.linalg.norm(v1)
                    v2_norm = np.linalg.norm(v2)
                    
                    if v1_norm > 1e-6 and v2_norm > 1e-6:
                        v1_normalized = v1 / v1_norm
                        v2_normalized = v2 / v2_norm
                        
                        # 回転軸（外積）
                        rotation_axis = np.cross(v1_normalized, v2_normalized)
                        axis_norm = np.linalg.norm(rotation_axis)
                        
                        if axis_norm > 1e-6:
                            rotation_axis = rotation_axis / axis_norm
                            
                            # 回転角（内積）
                            cos_angle = np.clip(np.dot(v1_normalized, v2_normalized), -1.0, 1.0)
                            angle = np.arccos(cos_angle)
                            
                            # Rodriguesの回転公式で回転行列を作成
                            K = np.array([
                                [0, -rotation_axis[2], rotation_axis[1]],
                                [rotation_axis[2], 0, -rotation_axis[0]],
                                [-rotation_axis[1], rotation_axis[0], 0]
                            ])
                            
                            rot_matrix = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
                            
                            # グループ全体を重心周りに回転
                            for atom_idx in move_group_dialog.group_atoms:
                                initial_pos = move_group_dialog._initial_positions[atom_idx]
                                # 重心からの相対座標
                                relative_pos = initial_pos - centroid
                                # 回転を適用
                                rotated_pos = rot_matrix @ relative_pos
                                # 絶対座標に戻す
                                new_pos = rotated_pos + centroid
                                
                                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                                mw.atom_positions_3d[atom_idx] = new_pos
                            
                            # 3D表示を更新
                            mw.draw_molecule_3d(mw.current_mol)
                            mw.update_chiral_labels()
                            move_group_dialog.show_atom_labels()
                            mw.push_undo_state()
                except Exception as e:
                    print(f"Error finalizing group rotation: {e}")
            
            # 状態をリセット
            move_group_dialog._is_rotating_group_vtk = False
            move_group_dialog._rotation_start_pos = None
            move_group_dialog._rotation_mouse_moved = False
            if hasattr(move_group_dialog, '_initial_positions'):
                delattr(move_group_dialog, '_initial_positions')
            if hasattr(move_group_dialog, '_group_centroid'):
                delattr(move_group_dialog, '_group_centroid')
            if hasattr(move_group_dialog, '_rotation_atom_idx'):
                delattr(move_group_dialog, '_rotation_atom_idx')
            
            try:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except Exception:
                pass
            return
        
        # 通常の右クリックリリース処理
        super().OnRightButtonUp()

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
        if self.plotter and self.plotter.renderer:
            bg_color = self.settings.get('background_color', '#919191')
            self.plotter.set_background(bg_color)
            self.apply_3d_settings()
        # Apply any CPK overrides here so both 2D and 3D use the configured colors
        try:
            self.update_cpk_colors_from_settings()
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


    def apply_3d_settings(self):
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


# --- 3D Editing Dialog Classes ---

class BondLengthDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 2:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Bond Length")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click two atoms in the 3D view to select a bond, then specify the new length.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current distance display
        self.distance_label = QLabel("")
        layout.addWidget(self.distance_label)
        
        # New distance input
        distance_layout = QHBoxLayout()
        distance_layout.addWidget(QLabel("New distance (Å):"))
        self.distance_input = QLineEdit()
        self.distance_input.setPlaceholderText("1.54")
        distance_layout.addWidget(self.distance_input)
        layout.addLayout(distance_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Movement Options:"))
        
        self.atom1_fix_group_radio = QRadioButton("Atom 1: Fixed, Atom 2: Move connected group")
        self.atom1_fix_group_radio.setChecked(True)
        group_layout.addWidget(self.atom1_fix_group_radio)

        self.atom1_fix_radio = QRadioButton("Atom 1: Fixed, Atom 2: Move atom only")
        group_layout.addWidget(self.atom1_fix_radio)
        
        self.both_groups_radio = QRadioButton("Both groups: Move towards center equally")
        group_layout.addWidget(self.both_groups_radio)
        
        layout.addWidget(group_box)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
        
        # 事前選択された原子がある場合は初期表示を更新
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_changes()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
    
    def clear_selection(self):
        """選択をクリア"""
        self.atom1_idx = None
        self.atom2_idx = None
        self.clear_selection_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx]
        labels = ["1st", "2nd"]
        colors = ["yellow", "yellow"]
        
        for i, atom_idx in enumerate(selected_atoms):
            if atom_idx is not None:
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"{labels[i]}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color=colors[i],
                    always_visible=True
                )
                self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def clear_selection_labels(self):
        """選択ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def add_selection_label(self, atom_idx, label_text):
        """選択された原子にラベルを追加"""
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        # 原子の位置を取得
        pos = self.main_window.atom_positions_3d[atom_idx]
        
        # ラベルを追加
        label_actor = self.main_window.plotter.add_point_labels(
            [pos], [label_text], 
            point_size=20, 
            font_size=12,
            text_color='yellow',
            always_visible=True
        )
        self.selection_labels.append(label_actor)
    
    def update_display(self):
        """表示を更新"""
        # 既存のラベルをクリア
        self.clear_selection_labels()
        
        if self.atom1_idx is None:
            self.selection_label.setText("No atoms selected")
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear distance input when no selection
            try:
                self.distance_input.clear()
            except Exception:
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            # Clear distance input while selection is incomplete
            try:
                self.distance_input.clear()
            except Exception:
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(f"Bond: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx})")
            
            # Calculate current distance
            conf = self.mol.GetConformer()
            pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
            pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
            current_distance = np.linalg.norm(pos2 - pos1)
            self.distance_label.setText(f"Current distance: {current_distance:.3f} Å")
            self.apply_button.setEnabled(True)
            # Update the distance input box to show current distance
            try:
                self.distance_input.setText(f"{current_distance:.3f}")
            except Exception:
                pass
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")
    
    def apply_changes(self):
        """変更を適用"""
        if self.atom1_idx is None or self.atom2_idx is None:
            return
        
        try:
            new_distance = float(self.distance_input.text())
            if new_distance <= 0:
                QMessageBox.warning(self, "Invalid Input", "Distance must be positive.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Undo状態を保存
        self.main_window.push_undo_state()
        
        # Apply the bond length change
        self.adjust_bond_length(new_distance)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()
    
    def adjust_bond_length(self, new_distance):
        """結合長を調整"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        
        # Direction vector from atom1 to atom2
        direction = pos2 - pos1
        current_distance = np.linalg.norm(direction)
        
        if current_distance == 0:
            return
        
        direction = direction / current_distance
        
        if self.both_groups_radio.isChecked():
            # Both groups move towards center equally
            bond_center = (pos1 + pos2) / 2
            half_distance = new_distance / 2
            
            # New positions for both atoms
            new_pos1 = bond_center - direction * half_distance
            new_pos2 = bond_center + direction * half_distance
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom1_idx, exclude=self.atom2_idx)
            group2_atoms = self.get_connected_group(self.atom2_idx, exclude=self.atom1_idx)
            
            # Calculate displacements
            displacement1 = new_pos1 - pos1
            displacement2 = new_pos2 - pos2
            
            # Move group 1
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement1
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Move group 2
            for atom_idx in group2_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement2
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.atom1_fix_radio.isChecked():
            # Move only the second atom
            new_pos2 = pos1 + direction * new_distance
            conf.SetAtomPosition(self.atom2_idx, new_pos2.tolist())
            self.main_window.atom_positions_3d[self.atom2_idx] = new_pos2
        else:
            # Move the connected group (default behavior)
            new_pos2 = pos1 + direction * new_distance
            atoms_to_move = self.get_connected_group(self.atom2_idx, exclude=self.atom1_idx)
            displacement = new_pos2 - pos2
            
            for atom_idx in atoms_to_move:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)
    
    def get_connected_group(self, start_atom, exclude=None):
        """指定された原子から連結されているグループを取得"""
        visited = set()
        to_visit = [start_atom]
        
        while to_visit:
            current = to_visit.pop()
            if current in visited or current == exclude:
                continue
            
            visited.add(current)
            
            # Get neighboring atoms
            atom = self.mol.GetAtomWithIdx(current)
            for bond in atom.GetBonds():
                other_idx = bond.GetOtherAtomIdx(current)
                if other_idx not in visited and other_idx != exclude:
                    to_visit.append(other_idx)
        
        return visited


class AngleDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 3:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # vertex
            self.atom3_idx = preselected_atoms[2]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Angle")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click three atoms in order: first-vertex-third. The angle around the vertex atom will be adjusted.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current angle display
        self.angle_label = QLabel("")
        layout.addWidget(self.angle_label)
        
        # New angle input
        angle_layout = QHBoxLayout()
        angle_layout.addWidget(QLabel("New angle (degrees):"))
        self.angle_input = QLineEdit()
        self.angle_input.setPlaceholderText("109.5")
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Rotation Options:"))
        
        self.rotate_group_radio = QRadioButton("Atom 1,2: Fixed, Atom 3: Rotate connected group")
        self.rotate_group_radio.setChecked(True)
        group_layout.addWidget(self.rotate_group_radio)

        self.rotate_atom_radio = QRadioButton("Atom 1,2: Fixed, Atom 3: Rotate atom only")
        group_layout.addWidget(self.rotate_atom_radio)
        
        self.both_groups_radio = QRadioButton("Vertex fixed: Both arms rotate equally")
        group_layout.addWidget(self.both_groups_radio)
    
        
        layout.addWidget(group_box)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker for AngleDialog
        self.picker_connection = None
        self.enable_picking()
        
        # 事前選択された原子がある場合は初期表示を更新
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_changes()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
    
    def clear_selection(self):
        """選択をクリア"""
        self.atom1_idx = None
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None
        self.clear_selection_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx]
        labels = ["1st", "2nd (vertex)", "3rd"]
        colors = ["yellow", "yellow", "yellow"]  # 全て黄色に統一
        
        for i, atom_idx in enumerate(selected_atoms):
            if atom_idx is not None:
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"{labels[i]}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color=colors[i],
                    always_visible=True
                )
                self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def clear_selection_labels(self):
        """選択ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def add_selection_label(self, atom_idx, label_text):
        """選択された原子にラベルを追加"""
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        # 原子の位置を取得
        pos = self.main_window.atom_positions_3d[atom_idx]
        
        # ラベルを追加
        label_actor = self.main_window.plotter.add_point_labels(
            [pos], [label_text], 
            point_size=20, 
            font_size=12,
            text_color='yellow',
            always_visible=True
        )
        self.selection_labels.append(label_actor)
    
    def update_display(self):
        """表示を更新"""
        # 既存のラベルをクリア
        self.clear_selection_labels()
        
        if self.atom1_idx is None:
            self.selection_label.setText("No atoms selected")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear angle input when no selection
            try:
                self.angle_input.clear()
            except Exception:
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.clear()
            except Exception:
                pass
        elif self.atom3_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - ?")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.clear()
            except Exception:
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            self.selection_label.setText(f"Angle: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - {symbol3}({self.atom3_idx})")
            
            # Calculate current angle
            current_angle = self.calculate_angle()
            self.angle_label.setText(f"Current angle: {current_angle:.2f}°")
            self.apply_button.setEnabled(True)
            # Update angle input box with current angle
            try:
                self.angle_input.setText(f"{current_angle:.2f}")
            except Exception:
                pass
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            self.add_selection_label(self.atom3_idx, "3")
    
    def calculate_angle(self):
        """現在の角度を計算"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))  # vertex
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)
    
    def apply_changes(self):
        """変更を適用"""
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return
        
        try:
            new_angle = float(self.angle_input.text())
            if new_angle < 0 or new_angle >= 360:
                QMessageBox.warning(self, "Invalid Input", "Angle must be between 0 and 360 degrees.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Undo状態を保存
        self.main_window.push_undo_state()
        
        # Apply the angle change
        self.adjust_angle(new_angle)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()
    
    def adjust_angle(self, new_angle_deg):
        """角度を調整（均等回転オプション付き）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))  # vertex
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        # Current angle
        current_angle_rad = np.arccos(np.clip(
            np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)), -1.0, 1.0))
        
        # Target angle
        target_angle_rad = np.radians(new_angle_deg)
        
        # Rotation axis (perpendicular to the plane containing vec1 and vec2)
        rotation_axis = np.cross(vec1, vec2)
        rotation_axis_norm = np.linalg.norm(rotation_axis)
        
        if rotation_axis_norm == 0:
            # Vectors are parallel, cannot rotate
            return
        
        rotation_axis = rotation_axis / rotation_axis_norm
        
        # Total rotation angle needed
        total_rotation_angle = target_angle_rad - current_angle_rad
        
        # Rodrigues' rotation formula
        def rotate_vector(v, axis, angle):
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)
        
        if self.both_groups_radio.isChecked():
            # Both arms rotate equally (half angle each in opposite directions)
            half_rotation = total_rotation_angle / 2
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom1_idx, exclude=self.atom2_idx)
            group3_atoms = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate group 1 by -half_rotation
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                relative_pos = current_pos - pos2
                rotated_pos = rotate_vector(relative_pos, rotation_axis, -half_rotation)
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Rotate group 3 by +half_rotation
            for atom_idx in group3_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                relative_pos = current_pos - pos2
                rotated_pos = rotate_vector(relative_pos, rotation_axis, half_rotation)
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.rotate_atom_radio.isChecked():
            # Move only the third atom
            new_vec2 = rotate_vector(vec2, rotation_axis, total_rotation_angle)
            new_pos3 = pos2 + new_vec2
            conf.SetAtomPosition(self.atom3_idx, new_pos3.tolist())
            self.main_window.atom_positions_3d[self.atom3_idx] = new_pos3
        else:
            # Rotate the connected group around atom2 (vertex) - default behavior
            atoms_to_move = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            for atom_idx in atoms_to_move:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                # Transform to coordinate system centered at atom2
                relative_pos = current_pos - pos2
                # Rotate around the rotation axis
                rotated_pos = rotate_vector(relative_pos, rotation_axis, total_rotation_angle)
                # Transform back to world coordinates
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)
    
    def get_connected_group(self, start_atom, exclude=None):
        """指定された原子から連結されているグループを取得"""
        visited = set()
        to_visit = [start_atom]
        
        while to_visit:
            current = to_visit.pop()
            if current in visited or current == exclude:
                continue
            
            visited.add(current)
            
            # Get neighboring atoms
            atom = self.mol.GetAtomWithIdx(current)
            for bond in atom.GetBonds():
                other_idx = bond.GetOtherAtomIdx(current)
                if other_idx not in visited and other_idx != exclude:
                    to_visit.append(other_idx)
        
        return visited


class DihedralDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 4:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # central bond start
            self.atom3_idx = preselected_atoms[2]  # central bond end
            self.atom4_idx = preselected_atoms[3]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Dihedral Angle")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click four atoms in order to define a dihedral angle. The rotation will be around the bond between the 2nd and 3rd atoms.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current dihedral angle display
        self.dihedral_label = QLabel("")
        layout.addWidget(self.dihedral_label)
        
        # New dihedral angle input
        dihedral_layout = QHBoxLayout()
        dihedral_layout.addWidget(QLabel("New dihedral angle (degrees):"))
        self.dihedral_input = QLineEdit()
        self.dihedral_input.setPlaceholderText("180.0")
        dihedral_layout.addWidget(self.dihedral_input)
        layout.addLayout(dihedral_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Move:"))
        
        self.move_group_radio = QRadioButton("Atom 1,2,3: Fixed, Atom 4 group: Rotate")
        self.move_group_radio.setChecked(True)
        group_layout.addWidget(self.move_group_radio)
        
        self.move_atom_radio = QRadioButton("Atom 1,2,3: Fixed, Atom 4: Rotate atom only")
        group_layout.addWidget(self.move_atom_radio)
        
        self.both_groups_radio = QRadioButton("Central bond fixed: Both groups rotate equally")
        group_layout.addWidget(self.both_groups_radio)
        
        layout.addWidget(group_box)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker for DihedralDialog
        self.picker_connection = None
        self.enable_picking()
        
        # 事前選択された原子がある場合は初期表示を更新
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
        elif self.atom4_idx is None:
            self.atom4_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self.atom4_idx = None
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_changes()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
    
    def clear_selection(self):
        """選択をクリア"""
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None
        self.clear_atom_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        labels = ["1st", "2nd (bond start)", "3rd (bond end)", "4th"]
        colors = ["yellow", "yellow", "yellow", "yellow"]  # 全て黄色に統一
        
        for i, atom_idx in enumerate(selected_atoms):
            if atom_idx is not None:
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"{labels[i]}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color=colors[i],
                    always_visible=True
                )
                self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def update_display(self):
        """表示を更新"""
        selected_count = sum(x is not None for x in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx])
        
        if selected_count == 0:
            self.selection_label.setText("No atoms selected")
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input when no selection
            try:
                self.dihedral_input.clear()
            except Exception:
                pass
        elif selected_count < 4:
            selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
            
            display_parts = []
            for atom_idx in selected_atoms:
                if atom_idx is not None:
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    display_parts.append(f"{symbol}({atom_idx})")
                else:
                    display_parts.append("?")
            
            self.selection_label.setText(" - ".join(display_parts))
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input while selection is incomplete
            try:
                self.dihedral_input.clear()
            except Exception:
                pass
        else:
            selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
            
            display_parts = []
            for atom_idx in selected_atoms:
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                display_parts.append(f"{symbol}({atom_idx})")
            
            self.selection_label.setText(" - ".join(display_parts))
            
            # Calculate current dihedral angle
            current_dihedral = self.calculate_dihedral()
            self.dihedral_label.setText(f"Current dihedral: {current_dihedral:.2f}°")
            self.apply_button.setEnabled(True)
            # Update dihedral input box with current dihedral
            try:
                self.dihedral_input.setText(f"{current_dihedral:.2f}")
            except Exception:
                pass
    
    def calculate_dihedral(self):
        """現在の二面角を計算（正しい公式を使用）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        pos4 = np.array(conf.GetAtomPosition(self.atom4_idx))
        
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
    
    def apply_changes(self):
        """変更を適用"""
        if any(idx is None for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]):
            return
        
        try:
            new_dihedral = float(self.dihedral_input.text())
            if new_dihedral < -180 or new_dihedral > 180:
                QMessageBox.warning(self, "Invalid Input", "Dihedral angle must be between -180 and 180 degrees.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Apply the dihedral angle change
        self.adjust_dihedral(new_dihedral)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()

        # Undo状態を保存
        self.main_window.push_undo_state()
    
    def adjust_dihedral(self, new_dihedral_deg):
        """二面角を調整（改善されたアルゴリズム）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        pos4 = np.array(conf.GetAtomPosition(self.atom4_idx))
        
        # Current dihedral angle
        current_dihedral = self.calculate_dihedral()
        
        # Calculate rotation angle needed
        rotation_angle_deg = new_dihedral_deg - current_dihedral
        
        # Handle angle wrapping for shortest rotation
        if rotation_angle_deg > 180:
            rotation_angle_deg -= 360
        elif rotation_angle_deg < -180:
            rotation_angle_deg += 360
        
        rotation_angle_rad = np.radians(rotation_angle_deg)
        
        # Skip if no rotation needed
        if abs(rotation_angle_rad) < 1e-6:
            return
        
        # Rotation axis is the bond between atom2 and atom3
        rotation_axis = pos3 - pos2
        axis_length = np.linalg.norm(rotation_axis)
        
        if axis_length == 0:
            return  # Atoms are at the same position
        
        rotation_axis = rotation_axis / axis_length
        
        # Rodrigues' rotation formula implementation
        def rotate_point_around_axis(point, axis_point, axis_direction, angle):
            """Rotate a point around an axis using Rodrigues' formula"""
            # Translate point so axis passes through origin
            translated_point = point - axis_point
            
            # Apply Rodrigues' rotation formula
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            
            rotated = (translated_point * cos_a + 
                      np.cross(axis_direction, translated_point) * sin_a + 
                      axis_direction * np.dot(axis_direction, translated_point) * (1 - cos_a))
            
            # Translate back to original coordinate system
            return rotated + axis_point
        
        if self.both_groups_radio.isChecked():
            # Both groups rotate equally around the central bond (half angle each in opposite directions)
            half_rotation = rotation_angle_rad / 2
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom2_idx, exclude=self.atom3_idx)
            group4_atoms = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate group1 (atom1 side) by -half_rotation
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, -half_rotation)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Rotate group4 (atom4 side) by +half_rotation
            for atom_idx in group4_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, half_rotation)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.move_group_radio.isChecked():
            # Move the connected group containing atom4
            # Find all atoms connected to atom3 (excluding atom2 side)
            atoms_to_rotate = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate all atoms in the group
            for atom_idx in atoms_to_rotate:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, rotation_angle_rad)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        else:
            # Move only atom4
            new_pos4 = rotate_point_around_axis(pos4, pos2, rotation_axis, rotation_angle_rad)
            conf.SetAtomPosition(self.atom4_idx, new_pos4.tolist())
            self.main_window.atom_positions_3d[self.atom4_idx] = new_pos4
        
        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)
    
    def get_connected_group(self, start_atom, exclude=None):
        """指定された原子から連結されているグループを取得"""
        visited = set()
        to_visit = [start_atom]
        
        while to_visit:
            current = to_visit.pop()
            if current in visited or current == exclude:
                continue
            
            visited.add(current)
            
            # Get neighboring atoms
            atom = self.mol.GetAtomWithIdx(current)
            for bond in atom.GetBonds():
                other_idx = bond.GetOtherAtomIdx(current)
                if other_idx not in visited and other_idx != exclude:
                    to_visit.append(other_idx)
        
        return visited

class ConstrainedOptimizationDialog(Dialog3DPickingMixin, QDialog):
    """制約付き最適化ダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = []  # 順序が重要なのでリストを使用
        self.constraints = []  # (type, atoms_indices, value)
        self.constraint_labels = [] # 3Dラベルアクター
        self.init_ui()
        self.enable_picking()

        # MainWindowから既存の制約を読み込む
        if self.main_window.constraints_3d:
            self.constraint_table.blockSignals(True) # 読み込み中のシグナルをブロック
            try:
                # self.constraints には (Type, (Idx...), Value, Force) のタプル形式で読み込む
                for const_data in self.main_window.constraints_3d:
                    # 後方互換性のため、3要素または4要素の制約に対応
                    if len(const_data) == 4:
                        const_type, atom_indices, value, force_const = const_data
                    else:
                        const_type, atom_indices, value = const_data
                        force_const = 1.0e5  # デフォルト値
                    
                    # タプル化して内部リストに追加
                    self.constraints.append((const_type, tuple(atom_indices), value, force_const))
                    
                    row_count = self.constraint_table.rowCount()
                    self.constraint_table.insertRow(row_count)
                    
                    value_str = ""
                    if const_type == "Distance":
                        value_str = f"{value:.3f}"
                    else:
                        value_str = f"{value:.2f}"

                    # カラム 0 (Type)
                    item_type = QTableWidgetItem(const_type)
                    item_type.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 0, item_type)

                    # カラム 1 (Atom Indices)
                    item_indices = QTableWidgetItem(str(atom_indices))
                    item_indices.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 1, item_indices)

                    # カラム 2 (Value)
                    item_value = QTableWidgetItem(value_str)
                    item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 2, item_value)
                    
                    # カラム 3 (Force)
                    item_force = QTableWidgetItem(f"{force_const:.2e}")
                    item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 3, item_force)
            finally:
                self.constraint_table.blockSignals(False)

            # <<< MainWindowの現在の最適化設定を読み込み、デフォルトにする >>>
        try:
            # (修正) None の場合に備えてフォールバックを追加
            current_method_str = self.main_window.optimization_method or "MMFF_RDKIT"
            current_method = current_method_str.upper()
            
            # (修正) 比較順序を厳密化
            
            # 1. UFF_RDKIT
            if current_method == "UFF_RDKIT":
                self.ff_combo.setCurrentText("UFF")
            
            # 2. MMFF94_RDKIT (MMFF94)
            elif current_method == "MMFF94_RDKIT":
                self.ff_combo.setCurrentText("MMFF94")

            # 3. MMFF_RDKIT (MMFF94s) - これがデフォルトでもある
            elif current_method == "MMFF_RDKIT":
                self.ff_combo.setCurrentText("MMFF94s")

            # 4. (古い設定ファイルなどからのフォールバック)
            elif "UFF" in current_method:
                self.ff_combo.setCurrentText("UFF")
            elif "MMFF94S" in current_method:
                self.ff_combo.setCurrentText("MMFF94s")
            elif "MMFF94" in current_method: # MMFF94_RDKITも含むが、先で処理済み
                 self.ff_combo.setCurrentText("MMFF94")

            # 5. デフォルト
            else:
                self.ff_combo.setCurrentText("MMFF94s")
                
        except Exception as e:
            print(f"Could not set default force field: {e}")

    def init_ui(self):
        self.setWindowTitle("Constrained Optimization")
        self.setModal(False)
        self.resize(450, 500)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        layout = QVBoxLayout(self)

        # 1. 説明
        instruction_label = QLabel("Select 2-4 atoms to add a constraint. Select constraints in the table to remove them.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # 2. 最適化方法とForce Constant
        form_layout = QFormLayout()
        self.ff_combo = QComboBox()
        self.ff_combo.addItems(["MMFF94s", "MMFF94", "UFF"])
        form_layout.addRow("Force Field:", self.ff_combo)
        
        # Force Constant設定
        self.force_const_input = QLineEdit("1.0e5")
        self.force_const_input.setToolTip("Force constant for constraints (default: 1.0e5)")
        form_layout.addRow("Force Constant:", self.force_const_input)
        
        layout.addLayout(form_layout)
        
        # 3. 選択中の原子
        self.selection_label = QLabel("Selected atoms: None")
        layout.addWidget(self.selection_label)

        # 4. 制約の表
        self.constraint_table = QTableWidget()
        self.constraint_table.setColumnCount(4)
        self.constraint_table.setHorizontalHeaderLabels(["Type", "Atom Indices", "Value (Å or °)", "Force"])
        self.constraint_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        # 編集トリガーをダブルクリックなどに変更
        self.constraint_table.setEditTriggers(QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.SelectedClicked | QTableWidget.EditTrigger.EditKeyPressed)
        self.constraint_table.itemSelectionChanged.connect(self.show_constraint_labels)
        self.constraint_table.cellChanged.connect(self.on_cell_changed)

        self.constraint_table.setStyleSheet("""
            QTableWidget QLineEdit {
                background-color: white;
                color: black;
                border: none;
            }
        """)

        layout.addWidget(self.constraint_table)

        # 5. ボタン (Add / Remove)
        button_layout = QHBoxLayout()
        self.add_button = QPushButton("Add Constraint")
        self.add_button.clicked.connect(self.add_constraint)
        self.add_button.setEnabled(False)
        button_layout.addWidget(self.add_button)
        
        self.remove_button = QPushButton("Remove Selected")
        self.remove_button.clicked.connect(self.remove_constraint)
        self.remove_button.setEnabled(False)
        button_layout.addWidget(self.remove_button)
        layout.addLayout(button_layout)

        self.remove_all_button = QPushButton("Remove All")
        self.remove_all_button.clicked.connect(self.remove_all_constraints)
        button_layout.addWidget(self.remove_all_button)
        
        # 6. メインボタン (Optimize / Close)
        main_buttons = QHBoxLayout()
        main_buttons.addStretch()
        self.optimize_button = QPushButton("Optimize")
        self.optimize_button.clicked.connect(self.apply_optimization)
        main_buttons.addWidget(self.optimize_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        main_buttons.addWidget(close_button)
        layout.addLayout(main_buttons)

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            if len(self.selected_atoms) >= 4:
                self.selected_atoms.pop(0)  # 4つまで
            self.selected_atoms.append(atom_idx)
        
        self.show_selection_labels()
        self.update_selection_display()

    def update_selection_display(self):
        self.show_selection_labels()
        n = len(self.selected_atoms)

        atom_str = ", ".join(map(str, self.selected_atoms))
        prefix = ""
        can_add = False

        if n == 0:
            prefix = "Selected atoms: None"
            atom_str = ""  # atom_str を空にする
        elif n == 1:
            prefix = "Selected atoms: "
        elif n == 2:
            prefix = "Selected atoms: <b>Distance</b> "
            can_add = True
        elif n == 3:
            prefix = "Selected atoms: <b>Angle</b> "
            can_add = True
        elif n == 4:
            prefix = "Selected atoms: <b>Torsion</b> "
            can_add = True
        else: # n > 4
            prefix = "Selected atoms (max 4): "

        # ラベルテキストを設定
        if n == 0:
            self.selection_label.setText(prefix)
        else:
            self.selection_label.setText(f"{prefix}[{atom_str}]")

        # ボタンのテキストは常に固定
        self.add_button.setText("Add Constraint")
        # ボタンの有効状態を設定
        self.add_button.setEnabled(can_add)

    def add_constraint(self):
        n = len(self.selected_atoms)
        conf = self.mol.GetConformer()
        
        # Force Constantを取得
        try:
            force_const = float(self.force_const_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid Force Constant. Using default 1.0e5.")
            force_const = 1.0e5
        
        if n == 2:
            constraint_type = "Distance"
            value = conf.GetAtomPosition(self.selected_atoms[0]).Distance(conf.GetAtomPosition(self.selected_atoms[1]))
            value_str = f"{value:.3f}"
        elif n == 3:
            constraint_type = "Angle"
            value = rdMolTransforms.GetAngleDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        elif n == 4:
            constraint_type = "Torsion"
            value = rdMolTransforms.GetDihedralDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        else:
            return

        atom_indices = tuple(self.selected_atoms)
        
        # 既存の制約と重複チェック (原子インデックスが同じもの)
        for const in self.constraints:
            if const[0] == constraint_type and const[1] == atom_indices:
                QMessageBox.warning(self, "Warning", "This exact constraint already exists.")
                return

        self.constraints.append((constraint_type, atom_indices, value, force_const))
        
        # 表を更新
        # 表を更新
        row_count = self.constraint_table.rowCount()
        self.constraint_table.insertRow(row_count)

        # --- カラム 0 (Type) ---
        item_type = QTableWidgetItem(constraint_type)
        # 編集不可フラグを設定 (ItemIsEnabled | ItemIsSelectable)
        item_type.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 0, item_type)

        # --- カラム 1 (Atom Indices) ---
        item_indices = QTableWidgetItem(str(atom_indices))
        # 編集不可フラグを設定
        item_indices.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 1, item_indices)

        # --- カラム 2 (Value) ---
        item_value = QTableWidgetItem(value_str)
        item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # 編集可能フラグはデフォルトで有効 (ItemIsEnabled | ItemIsSelectable | ItemIsEditable)
        self.constraint_table.setItem(row_count, 2, item_value)
        
        # --- カラム 3 (Force) ---
        item_force = QTableWidgetItem(f"{force_const:.2e}")
        item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # 編集可能
        self.constraint_table.setItem(row_count, 3, item_force)

        # 選択をクリア
        self.selected_atoms.clear()
        self.update_selection_display()

    def remove_constraint(self):
        selected_rows = sorted(list(set(index.row() for index in self.constraint_table.selectedIndexes())), reverse=True)
        if not selected_rows:
            return

        self.constraint_table.blockSignals(True) 
        
        for row in selected_rows:
            self.constraints.pop(row)
            self.constraint_table.removeRow(row)
            
        self.constraint_table.blockSignals(False) 
            
        self.clear_constraint_labels()

    def remove_all_constraints(self):
        """全ての制約をクリアする"""
        if not self.constraints:
            return
            
        # 内部リストをクリア
        self.constraints.clear()
        
        # テーブルの行を全て削除
        self.constraint_table.blockSignals(True) 
        self.constraint_table.setRowCount(0)
        self.constraint_table.blockSignals(False) 
            
        # 3Dラベルをクリア
        self.clear_constraint_labels()
        
        # 選択ボタンを無効化
        self.remove_button.setEnabled(False)

    def show_constraint_labels(self):
        self.clear_constraint_labels()
        selected_items = self.constraint_table.selectedItems()
        if not selected_items:
            self.remove_button.setEnabled(False)
            return
            
        self.remove_button.setEnabled(True)
        
        # 選択された行の制約を取得 (最初の選択行のみ)
        try:
            row = selected_items[0].row()
            constraint_type, atom_indices, value, force_const = self.constraints[row]
        except (IndexError, TypeError, ValueError):
            # 古い形式の制約の場合は3要素でunpack
            try:
                constraint_type, atom_indices, value = self.constraints[row]
            except (IndexError, TypeError):
                return
        
        labels = []
        if constraint_type == "Distance":
            labels = ["A1", "A2"]
        elif constraint_type == "Angle":
            labels = ["A1", "A2 (V)", "A3"]
        elif constraint_type == "Torsion":
            labels = ["A1", "A2", "A3", "A4"]
        
        positions = []
        texts = []
        for i, atom_idx in enumerate(atom_indices):
            positions.append(self.main_window.atom_positions_3d[atom_idx])
            texts.append(labels[i])
        
        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions, texts,
                point_size=20, font_size=12, text_color='cyan', always_visible=True
            )
            self.constraint_labels.append(label_actor)

    def clear_constraint_labels(self):
        for label_actor in self.constraint_labels:
            try:
                self.main_window.plotter.remove_actor(label_actor)
            except:
                pass
        self.constraint_labels = []

    def apply_optimization(self):
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No valid 3D molecule found.")
            return

        ff_name = self.ff_combo.currentText()
        conf = self.mol.GetConformer()
        
        try:
            if ff_name.startswith("MMFF"):
                props = AllChem.MMFFGetMoleculeProperties(self.mol, mmffVariant=ff_name)
                ff = AllChem.MMFFGetMoleculeForceField(self.mol, props, confId=0)
                add_dist_constraint = ff.MMFFAddDistanceConstraint
                add_angle_constraint = ff.MMFFAddAngleConstraint
                add_torsion_constraint = ff.MMFFAddTorsionConstraint
            else: # UFF
                ff = AllChem.UFFGetMoleculeForceField(self.mol, confId=0)
                add_dist_constraint = ff.UFFAddDistanceConstraint
                add_angle_constraint = ff.UFFAddAngleConstraint
                add_torsion_constraint = ff.UFFAddTorsionConstraint

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initialize force field {ff_name}: {e}")
            return

        # 制約を追加
        try:
            for constraint in self.constraints:
                # 後方互換性のため、4要素または3要素の制約に対応
                if len(constraint) == 4:
                    const_type, atoms, value, force_const = constraint
                else:
                    const_type, atoms, value = constraint
                    force_const = 1.0e5  # デフォルト値
                
                if const_type == "Distance":
                    # C++ signature: (self, idx1, idx2, bool relative, minLen, maxLen, forceConst)
                    add_dist_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        False, 
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
                elif const_type == "Angle":
                    # C++ signature: (self, idx1, idx2, idx3, bool relative, minDeg, maxDeg, forceConst)
                    add_angle_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        int(atoms[2]),
                        False,  
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
                elif const_type == "Torsion":
                    # C++ signature: (self, idx1, idx2, idx3, idx4, bool relative, minDeg, maxDeg, forceConst)
                    add_torsion_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        int(atoms[2]), 
                        int(atoms[3]),
                        False, 
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to add constraints: {e}")
            print(e)
            return

        # 最適化の実行
        try:
            self.main_window.statusBar().showMessage(f"Running constrained {ff_name} optimization...")
            ff.Minimize(maxIts=20000)
            
            # 最適化後の座標をメインウィンドウの numpy 配列に反映
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                self.main_window.atom_positions_3d[i] = [pos.x, pos.y, pos.z]
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            self.main_window.update_chiral_labels()
            self.main_window.push_undo_state()
            self.main_window.statusBar().showMessage("Constrained optimization finished.")

            try:
                constrained_method_name = f"Constrained_{ff_name}"
                self.main_window.last_successful_optimization_method = constrained_method_name
            except Exception as e:
                print(f"Failed to set last_successful_optimization_method: {e}")

            # (修正) 最適化成功時にも制約リストをMainWindowに保存 (reject と同じロジック)
            try:
                # JSON互換のため、タプルをリストに変換して保存
                json_safe_constraints = []
                for const in self.constraints:
                    # 4要素の制約（Type, Indices, Value, Force）
                    if len(const) == 4:
                        json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                    else:
                        # 古い形式の場合は3要素にデフォルトのForceを追加
                        json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
                
                # 変更があった場合のみ MainWindow を更新
                if self.main_window.constraints_3d != json_safe_constraints:
                    self.main_window.constraints_3d = json_safe_constraints
                    self.main_window.has_unsaved_changes = True # 制約の変更も「未保存」扱い
                    self.main_window.update_window_title()
                    
            except Exception as e:
                print(f"Failed to save constraints post-optimization: {e}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Optimization failed: {e}")

    def closeEvent(self, event):
        self.reject()
        event.accept()
    
    def reject(self):
        self.clear_constraint_labels()
        self.clear_selection_labels()
        self.disable_picking()

        # ダイアログを閉じる際に現在の制約リストをMainWindowに保存
        try:
            # JSON互換のため、タプルをリストに変換して保存
            json_safe_constraints = []
            for const in self.constraints:
                # (Type, (Idx...), Value, Force) -> [Type, [Idx...], Value, Force]
                if len(const) == 4:
                    json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                else:
                    # 古い形式の場合は3要素にデフォルトのForceを追加
                    json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
            
            # 変更があった場合のみ MainWindow を更新
            if self.main_window.constraints_3d != json_safe_constraints:
                self.main_window.constraints_3d = json_safe_constraints
                self.main_window.has_unsaved_changes = True # 制約の変更も「未保存」扱い
                self.main_window.update_window_title()
                
        except Exception as e:
            print(f"Failed to save constraints to main window: {e}")

        super().reject()

    def clear_selection(self):
        """選択をクリア (原子以外をクリックした時にMixinから呼ばれる)"""
        self.selected_atoms.clear()
        self.clear_selection_labels()
        self.update_selection_display()

    def show_selection_labels(self):
        """選択された原子にラベルを表示"""
        self.clear_selection_labels()
        
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        if not hasattr(self.main_window, 'atom_positions_3d') or self.main_window.atom_positions_3d is None:
            return  # 3D座標データがない場合は何もしない
            
        max_idx = len(self.main_window.atom_positions_3d) - 1
        positions = []
        texts = []
        
        for i, atom_idx in enumerate(self.selected_atoms):
            if atom_idx is not None and 0 <= atom_idx <= max_idx:
                positions.append(self.main_window.atom_positions_3d[atom_idx])
                texts.append(f"A{i+1}")
            elif atom_idx is not None:
                # インデックスが無効な場合はログ（デバッグ用）
                print(f"Warning: Invalid atom index {atom_idx} in show_selection_labels")
        
        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions, texts,
                point_size=20, font_size=12, text_color='yellow', always_visible=True
            )
            # add_point_labelsがリストを返す場合も考慮
            if isinstance(label_actor, list):
                self.selection_labels.extend(label_actor)
            else:
                self.selection_labels.append(label_actor)

    def clear_selection_labels(self):
        """選択ラベル(A1, A2...)をクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []

    def on_cell_changed(self, row, column):
        """テーブルのセルが編集されたときに内部データを更新する"""
        
        # "Value" 列 (カラムインデックス 2) と "Force" 列 (カラムインデックス 3) のみ対応
        if column not in [2, 3]:
            return

        try:
            # 変更されたアイテムからテキストを取得
            item = self.constraint_table.item(row, column)
            if not item:
                return
            
            new_value_str = item.text()
            new_value = float(new_value_str)
            
            # 内部の constraints リストを更新
            old_constraint = self.constraints[row]
            
            # 後方互換性のため、3要素または4要素の制約に対応
            if len(old_constraint) == 4:
                if column == 2:  # Value列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], new_value, old_constraint[3])
                elif column == 3:  # Force列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], old_constraint[2], new_value)
            else:
                # 古い3要素形式の場合
                if column == 2:  # Value列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], new_value, 1.0e5)
                elif column == 3:  # Force列（新規追加）
                    self.constraints[row] = (old_constraint[0], old_constraint[1], old_constraint[2], new_value)

        except (ValueError, TypeError):
            # 不正な値（数値以外）が入力された場合
            # 元の値をテーブルに戻す
            self.constraint_table.blockSignals(True)
            
            if column == 2:  # Value列
                old_value = self.constraints[row][2]
                if self.constraints[row][0] == "Distance":
                    item.setText(f"{old_value:.3f}")
                else:
                    item.setText(f"{old_value:.2f}")
            elif column == 3:  # Force列
                old_force = self.constraints[row][3] if len(self.constraints[row]) == 4 else 1.0e5
                item.setText(f"{old_force:.2e}")
            
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.constraint_table.blockSignals(False)
            
            QMessageBox.warning(self, "Invalid Value", "Please enter a valid floating-point number.")
        except IndexError:
            # constraints リストとテーブルが同期していない場合（通常発生しない）
            pass
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理 (Delete/Backspaceで制約を削除, Enterで最適化)"""
        key = event.key()
        
        # DeleteキーまたはBackspaceキーが押されたかチェック
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            # テーブルがフォーカスを持っているか、またはアイテムが選択されているか確認
            if self.constraint_table.hasFocus() or len(self.constraint_table.selectedIndexes()) > 0:
                self.remove_constraint()
                event.accept()
                return

        # Enter/Returnキーが押されたかチェック (最適化を実行)
        if key == Qt.Key.Key_Return or key == Qt.Key.Key_Enter:
            # テーブルが編集中でないことを確認（セルの編集中にEnterを押した場合）
            if self.constraint_table.state() != QAbstractItemView.State.EditingState:
                if self.optimize_button.isEnabled():
                    self.apply_optimization()
                event.accept()
                return
            
        # それ以外のキーはデフォルトの処理
        super().keyPressEvent(event)


# --- Application Execution ---
if __name__ == '__main__':
    main()

