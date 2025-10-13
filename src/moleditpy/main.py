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
VERSION = '1.8.0'

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

from collections import deque

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QSplitter, QGraphicsView, QGraphicsScene, QGraphicsItem,
    QToolBar, QStatusBar, QGraphicsTextItem, QGraphicsLineItem, QDialog, QGridLayout,
    QFileDialog, QSizePolicy, QLabel, QLineEdit, QToolButton, QMenu, QMessageBox, QInputDialog,
    QColorDialog, QCheckBox, QSlider, QFormLayout, QRadioButton, QComboBox, QListWidget, QListWidgetItem, QButtonGroup, QTabWidget
)

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QPainterPath, QPainterPathStroker, QFontMetrics, QFontMetricsF, QKeySequence, QTransform, QCursor, QPixmap, QIcon, QShortcut, QDesktopServices, QImage
)


from PyQt6.QtCore import Qt, QPointF, QRectF, QLineF, QObject, QThread, pyqtSignal, QEvent, QMimeData, QByteArray, QUrl, QTimer

from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


# Open Babel Python binding (optional; required for fallback)
from openbabel import pybel

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

pt = Chem.GetPeriodicTable()
VDW_RADII = {pt.GetElementSymbol(i): pt.GetRvdw(i) * 0.3 for i in range(1, 119)}

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
                                self.on_atom_picked(int(closest_atom_idx))
                                return True  # Consume event
                
                # 原子以外をクリックした場合は選択をクリア（Measurementモードと同じロジック）
                if hasattr(self, 'clear_selection'):
                    self.clear_selection()
                return True  # Consume event
                    
            except Exception as e:
                print(f"Error in eventFilter: {e}")
                
        return super().eventFilter(obj, event)
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, 'picking_enabled') and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
    
    def try_alternative_picking(self, x, y):
        """代替のピッキング方法（使用しない）"""
        pass

class TranslationDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.reference_atom_idx = None
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Translation")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click an atom as reference point, then enter target coordinates.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atom display
        self.selection_label = QLabel("No atom selected")
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
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Translation")
        self.apply_button.clicked.connect(self.apply_translation)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        self.reference_atom_idx = atom_idx
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
        if self.reference_atom_idx is None:
            self.selection_label.setText("No atom selected")
            self.apply_button.setEnabled(False)
        else:
            # 分子の有効性チェック
            if not self.mol or self.mol.GetNumConformers() == 0:
                self.selection_label.setText("Error: No valid molecule or conformer")
                self.apply_button.setEnabled(False)
                return
            
            try:
                symbol = self.mol.GetAtomWithIdx(self.reference_atom_idx).GetSymbol()
                conf = self.mol.GetConformer()
                pos = conf.GetAtomPosition(self.reference_atom_idx)
                self.selection_label.setText(f"Reference: {symbol}({self.reference_atom_idx}) at ({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})")
                self.apply_button.setEnabled(True)
            except Exception as e:
                self.selection_label.setText(f"Error accessing atom data: {str(e)}")
                self.apply_button.setEnabled(False)
    
    def apply_translation(self):
        """平行移動を適用"""
        if self.reference_atom_idx is None:
            QMessageBox.warning(self, "Warning", "Please select a reference atom.")
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
            # 参照原子の現在位置を取得
            conf = self.mol.GetConformer()
            current_pos = np.array(conf.GetAtomPosition(self.reference_atom_idx))
            target_pos = np.array([target_x, target_y, target_z])
            
            # 移動ベクトルを計算
            translation_vector = target_pos - current_pos
            
            # Undo状態を保存
            self.main_window.push_undo_state()
            
            # 全原子を平行移動
            for i in range(self.mol.GetNumAtoms()):
                atom_pos = np.array(conf.GetAtomPosition(i))
                new_pos = atom_pos + translation_vector
                conf.SetAtomPosition(i, new_pos.tolist())
                self.main_window.atom_positions_3d[i] = new_pos
            
            # 3D表示を更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # キラルラベルを更新
            self.main_window.update_chiral_labels()
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply translation: {str(e)}")
    
    def clear_selection(self):
        """選択をクリア"""
        self.reference_atom_idx = None
        self.clear_atom_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        if self.reference_atom_idx is not None:
            pos = self.main_window.atom_positions_3d[self.reference_atom_idx]
            label_text = "Ref"
            
            # ラベルを追加
            label_actor = self.main_window.plotter.add_point_labels(
                [pos], [label_text], 
                point_size=20, 
                font_size=12,
                text_color='cyan',
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

class SymmetrizeDialog(QDialog):
    """分子構造の対称化機能を提供するダイアログ"""
    """ The parameters have not been checked for accuracy. Temporary measure. Under Development"""

    # 黄金比 (正二十面体群の記述に使用)
    PHI = (1 + np.sqrt(5)) / 2

    POINT_GROUPS = {
        # ===============================================================
        # 1. 低対称性群 (Low Symmetry Groups)
        # ===============================================================
        "C1": {
            "name": "C1 (No symmetry)", 
            "operations": [np.eye(3)] # E
        },
        "Ci": {
            "name": "Ci (Inversion center)", 
            "operations": [
                np.eye(3),      # E
                -np.eye(3)      # i
            ]
        },
        "Cs": {
            "name": "Cs (Mirror plane)", 
            "operations": [
                np.eye(3),      # E
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])  # σh (xy)
            ]
        },

        # ===============================================================
        # 2. 単一軸を持つ群 (Groups with a single axis)
        # ===============================================================
        # Cn 群 (カイラル)
        "C2": {
            "name": "C2 (Rotation)", 
            "operations": [
                np.eye(3), 
                np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])  # C2(z)
            ]
        },
        "C3": {
            "name": "C3 (Rotation)", 
            "operations": [
                np.eye(3),
                np.array([[-0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -0.5, 0], [0, 0, 1]]), # C3
                np.array([[-0.5,  np.sqrt(3)/2, 0], [-np.sqrt(3)/2, -0.5, 0], [0, 0, 1]])  # C3^2
            ]
        },
        # Cnh 群
        "C2h": {
            "name": "C2h", 
            "operations": [
                np.eye(3),                                      # E
                np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),  # C2(z)
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),  # σh
                -np.eye(3)                                      # i
            ]
        },
        "C3h": {
            "name": "C3h",
            "operations": [
                np.array([[-0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -0.5, 0], [0, 0, 1]]), # C3
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),                          # σh
            ]
        },
        # Cnv 群
        "C2v": {
            "name": "C2v", 
            "operations": [
                np.eye(3),
                np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),  # C2(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),  # σv(xz)
                np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])   # σv(yz)
            ]
        },
        "C3v": {
            "name": "C3v", 
            "operations": [
                 np.array([[-0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -0.5, 0], [0, 0, 1]]), # C3
                 np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])                           # σv(xz)
            ]
        },

        # ===============================================================
        # 3. 二面体群 (Dihedral Groups)
        # ===============================================================
        # Dn 群 (カイラル)
        "D2": {
            "name": "D2", 
            "operations": [
                np.eye(3),
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # C2(x)
                np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]),  # C2(y)
                np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])   # C2(z)
            ]
        },
        "D3": {
            "name": "D3", 
            "operations": [
                np.array([[-0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -0.5, 0], [0, 0, 1]]), # C3(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])                           # C2(x)
            ]
        },
        # Dnh 群
        "D2h": {
            "name": "D2h", 
            "operations": [
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # C2(x)
                np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]),  # C2(y)
                -np.eye(3)                                      # i
            ]
        },
        "D3h": {
            "name": "D3h", 
            "operations": [
                np.array([[-0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, -0.5, 0], [0, 0, 1]]), # C3(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),                           # C2(x)
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])                            # σh
            ]
        },
        "D4h": {
            "name": "D4h",
            "operations": [
                np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),    # C4(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # C2(x)
                -np.eye(3)                                      # i
            ]
        },
        "D5h": {
            "name": "D5h",
            "operations": [
                # C5(z)
                np.array([[np.cos(2*np.pi/5), -np.sin(2*np.pi/5), 0], 
                          [np.sin(2*np.pi/5), np.cos(2*np.pi/5), 0], 
                          [0, 0, 1]]),
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # C2(x)
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),  # σh
            ]
        },
        "D6h": {
            "name": "D6h",
            "operations": [
                np.array([[0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, 0.5, 0], [0, 0, 1]]), # C6(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),                        # C2(x)
                -np.eye(3)                                                            # i
            ]
        },
        # Dnd 群
        "D2d": {
            "name": "D2d",
            "operations": [
                np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]]),   # S4(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])   # C2(x)
            ]
        },
        "D3d": {
            "name": "D3d",
            "operations": [
                np.array([[0.5, -np.sqrt(3)/2, 0], [np.sqrt(3)/2, 0.5, 0], [0, 0, -1]]), # S6(z)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])                          # C2(x)
            ]
        },

        # ===============================================================
        # 4. 高対称性群 (Cubic and Icosahedral Groups)
        # ===============================================================
        # Td 群 (正四面体)
        "Td": {
            "name": "Td (Tetrahedral)", 
            "operations": [
                np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),    # E (恒等操作)
                np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),    # C3(111)
                np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),    # C3(111)^2
                np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),   # S4(z)
                np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),   # S4(z)^2 = C2(z)
                np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]]),  # S4(z)^3
                np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]]),   # S4(y)
                np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]),   # S4(y)^3
                np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),   # S4(x)
                np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]]),   # S4(x)^3
                np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),  # C2(z)
                np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]]),  # C2(y)
                np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # C2(x)
            ]
        },
        # Oh 群 (正八面体 / 立方体)
        "Oh": {
            "name": "Oh (Octahedral)",
            "operations": [
                np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),    # C3(111)
                np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),   # C4(z, inv)
                -np.eye(3)                                      # i
            ]
        },
        # Ih 群 (正二十面体)
        "Ih": {
            "name": "Ih (Icosahedral)",
            "operations": [
                # C5 z-phi plane
                np.array([[ (PHI-1)/2, -PHI/2,  0.5],
                          [      PHI/2,  (PHI-1)/2, -0.5],
                          [       -0.5,      0.5,  PHI/2]]),
                # C3 111
                 np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
                 -np.eye(3) # i
            ]
        }
    }
    
    def __init__(self, mol, main_window, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.main_window = main_window
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Symmetrize Molecule")
        self.setModal(True)
        self.setFixedSize(450, 350)
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(
            "Select a point group and tolerance to automatically symmetrize the molecular structure. "
            "The algorithm will adjust atomic positions to enforce the selected symmetry."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        layout.addWidget(QLabel(""))  # Spacer
        
        # Point group selection
        point_group_layout = QFormLayout()
        
        # Point group selection with auto-detect button
        pg_selection_layout = QHBoxLayout()
        self.point_group_combo = QComboBox()
        for key, value in self.POINT_GROUPS.items():
            self.point_group_combo.addItem(value["name"], key)
        self.point_group_combo.setCurrentIndex(0)  # Default to C1
        pg_selection_layout.addWidget(self.point_group_combo)
        
        self.auto_detect_button = QPushButton("Auto-Detect")
        self.auto_detect_button.clicked.connect(self.auto_detect_symmetry)
        self.auto_detect_button.setToolTip("Automatically detect the most suitable point group for current structure")
        pg_selection_layout.addWidget(self.auto_detect_button)
        
        point_group_layout.addRow("Point Group:", pg_selection_layout)
        
        # Tolerance input
        self.tolerance_input = QLineEdit("0.1")
        self.tolerance_input.setToolTip("Maximum allowed displacement (Angstroms) when applying symmetry operations")
        point_group_layout.addRow("Tolerance (Å):", self.tolerance_input)
        
        layout.addLayout(point_group_layout)
        
        layout.addWidget(QLabel(""))  # Spacer
        
        # Preview information
        self.info_label = QLabel("Select parameters and click Apply to symmetrize the structure.")
        self.info_label.setWordWrap(True)
        self.info_label.setStyleSheet("QLabel { color: #666; font-style: italic; }")
        layout.addWidget(self.info_label)
        
        layout.addStretch()
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.preview_button = QPushButton("Preview Symmetry")
        self.preview_button.clicked.connect(self.preview_symmetry)
        self.preview_button.setToolTip("Analyze current structure and show symmetry information")
        button_layout.addWidget(self.preview_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Symmetrization")
        self.apply_button.clicked.connect(self.apply_symmetrization)
        self.apply_button.setDefault(True)
        button_layout.addWidget(self.apply_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def get_selected_point_group(self):
        """選択されたポイントグループの情報を取得"""
        key = self.point_group_combo.currentData()
        return key, self.POINT_GROUPS[key]
    
    def get_tolerance(self):
        """入力されたトレランス値を取得"""
        try:
            return float(self.tolerance_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter a valid tolerance value.")
            return None
    
    def preview_symmetry(self):
        """現在の構造の対称性を分析してプレビュー表示（高度な解析付き）"""
        tolerance = self.get_tolerance()
        if tolerance is None:
            return
        
        # 分子の有効性チェック
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Warning", "No valid molecule or conformer available.")
            return
        
        key, point_group = self.get_selected_point_group()
        
        try:
            # 現在の分子の座標を取得
            conf = self.mol.GetConformer()
            positions = np.array([conf.GetAtomPosition(i) for i in range(self.mol.GetNumAtoms())])
            
            # 高度な対称性分析を実行
            analysis_result = self.analyze_symmetry(positions, point_group["operations"], tolerance)
            
            # 等価原子グループ情報を取得
            try:
                symmetry_classes = self.get_molecular_symmetry_classes()
                equivalent_groups = self.group_equivalent_atoms(symmetry_classes)
                
                # グループ情報を整理
                group_info = []
                equivalent_count = 0
                
                for i, group in enumerate(equivalent_groups):
                    if len(group) > 1:
                        equivalent_count += 1
                        symbols = [self.mol.GetAtomWithIdx(idx).GetSymbol() for idx in group]
                        group_info.append(f"Group {equivalent_count}: {len(group)} {symbols[0]} atoms (indices: {group})")
                    else:
                        symbols = [self.mol.GetAtomWithIdx(idx).GetSymbol() for idx in group]
                        group_info.append(f"Single: {symbols[0]} atom (index: {group[0]})")
                
                if equivalent_count > 0:
                    group_text = f"Found {equivalent_count} equivalent atom groups:\n" + "\n".join(group_info)
                else:
                    group_text = "No equivalent atom groups found.\nAll atoms are chemically unique.\n" + "\n".join(group_info)
                
            except Exception as e:
                group_text = f"Symmetry analysis error: {str(e)}\nUsing fallback simple method."
            
            # 結果をダイアログに表示
            info_text = f"Point Group: {point_group['name']}\n"
            info_text += f"Tolerance: {tolerance} Å\n"
            info_text += f"Atoms to be moved: {analysis_result['atoms_to_move']}\n"
            info_text += f"Max displacement: {analysis_result['max_displacement']:.3f} Å\n\n"
            info_text += "Equivalent Atom Groups:\n" + group_text + "\n\n"
            
            if analysis_result['atoms_to_move'] == 0:
                info_text += "Structure already satisfies the selected symmetry."
            else:
                info_text += "Ready to apply symmetrization."
            
            self.info_label.setText(info_text)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to analyze symmetry: {str(e)}")
            print(f"Symmetry analysis error: {e}")
            import traceback
            traceback.print_exc()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            self.apply_symmetrization()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def apply_symmetrization(self):
        """対称化を適用"""
        tolerance = self.get_tolerance()
        if tolerance is None:
            return
        
        key, point_group = self.get_selected_point_group()
        
        try:
            # Undo状態を保存（操作前の状態のみ）
            self.main_window.push_undo_state()
            
            # 現在の分子の座標を取得
            conf = self.mol.GetConformer()
            original_positions = np.array([conf.GetAtomPosition(i) for i in range(self.mol.GetNumAtoms())])
            
            for i, pos in enumerate(original_positions):
                print(f"  Atom {i}: {pos}")
            
            # 対称化を適用
            new_positions = self.apply_symmetry_operations(
                original_positions, point_group["operations"], tolerance
            )
            
            # 新しい座標を分子に適用（3D座標のみ）
            for i, new_pos in enumerate(new_positions):
                conf.SetAtomPosition(i, new_pos.tolist())
                self.main_window.atom_positions_3d[i] = new_pos
            
            # 3D表示のみを更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # 成功メッセージ（ダイアログを閉じる前に表示）
            QMessageBox.information(
                self, "Success", 
                f"Molecular structure has been symmetrized according to {point_group['name']} point group."
            )
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply symmetrization: {str(e)}")
            print(f"Symmetrization error: {e}")
            import traceback
            traceback.print_exc()
    
    def analyze_symmetry(self, positions, operations, tolerance):
        """対称性を分析し、必要な変更を計算（高度なアルゴリズム）"""
        atoms_to_move = 0
        max_displacement = 0.0
        
        # RDKitを使用して分子の対称性と等価原子グループを取得
        try:
            # 分子の対称性解析
            symmetry_classes = self.get_molecular_symmetry_classes()
            equivalent_atom_groups = self.group_equivalent_atoms(symmetry_classes)
            
            # 各等価原子グループに対して対称化を適用
            centroid = np.mean(positions, axis=0)
            centered_positions = positions - centroid
            
            for group_atoms in equivalent_atom_groups:
                group_positions = centered_positions[group_atoms]
                
                # グループ内での理想的な対称位置を計算
                ideal_positions = self.calculate_ideal_symmetric_positions(
                    group_positions, operations, tolerance
                )
                
                # 各原子の移動距離を計算
                for i, atom_idx in enumerate(group_atoms):
                    displacement = np.linalg.norm(ideal_positions[i] - centered_positions[atom_idx])
                    if displacement > tolerance:
                        atoms_to_move += 1
                        max_displacement = max(max_displacement, displacement)
        
        except Exception as e:
            print(f"Advanced symmetry analysis failed, falling back to simple method: {e}")
            # フォールバック: 従来の単純な方法
            return self.analyze_symmetry_simple(positions, operations, tolerance)
        
        return {
            'atoms_to_move': atoms_to_move,
            'max_displacement': max_displacement
        }
    
    def get_molecular_symmetry_classes(self):
        """RDKitを使用して分子の対称性クラスを取得"""
        try:
            from rdkit.Chem import rdMolDescriptors
            
            # RDKitの正しいAPI名を試行
            try:
                # 新しいバージョンのRDKit
                canonical_ranks = list(rdMolDescriptors.GetAtomSymmetryClasses(self.mol))
                print(f"Debug: Using GetAtomSymmetryClasses - Symmetry classes: {canonical_ranks}")
                return canonical_ranks
            except AttributeError:
                # 代替手法: Canonical SMILES順序を利用
                try:
                    from rdkit.Chem import Descriptors
                    
                    # 分子のcanonical atom orderingを取得
                    mol_copy = Chem.Mol(self.mol)
                    canonical_order = tuple(mol_copy.GetPropsAsDict().get('_smilesAtomOutputOrder', range(mol_copy.GetNumAtoms())))
                    
                    # Morgan fingerprintを使用して原子の環境を比較
                    from rdkit.Chem import rdMolDescriptors
                    atom_invariants = []
                    
                    for atom in mol_copy.GetAtoms():
                        # 原子の化学環境を記述する不変量を計算
                        invariant = (
                            atom.GetAtomicNum(),
                            atom.GetDegree(),
                            atom.GetFormalCharge(),
                            atom.GetHybridization(),
                            atom.GetTotalNumHs(),
                            atom.IsInRing()
                        )
                        atom_invariants.append(invariant)
                    
                    # 同じ不変量を持つ原子に同じクラス番号を割り当て
                    unique_invariants = list(set(atom_invariants))
                    symmetry_classes = []
                    
                    for invariant in atom_invariants:
                        class_id = unique_invariants.index(invariant)
                        symmetry_classes.append(class_id)
                    
                    print(f"Debug: Using Morgan-based approach - Symmetry classes: {symmetry_classes}")
                    return symmetry_classes
                    
                except Exception as e2:
                    print(f"Morgan-based approach also failed: {e2}")
                    raise e
            
        except Exception as e:
            print(f"Failed to get molecular symmetry classes: {e}")
            print("Debug: Falling back to chemical-based symmetry analysis")
            
            # 手動でよくある分子パターンをチェック
            if self.is_methane_like():
                print("Debug: Detected methane-like molecule, applying manual grouping")
                return self.get_methane_symmetry_classes()
            elif self.is_water_like():
                print("Debug: Detected water-like molecule, applying manual grouping")
                return self.get_water_symmetry_classes()
            elif self.is_ammonia_like():
                print("Debug: Detected ammonia-like molecule, applying manual grouping")
                return self.get_ammonia_symmetry_classes()
            
            # フォールバック: 化学的知識に基づく推定
            return self.get_chemical_symmetry_classes()
    
    def is_methane_like(self):
        """メタン様分子（CH4, CCl4など）かどうか判定"""
        if self.mol.GetNumAtoms() != 5:
            return False
        
        # 中心原子（通常は炭素）と4つの等価な原子
        atoms = list(self.mol.GetAtoms())
        center_candidates = [atom for atom in atoms if atom.GetDegree() == 4]
        
        if len(center_candidates) != 1:
            return False
            
        center_atom = center_candidates[0]
        neighbors = [atom.GetSymbol() for atom in center_atom.GetNeighbors()]
        
        # 4つの隣接原子がすべて同じ元素かチェック
        return len(set(neighbors)) == 1 and len(neighbors) == 4
    
    def is_water_like(self):
        """水様分子（H2O）かどうか判定"""
        if self.mol.GetNumAtoms() != 3:
            return False
        
        atoms = list(self.mol.GetAtoms())
        center_candidates = [atom for atom in atoms if atom.GetDegree() == 2]
        
        if len(center_candidates) != 1:
            return False
            
        center_atom = center_candidates[0]
        neighbors = [atom.GetSymbol() for atom in center_atom.GetNeighbors()]
        
        return len(set(neighbors)) == 1 and len(neighbors) == 2
    
    def is_ammonia_like(self):
        """アンモニア様分子（NH3）かどうか判定"""
        if self.mol.GetNumAtoms() != 4:
            return False
        
        atoms = list(self.mol.GetAtoms())
        center_candidates = [atom for atom in atoms if atom.GetDegree() == 3]
        
        if len(center_candidates) != 1:
            return False
            
        center_atom = center_candidates[0]
        neighbors = [atom.GetSymbol() for atom in center_atom.GetNeighbors()]
        
        return len(set(neighbors)) == 1 and len(neighbors) == 3
    
    def get_methane_symmetry_classes(self):
        """メタン様分子の対称性クラス"""
        # 中心原子のインデックスを見つける
        center_idx = None
        for i, atom in enumerate(self.mol.GetAtoms()):
            if atom.GetDegree() == 4:
                center_idx = i
                break
        
        # 中心原子は独自のクラス、4つの隣接原子は同じクラス
        symmetry_classes = [1] * self.mol.GetNumAtoms()  # 隣接原子のクラス
        if center_idx is not None:
            symmetry_classes[center_idx] = 0  # 中心原子のクラス
        
        return symmetry_classes
    
    def get_water_symmetry_classes(self):
        """水様分子の対称性クラス"""
        center_idx = None
        for i, atom in enumerate(self.mol.GetAtoms()):
            if atom.GetDegree() == 2:
                center_idx = i
                break
        
        symmetry_classes = [1] * self.mol.GetNumAtoms()  # 隣接原子のクラス
        if center_idx is not None:
            symmetry_classes[center_idx] = 0  # 中心原子のクラス
        
        return symmetry_classes
    
    def get_ammonia_symmetry_classes(self):
        """アンモニア様分子の対称性クラス"""
        center_idx = None
        for i, atom in enumerate(self.mol.GetAtoms()):
            if atom.GetDegree() == 3:
                center_idx = i
                break
        
        symmetry_classes = [1] * self.mol.GetNumAtoms()  # 隣接原子のクラス
        if center_idx is not None:
            symmetry_classes[center_idx] = 0  # 中心原子のクラス
        
        return symmetry_classes
    
    def get_chemical_symmetry_classes(self):
        """化学的知識に基づく対称性クラス推定"""
        num_atoms = self.mol.GetNumAtoms()
        
        # 分子の3D座標を使用してより精密な解析
        try:
            conf = self.mol.GetConformer()
            positions = np.array([list(conf.GetAtomPosition(i)) for i in range(num_atoms)])
            
            # 距離行列ベースの解析
            symmetry_classes = self.analyze_by_distance_matrix(positions)
            if symmetry_classes:
                print(f"Debug: Distance-matrix-based symmetry classes: {symmetry_classes}")
                return symmetry_classes
        except Exception as e:
            print(f"Debug: Distance matrix analysis failed: {e}")
        
        # フォールバック: 各原子を元素と結合数で分類
        atom_signatures = []
        for atom in self.mol.GetAtoms():
            # より詳細な原子環境の記述
            neighbors = atom.GetNeighbors()
            neighbor_symbols = sorted([n.GetSymbol() for n in neighbors])
            
            signature = (
                atom.GetSymbol(),           # 元素記号
                atom.GetDegree(),          # 結合数
                atom.GetFormalCharge(),    # 電荷
                tuple(neighbor_symbols)    # 隣接原子の元素記号（ソート済み）
            )
            atom_signatures.append(signature)
        
        # 同じsignatureの原子に同じクラス番号を割り当て
        unique_signatures = list(set(atom_signatures))
        symmetry_classes = []
        
        for signature in atom_signatures:
            class_id = unique_signatures.index(signature)
            symmetry_classes.append(class_id)
        
        print(f"Debug: Chemical-based symmetry classes: {symmetry_classes}")
        print(f"Debug: Unique signatures: {unique_signatures}")
        return symmetry_classes
    
    def analyze_by_distance_matrix(self, positions):
        """距離行列を使用した対称性解析"""
        num_atoms = len(positions)
        
        # 各原子から他の全原子への距離ベクトルを計算
        distance_patterns = []
        
        for i in range(num_atoms):
            # 原子iから他の全原子への距離を計算
            distances = []
            for j in range(num_atoms):
                if i != j:
                    dist = np.linalg.norm(positions[i] - positions[j])
                    distances.append(round(dist, 3))  # 丸めて比較可能にする
            
            # 距離をソートして標準化
            distances.sort()
            
            # 元素記号と組み合わせて特徴ベクトルを作成
            atom_symbol = self.mol.GetAtomWithIdx(i).GetSymbol()
            pattern = (atom_symbol, tuple(distances))
            distance_patterns.append(pattern)
        
        # 同じパターンを持つ原子をグループ化
        unique_patterns = []
        symmetry_classes = []
        
        for pattern in distance_patterns:
            # 既存パターンとの近似比較
            matched_class = None
            for idx, unique_pattern in enumerate(unique_patterns):
                if self.patterns_are_similar(pattern, unique_pattern):
                    matched_class = idx
                    break
            
            if matched_class is not None:
                symmetry_classes.append(matched_class)
            else:
                unique_patterns.append(pattern)
                symmetry_classes.append(len(unique_patterns) - 1)
        
        return symmetry_classes
    
    def patterns_are_similar(self, pattern1, pattern2, tolerance=0.05):
        """2つの距離パターンが類似しているかチェック"""
        symbol1, distances1 = pattern1
        symbol2, distances2 = pattern2
        
        # 元素記号が違えば類似ではない
        if symbol1 != symbol2:
            return False
        
        # 距離数が違えば類似ではない
        if len(distances1) != len(distances2):
            return False
        
        # 各距離の差をチェック
        for d1, d2 in zip(distances1, distances2):
            if abs(d1 - d2) > tolerance:
                return False
        
        return True
    
    def auto_detect_symmetry(self):
        """分子構造から最適な点群を自動検出（複数候補を表示）"""
        try:
            conf = self.mol.GetConformer()
            positions = np.array([conf.GetAtomPosition(i) for i in range(self.mol.GetNumAtoms())])
            
            # 各点群に対する適合度を計算
            candidates = self.find_all_point_group_candidates(positions)
            
            if candidates:
                # 候補選択ダイアログを表示
                selected = self.show_symmetry_candidates_dialog(candidates)
                
                if selected:
                    # 選択された点群をコンボボックスで設定
                    for i in range(self.point_group_combo.count()):
                        if self.point_group_combo.itemData(i) == selected['key']:
                            self.point_group_combo.setCurrentIndex(i)
                            break
                    
                    # 結果を表示
                    info_text = f"Selected Point Group: {selected['name']}\n"
                    info_text += f"Confidence: {selected['confidence']:.1%}\n"
                    info_text += f"Reason: {selected['reason']}\n\n"
                    info_text += "Click 'Preview Symmetry' to see detailed analysis."
                    
                    self.info_label.setText(info_text)
            else:
                self.info_label.setText("Could not automatically detect suitable point group.\nManual selection recommended.")
                
        except Exception as e:
            QMessageBox.warning(self, "Auto-Detection Error", f"Failed to auto-detect symmetry: {str(e)}")
    
    def show_symmetry_candidates_dialog(self, candidates):
        """対称性候補選択ダイアログを表示"""
        dialog = SymmetryCandidatesDialog(candidates, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            return dialog.get_selected_candidate()
        return None
    
    def find_all_point_group_candidates(self, positions):
        """分子構造に適した点群候補をすべて検索"""
        tolerance = 0.2  # Auto-detection用の比較的緩い許容値
        
        # 候補点群のリスト
        candidates = []
        
        # 分子の基本情報を取得
        num_atoms = len(positions)
        atom_symbols = [self.mol.GetAtomWithIdx(i).GetSymbol() for i in range(num_atoms)]
        
        # 1. 特殊な分子パターンを直接検出
        special_match = self.detect_special_molecules(positions, atom_symbols)
        if special_match:
            candidates.append(special_match)
        
        # 2. 各点群の対称操作に対する適合度を計算
        for key, point_group in self.POINT_GROUPS.items():
            try:
                # 既に特殊検出で追加済みの場合はスキップ
                if special_match and key == special_match['key']:
                    continue
                
                # 対称性適合度を計算
                fit_score = self.calculate_symmetry_fit(positions, point_group["operations"], tolerance)
                
                # 最低閾値をクリアした候補のみ追加
                if fit_score >= 0.3:  # 30%以上の適合度
                    reason = f"Symmetry operations fit: {fit_score:.1%}"
                    if fit_score >= 0.8:
                        reason += " (Excellent match)"
                    elif fit_score >= 0.6:
                        reason += " (Good match)"
                    elif fit_score >= 0.4:
                        reason += " (Fair match)"
                    else:
                        reason += " (Poor match)"
                    
                    candidates.append({
                        'key': key,
                        'name': point_group['name'],
                        'confidence': fit_score,
                        'reason': reason,
                        'operations_count': len(point_group["operations"])
                    })
                
            except Exception as e:
                print(f"Error evaluating {key}: {e}")
                continue
        
        # 3. スコアでソート（高い方が良い）
        try:
            candidates.sort(key=lambda x: (-x['confidence'], -x.get('operations_count', 0)))
        except KeyError as e:
            print(f"Warning: Missing key in candidate sorting: {e}")
            # フォールバック: confidenceのみでソート
            candidates.sort(key=lambda x: -x['confidence'])
        
        # 4. 上位5候補まで
        return candidates[:5]
    
    def find_best_point_group(self, positions):
        """分子構造に最も適した点群を検索"""
        tolerance = 0.2  # Auto-detection用の比較的緩い許容値
        
        # 候補点群のリスト（優先度順）
        candidates = []
        
        # 分子の基本情報を取得
        num_atoms = len(positions)
        atom_symbols = [self.mol.GetAtomWithIdx(i).GetSymbol() for i in range(num_atoms)]
        
        # 1. 特殊な分子パターンを直接検出
        special_match = self.detect_special_molecules(positions, atom_symbols)
        if special_match:
            return special_match
        
        # 2. 各点群の対称操作に対する適合度を計算
        for key, point_group in self.POINT_GROUPS.items():
            try:
                # 対称性適合度を計算
                fit_score = self.calculate_symmetry_fit(positions, point_group["operations"], tolerance)
                
                candidates.append({
                    'key': key,
                    'name': point_group['name'],
                    'score': fit_score,
                    'operations_count': len(point_group["operations"])
                })
                
            except Exception as e:
                print(f"Error evaluating {key}: {e}")
                continue
        
        # 3. 最適な候補を選択
        if not candidates:
            return None
        
        # スコアでソート（高い方が良い）
        candidates.sort(key=lambda x: (-x['score'], -x['operations_count']))
        best = candidates[0]
        
        # 最低限の閾値をチェック
        if best['score'] < 0.5:  # 50%未満の適合度は却下
            return {
                'key': 'C1',
                'name': 'C1 (No symmetry)',
                'confidence': 0.9,
                'reason': 'Low symmetry detected, suggesting C1'
            }
        
        return {
            'key': best['key'],
            'name': best['name'],
            'confidence': best['score'],
            'reason': f'Best fit among {len(candidates)} point groups tested'
        }
    
    def detect_special_molecules(self, positions, atom_symbols):
        """特殊な分子パターンを直接検出"""
        num_atoms = len(positions)
        
        # メタン様分子 (AX4)
        if self.is_methane_like():
            return {
                'key': 'Td',
                'name': 'Td (Tetrahedral)',
                'confidence': 0.95,
                'reason': 'Tetrahedral molecule detected (e.g., CH4, CCl4)',
                'operations_count': len(self.POINT_GROUPS['Td']['operations'])
            }
        
        # 水様分子 (AX2)
        if self.is_water_like():
            return {
                'key': 'C2v',
                'name': 'C2v (C2 + 2 vertical mirrors)',
                'confidence': 0.9,
                'reason': 'Bent molecule detected (e.g., H2O)',
                'operations_count': len(self.POINT_GROUPS['C2v']['operations'])
            }
        
        # アンモニア様分子 (AX3)
        if self.is_ammonia_like():
            return {
                'key': 'C3v',
                'name': 'C3v (C3 + 3 vertical mirrors)',
                'confidence': 0.9,
                'reason': 'Trigonal pyramidal molecule detected (e.g., NH3)',
                'operations_count': len(self.POINT_GROUPS['C3v']['operations'])
            }
        
        # 直線分子 (2原子または3原子直線)
        if num_atoms == 2:
            return {
                'key': 'Ci',
                'name': 'Ci (Inversion center)',
                'confidence': 0.85,
                'reason': 'Diatomic molecule detected',
                'operations_count': len(self.POINT_GROUPS['Ci']['operations'])
            }
        
        # 平面分子の検出
        if self.is_planar_molecule(positions):
            if self.has_triangular_symmetry(positions, atom_symbols):
                return {
                    'key': 'D3h',
                    'name': 'D3h (Trigonal planar)',
                    'confidence': 0.85,
                    'reason': 'Triangular planar molecule detected (e.g., BF3)',
                    'operations_count': len(self.POINT_GROUPS['D3h']['operations'])
                }
            else:
                return {
                    'key': 'Cs',
                    'name': 'Cs (Mirror plane xy)',
                    'confidence': 0.75,
                    'reason': 'Planar molecule detected',
                    'operations_count': len(self.POINT_GROUPS['Cs']['operations'])
                }
        
        return None
    
    def calculate_symmetry_fit(self, positions, operations, tolerance):
        """対称操作に対する構造の適合度を計算（0-1のスコア）"""
        if len(operations) <= 1:  # 恒等操作のみ
            return 0.1
        
        centroid = np.mean(positions, axis=0)
        centered_positions = positions - centroid
        
        total_score = 0.0
        valid_operations = 0
        
        for operation in operations[1:]:  # 恒等操作を除く
            operation_score = 0.0
            
            for pos in centered_positions:
                transformed_pos = operation @ pos
                
                # 変換された位置に最も近い実際の原子を探す
                distances = [np.linalg.norm(transformed_pos - other_pos) 
                           for other_pos in centered_positions]
                min_distance = min(distances)
                
                # 距離に基づくスコア（近いほど高スコア）
                if min_distance < tolerance:
                    operation_score += 1.0 - (min_distance / tolerance)
                
            # 正規化（原子数で割る）
            operation_score /= len(positions)
            total_score += operation_score
            valid_operations += 1
        
        # 全対称操作での平均スコア
        if valid_operations > 0:
            return total_score / valid_operations
        else:
            return 0.0
    
    def is_planar_molecule(self, positions, tolerance=0.1):
        """分子が平面構造かどうか判定"""
        if len(positions) < 4:
            return True  # 3原子以下は常に平面
        
        # 主成分分析で平面性をチェック
        centroid = np.mean(positions, axis=0)
        centered = positions - centroid
        
        # 共分散行列の固有値を計算
        cov_matrix = np.cov(centered.T)
        eigenvalues = np.linalg.eigvals(cov_matrix)
        eigenvalues = np.sort(eigenvalues)
        
        # 最小固有値が小さければ平面的
        return eigenvalues[0] < tolerance
    
    def has_triangular_symmetry(self, positions, atom_symbols):
        """三角対称性を持つかどうか判定"""
        if len(positions) < 3:
            return False
        
        # 中心原子を特定
        center_candidates = []
        for i, atom in enumerate(self.mol.GetAtoms()):
            if atom.GetDegree() >= 3:
                center_candidates.append(i)
        
        if len(center_candidates) != 1:
            return False
        
        center_idx = center_candidates[0]
        center_atom = self.mol.GetAtomWithIdx(center_idx)
        
        # 隣接原子が3個で、すべて同じ元素かチェック
        neighbors = list(center_atom.GetNeighbors())
        if len(neighbors) != 3:
            return False
        
        neighbor_symbols = [atom.GetSymbol() for atom in neighbors]
        return len(set(neighbor_symbols)) == 1
    
    def group_equivalent_atoms(self, symmetry_classes):
        """対称性クラスから等価原子グループを作成"""
        from collections import defaultdict
        
        groups = defaultdict(list)
        for atom_idx, sym_class in enumerate(symmetry_classes):
            groups[sym_class].append(atom_idx)
        
        print(f"Debug: Raw groups from symmetry classes: {dict(groups)}")
        
        # 等価原子グループ（2個以上）と単独原子を分離
        equivalent_groups = []
        single_atoms = []
        
        for sym_class, atom_indices in groups.items():
            if len(atom_indices) > 1:
                equivalent_groups.append(atom_indices)
                print(f"Debug: Equivalent group found - Class {sym_class}: {atom_indices}")
            else:
                single_atoms.extend([[idx] for idx in atom_indices])
        
        # 結果を組み合わせ
        all_groups = equivalent_groups + single_atoms
        
        print(f"Debug: Final atom groups: {all_groups}")
        print(f"Debug: Number of equivalent groups (>1 atom): {len(equivalent_groups)}")
        
        return all_groups
    
    def calculate_ideal_symmetric_positions(self, group_positions, operations, tolerance):
        """等価原子グループの理想的な対称位置を計算"""
        if len(group_positions) == 1:
            # 単独原子の場合、全対称操作の平均位置を計算
            pos = group_positions[0]
            equivalent_positions = [op @ pos for op in operations]
            return [np.mean(equivalent_positions, axis=0)]
        
        # 複数原子グループの場合、グループ全体の対称性を考慮
        group_centroid = np.mean(group_positions, axis=0)
        ideal_positions = []
        
        for pos in group_positions:
            # 各原子に対して対称操作を適用し、理想位置を計算
            relative_pos = pos - group_centroid
            equivalent_relatives = [op @ relative_pos for op in operations]
            ideal_relative = np.mean(equivalent_relatives, axis=0)
            ideal_positions.append(group_centroid + ideal_relative)
        
        return ideal_positions
    
    def analyze_symmetry_simple(self, positions, operations, tolerance):
        """従来の単純な対称性分析（フォールバック）"""
        atoms_to_move = 0
        max_displacement = 0.0
        
        centroid = np.mean(positions, axis=0)
        centered_positions = positions - centroid
        
        for i, pos in enumerate(centered_positions):
            # 各対称操作を適用した等価位置を計算
            equivalent_positions = []
            for operation in operations:
                equivalent_positions.append(operation @ pos)
            
            # 等価位置の重心を計算（理想的な対称位置）
            ideal_pos = np.mean(equivalent_positions, axis=0)
            
            # 現在位置と理想位置の距離を計算
            displacement = np.linalg.norm(ideal_pos - pos)
            
            if displacement > tolerance:
                atoms_to_move += 1
                max_displacement = max(max_displacement, displacement)
        
        return {
            'atoms_to_move': atoms_to_move,
            'max_displacement': max_displacement
        }
    
    def apply_symmetry_operations(self, positions, operations, tolerance):
        """対称操作を適用して構造を対称化（安全なバージョン）"""
        try:
            # シンプルで直接的な方法を使用
            return self.apply_symmetry_direct(positions, operations, tolerance)
        except Exception as e:
            print(f"Symmetrization failed: {e}")
            return positions.copy()
    
    def apply_symmetry_direct(self, positions, operations, tolerance):
        """直接的で安全な対称化（分子タイプを考慮）"""
        # 対称操作の妥当性をチェック
        valid_operations = []
        for op in operations:
            try:
                if op.shape == (3, 3) and not np.any(np.isnan(op)) and not np.any(np.isinf(op)):
                    det = np.linalg.det(op)
                    if abs(abs(det) - 1.0) < 0.1:
                        valid_operations.append(op)
            except:
                continue
        
        print(f"Debug: Total operations: {len(operations)}, Valid operations: {len(valid_operations)}")
        
        if not valid_operations:
            print("警告: 有効な対称操作がありません。元の座標を返します。")
            return positions.copy()
        
        # 分子の中心を計算
        centroid = np.mean(positions, axis=0)
        print(f"Debug: Molecular centroid: {centroid}")
        
        # 等価原子グループを取得
        symmetry_classes = self.get_molecular_symmetry_classes()
        equivalent_atom_groups = self.group_equivalent_atoms(symmetry_classes)
        
        new_positions = positions.copy()
        
        # メタン分子の特別処理
        if len(positions) == 5 and len(equivalent_atom_groups) == 2:
            print("Debug: Detected methane-like molecule")
            
            # 中心原子（通常は炭素）を原点に移動
            center_atom_group = [group for group in equivalent_atom_groups if len(group) == 1][0]
            center_atom_idx = center_atom_group[0]
            
            # 周辺原子グループ（通常は水素）
            peripheral_atoms = [group for group in equivalent_atom_groups if len(group) > 1][0]
            
            print(f"Debug: Center atom: {center_atom_idx}, Peripheral atoms: {peripheral_atoms}")
            
            # 中心原子を分子の重心に配置（保守的なアプローチ）
            current_center = positions[center_atom_idx]
            if np.linalg.norm(current_center - centroid) > tolerance:
                new_positions[center_atom_idx] = centroid
                print(f"Debug: Moved center atom to centroid")
            
            # 周辺原子を中心からの相対位置で対称化（距離を保持）
            for atom_idx in peripheral_atoms:
                current_pos = positions[atom_idx]
                relative_pos = current_pos - centroid
                original_distance = np.linalg.norm(relative_pos)
                
                # 対称操作を適用（中心からの相対座標で）
                equivalent_relative_positions = []
                for op in valid_operations:
                    try:
                        transformed_relative = op @ relative_pos
                        equivalent_relative_positions.append(transformed_relative)
                    except:
                        continue
                
                if len(equivalent_relative_positions) > 1:
                    # 等価相対位置の平均
                    average_relative = np.mean(equivalent_relative_positions, axis=0)
                    
                    # 元の距離を保持
                    if np.linalg.norm(average_relative) > 1e-10:
                        average_relative = average_relative / np.linalg.norm(average_relative) * original_distance
                    
                    average_pos = centroid + average_relative
                    displacement = np.linalg.norm(average_pos - current_pos)
                    
                    print(f"Debug: Atom {atom_idx} - Original distance: {original_distance:.6f}, Displacement: {displacement:.6f}")
                    
                    # より保守的な条件：大きな変位は避ける
                    if displacement > tolerance and displacement < 2.0:  # 最大2Å未満の変位のみ許可
                        new_positions[atom_idx] = average_pos
                        print(f"Debug: Applied symmetrization to atom {atom_idx}")
                    elif displacement >= 2.0:
                        print(f"Debug: Skipped atom {atom_idx} - displacement too large: {displacement:.6f}")
        else:
            # 一般的な分子の処理
            for i, pos in enumerate(positions):
                # すべての対称操作を適用して等価位置を取得
                equivalent_positions = []
                for op in valid_operations:
                    try:
                        transformed_pos = op @ pos
                        equivalent_positions.append(transformed_pos)
                    except:
                        continue
                
                if len(equivalent_positions) > 1:
                    # 等価位置の重心を新しい位置とする
                    average_pos = np.mean(equivalent_positions, axis=0)
                    displacement = np.linalg.norm(average_pos - pos)
                    
                    print(f"Debug: Atom {i} - Displacement: {displacement}")
                    
                    # 保守的な条件：大きな変位は避ける
                    if displacement > tolerance and displacement < 2.0:  # 最大2Å未満の変位のみ許可
                        new_positions[i] = average_pos
                        print(f"Debug: Applied symmetrization to atom {i}")
                    elif displacement >= 2.0:
                        print(f"Debug: Skipped atom {i} - displacement too large: {displacement:.6f}")
        
        return new_positions
    
    def apply_symmetry_advanced(self, positions, operations, tolerance):
        """等価原子グループを考慮した高度な対称化"""
        # 対称操作の妥当性をチェック
        valid_operations = []
        for op in operations:
            try:
                # 3x3の行列であることを確認
                if op.shape == (3, 3) and not np.any(np.isnan(op)) and not np.any(np.isinf(op)):
                    # 行列式が±1に近いことを確認（回転・反射行列の条件）
                    det = np.linalg.det(op)
                    if abs(abs(det) - 1.0) < 0.1:  # より緩い条件
                        valid_operations.append(op)
            except:
                continue
        
        print(f"Debug: Total operations: {len(operations)}, Valid operations: {len(valid_operations)}")
        
        if not valid_operations:
            print("警告: 有効な対称操作がありません。元の座標を返します。")
            return positions.copy()
        
        # 分子の対称性クラスを取得
        symmetry_classes = self.get_molecular_symmetry_classes()
        equivalent_atom_groups = self.group_equivalent_atoms(symmetry_classes)
        
        centroid = np.mean(positions, axis=0)
        centered_positions = positions - centroid
        new_positions = centered_positions.copy()
        
        print(f"Debug: Original centroid: {centroid}")
        print(f"Debug: Original positions shape: {positions.shape}")
        print(f"Debug: Valid operations count: {len(valid_operations)}")
        
        # 各等価原子グループを個別に対称化
        max_iterations = 5  # 反復回数を制限
        for iteration in range(max_iterations):
            print(f"Debug: Iteration {iteration + 1}")
            moved_any = False
            
            for group_atoms in equivalent_atom_groups:
                if len(group_atoms) == 1:
                    # 単独原子の場合
                    atom_idx = group_atoms[0]
                    pos = new_positions[atom_idx]
                    
                    print(f"Debug: Processing single atom {atom_idx}")
                    
                    # 全ての対称操作で得られる等価位置の平均
                    equivalent_positions = []
                    for op in valid_operations:
                        try:
                            transformed_pos = op @ pos
                            equivalent_positions.append(transformed_pos)
                        except:
                            continue
                    
                        if equivalent_positions:
                            average_pos = np.mean(equivalent_positions, axis=0)
                            displacement = np.linalg.norm(average_pos - pos)
                            
                            print(f"Debug: Atom {atom_idx} - Displacement: {displacement}")
                            
                            # 数値精度を考慮したより緩い条件
                            if displacement > max(tolerance / 100, 1e-10):
                                new_positions[atom_idx] = average_pos
                                moved_any = True
                        
                else:
                    # 等価原子グループの場合、グループ全体の対称性を保持
                    group_positions = new_positions[group_atoms]
                    group_centroid = np.mean(group_positions, axis=0)
                    
                    print(f"Debug: Processing group {group_atoms}")
                    
                    # グループの重心を対称化
                    equivalent_centroids = []
                    for op in valid_operations:
                        try:
                            transformed_centroid = op @ group_centroid
                            equivalent_centroids.append(transformed_centroid)
                        except:
                            continue
                    
                    if equivalent_centroids:
                        ideal_centroid = np.mean(equivalent_centroids, axis=0)
                        centroid_shift = ideal_centroid - group_centroid
                        
                        print(f"Debug: Group centroid shift magnitude: {np.linalg.norm(centroid_shift)}")
                        
                        # グループ内の相対位置を保持しながら重心を移動
                        if np.linalg.norm(centroid_shift) > max(tolerance / 100, 1e-10):
                            for atom_idx in group_atoms:
                                new_positions[atom_idx] += centroid_shift
                            moved_any = True
                    
                    # グループ内の原子間の相対位置も対称化
                    for atom_idx in group_atoms:
                        pos = new_positions[atom_idx]
                        
                        # 全ての対称操作で得られる等価位置の平均
                        equivalent_positions = []
                        for op in valid_operations:
                            try:
                                transformed_pos = op @ pos
                                equivalent_positions.append(transformed_pos)
                            except:
                                continue
                        
                            if equivalent_positions:
                                average_pos = np.mean(equivalent_positions, axis=0)
                                displacement = np.linalg.norm(average_pos - pos)
                                
                                if displacement > max(tolerance / 100, 1e-10):
                                    new_positions[atom_idx] = average_pos
                                    moved_any = True
            
            if not moved_any:
                break
        
        # 重心を元に戻す
        final_positions = new_positions + centroid
        
        print(f"Debug: Final positions before return:")
        for i, pos in enumerate(final_positions):
            print(f"  Atom {i}: {pos}")
        
        # デバッグ: 座標チェック
        if np.any(np.isnan(final_positions)) or np.any(np.isinf(final_positions)):
            print("警告: 無効な座標が検出されました。元の座標を返します。")
            return positions.copy()
        
        return final_positions
    
    def apply_symmetry_simple(self, positions, operations, tolerance):
        """従来の単純な対称化（フォールバック）"""
        # 対称操作の妥当性をチェック
        valid_operations = []
        for op in operations:
            try:
                # 3x3の行列であることを確認
                if op.shape == (3, 3) and not np.any(np.isnan(op)) and not np.any(np.isinf(op)):
                    # 行列式が±1に近いことを確認（回転・反射行列の条件）
                    det = np.linalg.det(op)
                    if abs(abs(det) - 1.0) < 0.1:  # より緩い条件
                        valid_operations.append(op)
            except:
                continue
        
        if not valid_operations:
            print("警告: 有効な対称操作がありません。元の座標を返します。")
            return positions.copy()
        
        centroid = np.mean(positions, axis=0)
        centered_positions = positions - centroid
        new_positions = centered_positions.copy()
        
        # 反復的に対称化を適用
        max_iterations = 5  # 反復回数を制限
        for iteration in range(max_iterations):
            moved_any = False
            
            for i, pos in enumerate(new_positions):
                # 全ての対称操作で得られる等価位置の平均を計算
                equivalent_positions = []
                for operation in valid_operations:
                    try:
                        equivalent_positions.append(operation @ pos)
                    except:
                        continue
                
                if equivalent_positions:
                    # 等価位置の重心を新しい位置とする
                    average_pos = np.mean(equivalent_positions, axis=0)
                    displacement = np.linalg.norm(average_pos - pos)
                    
                    if displacement > max(tolerance / 100, 1e-10):  # 数値精度を考慮
                        new_positions[i] = average_pos
                        moved_any = True
            
            if not moved_any:
                break
        
        # 重心を元に戻す
        final_positions = new_positions + centroid
        
        # デバッグ: 座標チェック
        if np.any(np.isnan(final_positions)) or np.any(np.isinf(final_positions)):
            print("警告: 無効な座標が検出されました。元の座標を返します。")
            return positions.copy()
        
        return final_positions

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
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        button_layout.addWidget(apply_button)
        button_layout.addWidget(cancel_button)
        
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
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply mirror transformation: {str(e)}")

class PlanarizationDialog(Dialog3DPickingMixin, QDialog):
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
        self.setWindowTitle(f"Planarize to {plane_names[self.plane]} Plane")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(f"Click atoms in the 3D view to select them for planarization to the {plane_names[self.plane]} plane. At least 3 atoms are required.")
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
        
        self.apply_button = QPushButton("Apply Planarization")
        self.apply_button.clicked.connect(self.apply_planarization)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
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
                self.apply_planarization()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
    
    def update_display(self):
        """表示を更新"""
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select for planarization (minimum 3 required)")
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
    
    def apply_planarization(self):
        """平面化を適用（回転ベース）"""
        if len(self.selected_atoms) < 3:
            QMessageBox.warning(self, "Warning", "Please select at least 3 atoms for planarization.")
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
            
            self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply planarization: {str(e)}")
    
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
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
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
            
            QMessageBox.information(self, "Success", f"Alignment to {self.axis.upper()}-axis completed.")
            self.accept()
            
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
                import traceback
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
            import traceback
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
                    ay = pos.y() * ANGSTROM_PER_PIXEL
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
            total_dx = sum((b.atom2.pos().x() if b.atom1 is self else b.atom1.pos().x()) - my_pos_x for b in self.bonds)
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
            if self.charge == 1: charge_str = "+"
            elif self.charge == -1: charge_str = "-"
            else: charge_str = f"{self.charge:+}"
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
                total_dx = 0
                for bond in self.bonds:
                    other_atom = bond.atom1 if bond.atom2 is self else bond.atom2
                    total_dx += (other_atom.pos().x() - my_pos_x)

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
                if self.charge == 1: charge_str = "+"
                elif self.charge == -1: charge_str = "-"
                else: charge_str = f"{self.charge:+}"
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
                from PyQt6.QtWidgets import QGraphicsScene
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

    def set_geometry(self, points, is_aromatic=False):
        self.prepareGeometryChange()
        self.polygon = QPolygonF(points)
        self.is_aromatic = is_aromatic
        self.update()

    def boundingRect(self):
        return self.polygon.boundingRect().adjusted(-5, -5, 5, 5)

    def paint(self, painter, option, widget):
            
        painter.setPen(self.pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        if not self.polygon.isEmpty():
            painter.drawPolygon(self.polygon)
            if self.is_aromatic:
                center = self.polygon.boundingRect().center()
                radius = QLineF(center, self.polygon.first()).length() * 0.6
                painter.drawEllipse(center, radius, radius)

class MoleculeScene(QGraphicsScene):
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
        self.clear_all_problem_flags()
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
                        import traceback
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
            self.removeItem(self.temp_line)
            self.temp_line = None

        if self.mode.startswith('template') and is_click:
            if self.template_context and self.template_context.get('points'):
                context = self.template_context
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
                    import traceback
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
            import traceback
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
        """指定されたアイテムセット（原子・結合）を安全に削除する共通メソッド"""
        try:
            if not items_to_delete:
                return False

            atoms_to_delete = {item for item in items_to_delete if isinstance(item, AtomItem)}
            bonds_to_delete = {item for item in items_to_delete if isinstance(item, BondItem)}

            # 削除対象の原子に接続している結合も、すべて削除対象に加える
            for atom in atoms_to_delete:
                if hasattr(atom, 'bonds'):
                    bonds_to_delete.update(atom.bonds)

            # 影響を受ける（が削除はされない）原子を特定する
            atoms_to_update = set()
            for bond in bonds_to_delete:
                if bond.atom1 and bond.atom1 not in atoms_to_delete:
                    atoms_to_update.add(bond.atom1)
                if bond.atom2 and bond.atom2 not in atoms_to_delete:
                    atoms_to_update.add(bond.atom2)

            # --- データモデルからの削除 ---
            # 最初に原子をデータモデルから削除（関連する結合も内部で削除される）
            for atom in atoms_to_delete:
                if hasattr(atom, 'atom_id'):
                    self.data.remove_atom(atom.atom_id)
                    
            # 次に、明示的に選択された結合（まだ残っているもの）を削除
            for bond in bonds_to_delete:
                if bond.atom1 and bond.atom2 and hasattr(bond.atom1, 'atom_id') and hasattr(bond.atom2, 'atom_id'):
                    self.data.remove_bond(bond.atom1.atom_id, bond.atom2.atom_id)
            
            # --- シーンからのグラフィックアイテム削除（必ず結合を先に）---
            # まずremoveItem
            for bond in bonds_to_delete:
                if bond.scene(): 
                    self.removeItem(bond)
                    
            for atom in atoms_to_delete:
                if atom.scene(): 
                    self.removeItem(atom)
            
            # その後BondItem/AtomItemの参照をクリア
            for bond in bonds_to_delete:
                # atom1/atom2のbondsリストからこのbondを除去
                try:
                    if bond.atom1 and hasattr(bond.atom1, 'bonds') and bond in bond.atom1.bonds:
                        bond.atom1.bonds = [b for b in bond.atom1.bonds if b is not bond]
                    if bond.atom2 and hasattr(bond.atom2, 'bonds') and bond in bond.atom2.bonds:
                        bond.atom2.bonds = [b for b in bond.atom2.bonds if b is not bond]
                    # Clear bond references to prevent dangling pointers
                    bond.atom1 = None
                    bond.atom2 = None
                except Exception as bond_cleanup_error:
                    print(f"Error cleaning up bond references: {bond_cleanup_error}")
            
            for atom in atoms_to_delete:
                # Clear all bond references from deleted atoms
                try:
                    if hasattr(atom, 'bonds'):
                        atom.bonds.clear()
                except Exception as atom_cleanup_error:
                    print(f"Error cleaning up atom references: {atom_cleanup_error}")

            # --- 生き残った原子の内部参照とスタイルを更新 ---
            for atom in atoms_to_update:
                try:
                    if hasattr(atom, 'bonds'):
                        atom.bonds = [b for b in atom.bonds if b not in bonds_to_delete]
                    if hasattr(atom, 'update_style'):
                        atom.update_style()
                except Exception as update_error:
                    print(f"Error updating atom style: {update_error}")
                    
            return True
            
        except Exception as e:
            print(f"Error during delete_items operation: {e}")
            import traceback
            traceback.print_exc()
            return False

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
                self.removeItem(self.temp_line)
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
            import traceback
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
    finished=pyqtSignal(object); error=pyqtSignal(str)
    
    def run_calculation(self, mol_block):
        try:
            if not mol_block:
                raise ValueError("No atoms to convert.")
            
            self.status_update.emit("Creating 3D structure...")
            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None:
                raise ValueError("Failed to create molecule from MOL block.")

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
            
            self.status_update.emit("Embedding 3D coordinates...")
            
            # Try multiple times with different approaches if needed
            conf_id = -1
            
            # First attempt: Standard ETKDG with stereo enforcement
            try:
                conf_id = AllChem.EmbedMolecule(mol, params)
            except Exception as e:
                self.status_update.emit(f"Standard embedding failed: {e}")
            
            # Second attempt: Use constraint embedding if available
            if conf_id == -1:
                try:
                    # Create distance constraints for double bonds to enforce E/Z geometry
                    from rdkit.DistanceGeometry import DoTriangleSmoothing
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
                    self.status_update.emit("Constraint-based embedding succeeded")
                except Exception as e:
                    self.status_update.emit(f"Constraint embedding failed: {e}")
                    
            # Fallback: Try basic embedding
            if conf_id == -1:
                try:
                    basic_params = AllChem.ETKDGv2()
                    basic_params.randomSeed = 42
                    conf_id = AllChem.EmbedMolecule(mol, basic_params)
                except Exception:
                    pass
            '''
            if conf_id == -1:
                self.status_update.emit("Initial embedding failed, retrying with ignoreSmoothingFailures=True...")
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

            if conf_id != -1:
                # Success with RDKit: optimize and finish
                # CRITICAL: Restore original stereochemistry after embedding (explicit labels first)
                for bond_idx, stereo, stereo_atoms in original_stereo_info:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if len(stereo_atoms) == 2:
                        bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                    bond.SetStereo(stereo)
                
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except Exception:
                    # fallback to UFF if MMFF fails
                    try:
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
                self.finished.emit(mol)
                self.status_update.emit("RDKit 3D conversion succeeded.")
                return

            # ---------- RDKit failed: try Open Babel via pybel only (no CLI fallback) ----------
            self.status_update.emit("RDKit embedding failed. Attempting Open Babel fallback...")

            try:
                # pybel expects an input format; provide mol block
                # pybel.readstring accepts format strings like "mol" or "smi"
                ob_mol = pybel.readstring("mol", mol_block)
                # ensure hydrogens
                try:
                    ob_mol.addh()
                except Exception:
                    pass
                # build 3D coordinates
                ob_mol.make3D()
                try:
                    # まず第一候補であるMMFF94で最適化を試みる
                    self.status_update.emit("Optimizing with Open Babel (MMFF94)...")
                    ob_mol.localopt(forcefield='mmff94', steps=500)
                except Exception:
                    # MMFF94が失敗した場合、UFFにフォールバックして再試行
                    try:
                        self.status_update.emit("MMFF94 failed, falling back to UFF...")
                        ob_mol.localopt(forcefield='uff', steps=500)
                    except Exception:
                        # UFFも失敗した場合はスキップ
                        self.status_update.emit("UFF optimization also failed.")
                        pass
                # get molblock and convert to RDKit
                molblock_ob = ob_mol.write("mol")
                rd_mol = Chem.MolFromMolBlock(molblock_ob, removeHs=False)
                if rd_mol is None:
                    raise ValueError("Open Babel produced invalid MOL block.")
                # optimize in RDKit as a final step if possible
                rd_mol = Chem.AddHs(rd_mol)
                # Do NOT call AssignStereochemistry here to preserve original stereo info
                try:
                    AllChem.MMFFOptimizeMolecule(rd_mol)
                except Exception:
                    try:
                        AllChem.UFFOptimizeMolecule(rd_mol)
                    except Exception:
                        pass
                # Do NOT call AssignStereochemistry after optimization either
                self.status_update.emit("Open Babel embedding succeeded. Warning: Conformation accuracy may be limited.")
                self.finished.emit(rd_mol)
                return
            except Exception as ob_err:
                # pybel was available but failed
                raise RuntimeError(f"Open Babel 3D conversion failed: {ob_err}")

        except Exception as e:
            self.error.emit(str(e))


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
            ('Cs',6,1), ('Ba',6,2), ('La',6,3), ('Hf',6,4), ('Ta',6,5), ('W',6,6), ('Re',6,7), ('Os',6,8),
            ('Ir',6,9), ('Pt',6,10), ('Au',6,11), ('Hg',6,12), ('Tl',6,13), ('Pb',6,14), ('Bi',6,15), ('Po',6,16),
            ('At',6,17), ('Rn',6,18),
            ('Fr',7,1), ('Ra',7,2), ('Ac',7,3), ('Rf',7,4), ('Db',7,5), ('Sg',7,6), ('Bh',7,7), ('Hs',7,8),
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

            # CPK_COLORSから色を取得。見つからない場合はデフォルト色を使用
            q_color = CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT'])

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
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors

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
                from rdkit import Chem
                
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
                except:
                    inchi = "N/A"
                
                # 表示するプロパティを辞書にまとめる
                properties = {
                    "SMILES:": smiles,
                    "InChI:": inchi,
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
        self.default_settings = {
            'background_color': '#919191',
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
            # Multiple bond offset parameters
            'double_bond_offset_factor': 2.0,
            'triple_bond_offset_factor': 2.0,
            'double_bond_radius_factor': 0.8,
            'triple_bond_radius_factor': 0.7,
        }
        
        # --- 選択された色を管理する専用のインスタンス変数 ---
        self.current_bg_color = None

        # --- UI要素の作成 ---
        layout = QVBoxLayout(self)
        
        # タブウィジェットを作成
        self.tab_widget = QTabWidget()
        layout.addWidget(self.tab_widget)
        
        # 基本設定タブ
        self.create_general_tab()
        
        # Common設定タブ
        self.create_common_tab()
        
        # Ball and Stick設定タブ
        self.create_ball_stick_tab()
        
        # CPK設定タブ
        self.create_cpk_tab()
        
        # Wireframe設定タブ
        self.create_wireframe_tab()
        
        # Stick設定タブ
        self.create_stick_tab()

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
    
    def create_general_tab(self):
        """基本設定タブを作成"""
        general_widget = QWidget()
        form_layout = QFormLayout(general_widget)

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
        
        self.tab_widget.addTab(general_widget, "Scene")
    
    def create_common_tab(self):
        """Common設定タブを作成"""
        common_widget = QWidget()
        form_layout = QFormLayout(common_widget)
        
        # 二重結合オフセット倍率
        self.double_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.double_offset_slider.setRange(100, 400)  # 1.0 ~ 4.0
        self.double_offset_label = QLabel("1.50")
        self.double_offset_slider.valueChanged.connect(lambda v: self.double_offset_label.setText(f"{v/100:.2f}"))
        double_offset_layout = QHBoxLayout()
        double_offset_layout.addWidget(self.double_offset_slider)
        double_offset_layout.addWidget(self.double_offset_label)
        form_layout.addRow("Double Bond Offset:", double_offset_layout)
        
        # 三重結合オフセット倍率
        self.triple_offset_slider = QSlider(Qt.Orientation.Horizontal)
        self.triple_offset_slider.setRange(100, 400)  # 1.0 ~ 4.0
        self.triple_offset_label = QLabel("1.20")
        self.triple_offset_slider.valueChanged.connect(lambda v: self.triple_offset_label.setText(f"{v/100:.2f}"))
        triple_offset_layout = QHBoxLayout()
        triple_offset_layout.addWidget(self.triple_offset_slider)
        triple_offset_layout.addWidget(self.triple_offset_label)
        form_layout.addRow("Triple Bond Offset:", triple_offset_layout)
        
        # 二重結合半径倍率
        self.double_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.double_radius_slider.setRange(50, 100)  # 0.5 ~ 1.0
        self.double_radius_label = QLabel("0.80")
        self.double_radius_slider.valueChanged.connect(lambda v: self.double_radius_label.setText(f"{v/100:.2f}"))
        double_radius_layout = QHBoxLayout()
        double_radius_layout.addWidget(self.double_radius_slider)
        double_radius_layout.addWidget(self.double_radius_label)
        form_layout.addRow("Double Bond Thickness:", double_radius_layout)
        
        # 三重結合半径倍率
        self.triple_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.triple_radius_slider.setRange(50, 100)  # 0.5 ~ 1.0
        self.triple_radius_label = QLabel("0.70")
        self.triple_radius_slider.valueChanged.connect(lambda v: self.triple_radius_label.setText(f"{v/100:.2f}"))
        triple_radius_layout = QHBoxLayout()
        triple_radius_layout.addWidget(self.triple_radius_slider)
        triple_radius_layout.addWidget(self.triple_radius_label)
        form_layout.addRow("Triple Bond Thickness:", triple_radius_layout)
        
        self.tab_widget.addTab(common_widget, "Common")
    
    def create_ball_stick_tab(self):
        """Ball and Stick設定タブを作成"""
        ball_stick_widget = QWidget()
        form_layout = QFormLayout(ball_stick_widget)
        
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
        
        # 解像度
        self.bs_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.bs_resolution_slider.setRange(6, 32)
        self.bs_resolution_label = QLabel("16")
        self.bs_resolution_slider.valueChanged.connect(lambda v: self.bs_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.bs_resolution_slider)
        resolution_layout.addWidget(self.bs_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)
        
        self.tab_widget.addTab(ball_stick_widget, "Ball & Stick")
    
    def create_cpk_tab(self):
        """CPK設定タブを作成"""
        cpk_widget = QWidget()
        form_layout = QFormLayout(cpk_widget)
        
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
        
        info_label = QLabel("CPK model shows atoms as space-filling spheres using van der Waals radii.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow("", info_label)
        
        self.tab_widget.addTab(cpk_widget, "CPK (Space-filling)")
    
    def create_wireframe_tab(self):
        """Wireframe設定タブを作成"""
        wireframe_widget = QWidget()
        form_layout = QFormLayout(wireframe_widget)
        
        # ボンド半径
        self.wf_bond_radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_bond_radius_slider.setRange(1, 10)  # 0.01 ~ 0.1
        self.wf_bond_radius_label = QLabel("0.01")
        self.wf_bond_radius_slider.valueChanged.connect(lambda v: self.wf_bond_radius_label.setText(f"{v/100:.2f}"))
        bond_radius_layout = QHBoxLayout()
        bond_radius_layout.addWidget(self.wf_bond_radius_slider)
        bond_radius_layout.addWidget(self.wf_bond_radius_label)
        form_layout.addRow("Bond Radius:", bond_radius_layout)
        
        # 解像度
        self.wf_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.wf_resolution_slider.setRange(4, 16)
        self.wf_resolution_label = QLabel("6")
        self.wf_resolution_slider.valueChanged.connect(lambda v: self.wf_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.wf_resolution_slider)
        resolution_layout.addWidget(self.wf_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)
        
        info_label = QLabel("Wireframe model shows molecular structure with thin lines only (no atoms displayed).")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow("", info_label)
        
        self.tab_widget.addTab(wireframe_widget, "Wireframe")
    
    def create_stick_tab(self):
        """Stick設定タブを作成"""
        stick_widget = QWidget()
        form_layout = QFormLayout(stick_widget)
        
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
        
        # 解像度
        self.stick_resolution_slider = QSlider(Qt.Orientation.Horizontal)
        self.stick_resolution_slider.setRange(6, 32)
        self.stick_resolution_label = QLabel("16")
        self.stick_resolution_slider.valueChanged.connect(lambda v: self.stick_resolution_label.setText(str(v)))
        resolution_layout = QHBoxLayout()
        resolution_layout.addWidget(self.stick_resolution_slider)
        resolution_layout.addWidget(self.stick_resolution_label)
        form_layout.addRow("Resolution (Quality):", resolution_layout)
        
        info_label = QLabel("Stick model shows bonds as thick cylinders with atoms as small spheres.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic; margin-top: 10px;")
        form_layout.addRow("", info_label)
        
        self.tab_widget.addTab(stick_widget, "Stick")

    def reset_current_tab(self):
        """現在選択されているタブの設定のみをデフォルトに戻す"""
        current_tab_index = self.tab_widget.currentIndex()
        tab_name = self.tab_widget.tabText(current_tab_index)
        
        # 各タブの設定項目を定義
        tab_settings = {
            "Scene": {
                'background_color': self.default_settings['background_color'],
                'show_3d_axes': self.default_settings['show_3d_axes'],
                'lighting_enabled': self.default_settings['lighting_enabled'],
                'light_intensity': self.default_settings['light_intensity'],
                'specular': self.default_settings['specular'],
                'specular_power': self.default_settings['specular_power']
            },
            "Common": {
                'double_bond_offset_factor': self.default_settings['double_bond_offset_factor'],
                'triple_bond_offset_factor': self.default_settings['triple_bond_offset_factor'],
                'double_bond_radius_factor': self.default_settings['double_bond_radius_factor'],
                'triple_bond_radius_factor': self.default_settings['triple_bond_radius_factor']
            },
            "Ball & Stick": {
                'ball_stick_atom_scale': self.default_settings['ball_stick_atom_scale'],
                'ball_stick_bond_radius': self.default_settings['ball_stick_bond_radius'],
                'ball_stick_resolution': self.default_settings['ball_stick_resolution']
            },
            "CPK (Space-filling)": {
                'cpk_atom_scale': self.default_settings['cpk_atom_scale'],
                'cpk_resolution': self.default_settings['cpk_resolution']
            },
            "Wireframe": {
                'wireframe_bond_radius': self.default_settings['wireframe_bond_radius'],
                'wireframe_resolution': self.default_settings['wireframe_resolution']
            },
            "Stick": {
                'stick_atom_radius': self.default_settings['stick_atom_radius'],
                'stick_bond_radius': self.default_settings['stick_bond_radius'],
                'stick_resolution': self.default_settings['stick_resolution']
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
            self.update_ui_from_settings(self.default_settings)
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
            # Multi-bond settings
            'double_bond_offset_factor': self.double_offset_slider.value() / 100.0,
            'triple_bond_offset_factor': self.triple_offset_slider.value() / 100.0,
            'double_bond_radius_factor': self.double_radius_slider.value() / 100.0,
            'triple_bond_radius_factor': self.triple_radius_slider.value() / 100.0,
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
        
        # 多重結合設定
        double_offset = int(settings_dict.get('double_bond_offset_factor', self.default_settings['double_bond_offset_factor']) * 100)
        self.double_offset_slider.setValue(double_offset)
        self.double_offset_label.setText(f"{double_offset/100:.2f}")
        
        triple_offset = int(settings_dict.get('triple_bond_offset_factor', self.default_settings['triple_bond_offset_factor']) * 100)
        self.triple_offset_slider.setValue(triple_offset)
        self.triple_offset_label.setText(f"{triple_offset/100:.2f}")
        
        double_radius = int(settings_dict.get('double_bond_radius_factor', self.default_settings['double_bond_radius_factor']) * 100)
        self.double_radius_slider.setValue(double_radius)
        self.double_radius_label.setText(f"{double_radius/100:.2f}")
        
        triple_radius = int(settings_dict.get('triple_bond_radius_factor', self.default_settings['triple_bond_radius_factor']) * 100)
        self.triple_radius_slider.setValue(triple_radius)
        self.triple_radius_label.setText(f"{triple_radius/100:.2f}")
      
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
            # Multiple bond offset settings
            'double_bond_offset_factor': self.double_offset_slider.value() / 100.0,
            'triple_bond_offset_factor': self.triple_offset_slider.value() / 100.0,
            'double_bond_radius_factor': self.double_radius_slider.value() / 100.0,
            'triple_bond_radius_factor': self.triple_radius_slider.value() / 100.0,
        }

    def apply_settings(self):
        """設定を適用（ダイアログは開いたまま）"""
        # 親ウィンドウの設定を更新
        if self.parent_window:
            settings = self.get_settings()
            self.parent_window.settings.update(settings)
            self.parent_window.save_settings()
            # 3Dビューの設定を適用
            self.parent_window.apply_3d_settings()
            # 現在の分子を再描画（設定変更を反映）
            if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
            # ステータスバーに適用完了を表示
            self.parent_window.statusBar().showMessage("Settings applied successfully", 2000)

    def accept(self):
        """ダイアログの設定を適用してから閉じる"""
        # apply_settingsを呼び出して設定を適用
        self.apply_settings()
        super().accept()


class CustomQtInteractor(QtInteractor):
    def __init__(self, parent=None, main_window=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.main_window = main_window

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
        self.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.AddObserver("LeftButtonReleaseEvent", self.on_left_button_up)

    def on_left_button_down(self, obj, event):
        """
        クリック時の処理を振り分けます。
        原子を掴めた場合のみカスタム動作に入り、それ以外は親クラス（カメラ回転）に任せます。
        """
        mw = self.main_window
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
        
        # Ctrl+クリックでの原子選択（3D編集用）
        if is_ctrl_click and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

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
                            # 3D編集用の原子選択をトグル
                            mw.toggle_atom_selection_3d(int(closest_atom_idx))
                            return  # カメラ回転は無効

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

    def on_mouse_move(self, obj, event):
        """
        マウス移動時の処理。原子ドラッグ中か、それ以外（カメラ回転＋ホバー）かをハンドリングします。
        """
        mw = self.main_window
        interactor = self.GetInteractor()

        # マウス移動があったことを記録
        if self._mouse_press_pos is not None:
            current_pos = interactor.GetEventPosition()
            if abs(current_pos[0] - self._mouse_press_pos[0]) > 3 or abs(current_pos[1] - self._mouse_press_pos[1]) > 3:
                self._mouse_moved_during_drag = True

        if self._is_dragging_atom:
            # カスタムの原子ドラッグ処理
            self.is_dragging = True
            atom_id = mw.dragged_atom_info['id']
            conf = mw.current_mol.GetConformer()
            renderer = mw.plotter.renderer
            current_display_pos = interactor.GetEventPosition()
            pos_3d = conf.GetAtomPosition(atom_id)
            renderer.SetWorldPoint(pos_3d.x, pos_3d.y, pos_3d.z, 1.0)
            renderer.WorldToDisplay()
            display_coords = renderer.GetDisplayPoint()
            new_display_pos = (current_display_pos[0], current_display_pos[1], display_coords[2])
            renderer.SetDisplayPoint(new_display_pos[0], new_display_pos[1], new_display_pos[2])
            renderer.DisplayToWorld()
            new_world_coords_tuple = renderer.GetWorldPoint()
            new_world_coords = list(new_world_coords_tuple)[:3]
            mw.atom_positions_3d[atom_id] = new_world_coords
            mw.glyph_source.points = mw.atom_positions_3d
            mw.glyph_source.Modified()
            conf.SetAtomPosition(atom_id, new_world_coords)
            interactor.Render()
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
                if mw.current_mol:
                    mw.draw_molecule_3d(mw.current_mol)
                mw.push_undo_state()
            mw.dragged_atom_info = None
        else:
            # カメラ回転の後始末を親クラスに任せます
            super().OnLeftButtonUp()

        # 状態をリセット
        self._is_dragging_atom = False
        self.is_dragging = False # is_draggingもリセット
        self._mouse_press_pos = None  # マウスプレス位置もリセット
        
        # ピックリセットは測定モードで実際に問題が発生した場合のみ行う
        # （通常のドラッグ回転では行わない）
        
        # ボタンを離した後のカーソル表示を最新の状態に更新
        self.on_mouse_move(obj, event)

        # 2Dビューにフォーカスを戻し、ショートカットキーなどが使えるようにする
        if mw and mw.view_2d:
            mw.view_2d.setFocus()

class MainWindow(QMainWindow):

    start_calculation = pyqtSignal(str)
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
        self.axes_actor = None
        self.axes_widget = None
        self.undo_stack = []
        self.redo_stack = []
        self.mode_actions = {} 
        
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
        # 初期状態は_enable_3d_features(False)で統一的に設定
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

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
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

        # --- 結合ボタンのアイコンを生成するヘルパー関数 ---
        def create_bond_icon(bond_type, size=32):
            pixmap = QPixmap(size, size)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)

            p1 = QPointF(6, size / 2)
            p2 = QPointF(size - 6, size / 2)
            line = QLineF(p1, p2)

            painter.setPen(QPen(Qt.GlobalColor.black, 2))
            painter.setBrush(QBrush(Qt.GlobalColor.black))

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
                    painter.setPen(QPen(Qt.GlobalColor.black, 1.5))
                    painter.drawLine(start_pt - offset, start_pt + offset)

            elif bond_type == 'ez_toggle':
                # アイコン下部に二重結合を描画
                p1 = QPointF(6, size * 0.75)
                p2 = QPointF(size - 6, size * 0.75)
                line = QLineF(p1, p2)
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * 2.0
                painter.setPen(QPen(Qt.GlobalColor.black, 2))
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))
                # 上部に "Z⇌E" のテキストを描画
                painter.setPen(QPen(Qt.GlobalColor.black, 1))
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

        toolbar.addSeparator()
        toolbar.addWidget(QLabel(" Templates:"))
        
        # --- アイコンを生成するヘルパー関数 ---
        def create_template_icon(n, is_benzene=False):
            size = 32
            pixmap = QPixmap(size, size)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setPen(QPen(Qt.GlobalColor.black, 2))

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
            toolbar.addAction(action)
            self.tool_group.addAction(action)

        # 初期モードを'select'から'atom_C'（炭素原子描画モード）に変更
        self.set_mode('atom_C')
        # 対応するツールバーの'C'ボタンを選択状態にする
        if 'atom_C' in self.mode_actions:
            self.mode_actions['atom_C'].setChecked(True)

        # スペーサーを追加して、次のウィジェットを右端に配置する
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        toolbar.addWidget(spacer)

        # 測定機能ボタンを追加
        self.measurement_action = QAction("Measure", self, checkable=True)
        self.measurement_action.setToolTip("Enable distance, angle, and dihedral measurement in 3D view")
        # 初期状態は_enable_3d_features(False)で統一的に設定
        self.measurement_action.triggered.connect(self.toggle_measurement_mode)
        toolbar.addAction(self.measurement_action)

        self.edit_3d_action = QAction("3D Edit", self, checkable=True)
        self.edit_3d_action.setToolTip("Toggle 3D atom editing mode (Hold Alt for temporary mode)")
        self.edit_3d_action.setEnabled(False)
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
        
        # NEW メニュー（Clear allと同じ機能）
        new_action = QAction("&New", self)
        new_action.setShortcut("Ctrl+N")
        new_action.triggered.connect(self.clear_all)
        file_menu.addAction(new_action)
        
        file_menu.addSeparator()
        
        load_mol_action = QAction("Import MOL/SDF...", self); load_mol_action.triggered.connect(self.load_mol_file)
        file_menu.addAction(load_mol_action)
        import_smiles_action = QAction("Import SMILES...", self)
        import_smiles_action.triggered.connect(self.import_smiles_dialog)
        file_menu.addAction(import_smiles_action)
        import_inchi_action = QAction("Import InChI...", self)
        import_inchi_action.triggered.connect(self.import_inchi_dialog)
        file_menu.addAction(import_inchi_action)

        file_menu.addSeparator()
        load_3d_mol_action = QAction("Load 3D MOL/SDF (3D Only)...", self)
        load_3d_mol_action.triggered.connect(self.load_mol_file_for_3d_viewing)
        file_menu.addAction(load_3d_mol_action)
        
        load_3d_xyz_action = QAction("Load 3D XYZ (3D Only)...", self)
        load_3d_xyz_action.triggered.connect(self.load_xyz_for_3d_viewing)
        file_menu.addAction(load_3d_xyz_action)

        file_menu.addSeparator()
        save_mol_action = QAction("Save 2D as MOL...", self); save_mol_action.triggered.connect(self.save_as_mol)
        file_menu.addAction(save_mol_action)
        
        save_3d_mol_action = QAction("Save 3D as MOL...", self); save_3d_mol_action.triggered.connect(self.save_3d_as_mol)
        file_menu.addAction(save_3d_mol_action)
        
        save_xyz_action = QAction("Save 3D as XYZ...", self); save_xyz_action.triggered.connect(self.save_as_xyz)
        file_menu.addAction(save_xyz_action)
        file_menu.addSeparator()
        save_raw_action = QAction("Save Project...", self); save_raw_action.triggered.connect(self.save_raw_data)
        save_raw_action.setShortcut(QKeySequence.StandardKey.Save) 
        file_menu.addAction(save_raw_action)
        load_raw_action = QAction("Open Project...", self); load_raw_action.triggered.connect(self.load_raw_data)
        load_raw_action.setShortcut(QKeySequence.StandardKey.Open) 
        file_menu.addAction(load_raw_action)
        
        file_menu.addSeparator()
        
        export_2d_png_action = QAction("Export 2D as PNG...", self)
        export_2d_png_action.triggered.connect(self.export_2d_png)
        file_menu.addAction(export_2d_png_action)

        export_3d_png_action = QAction("Export 3D as PNG...", self)
        export_3d_png_action.triggered.connect(self.export_3d_png)
        file_menu.addAction(export_3d_png_action)
        
        # 3Dエクスポートサブメニュー
        export_3d_menu = file_menu.addMenu("Export 3D Files")
        
        export_mol_action = QAction("Export as MOL...", self)
        export_mol_action.triggered.connect(self.save_3d_as_mol)
        export_3d_menu.addAction(export_mol_action)

        export_xyz_action = QAction("Export as XYZ...", self)
        export_xyz_action.triggered.connect(self.save_as_xyz)
        export_3d_menu.addAction(export_xyz_action)
        
        export_3d_menu.addSeparator()
        
        export_stl_action = QAction("Export as STL...", self)
        export_stl_action.triggered.connect(self.export_stl)
        export_3d_menu.addAction(export_stl_action)
        
        export_obj_action = QAction("Export as OBJ/MTL (with colors)...", self)
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

        remove_hydrogen_action = QAction("Remove Hydrogen", self)
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

        # 3D Edit menu
        edit_3d_menu = menu_bar.addMenu("3D &Edit")
        
        # Translation action
        translation_action = QAction("Translation...", self)
        translation_action.triggered.connect(self.open_translation_dialog)
        translation_action.setEnabled(False)
        edit_3d_menu.addAction(translation_action)
        self.translation_action = translation_action
        
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
        
        # Plane alignment submenu (旧Planarization)
        plane_align_menu = align_menu.addMenu("Plane")
        
        planar_xy_action = QAction("XY-plane", self)
        planar_xy_action.triggered.connect(lambda: self.open_planarization_dialog('xy'))
        planar_xy_action.setEnabled(False)
        plane_align_menu.addAction(planar_xy_action)
        self.planar_xy_action = planar_xy_action
        
        planar_xz_action = QAction("XZ-plane", self)
        planar_xz_action.triggered.connect(lambda: self.open_planarization_dialog('xz'))
        planar_xz_action.setEnabled(False)
        plane_align_menu.addAction(planar_xz_action)
        self.planar_xz_action = planar_xz_action
        
        planar_yz_action = QAction("YZ-plane", self)
        planar_yz_action.triggered.connect(lambda: self.open_planarization_dialog('yz'))
        planar_yz_action.setEnabled(False)
        plane_align_menu.addAction(planar_yz_action)
        self.planar_yz_action = planar_yz_action

        edit_3d_menu.addSeparator()

        # Mirror action
        mirror_action = QAction("Mirror...", self)
        mirror_action.triggered.connect(self.open_mirror_dialog)
        mirror_action.setEnabled(False)
        edit_3d_menu.addAction(mirror_action)
        self.mirror_action = mirror_action
        
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
        
        #edit_3d_menu.addSeparator()
        
        # Symmetrize action
        #symmetrize_action = QAction("Symmetrize...", self)
        #symmetrize_action.triggered.connect(self.open_symmetrize_dialog)
        #symmetrize_action.setEnabled(False)
        #edit_3d_menu.addAction(symmetrize_action)
        #self.symmetrize_action = symmetrize_action
        
        

        analysis_menu = menu_bar.addMenu("&Analysis")
        self.analysis_action = QAction("Show Analysis...", self)
        self.analysis_action.triggered.connect(self.open_analysis_window)
        self.analysis_action.setEnabled(False)
        analysis_menu.addAction(self.analysis_action)

        settings_menu = menu_bar.addMenu("&Settings")
        view_settings_action = QAction("3D View Settings...", self)
        view_settings_action.triggered.connect(self.open_settings_dialog)
        settings_menu.addAction(view_settings_action)

        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("About", self)
        about_action.triggered.connect(lambda: QMessageBox.about(
            self,
            "About MoleditPy",
            f"MoleditPy Ver. {VERSION}\nAuthor: Hiromichi Yokoyama\nLicense: Apache-2.0"
        ))
        help_menu.addAction(about_action)

        github_action = QAction("GitHub", self)
        github_action.triggered.connect(
            lambda: QDesktopServices.openUrl(QUrl("https://github.com/HiroYokoyama/python_molecular_editor"))
        )
        help_menu.addAction(github_action)
        
        # 3D関連機能の初期状態を統一的に設定
        self._enable_3d_features(False)
        
    def init_worker_thread(self):
        self.thread=QThread();self.worker=CalculationWorker();self.worker.moveToThread(self.thread)
        self.start_calculation.connect(self.worker.run_calculation)
        self.worker.finished.connect(self.on_calculation_finished); self.worker.error.connect(self.on_calculation_error)
        self.worker.status_update.connect(self.update_status_bar)
        self.thread.start()

    def update_status_bar(self, message):
        """ワーカースレッドからのメッセージでステータスバーを更新するスロット"""
        self.statusBar().showMessage(message)

    def set_mode(self, mode_str):
        self.scene.mode = mode_str
        self.view_2d.setMouseTracking(True) 
        if not mode_str.startswith('template'):
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
        if mode_str in self.mode_actions:
            self.mode_actions[mode_str].setChecked(True)

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
            self.statusBar().showMessage(f"Copied {len(fragment_atoms)} atoms and {len(fragment_bonds)} bonds.", 2000)
            
        except Exception as e:
            print(f"Error during copy operation: {e}")
            import traceback
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during copy operation: {e}", 5000)

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
            import traceback
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during cut operation: {e}", 5000)

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
                self.statusBar().showMessage("Error: Invalid clipboard data format", 3000)
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
            import traceback
            traceback.print_exc()
            self.statusBar().showMessage(f"Error during paste operation: {e}", 5000)
        self.statusBar().showMessage(f"Pasted {len(new_atoms)} atoms.", 2000)
        self.activate_select_mode()

    def remove_hydrogen_atoms(self):
        """2Dビューで水素原子とその結合を削除する"""
        try:
            hydrogen_atoms = []
            
            # 水素原子を特定
            for atom_id, atom_data in self.data.atoms.items():
                if atom_data['symbol'] == 'H':
                    hydrogen_atoms.append(atom_data['item'])
            
            if not hydrogen_atoms:
                self.statusBar().showMessage("No hydrogen atoms found to remove.", 2000)
                return
            
            # 削除を実行
            if self.scene.delete_items(set(hydrogen_atoms)):
                self.push_undo_state()
                self.statusBar().showMessage(f"Removed {len(hydrogen_atoms)} hydrogen atoms.", 2000)
            else:
                self.statusBar().showMessage("Failed to remove hydrogen atoms.", 3000)
                
        except Exception as e:
            print(f"Error during hydrogen removal: {e}")
            import traceback
            traceback.print_exc()
            self.statusBar().showMessage(f"Error removing hydrogen atoms: {e}", 5000)

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


    def activate_select_mode(self):
        self.set_mode('select')
        if 'select' in self.mode_actions:
            self.mode_actions['select'].setChecked(True)

    def trigger_conversion(self):
        self.scene.clear_all_problem_flags()
        mol = self.data.to_rdkit_mol(use_2d_stereo=False)
        if not mol or mol.GetNumAtoms() == 0:
            # 3Dビューと関連データをクリア
            self.plotter.clear()
            self.current_mol = None
            self.analysis_action.setEnabled(False)
            self.statusBar().showMessage("3D view cleared.")
            self.view_2d.setFocus() 
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

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus() 
            return

        if len(Chem.GetMolFrags(mol)) > 1:
            self.statusBar().showMessage("Error: 3D conversion not supported for multiple molecules.")
            self.view_2d.setFocus() 
            return
            
        mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
        
        # CRITICAL FIX: Manually add stereo information to MOL block for E/Z bonds
        # This ensures the stereo info survives the MOL block round-trip
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
        
        self.convert_button.setEnabled(False)
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
        text_actor.GetTextProperty().SetOpacity(1)
        self.plotter.render()
        self.start_calculation.emit(mol_block)
        
        self.view_2d.setFocus()

    def optimize_3d_structure(self):
        """現在の3D分子構造を力場で最適化する"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to optimize.")
            return

        self.statusBar().showMessage("Optimizing 3D structure...")
        QApplication.processEvents() # UIの更新を確実に行う

        try:
            # MMFF力場での最適化を試みる
            AllChem.MMFFOptimizeMolecule(self.current_mol)
        except Exception:
            # MMFFが失敗した場合、UFF力場でフォールバック
            try:
                AllChem.UFFOptimizeMolecule(self.current_mol)
            except Exception as e:
                self.statusBar().showMessage(f"3D optimization failed: {e}")
                return
        
        # 最適化後の構造で3Dビューを再描画
        try:
            # キラル中心を3D座標から再計算（R/Sのみ）
            if self.current_mol.GetNumConformers() > 0:
                Chem.AssignAtomChiralTagsFromStructure(self.current_mol, confId=0)
            
            self.update_chiral_labels() # キラル中心のラベルも更新
        except Exception:
            pass
            
        self.draw_molecule_3d(self.current_mol)
        
        self.statusBar().showMessage("3D structure optimization successful.")
        self.push_undo_state() # Undo履歴に保存
        self.view_2d.setFocus()

    def on_calculation_finished(self, mol):
        self.dragged_atom_info = None
        self.current_mol = mol
        self.is_xyz_derived = False  # 2Dから生成した3D構造はXYZ由来ではない
        
        # 原子プロパティを復元（ワーカープロセスで失われたため）
        if hasattr(self, 'original_atom_properties'):
            for i, original_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(i)
                    atom.SetIntProp("_original_atom_id", original_id)
        
        # 原子IDマッピングを作成
        self.create_atom_id_mapping()
        
        # キラル中心を3D座標から再計算（R/Sのみ）
        try:
            if mol.GetNumConformers() > 0:
                Chem.AssignAtomChiralTagsFromStructure(mol, confId=0)
            
            self.update_chiral_labels()
        except Exception:
            # 念のためエラーを握り潰して UI を壊さない
            pass

        self.draw_molecule_3d(mol)

        #self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
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

    def on_calculation_error(self, error_message):
        self.plotter.clear()
        self.dragged_atom_info = None
        self.statusBar().showMessage(f"Error: {error_message}")
        self.cleanup_button.setEnabled(True)
        self.convert_button.setEnabled(True)
        self.analysis_action.setEnabled(False)
        self.edit_3d_action.setEnabled(False)
        # Disable 3D editing menu items
        self._enable_3d_edit_actions(False)
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
        
        self.update_implicit_hydrogens()
        self.update_realtime_info()
        self.update_undo_redo_actions()

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

    def clear_all(self):

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
        
        # 3D関連機能を統一的に無効化
        self._enable_3d_features(False)
        
        # Undo/Redoスタックをリセットする
        self.reset_undo_stack()
        
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
        if not self.data.atoms:
            return

        try:
            mol = self.data.to_rdkit_mol()
            if mol is None:
                # 構造が不正な場合、全原子の水素カウントを0に戻して再描画
                for atom_data in self.data.atoms.values():
                    if atom_data.get('item') and atom_data['item'].implicit_h_count != 0:
                        atom_data['item'].implicit_h_count = 0
                        atom_data['item'].update()
                return
            
            items_to_update = []
            for atom in mol.GetAtoms():
                if atom.HasProp("_original_atom_id"):
                    original_id = atom.GetIntProp("_original_atom_id")
                    if original_id in self.data.atoms:
                        item = self.data.atoms[original_id].get('item')
                        if item:
                            h_count = atom.GetNumImplicitHs()
                            if item.implicit_h_count != h_count:
                                item.prepareGeometryChange()
                                item.implicit_h_count = h_count
                                items_to_update.append(item)
            
            # カウントが変更されたアイテムのみ再描画をトリガー
            for item in items_to_update:
                item.update()
        except Exception:
            # 編集中に一時的に発生するエラーなどで計算が失敗してもアプリは継続
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
            QTimer.singleShot(0, self.fit_to_view)
            
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid SMILES: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading from SMILES: {e}")
            import traceback
            traceback.print_exc()

    def load_from_inchi(self, inchi_string):
        """InChI文字列から分子を読み込み、2Dエディタに表示する"""
        try:
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
            QTimer.singleShot(0, self.fit_to_view)
            
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid InChI: {e}")
        except Exception as e:
            self.statusBar().showMessage(f"Error loading from InChI: {e}")
            import traceback
            traceback.print_exc()

    def load_mol_file(self, file_path=None):
        if not file_path:
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getOpenFileName(self, "Import MOL File", "", "Chemical Files (*.mol *.sdf);;All Files (*)", options=options)
            if not file_path: 
                return

        try:
            self.dragged_atom_info = None
            suppl = Chem.SDMolSupplier(file_path, removeHs=False)
            mol = next(suppl, None)
            if mol is None: raise ValueError("Failed to read molecule from file.")

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
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid MOL file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading file: {e}")
            import traceback
            traceback.print_exc()
    
    def load_mol_for_3d_viewing(self):
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Load 3D MOL (View Only)", "", "Chemical Files (*.mol *.sdf);;All Files (*)", options=options)
        if not file_path:
            return

        try:
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

        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}", 5000)
            self.restore_ui_for_editing()
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid 3D MOL file: {e}", 5000)
            self.restore_ui_for_editing()
        except Exception as e:
            self.statusBar().showMessage(f"Error loading 3D file: {e}", 5000)
            self.restore_ui_for_editing()
            import traceback
            traceback.print_exc()

    def load_xyz_for_3d_viewing(self, file_path=None):
        """XYZファイルを読み込んで3Dビューアで表示する"""
        if not file_path:
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getOpenFileName(self, "Load 3D XYZ (View Only)", "", "XYZ Files (*.xyz);;All Files (*)", options=options)
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
            self.current_mol = mol
            self.is_xyz_derived = True  # XYZ由来であることを記録
            
            # XYZファイル読み込み時はマッピングをクリア（2D構造がないため）
            self.atom_id_to_rdkit_idx_map = {}
            
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

        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}", 5000)
            self.restore_ui_for_editing()
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid XYZ file: {e}", 5000)
            self.restore_ui_for_editing()
        except Exception as e:
            self.statusBar().showMessage(f"Error loading XYZ file: {e}", 5000)
            self.restore_ui_for_editing()
            import traceback
            traceback.print_exc()

    def load_xyz_file(self, file_path):
        """XYZファイルを読み込んでRDKitのMolオブジェクトを作成する"""
        import math
        from rdkit.Chem import rdGeometry
        
        try:
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
                except:
                    # 認識されない場合、最初の文字を大文字にして再試行
                    symbol = symbol.capitalize()
                    try:
                        test_atom = Chem.Atom(symbol)
                    except:
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
            
            # 結合を推定（距離ベース）
            self.estimate_bonds_from_distances(mol)
            
            # 分子を最終化
            try:
                mol = mol.GetMol()
                # 基本的な妥当性チェック
                if mol is None:
                    raise ValueError("Failed to create valid molecule object")
                Chem.SanitizeMol(mol)
            except Exception as e:
                # 化学的に不正な構造でも表示は可能にする
                mol = mol.GetMol()
                if mol is None:
                    raise ValueError("Failed to create molecule object")
            
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
        from rdkit.Chem import rdMolTransforms
        
        # 一般的な共有結合半径（Ångström）- より正確な値
        covalent_radii = {
            'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76,
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41,
            'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
            'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
            'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
            'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40
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


    def save_raw_data(self):
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            save_data = self.get_current_state()
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeraw'): 
                file_path += '.pmeraw'
                
            with open(file_path, 'wb') as f: 
                pickle.dump(save_data, f)
            self.statusBar().showMessage(f"Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving project file: {e}")
            import traceback
            traceback.print_exc()


    def load_raw_data(self, file_path=None):
        if not file_path:
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getOpenFileName(self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
            if not file_path: 
                return
        
        try:
            with open(file_path, 'rb') as f: 
                loaded_data = pickle.load(f)
            self.restore_ui_for_editing()
            self.set_state_from_data(loaded_data)
            self.statusBar().showMessage(f"Project loaded from {file_path}")
            self.reset_undo_stack()
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.UnpicklingError as e:
            self.statusBar().showMessage(f"Invalid project file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading project file: {e}")
            import traceback
            traceback.print_exc()

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
            
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(self, "Save 2D MOL File", "", "MOL Files (*.mol);;All Files (*)", options=options)
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
            import traceback
            traceback.print_exc()
            
    def save_3d_as_mol(self):
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        try:
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getSaveFileName(self, "Save 3D MOL File", "", "MOL Files (*.mol);;All Files (*)", options=options)
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
            import traceback
            traceback.print_exc()

    def save_as_xyz(self):
        if not self.current_mol: self.statusBar().showMessage("Error: Please generate a 3D structure first."); return
        options=QFileDialog.Option.DontUseNativeDialog
        file_path,_=QFileDialog.getSaveFileName(self,"Save 3D XYZ File","","XYZ Files (*.xyz);;All Files (*)",options=options)
        if file_path:
            if not file_path.lower().endswith('.xyz'): file_path += '.xyz'
            try:
                conf=self.current_mol.GetConformer(); num_atoms=self.current_mol.GetNumAtoms()
                xyz_lines=[str(num_atoms)]; smiles=Chem.MolToSmiles(Chem.RemoveHs(self.current_mol))
                xyz_lines.append(f"Generated by MoleditPy Ver. {VERSION}")
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
            
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as STL", "", "STL Files (*.stl);;All Files (*)", options=options
        )
        
        if not file_path:
            return
            
        try:
            import pyvista as pv
            
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
        """OBJ/MTLファイルとしてエクスポート（個別色付き）"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as OBJ/MTL (with colors)", "", "OBJ Files (*.obj);;All Files (*)", options=options
        )
        
        if not file_path:
            return
            
        try:
            import pyvista as pv
            
            # 3Dビューから個別のメッシュデータを取得（オブジェクトごとに分離）
            meshes_with_colors = self.export_individual_meshes_from_3d_view()
            
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

    def export_individual_meshes_from_3d_view(self):
        """3Dビューから個別のメッシュデータを色情報とともに取得（分子データベース）"""
        try:
            import pyvista as pv
            import numpy as np
            import vtk
            
            meshes_with_colors = []
            
            # 分子データから直接色情報を取得
            if self.current_mol is None:
                return []
            
            mol = self.current_mol
            conf = mol.GetConformer()
            
            # 原子ごとに個別のメッシュを作成
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                symbol = atom.GetSymbol()
                position = np.array(conf.GetAtomPosition(i))
                
                # CPKカラーから色を取得
                color_rgb = CPK_COLORS_PV.get(symbol, [0.5, 0.5, 0.5])
                color = [int(c * 255) for c in color_rgb]
                
                # 原子の半径を取得（設定値を反映）
                if self.current_3d_style == 'cpk':
                    base_radius = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                    radius = base_radius * self.settings.get('cpk_atom_scale', 1.0)
                elif self.current_3d_style == 'wireframe':
                    radius = 0.1  # Wireframeでは固定
                elif self.current_3d_style == 'stick':
                    radius = self.settings.get('stick_atom_radius', 0.15)
                else:  # ball_and_stick
                    base_radius = VDW_RADII.get(symbol, 0.4)
                    radius = base_radius * self.settings.get('ball_stick_atom_scale', 1.0)
                
                # Wireframeモードでは原子を描画しない
                if self.current_3d_style != 'wireframe':
                    # 原子球体の解像度を設定値に基づいて決定
                    if self.current_3d_style == 'cpk':
                        atom_resolution = self.settings.get('cpk_resolution', 32)
                    elif self.current_3d_style == 'stick':
                        atom_resolution = self.settings.get('stick_resolution', 16)
                    else:  # ball_and_stick
                        atom_resolution = self.settings.get('ball_stick_resolution', 16)
                    
                    # 個別の原子球体を作成
                    sphere = pv.Sphere(radius=radius, center=position, theta_resolution=atom_resolution, phi_resolution=atom_resolution)
                    
                    meshes_with_colors.append({
                        'mesh': sphere,
                        'color': color,
                        'name': f'atom_{i}_{symbol}',
                        'type': 'atom',
                        'atom_index': i
                    })
            
            # 結合ごとに個別のメッシュを作成
            if self.current_3d_style in ['ball_and_stick', 'wireframe', 'stick']:
                # スタイルに応じてボンドの太さを設定（設定値を反映）
                if self.current_3d_style == 'wireframe':
                    cyl_radius = self.settings.get('wireframe_bond_radius', 0.01)
                    bond_resolution = self.settings.get('wireframe_resolution', 6)
                elif self.current_3d_style == 'stick':
                    cyl_radius = self.settings.get('stick_bond_radius', 0.15)
                    bond_resolution = self.settings.get('stick_resolution', 16)
                else:  # ball_and_stick
                    cyl_radius = self.settings.get('ball_stick_bond_radius', 0.1)
                    bond_resolution = self.settings.get('ball_stick_resolution', 16)
                
                bond_counter = 0
                
                for bond in mol.GetBonds():
                    begin_atom_idx = bond.GetBeginAtomIdx()
                    end_atom_idx = bond.GetEndAtomIdx()
                    sp = np.array(conf.GetAtomPosition(begin_atom_idx))
                    ep = np.array(conf.GetAtomPosition(end_atom_idx))
                    bt = bond.GetBondType()
                    
                    # 結合の色を原子の色から決定
                    begin_symbol = mol.GetAtomWithIdx(begin_atom_idx).GetSymbol()
                    end_symbol = mol.GetAtomWithIdx(end_atom_idx).GetSymbol()
                    begin_color = CPK_COLORS_PV.get(begin_symbol, [0.5, 0.5, 0.5])
                    end_color = CPK_COLORS_PV.get(end_symbol, [0.5, 0.5, 0.5])
                    begin_color_rgb = [int(c * 255) for c in begin_color]
                    end_color_rgb = [int(c * 255) for c in end_color]
                    
                    c = (sp + ep) / 2
                    d = ep - sp
                    h = np.linalg.norm(d)
                    if h == 0: 
                        continue
                    
                    # 結合タイプに応じて描画
                    if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                        # 単結合・芳香族結合の処理
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickはグレー（単一の円柱）
                            cyl = pv.Cylinder(center=c, direction=d, radius=cyl_radius, height=h, resolution=bond_resolution)
                            meshes_with_colors.append({
                                'mesh': cyl,
                                'color': [127, 127, 127],
                                'name': f'bond_{bond_counter}_{begin_symbol}{begin_atom_idx}-{end_symbol}{end_atom_idx}',
                                'type': 'bond'
                            })
                        else:
                            # その他（stick, wireframe）は中央で色が変わる2つの円柱
                            mid_point = (sp + ep) / 2
                            
                            # 前半（開始原子の色）
                            cyl1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                            meshes_with_colors.append({
                                'mesh': cyl1,
                                'color': begin_color_rgb,
                                'name': f'bond_{bond_counter}_start_{begin_symbol}{begin_atom_idx}',
                                'type': 'bond_segment'
                            })
                            
                            # 後半（終了原子の色）
                            cyl2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                            meshes_with_colors.append({
                                'mesh': cyl2,
                                'color': end_color_rgb,
                                'name': f'bond_{bond_counter}_end_{end_symbol}{end_atom_idx}',
                                'type': 'bond_segment'
                            })
                    
                    elif bt == Chem.rdchem.BondType.DOUBLE:
                        # 二重結合の処理
                        # 結合ベクトルに垂直な方向を求める
                        d_norm = d / h
                        # 適当な方向ベクトルとの外積で垂直ベクトルを作成
                        if abs(d_norm[0]) < 0.9:
                            perp = np.cross(d_norm, [1, 0, 0])
                        else:
                            perp = np.cross(d_norm, [0, 1, 0])
                        perp = perp / np.linalg.norm(perp)
                        
                        # 2本の平行な円柱を作成（オフセット距離を設定値に基づいて調整）
                        offset_factor = self.settings.get('double_bond_offset_factor', 2.0)
                        offset_distance = cyl_radius * offset_factor
                        offset1 = c + perp * offset_distance
                        offset2 = c - perp * offset_distance
                        
                        # 各円柱の半径を設定値に基づいて調整
                        radius_factor = self.settings.get('double_bond_radius_factor', 0.8)
                        double_radius = cyl_radius * radius_factor
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickはグレー
                            cyl1 = pv.Cylinder(center=offset1, direction=d, radius=double_radius, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=offset2, direction=d, radius=double_radius, height=h, resolution=bond_resolution)
                            
                            meshes_with_colors.append({
                                'mesh': cyl1,
                                'color': [127, 127, 127],
                                'name': f'bond_{bond_counter}_double1_{begin_symbol}{begin_atom_idx}-{end_symbol}{end_atom_idx}',
                                'type': 'double_bond'
                            })
                            meshes_with_colors.append({
                                'mesh': cyl2,
                                'color': [127, 127, 127],
                                'name': f'bond_{bond_counter}_double2_{begin_symbol}{begin_atom_idx}-{end_symbol}{end_atom_idx}',
                                'type': 'double_bond'
                            })
                        else:
                            # その他は色分け
                            mid_point = (sp + ep) / 2
                            
                            # 第1本目
                            cyl1_start = pv.Cylinder(center=(sp + mid_point) / 2 + perp * offset_distance, direction=d, radius=double_radius, height=h/2, resolution=bond_resolution)
                            cyl1_end = pv.Cylinder(center=(mid_point + ep) / 2 + perp * offset_distance, direction=d, radius=double_radius, height=h/2, resolution=bond_resolution)
                            
                            # 第2本目
                            cyl2_start = pv.Cylinder(center=(sp + mid_point) / 2 - perp * offset_distance, direction=d, radius=double_radius, height=h/2, resolution=bond_resolution)
                            cyl2_end = pv.Cylinder(center=(mid_point + ep) / 2 - perp * offset_distance, direction=d, radius=double_radius, height=h/2, resolution=bond_resolution)
                            
                            meshes_with_colors.extend([
                                {'mesh': cyl1_start, 'color': begin_color_rgb, 'name': f'bond_{bond_counter}_double1_start_{begin_symbol}{begin_atom_idx}', 'type': 'double_bond_segment'},
                                {'mesh': cyl1_end, 'color': end_color_rgb, 'name': f'bond_{bond_counter}_double1_end_{end_symbol}{end_atom_idx}', 'type': 'double_bond_segment'},
                                {'mesh': cyl2_start, 'color': begin_color_rgb, 'name': f'bond_{bond_counter}_double2_start_{begin_symbol}{begin_atom_idx}', 'type': 'double_bond_segment'},
                                {'mesh': cyl2_end, 'color': end_color_rgb, 'name': f'bond_{bond_counter}_double2_end_{end_symbol}{end_atom_idx}', 'type': 'double_bond_segment'}
                            ])
                    
                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        # 三重結合の処理
                        # 結合ベクトルに垂直な2つの方向を求める
                        d_norm = d / h
                        if abs(d_norm[0]) < 0.9:
                            perp1 = np.cross(d_norm, [1, 0, 0])
                        else:
                            perp1 = np.cross(d_norm, [0, 1, 0])
                        perp1 = perp1 / np.linalg.norm(perp1)
                        perp2 = np.cross(d_norm, perp1)
                        perp2 = perp2 / np.linalg.norm(perp2)
                        
                        # 3本の円柱を120度間隔で配置（オフセット距離を設定値に基づいて調整）
                        offset_factor = self.settings.get('triple_bond_offset_factor', 2.0)
                        offset_distance = cyl_radius * offset_factor
                        angle_step = 2 * np.pi / 3  # 120度
                        
                        # 各円柱の半径を設定値に基づいて調整
                        radius_factor = self.settings.get('triple_bond_radius_factor', 0.7)
                        triple_radius = cyl_radius * radius_factor
                        
                        for i in range(3):
                            angle = i * angle_step
                            offset = perp1 * np.cos(angle) * offset_distance + perp2 * np.sin(angle) * offset_distance
                            
                            if self.current_3d_style == 'ball_and_stick':
                                # Ball and stickはグレー
                                cyl = pv.Cylinder(center=c + offset, direction=d, radius=triple_radius, height=h, resolution=bond_resolution)
                                meshes_with_colors.append({
                                    'mesh': cyl,
                                    'color': [127, 127, 127],
                                    'name': f'bond_{bond_counter}_triple{i+1}_{begin_symbol}{begin_atom_idx}-{end_symbol}{end_atom_idx}',
                                    'type': 'triple_bond'
                                })
                            else:
                                # その他は色分け
                                mid_point = (sp + ep) / 2
                                
                                cyl_start = pv.Cylinder(center=(sp + mid_point) / 2 + offset, direction=d, radius=triple_radius, height=h/2, resolution=bond_resolution)
                                cyl_end = pv.Cylinder(center=(mid_point + ep) / 2 + offset, direction=d, radius=triple_radius, height=h/2, resolution=bond_resolution)
                                
                                meshes_with_colors.extend([
                                    {'mesh': cyl_start, 'color': begin_color_rgb, 'name': f'bond_{bond_counter}_triple{i+1}_start_{begin_symbol}{begin_atom_idx}', 'type': 'triple_bond_segment'},
                                    {'mesh': cyl_end, 'color': end_color_rgb, 'name': f'bond_{bond_counter}_triple{i+1}_end_{end_symbol}{end_atom_idx}', 'type': 'triple_bond_segment'}
                                ])
                    
                    else:
                        # その他の結合タイプ（単結合として処理）
                        if self.current_3d_style == 'ball_and_stick':
                            cyl = pv.Cylinder(center=c, direction=d, radius=cyl_radius, height=h, resolution=bond_resolution)
                            meshes_with_colors.append({
                                'mesh': cyl,
                                'color': [127, 127, 127],
                                'name': f'bond_{bond_counter}_other_{begin_symbol}{begin_atom_idx}-{end_symbol}{end_atom_idx}',
                                'type': 'bond'
                            })
                        else:
                            mid_point = (sp + ep) / 2
                            
                            cyl1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                            
                            meshes_with_colors.extend([
                                {'mesh': cyl1, 'color': begin_color_rgb, 'name': f'bond_{bond_counter}_other_start_{begin_symbol}{begin_atom_idx}', 'type': 'bond_segment'},
                                {'mesh': cyl2, 'color': end_color_rgb, 'name': f'bond_{bond_counter}_other_end_{end_symbol}{end_atom_idx}', 'type': 'bond_segment'}
                            ])
                    
                    bond_counter += 1
            
            return meshes_with_colors
            
        except Exception as e:
            return []

    def create_multi_material_obj(self, meshes_with_colors, obj_path, mtl_path):
        """複数のマテリアルを持つOBJファイルとMTLファイルを作成（改良版）"""
        try:
            import os
            
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
            
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as Color STL", "", "STL Files (*.stl);;All Files (*)", options=options
        )
        
        if not file_path:
            return
            
        try:
            import pyvista as pv
            import numpy as np
            
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
            import pyvista as pv
            import numpy as np
            import vtk
            
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
            import pyvista as pv
            import numpy as np
            import vtk
            
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

    def export_2d_png(self):
        if not self.data.atoms:
            self.statusBar().showMessage("Nothing to export.", 2000)
            return

        options = QFileDialog.Option.DontUseNativeDialog
        filePath, _ = QFileDialog.getSaveFileName(self, "Export 2D as PNG", "", "PNG Files (*.png)", options=options)
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for a white background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.No)

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
                self.statusBar().showMessage("Error: Could not determine molecule bounds for export.", 5000)
                return

            if is_transparent:
                self.scene.setBackgroundBrush(QBrush(Qt.BrushStyle.NoBrush))
            else:
                self.scene.setBackgroundBrush(QBrush(QColor("#FFFFFF")))

            rect_to_render = molecule_bounds.adjusted(-20, -20, 20, 20)

            w = max(1, int(math.ceil(rect_to_render.width())))
            h = max(1, int(math.ceil(rect_to_render.height())))

            if w <= 0 or h <= 0:
                self.statusBar().showMessage("Error: Invalid image size calculated.", 5000)
                return

            image = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
            if is_transparent:
                image.fill(Qt.GlobalColor.transparent)
            else:
                image.fill(Qt.GlobalColor.white)

            painter = QPainter()
            ok = painter.begin(image)
            if not ok or not painter.isActive():
                self.statusBar().showMessage("Failed to start QPainter for image rendering.", 5000)
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
                self.statusBar().showMessage(f"2D view exported to {filePath}", 3000)
            else:
                self.statusBar().showMessage(f"Failed to save image. Check file path or permissions.", 5000)

        except Exception as e:
            self.statusBar().showMessage(f"An unexpected error occurred during 2D export: {e}", 5000)

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

        options = QFileDialog.Option.DontUseNativeDialog
        filePath, _ = QFileDialog.getSaveFileName(self, "Export 3D as PNG", "", "PNG Files (*.png)", options=options)
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for current background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.No)

        if reply == QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = (reply == QMessageBox.StandardButton.Yes)

        try:
            self.plotter.screenshot(filePath, transparent_background=is_transparent)
            self.statusBar().showMessage(f"3D view exported to {filePath}", 3000)
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting 3D PNG: {e}", 3000)


    def open_periodic_table_dialog(self):
        dialog=PeriodicTableDialog(self); dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action=self.tool_group.checkedAction()
        if checked_action: self.tool_group.setExclusive(False); checked_action.setChecked(False); self.tool_group.setExclusive(True)
        dialog.exec()

    def set_atom_from_periodic_table(self, symbol): 
        self.set_mode(f'atom_{symbol}')

   
    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")
        mol = self.data.to_rdkit_mol()
        if mol is None or mol.GetNumAtoms() == 0:
            self.statusBar().showMessage("Error: No atoms to optimize."); return

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
            for bond_data in self.data.bonds.values():
                if bond_data.get('item'):
                    bond_data['item'].update_position()

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
            if bond_data and 'item' in bond_data:
                bond_data['item'].update_position()
        
        # 重なり解消後に測定ラベルの位置を更新
        self.update_2d_measurement_labels()
        
        self.scene.update()
        self.push_undo_state()
        self.statusBar().showMessage("Resolved overlapping groups.", 2000)



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

                if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                    if self.current_3d_style == 'ball_and_stick':
                        # Ball and stickはグレー
                        cyl = pv.Cylinder(center=c, direction=d, radius=cyl_radius, height=h, resolution=bond_resolution)
                        actor = self.plotter.add_mesh(cyl, color='grey', **mesh_props)
                        self._3d_color_map[f'bond_{bond_counter}'] = [127, 127, 127]  # グレー
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
                    r = cyl_radius * 0.8
                    # 設定からオフセットファクターを取得
                    double_offset_factor = self.settings.get('double_bond_offset_factor', 2.0)
                    triple_offset_factor = self.settings.get('triple_bond_offset_factor', 2.0)
                    s = cyl_radius * 2.0  # デフォルト値
                    
                    if bt == Chem.rdchem.BondType.DOUBLE:
                        # 二重結合の場合、結合している原子の他の結合を考慮してオフセット方向を決定
                        off_dir = self._calculate_double_bond_offset(mol, bond, conf)
                        # 設定から二重結合のオフセットファクターを適用
                        s_double = cyl_radius * double_offset_factor
                        c1, c2 = c + off_dir * (s_double / 2), c - off_dir * (s_double / 2)
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickはグレー
                            cyl1 = pv.Cylinder(center=c1, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c2, direction=d, radius=r, height=h, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1, color='grey', **mesh_props)
                            self.plotter.add_mesh(cyl2, color='grey', **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1'] = [127, 127, 127]
                            self._3d_color_map[f'bond_{bond_counter}_2'] = [127, 127, 127]
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
                        # 三重結合
                        v_arb = np.array([0, 0, 1])
                        if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                        off_dir = np.cross(v1, v_arb)
                        off_dir /= np.linalg.norm(off_dir)
                        
                        # 設定から三重結合のオフセットファクターを適用
                        s_triple = cyl_radius * triple_offset_factor
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickはグレー
                            cyl1 = pv.Cylinder(center=c, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c + off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl3 = pv.Cylinder(center=c - off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1, color='grey', **mesh_props)
                            self.plotter.add_mesh(cyl2, color='grey', **mesh_props)
                            self.plotter.add_mesh(cyl3, color='grey', **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1'] = [127, 127, 127]
                            self._3d_color_map[f'bond_{bond_counter}_2'] = [127, 127, 127]
                            self._3d_color_map[f'bond_{bond_counter}_3'] = [127, 127, 127]
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
        
        # 各原子に対してラベルを作成
        texts = []
        positions = []
        
        for atom_idx in range(len(self.atom_positions_3d)):
            pos = self.atom_positions_3d[atom_idx]
            
            if self.atom_info_display_mode == 'id':
                # Original IDがある場合は優先表示、なければRDKitインデックス
                try:
                    if self.current_mol:
                        atom = self.current_mol.GetAtomWithIdx(atom_idx)
                        if atom.HasProp("_original_atom_id"):
                            original_id = atom.GetIntProp("_original_atom_id")
                            text = f"ID:{original_id}"
                        elif atom.HasProp("xyz_unique_id"):
                            unique_id = atom.GetIntProp("xyz_unique_id")
                            text = f"XYZ:{unique_id}"
                        else:
                            text = f"RDKit:{atom_idx}"
                    else:
                        text = f"RDKit:{atom_idx}"
                except Exception:
                    text = f"RDKit:{atom_idx}"
            elif self.atom_info_display_mode == 'rdkit_id':
                # RDKitで再生成された原子インデックスのみを表示
                text = f"RDKit:{atom_idx}"
            elif self.atom_info_display_mode == 'coords':
                text = f"({pos[0]:.2f},{pos[1]:.2f},{pos[2]:.2f})"
            elif self.atom_info_display_mode == 'symbol':
                if self.current_mol:
                    symbol = self.current_mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    text = f"{symbol}"
                else:
                    text = "?"
            else:
                continue
            
            texts.append(text)
            positions.append(pos)
        
        # すべてのラベルを一度に表示
        if texts and positions:
            try:
                label_actor = self.plotter.add_point_labels(
                    positions, texts,
                    point_size=12,
                    font_size=18,
                    text_color='black',
                    always_visible=True,
                    tolerance=0.01,
                    show_points=False
                )
                self.current_atom_info_labels = label_actor
            except Exception as e:
                print(f"Error adding atom info labels: {e}")

    def clear_all_atom_info_labels(self):
        """すべての原子情報ラベルをクリア"""
        if self.current_atom_info_labels is not None:
            try:
                self.plotter.remove_actor(self.current_atom_info_labels)
            except:
                pass
            self.current_atom_info_labels = None

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
        if self.settings != self.initial_settings:
            self.save_settings()
        reply = QMessageBox.question(self, 'Confirm Exit', 
                                     "Are you sure you want to exit?", 
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
                                     QMessageBox.StandardButton.No)

        if reply == QMessageBox.StandardButton.Yes:
            if self.scene and self.scene.template_preview:
                self.scene.template_preview.hide()

            self.thread.quit()
            self.thread.wait()
            
            event.accept()
        else:
            event.ignore()

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
            
        bounds = self.scene.itemsBoundingRect()
        visible_items_rect = QRectF()
        for item in self.scene.items():
            if item.isVisible() and not isinstance(item, TemplatePreviewItem):
                 if visible_items_rect.isEmpty():
                     visible_items_rect = item.sceneBoundingRect()
                 else:
                     visible_items_rect = visible_items_rect.united(item.sceneBoundingRect())
        
        if not visible_items_rect.isEmpty():
             self.view_2d.fitInView(visible_items_rect, Qt.AspectRatioMode.KeepAspectRatio)
             self.view_2d.scale(0.6, 0.6)
        else:
             self.reset_zoom()

    def toggle_3d_edit_mode(self, checked):
        """「3D Edit」ボタンの状態に応じて編集モードを切り替える"""
        if checked:
            # 3D Editモードをオンにする時は、Measurementモードを無効化
            if self.measurement_mode:
                self.measurement_action.setChecked(False)
                self.toggle_measurement_mode(False)
        
        self.is_3d_edit_mode = checked
        if checked:
            self.statusBar().showMessage("3D Edit Mode: ON.")
        else:
            self.statusBar().showMessage("3D Edit Mode: OFF.")
        self.view_2d.setFocus()

    def _setup_3d_picker(self):
        self.plotter.picker = vtk.vtkCellPicker()
        self.plotter.picker.SetTolerance(0.025)

        # 新しいカスタムスタイル（原子移動用）のインスタンスを作成
        style = CustomInteractorStyle(self)
        
        # 調査の結果、'style' プロパティへの代入が正しい設定方法と判明
        self.plotter.interactor.SetInteractorStyle(style)
        self.plotter.interactor.Initialize()
        
    def load_mol_file_for_3d_viewing(self, file_path=None):
        """MOL/SDFファイルを3Dビューアーで開く"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open MOL/SDF File", "", 
                "MOL/SDF Files (*.mol *.sdf);;All Files (*)"
            )
            if not file_path:
                return
        
        try:
            # RDKitでMOL/SDFファイルを読み込み
            if file_path.lower().endswith('.sdf'):
                suppl = Chem.SDMolSupplier(file_path)
                mol = next(suppl, None)
            else:
                mol = Chem.MolFromMolFile(file_path)
            
            if mol is None:
                self.statusBar().showMessage(f"Failed to load molecule from {file_path}")
                return
            
            # 3D座標がない場合は2Dから3D変換（最適化なし）
            if mol.GetNumConformers() == 0:
                self.statusBar().showMessage("No 3D coordinates found. Converting to 3D...")
                try:
                    AllChem.EmbedMolecule(mol)
                    # 最適化は実行しない
                except:
                    self.statusBar().showMessage("Failed to generate 3D coordinates")
                    return
            
            # 3Dビューアーに表示
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
        elif file_ext == 'pmeraw':
            self.load_raw_data(file_path=file_path)
        else:
            self.statusBar().showMessage(f"Unsupported file type: {file_ext}")
        
    def dragEnterEvent(self, event):
        """ウィンドウ全体で .pmeraw、.mol、.sdf、.xyz ファイルのドラッグを受け入れる"""
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if urls and urls[0].isLocalFile():
                file_path = urls[0].toLocalFile()
                if file_path.lower().endswith(('.pmeraw', '.mol', '.sdf', '.xyz')):
                    event.acceptProposedAction()
                    return
        event.ignore()

    def dropEvent(self, event):
        """ファイルがウィンドウ上でドロップされたときに呼び出される"""
        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            file_path = urls[0].toLocalFile()
            
            # ドロップ位置を取得
            drop_pos = event.position().toPoint()
            
            # 拡張子に応じて適切な読み込みメソッドを呼び出す
            if file_path.lower().endswith('.pmeraw'):
                self.load_raw_data(file_path=file_path)
                event.acceptProposedAction()
            elif file_path.lower().endswith(('.mol', '.sdf')):
                # MOL/SDFファイルは3Dビューアー領域にドロップされたかを判定
                plotter_widget = self.splitter.widget(1)  # 3Dビューアーウィジェット
                plotter_rect = plotter_widget.geometry()
                
                if plotter_rect.contains(drop_pos):
                    # 3Dビューアー領域にドロップ → 3Dビューアーで開く
                    self.load_mol_file_for_3d_viewing(file_path=file_path)
                else:
                    # 2D領域にドロップ → 従来の2D読み込み
                    if hasattr(self, 'load_mol_file'):
                        self.load_mol_file(file_path=file_path)
                    else:
                        self.statusBar().showMessage("MOL file import not implemented for 2D editor.")
                event.acceptProposedAction()
            elif file_path.lower().endswith('.xyz'):
                # XYZファイルは常に3Dビューアーで開く
                self.load_xyz_for_3d_viewing(file_path=file_path)
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
            'planar_xy_action', 
            'planar_xz_action',
            'planar_yz_action',
            'align_x_action',
            'align_y_action', 
            'align_z_action',
            'bond_length_action',
            'angle_action',
            'dihedral_action',
            'symmetrize_action',
            'mirror_action'
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
        # 基本的な3D機能
        basic_3d_actions = [
            'optimize_3d_button',
            'export_button', 
            'edit_3d_action',
            'analysis_action',
            'measurement_action'
        ]
        
        for action_name in basic_3d_actions:
            if hasattr(self, action_name):
                getattr(self, action_name).setEnabled(enabled)
        
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


    def apply_3d_settings(self):
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

    def load_settings(self):
        default_settings = {
            'background_color': '#919191',
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
            # Multiple bond offset parameters
            'double_bond_offset_factor': 2.0,
            'triple_bond_offset_factor': 2.0,
            'double_bond_radius_factor': 0.8,
            'triple_bond_radius_factor': 0.7,
        }

        try:
            if os.path.exists(self.settings_file):
                with open(self.settings_file, 'r') as f:
                    loaded_settings = json.load(f)
                
                for key, value in default_settings.items():
                    loaded_settings.setdefault(key, value)
                self.settings = loaded_settings
            
            else:
                self.settings = default_settings
        
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
            # 測定モードをオンにする時は、3D Editモードを無効化
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
        
        # 4つ以上選択された場合はクリア
        if len(self.selected_atoms_for_measurement) > 4:
            self.clear_measurement_selection()
            self.selected_atoms_for_measurement.append(atom_idx)
        
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
                    if label_item.scene():
                        self.scene.removeItem(label_item)
                except RuntimeError:
                    # オブジェクトが削除されている場合はスキップ
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

    # --- 3D Edit functionality ---
    
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
    
    def planarize_selection(self, plane):
        """選択された原子群を指定された平面に平面化する"""
        if not self.selected_atoms_3d or not self.current_mol:
            self.statusBar().showMessage("No atoms selected for planarization.")
            return
        
        if len(self.selected_atoms_3d) < 3:
            self.statusBar().showMessage("Please select at least 3 atoms for planarization.")
            return
        
        try:
            # 選択された原子の位置を取得
            selected_indices = list(self.selected_atoms_3d)
            selected_positions = self.atom_positions_3d[selected_indices].copy()
            
            # 重心を計算
            centroid = np.mean(selected_positions, axis=0)
            
            # 重心を原点に移動
            centered_positions = selected_positions - centroid
            
            if plane == 'xy':
                # XY平面に平面化（Z座標を0にする）
                centered_positions[:, 2] = 0
                plane_name = "XY"
            elif plane == 'xz':
                # XZ平面に平面化（Y座標を0にする）
                centered_positions[:, 1] = 0
                plane_name = "XZ"
            elif plane == 'yz':
                # YZ平面に平面化（X座標を0にする）
                centered_positions[:, 0] = 0
                plane_name = "YZ"
            else:
                self.statusBar().showMessage("Invalid plane specified.")
                return
            
            # 重心の位置に戻す
            new_positions = centered_positions + centroid
            
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
    
    def open_planarization_dialog(self, plane):
        """平面化ダイアログを開く"""
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
        
        dialog = PlanarizationDialog(self.current_mol, self, plane, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage(f"Atoms planarized to {plane.upper()} plane."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    
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
    
    def open_symmetrize_dialog(self):
        """対称化ダイアログを開く"""
        # Under Development メッセージを表示
        from PyQt6.QtWidgets import QMessageBox
        
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setWindowTitle("Symmetrize Function")
        msg.setText("Symmetrize Function - Under Development")
        msg.setInformativeText(
            "This function is under development and will be available in a future release."
        )
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
        
        # 元のコードは残しておく（コメントアウト）
        # # 測定モードを無効化
        # if self.measurement_mode:
        #     self.measurement_action.setChecked(False)
        #     self.toggle_measurement_mode(False)
        # 
        # dialog = SymmetrizeDialog(self.current_mol, self)
        # dialog.exec()  # モーダルダイアログとして表示
        # # 結果はダイアログ内で直接適用される

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
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
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
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
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
        
        self.accept()
    
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
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
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
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
        elif self.atom3_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - ?")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            self.selection_label.setText(f"Angle: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - {symbol3}({self.atom3_idx})")
            
            # Calculate current angle
            current_angle = self.calculate_angle()
            self.angle_label.setText(f"Current angle: {current_angle:.2f}°")
            self.apply_button.setEnabled(True)
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
            if new_angle <= 0 or new_angle >= 180:
                QMessageBox.warning(self, "Invalid Input", "Angle must be between 0 and 180 degrees.")
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
        
        self.accept()
    
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
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
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
        
        # Undo状態を保存
        self.main_window.push_undo_state()
        
        # Apply the dihedral angle change
        self.adjust_dihedral(new_dihedral)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()
        
        self.accept()
    
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


class SymmetryCandidatesDialog(QDialog):
    """対称性候補選択ダイアログ"""
    
    def __init__(self, candidates, parent=None):
        super().__init__(parent)
        self.candidates = candidates
        self.selected_candidate = None
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Symmetry Detection Results")
        self.setModal(True)
        self.setFixedSize(500, 400)
        layout = QVBoxLayout(self)
        
        # 説明ラベル
        instruction_label = QLabel(
            "Multiple point groups detected for your molecule. "
            "Select the most appropriate one based on confidence and chemical knowledge:"
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # 候補リスト
        self.candidates_list = QListWidget()
        self.candidates_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        
        for i, candidate in enumerate(self.candidates):
            # リストアイテムのテキストを作成
            confidence_bar = "★" * int(candidate['confidence'] * 5)  # 5段階評価
            item_text = f"{candidate['name']}\n"
            item_text += f"Confidence: {candidate['confidence']:.1%} {confidence_bar}\n"
            item_text += f"Reason: {candidate['reason']}"
            
            item = QListWidgetItem(item_text)
            item.setData(Qt.ItemDataRole.UserRole, candidate)
            
            # 信頼度に基づく色分け
            if candidate['confidence'] >= 0.8:
                item.setBackground(QColor(200, 255, 200))  # 淡い緑
            elif candidate['confidence'] >= 0.6:
                item.setBackground(QColor(255, 255, 200))  # 淡い黄色
            elif candidate['confidence'] >= 0.4:
                item.setBackground(QColor(255, 230, 200))  # 淡いオレンジ
            else:
                item.setBackground(QColor(255, 200, 200))  # 淡い赤
            
            self.candidates_list.addItem(item)
        
        # 最初の候補を選択
        if self.candidates:
            self.candidates_list.setCurrentRow(0)
        
        layout.addWidget(self.candidates_list)
        
        # 詳細情報表示エリア
        self.detail_label = QLabel()
        self.detail_label.setWordWrap(True)
        self.detail_label.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 10px; border: 1px solid #ccc; }")
        self.detail_label.setMaximumHeight(80)
        layout.addWidget(self.detail_label)
        
        # 選択が変更されたときの処理
        self.candidates_list.currentItemChanged.connect(self.on_selection_changed)
        
        # 初期表示
        self.on_selection_changed()
        
        # ボタン
        button_layout = QHBoxLayout()
        
        info_button = QPushButton("More Info")
        info_button.clicked.connect(self.show_more_info)
        button_layout.addWidget(info_button)
        
        button_layout.addStretch()
        
        self.select_button = QPushButton("Select This Point Group")
        self.select_button.clicked.connect(self.accept)
        self.select_button.setDefault(True)
        button_layout.addWidget(self.select_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.select_button.isEnabled():
                self.accept()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def on_selection_changed(self):
        """選択が変更されたときの処理"""
        current_item = self.candidates_list.currentItem()
        if current_item:
            candidate = current_item.data(Qt.ItemDataRole.UserRole)
            self.selected_candidate = candidate
            
            # 詳細情報を表示
            detail_text = f"Point Group: {candidate['name']}\n"
            
            # operations_countを安全に取得
            ops_count = candidate.get('operations_count', 'Unknown')
            if ops_count != 'Unknown':
                detail_text += f"Operations: {ops_count} symmetry operations\n"
            else:
                detail_text += "Operations: Computing...\n"
            
            # 化学的推奨情報
            recommendations = self.get_chemical_recommendations(candidate['key'])
            if recommendations:
                detail_text += f"Common molecules: {recommendations}"
            
            self.detail_label.setText(detail_text)
    
    def get_chemical_recommendations(self, point_group_key):
        """点群に対応する代表的な分子例を返す"""
        examples = {
            'Td': 'CH₄, CCl₄, SiF₄',
            'C3v': 'NH₃, PCl₃, SO₃²⁻',
            'C2v': 'H₂O, SO₂, NO₂⁻',
            'D3h': 'BF₃, CO₃²⁻, NO₃⁻',
            'D2h': 'C₂H₄, benzene (planar)',
            'Oh': 'SF₆, [Co(NH₃)₆]³⁺',
            'Cs': 'HOCl, planar molecules',
            'Ci': 'trans-compounds',
            'C2h': 'trans-N₂F₂',
            'C1': 'CHFClBr, most organic molecules'
        }
        return examples.get(point_group_key, '')
    
    def show_more_info(self):
        """詳細情報ダイアログを表示"""
        if not self.selected_candidate:
            return
        
        info_text = f"Point Group: {self.selected_candidate['name']}\n\n"
        info_text += f"Confidence Score: {self.selected_candidate['confidence']:.3f}\n"
        info_text += f"Detection Method: {self.selected_candidate['reason']}\n"
        info_text += f"Symmetry Operations: {self.selected_candidate.get('operations_count', 'N/A')}\n\n"
        
        examples = self.get_chemical_recommendations(self.selected_candidate['key'])
        if examples:
            info_text += f"Example Molecules:\n{examples}\n\n"
        
        info_text += "This point group describes the molecular symmetry and will be used to "
        info_text += "adjust atomic positions to achieve ideal symmetric geometry."
        
        QMessageBox.information(self, "Point Group Information", info_text)
    
    def get_selected_candidate(self):
        """選択された候補を返す"""
        return self.selected_candidate


# --- Application Execution ---
if __name__ == '__main__':
    main()

