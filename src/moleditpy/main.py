#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
moleditpy — simple molecule editor
Short: Atom/bond edit, template preview, keyboard shortcuts).
Author: HiroYokoyama
License: Apache-2.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
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

from collections import deque

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QPushButton, QSplitter, QGraphicsView, QGraphicsScene, QGraphicsItem,
    QToolBar, QStatusBar, QGraphicsTextItem, QGraphicsLineItem, QDialog, QGridLayout,
    QFileDialog, QSizePolicy, QLabel, QLineEdit, QToolButton, QMenu, QMessageBox
)
from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QPainterPath, QFontMetrics, QKeySequence, QTransform, QCursor, QPixmap, QIcon, QShortcut
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QLineF, QObject, QThread, pyqtSignal, QEvent, QMimeData, QByteArray


# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

# PyVista
import pyvista as pv
from pyvistaqt import QtInteractor

# --- Constants ---
ATOM_RADIUS = 18
BOND_OFFSET = 3.5
DEFAULT_BOND_LENGTH = 75 # テンプレートで使用する標準結合長
CLIPBOARD_MIME_TYPE = "application/x-moleditpy-fragment"

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

def main():
    # --- Windows タスクバーアイコンのための追加処理 ---
    if sys.platform == 'win32':
        myappid = 'hyoko.moleditpy.0.2' # アプリケーション固有のID（任意）
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
        if id1 > id2: id1, id2 = id2, id1
        bond_data = {'order': order, 'stereo': stereo, 'item': None}
        
        is_new_bond = (id1, id2) not in self.bonds
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
            # Safely get neighbors before deleting the atom's own entry
            neighbors = self.adjacency_list.get(atom_id, [])
            for neighbor_id in neighbors:
                if neighbor_id in self.adjacency_list and atom_id in self.adjacency_list[neighbor_id]:
                    self.adjacency_list[neighbor_id].remove(atom_id)

            # Now, safely delete the atom's own entry from the adjacency list
            if atom_id in self.adjacency_list:
                del self.adjacency_list[atom_id]

            del self.atoms[atom_id]
            bonds_to_remove = [key for key in self.bonds if atom_id in key]
            for key in bonds_to_remove:
                del self.bonds[key]

    def remove_bond(self, id1, id2):
        if id1 > id2: id1, id2 = id2, id1
        if (id1, id2) in self.bonds:
            if id1 in self.adjacency_list and id2 in self.adjacency_list[id1]:
                self.adjacency_list[id1].remove(id2)
            if id2 in self.adjacency_list and id1 in self.adjacency_list[id2]:
                self.adjacency_list[id2].remove(id1)
            del self.bonds[(id1, id2)]

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
        mol_block = "\n  PyQtEditor\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
        for old_id, atom in self.atoms.items():
            x, y, z, symbol = atom['item'].pos().x(), -atom['item'].pos().y(), 0.0, atom['symbol']
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
            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}  0  0  0  0\n"
        mol_block += "M  END\n"
        return mol_block


    def to_rdkit_mol(self):
        if not self.atoms: return None
        mol = Chem.RWMol()

        for atom_id, atom_data in self.atoms.items():
            atom = Chem.Atom(atom_data['symbol'])
            atom.SetFormalCharge(atom_data.get('charge', 0))
            atom.SetNumRadicalElectrons(atom_data.get('radical', 0))
            atom.SetIntProp("_original_atom_id", atom_id)
            mol.AddAtom(atom)

        atom_id_to_idx_map = {a.GetIntProp("_original_atom_id"): a.GetIdx() for a in mol.GetAtoms()}

        for (id1, id2), bond_data in self.bonds.items():
            if id1 not in atom_id_to_idx_map or id2 not in atom_id_to_idx_map:
                continue
            idx1 = atom_id_to_idx_map[id1]
            idx2 = atom_id_to_idx_map[id2]

            order_val = float(bond_data['order'])
            if order_val == 1.5:
                order = Chem.BondType.AROMATIC
            elif order_val == 2.0:
                order = Chem.BondType.DOUBLE
            elif order_val == 3.0:
                order = Chem.BondType.TRIPLE
            else:
                order = Chem.BondType.SINGLE

            bond_idx = mol.AddBond(idx1, idx2, order)
            bond = mol.GetBondWithIdx(bond_idx - 1)

            if bond_data.get('stereo', 0) != 0 and order == Chem.BondType.SINGLE:
                if bond_data['stereo'] == 1: # Wedge
                    bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif bond_data['stereo'] == 2: # Dash
                    bond.SetBondDir(Chem.BondDir.BEGINDASH)

        mol = mol.GetMol()
        try:
            Chem.SanitizeMol(mol)
        except Exception:
             pass
        return mol

class AtomItem(QGraphicsItem):
    def __init__(self, atom_id, symbol, pos, charge=0, radical=0):
        super().__init__()
        self.atom_id, self.symbol, self.charge, self.radical, self.bonds, self.chiral_label = atom_id, symbol, charge, radical, [], None
        self.setPos(pos)
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(1); self.font = QFont("Arial", 20, QFont.Weight.Bold); self.update_style()
    
    def boundingRect(self):
        r = globals().get('ATOM_RADIUS', 18)
        extra = 18.0 # 元の値(15.0)から拡大し、電荷・ラジカル表示領域を確保
        return QRectF(-r - extra, -r - extra, (r + extra) * 2.0, (r + extra) * 2.0)

    def shape(self):
        path = QPainterPath()
        hit_r = max(4.0, ATOM_RADIUS - 6.0) * 2
        path.addEllipse(QRectF(-hit_r, -hit_r, hit_r * 2.0, hit_r * 2.0))
        return path

    def paint(self, painter, option, widget):
        color = CPK_COLORS.get(self.symbol, CPK_COLORS['DEFAULT'])
        if self.is_visible:
            # 1. 原子記号（C, Nなど）を描画するための準備
            painter.setFont(self.font)
            fm = painter.fontMetrics()
            text_rect = fm.boundingRect(self.symbol)
            text_rect.moveCenter(QPointF(0, 0).toPoint())
            
            # 2. 原子記号の背景を白で塗りつぶす（結合線と重なっても見やすくするため）
            if self.scene():
                bg_brush = self.scene().backgroundBrush()
                bg_rect = text_rect.adjusted(-3, -3, 3, 3)
                painter.setBrush(bg_brush)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.drawRect(bg_rect)
            
            # 3. 原子記号自体を描画
            # 水素だけ黒文字、それ以外は CPK カラー
            if self.symbol == 'H':
                painter.setPen(QPen(Qt.GlobalColor.black))
            else:
                painter.setPen(QPen(color))
            painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, self.symbol)
            
            
            # --- ここから電荷とラジカルの描画 ---
            # 電荷の描画
            if self.charge != 0:
                if self.charge == 1: charge_str = "+"
                elif self.charge == -1: charge_str = "-"
                else: charge_str = f"{self.charge:+}"
                charge_font = QFont("Arial", 12, QFont.Weight.Bold)
                painter.setFont(charge_font)
                charge_rect = painter.fontMetrics().boundingRect(charge_str)
                charge_pos = QPointF(text_rect.right() + 2, text_rect.top() + charge_rect.height() - 2)
                painter.setPen(Qt.GlobalColor.black)
                painter.drawText(charge_pos, charge_str)
            
            # ラジカルの描画
            if self.radical > 0:
                painter.setBrush(QBrush(Qt.GlobalColor.black))
                painter.setPen(Qt.PenStyle.NoPen)
                radical_pos_y = text_rect.top() - 5
                if self.radical == 1:
                    painter.drawEllipse(QPointF(text_rect.center().x(), radical_pos_y), 3, 3)
                elif self.radical == 2:
                    painter.drawEllipse(QPointF(text_rect.center().x() - 5, radical_pos_y), 3, 3)
                    painter.drawEllipse(QPointF(text_rect.center().x() + 5, radical_pos_y), 3, 3)

        # 選択時のハイライト
        if self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.boundingRect())

    def update_style(self):
        self.is_visible = not (self.symbol == 'C' and len(self.bonds) > 0 and self.charge == 0 and self.radical == 0)
        self.update()

    # 約203行目 AtomItem クラス内

    def itemChange(self, change, value):
        res = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            if self.flags() & QGraphicsItem.GraphicsItemFlag.ItemIsMovable:
                 for bond in self.bonds: bond.update_position()
            
        return res

class BondItem(QGraphicsItem):
    def __init__(self, atom1_item, atom2_item, order=1, stereo=0):
        super().__init__()
        self.atom1, self.atom2, self.order, self.stereo = atom1_item, atom2_item, order, stereo
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        self.update_position()

    def get_line_in_local_coords(self):
        p2 = self.mapFromItem(self.atom2, 0, 0)
        return QLineF(QPointF(0, 0), p2)

    def boundingRect(self):
        try: line = self.get_line_in_local_coords()
        except Exception: line = QLineF(0, 0, 0, 0)
        bond_offset = globals().get('BOND_OFFSET', 2)
        extra = (getattr(self, 'order', 1) - 1) * bond_offset + 15 # extraを拡大
        return QRectF(line.p1(), line.p2()).normalized().adjusted(-extra, -extra, extra, extra)

    def shape(self):
        path = QPainterPath()
        try: line = self.get_line_in_local_coords()
        except Exception: return path
        if line.length() == 0: return path

        normal = line.normalVector()
        offset = QPointF(normal.dx(), normal.dy()) * 12.0 # 当たり判定を拡大
        poly = QPolygonF([line.p1() - offset, line.p1() + offset, line.p2() + offset, line.p2() - offset])
        path.addPolygon(poly)
        return path

    def paint(self, painter, option, widget):
        line = self.get_line_in_local_coords()
        if line.length() == 0: return

        current_pen = self.pen
        if self.isSelected():
            current_pen = QPen(QColor("blue"), 3)
        painter.setPen(current_pen)
        painter.setBrush(QBrush(Qt.GlobalColor.black))

        # --- Stereo (Wedge/Dash) Drawing ---
        if self.order == 1 and self.stereo in [1, 2]:
            vec = line.unitVector()
            normal = vec.normalVector()
            
            p1 = line.p1() + vec.p2() * 5 # 少し原子から離す
            p2 = line.p2() - vec.p2() * 5

            if self.stereo == 1: # Wedge
                offset = QPointF(normal.dx(), normal.dy()) * 6.0
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)
            
            elif self.stereo == 2: # Dash
                num_dashes = 8
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = 5.0 * (1 - t)
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.drawLine(start_pt - offset, start_pt + offset)
            return

        # --- Default (Multi-order) Drawing ---
        if self.order == 1:
            painter.drawLine(line)
        else:
            v = line.unitVector().normalVector()
            offset = QPointF(v.dx(), v.dy()) * BOND_OFFSET
            if self.order == 2:
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))
            elif self.order == 3:
                painter.drawLine(line)
                painter.drawLine(line.translated(offset))
                painter.drawLine(line.translated(-offset))

    def update_position(self):
        self.prepareGeometryChange()
        if self.atom1:
            self.setPos(self.atom1.pos())
        self.update()


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

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.initial_positions_in_event = {item: item.pos() for item in self.items() if isinstance(item, AtomItem)}

        if self.mode.startswith('template'):
            event.accept()
            return
        
        item = self.itemAt(self.press_pos, self.views()[0].transform())
        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != 'select':
                self.temp_line = QGraphicsLineItem(QLineF(self.start_atom.pos(), self.press_pos))
                self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
                self.addItem(self.temp_line)
            else: super().mousePressEvent(event)
        elif item is None and self.mode.startswith('atom'):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos)); self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine)); self.addItem(self.temp_line)
        else: super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.mode.startswith('template'):
            self.update_template_preview(event.scenePos()); return

        if not self.mouse_moved_since_press and self.press_pos:
            if (event.scenePos() - self.press_pos).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True
        
        if self.temp_line:
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
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
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

        released_item = self.itemAt(end_pos, self.views()[0].transform())
        
        # 1. 特殊モード（ラジカル/電荷）の処理
        if (self.mode == 'radical') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            # ラジカルの状態をトグル (0 -> 1 -> 2 -> 0)
            atom.radical = (atom.radical + 1) % 3 
            self.data.atoms[atom.atom_id]['radical'] = atom.radical
            atom.update_style()
            self.data_changed_in_event = True
        elif (self.mode == 'charge_plus' or self.mode == 'charge_minus') and is_click and isinstance(released_item, AtomItem):
            atom = released_item
            delta = 1 if self.mode == 'charge_plus' else -1
            atom.charge += delta
            self.data.atoms[atom.atom_id]['charge'] = atom.charge
            atom.update_style()
            self.data_changed_in_event = True

        # 2. 既存結合の様式更新処理 (全ての結合モードで動作)
        elif self.mode.startswith('bond') and is_click and isinstance(released_item, BondItem):
            b = released_item
            # 既存結合をクリックした場合は、現在のモードの結合様式に更新する
            b.order = self.bond_order
            b.stereo = self.bond_stereo
            id1, id2 = b.atom1.atom_id, b.atom2.atom_id
            if id1 > id2: id1, id2 = id2, id1
            self.data.bonds[(id1,id2)]['order'] = self.bond_order
            self.data.bonds[(id1,id2)]['stereo'] = self.bond_stereo
            b.update()
            self.data_changed_in_event = True

        # 3. 新規原子・結合の作成処理 (atom_* モード および すべての bond_* モードで許可)
        elif self.start_atom and (self.mode.startswith('atom') or self.mode.startswith('bond')):
            line = QLineF(self.start_atom.pos(), end_pos); end_item = self.itemAt(end_pos, self.views()[0].transform())

            # 使用する結合様式を決定
            # atomモードの場合は bond_order/stereo を None にして create_bond にデフォルト値(1, 0)を適用
            # bond_* モードの場合は現在の設定 (self.bond_order/stereo) を使用
            order_to_use = self.bond_order if self.mode.startswith('bond') else None # ★変更点★
            stereo_to_use = self.bond_stereo if self.mode.startswith('bond') else None # ★変更点★
    
            
            if line.length() < 10:
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


        moved_atoms = [item for item, old_pos in self.initial_positions_in_event.items() if item.scene() and item.pos() != old_pos]
        if moved_atoms:
            self.data_changed_in_event = True
            bonds_to_update = set()
            for atom in moved_atoms:
                self.data.atoms[atom.atom_id]['pos'] = atom.pos()
                bonds_to_update.update(atom.bonds)
            for bond in bonds_to_update: bond.update_position()
            if self.views(): self.views()[0].viewport().update()
        
        self.start_atom=None; self.start_pos = None; self.press_pos = None; self.temp_line = None
        self.template_context = {}
        if self.data_changed_in_event: self.window.push_undo_state()

    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self.data.add_atom(symbol, pos, charge=charge, radical=radical)
        atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
        self.data.atoms[atom_id]['item'] = atom_item; self.addItem(atom_item); return atom_id


    def create_bond(self, start_atom, end_atom, bond_order=None, bond_stereo=None):
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
            start_atom.bonds.append(bond_item)
            end_atom.bonds.append(bond_item)
            self.addItem(bond_item)
        
        start_atom.update_style()
        end_atom.update_style()

    def add_molecule_fragment(self, points, bonds_info, existing_items=None, symbol='C'):
        """
        add_molecule_fragment の最終確定版。
        - 既存の結合次数を変更しないポリシーを徹底（最重要）。
        - ベンゼン環テンプレートは、フューズされる既存結合の次数に基づき、
          「新規に作られる二重結合が2本になるように」回転を決定するロジックを適用（条件分岐あり）。
        """
        import math
    
        num_points = len(points)
        atom_items = [None] * num_points
    
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
    
        # --- 1) 既にクリックされた existing_items をテンプレート頂点にマップ (変更なし) ---
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
    
        # --- 2) シーン内既存原子を self.data.atoms から列挙してマップ (変更なし) ---
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
    
        # --- 3) 足りない頂点は新規作成 (変更なし) ---
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

                # --- ★フューズされた辺の数による条件分岐★ ---
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
    
        # --- 5) ボンド作成／更新 (変更なし) ---
        for id1_idx, id2_idx, order in template_bonds_to_use:
            if id1_idx < len(atom_items) and id2_idx < len(atom_items):
                a_item, b_item = atom_items[id1_idx], atom_items[id2_idx]
                if not a_item or not b_item or a_item is b_item: continue
    
                # 順序正規化
                id1, id2 = a_item.atom_id, b_item.atom_id
                if id1 > id2: id1, id2 = id2, id1
    
                exist_b = self.find_bond_between(a_item, b_item)
                if exist_b:
                    # ★最重要ポリシー（フューズ辺は次数を変更しない）★
                    continue 
                else:
                    # 新規ボンド作成
                    self.create_bond(a_item, b_item, bond_order=order)
    
        # --- 6) 表示更新（必要なら） (変更なし) ---
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
        # -------------------------------------------------------------

        points, bonds_info = [], []
        l = DEFAULT_BOND_LENGTH
        self.template_context = {}

        if isinstance(item, AtomItem):
            # 原子にスナップ
            p0 = item.pos()
            direction = QLineF(p0, pos).unitVector()
            p1 = p0 + direction.p2() * l if direction.length() > 0 else p0 + QPointF(l, 0)
            points = self._calculate_polygon_from_edge(p0, p1, n)
            self.template_context['items'] = [item]

        elif isinstance(item, BondItem):
            # 結合にスナップ
            p0, p1 = item.atom1.pos(), item.atom2.pos()
            points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=pos, use_existing_length=True)
            self.template_context['items'] = [item.atom1, item.atom2]

        else: # 空白領域に配置
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
            if cross_product_z > 0:
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
        if not items_to_delete:
            return False

        atoms_to_delete = {item for item in items_to_delete if isinstance(item, AtomItem)}
        bonds_to_delete = {item for item in items_to_delete if isinstance(item, BondItem)}

        # 削除対象の原子に接続している結合も、すべて削除対象に加える
        for atom in atoms_to_delete:
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
            self.data.remove_atom(atom.atom_id)
        # 次に、明示的に選択された結合（まだ残っているもの）を削除
        for bond in bonds_to_delete:
            if bond.atom1 and bond.atom2:
                self.data.remove_bond(bond.atom1.atom_id, bond.atom2.atom_id)
        
        # --- シーンからのグラフィックアイテム削除（必ず結合を先に）---
        for bond in bonds_to_delete:
            if bond.scene(): self.removeItem(bond)
        for atom in atoms_to_delete:
            if atom.scene(): self.removeItem(atom)

        # --- 生き残った原子の内部参照とスタイルを更新 ---
        for atom in atoms_to_update:
            atom.bonds = [b for b in atom.bonds if b not in bonds_to_delete]
            atom.update_style()
            
        return True

    def leaveEvent(self, event):
        self.template_preview.hide(); super().leaveEvent(event)

    def keyPressEvent(self, event):
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())
        key = event.key()
        modifiers = event.modifiers()
        
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
                item_at_cursor.symbol = new_symbol
                self.data.atoms[item_at_cursor.atom_id]['symbol'] = new_symbol
                item_at_cursor.update_style()
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
            needs_update = False
            for bond in target_bonds:
                id1, id2 = bond.atom1.atom_id, bond.atom2.atom_id
                if id1 > id2: id1, id2 = id2, id1

                if key == Qt.Key.Key_1 and (bond.order != 1 or bond.stereo != 0):
                    bond.order = 1; bond.stereo = 0; needs_update = True
                elif key == Qt.Key.Key_2 and bond.order != 2:
                    bond.order = 2; bond.stereo = 0; needs_update = True
                elif key == Qt.Key.Key_3 and bond.order != 3:
                    bond.order = 3; bond.stereo = 0; needs_update = True
                elif key == Qt.Key.Key_W and bond.stereo != 1: # Wedge
                    bond.order = 1; bond.stereo = 1; needs_update = True
                elif key == Qt.Key.Key_D and bond.stereo != 2: # Dash
                    bond.order = 1; bond.stereo = 2; needs_update = True

                if needs_update:
                    self.data.bonds[(id1, id2)]['order'] = bond.order
                    self.data.bonds[(id1, id2)]['stereo'] = bond.stereo
                    bond.update()

            if needs_update:
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
                    
                    SUM_TOLERANCE = 5.0 # 総和ベクトルのマンハッタン長がこの値以下の場合、ゼロとみなす
                    
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


                SNAP_DISTANCE = 14.0
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
        

    # In MoleculeScene class
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
    

class CalculationWorker(QObject):
    # --- 変更点1: 新しいシグナルを追加 ---
    status_update = pyqtSignal(str) 
    
    finished=pyqtSignal(object); error=pyqtSignal(str)
    def run_calculation(self, mol_block):
        try:
            if not mol_block: raise ValueError("No atoms to convert.")
            mol=Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None: raise ValueError("Failed to create molecule from MOL block.")
            
            mol = Chem.AddHs(mol)
            
            params = AllChem.ETKDGv2()
            params.randomSeed = 42
            conf_id = AllChem.EmbedMolecule(mol, params)
            
            if conf_id == -1:
                self.status_update.emit("Initial embedding failed, retrying with a random seed...")
                params.randomSeed = -1
                conf_id = AllChem.EmbedMolecule(mol, params)

            if conf_id == -1:
                self.status_update.emit("Random seed retry failed, attempting with random coordinates...")
                conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
            if conf_id == -1:
                raise ValueError("Failed to generate a 3D conformer after multiple attempts.")
            
            AllChem.MMFFOptimizeMolecule(mol)
            
            self.finished.emit(mol)
        except Exception as e: self.error.emit(str(e))

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
    def __init__(self, mol, parent=None):
        super().__init__(parent)
        self.mol = mol
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
            
            # 表示するプロパティを辞書にまとめる
            properties = {
                "SMILES:": smiles,
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

class MainWindow(QMainWindow):

    start_calculation = pyqtSignal(str)
    def __init__(self, initial_file=None):
        super().__init__()
        self.setWindowTitle("moleditpy -- Python Molecular Editor by HY"); self.setGeometry(100, 100, 1400, 800)
        self.data = MolecularData(); self.current_mol = None
        self.current_3d_style = 'ball_and_stick'
        self.show_chiral_labels = False
        self.undo_stack = []
        self.redo_stack = []
        self.mode_actions = {} 
        self.init_ui()
        self.init_worker_thread()
        self.reset_undo_stack()
        self.scene.selectionChanged.connect(self.update_edit_menu_actions)
        QApplication.clipboard().dataChanged.connect(self.update_edit_menu_actions)
        self.update_edit_menu_actions()

        if initial_file:
            self.load_raw_data(file_path=initial_file)

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

        splitter=QSplitter(Qt.Orientation.Horizontal)
        self.setCentralWidget(splitter)

        left_pane=QWidget()
        left_layout=QVBoxLayout(left_pane)

        self.scene=MoleculeScene(self.data,self)
        self.scene.setSceneRect(-4000,-4000,4000,4000)
        self.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.view_2d=QGraphicsView(self.scene)
        self.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.view_2d.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        left_layout.addWidget(self.view_2d)

        self.view_2d.scale(0.5, 0.5)

        self.cleanup_button=QPushButton("Optimize 2D")
        self.cleanup_button.clicked.connect(self.clean_up_2d_structure)
        left_layout.addWidget(self.cleanup_button)
        splitter.addWidget(left_pane)

        right_pane=QWidget()
        right_layout=QVBoxLayout(right_pane)
        self.plotter=QtInteractor(right_pane)
        self.plotter.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        right_layout.addWidget(self.plotter)

        self.plotter.installEventFilter(self)
        
        self.convert_button=QPushButton("Convert to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        right_layout.addWidget(self.convert_button)
        splitter.addWidget(right_pane)
        splitter.setSizes([600, 600])

        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        
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
                offset = QPointF(normal.dx(), normal.dy()) * 6.0
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)
            elif bond_type == 'dash':
                vec = line.unitVector()
                normal = vec.normalVector()

                num_dashes = 6
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = 5.0 * (1 - t)
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.setPen(QPen(Qt.GlobalColor.black, 1.5))
                    painter.drawLine(start_pt - offset, start_pt + offset)

            painter.end()
            return QIcon(pixmap)

        # --- 結合ボタンをツールバーに追加 ---
        bond_actions_data = [
            ("Single Bond", 'bond_1_0', '1', 'single'),
            ("Double Bond", 'bond_2_0', '2', 'double'),
            ("Triple Bond", 'bond_3_0', '3', 'triple'),
            ("Wedge Bond", 'bond_1_1', 'W', 'wedge'),
            ("Dash Bond", 'bond_1_2', 'D', 'dash'),
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

            if n in [8, 9]:
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

        quit_shortcut = QShortcut(QKeySequence("Ctrl+Q"), self)
        quit_shortcut.activated.connect(self.close)

        self.view_2d.setFocus()

    def init_menu_bar(self):
        menu_bar = self.menuBar()
        
        file_menu = menu_bar.addMenu("&File")
        load_mol_action = QAction("Open MOL/SDF...", self); load_mol_action.triggered.connect(self.load_mol_file)
        file_menu.addAction(load_mol_action)
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

        optimize_2d_action = QAction("Optimize 2D", self)
        optimize_2d_action.setShortcut(QKeySequence("Ctrl+K"))
        optimize_2d_action.triggered.connect(self.clean_up_2d_structure)
        edit_menu.addAction(optimize_2d_action)
        
        convert_3d_action = QAction("Convert to 3D", self)
        convert_3d_action.setShortcut(QKeySequence("Ctrl+L"))
        convert_3d_action.triggered.connect(self.trigger_conversion)
        edit_menu.addAction(convert_3d_action)

        edit_menu.addSeparator()
        
        select_all_action = QAction("Select All", self); select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.select_all); edit_menu.addAction(select_all_action)
        
        clear_all_action = QAction("Clear All", self)
        clear_all_action.setShortcut(QKeySequence("Ctrl+Shift+C"))
        clear_all_action.triggered.connect(self.clear_all); edit_menu.addAction(clear_all_action)

        view_menu = menu_bar.addMenu("&View")
        self.toggle_chiral_action = QAction("Show Chiral Labels", self, checkable=True)
        self.toggle_chiral_action.setChecked(self.show_chiral_labels)
        self.toggle_chiral_action.triggered.connect(self.toggle_chiral_labels_display)
        view_menu.addAction(self.toggle_chiral_action)

        analysis_menu = menu_bar.addMenu("&Analysis")
        self.analysis_action = QAction("Show Analysis...", self)
        self.analysis_action.triggered.connect(self.open_analysis_window)
        self.analysis_action.setEnabled(False)
        analysis_menu.addAction(self.analysis_action)
        
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
        
        if mode_str.startswith('template'):
            self.view_2d.setMouseTracking(True)
        else:
            self.view_2d.setMouseTracking(False)
            self.scene.template_preview.hide()

        # カーソル形状の設定
        if mode_str == 'select':
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)
        elif mode_str.startswith(('atom', 'bond', 'template')):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        elif mode_str.startswith(('charge', 'radical')):
            self.view_2d.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)

        if mode_str.startswith('atom'): 
            self.scene.current_atom_symbol = mode_str.split('_')[1]
            self.statusBar().showMessage(f"Mode: Draw Atom ({self.scene.current_atom_symbol})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
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
            return # スタイルが変更されていない場合は何もしない

        self.current_3d_style = style_name
        self.statusBar().showMessage(f"3D style set to: {style_name}")
        
        # 現在表示中の分子があれば、新しいスタイルで再描画する
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)

    def copy_selection(self):
        """選択された原子と結合をクリップボードにコピーする"""
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
                    'stereo': bond_data['stereo'],
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

    def cut_selection(self):
        """選択されたアイテムを切り取り（コピーしてから削除）"""
        selected_items = self.scene.selectedItems()
        if not selected_items:
            return
        
        # 最初にコピー処理を実行
        self.copy_selection()
        
        if self.scene.delete_items(set(selected_items)):
            self.push_undo_state()
            self.statusBar().showMessage("Cut selection.", 2000)

    def paste_from_clipboard(self):
        """クリップボードから分子フラグメントを貼り付け"""
        clipboard = QApplication.clipboard()
        mime_data = clipboard.mimeData()
        if not mime_data.hasFormat(CLIPBOARD_MIME_TYPE):
            return

        byte_array = mime_data.data(CLIPBOARD_MIME_TYPE)
        buffer = io.BytesIO(byte_array)
        try:
            fragment_data = pickle.load(buffer)
        except pickle.UnpicklingError:
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
                bond_stereo=bond_data.get('stereo', 0)
            )
        
        self.push_undo_state()
        self.statusBar().showMessage(f"Pasted {len(new_atoms)} atoms.", 2000)
        self.activate_select_mode()

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
        mol = self.data.to_rdkit_mol()
        if not mol or mol.GetNumAtoms() == 0:
            # 3Dビューと関連データをクリア
            self.plotter.clear()
            self.current_mol = None
            self.analysis_action.setEnabled(False)
            self.statusBar().showMessage("3D view cleared.")
            self.view_2d.setFocus() 
            return

        problems = Chem.DetectChemistryProblems(mol)
        if problems:
            self.statusBar().showMessage(f"Error: {len(problems)} chemistry problem(s) found.")
            self.scene.clearSelection()
            
            # 問題のある原子をハイライト
            for prob in problems:
                atom_idx = prob.GetAtomIdx()
                rdkit_atom = mol.GetAtomWithIdx(atom_idx)
                original_id = rdkit_atom.GetIntProp("_original_atom_id")
                if original_id in self.data.atoms:
                    item = self.data.atoms[original_id]['item']
                    item.setSelected(True)
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
        self.convert_button.setEnabled(False)
        self.statusBar().showMessage("Calculating 3D structure...")
        self.start_calculation.emit(mol_block)
        
        self.view_2d.setFocus()

    def on_calculation_finished(self, mol):
        self.current_mol = mol
        self.draw_molecule_3d(mol)
        # ここで最適化済みの current_mol を用いて R/S を再解析して表示を更新
        try:
            self.update_chiral_labels()
        except Exception:
            # 念のためエラーを握り潰して UI を壊さない
            pass
    
        self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
        self.analysis_action.setEnabled(True)
        self.push_undo_state()
        self.view_2d.setFocus()
        
    def on_calculation_error(self, error_message):
        self.statusBar().showMessage(f"Error: {error_message}")
        self.convert_button.setEnabled(True)
        self.analysis_action.setEnabled(False)
        self.view_2d.setFocus() 

    def eventFilter(self, obj, event):
        if obj is self.plotter and event.type() == QEvent.Type.MouseButtonPress:
            self.view_2d.setFocus()
        return super().eventFilter(obj, event)

    def get_current_state(self):
        atoms = {atom_id: {'symbol': data['symbol'],
                           'pos': (data['item'].pos().x(), data['item'].pos().y()),
                           'charge': data.get('charge', 0),
                           'radical': data.get('radical', 0)}  # <-- この行を追加
                 for atom_id, data in self.data.atoms.items()}
        bonds = {key: {'order': data['order'], 'stereo': data.get('stereo', 0)} for key, data in self.data.bonds.items()}
        state = {'atoms': atoms, 'bonds': bonds, '_next_atom_id': self.data._next_atom_id}
        
        if self.current_mol: state['mol_3d'] = self.current_mol.ToBinary()
            
        return state

    def set_state_from_data(self, state_data):
        self.clear_2d_editor(push_to_undo=False)
        
        loaded_data = copy.deepcopy(state_data)
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

        if 'mol_3d' in loaded_data:
            try:
                self.current_mol = Chem.Mol(loaded_data['mol_3d'])
                self.draw_molecule_3d(self.current_mol)
                self.analysis_action.setEnabled(True)
            except Exception as e:
                self.statusBar().showMessage(f"Could not load 3D model from project: {e}")
                self.current_mol = None; self.analysis_action.setEnabled(False)
        else:
            self.current_mol = None; self.plotter.clear(); self.analysis_action.setEnabled(False)

        self.update_chiral_labels()
        

    def push_undo_state(self):
        current_state_for_comparison = {
            'atoms': {k: (v['symbol'], v['item'].pos().x(), v['item'].pos().y(), v.get('charge', 0), v.get('radical', 0)) for k, v in self.data.atoms.items()},
            'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in self.data.bonds.items()},
            '_next_atom_id': self.data._next_atom_id,
            'has_3d': self.current_mol is not None
        }
        
        last_state_for_comparison = None
        if self.undo_stack:
            last_atoms = self.undo_stack[-1].get('atoms', {})
            last_bonds = self.undo_stack[-1].get('bonds', {})
            last_state_for_comparison = {
                'atoms': {k: (v['symbol'], v['pos'][0], v['pos'][1], v.get('charge', 0), v.get('radical', 0)) for k, v in last_atoms.items()},
                'bonds': {k: (v['order'], v.get('stereo', 0)) for k, v in last_bonds.items()},
                '_next_atom_id': self.undo_stack[-1].get('_next_atom_id'),
                'has_3d': 'mol_3d' in self.undo_stack[-1] 
            }

        if not last_state_for_comparison or current_state_for_comparison != last_state_for_comparison:
            state = self.get_current_state()
            self.undo_stack.append(state)
            self.redo_stack.clear()

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
        self.update_undo_redo_actions()
        self.update_realtime_info()
        self.view_2d.setFocus() 

    def redo(self):
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self.set_state_from_data(state)
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
                mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                num_atoms = mol.GetNumAtoms()
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
        # データが存在しない場合は何もしない
        if not self.data.atoms and self.current_mol is None:
            return
            
        # 2Dエディタをクリアする（Undoスタックにはプッシュしない）
        self.clear_2d_editor(push_to_undo=False)
        
        # 3Dモデルをクリアする
        self.current_mol = None
        self.plotter.clear()
        
        # 解析メニューを無効化する
        self.analysis_action.setEnabled(False)
        
        # Undo/Redoスタックをリセットする
        self.reset_undo_stack()
        
        # シーンとビューの明示的な更新
        self.scene.update()
        if self.view_2d:
            self.view_2d.viewport().update()
        
        # 3Dプロッターの再描画
        self.plotter.render()
        
        # アプリケーションのイベントループを強制的に処理し、画面の再描画を確実に行う
        QApplication.processEvents()
        
        self.statusBar().showMessage("Cleared all data.")
        
    def clear_2d_editor(self, push_to_undo=True):
        self.data = MolecularData()
        self.scene.data = self.data
        self.scene.clear()
        self.scene.reinitialize_items()
        if push_to_undo:
            self.push_undo_state()


    def load_mol_file(self):
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Open MOL/SDF File", "", "Chemical Files (*.mol *.sdf);;All Files (*)", options=options)
        if not file_path: return
        try:
            suppl = Chem.SDMolSupplier(file_path, removeHs=False)
            mol = next(suppl, None)
            if mol is None: raise ValueError("Failed to read molecule from file.")
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None; self.plotter.clear(); self.analysis_action.setEnabled(False)
            
            if mol.GetNumConformers() == 0: AllChem.Compute2DCoords(mol)

            conf = mol.GetConformer(); SCALE_FACTOR = 50.0
            rdkit_idx_to_my_id = {}
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            center_x = sum(p.x for p in positions)/len(positions) if positions else 0
            center_y = sum(p.y for p in positions)/len(positions) if positions else 0
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i); pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()
                scene_x=(pos.x-center_x)*SCALE_FACTOR; scene_y=-(pos.y-center_y)*SCALE_FACTOR
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge)
                rdkit_idx_to_my_id[i] = atom_id
            for bond in mol.GetBonds():
                # (以降のBond読み込み部分は変更なし)
                b_idx,e_idx=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble(); b_dir = bond.GetBondDir()
                stereo = 0
                if b_dir == Chem.BondDir.BEGINWEDGE: stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH: stereo = 2
                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item,a2_item=self.data.atoms[a1_id]['item'],self.data.atoms[a2_id]['item']

                self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo)

            self.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.reset_undo_stack()
        except Exception as e: self.statusBar().showMessage(f"Error loading file: {e}")

    def save_raw_data(self):
        if not self.data.atoms: self.statusBar().showMessage("Error: Nothing to save."); return
        save_data = self.get_current_state()
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
        if file_path:
            if not file_path.lower().endswith('.pmeraw'): file_path += '.pmeraw'
            try:
                with open(file_path, 'wb') as f: pickle.dump(save_data, f)
                self.statusBar().showMessage(f"Project saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving project file: {e}")

    def load_raw_data(self, file_path=None):
        if not file_path:
            options = QFileDialog.Option.DontUseNativeDialog
            file_path, _ = QFileDialog.getOpenFileName(self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
            if not file_path: return
        
        try:
            with open(file_path, 'rb') as f: loaded_data = pickle.load(f)
            self.set_state_from_data(loaded_data)
            self.statusBar().showMessage(f"Project loaded from {file_path}")
            self.reset_undo_stack()
        except Exception as e: self.statusBar().showMessage(f"Error loading project file: {e}")

    def save_as_mol(self):
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No 2D data to save."); return
        options=QFileDialog.Option.DontUseNativeDialog
        file_path,_=QFileDialog.getSaveFileName(self,"Save 2D MOL File","","MOL Files (*.mol);;All Files (*)",options=options)
        if file_path:
            if not file_path.lower().endswith('.mol'): file_path += '.mol'
            try:
                with open(file_path,'w') as f: f.write(mol_block)
                self.statusBar().showMessage(f"2D data saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")
            
    def save_3d_as_mol(self):
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(self, "Save 3D MOL File", "", "MOL Files (*.mol);;All Files (*)", options=options)
        if file_path:
            if not file_path.lower().endswith('.mol'):
                file_path += '.mol'
            try:
                mol_block = Chem.MolToMolBlock(self.current_mol, includeStereo=True)
                with open(file_path, 'w') as f:
                    f.write(mol_block)
                self.statusBar().showMessage(f"3D data saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving 3D MOL file: {e}")

    def save_as_xyz(self):
        if not self.current_mol: self.statusBar().showMessage("Error: Please generate a 3D structure first."); return
        options=QFileDialog.Option.DontUseNativeDialog
        file_path,_=QFileDialog.getSaveFileName(self,"Save 3D XYZ File","","XYZ Files (*.xyz);;All Files (*)",options=options)
        if file_path:
            if not file_path.lower().endswith('.xyz'): file_path += '.xyz'
            try:
                conf=self.current_mol.GetConformer(); num_atoms=self.current_mol.GetNumAtoms()
                xyz_lines=[str(num_atoms)]; smiles=Chem.MolToSmiles(Chem.RemoveHs(self.current_mol))
                xyz_lines.append(f"Generated by moleditpy. SMILES: {smiles}")
                for i in range(num_atoms):
                    pos=conf.GetAtomPosition(i); symbol=self.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
                with open(file_path,'w') as f: f.write("\n".join(xyz_lines))
                self.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")

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
            original_center_x = sum(item.pos().x() for item in target_atom_items) / len(target_atom_items)
            original_center_y = sum(item.pos().y() for item in target_atom_items) / len(target_atom_items)

            positions = list(new_positions_map.values())
            rdkit_cx = sum(p.x for p in positions) / len(positions)
            rdkit_cy = sum(p.y for p in positions) / len(positions)

            SCALE = 50.0

            # 新しい座標を適用
            for atom_id, rdkit_pos in new_positions_map.items():
                if atom_id in self.data.atoms:
                    item = self.data.atoms[atom_id]['item']
                    sx = ((rdkit_pos.x - rdkit_cx) * SCALE) + original_center_x
                    sy = (-(rdkit_pos.y - rdkit_cy) * SCALE) + original_center_y
                    new_scene_pos = QPointF(sx, sy)
                    item.setPos(new_scene_pos)
                    self.data.atoms[atom_id]['pos'] = new_scene_pos

            # 最終的な座標に基づき、全ての結合表示を一度に更新
            for bond_data in self.data.bonds.values():
                if bond_data.get('item'):
                    bond_data['item'].update_position()

            # 重なり解消ロジックを実行
            self. resolve_overlapping_groups()
            
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
        self.scene.update()
        self.push_undo_state()
        self.statusBar().showMessage("Resolved overlapping groups.", 2000)


    def draw_molecule_3d(self, mol):
        self.plotter.clear()
        conf = mol.GetConformer()
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        sym = [a.GetSymbol() for a in mol.GetAtoms()]
        col = np.array([CPK_COLORS_PV.get(s, [0.5, 0.5, 0.5]) for s in sym])
                
        # スタイルに応じて原子半径を決定
        if self.current_3d_style == 'cpk':
            # CPKモデル: Van der Waals半径をそのまま使う（見栄えのため少しだけ縮小なら0.8くらい）
            rad = np.array([pt.GetRvdw(pt.GetAtomicNumber(s)) * 1.0 for s in sym])
        else: # ball_and_stick (デフォルト)
            # Ball & Stickモデル: 元のコード通り半径を小さくする
            rad = np.array([VDW_RADII.get(s, 0.4) for s in sym])

        poly = pv.PolyData(pos)
        poly['colors'] = col
        poly['radii'] = rad
        glyphs = poly.glyph(scale='radii', geom=pv.Sphere(radius=1.0), orient=False)
        
        edge_color = '#505050'
        self.plotter.add_mesh(glyphs, scalars='colors', rgb=True, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3, line_width=0.1)
        
        # Ball & Stick スタイルの時のみ結合を描画
        if self.current_3d_style == 'ball_and_stick':
            for bond in mol.GetBonds():
                sp = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
                ep = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
                bt = bond.GetBondType()
                c = (sp + ep) / 2
                d = ep - sp
                h = np.linalg.norm(d)
                if h == 0: continue
                color = [0.5, 0.5, 0.5]
                if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                    cyl = pv.Cylinder(center=c, direction=d, radius=0.1, height=h)
                    self.plotter.add_mesh(cyl, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                else:
                    v1 = d / h
                    v_arb = np.array([0, 0, 1])
                    if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                    off_dir = np.cross(v1, v_arb)
                    off_dir /= np.linalg.norm(off_dir)
                    if bt == Chem.rdchem.BondType.DOUBLE:
                        r = 0.08
                        s = 0.20
                        c1 = c + off_dir * (s / 2)
                        c2 = c - off_dir * (s / 2)
                        cyl1 = pv.Cylinder(center=c1, direction=d, radius=r, height=h)
                        cyl2 = pv.Cylinder(center=c2, direction=d, radius=r, height=h)
                        self.plotter.add_mesh(cyl1, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                        self.plotter.add_mesh(cyl2, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        r = 0.08
                        s = 0.20
                        cc = pv.Cylinder(center=c, direction=d, radius=r, height=h)
                        c1 = pv.Cylinder(center=c + off_dir * s, direction=d, radius=r, height=h)
                        c2 = pv.Cylinder(center=c - off_dir * s, direction=d, radius=r, height=h)
                        self.plotter.add_mesh(cc, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                        self.plotter.add_mesh(c1, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                        self.plotter.add_mesh(c2, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)

        # chiral labels: 3D コンフォマーを持つ場合はここで R/S を表示
        if getattr(self, 'show_chiral_labels', False):
            try:
                # RDKit 側でキラル中心を取得（3D コンフォマーがあればそれを使う）
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if chiral_centers:
                    # ラベル描画用の点配列と文字列リストを作成
                    # pos は numpy array (N,3) なので、index で取れば座標が得られる
                    pts = []
                    labels = []
                    # 少しラベルが球と重ならないよう z 方向にオフセット（全体のサイズに合わせて調整）
                    z_off = 0
                    for idx, lbl in chiral_centers:
                        coord = pos[idx].copy()
                        coord[2] += z_off
                        pts.append(coord)
                        labels.append(lbl if lbl is not None else '?')

                    # 既存のラベルを削除し（存在すれば）、再描画する
                    try:
                        # actor 名 'chiral_labels' をキーとして remove を試みる（無ければ ignore）
                        self.plotter.remove_actor('chiral_labels')
                    except Exception:
                        pass

                    # PyVista のラベル描画（font_size 等はお好みで調整）
                    self.plotter.add_point_labels(
                        np.array(pts),
                        labels,
                        font_size=20,
                        point_size=0,
                        text_color='k',
                        name='chiral_labels',
                        always_visible=True,   # ← ここを追加（ラベルの可視性フィルタをスキップ）
                        tolerance=0.01,        # ← 必要なら少し大きめにして消えにくくする（デフォルトは 0.001）
                        show_points=False      # ← ポイント自体の小さなドット表示が不要ならオフに（任意）
                    )
            except Exception as e: 
                self.statusBar().showMessage(f"3D chiral label drawing error: {e}")


    def toggle_chiral_labels_display(self, checked):
        """Viewメニューのアクションに応じてキラルラベル表示を切り替える"""
        self.show_chiral_labels = checked
        if checked:
            draw_molecule_3d()
            self.statusBar().showMessage("Chiral labels: will be (re)computed after Convert→3D.")
        else:
            draw_molecule_3d()
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
            Chem.AssignStereochemistry(mol_for_chirality, cleanIt=True, force=True, flagPossibleStereoCenters=True)

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


    def open_analysis_window(self):
        if self.current_mol:
            dialog = AnalysisWindow(self.current_mol, self)
            dialog.exec()
        else:
            self.statusBar().showMessage("Please generate a 3D structure first to show analysis.")

    def closeEvent(self, event):
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

# --- Application Execution ---
if __name__ == '__main__':
    main()
