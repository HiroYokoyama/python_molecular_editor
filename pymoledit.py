#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pymoledit.py — simple molecule editor
Short: Atom/bond edit, template preview, keyboard shortcuts (space, 1/2/3, Del etc.).
Author: HiroYokoyama
License: Apache-2.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
"""

import sys
import numpy as np
import pickle
import copy
import math

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QPushButton, QSplitter, QGraphicsView, QGraphicsScene, QGraphicsItem,
    QToolBar, QStatusBar, QGraphicsTextItem, QGraphicsLineItem, QDialog, QGridLayout,
    QFileDialog
)
from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QPainterPath, QFontMetrics, QKeySequence, QTransform, QCursor
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QLineF, QObject, QThread, pyqtSignal, QEvent

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# PyVista
import pyvista as pv
from pyvistaqt import QtInteractor

# --- Constants ---
ATOM_RADIUS = 18
BOND_OFFSET = 3.5
DEFAULT_BOND_LENGTH = 50.0 # テンプレートで使用する標準結合長
CPK_COLORS = {
    'H': QColor("white"), 'C': QColor("#333333"), 'N': QColor("#0000FF"),
    'O': QColor("#FF0000"), 'F': QColor("#00FF00"), 'S': QColor("#FFC000"),
    'Cl': QColor("#00FF00"), 'Br': QColor("#A52A2A"), 'I': QColor("#9400D3"),
    'P': QColor("#FFA500"), 'Si': QColor("#DAA520"), 'B': QColor("#FA8072"),
}
CPK_COLORS_PV = {
    'H': [0.9, 0.9, 0.9], 'C': [0.2, 0.2, 0.2], 'N': [0.0, 0.0, 1.0],
    'O': [1.0, 0.0, 0.0], 'F': [0.0, 1.0, 0.0], 'S': [1.0, 0.8, 0.0],
    'Cl': [0.0, 1.0, 0.0], 'Br': [0.6, 0.2, 0.0], 'I': [0.4, 0.0, 0.6],
    'P': [1.0, 0.65, 0.0], 'Si': [0.85, 0.65, 0.125], 'B': [0.98, 0.5, 0.45],
}
pt = Chem.GetPeriodicTable()
VDW_RADII = {pt.GetElementSymbol(i): pt.GetRvdw(i) * 0.3 for i in range(1, 119)}


# --- Data Model ---
class MolecularData:
    def __init__(self):
        self.atoms = {}; self.bonds = {}; self._next_atom_id = 0
    def add_atom(self, symbol, pos):
        atom_id = self._next_atom_id
        self.atoms[atom_id] = {'symbol': symbol, 'pos': pos, 'item': None}
        self._next_atom_id += 1; return atom_id
    def add_bond(self, id1, id2, order=1):
        if id1 > id2: id1, id2 = id2, id1
        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)]['order'] = order; return (id1, id2), 'updated'
        else:
            self.bonds[(id1, id2)] = {'order': order, 'item': None}; return (id1, id2), 'created'
    def remove_atom(self, atom_id):
        if atom_id in self.atoms:
            del self.atoms[atom_id]
            bonds_to_remove = [key for key in self.bonds if atom_id in key]
            for key in bonds_to_remove: del self.bonds[key]
    def remove_bond(self, id1, id2):
        if id1 > id2: id1, id2 = id2, id1
        if (id1, id2) in self.bonds: del self.bonds[(id1, id2)]
    def to_mol_block(self):
        if not self.atoms: return None
        atom_map = {old_id: new_id for new_id, old_id in enumerate(self.atoms.keys())}
        num_atoms, num_bonds = len(self.atoms), len(self.bonds)
        mol_block = "\n  PyQtEditor\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
        for old_id, atom in self.atoms.items():
            x, y, z, symbol = atom['item'].pos().x(), -atom['item'].pos().y(), 0.0, atom['symbol']
            mol_block += f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
        for (id1, id2), bond in self.bonds.items():
            idx1, idx2, order = atom_map[id1] + 1, atom_map[id2] + 1, bond['order']
            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}  0  0  0  0\n"
        mol_block += "M  END\n"
        return mol_block


# --- Custom 2D Graphics Items ---
class AtomItem(QGraphicsItem):
    def __init__(self, atom_id, symbol, pos):
        super().__init__()
        self.atom_id, self.symbol, self.bonds = atom_id, symbol, []
        self.setPos(pos)
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(1); self.font = QFont("Arial", 16, QFont.Weight.Bold); self.update_style()
    def boundingRect(self):
        # 原子の表示半径（クラス属性かグローバル定義を参照、なければ 10 を使う）
        r = getattr(self, 'radius', None)
        if r is None:
            r = globals().get('ATOM_RADIUS', 10)
        extra = 4.0  # 描画余白
        return QRectF(-r - extra, -r - extra, (r + extra) * 2.0, (r + extra) * 2.0)


    def shape(self):
        """
        原子ラベルの当たり判定を少し小さくする（見た目のラベルサイズより余裕を縮める）。
        ATOM_RADIUS をベースに小さめの円を返す（最低 4px）。
        """
        path = QPainterPath()
        hit_r = max(4.0, ATOM_RADIUS - 6.0)   # ここを調整すれば当たり判定の大きさを変えられます
        path.addEllipse(QRectF(-hit_r, -hit_r, hit_r * 2.0, hit_r * 2.0))
        return path

    def paint(self, painter, option, widget):
        if self.is_visible:
            painter.setFont(self.font); fm = painter.fontMetrics()
            text_rect = fm.boundingRect(self.symbol); text_rect.moveCenter(QPointF(0, 0).toPoint())
            if self.scene():
                bg_brush = self.scene().backgroundBrush(); bg_rect = text_rect.adjusted(-3, -3, 3, 3)
                painter.setBrush(bg_brush); painter.setPen(Qt.PenStyle.NoPen); painter.drawRect(bg_rect)
            if self.symbol == 'H':
                path = QPainterPath(); path.addText(text_rect.left(), text_rect.bottom(), self.font, self.symbol)
                painter.setPen(QPen(Qt.GlobalColor.black, 1)); painter.setBrush(QBrush(CPK_COLORS.get(self.symbol)))
                painter.drawPath(path)
            else:
                painter.setPen(QPen(CPK_COLORS.get(self.symbol, QColor("pink"))))
                painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, self.symbol)
        if self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush); painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.boundingRect())
    def update_style(self):
        self.is_visible = not (self.symbol == 'C' and len(self.bonds) > 0); self.update()
    def itemChange(self, change, value):
        res = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            for bond in self.bonds:
                bond.update_position()
        return res


class BondItem(QGraphicsItem):
    def __init__(self, atom1_item, atom2_item, order=1):
        super().__init__()
        self.atom1, self.atom2, self.order = atom1_item, atom2_item, order
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        # Set initial position and geometry
        self.update_position()

    def get_line_in_local_coords(self):
        """Returns the bond's line in its own local coordinates."""
        # The origin (0,0) of this item is atom1's position.
        # We need the position of atom2 relative to this item.
        p2 = self.mapFromItem(self.atom2, 0, 0)
        return QLineF(QPointF(0, 0), p2)

    def boundingRect(self):
        # 可能なら BondItem 側のユーティリティで線を取得、なければ atom の位置から作る
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            # get_line_in_local_coords が無ければ atom の位置から計算する（安全策）
            try:
                p1 = self.atom1.pos()
                p2 = self.atom2.pos()
                line = QLineF(p1, p2)
            except Exception:
                # 最終フォールバック：0 長さの線
                line = QLineF(0, 0, 0, 0)

        # 当たり判定余白（order によるオフセット + 10px の拡張）
        bond_offset = globals().get('BOND_OFFSET', 2)
        extra = (getattr(self, 'order', 1) - 1) * bond_offset + 10
        return QRectF(line.p1(), line.p2()).normalized().adjusted(-extra, -extra, extra, extra)


    def shape(self):
        path = QPainterPath()
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            try:
                p1 = self.atom1.pos()
                p2 = self.atom2.pos()
                line = QLineF(p1, p2)
            except Exception:
                return path

        if line.length() == 0:
            return path

        # 法線ベクトルを使って線の周りに幅を作る（ここで幅を約 10px に）
        normal = line.normalVector()
        offset = QPointF(normal.dx(), normal.dy()) * 10.0
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
        # This function is called when either atom moves.
        # We must always update the geometry and position.
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


# --- 2D Editor Scene Class ---
class MoleculeScene(QGraphicsScene):
    def __init__(self, data, window):
        super().__init__()
        self.data, self.window = data, window
        self.mode, self.current_atom_symbol, self.bond_order = 'select', 'C', 1
        self.start_atom, self.temp_line, self.start_pos = None, None, None; self.press_pos = None
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        
        # ホバー中の原子置換用マッピング
        self.key_to_symbol_map = {
            Qt.Key.Key_C: 'C', Qt.Key.Key_N: 'N', Qt.Key.Key_O: 'O', Qt.Key.Key_S: 'S',
            Qt.Key.Key_F: 'F', Qt.Key.Key_B: 'B', Qt.Key.Key_I: 'I', Qt.Key.Key_H: 'H',
        }
        self.key_to_symbol_map_shift = {
            Qt.Key.Key_C: 'Cl', Qt.Key.Key_B: 'Br',
        }

        self.reinitialize_items()

    def reinitialize_items(self):
        """シーンクリア後に再作成が必要なアイテムを初期化する"""
        self.template_preview = TemplatePreviewItem()
        self.addItem(self.template_preview)
        self.template_preview.hide()
        self.template_preview_points = []
        self.template_context = {}

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.initial_positions_in_event = {item: item.pos() for item in self.items() if isinstance(item, AtomItem)}


        if self.mode.startswith('template'):
            super().mousePressEvent(event)
            return
        
        item = self.itemAt(self.press_pos, self.views()[0].transform())
        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != 'select':
                self.temp_line=QGraphicsLineItem(QLineF(self.start_atom.pos(),self.press_pos)); self.temp_line.setPen(QPen(Qt.GlobalColor.red,2,Qt.PenStyle.DotLine)); self.addItem(self.temp_line)
            else: super().mousePressEvent(event)
        elif item is None and self.mode.startswith('atom'):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos))
            self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
            self.addItem(self.temp_line)
        else: super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.mode.startswith('template'):
            self.update_template_preview(event.scenePos())
            return

        if not self.mouse_moved_since_press and self.press_pos:
            if (event.scenePos() - self.press_pos).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True
        
        if self.temp_line:
            start_point = self.start_atom.pos() if self.start_atom else self.start_pos
            if start_point:
                self.temp_line.setLine(QLineF(start_point, event.scenePos()))
        else:
            super().mouseMoveEvent(event)


    def mouseReleaseEvent(self, event):
        end_pos = event.scenePos()
        
        is_click = self.press_pos and (end_pos - self.press_pos).manhattanLength() < QApplication.startDragDistance()


        if self.mode.startswith('template') and is_click:
            if self.template_context and self.template_context.get('points'):
                context = self.template_context
                self.add_molecule_fragment(
                    context['points'],
                    context['bonds_info'],
                    existing_items=context.get('items', [])
                )
                self.data_changed_in_event = True

        released_item = self.itemAt(end_pos, self.views()[0].transform())

        if self.mode.startswith('atom') and is_click and isinstance(released_item, BondItem):
            b=released_item; new_order=(b.order%3)+1; id1,id2=b.atom1.atom_id,b.atom2.atom_id
            if id1>id2: id1,id2=id2,id1
            self.data.bonds[(id1,id2)]['order']=new_order; b.order=new_order; b.update()
            self.data_changed_in_event = True
            if self.temp_line: self.removeItem(self.temp_line); self.temp_line=None
        elif self.start_atom and self.temp_line:
            self.removeItem(self.temp_line); self.temp_line = None
            line = QLineF(self.start_atom.pos(), end_pos); end_item = self.itemAt(end_pos, self.views()[0].transform())
            if line.length() < 10:
                if self.mode.startswith('atom'):
                    if self.start_atom.symbol != self.current_atom_symbol:
                        self.start_atom.symbol=self.current_atom_symbol; self.data.atoms[self.start_atom.atom_id]['symbol']=self.current_atom_symbol; self.start_atom.update_style()
                        self.data_changed_in_event = True
                elif self.mode == 'select':
                    pass
            else:
                if isinstance(end_item, AtomItem) and self.start_atom!=end_item: self.create_bond(self.start_atom, end_item)
                else:
                    new_id = self.create_atom(self.current_atom_symbol, end_pos); new_item = self.data.atoms[new_id]['item']
                    self.create_bond(self.start_atom, new_item)
                self.data_changed_in_event = True
        elif self.start_pos:
            if self.temp_line:
                self.removeItem(self.temp_line); self.temp_line = None
            line = QLineF(self.start_pos, end_pos)
            if line.length() < 10:
                self.create_atom(self.current_atom_symbol, end_pos)
                self.data_changed_in_event = True
            else:
                start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                end_id = self.create_atom(self.current_atom_symbol, end_pos)
                self.create_bond(self.data.atoms[start_id]['item'], self.data.atoms[end_id]['item'], bond_order=1)
                self.data_changed_in_event = True
        else:
            super().mouseReleaseEvent(event)

        # ★★★ 変更点 ★★★
        # itemChangeが全ての移動アイテムに発行されない問題に対処するため、
        # ここで明示的に移動した全原子と、関連する結合を更新する。
        moved_atoms = []
        for item, old_pos in self.initial_positions_in_event.items():
            if item.scene() and item.pos() != old_pos:
                moved_atoms.append(item)
        
        if moved_atoms:
            self.data_changed_in_event = True
            bonds_to_update = set()
            for atom in moved_atoms:
                # データモデル内の位置情報を更新
                self.data.atoms[atom.atom_id]['pos'] = atom.pos()
                # 更新が必要な結合を収集
                for bond in atom.bonds:
                    bonds_to_update.add(bond)
            
            # 収集したユニークな結合をすべて更新
            for bond in bonds_to_update:
                bond.update_position()

            # 移動完了後にビューポートを強制的に再描画し、表示の不整合を防ぐ
            if self.views():
                self.views()[0].viewport().update()
        
        self.start_atom=None; self.start_pos = None; self.press_pos = None
        self.template_context = {}

        if self.data_changed_in_event:
            self.window.push_undo_state()

    def create_atom(self, symbol, pos):
        atom_id = self.data.add_atom(symbol, pos); atom_item = AtomItem(atom_id, symbol, pos)
        self.data.atoms[atom_id]['item'] = atom_item; self.addItem(atom_item); return atom_id

    def create_bond(self, start_atom, end_atom, bond_order=None):
        order = bond_order if bond_order is not None else self.bond_order
        key, status = self.data.add_bond(start_atom.atom_id, end_atom.atom_id, order)
        if status == 'created':
            bond_item=BondItem(start_atom, end_atom, order); self.data.bonds[key]['item']=bond_item
            start_atom.bonds.append(bond_item); end_atom.bonds.append(bond_item); self.addItem(bond_item)
        else: 
            bond_item = self.data.bonds[key]['item']
            bond_item.order=order; bond_item.update()
        start_atom.update_style(); end_atom.update_style()

    def add_molecule_fragment(self, points, bonds_info, existing_items, symbol='C'):
        num_points = len(points)
        atom_items = [None] * num_points

        # Step 1: 既存のアイテム（基準となる原子や結合）を新しい原子リストに配置する
        if len(existing_items) == 2 and all(isinstance(i, AtomItem) for i in existing_items):
            # 基準が結合の場合
            atom_items[0] = existing_items[0]
            atom_items[1] = existing_items[1]
        elif len(existing_items) == 1 and isinstance(existing_items[0], AtomItem):
            # 基準が原子の場合
            atom_items[0] = existing_items[0]
        
        # Step 2: まだ配置されていない頂点に対して新しい原子を作成する
        for i, pos in enumerate(points):
            if atom_items[i] is None:
                atom_id = self.create_atom(symbol, pos)
                atom_items[i] = self.data.atoms[atom_id]['item']

        # Step 3: 原子間に結合を作成する
        for id1_idx, id2_idx, order in bonds_info:
            if atom_items[id1_idx] and atom_items[id2_idx]:
                self.create_bond(atom_items[id1_idx], atom_items[id2_idx], bond_order=order)


    def update_template_preview(self, pos):
        mode_parts = self.mode.split('_')
        is_aromatic = False
        if mode_parts[1] == 'benzene':
            n = 6
            is_aromatic = True
        else:
            try: n = int(mode_parts[1])
            except ValueError: return

        # --- 変更箇所：item の取得を「template_preview を無視して下位の Atom/Bond を探す」方式にする ---
        # NOTE: self.items(pos) は Z-order（上位→下位）のリストを返すので、
        # テンプレート自体 (TemplatePreviewItem) に当たっても、その下にある Atom/Bond を見つけられます。
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
            points = self._calculate_polygon_from_edge(p0, p1, n, cursor_pos=pos)
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
            # ビューポート再描画を強制して更新漏れを防ぐ
            if self.views():
                self.views()[0].viewport().update()
        else:
            self.template_preview.hide()
            if self.views():
                self.views()[0].viewport().update()


    def _calculate_polygon_from_edge(self, p0, p1, n, cursor_pos=None):
        if n < 3: return []
        v_edge = p1 - p0
        edge_length = math.sqrt(v_edge.x()**2 + v_edge.y()**2)
        if edge_length == 0: return []
        
        # --- ▼▼▼ ここから変更 ▼▼▼ ---
        # テンプレートのサイズが常に一定になるように、基準辺の長さを標準結合長に正規化します。
        # これにより、プレビューが小さすぎたり大きすぎたりして見えなくなる問題を防ぎます。
        v_edge = (v_edge / edge_length) * DEFAULT_BOND_LENGTH
        # 正規化後のベクトルに基づいて、2番目の頂点p1を再計算します。
        p1 = p0 + v_edge
        # --- ▲▲▲ ここまで変更 ▲▲▲ ---

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

    def leaveEvent(self, event):
        self.template_preview.hide()
        super().leaveEvent(event)

    def keyPressEvent(self, event):
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())

        # --- Atom hover: element change (既存) ---
        if isinstance(item_at_cursor, AtomItem):
            key = event.key()
            modifiers = event.modifiers()
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

        # --- 追加：選択モードで Space を押したら全選択 ---
        if event.key() == Qt.Key.Key_Space and self.mode == 'select':
            self.window.select_all()
            event.accept()
            return

        # --- 追加：ホバー／選択されている結合に対して 1/2/3 で結合次数を設定 ---
        if event.key() in (Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3):
            # numeric -> desired bond order
            order = {Qt.Key.Key_1: 1, Qt.Key.Key_2: 2, Qt.Key.Key_3: 3}[event.key()]

            target_bonds = []
            # 1) カーソル下の結合があればそれを優先
            if isinstance(item_at_cursor, BondItem):
                target_bonds = [item_at_cursor]
            else:
                # 2) 選択中の結合を対象
                target_bonds = [it for it in self.selectedItems() if isinstance(it, BondItem)]

            if target_bonds:
                for b in target_bonds:
                    # 更新：BondItem とデータモデル双方を更新
                    b.order = order
                    # データ側のキー (id1,id2) は小さい方を先にする
                    id1, id2 = b.atom1.atom_id, b.atom2.atom_id
                    if id1 > id2: id1, id2 = id2, id1
                    if (id1, id2) in self.data.bonds:
                        self.data.bonds[(id1, id2)]['order'] = order
                    b.update()
                self.window.push_undo_state()
                event.accept()
                return
            # ※ 該当する結合がなければ下へフォールバック（例：既存の '1' キーの原子追加処理）
        
        # --- ▼▼▼ 以下は既存の '1' キーで原子を追加する処理（フォールバック）--- 
        # （あなたの元コードの '1' による原子追加処理をそのまま残します）
        if event.key() == Qt.Key.Key_1:
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
                new_pos_offset = QPointF(0, -l) # デフォルトのオフセット

                # --- 60度にするためのロジック ---
                if len(start_atom.bonds) == 1:
                    bond = start_atom.bonds[0]
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2

                    existing_bond_vector = start_pos - other_atom.pos()
                    angle_rad = math.radians(60)
                    cos_a = math.cos(angle_rad)
                    sin_a = math.sin(angle_rad)
                    vx = existing_bond_vector.x()
                    vy = existing_bond_vector.y()
                    new_vx = vx * cos_a - vy * sin_a
                    new_vy = vx * sin_a + vy * cos_a
                    rotated_vector = QPointF(new_vx, new_vy)
                    line = QLineF(QPointF(0, 0), rotated_vector)
                    line.setLength(l)
                    new_pos_offset = line.p2()

                elif start_atom.bonds:
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

                new_atom_id = self.create_atom('C', start_pos + new_pos_offset)
                new_atom_item = self.data.atoms[new_atom_id]['item']
                self.create_bond(start_atom, new_atom_item, bond_order=1)

                self.clearSelection()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 既存の削除処理等はそのまま継続 ---
        if event.key() == Qt.Key.Key_Delete or event.key() == Qt.Key.Key_Backspace:
            selected_items = self.selectedItems()
            if not selected_items: return

            atoms_to_delete = {item for item in selected_items if isinstance(item, AtomItem)}
            bonds_to_delete = {item for item in selected_items if isinstance(item, BondItem)}

            # 削除対象の原子に接続している結合もすべて削除リストに追加
            for atom in atoms_to_delete:
                bonds_to_delete.update(atom.bonds)

            # 更新が必要な（削除されない）原子を特定
            atoms_to_update = set()
            for bond in bonds_to_delete:
                if bond.atom1 not in atoms_to_delete:
                    atoms_to_update.add(bond.atom1)
                if bond.atom2 not in atoms_to_delete:
                    atoms_to_update.add(bond.atom2)

            # --- 安全な削除処理 ---
            for bond in list(bonds_to_delete):
                if bond.atom1 in atoms_to_update and bond in bond.atom1.bonds:
                    bond.atom1.bonds.remove(bond)
                if bond.atom2 in atoms_to_update and bond in bond.atom2.bonds:
                    bond.atom2.bonds.remove(bond)

                # データモデルとシーンから結合を削除
                self.data.remove_bond(bond.atom1.atom_id, bond.atom2.atom_id)
                self.removeItem(bond)

            for atom in list(atoms_to_delete):
                # データモデルとシーンから原子を削除
                self.data.remove_atom(atom.atom_id)
                self.removeItem(atom)

            # 残った原子の表示スタイルを更新
            for atom in atoms_to_update:
                atom.update_style()

            # undo スタックに追加
            self.window.push_undo_state()

            # 表示の残り（描画のほんのわずかなゴミ）を消すために強制再描画
            if self.views():
                self.views()[0].viewport().update()
                # さらにシーン自体も update() しておく（安全措置）
                self.update()

        else:
            super().keyPressEvent(event)



# --- Worker Thread for Calculations ---
class CalculationWorker(QObject):
    finished=pyqtSignal(object); error=pyqtSignal(str)
    def run_calculation(self, mol_block):
        try:
            if not mol_block: raise ValueError("No atoms to convert.")
            mol=Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None: raise ValueError("Failed to create molecule from MOL block.")
            mol=Chem.AddHs(mol); AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.MMFFOptimizeMolecule(mol)
            self.finished.emit(mol)
        except Exception as e: self.error.emit(str(e))


# --- Periodic Table Dialog ---
class PeriodicTableDialog(QDialog):
    element_selected = pyqtSignal(str)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select an Element"); layout = QGridLayout(self); self.setLayout(layout)
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
            ('Ce',8,4), ('Pr',8,5), ('Nd',8,6), ('Pm',8,7), ('Sm',8,8), ('Eu',8,9), ('Gd',8,10), ('Tb',8,11),
            ('Dy',8,12), ('Ho',8,13), ('Er',8,14), ('Tm',8,15), ('Yb',8,16), ('Lu',8,17),
            # Actinides (separate row)
            ('Th',9,4), ('Pa',9,5), ('U',9,6), ('Np',9,7), ('Pu',9,8), ('Am',9,9), ('Cm',9,10), ('Bk',9,11),
            ('Cf',9,12), ('Es',9,13), ('Fm',9,14), ('Md',9,15), ('No',9,16), ('Lr',9,17),
        ]
        for symbol, row, col in elements:
            b=QPushButton(symbol); b.setFixedSize(40,40); b.clicked.connect(self.on_button_clicked); layout.addWidget(b, row, col)
    def on_button_clicked(self):
        b=self.sender(); self.element_selected.emit(b.text()); self.accept()


# --- Main Window ---
class MainWindow(QMainWindow):
    start_calculation = pyqtSignal(str)
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Python Molecular Editor by HY"); self.setGeometry(100, 100, 1600, 900)
        self.data = MolecularData(); self.current_mol = None
        self.undo_stack = []
        self.redo_stack = []
        self.mode_actions = {} 
        self.init_ui()
        self.init_worker_thread()
        self.reset_undo_stack()

    def init_ui(self):
        self.init_menu_bar()

        splitter=QSplitter(Qt.Orientation.Horizontal)
        self.setCentralWidget(splitter)

        left_pane=QWidget()
        left_layout=QVBoxLayout(left_pane)
        left_layout.setContentsMargins(0,0,0,0)

        self.scene=MoleculeScene(self.data,self)
        self.scene.setSceneRect(-800,-800,1600,1600)
        self.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.view_2d=QGraphicsView(self.scene)
        self.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        left_layout.addWidget(self.view_2d)

        self.cleanup_button=QPushButton("Optimize 2D")
        self.cleanup_button.clicked.connect(self.clean_up_2d_structure)
        left_layout.addWidget(self.cleanup_button)
        splitter.addWidget(left_pane)

        right_pane=QWidget()
        right_layout=QVBoxLayout(right_pane)
        self.plotter=QtInteractor(right_pane)
        right_layout.addWidget(self.plotter)

        # --- ▼▼▼ ここから追加 ▼▼▼ ---
        # 3D画面(plotter)にイベントフィルターを設定
        self.plotter.installEventFilter(self)
        # --- ▲▲▲ ここまで追加 ▲▲▲ ---

        self.convert_button=QPushButton("Convert to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        right_layout.addWidget(self.convert_button)
        splitter.addWidget(right_pane)
        splitter.setSizes([800, 800])

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        self.tool_group = QActionGroup(self)
        self.tool_group.setExclusive(True)

        actions_data = [
            ("Select", 'select', 'Space'), ("C", 'atom_C', 'c'), ("H", 'atom_H', 'h'), ("B", 'atom_B', 'b'),
            ("N", 'atom_N', 'n'), ("O", 'atom_O', 'o'), ("S", 'atom_S', 's'), ("F", 'atom_F', 'f'),
            ("Cl", 'atom_Cl', 'Shift+C'), ("Br", 'atom_Br', 'Shift+B'), ("I", 'atom_I', 'i'), 
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
                action.setChecked(True); self.set_mode('select')
        
        toolbar.addSeparator()
        toolbar.addWidget(QPushButton("Templates:"))
        
        templates = [("Benzene", "template_benzene")] + [(f"{i}-Ring", f"template_{i}") for i in range(3, 10)]
        for text, mode in templates:
            action = QAction(text, self, checkable=True)
            action.setToolTip(f"{text} Template")
            action.triggered.connect(lambda c, m=mode: self.set_mode(m))
            self.mode_actions[mode] = action
            toolbar.addAction(action)
            self.tool_group.addAction(action)
    # MainWindowクラスにこのメソッドを丸ごと追加してください
    def eventFilter(self, obj, event):
        """イベントフィルター：3Dビューがクリックされた際の処理"""
        if obj is self.plotter and event.type() == QEvent.Type.MouseButtonPress:
            # 3Dビューがクリックされたら、即座に2Dビューにフォーカスを戻す
            self.view_2d.setFocus()
        return super().eventFilter(obj, event)
    # --- ▲▲▲ ここまで追加 ▲▲▲ ---

    def keyPressEvent(self, event):
        key = event.key()
        modifiers = event.modifiers()

        if key == Qt.Key.Key_Space:
            if self.scene.mode == 'select':
                self.select_all()
                event.accept()
                return
            else:
                self.set_mode('select')
                if 'select' in self.mode_actions:
                    self.mode_actions['select'].setChecked(True)
                event.accept()
                return

        mode_map = {
            (Qt.Key.Key_Space, Qt.KeyboardModifier.NoModifier): 'select',
            (Qt.Key.Key_C, Qt.KeyboardModifier.NoModifier): 'atom_C',
            (Qt.Key.Key_N, Qt.KeyboardModifier.NoModifier): 'atom_N',
            (Qt.Key.Key_O, Qt.KeyboardModifier.NoModifier): 'atom_O',
            (Qt.Key.Key_S, Qt.KeyboardModifier.NoModifier): 'atom_S',
            (Qt.Key.Key_F, Qt.KeyboardModifier.NoModifier): 'atom_F',
            (Qt.Key.Key_H, Qt.KeyboardModifier.NoModifier): 'atom_H',
            (Qt.Key.Key_B, Qt.KeyboardModifier.NoModifier): 'atom_B',
            (Qt.Key.Key_I, Qt.KeyboardModifier.NoModifier): 'atom_I',
            (Qt.Key.Key_C, Qt.KeyboardModifier.ShiftModifier): 'atom_Cl',
            (Qt.Key.Key_B, Qt.KeyboardModifier.ShiftModifier): 'atom_Br',
        }

        target_mode = mode_map.get((key, modifiers))

        if target_mode:
            self.set_mode(target_mode)
            if target_mode in self.mode_actions:
                self.mode_actions[target_mode].setChecked(True)
            event.accept() 
        else:
            super().keyPressEvent(event)


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
        file_menu.addAction(save_raw_action)
        load_raw_action = QAction("Open Project...", self); load_raw_action.triggered.connect(self.load_raw_data)
        file_menu.addAction(load_raw_action)
        
        edit_menu = menu_bar.addMenu("&Edit")
        self.undo_action = QAction("Undo", self); self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.triggered.connect(self.undo); edit_menu.addAction(self.undo_action)
        
        self.redo_action = QAction("Redo", self); self.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.redo_action.triggered.connect(self.redo); edit_menu.addAction(self.redo_action)
        
        edit_menu.addSeparator()
        
        select_all_action = QAction("Select All", self); select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.select_all); edit_menu.addAction(select_all_action)
        
        clear_all_action = QAction("Clear All", self)
        clear_all_action.triggered.connect(self.clear_all); edit_menu.addAction(clear_all_action)

    def get_current_state(self):
        atoms = {atom_id: {'symbol': data['symbol'], 'pos': (data['item'].pos().x(), data['item'].pos().y())}
                 for atom_id, data in self.data.atoms.items()}
        bonds = {key: {'order': data['order']} for key, data in self.data.bonds.items()}
        state = {'atoms': atoms, 'bonds': bonds, '_next_atom_id': self.data._next_atom_id}
        
        if self.current_mol:
            state['mol_3d'] = self.current_mol.ToBinary()
            
        return state

    def set_state_from_data(self, state_data):
        self.clear_2d_editor(push_to_undo=False)
        
        loaded_data = copy.deepcopy(state_data)
        raw_atoms = loaded_data.get('atoms', {})
        raw_bonds = loaded_data.get('bonds', {})

        for atom_id_str, data in raw_atoms.items():
            atom_id = int(atom_id_str)
            pos = QPointF(data['pos'][0], data['pos'][1])
            atom_item = AtomItem(atom_id, data['symbol'], pos)
            self.data.atoms[atom_id] = {'symbol': data['symbol'], 'pos': pos, 'item': atom_item}
            self.scene.addItem(atom_item)
        
        self.data._next_atom_id = loaded_data.get('_next_atom_id', max(self.data.atoms.keys()) + 1 if self.data.atoms else 0)

        for key_tuple, data in raw_bonds.items():
            id1, id2 = int(key_tuple[0]), int(key_tuple[1])
            if id1 in self.data.atoms and id2 in self.data.atoms:
                atom1_item = self.data.atoms[id1]['item']
                atom2_item = self.data.atoms[id2]['item']
                self.scene.create_bond(atom1_item, atom2_item, bond_order=data['order'])

        for atom_data in self.data.atoms.values():
            if atom_data['item']:
                atom_data['item'].update_style()
        self.scene.update()

        if 'mol_3d' in loaded_data:
            try:
                self.current_mol = Chem.Mol(loaded_data['mol_3d'])
                self.draw_molecule_3d(self.current_mol)
            except Exception as e:
                self.statusBar().showMessage(f"Could not load 3D model from project: {e}")
                self.current_mol = None
        else:
            self.current_mol = None
            self.plotter.clear()

    def push_undo_state(self):
        current_state_for_comparison = {
            'atoms': {k: (v['symbol'], v['item'].pos().x(), v['item'].pos().y()) for k, v in self.data.atoms.items()},
            'bonds': self.data.bonds,
            '_next_atom_id': self.data._next_atom_id
        }
        
        last_state_for_comparison = None
        if self.undo_stack:
            last_atoms = self.undo_stack[-1].get('atoms', {})
            last_state_for_comparison = {
                'atoms': {k: (v['symbol'], v['pos'][0], v['pos'][1]) for k, v in last_atoms.items()},
                'bonds': self.undo_stack[-1].get('bonds', {}),
                '_next_atom_id': self.undo_stack[-1].get('_next_atom_id')
            }

        if not last_state_for_comparison or current_state_for_comparison != last_state_for_comparison:
            state = self.get_current_state()
            self.undo_stack.append(state)
            self.redo_stack.clear()
        
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

    def redo(self):
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self.set_state_from_data(state)
        self.update_undo_redo_actions()
        
    def update_undo_redo_actions(self):
        self.undo_action.setEnabled(len(self.undo_stack) > 1)
        self.redo_action.setEnabled(len(self.redo_stack) > 0)

    def select_all(self):
        for item in self.scene.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.setSelected(True)


    def clear_all(self):
        if not self.data.atoms and self.current_mol is None: return
        self.clear_2d_editor(push_to_undo=False)
        self.current_mol = None
        self.plotter.clear()
        self.reset_undo_stack()
        
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
            self.current_mol = None
            self.plotter.clear()
            
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)

            conf = mol.GetConformer(); SCALE_FACTOR = 50.0
            rdkit_idx_to_my_id = {}
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            center_x = sum(p.x for p in positions)/len(positions) if positions else 0
            center_y = sum(p.y for p in positions)/len(positions) if positions else 0
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i); pos = conf.GetAtomPosition(i)
                scene_x=(pos.x-center_x)*SCALE_FACTOR; scene_y=-(pos.y-center_y)*SCALE_FACTOR
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y))
                rdkit_idx_to_my_id[i] = atom_id
            for bond in mol.GetBonds():
                b_idx,e_idx,b_type=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(),bond.GetBondTypeAsDouble()
                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item,a2_item=self.data.atoms[a1_id]['item'],self.data.atoms[a2_id]['item']
                self.scene.create_bond(a1_item, a2_item, bond_order=int(b_type))
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

    def load_raw_data(self):
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
                mol_block = Chem.MolToMolBlock(self.current_mol)
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
                xyz_lines.append(f"Generated by Python Molecular Editor coded by HY. SMILES: {smiles}")
                for i in range(num_atoms):
                    pos=conf.GetAtomPosition(i); symbol=self.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
                with open(file_path,'w') as f: f.write("\n".join(xyz_lines))
                self.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")

    def init_worker_thread(self):
        self.thread=QThread();self.worker=CalculationWorker();self.worker.moveToThread(self.thread)
        self.start_calculation.connect(self.worker.run_calculation)
        self.worker.finished.connect(self.on_calculation_finished); self.worker.error.connect(self.on_calculation_error)
        self.thread.start()
    
    def open_periodic_table_dialog(self):
        dialog=PeriodicTableDialog(self); dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action=self.tool_group.checkedAction()
        if checked_action: self.tool_group.setExclusive(False); checked_action.setChecked(False); self.tool_group.setExclusive(True)
        dialog.exec()

    def set_atom_from_periodic_table(self, symbol): self.set_mode(f'atom_{symbol}')

    def set_mode(self, mode_str):
        self.scene.mode = mode_str
        
        if mode_str.startswith('template'):
            self.view_2d.setMouseTracking(True)
            self.view_2d.viewport().setMouseTracking(True)
        else:
            self.view_2d.setMouseTracking(False)
            self.view_2d.viewport().setMouseTracking(False)
            self.scene.template_preview.hide()

        if mode_str.startswith('atom'): 
            parts = mode_str.split('_')
            self.scene.current_atom_symbol = parts[1]
            self.statusBar().showMessage(f"Mode: Draw Atom ({parts[1]})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str.startswith('bond'): 
            parts = mode_str.split('_')
            self.scene.bond_order = int(parts[1])
            self.statusBar().showMessage(f"Mode: Draw Bond (Order: {parts[1]})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str.startswith('template'):
            name = mode_str.split('_')[1]
            self.statusBar().showMessage(f"Mode: {name.capitalize()} Template")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        else: # Select mode
            self.scene.bond_order = 1
            self.statusBar().showMessage("Mode: Select")
            self.view_2d.setDragMode(QGraphicsView.DragMode.RubberBandDrag)

    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No atoms to optimize."); return
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if mol is None: self.statusBar().showMessage("Error: Failed to create molecule for optimization."); return
        try:
            AllChem.Compute2DCoords(mol); conf=mol.GetConformer(); original_ids=list(self.data.atoms.keys()); SCALE=50.0
            positions=[conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if not positions: self.statusBar().showMessage("Optimization complete."); return
            cx=sum(p.x for p in positions)/len(positions); cy=sum(p.y for p in positions)/len(positions)
            for i, atom_id in enumerate(original_ids):
                if i < mol.GetNumAtoms() and atom_id in self.data.atoms:
                    item=self.data.atoms[atom_id]['item']
                    new_pos=conf.GetAtomPosition(i)
                    sx=(new_pos.x-cx)*SCALE; sy=-(new_pos.y-cy)*SCALE; item.setPos(sx,sy)
            for bond_data in self.data.bonds.values():
                if bond_data.get('item'): bond_data['item'].update_position()
            self.statusBar().showMessage("2D structure optimization successful.")
            self.push_undo_state()
        except Exception as e: self.statusBar().showMessage(f"Error during 2D optimization: {e}")

    def trigger_conversion(self):
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No atoms to convert."); return
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
        if mol is None: self.statusBar().showMessage("Error: Invalid chemical structure."); return
        if len(Chem.GetMolFrags(mol)) > 1: self.statusBar().showMessage("Error: 3D conversion not supported for multiple molecules."); return
        self.convert_button.setEnabled(False)
        self.statusBar().showMessage("Calculating 3D structure...")
        self.start_calculation.emit(mol_block)
        
    def on_calculation_finished(self, mol):
        self.current_mol=mol
        self.draw_molecule_3d(mol)
        self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
        self.push_undo_state()
        
    def on_calculation_error(self, error_message):
        self.statusBar().showMessage(f"Error: {error_message}")
        self.convert_button.setEnabled(True)
        
    def draw_molecule_3d(self, mol):
        self.plotter.clear(); conf = mol.GetConformer()
        pos=np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        sym=[a.GetSymbol() for a in mol.GetAtoms()]
        rad=np.array([VDW_RADII.get(s,0.4) for s in sym]); col=np.array([CPK_COLORS_PV.get(s,[0.5,0.5,0.5]) for s in sym])
        poly=pv.PolyData(pos); poly['colors']=col; poly['radii']=rad
        glyphs=poly.glyph(scale='radii',geom=pv.Sphere(radius=1.0),orient=False)
        
        edge_color = '#505050'
        self.plotter.add_mesh(glyphs, scalars='colors', rgb=True, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3, line_width=0.1)
        
        for bond in mol.GetBonds():
            sp=np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx())); ep=np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
            bt=bond.GetBondType(); c=(sp+ep)/2; d=ep-sp; h=np.linalg.norm(d)
            if h==0: continue
            color=[0.5, 0.5, 0.5]
            if bt==Chem.rdchem.BondType.SINGLE or bt==Chem.rdchem.BondType.AROMATIC:
                cyl=pv.Cylinder(center=c,direction=d,radius=0.1,height=h)
                self.plotter.add_mesh(cyl, color=color, smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
            else:
                v1=d/h; v_arb=np.array([0,0,1])
                if np.allclose(np.abs(np.dot(v1,v_arb)),1.0): v_arb=np.array([0,1,0])
                off_dir=np.cross(v1,v_arb); off_dir/=np.linalg.norm(off_dir)
                if bt==Chem.rdchem.BondType.DOUBLE:
                    r=0.09; s=0.15; c1=c+off_dir*(s/2); c2=c-off_dir*(s/2)
                    cyl1=pv.Cylinder(center=c1,direction=d,radius=r,height=h); cyl2=pv.Cylinder(center=c2,direction=d,radius=r,height=h)
                    self.plotter.add_mesh(cyl1,color=color,smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                    self.plotter.add_mesh(cyl2,color=color,smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                elif bt==Chem.rdchem.BondType.TRIPLE:
                    r=0.08; s=0.18; cc=pv.Cylinder(center=c,direction=d,radius=r,height=h)
                    c1=pv.Cylinder(center=c+off_dir*s,direction=d,radius=r,height=h); c2=pv.Cylinder(center=c-off_dir*s,direction=d,radius=r,height=h)
                    self.plotter.add_mesh(cc,color=color,smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                    self.plotter.add_mesh(c1,color=color,smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
                    self.plotter.add_mesh(c2,color=color,smooth_shading=True, show_edges=True, edge_color=edge_color, edge_opacity=0.3)
        self.plotter.reset_camera()
        
    def closeEvent(self, event):
        self.thread.quit(); self.thread.wait(); super().closeEvent(event)

# --- Application Execution ---
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


