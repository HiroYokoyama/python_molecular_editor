#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import QGraphicsItem

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QFont, QPainterPath, QFontMetricsF
)

from PyQt6.QtCore import (
    Qt, QPointF, QRectF
)

try:
    from .constants import (
        ATOM_RADIUS, DESIRED_ATOM_PIXEL_RADIUS,
        FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD,
        CPK_COLORS,
    )
except Exception:
    from modules.constants import (
        ATOM_RADIUS, DESIRED_ATOM_PIXEL_RADIUS,
        FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD,
        CPK_COLORS,
    )
    
try:
    from . import sip_isdeleted_safe
except Exception:
    from modules import sip_isdeleted_safe

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
