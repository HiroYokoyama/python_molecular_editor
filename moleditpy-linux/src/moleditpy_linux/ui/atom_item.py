#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations
import logging  # [REPORT ERROR MISSING ATTRIBUTE]
from typing import Any, List, Optional

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontMetricsF,
    QPainter,
    QPainterPath,
    QPen,
)
from PyQt6.QtWidgets import QGraphicsItem, QWidget

try:
    from ..utils.constants import (
        ATOM_RADIUS,
        CPK_COLORS,
        DESIRED_ATOM_PIXEL_RADIUS,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
    )
except ImportError:
    from moleditpy_linux.utils.constants import (
        ATOM_RADIUS,
        CPK_COLORS,
        DESIRED_ATOM_PIXEL_RADIUS,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
    )
from PyQt6 import sip


def sip_isdeleted_safe(obj: Any) -> bool:
    try:
        return sip.isdeleted(obj)
    except (AttributeError, TypeError, ImportError, RuntimeError):
        # If the object does not support sip.isdeleted or is already in a state
        # where the check fails, we assume it's unsafe or "deleted" for our purposes.
        return True


class AtomItem(QGraphicsItem):
    def __init__(
        self, atom_id: int, symbol: str, pos: QPointF, charge: int = 0, radical: int = 0
    ) -> None:
        super().__init__()
        self.atom_id: int = atom_id
        self.symbol: str = symbol
        self.charge: int = charge
        self.radical: int = radical
        self.bonds: List[Any] = []
        self.chiral_label: Optional[str] = None

        self.setPos(pos)
        self.implicit_h_count: int = 0
        self.setFlags(
            QGraphicsItem.GraphicsItemFlag.ItemIsMovable
            | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
        )
        self.setZValue(1)
        self.update_style()
        self.setAcceptHoverEvents(True)
        self.hovered: bool = False
        self.has_problem: bool = False
        self.is_visible: bool = True
        self.font: QFont = QFont(FONT_FAMILY, 20, FONT_WEIGHT_BOLD)

    def update_style(self) -> None:
        if sip_isdeleted_safe(self):
            return
        # Allow updating font preference dynamically
        font_size = 20
        font_family = FONT_FAMILY

        scene = self.scene()
        if scene is not None:
            if hasattr(scene, "get_setting"):
                font_size = scene.get_setting("atom_font_size_2d", 20)
                font_family = scene.get_setting("atom_font_family_2d", FONT_FAMILY)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                )

        self.font = QFont(font_family, font_size, FONT_WEIGHT_BOLD)
        self.prepareGeometryChange()

        self.is_visible = not (
            self.symbol == "C"
            and len(self.bonds) > 0
            and self.charge == 0
            and self.radical == 0
        )
        self.update()

    def boundingRect(self) -> QRectF:
        """Calculate the bounding rectangle for the atom item."""
        # --- Calculate text position and size using logic matching paint() ---
        # Get dynamic font size and family
        font_size = 20
        font_family = FONT_FAMILY
        scene = self.scene()
        if scene is not None:
            if hasattr(scene, "get_setting"):
                font_size = scene.get_setting("atom_font_size_2d", 20)
                font_family = scene.get_setting("atom_font_family_2d", FONT_FAMILY)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                )

        font = QFont(font_family, font_size, FONT_WEIGHT_BOLD)
        fm = QFontMetricsF(font)

        hydrogen_part = ""
        if self.implicit_h_count > 0:
            is_skeletal_carbon = (
                self.symbol == "C"
                and self.charge == 0
                and self.radical == 0
                and len(self.bonds) > 0
            )
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
                except (AttributeError, RuntimeError, TypeError):
                    # Skip any bond that raises while inspecting; keep UI tolerant.
                    # This happens if the underlying C++ object is being destroyed.
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

        # 1. Calculate the background rectangle (bg_rect) used in paint()
        bg_rect = text_rect.adjusted(-5, -8, 5, 8)

        # 2. Construct the full visual rectangle relative to bg_rect
        full_visual_rect = QRectF(bg_rect)

        # Include charge symbol area in calculation
        if self.charge != 0:
            # Chemical convention: single charge as "+"/"-", multiple as "2+"/"2-"
            if self.charge == 1:
                charge_str = "+"
            elif self.charge == -1:
                charge_str = "-"
            else:
                sign = "+" if self.charge > 0 else "-"
                charge_str = f"{abs(self.charge)}{sign}"
            charge_font = QFont("Arial", 12, QFont.Weight.Bold)
            charge_fm = QFontMetricsF(charge_font)
            charge_rect = charge_fm.boundingRect(charge_str)

            if flip_text:
                charge_pos = QPointF(
                    text_rect.left() - charge_rect.width() - 2, text_rect.top()
                )
            else:
                charge_pos = QPointF(text_rect.right() + 2, text_rect.top())
            charge_rect.moveTopLeft(charge_pos)
            full_visual_rect = full_visual_rect.united(charge_rect)

        # Include radical symbol area in calculation
        if self.radical > 0:
            radical_area = QRectF(
                text_rect.center().x() - 8, text_rect.top() - 8, 16, 8
            )
            full_visual_rect = full_visual_rect.united(radical_area)

        # 3. Add final margins for selection highlights, etc.
        return full_visual_rect.adjusted(-3, -3, 3, 3)

    def shape(self) -> QPainterPath:
        """Define the shape of the atom item for collision detection."""
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

    def paint(
        self,
        painter: Optional[QPainter],
        option: Any,
        widget: Optional[QWidget] = None,
    ) -> None:
        """Paint the atom symbol and its associated labels (charge, radical)."""
        if painter is None:
            return
        # Color logic: check if we should use bond color (uniform) or CPK (element-specific)
        color = CPK_COLORS.get(self.symbol, CPK_COLORS["DEFAULT"])
        # Use bond color if specified in settings
        scene = self.scene()
        if hasattr(scene, "get_setting") and (
            self.symbol == "H" or scene.get_setting("atom_use_bond_color_2d", False)
        ):
            custom_color = scene.get_setting("bond_color_2d", "#222222")
            if isinstance(custom_color, str):
                color = QColor(custom_color)

        if self.is_visible:
            # 1. Preparation for painting
            # Ensure correct font is used (self.font should be updated by update_style)
            painter.setFont(self.font)
            fm = painter.fontMetrics()

            # --- Create text for the hydrogen part ---
            hydrogen_part = ""
            if self.implicit_h_count > 0:
                is_skeletal_carbon = (
                    self.symbol == "C"
                    and self.charge == 0
                    and self.radical == 0
                    and len(self.bonds) > 0
                )
                if not is_skeletal_carbon:
                    hydrogen_part = "H"
                    if self.implicit_h_count > 1:
                        subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                        hydrogen_part += str(self.implicit_h_count).translate(
                            subscript_map
                        )

            # --- Determine if the text should be flipped ---
            flip_text = False
            # Consider flipping only if H-label exists and there are one or more bonds
            if hydrogen_part and self.bonds:
                # Determine bias (left/right) based on relative X-coordinates
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
                        except (AttributeError, RuntimeError, TypeError):
                            # If sip check fails, continue defensively.
                            # This usually means the object is in an inconsistent state.
                            pass  # Silent failure for non-critical partner state check

                        other_pos = None
                        try:
                            other_pos = other_atom.pos()
                        except (AttributeError, RuntimeError, TypeError):
                            # Accessing .pos() may raise if the C++ object was destroyed
                            other_pos = None

                        if other_pos is None:
                            continue

                        total_dx += other_pos.x() - my_pos_x
                    except (AttributeError, RuntimeError, TypeError):
                        # Skip any problematic bond/partner rather than crashing the paint
                        continue

                # Flip text if bonds are primarily on the right
                if total_dx > 0:
                    flip_text = True

            # --- Finalize display text and alignment ---
            if flip_text:
                display_text = hydrogen_part + self.symbol
                alignment_flag = (
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
            else:
                display_text = self.symbol + hydrogen_part
                alignment_flag = (
                    Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
                )

            text_rect = fm.boundingRect(display_text)
            text_rect.adjust(-2, -2, 2, 2)
            symbol_rect = fm.boundingRect(self.symbol)

            # --- Determine text drawing position ---
            # No hydrogen label (center align as before)
            if not hydrogen_part:
                alignment_flag = Qt.AlignmentFlag.AlignCenter
                text_rect.moveCenter(QPointF(0, 0).toPoint())
            # Hydrogen label exists and is flipped (right align)
            elif flip_text:
                # Adjust right edge to center the main element
                offset_x = symbol_rect.width() // 2
                text_rect.moveTo(offset_x - text_rect.width(), -text_rect.height() // 2)
            # Non-flipped H-label (left align)
            else:
                # Adjust left edge to center the main element
                offset_x = -symbol_rect.width() // 2
                text_rect.moveTo(offset_x, -text_rect.height() // 2)

            # 2. Handle background (fill with white or clear if transparent)
            if self.scene():
                bg_brush = self.scene().backgroundBrush()
                bg_rect = text_rect.adjusted(-5, -8, 5, 8)

                if bg_brush.style() == Qt.BrushStyle.NoBrush:
                    # Use CompositionMode_Clear to erase overlapping bond lines
                    painter.save()
                    painter.setCompositionMode(
                        QPainter.CompositionMode.CompositionMode_Clear
                    )
                    painter.setBrush(
                        QColor(0, 0, 0, 255)
                    )  # Color doesn't matter (alpha is key)
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.drawEllipse(bg_rect)
                    painter.restore()
                else:
                    # Fill with background color if it exists
                    painter.setBrush(bg_brush)
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.drawEllipse(bg_rect)

            # 3. Draw the atom symbol itself
            # Color is already determined above
            painter.setPen(QPen(color))
            painter.drawText(text_rect, int(alignment_flag), display_text)

            # --- Draw charge and radical ---
            if self.charge != 0:
                # Chemical convention: single charge as "+"/"-", multiple as "2+"/"2-"
                if self.charge == 1:
                    charge_str = "+"
                elif self.charge == -1:
                    charge_str = "-"
                else:
                    sign = "+" if self.charge > 0 else "-"
                    charge_str = f"{abs(self.charge)}{sign}"
                charge_font = QFont("Arial", 12, QFont.Weight.Bold)
                painter.setFont(charge_font)
                charge_rect = painter.fontMetrics().boundingRect(charge_str)
                # Charge position also supports flipping
                if flip_text:
                    charge_pos = QPointF(
                        text_rect.left() - charge_rect.width() - 2,
                        text_rect.top() + charge_rect.height() - 2,
                    )
                else:
                    charge_pos = QPointF(
                        text_rect.right() + 2,
                        text_rect.top() + charge_rect.height() - 2,
                    )
                painter.setPen(Qt.GlobalColor.black)
                painter.drawText(charge_pos, charge_str)

            if self.radical > 0:
                painter.setBrush(QBrush(Qt.GlobalColor.black))
                painter.setPen(Qt.PenStyle.NoPen)
                radical_pos_y = text_rect.top() - 5
                if self.radical == 1:
                    painter.drawEllipse(
                        QPointF(text_rect.center().x(), radical_pos_y), 3, 3
                    )
                elif self.radical == 2:
                    painter.drawEllipse(
                        QPointF(text_rect.center().x() - 5, radical_pos_y), 3, 3
                    )
                    painter.drawEllipse(
                        QPointF(text_rect.center().x() + 5, radical_pos_y), 3, 3
                    )

        # --- Selection highlights etc. ---
        if self.has_problem:
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(255, 0, 0, 200), 4))
            painter.drawRect(self.boundingRect())
        elif self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.boundingRect())
        if (not self.isSelected()) and getattr(self, "hovered", False):
            pen = QPen(QColor(144, 238, 144, 200), 3)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(pen)
            painter.drawRect(self.boundingRect())

    def itemChange(self, change: QGraphicsItem.GraphicsItemChange, value: Any) -> Any:
        res = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            if self.flags() & QGraphicsItem.GraphicsItemFlag.ItemIsMovable:
                for bond in self.bonds:
                    if bond.scene():
                        bond.update_position()

        return res

    def hoverEnterEvent(self, event: Any) -> None:
        # Enable highlight on hover regardless of scene mode
        self.hovered = True
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Any) -> None:
        if self.hovered:
            self.hovered = False
            self.update()
        super().hoverLeaveEvent(event)
