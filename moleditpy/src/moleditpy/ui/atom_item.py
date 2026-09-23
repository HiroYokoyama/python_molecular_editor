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
import logging
from typing import TYPE_CHECKING, Any, List, Optional

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
from PyQt6.QtWidgets import (
    QGraphicsItem,
    QGraphicsSceneHoverEvent,
    QStyleOptionGraphicsItem,
    QWidget,
)

from ..utils.constants import (
    ATOM_RADIUS,
    CPK_COLORS,
    FONT_FAMILY,
    FONT_WEIGHT_BOLD,
)
from ..utils.hit_radius import scene_hit_radius
from ..utils.sip_isdeleted_safe import sip_isdeleted_safe

if TYPE_CHECKING:
    from .bond_item import BondItem

SUBSCRIPT_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

# Every subscript digit, so one warm-up covers any implicit-H count.
_WARM_UP_TEXT = "CH" + "0123456789".translate(SUBSCRIPT_MAP)


def warm_font_cache(font: Optional[QFont] = None) -> None:
    """Pay Qt's one-off font costs up front instead of on the first heteroatom.

    Two process-wide costs hide in QFontMetricsF: the first text shaping in a
    family (~10-120 ms, worse on a cold Windows font cache) and the first
    lookup of a subscript digit, which Arial lacks — Qt answers that by walking
    the whole system font database (~200-330 ms, measured).

    Drawing a carbon skeleton pays neither: a bonded neutral carbon is
    invisible, so it renders no glyph at all. The first O or N a user places
    gets "OH₂" and pays both at once, mid-click, as a visible freeze.

    Call once on the GUI thread after the window is up. QFontMetricsF is not
    usable off the GUI thread, so this cannot move to a worker.
    """
    try:
        QFontMetricsF(font or QFont(FONT_FAMILY, 20, FONT_WEIGHT_BOLD)).boundingRect(
            _WARM_UP_TEXT
        )
    except (RuntimeError, TypeError, ValueError):
        logging.debug("Font cache warm-up skipped", exc_info=True)


class AtomItem(QGraphicsItem):
    """2D scene item representing a single atom in the molecule editor."""

    def __init__(
        self, atom_id: int, symbol: str, pos: QPointF, charge: int = 0, radical: int = 0
    ) -> None:
        """Initialize 2D atom graphics item."""
        super().__init__()
        self.atom_id: int = atom_id
        self.symbol: str = symbol
        self.charge: int = charge
        self.radical: int = radical
        self.bonds: List[BondItem] = []
        self.chiral_label: Optional[str] = None

        self.setPos(pos)
        self.implicit_h_count: int = 0
        self.setFlags(
            QGraphicsItem.GraphicsItemFlag.ItemIsMovable
            | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable
            # Without this Qt never calls itemChange for moves, so bonds lag.
            | QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges
        )
        self.setZValue(1)
        self.font: QFont = QFont(FONT_FAMILY, 22, FONT_WEIGHT_BOLD)
        self.update_style()
        self.setAcceptHoverEvents(True)
        self.hovered: bool = False
        self.has_problem: bool = False
        self.is_visible: bool = True

    def update_style(self) -> None:
        """Refresh font, color, and visibility based on current scene settings."""
        if sip_isdeleted_safe(self):
            return
        self.font = self._label_font()
        self.prepareGeometryChange()

        self.is_visible = not (
            self.symbol == "C"
            and len(self.bonds) > 0
            and self.charge == 0
            and self.radical == 0
        )
        self.update()

    def _label_font(self) -> QFont:
        """The atom label font from the scene's 2D settings (defaults off-scene)."""
        font_size = 22
        font_family = FONT_FAMILY
        font_bold = True
        font_italic = False
        font_underline = False
        scene: Any = self.scene()
        if scene is not None:
            font_size = scene.get_setting("atom_font_size_2d", 22)
            font_family = scene.get_setting("atom_font_family_2d", FONT_FAMILY)
            font_bold = scene.get_setting("atom_font_bold_2d", True)
            font_italic = scene.get_setting("atom_font_italic_2d", False)
            font_underline = scene.get_setting("atom_font_underline_2d", False)

        weight = QFont.Weight.Bold if font_bold else QFont.Weight.Normal
        font = QFont(font_family, font_size, weight)
        font.setItalic(font_italic)
        font.setUnderline(font_underline)
        return font

    def _hydrogen_part(self) -> str:
        """Implicit-H suffix such as "H" or "H₂"; empty for a skeletal carbon."""
        if self.implicit_h_count <= 0:
            return ""
        is_skeletal_carbon = (
            self.symbol == "C"
            and self.charge == 0
            and self.radical == 0
            and len(self.bonds) > 0
        )
        if is_skeletal_carbon:
            return ""
        hydrogen_part = "H"
        if self.implicit_h_count > 1:
            hydrogen_part += str(self.implicit_h_count).translate(SUBSCRIPT_MAP)
        return hydrogen_part

    def _label_flipped(self, hydrogen_part: str) -> bool:
        """Whether to write the label H-first ("H₂O"), away from the bonds.

        True when the bonds lie mostly to the right. Partners with a missing
        or deleted wrapper are skipped.
        """
        if not hydrogen_part or not self.bonds:
            return False
        my_pos_x = self.pos().x()
        total_dx = 0.0
        for bond in self.bonds:
            try:
                partner = bond.atom2 if bond.atom1 is self else bond.atom1
                if partner is None or sip_isdeleted_safe(partner):
                    continue
                partner_pos = partner.pos()
                if partner_pos is None:
                    continue
                total_dx += partner_pos.x() - my_pos_x
            except (AttributeError, RuntimeError, TypeError, ValueError):
                # The partner's C++ object can be mid-destruction.
                continue
        return total_dx > 0

    def _charge_text(self) -> str:
        """Charge label by chemical convention: "+"/"-" for ±1, "2+"/"2-" beyond."""
        if self.charge == 1:
            return "+"
        if self.charge == -1:
            return "-"
        sign = "+" if self.charge > 0 else "-"
        return f"{abs(self.charge)}{sign}"

    def _label_text_rect(self, hydrogen_part: str, flip_text: bool) -> QRectF:
        """Rectangle of the label text, element symbol centred on the atom."""
        fm = QFontMetricsF(self._label_font())
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
        return text_rect

    def visual_rect(self) -> QRectF:
        """Return the rectangle the atom draws into; highlights use this."""
        hydrogen_part = self._hydrogen_part()
        flip_text = self._label_flipped(hydrogen_part)
        text_rect = self._label_text_rect(hydrogen_part, flip_text)

        # 1. Calculate the background rectangle (bg_rect) used in paint()
        bg_rect = text_rect.adjusted(-5, -8, 5, 8)

        # 2. Construct the full visual rectangle relative to bg_rect
        full_visual_rect = QRectF(bg_rect)

        # Include charge symbol area in calculation
        if self.charge != 0:
            charge_font = QFont("Arial", 12, QFont.Weight.Bold)
            charge_rect = QFontMetricsF(charge_font).boundingRect(self._charge_text())

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

    def boundingRect(self) -> QRectF:
        """Return the drawn rect plus the hit shape, which Qt requires it to cover."""
        hit_r = self.hit_radius()
        return self.visual_rect().united(
            QRectF(-hit_r, -hit_r, hit_r * 2.0, hit_r * 2.0)
        )

    def get_bg_ellipse_path(self) -> QPainterPath:
        """Return the elliptical background path used for bond endpoint clipping."""
        path = QPainterPath()
        if not self.is_visible:
            return path

        hydrogen_part = self._hydrogen_part()
        text_rect = self._label_text_rect(
            hydrogen_part, self._label_flipped(hydrogen_part)
        )
        path.addEllipse(text_rect.adjusted(-5, -8, 5, 8))
        return path

    def hit_radius(self) -> float:
        """Return the hit/snap radius in scene units for the current zoom."""
        return scene_hit_radius(self.scene(), fallback=max(4.0, ATOM_RADIUS - 6.0) * 2)

    def shape(self) -> QPainterPath:
        """Define the collision area: a fixed pixel radius, never the drawn label."""
        scene_radius = self.hit_radius()

        path = QPainterPath()
        path.addEllipse(QPointF(0, 0), scene_radius, scene_radius)
        return path

    def paint(
        self,
        painter: Optional[QPainter],
        option: QStyleOptionGraphicsItem,  # type: ignore[override]
        widget: Optional[QWidget] = None,
    ) -> None:
        """Paint the atom symbol and its associated labels (charge, radical)."""
        if painter is None:
            return
        # Color logic: check if we should use bond color (uniform) or CPK (element-specific)
        color = CPK_COLORS.get(self.symbol, CPK_COLORS["DEFAULT"])
        # Use bond color if specified in settings
        scene: Any = self.scene()
        if scene is not None and (
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

            hydrogen_part = self._hydrogen_part()
            flip_text = self._label_flipped(hydrogen_part)

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

            # 3. Draw the atom symbol itself
            # Color is already determined above
            painter.setPen(QPen(color))
            painter.drawText(text_rect, int(alignment_flag), display_text)

            # --- Draw charge and radical ---
            if self.charge != 0:
                charge_str = self._charge_text()
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
            painter.drawRect(self.visual_rect())
        elif self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.visual_rect())
        if (not self.isSelected()) and getattr(self, "hovered", False):
            pen = QPen(QColor(144, 238, 144, 200), 3)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(pen)
            painter.drawRect(self.visual_rect())

    def itemChange(self, change: QGraphicsItem.GraphicsItemChange, value: Any) -> Any:
        """Propagate position and scene changes to connected bond items."""
        res = super().itemChange(change, value)
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            if self.flags() & QGraphicsItem.GraphicsItemFlag.ItemIsMovable:
                for bond in self.bonds:
                    if bond.scene():
                        bond.update_position()
        elif change == QGraphicsItem.GraphicsItemChange.ItemSceneHasChanged:
            if self.scene() is not None:
                self.update_style()
        return res

    def hoverEnterEvent(self, event: QGraphicsSceneHoverEvent) -> None:  # type: ignore[override]
        """Highlight the atom on mouse hover."""
        # Enable highlight on hover regardless of scene mode
        self.hovered = True
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: QGraphicsSceneHoverEvent) -> None:  # type: ignore[override]
        """Remove hover highlight when the mouse leaves."""
        if self.hovered:
            self.hovered = False
            self.update()
        super().hoverLeaveEvent(event)
