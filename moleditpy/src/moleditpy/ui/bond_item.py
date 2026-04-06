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
from typing import Any, Optional, Tuple, Union

from PyQt6.QtCore import QLineF, QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontMetricsF,
    QPainter,
    QPainterPath,
    QPainterPathStroker,
    QPen,
    QPolygonF,
)
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene, QWidget

try:
    from ..utils.constants import (
        DESIRED_BOND_PIXEL_WIDTH,
        EZ_LABEL_BOX_SIZE,
        EZ_LABEL_MARGIN,
        EZ_LABEL_TEXT_OUTLINE,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH,
    )
except ImportError:
    from moleditpy.utils.constants import (
        DESIRED_BOND_PIXEL_WIDTH,
        EZ_LABEL_BOX_SIZE,
        EZ_LABEL_MARGIN,
        EZ_LABEL_TEXT_OUTLINE,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH,
    )


class BondItem(QGraphicsItem):
    """Visual representation of a molecular bond in the 2D scene."""

    def get_ez_label_rect(self) -> Optional[QRectF]:
        """Returns the drawing range for E/Z labels (scene coords). Returns None if no label."""
        if self.order != 2 or self.stereo not in [3, 4]:
            return None
        line = self.get_line_in_local_coords()
        center = line.center()
        label_width = EZ_LABEL_BOX_SIZE
        label_height = EZ_LABEL_BOX_SIZE
        label_rect = QRectF(
            center.x() - label_width / 2,
            center.y() - label_height / 2,
            label_width,
            label_height,
        )
        # Convert to scene coordinates
        return self.mapToScene(label_rect).boundingRect()

    def set_stereo(self, new_stereo: int) -> None:
        """Update the stereochemical state of the bond."""
        try:
            # Invalidate scene area when removing label
            if new_stereo == 0 and self.stereo in [3, 4] and self.scene():
                rect = self.mapToScene(self.boundingRect()).boundingRect()
                self.scene().invalidate(
                    rect,
                    QGraphicsScene.SceneLayer.BackgroundLayer
                    | QGraphicsScene.SceneLayer.ForegroundLayer,
                )

            self.prepareGeometryChange()
            self.stereo = new_stereo
            self.update()

            if self.scene() and self.scene().views():
                view = self.scene().views()[0]
                if view and view.viewport():
                    view.viewport().update()

        except (AttributeError, RuntimeError, TypeError) as e:
            print(f"Error in BondItem.set_stereo: {e}")
            # Continue without crashing
            self.stereo = new_stereo

    def set_order(self, new_order: int) -> None:
        """Update the bond order."""
        self.prepareGeometryChange()
        self.order = new_order
        self.update()
        if self.scene() and self.scene().views():
            self.scene().views()[0].viewport().update()

    def update_style(self) -> None:
        """Force internal state refresh and redraw (primarily from settings)."""
        self.prepareGeometryChange()
        self.update()

    def __init__(
        self, atom1_item: Any, atom2_item: Any, order: int = 1, stereo: int = 0
    ) -> None:
        super().__init__()
        # Validate input parameters
        if atom1_item is None or atom2_item is None:
            raise ValueError("BondItem requires non-None atom items")
        self.atom1: Any = atom1_item
        self.atom2: Any = atom2_item
        self.order: int = order
        self.stereo: int = stereo

        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen: QPen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        self.update_position()
        self.setAcceptHoverEvents(True)
        self.hovered: bool = False
        self.is_in_ring: bool = False
        self.ring_center: Optional[Union[QPointF, Tuple[float, float]]] = None

    def get_line_in_local_coords(self) -> QLineF:
        if self.atom1 is None or self.atom2 is None:
            return QLineF(0, 0, 0, 0)
        try:
            # Use pos() directly - assuming items are in scene coords (no parent)
            # This is robust and efficient.
            p1 = self.atom1.pos()
            p2 = self.atom2.pos()
            return QLineF(QPointF(0, 0), p2 - p1)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            # Fallback for inconsistent/deleted atom references
            # return zero line to prevent downstream crashes.
            return QLineF(0, 0, 0, 0)

    def boundingRect(self) -> QRectF:
        try:
            line = self.get_line_in_local_coords()
        except (AttributeError, RuntimeError, ValueError, TypeError):
            line = QLineF(0, 0, 0, 0)

        # Get dynamic bond offset (spacing)
        scene = self.scene()
        bond_offset = 3.5
        wedge_width = 6.0
        if scene is not None:
            if hasattr(scene, "get_setting"):
                # Use specific spacing based on bond order
                key = (
                    "bond_spacing_triple_2d"
                    if getattr(self, "order", 1) == 3
                    else "bond_spacing_double_2d"
                )
                bond_offset = scene.get_setting(key, 3.5)
                wedge_width = scene.get_setting("bond_wedge_width_2d", 6.0)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                )

        extra = (getattr(self, "order", 1) - 1) * bond_offset + 50 + wedge_width
        rect = (
            QRectF(line.p1(), line.p2())
            .normalized()
            .adjusted(-extra, -extra, extra, extra)
        )

        # Expand bounding rect for E/Z labels (calculating precisely with QFontMetricsF)
        if self.order == 2 and self.stereo in [3, 4]:
            font_size = 20
            font_family = FONT_FAMILY
            if scene is not None:
                if hasattr(scene, "get_setting"):
                    font_size = scene.get_setting("atom_font_size_2d", 20)
                    font_family = scene.get_setting("atom_font_family_2d", FONT_FAMILY)
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error(
                        f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                    )

            font = QFont(font_family, font_size, FONT_WEIGHT_BOLD)
            font.setItalic(True)
            text = "Z" if self.stereo == 3 else "E"
            fm = QFontMetricsF(font)
            text_rect = fm.boundingRect(text)
            outline = EZ_LABEL_TEXT_OUTLINE  # Outline thickness
            margin = EZ_LABEL_MARGIN  # Additional margin
            center = line.center()
            label_rect = QRectF(
                center.x() - text_rect.width() / 2 - outline - margin,
                center.y() - text_rect.height() / 2 - outline - margin,
                text_rect.width() + 2 * outline + 2 * margin,
                text_rect.height() + 2 * outline + 2 * margin,
            )
            rect = rect.united(label_rect)
        return rect

    def shape(self) -> QPainterPath:
        """Define the precise collision/selection area."""
        path = QPainterPath()
        try:
            line = self.get_line_in_local_coords()
            # Create a simple path along the bond line
            path.moveTo(line.p1())
            path.lineTo(line.p2())

            # Stroke it to give it some width (e.g., 10px or dynamic based on settings) generally easier to click
            # even if the visual width is smaller.
            stroker = QPainterPathStroker()
            stroker.setWidth(DESIRED_BOND_PIXEL_WIDTH)  # Use constant (20.0)
            path = stroker.createStroke(path)

            # If there's an E/Z label, add its rect to the selection shape
            label_rect = self.get_ez_label_local_rect()
            if label_rect:
                path.addRect(label_rect)

        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Fallback to a small rect around the origin if calculation fails
            # This is non-critical for collision detection if points are missing.
            path.addRect(QRectF(-5, -5, 10, 10))

        return path

    def get_ez_label_local_rect(self) -> Optional[QRectF]:
        """Helper to get E/Z label rect in local coordinates."""
        if self.order != 2 or self.stereo not in [3, 4]:
            return None
        try:
            line = self.get_line_in_local_coords()
            center = line.center()

            # Logic similar to boundingRect but returning just the label box
            # ... (Simpler logic: just return a box around center)
            # Standard size estimate
            box_size = 30
            return QRectF(
                center.x() - box_size / 2, center.y() - box_size / 2, box_size, box_size
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return None

    def paint(
        self,
        painter: QPainter,
        option: Any,
        widget: Optional[QWidget] = None,
    ) -> None:
        if self.atom1 is None or self.atom2 is None:
            return
        line = self.get_line_in_local_coords()
        if line.length() == 0:
            return

        # Default values
        wedge_width_half = 6.0
        num_dashes = 8
        bond_color = QColor("#222222")

        try:
            scene = self.scene()
            if scene is not None:
                if hasattr(scene, "get_setting"):
                    # Color
                    if self.isSelected():
                        bond_color = QColor("blue")
                    else:
                        color_hex = scene.get_setting("bond_color_2d", "#222222")
                        bond_color = QColor(color_hex)

                    # Width and Cap Style
                    bond_width = scene.get_setting("bond_width_2d", 2.0)
                    cap_style_str = scene.get_setting("bond_cap_style_2d", "Round")

                    # Wedge/Dash settings
                    wedge_width_half = scene.get_setting("bond_wedge_width_2d", 6.0)
                    num_dashes = int(scene.get_setting("bond_dash_count_2d", 8))

                    # Cap Style logic
                    cap_style = Qt.PenCapStyle.RoundCap
                    if cap_style_str == "Flat":
                        cap_style = Qt.PenCapStyle.FlatCap
                    elif cap_style_str == "Square":
                        cap_style = Qt.PenCapStyle.SquareCap

                    pen = QPen(bond_color, bond_width)
                    pen.setCapStyle(cap_style)
                    painter.setPen(pen)
                    painter.setBrush(QBrush(bond_color))
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error(
                        f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                    )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Fallback
            painter.setPen(self.pen)
            painter.setBrush(QBrush(Qt.GlobalColor.black))
            # Default fallback width for wedge
            wedge_width_half = 6.0
            num_dashes = 8

        # --- Draw Stereochemistry (Wedge/Dash) ---
        if self.order == 1 and self.stereo in [1, 2]:
            vec = line.unitVector()
            normal = vec.normalVector()
            p1 = line.p1() + vec.p2() * 5
            p2 = line.p2() - vec.p2() * 5

            if self.stereo == 1:  # Wedge
                offset = QPointF(normal.dx(), normal.dy()) * wedge_width_half
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)

            elif self.stereo == 2:  # Dash
                painter.save()
                if not self.isSelected():
                    pen = painter.pen()
                    pen.setWidthF(2.5)
                    painter.setPen(pen)

                # Use configured number of dashes (default 8)
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = (wedge_width_half * 2.0) * t
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.drawLine(start_pt - offset, start_pt + offset)
                painter.restore()

        # --- Draw Regular Bonds (Single/Double/Triple) ---
        else:
            if self.order == 1:
                painter.drawLine(line)
            else:
                v = line.unitVector().normalVector()
                # Use dynamic offset
                bond_offset = 3.5
                scene = self.scene()
                if scene is not None:
                    if hasattr(scene, "get_setting"):
                        key = (
                            "bond_spacing_triple_2d"
                            if self.order == 3
                            else "bond_spacing_double_2d"
                        )
                        bond_offset = scene.get_setting(key, 3.5)
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error(
                            f"REPORT ERROR: Missing attribute 'get_setting' on scene of type {type(scene)}"
                        )

                offset = QPointF(v.dx(), v.dy()) * bond_offset

                if self.order == 2:
                    # Use cached ring info
                    is_in_ring = self.is_in_ring
                    ring_center = self.ring_center

                    v = line.unitVector().normalVector()
                    # Re-calculate offset in case loop variable scope issue, though strictly not needed if offset defined above works
                    offset = QPointF(v.dx(), v.dy()) * bond_offset

                    if is_in_ring and ring_center:
                        # Ring structure: 1 central line (single bond pos) + 1 short inner line
                        # Ring center direction in local coords
                        # ring_center may be QPointF (from old/legacy code) or tuple (from new decoupled core)
                        if isinstance(ring_center, (list, tuple)):
                            # Use iterator to satisfy strict linters about subscripting/unpacking
                            it = iter(ring_center)
                            rc_pt = QPointF(float(next(it)), float(next(it)))
                        else:
                            rc_pt = ring_center
                        local_ring_center = self.mapFromScene(rc_pt)
                        local_bond_center = line.center()
                        inward_vec = local_ring_center - local_bond_center

                        # Use dot product to determine inner side
                        if QPointF.dotProduct(offset, inward_vec) > 0:
                            # Offset is inward (double offset)
                            inner_offset = offset * 2
                        else:
                            # Negative offset is inward (double offset)
                            inner_offset = -offset * 2

                        # Draw central line (same position as single bond)
                        painter.drawLine(line)

                        # Draw short inner line (80% length)
                        inner_line = line.translated(inner_offset)
                        shorten_factor = 0.8
                        p1 = inner_line.p1()
                        p2 = inner_line.p2()
                        center = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
                        shortened_p1 = center + (p1 - center) * shorten_factor
                        shortened_p2 = center + (p2 - center) * shorten_factor
                        painter.drawLine(QLineF(shortened_p1, shortened_p2))
                    else:
                        # Non-ring structure: parallel lines
                        line1 = line.translated(offset)
                        line2 = line.translated(-offset)
                        painter.drawLine(line1)
                        painter.drawLine(line2)

                    # E/Z Label Drawing
                    if self.stereo in [3, 4]:
                        painter.save()  # Save current painter state

                        # --- Label Settings ---
                        font_size = 20
                        font_family = FONT_FAMILY
                        try:
                            if self.scene() and self.scene().views():
                                win = self.scene().views()[0].window()
                                if win and hasattr(win, "settings"):
                                    font_size = win.settings.get(
                                        "atom_font_size_2d", 20
                                    )
                                    font_family = win.settings.get(
                                        "atom_font_family_2d", FONT_FAMILY
                                    )
                        except (AttributeError, RuntimeError, TypeError, ValueError):
                            # Silent failure for non-critical 2D atom font setting
                            # If we can't get custom settings, we just use defaults.
                            pass

                        font = QFont(font_family, font_size, FONT_WEIGHT_BOLD)
                        font.setItalic(True)
                        text_color = QColor("gray")
                        # Match outline color to background (safe fallback if scene() is None)
                        outline_color = None
                        try:
                            sc = self.scene()
                            if sc is not None:
                                outline_color = sc.backgroundBrush().color()
                        except (AttributeError, RuntimeError, TypeError):
                            outline_color = None
                        if outline_color is None:
                            # Default to white background outline for visibility
                            outline_color = QColor(255, 255, 255)

                        # --- Create Drawing Paths ---
                        text = "Z" if self.stereo == 3 else "E"
                        path = QPainterPath()

                        # Center text precisely
                        fm = QFontMetricsF(font)
                        text_rect = fm.boundingRect(text)
                        text_rect.moveCenter(line.center())
                        path.addText(text_rect.topLeft(), font, text)

                        # --- Draw Outline ---
                        stroker = QPainterPathStroker()
                        stroker.setWidth(EZ_LABEL_TEXT_OUTLINE)  # Outline width
                        outline_path = stroker.createStroke(path)

                        painter.setBrush(outline_color)
                        painter.setPen(Qt.PenStyle.NoPen)
                        painter.drawPath(outline_path)

                        # --- Draw Text Body ---
                        painter.setBrush(text_color)
                        painter.setPen(text_color)
                        painter.drawPath(path)

                        painter.restore()  # Restore painter state

                elif self.order == 3:
                    painter.drawLine(line)
                    painter.drawLine(line.translated(offset))
                    painter.drawLine(line.translated(-offset))

        # --- 2. Draw hover effects on top ---
        if (not self.isSelected()) and getattr(self, "hovered", False):
            try:
                # Draw highlight as thick semi-transparent line
                hover_pen = QPen(
                    QColor(144, 238, 144, 180), HOVER_PEN_WIDTH
                )  # LightGreen, semi-transparent
                hover_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                painter.setPen(hover_pen)
                painter.drawLine(line)
            except (AttributeError, RuntimeError, TypeError, ValueError):
                # Silent failure for non-critical hover highlight drawing.
                # If highlight fails, it's just a visual artifact.
                pass

    def update_position(self, notify: bool = True) -> None:
        try:
            if notify:
                self.prepareGeometryChange()
            if self.atom1:
                self.setPos(self.atom1.pos())
            self.update()
        except (AttributeError, RuntimeError) as e:
            print(f"Error updating bond position: {e}")
            # Continue without crashing

    def hoverEnterEvent(self, event: Any) -> None:
        self.hovered = True
        self.update()
        if self.scene():
            self.scene().set_hovered_item(self)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event: Any) -> None:
        if self.hovered:
            self.hovered = False
            self.update()
        if self.scene():
            self.scene().set_hovered_item(None)
        super().hoverLeaveEvent(event)
