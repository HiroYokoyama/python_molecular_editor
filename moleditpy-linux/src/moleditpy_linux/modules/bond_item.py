#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import QLineF, QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontMetricsF,
    QPainterPath,
    QPainterPathStroker,
    QPen,
    QPolygonF,
)
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene

try:
    from .constants import (
        DESIRED_BOND_PIXEL_WIDTH,
        EZ_LABEL_BOX_SIZE,
        EZ_LABEL_MARGIN,
        EZ_LABEL_TEXT_OUTLINE,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH,
    )
except ImportError:
    from modules.constants import (
        DESIRED_BOND_PIXEL_WIDTH,
        EZ_LABEL_BOX_SIZE,
        EZ_LABEL_MARGIN,
        EZ_LABEL_TEXT_OUTLINE,
        FONT_FAMILY,
        FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH,
    )


class BondItem(QGraphicsItem):
    def get_ez_label_rect(self):
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

    def set_stereo(self, new_stereo):
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
                try:
                    self.scene().views()[0].viewport().update()
                except (IndexError, RuntimeError):
                    # Handle case where views are being destroyed
                    pass

        except (AttributeError, RuntimeError, TypeError) as e:
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
        self.atom1, self.atom2, self.order, self.stereo = (
            atom1_item,
            atom2_item,
            order,
            stereo,
        )
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        self.update_position()
        self.setAcceptHoverEvents(True)
        self.setAcceptHoverEvents(True)
        self.hovered = False
        self.order = order
        self.stereo = stereo

    def get_line_in_local_coords(self):
        if self.atom1 is None or self.atom2 is None:
            return QLineF(0, 0, 0, 0)
        try:
            # Use pos() directly - assuming items are in scene coords (no parent)
            # This is robust and efficient.
            p1 = self.atom1.pos()
            p2 = self.atom2.pos()
            return QLineF(QPointF(0, 0), p2 - p1)
        except (AttributeError, RuntimeError):
            return QLineF(0, 0, 0, 0)

    def boundingRect(self):
        try:
            line = self.get_line_in_local_coords()
        except (AttributeError, RuntimeError):
            line = QLineF(0, 0, 0, 0)

        # Get dynamic bond offset (spacing)
        bond_offset = 3.5
        try:
            if self.scene() and hasattr(self.scene(), "views") and self.scene().views():
                win = self.scene().views()[0].window()
                if win and hasattr(win, "settings"):
                    # Use specific spacing based on bond order
                    if getattr(self, "order", 1) == 3:
                        bond_offset = win.settings.get("bond_spacing_triple_2d", 3.5)
                    else:
                        bond_offset = win.settings.get("bond_spacing_double_2d", 3.5)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            bond_offset = globals().get("BOND_OFFSET", 3.5)

        # Get dynamic wedge width
        wedge_width = 6.0
        try:
            if self.scene() and self.scene().views():
                win = self.scene().views()[0].window()
                if win and hasattr(win, "settings"):
                    wedge_width = win.settings.get("bond_wedge_width_2d", 6.0)
        except (AttributeError, RuntimeError, TypeError, ValueError):  # pragma: no cover
            import traceback
            traceback.print_exc()

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
            try:
                if self.scene() and self.scene().views():
                    win = self.scene().views()[0].window()
                    if win and hasattr(win, "settings"):
                        font_size = win.settings.get("atom_font_size_2d", 20)
                        font_family = win.settings.get(
                            "atom_font_family_2d", FONT_FAMILY
                        )
            except (AttributeError, RuntimeError, TypeError, ValueError):  # pragma: no cover
                import traceback
                traceback.print_exc()

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

    def shape(self):
        """Define the precise collision/selection area, separate from the drawing area (boundingRect)."""
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
            path.addRect(QRectF(-5, -5, 10, 10))

        return path

    def get_ez_label_local_rect(self):
        """Helper to get E/Z label rect in local coordinates."""
        if self.order != 2 or self.stereo not in [3, 4]:
            return None
        try:
            line = self.get_line_in_local_coords()
            center = line.center()

            # Logic similar to boundingRect but returning just the label box
            font_size = 20
            # ... (Simpler logic: just return a box around center)
            # Standard size estimate
            box_size = 30
            return QRectF(
                center.x() - box_size / 2, center.y() - box_size / 2, box_size, box_size
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return None

    def paint(self, painter, option, widget):
        if self.atom1 is None or self.atom2 is None:
            return
        line = self.get_line_in_local_coords()
        if line.length() == 0:
            return

        # Default values
        width_2d = 2.0
        wedge_width_half = 6.0
        num_dashes = 8
        bond_color = QColor("#222222")

        try:
            sc = self.scene()
            if sc is not None and hasattr(sc, "window") and sc.window is not None:
                # Get settings
                settings = sc.window.settings

                # Width
                width_2d = settings.get("bond_width_2d", 2.0)

                # Cap Style logic
                cap_style_str = settings.get("bond_cap_style_2d", "Round")
                cap_style = Qt.PenCapStyle.RoundCap  # Default

                if cap_style_str == "Flat":
                    cap_style = Qt.PenCapStyle.FlatCap
                elif cap_style_str == "Square":
                    cap_style = Qt.PenCapStyle.SquareCap

                # Color
                if self.isSelected():
                    bond_color = QColor("blue")  # Selection color
                else:
                    bond_hex = settings.get("bond_color_2d", "#222222")
                    bond_color = QColor(bond_hex)

                pen = QPen(bond_color, width_2d)
                pen.setCapStyle(cap_style)
                painter.setPen(pen)

                # Wedge/Dash Specific Settings
                wedge_width_half = settings.get("bond_wedge_width_2d", 6.0)
                num_dashes = int(settings.get("bond_dash_count_2d", 8))

            # Use bond color for fill
            painter.setBrush(QBrush(bond_color))
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
                try:
                    sc = self.scene()
                    if (
                        sc
                        and sc.views()
                        and hasattr(sc.views()[0].window(), "settings")
                    ):
                        if self.order == 3:
                            bond_offset = (
                                sc.views()[0]
                                .window()
                                .settings.get("bond_spacing_triple_2d", 3.5)
                            )
                        else:
                            bond_offset = (
                                sc.views()[0]
                                .window()
                                .settings.get("bond_spacing_double_2d", 3.5)
                            )
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    bond_offset = globals().get("BOND_OFFSET", 3.5)

                offset = QPointF(v.dx(), v.dy()) * bond_offset

                if self.order == 2:
                    # Determine if part of a ring structure to adjust drawing style
                    is_in_ring = False
                    ring_center = None

                    try:
                        # Get RDKit molecule from scene
                        sc = self.scene()
                        if sc and hasattr(sc, "window") and sc.window:
                            # Generate RDKit molecule from 2D data
                            mol = sc.window.data.to_rdkit_mol(use_2d_stereo=False)
                            if mol:
                                # Find RDKit bond corresponding to this editor bond
                                atom1_id = self.atom1.atom_id
                                atom2_id = self.atom2.atom_id

                                # Get RDKit indices
                                rdkit_idx1 = None
                                rdkit_idx2 = None
                                for atom in mol.GetAtoms():
                                    if atom.HasProp("_original_atom_id"):
                                        orig_id = atom.GetIntProp("_original_atom_id")
                                        if orig_id == atom1_id:
                                            rdkit_idx1 = atom.GetIdx()
                                        elif orig_id == atom2_id:
                                            rdkit_idx2 = atom.GetIdx()

                                if rdkit_idx1 is not None and rdkit_idx2 is not None:
                                    bond = mol.GetBondBetweenAtoms(
                                        rdkit_idx1, rdkit_idx2
                                    )
                                    if bond and bond.IsInRing():
                                        is_in_ring = True
                                        # Calculate ring center (smallest ring containing this bond)
                                        ring_info = mol.GetRingInfo()
                                        for ring in ring_info.AtomRings():
                                            if (
                                                rdkit_idx1 in ring
                                                and rdkit_idx2 in ring
                                            ):
                                                # Calculate average position of atoms in ring
                                                ring_positions = []
                                                for atom_idx in ring:
                                                    # Find corresponding atom in editor
                                                    rdkit_atom = mol.GetAtomWithIdx(
                                                        atom_idx
                                                    )
                                                    if rdkit_atom.HasProp(
                                                        "_original_atom_id"
                                                    ):
                                                        editor_atom_id = (
                                                            rdkit_atom.GetIntProp(
                                                                "_original_atom_id"
                                                            )
                                                        )
                                                        if (
                                                            editor_atom_id
                                                            in sc.window.data.atoms
                                                        ):
                                                            atom_item = (
                                                                sc.window.data.atoms[
                                                                    editor_atom_id
                                                                ]["item"]
                                                            )
                                                            if atom_item:
                                                                ring_positions.append(
                                                                    atom_item.pos()
                                                                )

                                                if ring_positions:
                                                    # Calculate ring center
                                                    center_x = sum(
                                                        p.x() for p in ring_positions
                                                    ) / len(ring_positions)
                                                    center_y = sum(
                                                        p.y() for p in ring_positions
                                                    ) / len(ring_positions)
                                                    ring_center = QPointF(
                                                        center_x, center_y
                                                    )
                                                    break
                    except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                        # Fallback to default drawing on error
                        is_in_ring = False

                    v = line.unitVector().normalVector()
                    # Re-calculate offset in case loop variable scope issue, though strictly not needed if offset defined above works
                    offset = QPointF(v.dx(), v.dy()) * bond_offset

                    if is_in_ring and ring_center:
                        # Ring structure: 1 central line (single bond pos) + 1 short inner line
                        # Calculate direction from bond center to ring center
                        bond_center = line.center()

                        # Ring center direction in local coords
                        local_ring_center = self.mapFromScene(ring_center)
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
                        except (AttributeError, RuntimeError, TypeError, ValueError):  # pragma: no cover
                            import traceback
                            traceback.print_exc()

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
            except (AttributeError, RuntimeError, TypeError, ValueError):  # pragma: no cover
                import traceback
                traceback.print_exc()

    def update_position(self, notify=True):
        try:
            if notify:
                self.prepareGeometryChange()
            if self.atom1:
                self.setPos(self.atom1.pos())
            self.update()
        except (AttributeError, RuntimeError) as e:
            print(f"Error updating bond position: {e}")
            # Continue without crashing

    def hoverEnterEvent(self, event):
        scene = self.scene()
        mode = getattr(scene, "mode", "")
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
