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

import math
from typing import Any, Dict, List, Optional, Set

from PyQt6.QtCore import QLineF, QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QPainter,
    QPainterPath,
    QPen,
    QPolygonF,
)
from PyQt6.QtWidgets import QGraphicsItem, QStyleOptionGraphicsItem, QWidget

from .atom_item import AtomItem
from .bond_item import BondItem
from .preview_molecule import PreviewScene, build_preview_items

GHOST_OPACITY = 0.55
GHOST_COLOR = QColor(80, 80, 80, 180)
LABEL_BACKING = QColor(255, 255, 255, 240)
FIRST_ATOM_FILL = QColor(175, 215, 255, 85)
FIRST_ATOM_BORDER = QColor(140, 195, 245, 120)
FIRST_ATOM_RADIUS = 15.0
SNAP_CELL = 4.0
SNAP_TOLERANCE = 2.0


class TemplatePreviewItem(QGraphicsItem):
    """Ghost overlay item that previews ring or user-template placement."""

    def __init__(self) -> None:
        super().__init__()
        self.setZValue(2)
        self.pen = QPen(GHOST_COLOR, 2)
        self.polygon = QPolygonF()
        self.is_aromatic = False
        self.user_template_points: List[QPointF] = []
        self.user_template_bonds: List[Any] = []
        self.user_template_atoms: List[Any] = []
        self.is_user_template = False
        self.chain_points: List[QPointF] = []
        self.chain_label = ""
        self.chain_label_pos = QPointF()
        self.is_chain = False
        self.ghost_scene = PreviewScene()
        self.ghost_atoms: List[AtomItem] = []
        self.ghost_bonds: List[BondItem] = []
        self.existing_indices: Set[int] = set()
        self.existing_h_counts: Dict[int, int] = {}
        self.replaced_label_path = QPainterPath()
        self.mark_first_atom = False
        self.ghost_signature: Any = None
        self.ring_atoms: Dict[int, Any] = {}
        self.cached_rect = QRectF()

    def set_geometry(self, points: list[QPointF], is_aromatic: bool = False) -> None:
        """Set polygon points for a standard ring template preview."""
        self.prepareGeometryChange()
        self.polygon = QPolygonF(points)
        self.is_aromatic = is_aromatic
        self.is_user_template = False
        self.is_chain = False
        self.mark_first_atom = False

        # add_molecule_fragment lays aromatic rings out as alternating Kekulé bonds
        n = len(points)
        bonds_info = [
            (i, (i + 1) % n, 2 if (is_aromatic and i % 2 == 0) else 1, 0)
            for i in range(n)
        ]
        self._rebuild_ghost(points, bonds_info, [{"symbol": "C"} for _ in points])
        self.update()

    def set_chain_geometry(
        self, points: list[QPointF], label: str, label_pos: QPointF
    ) -> None:
        """Set the zigzag vertices and the repeat-count label for a chain preview."""
        self.prepareGeometryChange()
        self.chain_points = list(points)
        self.chain_label = label
        self.chain_label_pos = QPointF(label_pos)
        self.is_chain = True
        self.is_user_template = False
        self.is_aromatic = False
        self.mark_first_atom = False
        self.polygon = QPolygonF()
        self._clear_ghost()
        self.update()

    def set_user_template_geometry(
        self,
        points: list[QPointF],
        bonds_info: list[Any],
        atoms_data: list[dict[str, Any]],
    ) -> None:
        """Set point and bond data for a user-defined template preview."""
        self.prepareGeometryChange()
        self.user_template_points = points
        self.user_template_bonds = bonds_info
        self.user_template_atoms = atoms_data
        self.is_user_template = True
        self.is_chain = False
        self.is_aromatic = False
        self.polygon = QPolygonF()
        # The first template atom is what a clicked existing atom turns into
        self.mark_first_atom = True
        self._rebuild_ghost(points, bonds_info, atoms_data)
        self.update()

    # ------------------------------------------------------------------
    # Ghost molecule: real AtomItem/BondItem instances, off-scene
    # ------------------------------------------------------------------

    def _clear_ghost_items(self) -> None:
        """Drop the ghost atom and bond items."""
        for item in list(self.ghost_scene.items()):
            self.ghost_scene.removeItem(item)
        self.ghost_atoms = []
        self.ghost_bonds = []
        self.ring_atoms = {}
        self.ghost_signature = None
        self.cached_rect = QRectF()

    def _clear_ghost(self) -> None:
        """Drop the ghost molecule and everything derived from the editor scene."""
        self._clear_ghost_items()
        self.existing_indices = set()
        self.existing_h_counts = {}
        self.replaced_label_path = QPainterPath()

    def _rebuild_ghost(
        self,
        points: list[QPointF],
        bonds_info: list[Any],
        atoms_data: list[dict[str, Any]],
    ) -> None:
        """Rebuild the preview as real atom/bond items so it renders like the result."""
        self.ghost_scene.host = self.scene()
        if not points:
            self._clear_ghost()
            return

        existing_indices: Set[int] = set()
        existing_h_counts: Dict[int, int] = {}
        replaced_label_path = QPainterPath()
        preview_atoms: List[Dict[str, Any]] = []
        occupied = self._occupied_cells()
        for i, point in enumerate(points):
            atom_data: Dict[str, Any] = (
                dict(atoms_data[i]) if i < len(atoms_data) else {"symbol": "C"}
            )
            existing = self._existing_atom_at(point, occupied)
            replaces_existing = i == 0 and self.mark_first_atom
            if existing is not None and replaces_existing:
                # Its label has to be covered, or the ghost reads as "N over OH"
                replaced_label_path = existing.get_bg_ellipse_path().translated(
                    existing.pos()
                )
            elif existing is not None:
                # A fused vertex keeps the atom already drawn there
                existing_indices.add(i)
                # Its real H count comes from bonds the ghost molecule cannot see
                existing_h_counts[i] = getattr(existing, "implicit_h_count", 0)
                atom_data = {
                    "symbol": existing.symbol,
                    "charge": existing.charge,
                    "radical": existing.radical,
                }
            atom_data["pos"] = QPointF(point)
            preview_atoms.append(atom_data)

        signature = self._signature(preview_atoms, bonds_info, existing_h_counts)
        self.existing_indices = existing_indices
        self.existing_h_counts = existing_h_counts
        self.replaced_label_path = replaced_label_path

        if signature == self.ghost_signature:
            # Same molecule, new cursor position: moving it beats rebuilding it
            self._move_ghost(points)
            self._refresh_bounding_rect()
            self.prepareGeometryChange()
            return

        self._clear_ghost_items()
        self.ghost_atoms, self.ghost_bonds, self.ring_atoms = build_preview_items(
            preview_atoms, bonds_info
        )
        self.ghost_signature = signature
        for item in self.ghost_bonds:
            self.ghost_scene.addItem(item)
        for atom in self.ghost_atoms:
            self.ghost_scene.addItem(atom)

        for index, h_count in existing_h_counts.items():
            if index < len(self.ghost_atoms):
                self.ghost_atoms[index].implicit_h_count = h_count
                self.ghost_atoms[index].update_style()

        self._refresh_bounding_rect()
        # find_atom_near() above made the scene cache an empty boundingRect
        self.prepareGeometryChange()

    @staticmethod
    def _signature(
        preview_atoms: List[Dict[str, Any]],
        bonds_info: list[Any],
        existing_h_counts: Dict[int, int],
    ) -> Any:
        """Identify the ghost molecule, ignoring where the cursor put it.

        The fused atoms' hydrogen counts belong here: they come from the editor,
        not from the template, so a ghost reused across cursor moves would keep
        the label of whichever atom it was fused to before.
        """
        return (
            tuple(
                (
                    atom.get("symbol", "C"),
                    atom.get("charge", 0),
                    atom.get("radical", 0),
                )
                for atom in preview_atoms
            ),
            tuple(tuple(bond) for bond in bonds_info),
            tuple(sorted(existing_h_counts.items())),
        )

    def _move_ghost(self, points: list[QPointF]) -> None:
        """Move an unchanged ghost molecule to the new cursor position."""
        for atom, point in zip(self.ghost_atoms, points):
            atom.setPos(point)
        for bond in self.ghost_bonds:
            bond.update_position()
        for bond_index, ring in self.ring_atoms.items():
            positions = [self.ghost_atoms[i].pos() for i in ring]
            self.ghost_bonds[bond_index].ring_center = (
                sum(p.x() for p in positions) / len(positions),
                sum(p.y() for p in positions) / len(positions),
            )

    def _occupied_cells(self) -> Dict[Any, List[AtomItem]]:
        """Bucket the editor's atoms by position.

        One pass per preview update; asking the scene per template atom instead
        makes every vertex of a big template re-walk the whole scene.
        """
        scene: Any = self.scene()
        atom_items = getattr(scene, "atom_items", None)
        if not atom_items:
            return {}
        cells: Dict[Any, List[AtomItem]] = {}
        for atom in atom_items.values():
            try:
                pos = atom.pos()
            except (AttributeError, RuntimeError, TypeError):
                continue
            key = (int(pos.x() // SNAP_CELL), int(pos.y() // SNAP_CELL))
            cells.setdefault(key, []).append(atom)
        return cells

    @staticmethod
    def _existing_atom_at(
        point: QPointF, cells: Dict[Any, List[AtomItem]]
    ) -> Optional[AtomItem]:
        """Return the editor atom sitting on point, if the preview snapped onto one."""
        cell_x, cell_y = int(point.x() // SNAP_CELL), int(point.y() // SNAP_CELL)
        best: Optional[AtomItem] = None
        best_distance = SNAP_TOLERANCE
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for atom in cells.get((cell_x + dx, cell_y + dy), ()):
                    pos = atom.pos()
                    distance = math.hypot(pos.x() - point.x(), pos.y() - point.y())
                    if distance <= best_distance:
                        best, best_distance = atom, distance
        return best

    # ------------------------------------------------------------------
    # Painting
    # ------------------------------------------------------------------

    def _refresh_bounding_rect(self) -> None:
        """Recompute the cached rect; Qt asks for it far too often to do it live."""
        rect = QRectF()
        for atom in self.ghost_atoms:
            rect = rect.united(atom.boundingRect().translated(atom.pos()))
        for bond in self.ghost_bonds:
            rect = rect.united(bond.boundingRect().translated(bond.pos()))
        self.cached_rect = (
            rect.adjusted(-10, -10, 10, 10)
            if not rect.isNull()
            else self.polygon.boundingRect().adjusted(-5, -5, 5, 5)
        )

    def boundingRect(self) -> QRectF:
        """Return the bounding rect encompassing the preview geometry."""
        if self.is_chain and self.chain_points:
            xs = [p.x() for p in self.chain_points] + [self.chain_label_pos.x()]
            ys = [p.y() for p in self.chain_points] + [self.chain_label_pos.y()]
            return QRectF(
                min(xs) - 10,
                min(ys) - 20,
                max(xs) - min(xs) + 90,
                max(ys) - min(ys) + 40,
            )
        return QRectF(self.cached_rect)

    def paint(
        self,
        painter: Optional[QPainter],
        option: Optional[QStyleOptionGraphicsItem],
        widget: Optional[QWidget] = None,
    ) -> None:
        """Dispatch to the chain or ghost-molecule painter."""
        if painter is None:
            return
        if self.is_chain:
            self.paint_chain(painter)
        else:
            self.paint_ghost(painter, option)

    def paint_chain(self, painter: QPainter) -> None:
        """Paint the ghost zigzag chain plus the repeat count next to the cursor."""
        if len(self.chain_points) < 2:
            return

        painter.setPen(QPen(GHOST_COLOR, 2.5))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for i in range(len(self.chain_points) - 1):
            painter.drawLine(QLineF(self.chain_points[i], self.chain_points[i + 1]))

        dot = QColor(80, 80, 80, 140)
        painter.setBrush(QBrush(dot))
        painter.setPen(QPen(dot, 1))
        for p in self.chain_points[1:]:
            painter.drawEllipse(p, 2.5, 2.5)

        if not self.chain_label:
            return
        font = QFont("Arial", 11, QFont.Weight.Bold)
        painter.setFont(font)
        metrics = painter.fontMetrics()
        rect = QRectF(metrics.boundingRect(self.chain_label))
        rect.moveTopLeft(self.chain_label_pos)
        rect = rect.adjusted(-5, -3, 5, 3)
        painter.setPen(QPen(GHOST_COLOR, 1))
        painter.setBrush(QBrush(QColor(255, 255, 255, 225)))
        painter.drawRoundedRect(rect, 4, 4)
        painter.setPen(QPen(QColor(40, 40, 40)))
        painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, self.chain_label)

    def paint_regular_template(self, painter: QPainter) -> None:
        """Paint a ring template; rings and user templates share one renderer."""
        self.paint_ghost(painter)

    def paint_user_template(self, painter: QPainter) -> None:
        """Paint a user template; rings and user templates share one renderer."""
        self.paint_ghost(painter)

    def paint_ghost(
        self, painter: QPainter, option: Optional[QStyleOptionGraphicsItem] = None
    ) -> None:
        """Paint the ghost molecule using the real atom and bond renderers."""
        if not self.ghost_atoms:
            return

        self._paint_label_backings(painter)

        painter.save()
        painter.setOpacity(GHOST_OPACITY)
        for bond in self.ghost_bonds:
            painter.save()
            painter.translate(bond.pos())
            bond.paint(painter, option)  # type: ignore[arg-type]
            painter.restore()
        for i, atom in enumerate(self.ghost_atoms):
            if i in self.existing_indices:
                # Already drawn by the editor itself; do not double-paint it
                continue
            painter.save()
            painter.translate(atom.pos())
            atom.paint(painter, option)  # type: ignore[arg-type]
            painter.restore()
        painter.restore()

    def _paint_label_backings(self, painter: QPainter) -> None:
        """Cover the label the first atom replaces and mark that atom.

        Ghost labels themselves need no backing: the bonds are shortened around
        them, exactly as in the editor.
        """
        painter.save()
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(LABEL_BACKING))
        if not self.replaced_label_path.isEmpty():
            painter.drawPath(self.replaced_label_path)

        if self.mark_first_atom and self.ghost_atoms:
            first = self.ghost_atoms[0]
            path = first.get_bg_ellipse_path()
            if path.isEmpty():
                rect = QRectF(
                    -FIRST_ATOM_RADIUS,
                    -FIRST_ATOM_RADIUS,
                    FIRST_ATOM_RADIUS * 2,
                    FIRST_ATOM_RADIUS * 2,
                )
            else:
                rect = path.boundingRect().adjusted(-3, -3, 3, 3)
            painter.setPen(QPen(FIRST_ATOM_BORDER, 1.5))
            painter.setBrush(QBrush(FIRST_ATOM_FILL))
            painter.drawEllipse(rect.translated(first.pos()))
        painter.restore()
