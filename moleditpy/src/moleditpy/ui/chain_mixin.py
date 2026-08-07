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
from typing import Any, List, Optional, Tuple

from PyQt6.QtCore import QPointF

from .atom_item import AtomItem
from ..utils.constants import DEFAULT_BOND_LENGTH

# Half of the 120 deg zigzag angle, measured from the chain axis.
CHAIN_HALF_ANGLE_DEG = 30.0
CHAIN_ANGLE_SNAP_DEG = 15.0
MAX_CHAIN_LENGTH = 100
# Reject an end snap that would collapse the tail bond onto its neighbour.
MIN_TAIL_BOND_RATIO = 0.4


class ChainMixin:
    """Drag-to-draw tool that grows a fixed-bond-length zigzag alkyl chain."""

    add_molecule_fragment: Any
    current_atom_symbol: Any
    find_atom_near: Any
    get_setting: Any
    template_preview: Any

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.chain_anchor: Optional[QPointF] = None
        self.chain_start_atom: Optional[AtomItem] = None
        self.chain_end_atom: Optional[AtomItem] = None
        self.chain_points: List[QPointF] = []
        self.chain_active: bool = False
        super().__init__(*args, **kwargs)

    def begin_chain(self, pos: QPointF, start_atom: Optional[AtomItem] = None) -> None:
        """Anchor a new chain drag at ``pos`` (or at ``start_atom`` if given)."""
        self.chain_start_atom = start_atom
        self.chain_end_atom = None
        self.chain_anchor = start_atom.pos() if start_atom is not None else QPointF(pos)
        self.chain_points = []
        self.chain_active = True

    def chain_geometry(
        self, cursor: QPointF
    ) -> Tuple[List[QPointF], Optional[AtomItem]]:
        """Zigzag for this cursor, with the last vertex pulled onto a nearby atom.

        Every other bond keeps the fixed length, so the tail bond absorbs the
        whole difference and may end up longer or shorter. Both ends of the drag
        use the shared 2D bond snapping distance, as bond drawing does.
        """
        if self.chain_anchor is None:
            return [], None

        points = self.calculate_chain_points(self.chain_anchor, cursor)
        if len(points) < 2:
            return points, None

        snap_dist = self.get_setting("bond_snapping_distance_2d", 14.0)
        target = self.find_atom_near(cursor, tol=snap_dist)
        if target is None or target is self.chain_start_atom:
            return points, None

        end = target.pos()
        prev = points[-2]
        tail = math.hypot(end.x() - prev.x(), end.y() - prev.y())
        if tail < MIN_TAIL_BOND_RATIO * DEFAULT_BOND_LENGTH:
            return points, None

        points[-1] = end
        return points, target

    def calculate_chain_points(self, anchor: QPointF, cursor: QPointF) -> List[QPointF]:
        """Return the zigzag vertices from ``anchor`` toward ``cursor``.

        The chain axis is snapped to 15 deg steps, every bond keeps the standard
        2D bond length, and the number of bonds follows the drag distance
        projected onto the axis.
        """
        vx = cursor.x() - anchor.x()
        vy = cursor.y() - anchor.y()

        snap = math.radians(CHAIN_ANGLE_SNAP_DEG)
        axis_angle = round(math.atan2(vy, vx) / snap) * snap
        ax, ay = math.cos(axis_angle), math.sin(axis_angle)

        half = math.radians(CHAIN_HALF_ANGLE_DEG)
        step = DEFAULT_BOND_LENGTH * math.cos(half)
        offset = DEFAULT_BOND_LENGTH * math.sin(half)

        count = int(round((vx * ax + vy * ay) / step))
        if count < 1:
            return []
        count = min(count, MAX_CHAIN_LENGTH)

        # Zigzag bulges toward the side of the axis the cursor drifted to.
        side = -1.0 if (ax * vy - ay * vx) < 0 else 1.0
        px, py = -ay * offset * side, ax * offset * side

        points = [QPointF(anchor)]
        for k in range(1, count + 1):
            wobble = 1.0 if k % 2 else 0.0
            points.append(
                QPointF(
                    anchor.x() + ax * step * k + px * wobble,
                    anchor.y() + ay * step * k + py * wobble,
                )
            )
        return points

    def update_chain_preview(self, pos: QPointF) -> None:
        """Refresh the ghost chain and the repeat-count label near the cursor."""
        if not self.chain_active or self.chain_anchor is None:
            return

        self.chain_points, self.chain_end_atom = self.chain_geometry(pos)
        if len(self.chain_points) < 2:
            self.template_preview.hide()
            return

        self.template_preview.set_chain_geometry(
            self.chain_points,
            f"n = {self.new_atom_count()}",
            QPointF(pos.x() + 16, pos.y() - 12),
        )
        self.template_preview.show()

    def new_atom_count(self) -> int:
        """Atoms the previewed chain would add; reused end atoms do not count."""
        if not self.chain_points:
            return 0
        reused = sum(
            1
            for atom in (self.chain_start_atom, self.chain_end_atom)
            if atom is not None
        )
        return len(self.chain_points) - reused

    def clear_chain_preview(self) -> None:
        """Hide the ghost chain and drop the drag state."""
        self.chain_anchor = None
        self.chain_start_atom = None
        self.chain_end_atom = None
        self.chain_points = []
        self.chain_active = False
        self.template_preview.hide()

    def commit_chain(self, pos: QPointF) -> bool:
        """Materialize the previewed chain as atoms and bonds. True if anything was added."""
        if self.chain_anchor is None:
            return False

        points, end_atom = self.chain_geometry(pos)
        if len(points) < 2:
            return False

        bonds_info = [(i, i + 1, 1) for i in range(len(points) - 1)]
        # Naming both reused atoms fuses them even with template fusing turned off.
        existing = [
            atom for atom in (self.chain_start_atom, end_atom) if atom is not None
        ]
        self.add_molecule_fragment(
            points, bonds_info, existing_items=existing, symbol=self.current_atom_symbol
        )
        return True
