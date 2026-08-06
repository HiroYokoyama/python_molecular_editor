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
import math
from typing import Any, List, Optional

from PyQt6.QtCore import QPointF

from ..utils.constants import DEFAULT_BOND_LENGTH

# Half of the 120 deg zigzag angle, measured from the chain axis.
CHAIN_HALF_ANGLE_DEG = 30.0
CHAIN_ANGLE_SNAP_DEG = 15.0
MAX_CHAIN_LENGTH = 100


class ChainMixin:
    """Drag-to-draw tool that grows a fixed-bond-length zigzag alkyl chain."""

    atom_items: Any
    add_molecule_fragment: Any
    current_atom_symbol: Any
    template_preview: Any
    views: Any

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.chain_anchor: Optional[QPointF] = None
        self.chain_start_atom: Optional[Any] = None
        self.chain_points: List[QPointF] = []
        self.chain_active: bool = False
        super().__init__(*args, **kwargs)

    def begin_chain(self, pos: QPointF, start_atom: Optional[Any] = None) -> None:
        """Anchor a new chain drag at ``pos`` (or at ``start_atom`` if given)."""
        self.chain_start_atom = start_atom
        self.chain_anchor = start_atom.pos() if start_atom is not None else QPointF(pos)
        self.chain_points = []
        self.chain_active = True

    def calculate_chain_points(self, anchor: QPointF, cursor: QPointF) -> List[QPointF]:
        """Return the zigzag vertices from ``anchor`` toward ``cursor``.

        The chain axis is snapped to 15 deg steps, every bond keeps the standard
        2D bond length, and the number of bonds follows the drag distance
        projected onto the axis.
        """
        vx = cursor.x() - anchor.x()
        vy = cursor.y() - anchor.y()
        if math.hypot(vx, vy) < 1e-6:
            return []

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

        self.chain_points = self.calculate_chain_points(self.chain_anchor, pos)
        if len(self.chain_points) < 2:
            self.template_preview.hide()
        else:
            self.template_preview.set_chain_geometry(
                self.chain_points,
                f"n = {len(self.chain_points) - 1}",
                QPointF(pos.x() + 16, pos.y() - 12),
            )
            self.template_preview.show()

        if self.views():
            self.views()[0].viewport().update()

    def clear_chain_preview(self) -> None:
        """Hide the ghost chain and drop the drag state."""
        self.chain_anchor = None
        self.chain_start_atom = None
        self.chain_points = []
        self.chain_active = False
        try:
            self.template_preview.hide()
        except (AttributeError, RuntimeError) as e:
            logging.debug(f"Could not hide chain preview: {e}")
        if self.views():
            self.views()[0].viewport().update()

    def commit_chain(self, pos: QPointF) -> bool:
        """Materialize the previewed chain as atoms and bonds. True if anything was added."""
        if self.chain_anchor is None:
            return False

        points = self.calculate_chain_points(self.chain_anchor, pos)
        if len(points) < 2:
            return False

        bonds_info = [(i, i + 1, 1) for i in range(len(points) - 1)]
        existing = [self.chain_start_atom] if self.chain_start_atom is not None else []
        symbol = getattr(self, "current_atom_symbol", "C") or "C"
        self.add_molecule_fragment(
            points, bonds_info, existing_items=existing, symbol=symbol
        )
        return True
