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

from typing import Any, Optional

import numpy as np

try:
    from ..utils.constants import VDW_RADII, pt
except ImportError:
    from moleditpy_linux.utils.constants import VDW_RADII, pt


def _world_to_display(renderer: Any, pos: Any) -> Optional[tuple[float, float, float]]:
    try:
        renderer.SetWorldPoint(float(pos[0]), float(pos[1]), float(pos[2]), 1.0)
        renderer.WorldToDisplay()
        display = renderer.GetDisplayPoint()
        return (float(display[0]), float(display[1]), float(display[2]))
    except (AttributeError, RuntimeError, TypeError, ValueError, IndexError):
        return None


def _atom_world_radius(view_3d_manager: Any, mol: Any, atom_idx: int) -> float:
    try:
        atom = mol.GetAtomWithIdx(int(atom_idx))
        symbol = atom.GetSymbol()
    except (AttributeError, RuntimeError, TypeError, ValueError):
        symbol = "C"

    settings = {}
    try:
        settings = view_3d_manager.host.init_manager.settings
    except (AttributeError, RuntimeError, TypeError):
        pass

    style = str(getattr(view_3d_manager, "current_3d_style", "ball_and_stick"))
    style = style.lower().replace(" ", "_")

    if style == "cpk":
        scale = settings.get("cpk_atom_scale", 1.0)
        try:
            radius = pt.GetRvdw(pt.GetAtomicNumber(symbol))
            return float(radius if radius > 0.1 else 1.5) * float(scale)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return 1.5 * float(scale)

    if style == "stick":
        return float(settings.get("stick_bond_radius", 0.15))

    if style == "wireframe":
        return 0.01

    scale = settings.get("ball_stick_atom_scale", 1.0)
    return float(VDW_RADII.get(symbol, 0.4)) * float(scale)


def _projected_radius_px(
    renderer: Any, center: Any, world_radius: float
) -> Optional[float]:
    center_display = _world_to_display(renderer, center)
    if center_display is None:
        return None

    offsets = (
        (world_radius, 0.0, 0.0),
        (0.0, world_radius, 0.0),
        (0.0, 0.0, world_radius),
    )
    radius_px = 0.0
    for offset in offsets:
        edge = (
            float(center[0]) + offset[0],
            float(center[1]) + offset[1],
            float(center[2]) + offset[2],
        )
        edge_display = _world_to_display(renderer, edge)
        if edge_display is None:
            continue
        radius_px = max(
            radius_px,
            float(
                np.hypot(
                    edge_display[0] - center_display[0],
                    edge_display[1] - center_display[1],
                )
            ),
        )

    return radius_px


def pick_atom_index_from_screen(
    view_3d_manager: Any,
    click_pos: tuple[int, int],
    mol: Optional[Any] = None,
    padding_px: float = 8.0,
    min_radius_px: float = 14.0,
    max_radius_px: float = 96.0,
) -> Optional[int]:
    """Return the atom nearest a screen click without invoking VTK cell picking."""
    try:
        plotter = view_3d_manager.plotter
        renderer = plotter.renderer
        positions = view_3d_manager.atom_positions_3d
    except (AttributeError, RuntimeError, TypeError):
        return None

    if positions is None:
        return None

    try:
        positions_array = np.asarray(positions, dtype=float)
    except (TypeError, ValueError):
        return None

    if positions_array.ndim != 2 or positions_array.shape[1] < 3:
        return None

    if mol is None:
        mol = getattr(view_3d_manager, "current_mol", None)

    try:
        atom_count = int(mol.GetNumAtoms()) if mol is not None else len(positions_array)
    except (AttributeError, RuntimeError, TypeError, ValueError):
        atom_count = len(positions_array)

    best_idx: Optional[int] = None
    best_score: Optional[tuple[float, float]] = None

    for atom_idx in range(min(atom_count, len(positions_array))):
        center = positions_array[atom_idx]
        display = _world_to_display(renderer, center)
        if display is None or not np.all(np.isfinite(display[:2])):
            continue

        world_radius = _atom_world_radius(view_3d_manager, mol, atom_idx)
        projected_radius = _projected_radius_px(renderer, center, world_radius)
        hit_radius = max(
            float(min_radius_px),
            min(float(max_radius_px), float(projected_radius or 0.0) + padding_px),
        )
        distance = float(np.hypot(display[0] - click_pos[0], display[1] - click_pos[1]))
        if distance > hit_radius:
            continue

        score = (distance / hit_radius, distance)
        if best_score is None or score < best_score:
            best_idx = atom_idx
            best_score = score

    return best_idx
