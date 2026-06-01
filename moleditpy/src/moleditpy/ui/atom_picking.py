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
import logging

import numpy as np

try:
    from ..utils.constants import VDW_RADII, pt
except ImportError:
    from moleditpy.utils.constants import VDW_RADII, pt


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


def pick_atom_index_from_screen_sequential(
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


def pick_atom_index_from_screen_vectorized(
    view_3d_manager: Any,
    click_pos: tuple[int, int],
    mol: Optional[Any] = None,
    padding_px: float = 8.0,
    min_radius_px: float = 14.0,
    max_radius_px: float = 96.0,
) -> Optional[int]:
    """
    Vectorized atom picking using the camera's projection matrix.
    Eliminates the O(N) Python loop and VTK C++ boundary calls.
    """
    try:
        plotter = view_3d_manager.plotter
        renderer = plotter.renderer
        positions = view_3d_manager.atom_positions_3d
    except (AttributeError, RuntimeError, TypeError):
        return None

    if positions is None or len(positions) == 0:
        return None

    try:
        positions_array = np.asarray(positions, dtype=float)  # Shape: (N, 3)
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

    active_atoms = min(atom_count, len(positions_array))
    if active_atoms == 0:
        return None

    positions_array = positions_array[:active_atoms]

    # 1. Retrieve the View-Projection (Composite) Matrix from the active camera
    try:
        camera = renderer.GetActiveCamera()
        aspect_ratio = renderer.GetTiledAspectRatio()
        # Get the 4x4 composite projection matrix (vtkMatrix4x4)
        vtk_matrix = camera.GetCompositeProjectionTransformMatrix(aspect_ratio, -1, 1)

        # Convert vtkMatrix4x4 to a NumPy 4x4 array
        matrix = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                matrix[i, j] = vtk_matrix.GetElement(i, j)
    except Exception:
        return None

    # 2. Convert all N world coordinates to homogeneous coordinates (N, 4)
    homogeneous_coords = np.hstack([positions_array, np.ones((active_atoms, 1))])

    # 3. Apply Matrix Multiplication: (N, 4) x (4, 4).T -> Clip Space Coordinates
    clip_coords = homogeneous_coords @ matrix.T

    # 4. Perform Perspective Divide -> Normalized Device Coordinates (NDC)
    w = clip_coords[:, 3:4]
    w_copy = np.copy(w)
    w_copy[np.abs(w_copy) < 1e-5] = 1.0
    ndc_coords = clip_coords[:, :3] / w_copy

    # 5. Transform NDC to Screen/Display Coordinates
    try:
        size = renderer.GetSize()  # (width, height)
    except Exception:
        return None

    # VTK display space coordinates: X: [0, W], Y: [0, H]
    display_coords = np.zeros((active_atoms, 2))
    display_coords[:, 0] = (ndc_coords[:, 0] + 1.0) * 0.5 * size[0]
    display_coords[:, 1] = (ndc_coords[:, 1] + 1.0) * 0.5 * size[1]

    # 6. Vectorized Distance and Hit Radius Calculation
    dx = display_coords[:, 0] - click_pos[0]
    dy = display_coords[:, 1] - click_pos[1]
    distances = np.hypot(dx, dy)

    try:
        if camera.GetParallelProjection():
            pixel_scale = size[1] / (2.0 * camera.GetParallelScale())
        else:
            view_angle_rad = np.radians(camera.GetViewAngle())
            pixel_scale = size[1] / (
                2.0 * np.abs(w.flatten()) * np.tan(view_angle_rad / 2.0)
            )
    except Exception:
        pixel_scale = 20.0  # Safe fallback scale

    # Pre-calculate world radii for all atoms
    world_radii = np.array(
        [_atom_world_radius(view_3d_manager, mol, idx) for idx in range(active_atoms)]
    )

    projected_radii = world_radii * pixel_scale
    hit_radii = np.maximum(
        min_radius_px, np.minimum(max_radius_px, projected_radii + padding_px)
    )

    # 7. Mask out atoms that are further than their hit radius
    valid_mask = distances <= hit_radii
    if not np.any(valid_mask):
        return None

    # 8. Score calculation and find the best index
    # Tie-breaking logic: (ratio) + (distances * 1e-8)
    scores = (distances / hit_radii) + (distances * 1e-8)

    scores[~valid_mask] = np.inf
    best_idx = int(np.argmin(scores))

    return best_idx if scores[best_idx] != np.inf else None


def pick_atom_index_from_screen(
    view_3d_manager: Any,
    click_pos: tuple[int, int],
    mol: Optional[Any] = None,
    padding_px: float = 8.0,
    min_radius_px: float = 14.0,
    max_radius_px: float = 96.0,
) -> Optional[int]:
    """Return the atom nearest a screen click, trying vectorized first, falling back to sequential."""
    try:
        best_idx = pick_atom_index_from_screen_vectorized(
            view_3d_manager, click_pos, mol, padding_px, min_radius_px, max_radius_px
        )
        if best_idx is not None:
            return best_idx
    except Exception as e:
        logging.debug("Vectorized picking failed, falling back to sequential: %s", e)

    return pick_atom_index_from_screen_sequential(
        view_3d_manager, click_pos, mol, padding_px, min_radius_px, max_radius_px
    )
