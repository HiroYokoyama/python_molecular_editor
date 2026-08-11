#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from typing import Any

DEFAULT_HIT_PX = 14.0


def scene_hit_radius(scene: Any, fallback: float = DEFAULT_HIT_PX) -> float:
    """Return the 2D hit/snap radius in scene units for the current zoom.

    The ``bond_snapping_distance_2d`` setting is in screen pixels, so it is
    divided by the view scale.  AtomItem.shape(), BondItem.shape() and
    ``find_atom_near()`` all resolve their hit zone through this, which is what
    keeps hover highlight, clicking and snapping on the same radius.
    """
    if not scene or not scene.views():
        return fallback

    scale = scene.views()[0].transform().m11()
    if scale <= 0.0:
        scale = 1.0

    hit_px = DEFAULT_HIT_PX
    if hasattr(scene, "get_setting"):
        hit_px = scene.get_setting("bond_snapping_distance_2d", DEFAULT_HIT_PX)
    return float(hit_px / scale)
