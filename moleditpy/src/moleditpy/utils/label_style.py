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

import re
from collections.abc import Mapping
from typing import Any, Dict, Tuple

# (kind, display name, default text color, base font size). The defaults
# reproduce how each 3D label looked before it could be styled. Grouped as
# LABEL_SECTIONS below lists them.
LABEL_KINDS: Tuple[Tuple[str, str, str, int], ...] = (
    ("index", "Index", "#003366", 18),
    ("original_id", "Original ID", "#009000", 18),
    ("xyz_index", "XYZ Index", "#8B0000", 18),
    ("symbol", "Element Symbol", "#000000", 18),
    ("coords", "Coordinates", "#000000", 18),
    ("chiral", "Chiral (R/S)", "#0000FF", 20),
    ("ez", "E/Z", "#006400", 18),
    ("selection", "Selection", "#FFFF00", 12),
    ("measurement", "Measurement", "#FF0000", 16),
    ("constraint", "Constraint", "#00FFFF", 12),
)

# Settings-tab subsections of the label colors: (title, kinds in display order).
LABEL_SECTIONS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("Atom Info", ("index", "original_id", "xyz_index", "symbol", "coords")),
    ("Stereo", ("chiral", "ez")),
    ("Tool", ("selection", "measurement", "constraint")),
)

BACKGROUND_COLOR_KEY = "label_background_color_3d"
BACKGROUND_OPACITY_KEY = "label_background_opacity_3d"
FONT_SIZE_RANGE = (6, 72)
FONT_FAMILY_KEY = "label_font_family_3d"
FONT_BOLD_KEY = "label_font_bold_3d"
FONT_ITALIC_KEY = "label_font_italic_3d"
# The fonts VTK renders labels with; (setting value, display name).
FONT_FAMILIES: Tuple[Tuple[str, str], ...] = (
    ("arial", "Arial"),
    ("courier", "Courier"),
    ("times", "Times"),
)


def color_key(kind: str) -> str:
    """Settings key of one label kind's text color."""
    return f"label_color_{kind}_3d"


def size_key(kind: str) -> str:
    """Settings key of one label kind's font size."""
    return f"label_font_size_{kind}_3d"


DEFAULT_LABEL_SETTINGS: Dict[str, Any] = {
    **{color_key(kind): color for kind, _, color, _ in LABEL_KINDS},
    **{size_key(kind): size for kind, _, _, size in LABEL_KINDS},
    BACKGROUND_COLOR_KEY: "#808080",
    BACKGROUND_OPACITY_KEY: 0.5,
    FONT_FAMILY_KEY: "arial",
    FONT_BOLD_KEY: True,
    FONT_ITALIC_KEY: False,
}

_HEX_COLOR = re.compile(r"^#[0-9A-Fa-f]{6}$")


def _get(settings: Any, key: str) -> Any:
    """Read *key* from a settings mapping; anything else yields None."""
    return settings.get(key) if isinstance(settings, Mapping) else None


def _color(settings: Any, key: str) -> str:
    """A stored color if it is a valid #RRGGBB, else the default."""
    value = _get(settings, key)
    if isinstance(value, str) and _HEX_COLOR.match(value):
        return value
    return str(DEFAULT_LABEL_SETTINGS[key])


def _number(settings: Any, key: str, low: float, high: float) -> float:
    """A stored number clamped to [low, high], else the default."""
    value = _get(settings, key)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return min(max(float(value), low), high)
    return float(DEFAULT_LABEL_SETTINGS[key])


def _flag(settings: Any, key: str) -> bool:
    """A stored boolean, else the default."""
    value = _get(settings, key)
    return value if isinstance(value, bool) else bool(DEFAULT_LABEL_SETTINGS[key])


def _family(settings: Any) -> str:
    """A stored VTK font family, else the default."""
    value = _get(settings, FONT_FAMILY_KEY)
    if value in {family for family, _ in FONT_FAMILIES}:
        return str(value)
    return str(DEFAULT_LABEL_SETTINGS[FONT_FAMILY_KEY])


def label_kwargs(settings: Any, kind: str) -> Dict[str, Any]:
    """PyVista ``add_point_labels`` style arguments for one label kind.

    *settings* may be missing, partial, hand-edited or not a mapping at all
    (test fakes); every value falls back to its default independently.
    """
    return {
        "text_color": _color(settings, color_key(kind)),
        "shape_color": _color(settings, BACKGROUND_COLOR_KEY),
        "shape_opacity": _number(settings, BACKGROUND_OPACITY_KEY, 0.0, 1.0),
        "font_size": int(round(_number(settings, size_key(kind), *FONT_SIZE_RANGE))),
        "font_family": _family(settings),
        "bold": _flag(settings, FONT_BOLD_KEY),
        "italic": _flag(settings, FONT_ITALIC_KEY),
    }
