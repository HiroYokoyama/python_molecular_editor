#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import math


def finite_float(text: str) -> float:
    """Parse *text* as a finite number, raising ValueError otherwise.

    ``float()`` accepts "nan" and "inf"; typed into a geometry dialog they
    would move atoms to NaN or infinite coordinates. Callers already turn a
    ValueError into an "invalid number" message.
    """
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"not a finite number: {text!r}")
    return value
