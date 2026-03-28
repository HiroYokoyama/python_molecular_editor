#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

# Import OBABEL_AVAILABLE from the top-level package and re-export it
try:
    from moleditpy_linux import OBABEL_AVAILABLE
except ImportError:
    OBABEL_AVAILABLE = False

# Re-export core UI utilities
from .sip_isdeleted_safe import sip_isdeleted_safe

__all__ = ["OBABEL_AVAILABLE", "sip_isdeleted_safe"]
