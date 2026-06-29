#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import importlib.util  # noqa: F401
from .utils.constants import VERSION as __version__  # noqa: F401

try:
    OBABEL_AVAILABLE = importlib.util.find_spec("openbabel") is not None
except ImportError:
    OBABEL_AVAILABLE = False
