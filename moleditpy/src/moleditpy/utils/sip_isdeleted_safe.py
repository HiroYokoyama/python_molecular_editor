#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: LGPL-3.0-or-later
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from shiboken6 import isValid


def sip_isdeleted_safe(obj: object) -> bool:
    """Check PySide wrapper lifetime; retain the historical name for plugins."""
    if obj is None:
        return True
    try:
        return not isValid(obj)
    except (AttributeError, RuntimeError, TypeError):
        return False
