#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

try:  # pragma: no cover
    import importlib.util

    OBABEL_AVAILABLE = importlib.util.find_spec("openbabel") is not None
except Exception:  # pragma: no cover
    OBABEL_AVAILABLE = False

try:  # pragma: no cover
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except Exception:  # pragma: no cover
    _sip = None
    _sip_isdeleted = None


def sip_isdeleted_safe(obj) -> bool:
    """Return True if sip reports the given wrapper object as deleted.

    This function is conservative: if SIP isn't available or any error
    occurs while checking, it returns False (i.e. not deleted) so that the
    caller can continue other lightweight guards (like checking scene()).
    """
    try:  # pragma: no cover
        if _sip_isdeleted is None:
            return False
        return bool(_sip_isdeleted(obj))
    except Exception:  # pragma: no cover
        return False
