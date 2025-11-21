#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""


# Open Babel Python binding (optional; required for fallback)
# Do not import `pybel` at module import time. Only expose the presence
# of Open Babel via `OBABEL_AVAILABLE`. Modules should import `pybel`
# lazily if and when they need it.
try:
    import importlib.util
    OBABEL_AVAILABLE = importlib.util.find_spec("openbabel") is not None
except Exception:
    OBABEL_AVAILABLE = False

# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

def sip_isdeleted_safe(obj) -> bool:
    """Return True if sip reports the given wrapper object as deleted.

    This function is conservative: if SIP isn't available or any error
    occurs while checking, it returns False (i.e. not deleted) so that the
    caller can continue other lightweight guards (like checking scene()).
    """
    try:
        if _sip_isdeleted is None:
            return False
        return bool(_sip_isdeleted(obj))
    except Exception:
        return False

