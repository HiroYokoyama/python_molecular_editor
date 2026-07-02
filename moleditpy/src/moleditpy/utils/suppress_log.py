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

import contextlib
import logging
from typing import Any, Iterator


@contextlib.contextmanager
def suppress_log(*exceptions: Any, note: str = "") -> Iterator[None]:
    """Like contextlib.suppress, but logs the suppressed exception at DEBUG level.

    Drop-in replacement for best-effort blocks: same suppression semantics,
    but never silent — the full traceback lands in the debug log so hidden
    failures (e.g. calls to renamed methods) remain diagnosable.
    """
    try:
        yield
    except exceptions:
        logging.debug("Suppressed%s", f" ({note})" if note else "", exc_info=True)
