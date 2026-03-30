#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
moleditpy.ui.molecular_parsers -- backward-compat stub.

The MolecularParsers functionality was merged into IOManager (io_logic.py).
Import from there:
    from moleditpy.ui.io_logic import IOManager
"""
from moleditpy.ui.io_logic import IOManager

# Backward-compat alias
IOManager = IOManager

__all__ = ["IOManager", "IOManager"]
