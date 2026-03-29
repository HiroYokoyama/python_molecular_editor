#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
moleditpy.ui.project_io — backward-compat stub.

The ProjectIO functionality was merged into IOManager (io_logic.py).
Import from there:
    from moleditpy.ui.io_logic import IOManager
"""
from moleditpy.ui.io_logic import IOManager

# Backward-compat alias
MainWindowProjectIo = IOManager

__all__ = ["IOManager", "MainWindowProjectIo"]
