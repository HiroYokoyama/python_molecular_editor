#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: LGPL-3.0-or-later
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import os

# qtpy/pyvistaqt must use the same binding as the application's native widgets,
# including when an older installation left QT_API=pyqt6 in the environment.
os.environ["QT_API"] = "pyside6"

from .utils.constants import VERSION as __version__  # noqa: E402, F401
