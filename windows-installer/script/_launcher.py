#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PyInstaller entry-point.  Imports the package by absolute name so that all
# relative imports inside moleditpy resolve correctly when frozen.

import sys

if __name__ == "__main__":
    if sys.argv[1:] == ["--openbabel-make3d"]:
        from moleditpy.utils.openbabel_worker import make3d_main

        make3d_main()
    else:
        from moleditpy.main import main

        main()
