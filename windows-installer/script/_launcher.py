#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PyInstaller entry-point.  Imports the package by absolute name so that all
# relative imports inside moleditpy resolve correctly when frozen.


if __name__ == "__main__":
    from moleditpy.main import main

    main()
