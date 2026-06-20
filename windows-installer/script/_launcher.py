#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PyInstaller entry-point.  Imports the package by absolute name so that all
# relative imports inside moleditpy_linux resolve correctly when frozen.

print("-----------------------------------------------------")
print("MoleditPy - A Python-based molecular editing software")
print("-----------------------------------------------------\n")

from moleditpy_linux.main import main  # noqa: E402

if __name__ == "__main__":
    main()
