#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PyInstaller entry-point for the macOS .app bundle.  Imports the package by
# absolute name so that all relative imports inside moleditpy_linux resolve
# correctly when frozen.  Nothing is printed here: the bundle is windowed, so
# stdout is not attached to a terminal (see macos-installer/macos_installer.md).

from moleditpy_linux.main import main

if __name__ == "__main__":
    main()
