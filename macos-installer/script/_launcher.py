#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PyInstaller entry-point for the macOS .app bundle.  Imports the package by
# absolute name so that all relative imports inside moleditpy resolve
# correctly when frozen.  Nothing is printed here: the bundle is windowed, so
# stdout is not attached to a terminal (see macos-installer/macos_installer.md).


if __name__ == "__main__":
    from moleditpy.main import main

    main()
