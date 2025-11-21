#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""


print("-----------------------------------------------------")
print("MoleditPy — A Python-based molecular editing software")
print("-----------------------------------------------------\n")

try:
    # Preferred when running as a package: python -m moleditpy
    from .main import main
except Exception:
    # Fallback when running the file directly: python __main__.py
    # This will import the top-level `main` module in the same folder.
    from main import main

# --- Application Execution ---
if __name__ == '__main__':
    main()

