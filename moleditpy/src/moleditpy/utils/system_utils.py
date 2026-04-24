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
import platform
import subprocess
from typing import Optional

try:
    import winreg
except ImportError:
    winreg = None


def detect_system_dark_mode() -> Optional[bool]:
    """Return True if the OS prefers dark app theme, False if light, or None if unknown."""
    theme = detect_system_theme()
    if theme == "dark":
        return True
    if theme == "light":
        return False
    return None


def detect_system_theme() -> Optional[str]:
    """Return the OS's preferred theme setting as 'dark', 'light', or None."""
    with contextlib.suppress(AttributeError, RuntimeError, OSError, FileNotFoundError):
        # Windows
        if platform.system() == "Windows" and winreg is not None:
            with winreg.OpenKey(
                winreg.HKEY_CURRENT_USER,
                r"Software\Microsoft\Windows\CurrentVersion\Themes\Personalize",
            ) as k:
                val, _ = winreg.QueryValueEx(k, "AppsUseLightTheme")
                return "dark" if int(val) == 0 else "light"

        # macOS fallback
        if platform.system() == "Darwin":
            return "light"

        # Linux / GNOME
        if platform.system() == "Linux":
            p = subprocess.run(
                ["gsettings", "get", "org.gnome.desktop.interface", "color-scheme"],
                capture_output=True,
                text=True,
            )
            if p.returncode == 0:
                out = p.stdout.strip().strip("'\n ")
                if "dark" in out.lower():
                    return "dark"
                if "light" in out.lower():
                    return "light"

            p = subprocess.run(
                ["gsettings", "get", "org.gnome.desktop.interface", "gtk-theme"],
                capture_output=True,
                text=True,
            )
            if p.returncode == 0 and "-dark" in p.stdout.lower():
                return "dark"

    return None
