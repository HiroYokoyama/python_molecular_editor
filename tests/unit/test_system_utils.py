"""
Unit tests for utils/system_utils.py.

Covers:
  - detect_system_theme: Windows dark/light via winreg, macOS always light,
    Linux via gsettings color-scheme and gtk-theme, unknown platform → None
  - Exception suppression (contextlib.suppress)
  - detect_system_dark_mode: wraps detect_system_theme correctly
"""

import sys
import types
import subprocess
from unittest.mock import MagicMock, patch, call

import pytest

# ---------------------------------------------------------------------------
# Path setup happens via conftest; just import directly.
# ---------------------------------------------------------------------------

from moleditpy.utils.system_utils import detect_system_theme, detect_system_dark_mode


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fake_winreg(val):
    """Return a fake winreg module that yields *val* for AppsUseLightTheme."""
    mod = types.ModuleType("winreg")
    mod.HKEY_CURRENT_USER = 0

    class _Key:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

    mod.OpenKey = MagicMock(return_value=_Key())
    mod.QueryValueEx = MagicMock(return_value=(val, 4))  # (value, type)
    return mod


# ---------------------------------------------------------------------------
# Windows
# ---------------------------------------------------------------------------


class TestWindowsTheme:
    def test_windows_dark_when_val_zero(self):
        fake_reg = _fake_winreg(0)
        with (
            patch("platform.system", return_value="Windows"),
            patch("moleditpy.utils.system_utils.winreg", fake_reg),
        ):
            assert detect_system_theme() == "dark"

    def test_windows_light_when_val_one(self):
        fake_reg = _fake_winreg(1)
        with (
            patch("platform.system", return_value="Windows"),
            patch("moleditpy.utils.system_utils.winreg", fake_reg),
        ):
            assert detect_system_theme() == "light"

    def test_windows_none_when_winreg_unavailable(self):
        """If winreg is None (non-Windows build), Windows branch is skipped."""
        with (
            patch("platform.system", return_value="Windows"),
            patch("moleditpy.utils.system_utils.winreg", None),
        ):
            assert detect_system_theme() is None

    def test_windows_oserror_returns_none(self):
        """OSError in winreg must be suppressed and return None."""
        fake_reg = types.ModuleType("winreg")
        fake_reg.HKEY_CURRENT_USER = 0
        fake_reg.OpenKey = MagicMock(side_effect=OSError("access denied"))
        with (
            patch("platform.system", return_value="Windows"),
            patch("moleditpy.utils.system_utils.winreg", fake_reg),
        ):
            assert detect_system_theme() is None


# ---------------------------------------------------------------------------
# macOS
# ---------------------------------------------------------------------------


class TestMacOSTheme:
    def test_macos_always_light(self):
        with patch("platform.system", return_value="Darwin"):
            assert detect_system_theme() == "light"


# ---------------------------------------------------------------------------
# Linux
# ---------------------------------------------------------------------------


def _proc(returncode, stdout):
    p = MagicMock()
    p.returncode = returncode
    p.stdout = stdout
    return p


class TestLinuxTheme:
    def test_linux_gnome_color_scheme_dark(self):
        responses = [_proc(0, "'prefer-dark'\n"), _proc(1, "")]
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=responses),
        ):
            assert detect_system_theme() == "dark"

    def test_linux_gnome_color_scheme_light(self):
        responses = [_proc(0, "'prefer-light'\n"), _proc(1, "")]
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=responses),
        ):
            assert detect_system_theme() == "light"

    def test_linux_gtk_theme_dark_fallback(self):
        """When color-scheme returns nothing useful, fall back to gtk-theme name."""
        responses = [_proc(0, "'default'\n"), _proc(0, "Adwaita-dark")]
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=responses),
        ):
            assert detect_system_theme() == "dark"

    def test_linux_gtk_theme_no_dark_keyword(self):
        """gtk-theme without '-dark' should not return dark."""
        responses = [_proc(0, "'default'\n"), _proc(0, "Adwaita")]
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=responses),
        ):
            assert detect_system_theme() is None

    def test_linux_gsettings_not_found_returns_none(self):
        """FileNotFoundError (gsettings not installed) is suppressed."""
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=FileNotFoundError),
        ):
            assert detect_system_theme() is None

    def test_linux_gsettings_both_fail_returns_none(self):
        responses = [_proc(1, ""), _proc(1, "")]
        with (
            patch("platform.system", return_value="Linux"),
            patch("subprocess.run", side_effect=responses),
        ):
            assert detect_system_theme() is None


# ---------------------------------------------------------------------------
# Unknown platform
# ---------------------------------------------------------------------------


class TestUnknownPlatform:
    def test_unknown_os_returns_none(self):
        with patch("platform.system", return_value="FreeBSD"):
            assert detect_system_theme() is None


# ---------------------------------------------------------------------------
# detect_system_dark_mode
# ---------------------------------------------------------------------------


class TestDetectSystemDarkMode:
    def test_dark_theme_returns_true(self):
        with patch(
            "moleditpy.utils.system_utils.detect_system_theme", return_value="dark"
        ):
            assert detect_system_dark_mode() is True

    def test_light_theme_returns_false(self):
        with patch(
            "moleditpy.utils.system_utils.detect_system_theme", return_value="light"
        ):
            assert detect_system_dark_mode() is False

    def test_none_theme_returns_none(self):
        with patch(
            "moleditpy.utils.system_utils.detect_system_theme", return_value=None
        ):
            assert detect_system_dark_mode() is None
