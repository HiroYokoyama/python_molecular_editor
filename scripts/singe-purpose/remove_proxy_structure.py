#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
remove_proxy_structure.py
=========================
Removes the backward-compatibility proxy structure from the moleditpy codebase.

What this removes
-----------------
1. The _CompatibilityProxy system in  moleditpy/modules/__init__.py
2. In MainWindow (main_window.py):
   - The five @property proxy accessors that forward to managers
   - _delegate_manager_methods() and its call-sites in __init__
   - The "delegations" backward-compat dict + setattr loop in __init__
   - __getattr__ (dynamic manager delegation)
3. Phase-2 (--phase2): replaces every mixin's `self.METHOD()` call whose
   METHOD lives on a manager with the explicit `self.MANAGER.METHOD()` form,
   then removes __getattr__ from main_window.py.

What is NOT touched
-------------------
- Manager-internal uses of their OWN attributes (e.g. self.current_3d_style
  inside View3DManager — that is a real attribute there, not a proxy).
- Manager.__getattr__ delegation back to host (that is a different pattern).

Usage
-----
    # Preview all changes without writing anything
    python scripts/remove_proxy_structure.py

    # Apply Phase-1 changes (proxy attrs + main_window.py boilerplate)
    python scripts/remove_proxy_structure.py --apply

    # Apply Phase-1 + Phase-2 (manager-method delegation in mixins)
    python scripts/remove_proxy_structure.py --apply --phase2
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "moleditpy" / "src" / "moleditpy"
UI_DIR = SRC_ROOT / "ui"

# Proxy attribute  →  manager attribute on MainWindow
PROXY_ATTRS: Dict[str, str] = {
    "current_3d_style":       "view_3d_manager",
    "atom_info_display_mode": "view_3d_manager",
    "atom_positions_3d":      "view_3d_manager",
    "show_chiral_labels":     "view_3d_manager",
    "is_3d_edit_mode":        "edit_3d_manager",
}

# ── Mixin files ─────────────────────────────────────────────────────────────
# `self` IS the MainWindow instance inside these classes.
# Fix:  self.ATTR  →  self.MANAGER.ATTR
MIXIN_FILES = [
    UI_DIR / "app_state.py",
    UI_DIR / "main_window_init.py",
    UI_DIR / "ui_manager.py",
    # compute_engine.py (MainWindowCompute) is NOT a real mixin of MainWindow —
    # it is the legacy standalone class used only in tests. Phase-2 must not
    # re-qualify its own internal method calls with compute_manager.*.
    UI_DIR / "molecular_parsers.py",
    UI_DIR / "project_io.py",
    UI_DIR / "string_importers.py",
    UI_DIR / "view_loaders.py",
    UI_DIR / "dialog_manager.py",
]

# ── Manager files ────────────────────────────────────────────────────────────
# `self` IS the manager (has `self.host` = MainWindow).
# Proxy attrs NOT owned by that manager are accessed via __getattr__ → host.
# Fix:  self.ATTR  →  self.host.MANAGER.ATTR
#
# Map: manager file → set of attr names it owns (do NOT rewrite these).
MANAGER_OWNED: Dict[Path, set] = {
    UI_DIR / "view_3d_logic.py": {
        "current_3d_style", "atom_positions_3d",
        "atom_info_display_mode", "show_chiral_labels",
    },
    UI_DIR / "edit_3d_logic.py": {
        "is_3d_edit_mode",
    },
    UI_DIR / "edit_actions_logic.py": set(),
    UI_DIR / "export_logic.py":       set(),
    UI_DIR / "compute_logic.py":      set(),
    UI_DIR / "dialog_logic.py":       set(),
}

# ── Dialog / interactor files ────────────────────────────────────────────────
# These hold a reference to MainWindow under different variable names.
# Fix:  VAR.ATTR  →  VAR.MANAGER.ATTR
DIALOG_MAIN_WINDOW_VARS = ["self.main_window"]   # most dialogs
INTERACTOR_MAIN_WINDOW_VARS = ["mw"]             # custom_interactor_style.py


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _write(path: Path, text: str, dry_run: bool) -> None:
    if dry_run:
        return
    path.write_text(text, encoding="utf-8")


def _show_diff(path: Path, original: str, updated: str) -> None:
    import difflib, sys
    if original == updated:
        return
    rel = path.relative_to(REPO_ROOT)
    lines_orig = original.splitlines(keepends=True)
    lines_new  = updated.splitlines(keepends=True)
    unified = list(difflib.unified_diff(
        lines_orig, lines_new,
        fromfile=str(rel), tofile=str(rel),
        n=0,
    ))
    added   = sum(1 for l in unified if l.startswith("+") and not l.startswith("+++"))
    removed = sum(1 for l in unified if l.startswith("-") and not l.startswith("---"))
    out = sys.stdout
    out.write(f"\n  [{rel}]  +{added}/-{removed} lines\n")
    for line in unified[:30]:
        # sanitise for cp932 / other narrow encodings
        safe = line.rstrip("\n").encode("ascii", errors="replace").decode("ascii")
        out.write(f"    {safe}\n")
    if len(unified) > 30:
        out.write(f"    ... ({len(unified) - 30} more diff lines)\n")


# ---------------------------------------------------------------------------
# Phase-1 helpers
# ---------------------------------------------------------------------------

def _fix_proxy_attrs_in_text(
    text: str,
    access_prefix: str,
    owned_attrs: Optional[set] = None,
    host_prefix: Optional[str] = None,
) -> str:
    """
    Replace  <access_prefix>.<ATTR>  with  <access_prefix>.<MANAGER>.<ATTR>

    If *host_prefix* is given, the replacement becomes
    <host_prefix>.<MANAGER>.<ATTR> (used for managers where the right-hand
    side is self.host.view_3d_manager.* rather than self.view_3d_manager.*).

    *owned_attrs*: attrs to skip (the manager owns them directly).
    """
    owned_attrs = owned_attrs or set()
    for attr, mgr in PROXY_ATTRS.items():
        if attr in owned_attrs:
            continue
        old = f"{access_prefix}.{attr}"
        if host_prefix:
            new = f"{host_prefix}.{mgr}.{attr}"
        else:
            new = f"{access_prefix}.{mgr}.{attr}"
        text = text.replace(old, new)
    return text


def fix_dialog_files(dry_run: bool) -> None:
    """Fix all dialog/ui files that hold a MainWindow reference."""
    skip_files = {
        UI_DIR / "main_window.py",          # handled separately
        UI_DIR / "view_3d_logic.py",        # manager file
        UI_DIR / "edit_3d_logic.py",        # manager file
        UI_DIR / "edit_actions_logic.py",   # manager file
        UI_DIR / "export_logic.py",         # manager file
        UI_DIR / "compute_logic.py",        # manager file
        UI_DIR / "dialog_logic.py",         # manager file
        *MIXIN_FILES,                       # handled by fix_mixin_files
    }
    print("\n=== Phase-1 : dialog / interactor files ===")
    all_py = list(UI_DIR.glob("*.py")) + list((UI_DIR / "settings_tabs").glob("*.py"))
    for fpath in sorted(all_py):
        if fpath in skip_files:
            continue
        original = _read(fpath)
        text = original

        # self.main_window.ATTR  →  self.main_window.MANAGER.ATTR
        for prefix in DIALOG_MAIN_WINDOW_VARS:
            text = _fix_proxy_attrs_in_text(text, prefix)

        # mw.ATTR  →  mw.MANAGER.ATTR  (custom_interactor_style.py etc.)
        if fpath.name == "custom_interactor_style.py":
            for prefix in INTERACTOR_MAIN_WINDOW_VARS:
                text = _fix_proxy_attrs_in_text(text, prefix)

        # getattr(self.main_window, "atom_positions_3d", ...)
        #  → getattr(self.main_window.view_3d_manager, "atom_positions_3d", ...)
        for attr, mgr in PROXY_ATTRS.items():
            text = text.replace(
                f'getattr(self.main_window, "{attr}"',
                f'getattr(self.main_window.{mgr}, "{attr}"',
            )
            text = text.replace(
                f"getattr(self.main_window, '{attr}'",
                f"getattr(self.main_window.{mgr}, '{attr}'",
            )
            # hasattr(self.main_window, "atom_positions_3d")
            #  → self.main_window.view_3d_manager.atom_positions_3d is not None
            # (only for atom_positions_3d since it's the one guarded this way)
            text = text.replace(
                f'not hasattr(self.main_window, "{attr}")',
                f'self.main_window.{mgr}.{attr} is None',
            )
            text = text.replace(
                f'hasattr(self.main_window, "{attr}")',
                f'self.main_window.{mgr}.{attr} is not None',
            )

        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)


def fix_mixin_files(dry_run: bool) -> None:
    """Fix mixin files where self is the MainWindow."""
    print("\n=== Phase-1 : mixin files ===")
    for fpath in MIXIN_FILES:
        if not fpath.exists():
            continue
        original = _read(fpath)
        text = _fix_proxy_attrs_in_text(text=original, access_prefix="self")
        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)


def fix_manager_files(dry_run: bool) -> None:
    """Fix manager files where proxy attrs go through self.__getattr__ → self.host."""
    print("\n=== Phase-1 : manager files ===")
    for fpath, owned in MANAGER_OWNED.items():
        if not fpath.exists():
            continue
        original = _read(fpath)
        # In a manager, not-owned attrs accessed as self.X come from __getattr__
        # which delegates to self.host.X (the proxy on MainWindow).
        # After removal, rewrite to self.host.MANAGER.X
        text = _fix_proxy_attrs_in_text(
            text=original,
            access_prefix="self",
            owned_attrs=owned,
            host_prefix="self.host",
        )
        # compute_logic.py already uses self.host.is_3d_edit_mode (explicit host ref)
        # fix that too
        for attr, mgr in PROXY_ATTRS.items():
            if attr not in owned:
                text = text.replace(
                    f"self.host.{attr}",
                    f"self.host.{mgr}.{attr}",
                )
        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)


# ---------------------------------------------------------------------------
# main_window.py structural surgery
# ---------------------------------------------------------------------------

_MW_DELEGATE_CALLS = """\
        self._delegate_manager_methods(self.export_manager)
        self._delegate_manager_methods(self.view_3d_manager)
        self._delegate_manager_methods(self.edit_3d_manager)
        self._delegate_manager_methods(self.edit_actions_manager)
        self._delegate_manager_methods(self.compute_manager)
        self._delegate_manager_methods(self.dialog_manager)

        # Initialize features via Mixins"""

_MW_DELEGATE_CALLS_NEW = """\
        # Initialize features via Mixins"""

_MW_DELEGATIONS_BLOCK = """\
        # Backwards compatibility for legacy plugins using delegation attributes
        delegations = {
            "main_window_app_state": self,
            "main_window_compute": self,
            "main_window_dialog_manager": self,
            "main_window_edit_3d": self.edit_3d_manager,
            "main_window_edit_actions": self.edit_actions_manager,
            "main_window_export": self.export_manager,
            "main_window_main_init": self,
            "main_window_molecular_parsers": self,
            "main_window_project_io": self,
            "main_window_string_importers": self,
            "main_window_ui_manager": self,
            "main_window_view_3d": self.view_3d_manager,
            "main_window_view_loaders": self,
        }
        for attr, target in delegations.items():
            setattr(self, attr, target)"""

_MW_DELEGATE_METHOD = """\
    def _delegate_manager_methods(self, manager):
        \"\"\"Bind public methods from manager to this instance for backwards compatibility.\"\"\"
        for name in dir(manager):
            # Also delegate _compute_h_counts explicitly for internal use by other mixins
            if name == "_compute_h_counts":
                setattr(self, name, getattr(manager, name))
            elif not name.startswith("_") and callable(getattr(manager, name)):
                if not hasattr(self, name):
                    setattr(self, name, getattr(manager, name))"""

_MW_PROXY_PROPERTIES = """\
    # Proxy properties for View3DManager state synchronization
    @property
    def current_3d_style(self):
        return self.view_3d_manager.current_3d_style

    @current_3d_style.setter
    def current_3d_style(self, value):
        self.view_3d_manager.current_3d_style = value

    @property
    def atom_info_display_mode(self):
        return self.view_3d_manager.atom_info_display_mode

    @atom_info_display_mode.setter
    def atom_info_display_mode(self, value):
        self.view_3d_manager.atom_info_display_mode = value

    @property
    def atom_positions_3d(self):
        return self.view_3d_manager.atom_positions_3d

    @atom_positions_3d.setter
    def atom_positions_3d(self, value):
        self.view_3d_manager.atom_positions_3d = value

    @property
    def show_chiral_labels(self):
        return self.view_3d_manager.show_chiral_labels

    @show_chiral_labels.setter
    def show_chiral_labels(self, value):
        self.view_3d_manager.show_chiral_labels = value

    @property
    def is_3d_edit_mode(self):
        return self.edit_3d_manager.is_3d_edit_mode

    @is_3d_edit_mode.setter
    def is_3d_edit_mode(self, value):
        self.edit_3d_manager.is_3d_edit_mode = value"""

_MW_GETATTR = """\
    def __getattr__(self, name):
        \"\"\"Dynamic delegation to managers for backwards compatibility.\"\"\"
        # Prevent recursion for manager attributes themselves
        if name in [
            "export_manager",
            "view_3d_manager",
            "edit_3d_manager",
            "edit_actions_manager",
        ]:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

        # Core attributes that MUST be on host (skip list)
        core_attrs = {
            "settings",
            "plotter",
            "statusBar",
            "current_mol",
            "scene",
            "view_2d",
            "data",
            "plugin_manager",
            "current_file_path",
        }
        if name in core_attrs or name.startswith("_"):
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

        # Try to find the attribute in any of the injected managers
        for mgr_attr in [
            "export_manager",
            "view_3d_manager",
            "edit_3d_manager",
            "edit_actions_manager",
        ]:
            manager = getattr(self, mgr_attr, None)
            if manager:
                # Avoid hasattr() as it triggers mgr.__getattr__
                # Check if name is in manager's class hierarchy or instance dict
                mgr_type = type(manager)
                if (
                    any(name in c.__dict__ for c in mgr_type.mro())
                    or name in manager.__dict__
                ):
                    return getattr(manager, name)

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )"""


def fix_main_window(dry_run: bool, remove_getattr: bool) -> None:
    """Remove proxy boilerplate from main_window.py."""
    fpath = UI_DIR / "main_window.py"
    print("\n=== Phase-1 : main_window.py ===")
    original = _read(fpath)
    text = original

    # 1. Remove _delegate_manager_methods call-sites in __init__
    text = text.replace(_MW_DELEGATE_CALLS, _MW_DELEGATE_CALLS_NEW)

    # 2. Remove delegations backward-compat block
    text = text.replace("\n" + _MW_DELEGATIONS_BLOCK, "")

    # 3. Remove _delegate_manager_methods method definition
    text = text.replace("\n" + _MW_DELEGATE_METHOD, "")

    # 4. Remove proxy @property definitions
    text = text.replace("\n" + _MW_PROXY_PROPERTIES, "")

    # 5. Optionally remove __getattr__  (Phase-2 flag gates this)
    if remove_getattr:
        text = text.replace("\n" + _MW_GETATTR, "")

    _show_diff(fpath, original, text)
    _write(fpath, text, dry_run)


# ---------------------------------------------------------------------------
# modules/__init__.py — replace proxy with direct re-export stub
# ---------------------------------------------------------------------------

_MODULES_STUB = '''\
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
moleditpy.modules - legacy namespace (backward-compat shim removed).

Import from the canonical locations instead, e.g.:
    from moleditpy.core.molecular_data import MolecularData
    from moleditpy.ui.main_window import MainWindow
"""
'''


def fix_modules_init(dry_run: bool) -> None:
    fpath = SRC_ROOT / "modules" / "__init__.py"
    print("\n=== Phase-1 : modules/__init__.py ===")
    original = _read(fpath)
    _show_diff(fpath, original, _MODULES_STUB)
    _write(fpath, _MODULES_STUB, dry_run)


# ---------------------------------------------------------------------------
# Phase-2 : manager-method delegation in mixin files
# ---------------------------------------------------------------------------

def _collect_manager_methods() -> Dict[str, str]:
    """
    Parse each manager file with regex and return
    { method_name: manager_attribute } for every public method.

    Only methods that are DEFINED in the manager class body are returned.
    """
    manager_files: Dict[Path, str] = {
        UI_DIR / "view_3d_logic.py":      "view_3d_manager",
        UI_DIR / "edit_3d_logic.py":      "edit_3d_manager",
        UI_DIR / "edit_actions_logic.py": "edit_actions_manager",
        UI_DIR / "export_logic.py":       "export_manager",
        UI_DIR / "compute_logic.py":      "compute_manager",
        UI_DIR / "dialog_logic.py":       "dialog_manager",
    }
    methods: Dict[str, str] = {}
    for fpath, mgr_attr in manager_files.items():
        if not fpath.exists():
            continue
        src = _read(fpath)
        for m in re.finditer(r"^\s{4}def ([a-z][a-zA-Z0-9_]*)\s*\(self", src, re.MULTILINE):
            name = m.group(1)
            if not name.startswith("_"):
                # Later manager wins if there is a duplicate name — flag it
                if name in methods and methods[name] != mgr_attr:
                    print(f"  [WARN] method '{name}' found in both "
                          f"{methods[name]} and {mgr_attr} — skipping auto-fix")
                    methods[name] = "__AMBIGUOUS__"
                else:
                    methods[name] = mgr_attr
    return {k: v for k, v in methods.items() if v != "__AMBIGUOUS__"}


def _fix_mixin_method_calls(text: str, method_map: Dict[str, str]) -> str:
    """Replace self.METHOD and self.METHOD(...) with self.MANAGER.METHOD in mixin code.

    Handles both call-sites  self.foo(...)  and callback refs  self.foo  (no parens).
    """
    for method, mgr in method_map.items():
        # Match  self.METHOD  but NOT  self.manager.METHOD  or  self.METHOD_something
        pattern = rf"\bself\.(?!{re.escape(mgr)}\.)({re.escape(method)})\b"
        replacement = rf"self.{mgr}.\1"
        text = re.sub(pattern, replacement, text)
    return text


def fix_mixin_method_delegation(dry_run: bool) -> None:
    print("\n=== Phase-2 : manager-method delegation in mixin files ===")
    method_map = _collect_manager_methods()
    print(f"  Collected {len(method_map)} unique public manager methods.")

    for fpath in MIXIN_FILES:
        if not fpath.exists():
            continue
        original = _read(fpath)
        text = _fix_mixin_method_calls(original, method_map)
        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)


# ---------------------------------------------------------------------------
# Entry-point
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply",   action="store_true",
                    help="Write changes to disk (default: dry-run, preview only)")
    ap.add_argument("--phase2",  action="store_true",
                    help="Also fix manager-method delegation in mixin files "
                         "and remove MainWindow.__getattr__")
    args = ap.parse_args()

    dry_run = not args.apply

    if dry_run:
        print("DRY-RUN mode - no files will be written.")
        print("Pass --apply to apply changes.\n")
    else:
        print("APPLY mode - files will be modified.\n")

    # Phase-1
    fix_mixin_files(dry_run)
    fix_manager_files(dry_run)
    fix_dialog_files(dry_run)
    fix_main_window(dry_run, remove_getattr=args.phase2)
    fix_modules_init(dry_run)

    if args.phase2:
        fix_mixin_method_delegation(dry_run)

    print("\n" + ("=" * 60))
    if dry_run:
        print("Preview complete.  Run with --apply to write changes.")
    else:
        print("Done.")
        if not args.phase2:
            print("\nNote: MainWindow.__getattr__ is still in place.")
            print("Run with --apply --phase2 to also remove it and fix")
            print("manager-method delegation in mixin files.")


if __name__ == "__main__":
    main()
