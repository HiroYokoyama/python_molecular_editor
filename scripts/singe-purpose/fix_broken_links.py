#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fix_broken_links.py
===================
Reconnects all broken links after the mixin→manager refactor, without proxies.

What this fixes
---------------
1. Broken self-references in managers:
   - self.host.METHOD() where METHOD is defined IN THE SAME manager class
     → replace with self.METHOD()
2. __getattr__ proxy methods in all manager files:
   - def __getattr__(self, name): return getattr(self.host, name)
     → remove entirely (callers must use explicit self.host.X)
3. Unused 'Any' import in managers after __getattr__ removal
4. Tests: update imports from old mixin names to new manager names
   and add backward-compat aliases where needed

Broken links found (static analysis):
  io_logic.py:178   self.host.fix_mol_block(raw)        → self.fix_mol_block(raw)
  io_logic.py:265   self.host.load_xyz_file(file_path)  → self.load_xyz_file(file_path)
  io_logic.py:296   self.host.fix_mol_block(raw)        → self.fix_mol_block(raw)

__getattr__ proxy in:
  io_logic.py, compute_logic.py, dialog_logic.py, edit_3d_logic.py,
  edit_actions_logic.py, export_logic.py, view_3d_logic.py

Usage
-----
    python scripts/fix_broken_links.py            # dry-run (preview only)
    python scripts/fix_broken_links.py --apply    # write changes
    python scripts/fix_broken_links.py --apply --tests  # also fix test imports
"""

from __future__ import annotations

import argparse
import ast
import difflib
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT  = REPO_ROOT / "moleditpy" / "src" / "moleditpy"
UI_DIR    = SRC_ROOT / "ui"
TESTS_UNIT = REPO_ROOT / "tests" / "unit"

# ---------------------------------------------------------------------------
# Manager files that HAVE __getattr__ proxy (all 7)
# ---------------------------------------------------------------------------
MANAGER_FILES_WITH_PROXY: List[Path] = [
    UI_DIR / "io_logic.py",
    UI_DIR / "compute_logic.py",
    UI_DIR / "dialog_logic.py",
    UI_DIR / "edit_3d_logic.py",
    UI_DIR / "edit_actions_logic.py",
    UI_DIR / "export_logic.py",
    UI_DIR / "view_3d_logic.py",
]

# ALL manager files (for self-reference scan)
ALL_MANAGER_FILES: List[Path] = MANAGER_FILES_WITH_PROXY + [
    UI_DIR / "app_state.py",
    UI_DIR / "string_importers.py",
    UI_DIR / "ui_manager.py",
    UI_DIR / "main_window_init.py",
]

# ---------------------------------------------------------------------------
# Host (MainWindow) attribute names that must be accessed via self.host.X
# in manager files.  Derived from MW_ATTRS in the old refactor script.
# ---------------------------------------------------------------------------
MW_HOST_ATTRS: List[str] = [
    # Data / scene
    "data", "scene", "view_2d", "current_mol", "plotter",
    # Methods that are bound to host
    "statusBar",
    # Settings / state
    "settings", "current_file_path", "has_unsaved_changes", "plugin_manager",
    "settings_file", "settings_dir", "initial_settings", "settings_dirty",
    "initialization_complete", "_ih_update_counter", "_active_calc_threads",
    "_is_restoring_state",
    # UI state flags
    "is_2d_editable", "is_xyz_derived", "chem_check_tried", "chem_check_failed",
    # Actions & buttons
    "measurement_action", "edit_3d_action", "analysis_action",
    "style_button", "optimize_3d_button", "export_button", "cleanup_button",
    "convert_button", "formula_label", "tool_group", "splitter", "mode_actions",
    "undo_action", "redo_action", "cut_action", "copy_action", "paste_action",
    # Layout helpers
    "minimize_2d_panel", "restore_2d_panel",
    # Manager UI methods (on host but delegate to UIManager)
    "update_window_title", "_enter_3d_viewer_ui_mode", "restore_ui_for_editing",
    "_enable_3d_features", "_enable_3d_edit_actions",
    # 3D view state
    "show_chiral_labels",
]

# Per-manager: attributes OWNED by each manager's own __init__ (NOT host).
# These must NOT be prefixed with self.host.
MANAGER_OWNED_ATTRS: Dict[str, Set[str]] = {
    "view_3d_logic.py": {
        "host", "current_3d_style", "atom_positions_3d",
        "atom_info_display_mode", "show_chiral_labels",
        "_camera_initialized", "atom_actor_original_opacity",
    },
    "edit_3d_logic.py": {
        "host", "is_3d_edit_mode", "dragged_atom_info",
        "constraints_3d", "active_3d_dialogs",
    },
    "export_logic.py": {"host"},
    "compute_logic.py": {"host", "_active_calc_threads", "halt_ids"},
    "dialog_logic.py": {"host"},
    "io_logic.py": {"host"},
    "edit_actions_logic.py": {
        "host", "undo_stack", "redo_stack",
    },
    "app_state.py": {"host", "_cls"},
    "string_importers.py": {"host"},
    "ui_manager.py": {"host"},
    "main_window_init.py": {"host"},
}

# Files (non-manager) that access MainWindow via a different variable
# (e.g. self.window, self.main_window, mw)
SCENE_WINDOW_METHODS_REMAP: Dict[str, str] = {
    # molecule_scene.py uses self.window.METHOD() where METHOD is on a manager
    "push_undo_state":              "edit_actions_manager.push_undo_state",
    "update_2d_measurement_labels": "edit_3d_manager.update_2d_measurement_labels",
}

# ---------------------------------------------------------------------------
# Test import aliases: old mixin name → (module_path, new_class_name)
# ---------------------------------------------------------------------------
TEST_IMPORT_FIXES: Dict[str, Tuple[str, str]] = {
    "MainWindowAppState":        ("moleditpy.ui.app_state",       "StateManager"),
    "MainWindowStringImporters": ("moleditpy.ui.string_importers", "StringImporterManager"),
    "MainWindowUiManager":       ("moleditpy.ui.ui_manager",       "UIManager"),
}

# Old modules that no longer exist → redirect to new location
TEST_MODULE_REDIRECTS: Dict[str, str] = {
    "moleditpy.ui.project_io":       "moleditpy.ui.io_logic",
    "moleditpy.ui.molecular_parsers":"moleditpy.ui.io_logic",
    "moleditpy.ui.view_loaders":     "moleditpy.ui.io_logic",
}

# Old class → new class in tests
TEST_CLASS_RENAMES: Dict[str, str] = {
    "MainWindowProjectIo":       "IOManager",
    "MainWindowMolecularParsers":"IOManager",
    "MainWindowViewLoaders":     "IOManager",
    "MainWindowAppState":        "StateManager",
    "MainWindowStringImporters": "StringImporterManager",
    "MainWindowUiManager":       "UIManager",
}

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
    if original == updated:
        return
    rel = path.relative_to(REPO_ROOT)
    lines_orig = original.splitlines(keepends=True)
    lines_new  = updated.splitlines(keepends=True)
    unified = list(difflib.unified_diff(
        lines_orig, lines_new,
        fromfile=str(rel), tofile=str(rel), n=2,
    ))
    added   = sum(1 for l in unified if l.startswith("+") and not l.startswith("+++"))
    removed = sum(1 for l in unified if l.startswith("-") and not l.startswith("---"))
    print(f"\n  [{rel}]  +{added}/-{removed} lines")
    for line in unified[:40]:
        safe = line.rstrip("\n").encode("ascii", errors="replace").decode("ascii")
        print(f"    {safe}")
    if len(unified) > 40:
        print(f"    ... ({len(unified) - 40} more diff lines)")


# ---------------------------------------------------------------------------
# Step 1: Collect all public methods defined in each manager class
# ---------------------------------------------------------------------------

def collect_methods_in_file(fpath: Path) -> Set[str]:
    """Return set of method names defined at class-body indent level."""
    methods: Set[str] = set()
    try:
        src = _read(fpath)
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                for item in node.body:
                    if isinstance(item, ast.FunctionDef):
                        methods.add(item.name)
    except SyntaxError:
        # Fallback: regex
        for m in re.finditer(r"^\s{4}def ([a-zA-Z_]\w*)\s*\(self", _read(fpath), re.MULTILINE):
            methods.add(m.group(1))
    return methods


# ---------------------------------------------------------------------------
# Step 2: Fix broken self.host.METHOD() → self.METHOD()
# ---------------------------------------------------------------------------

def fix_broken_self_refs(fpath: Path, dry_run: bool) -> None:
    """
    In manager files, replace self.host.METHOD(...) with self.METHOD(...)
    when METHOD is defined in the same manager class.
    """
    local_methods = collect_methods_in_file(fpath)
    # Remove dunder and private methods from consideration
    local_public = {m for m in local_methods if not m.startswith("_")}

    original = _read(fpath)
    text = original

    fixed = []
    for method in sorted(local_public):
        # Match  self.host.method(  or  self.host.method  (no parens)
        pattern = rf"\bself\.host\.({re.escape(method)})\b"
        new_text = re.sub(pattern, rf"self.\1", text)
        if new_text != text:
            fixed.append(method)
            text = new_text

    if fixed:
        print(f"  Broken self-refs fixed in {fpath.name}: {fixed}")
        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)
    else:
        print(f"  No broken self-refs in {fpath.name}")


# ---------------------------------------------------------------------------
# Step 3: Remove __getattr__ proxy methods
# ---------------------------------------------------------------------------

def remove_getattr_proxy(fpath: Path, dry_run: bool) -> None:
    """Remove the __getattr__ delegation method (simple or guarded form)."""
    original = _read(fpath)
    text = original

    # Matches the full method body up to (but not including) the next method/class/blank.
    # Handles two forms:
    #   Form A (simple):
    #       def __getattr__(self, name[: str]) [-> Any]:
    #           return getattr(self.host, name)
    #
    #   Form B (guarded):
    #       def __getattr__(self, name[: str]) [-> Any]:
    #           """..."""
    #           if name == "host":
    #               raise AttributeError(name)
    #           return getattr(self.host, name)
    pattern = (
        r"\n    def __getattr__\(self, name[^)]*\)[^:]*:\n"       # def line
        r"(?:        \"\"\"[^\"]*\"\"\"\n)?"                       # optional docstring (single line)
        r"(?:        if name == \"host\":\n"                       # optional guard
        r"            raise AttributeError\(name\)\n)?"
        r"        return getattr\(self\.host, name\)\n"            # delegation
    )

    new_text = re.sub(pattern, "\n", text)

    # Clean up triple blank lines
    new_text = re.sub(r"\n{3,}", "\n\n", new_text)

    if new_text != original:
        print(f"  Removed __getattr__ proxy from {fpath.name}")
        _show_diff(fpath, original, new_text)
        _write(fpath, new_text, dry_run)
    else:
        print(f"  No proxy __getattr__ pattern matched in {fpath.name} (may need manual check)")


# ---------------------------------------------------------------------------
# Step 4: Remove 'Any' from typing import if it was only used by __getattr__
# ---------------------------------------------------------------------------

def cleanup_any_import(fpath: Path, dry_run: bool) -> None:
    """
    Remove 'Any' from 'from typing import ...' if 'Any' no longer appears in file.
    """
    text = _read(fpath)
    # Count occurrences of ' Any' as a type annotation (excluding import line)
    lines = text.splitlines()
    import_line_idx = None
    for i, line in enumerate(lines):
        if re.match(r"from typing import .*\bAny\b", line):
            import_line_idx = i
            break

    if import_line_idx is None:
        return

    # Count uses of Any outside the import line
    rest = "\n".join(l for i, l in enumerate(lines) if i != import_line_idx)
    if re.search(r"\bAny\b", rest):
        return  # Still used elsewhere

    # Remove Any from the import
    original = text
    old_line = lines[import_line_idx]
    new_line = re.sub(r",\s*Any\b", "", old_line)
    new_line = re.sub(r"\bAny\b,?\s*", "", new_line)
    new_line = re.sub(r",\s*$", "", new_line)
    new_line = re.sub(r"from typing import\s*$", "", new_line)
    # If import line is now empty, remove it
    if re.match(r"\s*$", new_line):
        lines.pop(import_line_idx)
    else:
        lines[import_line_idx] = new_line

    new_text = "\n".join(lines)
    if new_text != original:
        print(f"  Cleaned 'Any' import from {fpath.name}")
        _write(fpath, new_text, dry_run)


# ---------------------------------------------------------------------------
# Step 5: Fix test imports
# ---------------------------------------------------------------------------

def fix_test_file(fpath: Path, dry_run: bool) -> None:
    """
    Fix a test file:
    - Redirect imports from removed modules to new locations
    - Rename old mixin class references to new manager class names
    """
    original = _read(fpath)
    text = original

    # Fix module paths in import statements
    for old_mod, new_mod in TEST_MODULE_REDIRECTS.items():
        text = text.replace(f"from {old_mod} import", f"from {new_mod} import")

    # Rename class names throughout the file
    for old_cls, new_cls in TEST_CLASS_RENAMES.items():
        # Only replace as whole words (not partial matches)
        text = re.sub(rf"\b{re.escape(old_cls)}\b", new_cls, text)

    if text != original:
        print(f"  Fixed test imports in {fpath.name}")
        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)


def add_compat_aliases(dry_run: bool) -> None:
    """
    Add backward-compat aliases at the bottom of source files so old test
    imports still resolve even without rewriting the tests.
    e.g.  MainWindowAppState = StateManager
    """
    aliases = {
        UI_DIR / "app_state.py": [
            ("StateManager", "MainWindowAppState"),
        ],
        UI_DIR / "string_importers.py": [
            ("StringImporterManager", "MainWindowStringImporters"),
        ],
        UI_DIR / "ui_manager.py": [
            ("UIManager", "MainWindowUiManager"),
        ],
        UI_DIR / "io_logic.py": [
            ("IOManager", "MainWindowProjectIo"),
            ("IOManager", "MainWindowMolecularParsers"),
            ("IOManager", "MainWindowViewLoaders"),
        ],
    }

    print("\n=== Adding backward-compat aliases ===")
    for fpath, pairs in aliases.items():
        original = _read(fpath)
        text = original
        appended = []
        for new_cls, old_name in pairs:
            alias_line = f"{old_name} = {new_cls}"
            if alias_line not in text:
                appended.append(alias_line)

        if appended:
            # Check if there's already an alias section
            if "# Backward-compat aliases" not in text:
                text = text.rstrip("\n") + "\n\n\n# Backward-compat aliases\n"
            for line in appended:
                text += line + "\n"
            print(f"  Added aliases to {fpath.name}: {[p[1] for p in pairs]}")
            _show_diff(fpath, original, text)
            _write(fpath, text, dry_run)
        else:
            print(f"  Aliases already present in {fpath.name}")


# ---------------------------------------------------------------------------
# Step 6: Scan for any remaining self.host.METHOD where METHOD is on a
#         DIFFERENT manager (cross-manager proxy leak)
# ---------------------------------------------------------------------------

def build_all_method_maps() -> Dict[str, str]:
    """Build {method_name: manager_attr} across all managers."""
    # manager attr → file
    mgr_files = {
        "state_manager":          UI_DIR / "app_state.py",
        "string_importer_manager":UI_DIR / "string_importers.py",
        "ui_manager":             UI_DIR / "ui_manager.py",
        "init_manager":           UI_DIR / "main_window_init.py",
        "io_manager":             UI_DIR / "io_logic.py",
        "compute_manager":        UI_DIR / "compute_logic.py",
        "dialog_manager":         UI_DIR / "dialog_logic.py",
        "edit_3d_manager":        UI_DIR / "edit_3d_logic.py",
        "edit_actions_manager":   UI_DIR / "edit_actions_logic.py",
        "export_manager":         UI_DIR / "export_logic.py",
        "view_3d_manager":        UI_DIR / "view_3d_logic.py",
    }
    method_map: Dict[str, str] = {}
    for mgr_attr, fpath in mgr_files.items():
        if not fpath.exists():
            continue
        methods = collect_methods_in_file(fpath)
        for m in methods:
            if not m.startswith("_"):
                method_map[m] = mgr_attr
    return method_map


def fix_host_attrs_comprehensive(dry_run: bool) -> None:
    """
    For each manager file that previously relied on __getattr__ to delegate
    self.X → self.host.X, replace every remaining self.ATTR with self.host.ATTR
    for all attrs in MW_HOST_ATTRS that the manager does NOT own itself.

    This is the main fix for the 'data', 'scene', 'plotter', 'current_mol',
    'settings', etc. accesses still using bare self.X after __getattr__ removal.
    """
    print("\n=== Step 4: Comprehensive host-attr fix (self.X -> self.host.X) ===")
    for fpath in MANAGER_FILES_WITH_PROXY:
        if not fpath.exists():
            continue
        owned = MANAGER_OWNED_ATTRS.get(fpath.name, set())
        original = text = _read(fpath)
        fixed_attrs = []
        for attr in MW_HOST_ATTRS:
            if attr in owned:
                continue
            # Match  self.ATTR  but NOT  self.host.ATTR  or  self.something_else.ATTR
            # Also NOT when ATTR is being defined as a function parameter or docstring word
            pattern = (
                r"(?<!\w)"                    # not preceded by word char
                r"\bself\."                   # self.
                r"(?!host\.)"                 # not already self.host.
                rf"({re.escape(attr)})"       # the attribute name
                r"\b"                         # word boundary
                r"(?!\w)"                     # not followed by word char (handled by \b)
            )
            new_text = re.sub(pattern, r"self.host.\1", text)
            if new_text != text:
                fixed_attrs.append(attr)
                text = new_text
        if text != original:
            count = sum(text.count(f"self.host.{a}") - original.count(f"self.host.{a}")
                        for a in fixed_attrs)
            print(f"  {fpath.name}: fixed {count} refs to {fixed_attrs}")
            _show_diff(fpath, original, text)
            _write(fpath, text, dry_run)
        else:
            print(f"  {fpath.name}: no host-attr fixes needed")


def fix_scene_window_calls(dry_run: bool) -> None:
    """
    Fix molecule_scene.py calls: self.window.METHOD() → self.window.MANAGER.METHOD()
    where METHOD was previously on MainWindow via __getattr__ delegation but now lives
    on a specific manager.
    """
    scene_file = UI_DIR / "molecule_scene.py"
    if not scene_file.exists():
        return
    print("\n=== Fixing molecule_scene.py window.METHOD() calls ===")
    original = text = _read(scene_file)
    for method, new_target in SCENE_WINDOW_METHODS_REMAP.items():
        old = f"self.window.{method}("
        new = f"self.window.{new_target}("
        count = text.count(old)
        if count:
            text = text.replace(old, new)
            print(f"  Fixed {count}x self.window.{method}() -> self.window.{new_target}()")
    if text != original:
        _show_diff(scene_file, original, text)
        _write(scene_file, text, dry_run)


def fix_host_attr_accesses(dry_run: bool) -> None:
    """
    Fix self.HOST_ATTR → self.host.HOST_ATTR in manager files.

    After __getattr__ removal, any manager method that accesses a host-side
    attribute (QAction, QPushButton, host state flags, etc.) via bare 'self.X'
    will raise AttributeError.  We fix each known case explicitly.

    HOST_ATTRS are identified from main_window_init.py where they are set as
    self.host.X = ...  (not self.X = ...).
    """

    # (file, old_pattern, new_replacement)
    # Use str.replace for literal strings.
    FIXES: List[Tuple[Path, str, str]] = [
        # ── main_window_init.py ──────────────────────────────────────────────
        # analysis_action is set on init_manager (self) instead of host
        (UI_DIR / "main_window_init.py",
         "        self.analysis_action = QAction(",
         "        self.host.analysis_action = QAction("),
        (UI_DIR / "main_window_init.py",
         "        self.analysis_action.triggered.connect(",
         "        self.host.analysis_action.triggered.connect("),
        (UI_DIR / "main_window_init.py",
         "        self.analysis_action.setEnabled(",
         "        self.host.analysis_action.setEnabled("),
        (UI_DIR / "main_window_init.py",
         "        analysis_menu.addAction(self.analysis_action)",
         "        analysis_menu.addAction(self.host.analysis_action)"),

        # ── app_state.py ─────────────────────────────────────────────────────
        (UI_DIR / "app_state.py",
         "            self.analysis_action.setEnabled(",
         "            self.host.analysis_action.setEnabled("),

        # ── string_importers.py ──────────────────────────────────────────────
        (UI_DIR / "string_importers.py",
         "self.analysis_action.setEnabled(",
         "self.host.analysis_action.setEnabled("),

        # ── edit_actions_logic.py ────────────────────────────────────────────
        (UI_DIR / "edit_actions_logic.py",
         "            self.cut_action.setEnabled(",
         "            self.host.cut_action.setEnabled("),
        (UI_DIR / "edit_actions_logic.py",
         "            self.copy_action.setEnabled(",
         "            self.host.copy_action.setEnabled("),
        (UI_DIR / "edit_actions_logic.py",
         "            self.paste_action.setEnabled(",
         "            self.host.paste_action.setEnabled("),
        (UI_DIR / "edit_actions_logic.py",
         "            self.edit_3d_action.setChecked(",
         "            self.host.edit_3d_action.setChecked("),
        (UI_DIR / "edit_actions_logic.py",
         "                    self.optimize_3d_button.setEnabled(",
         "                    self.host.optimize_3d_button.setEnabled("),
        # chem_check / xyz state flags live on host (set during init)
        (UI_DIR / "edit_actions_logic.py",
         "        self.chem_check_tried = False",
         "        self.host.chem_check_tried = False"),
        (UI_DIR / "edit_actions_logic.py",
         "        self.chem_check_failed = False",
         "        self.host.chem_check_failed = False"),
        (UI_DIR / "edit_actions_logic.py",
         "            self.chem_check_tried = True",
         "            self.host.chem_check_tried = True"),
        (UI_DIR / "edit_actions_logic.py",
         "            self.chem_check_failed = False",
         "            self.host.chem_check_failed = False"),
        (UI_DIR / "edit_actions_logic.py",
         "            self.chem_check_failed = True",
         "            self.host.chem_check_failed = True"),
        (UI_DIR / "edit_actions_logic.py",
         "        self.is_xyz_derived = False",
         "        self.host.is_xyz_derived = False"),
        # optimize_3d_button accesses without leading spaces (line 1450)
        (UI_DIR / "edit_actions_logic.py",
         "                self.optimize_3d_button.setEnabled(",
         "                self.host.optimize_3d_button.setEnabled("),

        # ── view_3d_logic.py ─────────────────────────────────────────────────
        (UI_DIR / "view_3d_logic.py",
         "            self.measurement_action.setChecked(",
         "            self.host.measurement_action.setChecked("),
        (UI_DIR / "view_3d_logic.py",
         "            self.edit_3d_action.setChecked(",
         "            self.host.edit_3d_action.setChecked("),
    ]

    print("\n=== Fixing host-attribute accesses in manager files ===")
    changed: Dict[str, int] = {}
    for fpath, old, new in FIXES:
        if not fpath.exists():
            continue
        original = _read(fpath)
        count = original.count(old)
        if count:
            text = original.replace(old, new)
            changed.setdefault(fpath.name, 0)
            changed[fpath.name] += count
            _write(fpath, text, dry_run)
            if dry_run:
                _show_diff(fpath, original, text)
            print(f"  Fixed {count}x in {fpath.name}: '{old.strip()}'")

    if not changed:
        print("  All host-attr accesses already fixed.")


def fix_known_cross_manager_leaks(dry_run: bool) -> None:
    """
    Fix known cross-manager proxy leaks found by static analysis.
    These are self.host.METHOD() calls where METHOD belongs to a different manager.
    """
    fixes: List[Tuple[Path, str, str]] = [
        # dialog_logic.py calls remove_dialog_from_list which is on edit_3d_manager
        (
            UI_DIR / "dialog_logic.py",
            "self.host.remove_dialog_from_list(",
            "self.host.edit_3d_manager.remove_dialog_from_list(",
        ),
    ]

    print("\n=== Fixing known cross-manager proxy leaks ===")
    for fpath, old, new in fixes:
        if not fpath.exists():
            continue
        original = _read(fpath)
        if old in original:
            count = original.count(old)
            text = original.replace(old, new)
            print(f"  Fixed {count}x '{old}' -> '{new}' in {fpath.name}")
            _show_diff(fpath, original, text)
            _write(fpath, text, dry_run)
        else:
            print(f"  Already fixed or not found: '{old}' in {fpath.name}")


def scan_host_method_leaks(dry_run: bool) -> None:
    """
    In each manager file, scan for self.host.METHOD() where METHOD lives on
    a DIFFERENT manager.  These should be self.host.MANAGER.METHOD().
    Prints findings; does NOT auto-fix (too risky without context).
    """
    method_map = build_all_method_maps()
    print("\n=== Scanning for cross-manager proxy leaks ===")
    found_any = False
    for fpath in ALL_MANAGER_FILES:
        if not fpath.exists():
            continue
        text = _read(fpath)
        local_methods = collect_methods_in_file(fpath)
        for match in re.finditer(r"\bself\.host\.([a-z][a-zA-Z0-9_]*)\b", text):
            method = match.group(1)
            if method in method_map and method not in local_methods:
                correct_mgr = method_map[method]
                lineno = text[:match.start()].count("\n") + 1
                print(f"  WARN {fpath.name}:{lineno}  "
                      f"self.host.{method}()  →  self.host.{correct_mgr}.{method}()")
                found_any = True
    if not found_any:
        print("  No cross-manager proxy leaks found.")


# ---------------------------------------------------------------------------
# Entry-point
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--apply",  action="store_true",
                    help="Write changes to disk (default: dry-run)")
    ap.add_argument("--tests",  action="store_true",
                    help="Also fix test file imports")
    ap.add_argument("--aliases", action="store_true",
                    help="Add backward-compat aliases instead of rewriting tests")
    args = ap.parse_args()

    dry_run = not args.apply

    if dry_run:
        print("DRY-RUN mode -- no files written.  Pass --apply to apply.\n")
    else:
        print("APPLY mode -- files will be modified.\n")

    # --- Step 1 + 2: Fix broken self-references in all manager files ---
    print("\n=== Step 1: Fix broken self.host.METHOD() references ===")
    for fpath in ALL_MANAGER_FILES:
        if fpath.exists():
            fix_broken_self_refs(fpath, dry_run)

    # --- Step 3: Remove __getattr__ proxy methods ---
    print("\n=== Step 2: Remove __getattr__ proxy methods ===")
    for fpath in MANAGER_FILES_WITH_PROXY:
        if fpath.exists():
            remove_getattr_proxy(fpath, dry_run)

    # --- Step 4: Cleanup Any imports ---
    print("\n=== Step 3: Cleanup unused 'Any' imports ===")
    for fpath in MANAGER_FILES_WITH_PROXY:
        if fpath.exists():
            cleanup_any_import(fpath, dry_run)

    # --- Step 4: Comprehensive host-attr fix (self.X -> self.host.X) ---
    fix_host_attrs_comprehensive(dry_run)

    # --- Step 4a: Fix specific known host-attribute accesses ---
    fix_host_attr_accesses(dry_run)

    # --- Step 4c: Fix molecule_scene.py window method calls ---
    fix_scene_window_calls(dry_run)

    # --- Step 4b: Fix known cross-manager proxy leaks ---
    fix_known_cross_manager_leaks(dry_run)

    # --- Step 5: Scan for remaining cross-manager leaks ---
    scan_host_method_leaks(dry_run)

    # --- Step 6: Fix tests or add aliases ---
    if args.aliases:
        add_compat_aliases(dry_run)

    if args.tests:
        print("\n=== Step 4: Fix test file imports ===")
        test_files_to_fix = [
            TESTS_UNIT / "test_app_state.py",
            TESTS_UNIT / "test_app_state_persistence.py",
            TESTS_UNIT / "test_string_importers.py",
            TESTS_UNIT / "test_ui_manager_robustness.py",
            TESTS_UNIT / "test_io.py",
            TESTS_UNIT / "test_parsers.py",
            TESTS_UNIT / "test_parser_robustness.py",
            TESTS_UNIT / "test_project_io_raw.py",
            TESTS_UNIT / "test_project_io_extended.py",
        ]
        for fpath in test_files_to_fix:
            if fpath.exists():
                fix_test_file(fpath, dry_run)
            else:
                print(f"  (not found) {fpath.name}")

    print("\n" + "=" * 60)
    if dry_run:
        print("Preview complete. Run with --apply to write changes.")
        print("Add --aliases to add compat aliases (easiest test fix).")
        print("Add --tests to rewrite test imports instead.")
    else:
        print("Done.")


if __name__ == "__main__":
    main()
