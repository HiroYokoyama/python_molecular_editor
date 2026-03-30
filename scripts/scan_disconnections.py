#!/usr/bin/env python3
"""
scan_disconnections.py  --  AST-based scanner for broken self.host.X connections.

Finds:
  • self.host.METHOD() calls where METHOD is not on MainWindow nor on any manager
  • self.host.ATTR accesses where ATTR appears to be missing
  • hasattr(self.host, "X") / hasattr(self, "X") checks — marks "soft" connections
  • self.host.MANAGER.METHOD() where MANAGER or METHOD is undefined

Usage:
    python scripts/scan_disconnections.py          # report only
    python scripts/scan_disconnections.py --fix    # auto-fix single-occurrence issues
    python scripts/scan_disconnections.py --verbose  # include hasattr/soft hits too
"""

import ast
import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "moleditpy" / "src" / "moleditpy" / "ui"

# ---------------------------------------------------------------------------
# Manager attribute name on MainWindow → Python source file
# ---------------------------------------------------------------------------
MANAGER_FILES = {
    "state_manager":             SRC / "app_state.py",
    "edit_actions_manager":      SRC / "edit_actions_logic.py",
    "ui_manager":                SRC / "ui_manager.py",
    "view_3d_manager":           SRC / "view_3d_logic.py",
    "edit_3d_manager":           SRC / "edit_3d_logic.py",
    "dialog_manager":            SRC / "dialog_logic.py",
    "compute_manager":           SRC / "compute_logic.py",
    "export_manager":            SRC / "export_logic.py",
    "io_manager":                SRC / "io_logic.py",
    "string_importer_manager":   SRC / "string_importers.py",
    "init_manager":              SRC / "main_window_init.py",
}

# Files to scan for disconnections (Now expanded to check all logic and dialog files)
SCAN_FILES = [f for f in SRC.glob("*.py") if f.name not in ["__init__.py", "molecular_parsers.py", "sip_isdeleted_safe.py", "project_io.py"]]

# MainWindow's direct attributes (not delegated to managers)
MAINWINDOW_DIRECT_ATTRS = {
    # Qt base methods
    "statusBar", "menuBar", "centralWidget", "close", "show", "hide",
    "update", "repaint", "resize", "move", "setWindowTitle", "windowTitle",
    "geometry", "rect", "width", "height", "setEnabled", "setVisible",
    "setStyleSheet", "setWindowIcon", "raise_", "activateWindow",
    "closeEvent", "resizeEvent", "keyPressEvent", "mousePressEvent",
    "focusInEvent", "focusOutEvent", "wheelEvent",
    # MainWindow data attributes
    "data", "scene", "view_2d", "view_3d", "plotter",
    "settings", "current_mol", "current_file_path", "has_unsaved_changes",
    "is_2d_editable", "constraints_3d", "is_xyz_derived",
    "active_worker_ids", "halt_ids", "_ih_update_counter",
    "chem_check_tried", "chem_check_failed",
    "dragged_atom_info", "_saved_state",
    # UI elements on MainWindow
    "undo_action", "redo_action", "cut_action", "copy_action", "paste_action",
    "edit_3d_action", "measurement_action", "analysis_action",
    "optimize_3d_button", "atom_type_group",
    "left_panel", "right_panel", "splitter",
    "toolbar", "toolbar_bottom",
    # UI buttons/widgets set on host in main_window_init.py
    "formula_label", "cleanup_button", "convert_button", "export_button",
    "tool_group", "mode_actions", "style_button", "plugin_manager",
    "plugin_toolbar", "import_menu", "_template_dialog",
    # Per-manager attrs stored on host
    "opt3d_actions", "opt3d_method_labels", "settings_dirty",
    "undo_stack", "redo_stack",
    # Unified UI actions on host
    "other_atom_action", "redraw_menu_action", "show_atom_id_action",
    "show_atom_coords_action", "show_atom_symbol_action",
    "intermolecular_rdkit_action",
    # Manager attrs (correct routing)
    "state_manager", "edit_actions_manager", "ui_manager", "view_3d_manager",
    "edit_3d_manager", "dialog_manager", "compute_manager", "export_manager",
    "io_manager", "string_importer_manager", "init_manager",
}

# Methods that IOManager intentionally exposes as delegation wrappers
DELEGATION_METHODS = {
    "create_json_data", "load_from_json_data", "set_state_from_data",
    "get_current_state", "update_window_title", "reset_undo_stack",
    "restore_ui_for_editing", "clear_all", "fit_to_view", "check_unsaved_changes",
}

# ---------------------------------------------------------------------------
# Step 1 — Build method maps via AST
# ---------------------------------------------------------------------------

def ast_collect_methods(fpath: Path) -> set:
    """Return set of method names (not dunder) defined at class level in fpath."""
    try:
        tree = ast.parse(fpath.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError:
        return set()
    methods = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.ClassDef,)):
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and not item.name.startswith("__"):
                    methods.add(item.name)
    return methods


def build_full_method_map() -> dict:
    """Build {method_name: manager_attr} across all manager files.
    IOManager delegation methods are excluded (they are intentional bridges)."""
    full_map: dict = {}
    for mgr_attr, fpath in MANAGER_FILES.items():
        if not fpath.exists():
            continue
        for name in ast_collect_methods(fpath):
            if name not in DELEGATION_METHODS:
                full_map.setdefault(name, mgr_attr)
    return full_map


def build_per_manager_methods() -> dict:
    """Return {manager_attr: set(method_names)}."""
    result = {}
    for mgr_attr, fpath in MANAGER_FILES.items():
        result[mgr_attr] = ast_collect_methods(fpath) if fpath.exists() else set()
    return result


# ---------------------------------------------------------------------------
# Step 2 — Regex-based scan (handles both dot access and hasattr)
# ---------------------------------------------------------------------------

SELF_HOST_ATTR_RE = re.compile(
    r"\bself\.host\."               # self.host.
    r"([a-zA-Z_]\w*)"               # ATTR
)

HASATTR_RE = re.compile(
    r'\bhasattr\s*\(\s*self(?:\.host)?\s*,\s*["\']([a-zA-Z_]\w*)["\']'
)

SELF_HOST_MANAGER_METHOD_RE = re.compile(
    r"\bself\.host\.([a-zA-Z_]\w+)\.([a-zA-Z_]\w+)"
)


def scan_file(fpath: Path, full_method_map: dict, per_mgr: dict,
              verbose: bool) -> list:
    """Return list of finding dicts for fpath."""
    if not fpath.exists():
        return []

    text   = fpath.read_text(encoding="utf-8", errors="replace")
    lines  = text.splitlines()

    # Which methods are locally defined in THIS file?
    local_methods = ast_collect_methods(fpath)

    # Which manager does THIS file belong to?
    this_mgr_attr = next(
        (k for k, v in MANAGER_FILES.items() if v == fpath), None
    )
    this_mgr_methods = per_mgr.get(this_mgr_attr, set()) if this_mgr_attr else set()

    findings = []

    for lineno, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue

        # ---- A) self.host.MANAGER.METHOD() — check both parts exist ----
        for m in SELF_HOST_MANAGER_METHOD_RE.finditer(line):
            mgr_attr   = m.group(1)
            method     = m.group(2)
            # Skip dunder attrs like __dict__, __class__, etc.
            if mgr_attr.startswith("__"):
                continue
            # Skip if first part is a direct MainWindow attribute (e.g. data.atoms)
            if mgr_attr in MAINWINDOW_DIRECT_ATTRS:
                continue
            if mgr_attr not in MANAGER_FILES:
                findings.append({
                    "kind":    "UNKNOWN_MANAGER",
                    "lineno":  lineno,
                    "attr":    mgr_attr,
                    "extra":   method,
                    "line":    stripped,
                    "soft":    False,
                })
                continue
            known_methods = per_mgr.get(mgr_attr, set())
            if method not in known_methods and not method.startswith("_"):
                findings.append({
                    "kind":    "MISSING_MANAGER_METHOD",
                    "lineno":  lineno,
                    "attr":    f"{mgr_attr}.{method}",
                    "extra":   "",
                    "line":    stripped,
                    "soft":    False,
                })

        # ---- B) self.host.ATTR (single level) — check if attr is known ----
        for m in SELF_HOST_ATTR_RE.finditer(line):
            attr = m.group(1)
            # Skip if this is actually "self.host.MANAGER.METHOD" (handled above)
            rest_after = line[m.end():]
            if rest_after.startswith("."):
                continue
            # Skip manager attrs, direct attrs, delegation methods, dunder
            if attr in MAINWINDOW_DIRECT_ATTRS:
                continue
            if attr in DELEGATION_METHODS:
                continue
            if attr.startswith("__"):
                continue
            # Skip if it's a local method on THIS manager (self-call, not host-call)
            # e.g. self.create_json_data() is an IOManager method
            if attr in local_methods or attr in this_mgr_methods:
                continue
            # If it's a method belonging to another manager → suggest routing
            if attr in full_method_map:
                mgr = full_method_map[attr]
                findings.append({
                    "kind":    "SHOULD_ROUTE_VIA_MANAGER",
                    "lineno":  lineno,
                    "attr":    attr,
                    "extra":   mgr,
                    "line":    stripped,
                    "soft":    False,
                    "old":     f"self.host.{attr}" if "self.host." in line[max(0, m.start()-5):m.end()+10] else f"self.{attr}",
                    "new":     f"self.host.{mgr}.{attr}",
                    "count":   text.count(f"self.host.{attr}"),
                })
            # else: unknown attr — skip unless verbose (could be a Qt attribute etc.)

        # ---- C) hasattr(self[.host], "ATTR") — soft connections ----
        if verbose:
            for m in HASATTR_RE.finditer(line):
                attr = m.group(1)
                # Report if attr is a manager method but accessed without routing
                if attr in full_method_map:
                    mgr = full_method_map[attr]
                    findings.append({
                        "kind":   "HASATTR_SOFT",
                        "lineno": lineno,
                        "attr":   attr,
                        "extra":  f"(belongs to {mgr})",
                        "line":   stripped,
                        "soft":   True,
                    })

    return findings


# ---------------------------------------------------------------------------
# Step 3 — Auto-fix
# ---------------------------------------------------------------------------

def apply_fix(fpath: Path, old: str, new: str) -> bool:
    text = fpath.read_text(encoding="utf-8", errors="replace")
    if old not in text:
        return False
    fpath.write_text(text.replace(old, new), encoding="utf-8")
    return True


# ---------------------------------------------------------------------------
# Entry-point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fix",     action="store_true",
                    help="Auto-apply single-occurrence SHOULD_ROUTE fixes")
    ap.add_argument("--verbose", action="store_true",
                    help="Also show hasattr/soft connection checks")
    args = ap.parse_args()

    print("Building method map from AST...")
    full_method_map = build_full_method_map()
    per_mgr         = build_per_manager_methods()
    print(f"  Known manager methods: {len(full_method_map)}")

    total_hard   = 0
    total_soft   = 0
    auto_fixed   = 0

    for fpath in SCAN_FILES:
        findings = scan_file(fpath, full_method_map, per_mgr, args.verbose)
        hard = [f for f in findings if not f["soft"]]
        soft = [f for f in findings if f["soft"]]

        if not hard and not (args.verbose and soft):
            continue

        rel = fpath.relative_to(ROOT)
        print(f"\n{'='*72}")
        print(f"  {rel}")
        print(f"{'='*72}")

        for f in hard:
            total_hard += 1
            kind = f["kind"]
            print(f"  [{kind}] line {f['lineno']}")
            if kind == "SHOULD_ROUTE_VIA_MANAGER":
                safe = "(single)" if f.get("count", 0) == 1 else f"(occurs {f.get('count',0)}x)"
                print(f"    Found:   {f['old']}")
                print(f"    Fix to:  {f['new']}  {safe}")
                if args.fix and f.get("count", 0) == 1:
                    if apply_fix(fpath, f["old"], f["new"]):
                        print(f"    [FIXED]")
                        auto_fixed += 1
            elif kind == "MISSING_MANAGER_METHOD":
                print(f"    self.host.{f['attr']}()  -- method not found in manager")
            elif kind == "UNKNOWN_MANAGER":
                print(f"    self.host.{f['attr']}.{f['extra']}  -- unknown manager attr '{f['attr']}'")
            print(f"    Context: {f['line'][:90]}")

        if args.verbose and soft:
            for f in soft:
                total_soft += 1
                print(f"  [SOFT/hasattr] line {f['lineno']}  {f['attr']}  {f['extra']}")
                print(f"    Context: {f['line'][:90]}")

    print(f"\n{'='*72}")
    print(f"Hard disconnections found: {total_hard}")
    if args.verbose:
        print(f"Soft (hasattr) checks:     {total_soft}")
    if args.fix:
        print(f"Auto-fixed:                {auto_fixed}")
        remaining = total_hard - auto_fixed
        if remaining:
            print(f"Remaining (manual):        {remaining}")
    else:
        print("Run with --fix to auto-apply single-occurrence fixes.")


if __name__ == "__main__":
    main()
