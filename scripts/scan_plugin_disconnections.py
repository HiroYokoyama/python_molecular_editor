#!/usr/bin/env python3
"""
scan_plugin_disconnections.py — Finds broken MainWindow API calls inside tmp/plugins.

Reuses the manager method maps from scan_disconnections_v2.py and scans every
plugin .py file for patterns that won't work with the V3 manager architecture:

  [PLUGIN_DIRECT_MV]    mw.METHOD or self.mw.METHOD   — METHOD is on a manager
  [PLUGIN_MGR_ATTR]     mw.MANAGER.METHOD              — METHOD not in that manager
  [PLUGIN_BAD_IMPORT]   from moleditpy.modules.*       — old module path
  [PLUGIN_NO_RUN]       has initialize() but no run()  — legacy loader logs error
  [PLUGIN_OLD_CONTEXT]  context.register_menu_action   — old 3-arg style (compat, not error)

Usage:
    python scripts/scan_plugin_disconnections.py
    python scripts/scan_plugin_disconnections.py --fix-report   # also write fix_plugin_api.py candidates
"""

import re
import sys
import argparse
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

ROOT        = Path(__file__).resolve().parent.parent
PLUGIN_DIR  = ROOT / "tmp" / "plugins"
SRC_UI      = ROOT / "moleditpy" / "src" / "moleditpy" / "ui"

# ---------------------------------------------------------------------------
# Re-use manager map builder from scan_disconnections_v2
# ---------------------------------------------------------------------------
sys.path.insert(0, str(ROOT / "scripts"))
from scan_disconnections_v2 import (
    MANAGER_FILES,
    ALL_MANAGER_NAMES,
    QT_METHODS,
    build_method_maps,
    MAIN_WINDOW_FILE,
)

# ---------------------------------------------------------------------------
# Patterns
# ---------------------------------------------------------------------------

# bare  mw.METHOD  or  self.mw.METHOD  (not followed by another dot)
RE_BARE_MW = re.compile(
    r"\b(?:self\.)?mw\.([a-zA-Z_]\w+)\b"
)

# mw.MANAGER.METHOD  (or self.mw.MANAGER.METHOD)
RE_MW_MGR_METHOD = re.compile(
    r"\b(?:self\.)?mw\.([a-zA-Z_]\w+)\.([a-zA-Z_]\w+)\b"
)

# old moleditpy.modules.* import
RE_OLD_IMPORT = re.compile(
    r"from\s+moleditpy\.modules[\.\w]*\s+import"
)

# def initialize / def run at module level
RE_INITIALIZE = re.compile(r"^def initialize\s*\(", re.MULTILINE)
RE_RUN        = re.compile(r"^def run\s*\(", re.MULTILINE)

# old 3-arg register_menu_action: context.register_menu_action(path, "text", callback)
RE_OLD_REGISTER = re.compile(
    r"context\.register_menu_action\s*\(\s*['\"][^'\"]+['\"]\s*,\s*['\"]"
)

# Additional valid attrs directly on MainWindow (dynamic attrs, proxies, Qt widgets etc.)
MW_EXTRA_VALID = {
    # core Qt/window
    "statusBar", "menuBar", "centralWidget", "addDockWidget", "close", "show",
    "setWindowTitle", "update", "resize", "geometry", "isVisible", "isEnabled",
    "findChild", "findChildren", "objectName", "winId", "screen", "style",
    # MainWindow-owned data attrs (set by managers on host)
    "current_mol", "plotter", "scene", "splitter",
    "is_xyz_derived", "plugin_manager", "settings",
    "view_2d", "align_x_action", "start_calculation",
    # proxy methods kept on MainWindow
    "trigger_conversion", "convert_button", "has_unsaved_changes",
    "on_calculation_finished", "push_undo_state",
    "update_window_title", "update_realtime_info", "update_undo_redo_actions",
    "check_chemistry_problems_fallback",
    # manager handles
} | ALL_MANAGER_NAMES | QT_METHODS


def scan_plugin_file(fpath: Path, per_mgr: dict, method_to_mgr: dict) -> list:
    if not fpath.exists():
        return []
    text  = fpath.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    findings = []

    is_entry = fpath.name == "__init__.py" or fpath.parent == PLUGIN_DIR
    has_init = bool(RE_INITIALIZE.search(text))
    has_run  = bool(RE_RUN.search(text))

    # [PLUGIN_NO_RUN] — only for entry-point files
    if is_entry and has_init and not has_run:
        findings.append({
            "kind": "PLUGIN_NO_RUN",
            "lineno": 0,
            "expr": "def run(...) missing",
            "detail": "Add a no-op run(mw) stub to silence legacy loader error",
            "line": "",
        })

    for lineno, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue

        # [PLUGIN_BAD_IMPORT]
        if RE_OLD_IMPORT.search(raw):
            findings.append({
                "kind": "PLUGIN_BAD_IMPORT",
                "lineno": lineno,
                "expr": stripped[:80],
                "detail": "moleditpy.modules.* → moleditpy.ui.* or moleditpy.utils.* or moleditpy.core.*",
                "line": stripped,
            })

        # [PLUGIN_OLD_CONTEXT]
        if RE_OLD_REGISTER.search(raw):
            findings.append({
                "kind": "PLUGIN_OLD_CONTEXT",
                "lineno": lineno,
                "expr": stripped[:80],
                "detail": "Old 3-arg register_menu_action(path, text, cb) — use add_menu_action(path, cb, text)",
                "line": stripped,
            })

        # [PLUGIN_MGR_ATTR] — mw.MANAGER.METHOD
        for m in RE_MW_MGR_METHOD.finditer(raw):
            mgr_name = m.group(1)
            method   = m.group(2)
            if mgr_name not in ALL_MANAGER_NAMES:
                continue  # not a manager, skip (handled by DIRECT_MW below)
            # next char must NOT be another dot (triple chain)
            end = m.end()
            if end < len(raw) and raw[end] == ".":
                continue
            known = per_mgr.get(mgr_name, set())
            if method not in known and not method.startswith("_"):
                findings.append({
                    "kind": "PLUGIN_MGR_ATTR",
                    "lineno": lineno,
                    "expr": m.group(0),
                    "detail": f"'{method}' not found in {mgr_name}",
                    "line": stripped,
                })

        # [PLUGIN_DIRECT_MW] — bare mw.METHOD (no second dot)
        for m in RE_BARE_MW.finditer(raw):
            attr = m.group(1)
            end  = m.end()
            # if followed by a dot → this is mw.MANAGER.X, handled above
            if end < len(raw) and raw[end] == ".":
                continue
            if attr in MW_EXTRA_VALID or attr.startswith("_"):
                continue
            if attr in method_to_mgr:
                mgr = method_to_mgr[attr]
                findings.append({
                    "kind": "PLUGIN_DIRECT_MW",
                    "lineno": lineno,
                    "expr": m.group(0),
                    "detail": f"→ mw.{mgr}.{attr}",
                    "line": stripped,
                })

    return findings


KIND_LABEL = {
    "PLUGIN_DIRECT_MW":   "DIRECT_MW_CALL  ",
    "PLUGIN_MGR_ATTR":    "MGR_ATTR_MISSING",
    "PLUGIN_BAD_IMPORT":  "BAD_IMPORT      ",
    "PLUGIN_NO_RUN":      "NO_RUN_STUB     ",
    "PLUGIN_OLD_CONTEXT": "OLD_CONTEXT_API ",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fix-report", action="store_true",
                    help="Print a machine-readable list of (file, lineno, pattern, fix)")
    args = ap.parse_args()

    # Build manager maps from core source
    manager_files = [v for v in MANAGER_FILES.values() if v.exists()]
    all_scan = manager_files[:]
    if MAIN_WINDOW_FILE.exists():
        all_scan.append(MAIN_WINDOW_FILE)

    print("Building manager method maps from AST...")
    per_mgr, method_to_mgr, _ = build_method_maps(all_scan)
    print(f"Managers indexed: {len(per_mgr)}  |  Methods indexed: {len(method_to_mgr)}\n")

    py_files = sorted(PLUGIN_DIR.rglob("*.py"))
    counter  = 0
    files_with_issues = 0

    for fpath in py_files:
        findings = scan_plugin_file(fpath, per_mgr, method_to_mgr)
        if not findings:
            continue
        files_with_issues += 1
        rel = fpath.relative_to(ROOT)
        print(f"\n{'='*72}\n  {rel}\n{'='*72}")
        for f in findings:
            counter += 1
            label = KIND_LABEL.get(f["kind"], f["kind"])
            loc   = f"line {f['lineno']}" if f["lineno"] else "file-level"
            print(f"  [{label}] {loc}")
            print(f"    Expr:  {f['expr']}")
            print(f"    Fix:   {f['detail']}")
            if f["line"]:
                print(f"    Code:  {f['line'][:88]}")

    print(f"\n{'='*72}")
    print(f"Plugin files scanned:  {len(py_files)}")
    print(f"Files with issues:     {files_with_issues}")
    print(f"Total issues found:    {counter}")


if __name__ == "__main__":
    main()
