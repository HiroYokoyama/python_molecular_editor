#!/usr/bin/env python3
"""
scan_disconnections_v2.py — Comprehensive scanner for broken manager connections.

Architecture: MainWindow mixins → separate manager classes.
  MainWindow.view_3d_manager  → View3DManager  (view_3d_logic.py)
  MainWindow.edit_actions_manager → EditActionsManager  (edit_actions_logic.py)
  ... etc.

Finds ALL of these broken patterns across manager files, main_window.py, and dialogs:

  [MISSING_MANAGER_METHOD]   self.host.MANAGER.METHOD  — METHOD not in MANAGER
  [SHOULD_ROUTE_VIA_MANAGER] self.host.METHOD or self.METHOD — METHOD moved to a manager
  [MAINWINDOW_SELF_METHOD]   self.METHOD in main_window.py — METHOD now on a manager
  [DIALOG_DIRECT_CALL]       self.main_window.METHOD  — METHOD now on a manager

Usage:
    python scripts/scan_disconnections_v2.py
    python scripts/scan_disconnections_v2.py --verbose   # also show hasattr checks
    python scripts/scan_disconnections_v2.py --dialogs   # also scan dialog files
"""

import ast
import re
import sys
import argparse
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8')

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "moleditpy" / "src" / "moleditpy" / "ui"

# ---------------------------------------------------------------------------
# Manager attribute name on MainWindow → source file
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

MAIN_WINDOW_FILE  = SRC / "main_window.py"
ALL_MANAGER_NAMES = set(MANAGER_FILES.keys())

# ---------------------------------------------------------------------------
# Attributes that are VALID to access directly on MainWindow (self.host.X)
# These are pure Qt base class methods + MainWindow-owned data attrs.
# Do NOT put manager-owned methods here.
# ---------------------------------------------------------------------------
QT_METHODS = {
    # QMainWindow / QWidget
    "statusBar", "menuBar", "centralWidget", "close", "show", "hide",
    "update", "repaint", "resize", "move", "setWindowTitle", "windowTitle",
    "geometry", "rect", "width", "height", "setEnabled", "setVisible",
    "setStyleSheet", "setWindowIcon", "raise_", "activateWindow",
    "closeEvent", "resizeEvent", "keyPressEvent", "mousePressEvent",
    "focusInEvent", "focusOutEvent", "wheelEvent",
    "parent", "layout", "setLayout", "sizePolicy", "setSizePolicy",
    "minimumSize", "maximumSize", "setMinimumSize", "setMaximumSize",
    "size", "pos", "x", "y", "isActiveWindow", "isMinimized",
    "isMaximized", "isFullScreen", "isHidden", "isVisible", "isEnabled",
    "hasFocus", "setFocus", "clearFocus", "addToolBar", "removeToolBar",
    "insertToolBar", "addDockWidget", "removeDockWidget", "setCentralWidget",
    "setMenuBar", "setStatusBar", "sender", "connect", "disconnect",
    "emit", "blockSignals", "signalsBlocked", "deleteLater", "findChild",
    "findChildren", "objectName", "setObjectName", "inherits", "style",
    "setStyle", "cursor", "setCursor", "font", "setFont", "palette",
    "setPalette", "toolTip", "setToolTip", "window", "showMaximized",
    "showMinimized", "showNormal", "showFullScreen", "timerEvent", "metaObject",
    "childAt", "setAcceptDrops", "setGeometry", "addToolBarBreak",
    "installEventFilter", "removeEventFilter", "eventFilter",
    "mapToGlobal", "mapFromGlobal", "mapToParent", "mapFromParent",
    "nativeParentWidget", "effectiveWinId", "winId", "screen",
    "grab", "grabKeyboard", "releaseKeyboard", "grabMouse", "releaseMouse",
    "setWindowFlags", "windowFlags", "setWindowState", "windowState",
    "setAttribute", "testAttribute", "clearMask", "setMask", "mask",
    "adjustSize", "updateGeometry", "setFixedSize", "setFixedWidth", "setFixedHeight",
    # Python internals
    "__dict__", "__class__", "__module__",
    # App-level methods still on MainWindow itself
    "start_calculation", "_is_restoring_state", "dragEnterEvent", "dropEvent",
    "check_unsaved_changes", "warning_message_box",
    # Manager handles — valid to reference directly
    "host",
} | ALL_MANAGER_NAMES

# ---------------------------------------------------------------------------
# Step 1 — AST-based name collection
# ---------------------------------------------------------------------------

def ast_collect_names(fpath: Path) -> set:
    """Collect method names + self.X assignments (including __init__)."""
    try:
        tree = ast.parse(fpath.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError:
        return set()
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if not isinstance(item, ast.FunctionDef):
                    continue
                if not item.name.startswith("__"):
                    names.add(item.name)
                # Collect self.X = ... in ALL methods (incl. __init__)
                for sub in ast.walk(item):
                    if isinstance(sub, ast.Assign):
                        for tgt in sub.targets:
                            if (isinstance(tgt, ast.Attribute)
                                    and isinstance(tgt.value, ast.Name)
                                    and tgt.value.id == "self"):
                                names.add(tgt.attr)
                    elif isinstance(sub, ast.AnnAssign):
                        if (isinstance(sub.target, ast.Attribute)
                                and isinstance(sub.target.value, ast.Name)
                                and sub.target.value.id == "self"):
                            names.add(sub.target.attr)
    return names


def collect_cross_assigned_attrs(all_files: list) -> dict:
    """
    Scan for 'self.host.MANAGER.ATTR = X' across all files.
    These attrs are valid members of MANAGER even if not in its source.
    Example: self.host.view_3d_manager.current_mol = mol
    """
    RE = re.compile(r"\bself\.host\.([a-zA-Z_]\w+)\.([a-zA-Z_]\w+)\s*=")
    extras: dict = {mgr: set() for mgr in ALL_MANAGER_NAMES}
    for fpath in all_files:
        if not fpath.exists():
            continue
        for line in fpath.read_text(encoding="utf-8", errors="replace").splitlines():
            for m in RE.finditer(line):
                mgr, attr = m.group(1), m.group(2)
                if mgr in extras:
                    extras[mgr].add(attr)
    return extras


def collect_host_direct_attrs(all_files: list) -> set:
    """
    Scan for 'self.host.ATTR = X' (single-level) across all files.
    These are valid MainWindow attributes set by init_manager / other managers.
    Example: self.host.initial_settings = ...
             self.host.align_x_action = ...
    """
    RE_SINGLE = re.compile(r"\bself\.host\.([a-zA-Z_]\w+)\s*=")
    RE_TWO    = re.compile(r"\bself\.host\.[a-zA-Z_]\w+\.[a-zA-Z_]\w+")
    known: set = set()
    for fpath in all_files:
        if not fpath.exists():
            continue
        for line in fpath.read_text(encoding="utf-8", errors="replace").splitlines():
            for m in RE_SINGLE.finditer(line):
                tail = line[m.start():]
                if not RE_TWO.match(tail):
                    known.add(m.group(1))
    return known


def build_method_maps(all_scan_files: list) -> tuple:
    """
    Returns:
      per_mgr         — {manager_name: set_of_names}
      method_to_mgr   — {method_name: manager_name}  (reverse lookup)
      host_direct     — set of valid single-level self.host.X attrs
    """
    cross = collect_cross_assigned_attrs(all_scan_files)
    host_direct = collect_host_direct_attrs(all_scan_files)

    per_mgr: dict = {}
    for mgr, fpath in MANAGER_FILES.items():
        names = ast_collect_names(fpath) if fpath.exists() else set()
        names |= cross.get(mgr, set())
        per_mgr[mgr] = names

    # Reverse map — first manager that defines a name wins
    method_to_mgr: dict = {}
    for mgr, names in per_mgr.items():
        for name in names:
            if name not in method_to_mgr:
                method_to_mgr[name] = mgr

    return per_mgr, method_to_mgr, host_direct

# ---------------------------------------------------------------------------
# Step 2 — Regex patterns (separate, no backtracking ambiguity)
# ---------------------------------------------------------------------------

# self.host.MANAGER.METHOD
RE_HOST_MGR_METHOD = re.compile(r"\bself\.host\.([a-zA-Z_]\w+)\.([a-zA-Z_]\w+)")
# self.host.ATTR  (single level, dot notation)
RE_HOST_ATTR       = re.compile(r"\bself\.host\.([a-zA-Z_]\w+)\b")
# getattr(self.host, "ATTR", ...) — catches getattr-based access bypassing dot notation
RE_GETATTR_HOST    = re.compile(r'\bgetattr\s*\(\s*self\.host\s*,\s*["\']([a-zA-Z_]\w+)["\']')
# self.ATTR  (for main_window.py scanning)
RE_SELF_ATTR       = re.compile(r"\bself\.([a-zA-Z_]\w+)\b")
# self.main_window.METHOD or self.parent_window.METHOD (dialogs)
RE_DIALOG_CALL     = re.compile(r"\bself\.(?:main_window|parent_window|mw)\.([a-zA-Z_]\w+)\b")
# hasattr
RE_HASATTR         = re.compile(r'\bhasattr\s*\(\s*self(?:\.host)?\s*,\s*["\']([a-zA-Z_]\w*)["\']')

# ---------------------------------------------------------------------------
# Step 3 — Scan a single file
# ---------------------------------------------------------------------------

def scan_manager_file(fpath: Path, per_mgr: dict, method_to_mgr: dict,
                      host_direct: set) -> list:
    """Scan a manager file for broken host.MANAGER.METHOD and misrouted host.METHOD."""
    if not fpath.exists():
        return []
    text   = fpath.read_text(encoding="utf-8", errors="replace")
    lines  = text.splitlines()
    findings = []

    for lineno, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue

        # -- Check A: self.host.MANAGER.METHOD ---------------------
        for m in RE_HOST_MGR_METHOD.finditer(raw):
            mgr_attr = m.group(1)
            method   = m.group(2)
            # skip if MANAGER is a known Qt/MainWindow direct attr (widget, action, etc.)
            if mgr_attr in QT_METHODS and mgr_attr not in ALL_MANAGER_NAMES:
                continue
            # skip if MANAGER is a dynamically-set host attr (plugin_manager, _template_dialog, etc.)
            if mgr_attr in host_direct and mgr_attr not in ALL_MANAGER_NAMES:
                continue
            if mgr_attr not in ALL_MANAGER_NAMES:
                findings.append({"kind": "INVALID_HOST_MANAGER",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"'{mgr_attr}' is not a known manager",
                                  "line": stripped})
                continue
            known = per_mgr.get(mgr_attr, set())
            if method not in known and not method.startswith("_"):
                findings.append({"kind": "MISSING_MANAGER_METHOD",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"'{method}' not found in {mgr_attr}",
                                  "line": stripped})

        # -- Check B: self.host.ATTR (single level, not followed by dot) --
        for m in RE_HOST_ATTR.finditer(raw):
            attr = m.group(1)
            # skip if another dot follows (that's Check A territory)
            if m.end() < len(raw) and raw[m.end()] == ".":
                continue
            if attr in QT_METHODS or attr in host_direct:
                continue
            # Is this attr something that belongs to a manager?
            if attr in method_to_mgr:
                mgr = method_to_mgr[attr]
                findings.append({"kind": "SHOULD_ROUTE_VIA_MANAGER",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"→ self.host.{mgr}.{attr}",
                                  "line": stripped})

        # -- Check B2: getattr(self.host, "ATTR") bypasses dot-notation regex --
        for m in RE_GETATTR_HOST.finditer(raw):
            attr = m.group(1)
            if attr in QT_METHODS or attr in host_direct:
                continue
            if attr in method_to_mgr:
                mgr = method_to_mgr[attr]
                findings.append({"kind": "SHOULD_ROUTE_VIA_MANAGER",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"→ getattr(self.host.{mgr}, '{attr}')",
                                  "line": stripped})

    return findings


def scan_main_window(fpath: Path, per_mgr: dict, method_to_mgr: dict) -> list:
    """Scan main_window.py for self.METHOD() that now lives on a manager."""
    if not fpath.exists():
        return []
    text   = fpath.read_text(encoding="utf-8", errors="replace")
    lines  = text.splitlines()

    # Collect names defined directly on MainWindow class
    mw_local = ast_collect_names(fpath)

    findings = []
    for lineno, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue
        for m in RE_SELF_ATTR.finditer(raw):
            attr = m.group(1)
            # skip if immediately followed by another dot (self.MANAGER.X handled elsewhere)
            if m.end() < len(raw) and raw[m.end()] == ".":
                continue
            if attr in QT_METHODS or attr in mw_local or attr.startswith("_"):
                continue
            if attr in method_to_mgr:
                mgr = method_to_mgr[attr]
                findings.append({"kind": "MAINWINDOW_SELF_METHOD",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"→ self.{mgr}.{attr}",
                                  "line": stripped})
    return findings


def scan_dialog_file(fpath: Path, method_to_mgr: dict) -> list:
    """
    Scan a dialog/widget file for self.main_window.METHOD where METHOD
    now lives on a manager (should be self.main_window.MANAGER.METHOD).
    """
    if not fpath.exists():
        return []
    text   = fpath.read_text(encoding="utf-8", errors="replace")
    lines  = text.splitlines()
    findings = []
    for lineno, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue
        for m in RE_DIALOG_CALL.finditer(raw):
            attr = m.group(1)
            # skip if another dot follows (self.main_window.MANAGER.X — already correct)
            if m.end() < len(raw) and raw[m.end()] == ".":
                continue
            if attr in QT_METHODS or attr.startswith("_"):
                continue
            if attr in method_to_mgr:
                mgr = method_to_mgr[attr]
                prefix = re.search(r"self\.(main_window|parent_window|mw)", raw)
                pname  = prefix.group(0) if prefix else "self.main_window"
                findings.append({"kind": "DIALOG_DIRECT_CALL",
                                  "lineno": lineno, "expr": m.group(0),
                                  "detail": f"→ {pname}.{mgr}.{attr}",
                                  "line": stripped})
    return findings

# ---------------------------------------------------------------------------
# Step 4 — Entry-point
# ---------------------------------------------------------------------------

KIND_LABEL = {
    "MISSING_MANAGER_METHOD":   "MISSING_MANAGER_METHOD",
    "INVALID_HOST_MANAGER":     "INVALID_HOST_MANAGER ",
    "SHOULD_ROUTE_VIA_MANAGER": "SHOULD_ROUTE_VIA_MGR ",
    "MAINWINDOW_SELF_METHOD":   "MAINWINDOW_SELF_METHOD",
    "DIALOG_DIRECT_CALL":       "DIALOG_DIRECT_CALL   ",
}

def print_findings(findings: list, rel: Path, counter: list):
    if not findings:
        return
    print(f"\n{'='*72}\n  {rel}\n{'='*72}")
    for f in findings:
        counter[0] += 1
        label = KIND_LABEL.get(f["kind"], f["kind"])
        print(f"  [{label}] line {f['lineno']}")
        print(f"    Expr:    {f['expr']}")
        print(f"    Fix:     {f['detail']}")
        print(f"    Context: {f['line'][:88]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true",
                    help="Also show hasattr/soft checks")
    ap.add_argument("--dialogs", action="store_true",
                    help="Also scan dialog files for self.main_window.METHOD patterns")
    args = ap.parse_args()

    # All files to pre-scan for attribute maps
    manager_files = [v for v in MANAGER_FILES.values() if v.exists()]
    all_scan = manager_files[:]
    if MAIN_WINDOW_FILE.exists():
        all_scan.append(MAIN_WINDOW_FILE)

    print("Building manager method maps from AST...")
    per_mgr, method_to_mgr, host_direct = build_method_maps(all_scan)

    counter = [0]

    # --- Scan manager files ---
    for fpath in manager_files:
        findings = scan_manager_file(fpath, per_mgr, method_to_mgr, host_direct)
        print_findings(findings, fpath.relative_to(ROOT), counter)

    # --- Scan main_window.py ---
    if MAIN_WINDOW_FILE.exists():
        findings = scan_main_window(MAIN_WINDOW_FILE, per_mgr, method_to_mgr)
        print_findings(findings, MAIN_WINDOW_FILE.relative_to(ROOT), counter)

    # --- Scan dialog files (optional) ---
    if args.dialogs:
        dialog_files = [
            f for f in SRC.glob("*.py")
            if f not in all_scan
            and f.name not in {"__init__.py", "sip_isdeleted_safe.py"}
            and "main_window" in f.read_text(encoding="utf-8", errors="replace")
        ]
        for fpath in sorted(dialog_files):
            findings = scan_dialog_file(fpath, method_to_mgr)
            print_findings(findings, fpath.relative_to(ROOT), counter)

    print(f"\n{'='*72}")
    print(f"Files scanned:        {len(all_scan)}" +
          (" + dialogs" if args.dialogs else ""))
    print(f"Disconnections found: {counter[0]}")


if __name__ == "__main__":
    main()
