#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
refactor_mixins_to_managers.py
==============================
Refactors the remaining MainWindow mixins into composition-based Managers.

This script performs the transformation in several steps:
1.  Collects all public methods from the mixins to build a re-qualification map.
2.  Updates MainWindow (main_window.py) to remove inheritance and add manager instances.
3.  Transforms Mixin classes into Manager classes (renaming, adding __init__(self, host)).
4.  Globally updates all call-sites from self.METHOD to self.MANAGER.METHOD.
"""

import argparse
import os
import re
import sys
import difflib
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "moleditpy" / "src" / "moleditpy"
UI_DIR = SRC_ROOT / "ui"

# Mixin File -> (Mixin Class Name, New Manager Class Name, Manager Attr Name on MW)
MIXIN_CONFIG = {
    UI_DIR / "app_state.py":        ("MainWindowAppState",        "StateManager",          "state_manager"),
    UI_DIR / "main_window_init.py": ("MainWindowMainInit",        "MainInitManager",       "init_manager"),
    UI_DIR / "string_importers.py": ("MainWindowStringImporters", "StringImporterManager", "string_importer_manager"),
    UI_DIR / "ui_manager.py":       ("MainWindowUiManager",       "UIManager",             "ui_manager"),
}

# Known MainWindow attributes that need 'self.host.' prefix in managers
MW_ATTRS = [
    "data", "scene", "view_2d", "current_mol", "plotter", "statusBar", "settings",
    "current_file_path", "has_unsaved_changes", "plugin_manager", "undo_stack",
    "redo_stack", "constraints_3d", "splitter", "mode_actions", "view_3d_manager",
    "edit_3d_manager", "edit_actions_manager", "compute_manager", "dialog_manager",
    "io_manager", "settings_file", "settings_dir", "initial_settings", "settings_dirty",
    "initialization_complete", "_ih_update_counter", "_active_calc_threads",
    "formula_label", "cleanup_button", "convert_button", "tool_group", "is_2d_editable",
    "minimize_2d_panel", "restore_2d_panel", "update_window_title",
    "_enter_3d_viewer_ui_mode", "restore_ui_for_editing", "_enable_3d_features",
    "_enable_3d_edit_actions", "show_chiral_labels", "is_xyz_derived",
    "chem_check_tried", "chem_check_failed", "measurement_action", "edit_3d_action",
    "style_button", "optimize_3d_button", "export_button", "undo_action", "redo_action",
    "cut_action", "copy_action", "paste_action", "active_3d_dialogs", "_is_restoring_state"
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")

def _write(path: Path, text: str, dry_run: bool) -> None:
    if not dry_run:
        path.write_text(text, encoding="utf-8")

def _show_diff(path: Path, original: str, updated: str) -> None:
    if original == updated:
        return
    rel = path.relative_to(REPO_ROOT)
    unified = list(difflib.unified_diff(
        original.splitlines(keepends=True),
        updated.splitlines(keepends=True),
        fromfile=str(rel), tofile=str(rel),
        n=1,
    ))
    added = sum(1 for l in unified if l.startswith("+") and not l.startswith("+++"))
    removed = sum(1 for l in unified if l.startswith("-") and not l.startswith("---"))
    print(f"\n  [{rel}]  +{added}/-{removed} lines")
    for line in unified[:15]:
        # Handle encoding issues for console
        safe = line.rstrip("\n").encode("ascii", errors="replace").decode("ascii")
        print(f"    {safe}")
    if len(unified) > 15:
        print(f"    ... ({len(unified) - 15} more diff lines)")

# ---------------------------------------------------------------------------
# Core Logic
# ---------------------------------------------------------------------------

def collect_methods() -> Dict[str, str]:
    """Build a mapping of public method name -> manager_attr_on_mw."""
    method_to_mgr = {}
    for fpath, (old_cls, new_cls, mgr_attr) in MIXIN_CONFIG.items():
        if not fpath.exists(): continue
        src = _read(fpath)
        # Find def method(self, ...) inside the class body
        # Simple indentation-based heuristic
        matches = re.finditer(r"^\s{4}def ([a-zA-Z0-9_]+)\s*\(self", src, re.MULTILINE)
        for m in matches:
            name = m.group(1)
            if name == "__init__": continue
            if name in method_to_mgr:
                print(f"  [WARN] Duplicate method '{name}' found in {mgr_attr} and {method_to_mgr[name]}")
            else:
                method_to_mgr[name] = mgr_attr
    return method_to_mgr

def fix_main_window(dry_run: bool):
    print("\n=== Phase: Refactoring MainWindow.py ===")
    fpath = UI_DIR / "main_window.py"
    if not fpath.exists(): return
    original = _read(fpath)
    text = original
    
    # 1. Update inheritance and imports
    for old_cls, new_cls, mgr_attr in MIXIN_CONFIG.values():
        text = text.replace(f"    {old_cls},", "") # Remove from class def
        text = text.replace(f"from .app_state import {old_cls}", f"from .app_state import {new_cls}")
        text = text.replace(f"from .main_window_init import {old_cls}", f"from .main_window_init import {new_cls}")
        text = text.replace(f"from .string_importers import {old_cls}", f"from .string_importers import {new_cls}")
        text = text.replace(f"from .ui_manager import {old_cls}", f"from .ui_manager import {new_cls}")

    # Remove extra commas/newlines in class def if any
    text = re.sub(r"class MainWindow\(\s*,\s*QMainWindow\s*\):", r"class MainWindow(QMainWindow):", text)

    # 2. Add manager initializations
    init_pos = text.find("self.io_manager = IOManager(self)")
    if init_pos != -1:
        insertion_point = text.find("\n", init_pos) + 1
        mgr_inits = [f"        self.{mgr_attr} = {new_cls}(self)" for _, (_, new_cls, mgr_attr) in MIXIN_CONFIG.items()]
        text = text[:insertion_point] + "\n".join(mgr_inits) + "\n" + text[insertion_point:]

    # 3. Update __init__ re-calls
    text = text.replace("MainWindowMainInit.__init__(self,", "self.init_manager.__init__(self,")
    text = text.replace("MainWindowAppState.__init__(self)", "self.state_manager.__init__(self)")

    # 4. Add bridge methods for eventFilter/closeEvent if they were in UiManager
    if "def closeEvent(self, event):" not in text:
        text = text.rstrip() + "\n\n    def closeEvent(self, event):\n        self.ui_manager.closeEvent(event)\n\n    def eventFilter(self, obj, event):\n        return self.ui_manager.eventFilter(obj, event)\n"

    _show_diff(fpath, original, text)
    _write(fpath, text, dry_run)

def fix_mixin_to_manager(dry_run: bool):
    print("\n=== Phase: Converting Mixins to Managers ===")
    for fpath, (old_cls, new_cls, mgr_attr) in MIXIN_CONFIG.items():
        if not fpath.exists(): continue
        original = _read(fpath)
        text = original
        
        # 1. Rename class
        text = re.sub(rf"\bclass {old_cls}\b", f"class {new_cls}", text)
        
        # 2. Add self.host = host to __init__
        init_match = re.search(r"def __init__\(self(.*?)\):", text)
        if init_match:
            args = init_match.group(1)
            if "host" not in args:
                text = text.replace(init_match.group(0), f"def __init__(self, host{args}):")
                insertion_point = text.find(f"def __init__(self, host{args}):") + len(f"def __init__(self, host{args}):")
                next_line = text.find("\n", insertion_point) + 1
                text = text[:next_line] + "        self.host = host\n" + text[next_line:]
        else:
            # Inline insertion after docstring or class def
            insert_pos = text.find(f"class {new_cls}")
            insert_pos = text.find("\n", insert_pos) + 1
            if text[insert_pos:].strip().startswith('"""'):
                doc_end = text.find('"""', insert_pos + 3) + 3
                insert_pos = text.find("\n", doc_end) + 1
            text = text[:insert_pos] + "    def __init__(self, host):\n        self.host = host\n\n" + text[insert_pos:]

    # 3. Hard-re-link known MW attributes (self.data -> self.host.data)
        for attr in MW_ATTRS:
            # self.attr -> self.host.attr (not self.manager.attr)
            text = re.sub(rf"\bself\.{re.escape(attr)}\b", f"self.host.{attr}", text)
        
        # Fix double nesting self.host.host
        text = text.replace("self.host.host", "self.host")

        _show_diff(fpath, original, text)
        _write(fpath, text, dry_run)

def refactor_call_sites(dry_run: bool, method_map: Dict[str, str]):
    print("\n=== Phase: Re-qualifying Method Calls ===")
    all_py = sorted(list(UI_DIR.glob("*.py")))
    
    for fpath in all_py:
        original = _read(fpath)
        text = original
        
        # Skip certain files where re-qualification might be dangerous or handled differently
        is_mainwindow = fpath.name == "main_window.py"
        is_manager = any(fpath == p for p in MIXIN_CONFIG.keys()) or "_logic.py" in fpath.name
        
        for method, mgr in method_map.items():
            # Pattern A: self.METHOD -> self.MGR.METHOD (in MainWindow or Mixins)
            # Pattern B: self.host.METHOD -> self.host.MGR.METHOD (in other Managers)
            # Pattern C: VAR.METHOD -> VAR.MGR.METHOD (in Dialogs) - too risky for generic script?
            
            if is_mainwindow:
                # self.init_ui -> self.init_manager.init_ui
                pattern = rf"\bself\.(?!{re.escape(mgr)}\.)({re.escape(method)})\b"
                text = re.sub(pattern, rf"self.{mgr}.\1", text)
            elif is_manager:
                # self.host.init_ui -> self.host.init_manager.init_ui
                pattern = rf"\bself\.host\.(?!{re.escape(mgr)}\.)({re.escape(method)})\b"
                text = re.sub(pattern, rf"self.host.{mgr}.\1", text)
                
                # Special: in the mixin-being-converted, self.method stays self.method
                # (because it's now in the manager class itself).
                # But if it's calling a method on ANOTHER manager, it must be re-qualified via host.
                if fpath in MIXIN_CONFIG and MIXIN_CONFIG[fpath][2] != mgr:
                    pattern = rf"\bself\.(?!{re.escape(mgr)}\.)({re.escape(method)})\b"
                    text = re.sub(pattern, rf"self.host.{mgr}.\1", text)
            else:
                # In Dialogs: self.main_window.init_ui -> self.main_window.init_manager.init_ui
                pattern = rf"\bself\.main_window\.(?!{re.escape(mgr)}\.)({re.escape(method)})\b"
                text = re.sub(pattern, rf"self.main_window.{mgr}.\1", text)

        if text != original:
            _show_diff(fpath, original, text)
            _write(fpath, text, dry_run)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    dry_run = not args.apply
    
    if dry_run:
        print("DRY RUN - no changes will be applied.")
    
    method_map = collect_methods()
    
    # 1. Transform structure
    fix_mixin_to_manager(dry_run)
    fix_main_window(dry_run)
    
    # 2. Transform calls
    refactor_call_sites(dry_run, method_map)

if __name__ == "__main__":
    main()
