#!/usr/bin/env python3
"""
refactor_gui_tests.py -- Systematic refactoring of GUI tests for manager-based composition.

This script updates:
1. Class renames in conftest.py (e.g., MainWindowMainInit -> MainInitManager)
2. Attribute/method redirections in test_*.py (e.g., window.set_mode -> window.ui_manager.set_mode)

Usage:
    python scripts/refactor_gui_tests.py           # dry-run
    python scripts/refactor_gui_tests.py --apply   # apply changes
"""

import re
import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS_GUI = ROOT / "tests" / "gui"

# --- 1. Module/Class Rename Rules (Mostly for conftest.py) ---
# Format: (old_string, new_string)
# NOTE: Order from most specific to least specific to avoid double-replacement.
MODULE_RENAMES = [
    ("moleditpy.ui.main_window_init.MainWindowMainInit", "moleditpy.ui.main_window_init.MainInitManager"),
    ("moleditpy.ui.ui_manager.MainWindowUiManager", "moleditpy.ui.ui_manager.UIManager"),
    ("moleditpy.ui.compute_engine.MainWindowCompute", "moleditpy.ui.compute_logic.ComputeManager"),
    ("moleditpy.ui.compute_engine", "moleditpy.ui.compute_logic"), # Module name change
    ("MainWindowProjectIo", "IOManager"),
    ("MainWindowMolecularParsers", "IOManager"),
    ("MainWindowViewLoaders", "IOManager"),
    ("MainWindowAnalysis", "ComputeManager"),
    ("MainWindowTemplate", "TemplateManager"),
    ("MainWindowDialogs", "DialogManager"),
    ("_mwcomp.MainWindowCompute", "_mwcomp.ComputeManager"),
    ("_mwcomp2.MainWindowCompute", "_mwcomp2.ComputeManager"),
]

# --- 2. Attribute/Method Redirections (Mostly for test_*.py) ---
# Format: (regex_pattern, replacement, description)
REDIRECTS = [
    # Edit Actions Manager
    (r"\bwindow\.clear_2d_editor\b", "window.edit_actions_manager.clear_2d_editor", "clear_2d_editor -> edit_actions_manager"),
    (r"\bwindow\.select_all\b", "window.edit_actions_manager.select_all", "select_all -> edit_actions_manager"),
    (r"\bwindow\.copy_selection\b", "window.edit_actions_manager.copy_selection", "copy_selection -> edit_actions_manager"),
    (r"\bwindow\.paste_from_clipboard\b", "window.edit_actions_manager.paste_from_clipboard", "paste_from_clipboard -> edit_actions_manager"),
    (r"\bwindow\.clear_all\b", "window.edit_actions_manager.clear_all", "clear_all -> edit_actions_manager"),
    
    # String Importer Manager
    (r"\bwindow\.import_smiles_dialog\b", "window.string_importer_manager.import_smiles_dialog", "import_smiles -> string_importer_manager"),
    (r"\bwindow\.import_inchi_dialog\b", "window.string_importer_manager.import_inchi_dialog", "import_inchi -> string_importer_manager"),
    (r"\bwindow\.import_inchikey_dialog\b", "window.string_importer_manager.import_inchikey_dialog", "import_inchikey -> string_importer_manager"),
    
    # View 3D Manager
    (r"\bwindow\.current_3d_style\b", "window.view_3d_manager.current_3d_style", "current_3d_style -> view_3d_manager"),
    (r"\bwindow\.view_3d_manager\.draw_molecule_3d\b", "window.view_3d_manager.draw_molecule_3d", "draw_3d -> view_3d_manager"),
    
    # IO Manager
    (r"\bwindow\.load_mol_file_for_3d_viewing\b", "window.io_manager.load_mol_file_for_3d_viewing", "load_mol_3d -> io_manager"),
    (r"\bwindow\.load_xyz_for_3d_viewing\b", "window.io_manager.load_xyz_for_3d_viewing", "load_xyz_3d -> io_manager"),
    (r"\bwindow\.open_project_file\b", "window.io_manager.open_project_file", "open_project -> io_manager"),
    
    # Init Manager (UI Elements & Init Logic)
    (r"\bwindow\.undo_action\b", "window.init_manager.undo_action", "undo_action -> init_manager"),
    (r"\bwindow\.redo_action\b", "window.init_manager.redo_action", "redo_action -> init_manager"),
    (r"\bwindow\.toolbar\b", "window.init_manager.toolbar", "toolbar -> init_manager"),
    (r"\bwindow\.toolbar_bottom\b", "window.init_manager.toolbar_bottom", "toolbar_bottom -> init_manager"),
    (r"\bwindow\.style_button\b", "window.init_manager.style_button", "style_button -> init_manager"),
    (r"\bwindow\.convert_button\b", "window.init_manager.convert_button", "convert_button -> init_manager"),
    (r"\bwindow\.optimize_3d_button\b", "window.init_manager.optimize_3d_button", "optimize_3d_button -> init_manager"),
    (r"\bwindow\.export_button\b", "window.init_manager.export_button", "export_button -> init_manager"),
    (r"\bwindow\.analysis_action\b", "window.init_manager.analysis_action", "analysis_action -> init_manager"),
    (r"\bwindow\.load_command_line_file\b", "window.init_manager.load_command_line_file", "load_cmd -> init_manager"),
    (r"\bwindow\.apply_initial_settings\b", "window.init_manager.apply_initial_settings", "apply_settings -> init_manager"),
    (r"\bwindow\.reset_all_settings_menu\b", "window.init_manager.reset_all_settings_menu", "reset_settings_menu -> init_manager"),
    (r"\bwindow\._perform_settings_reset\b", "window.init_manager._perform_settings_reset", "perform_reset -> init_manager"),
    (r"\bwindow\._refresh_ui_after_reset\b", "window.init_manager._refresh_ui_after_reset", "refresh_ui -> init_manager"),
    (r"\bwindow\.load_settings\b", "window.init_manager.load_settings", "load_settings -> init_manager"),
    
    # Patch targets in conftest.py
    (r"main_window\.toolbar =", "main_window.init_manager.toolbar =", "toolbar_init -> init_manager"),
    (r"main_window\.toolbar_bottom =", "main_window.init_manager.toolbar_bottom =", "toolbar_bottom_init -> init_manager"),
]

def refactor_file(fpath: Path, apply: bool):
    try:
        text = fpath.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = fpath.read_text(encoding="shift-jis", errors="replace")
    
    original_text = text
    changes = []

    # 1. Apply Module/Class Renames
    # Sort by length descending to avoid partial matches
    sorted_renames = sorted(MODULE_RENAMES, key=lambda x: len(x[0]), reverse=True)
    for old_name, new_name in sorted_renames:
        if old_name in text:
            pattern = re.compile(rf"\b{re.escape(old_name)}\b")
            new_text = pattern.sub(new_name, text)
            if new_text != text:
                text = new_text
                changes.append(f"Rename: {old_name} -> {new_name}")

    # 2. Apply Attribute Redirections
    for pattern_str, replacement, desc in REDIRECTS:
        pattern = re.compile(pattern_str)
        if pattern.search(text):
            text = pattern.sub(replacement, text)
            changes.append(f"Redirect: {desc}")

    if text != original_text and apply:
        fpath.write_text(text, encoding="utf-8")

    return changes

def main():
    parser = argparse.ArgumentParser(description="Refactor GUI tests for manager composition.")
    parser.add_argument("--apply", action="store_true", help="Apply changes to files")
    args = parser.parse_args()

    target_files = sorted(TESTS_GUI.glob("*.py"))
    total_files = 0
    total_changes = 0

    for fpath in target_files:
        if fpath.name == "sitecustomize.py": continue
        
        changes = refactor_file(fpath, args.apply)
        if changes:
            total_files += 1
            total_changes += len(changes)
            rel = fpath.relative_to(ROOT)
            status = "[FIXED]" if args.apply else "[DRY-RUN]"
            print(f"\n{rel}:")
            for c in changes:
                print(f"  {status} {c}")

    print(f"\nRefactored {total_files} files with {total_changes} changes.")
    if not args.apply and total_changes > 0:
        print("Run with --apply to save changes.")

if __name__ == "__main__":
    main()
