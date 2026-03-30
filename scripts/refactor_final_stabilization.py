#!/usr/bin/env python3
"""
refactor_composition_architecture.py -- Systematic refactoring for manager-based composition.

This script updates:
1. Class renames in conftest.py
2. Attribute/method redirections in both tests AND managers.

Usage:
    python scripts/refactor_composition_architecture.py           # dry-run
    python scripts/refactor_composition_architecture.py --apply   # apply changes
"""

import re
import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS_GUI = ROOT / "tests" / "gui"
SRC_UI = ROOT / "moleditpy" / "src" / "moleditpy" / "ui"

# --- 1. Module/Class Rename Rules ---
MODULE_RENAMES = [
    ("moleditpy.ui.main_window_init.MainWindowMainInit", "moleditpy.ui.main_window_init.MainInitManager"),
    ("moleditpy.ui.ui_manager.MainWindowUiManager", "moleditpy.ui.ui_manager.UIManager"),
    ("moleditpy.ui.compute_engine.MainWindowCompute", "moleditpy.ui.compute_logic.ComputeManager"),
    ("moleditpy.ui.compute_engine", "moleditpy.ui.compute_logic"),
    ("MainWindowProjectIo", "IOManager"),
    ("MainWindowMolecularParsers", "IOManager"),
    ("MainWindowViewLoaders", "IOManager"),
    ("MainWindowAnalysis", "ComputeManager"),
    ("MainWindowTemplate", "TemplateManager"),
    ("MainWindowDialogs", "DialogManager"),
]

# --- 2. Attribute/Method Redirections ---
# Format: (regex_pattern, replacement, description)
REDIRECTS = [
    # 1. Host-based redirections (e.g., self.host.data -> self.host.state_manager.data)
    (r"\b(host|main_window|window|mw)\.data\b", r"\1.state_manager.data", "data -> state_manager.data"),
    (r"\b(host|main_window|window|mw)\.scene\b", r"\1.init_manager.scene", "scene -> init_manager.scene"),
    (r"\b(host|main_window|window|mw)\.view_2d\b", r"\1.init_manager.view_2d", "view_2d -> init_manager.view_2d"),
    (r"\b(host|main_window|window|mw)\.view_3d\b", r"\1.view_3d_manager.view_3d", "view_3d -> view_3d_manager.view_3d"),
    (r"\b(host|main_window|window|mw)\.plotter\b", r"\1.view_3d_manager.plotter", "plotter -> view_3d_manager.plotter"),
    (r"\b(host|main_window|window|mw)\.current_mol\b", r"\1.view_3d_manager.current_mol", "current_mol -> view_3d_manager.current_mol"),
    (r"\b(host|main_window|window|mw)\.atom_positions_3d\b", r"\1.view_3d_manager.atom_positions_3d", "atom_pos_3d -> view_3d_manager"),
    (r"\b(host|main_window|window|mw)\.constraints_3d\b", r"\1.edit_3d_manager.constraints_3d", "constraints_3d -> edit_3d_manager"),
    (r"\b(host|main_window|window|mw)\.constraints_2d\b", r"\1.compute_manager.constraints_2d", "constraints_2d -> compute_manager"),
    (r"\b(host|main_window|window|mw)\.settings\b", r"\1.init_manager.settings", "settings -> init_manager.settings"),
    (r"\b(host|main_window|window|mw)\.settings_file\b", r"\1.init_manager.settings_file", "settings_file -> init_manager"),
    (r"\b(host|main_window|window|mw)\.settings_dirty\b", r"\1.init_manager.settings_dirty", "settings_dirty -> init_manager"),
    (r"\b(host|main_window|window|mw)\.optimization_method\b", r"\1.init_manager.optimization_method", "opt_method -> init_manager"),
    (r"\b(host|main_window|window|mw)\.last_successful_optimization_method\b", r"\1.compute_manager.last_successful_optimization_method", "last_opt -> compute_manager"),
    (r"\b(host|main_window|window|mw)\.optimization_results\b", r"\1.compute_manager.optimization_results", "opt_res -> compute_manager"),
    (r"\b(host|main_window|window|mw)\._calculating_text_actor\b", r"\1.compute_manager._calculating_text_actor", "_calc_text -> compute_manager"),
    (r"\b(host|main_window|window|mw)\.halt_ids\b", r"\1.compute_manager.halt_ids", "halt_ids -> compute_manager"),
    (r"\b(host|main_window|window|mw)\.active_worker_ids\b", r"\1.compute_manager.active_worker_ids", "active_workers -> compute_manager"),
    (r"\b(host|main_window|window|mw)\.next_conversion_id\b", r"\1.compute_manager.next_conversion_id", "next_conv_id -> compute_manager"),
    (r"\b(host|main_window|window|mw)\._active_calc_threads\b", r"\1.compute_manager._active_calc_threads", "threads -> compute_manager"),
    (r"\b(host|main_window|window|mw|mock_parser_host|dummy_host)\.current_file_path\b", r"\1.init_manager.current_file_path", "current_file_path -> init_manager"),
    (r"\b(host|main_window|window|mw|mock_parser_host|dummy_host)\.has_unsaved_changes\b", r"\1.state_manager.has_unsaved_changes", "has_unsaved_changes -> state_manager"),
    (r"\b(host|main_window|window|mw|mock_parser_host|dummy_host)\._saved_state\b", r"\1.state_manager._saved_state", "_saved_state -> state_manager"),
    (r"\b(io|self)\.current_file_path\b", r"\1.host.init_manager.current_file_path", "io.current_file_path -> io.host.init_manager"),
    (r"\b(io|self)\.has_unsaved_changes\b", r"\1.host.state_manager.has_unsaved_changes", "io.has_unsaved_changes -> io.host.state_manager"),

    # 2. Method-based redirections
    # --- Redirection Patterns (Regex) ---
    # (Pattern, Replacement, Description)

    # 1. Host Managers
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.data\b", r"\1.state_manager.data", "data -> state_manager.data"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.undo\b", r"\1.state_manager.undo", "undo -> state_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.scene\b", r"\1.init_manager.scene", "scene -> init_manager.scene"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.plotter\b", r"\1.view_3d_manager.plotter", "plotter -> view_3d_manager.plotter"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.view_2d\b", r"\1.init_manager.view_2d", "view_2d -> init_manager.view_2d"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.settings\b", r"\1.init_manager.settings", "settings -> init_manager.settings"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.current_mol\b", r"\1.view_3d_manager.current_mol", "mol -> view_3d_manager.current_mol"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.is_2d_editable\b", r"\1.ui_manager.is_2d_editable", "editable -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.constraints_3d\b", r"\1.edit_3d_manager.constraints_3d", "constraints -> edit_3d_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.atom_pos_3d\b", r"\1.view_3d_manager.atom_pos_3d", "atom_pos -> view_3d_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.active_worker_ids\b", r"\1.compute_manager.active_worker_ids", "workers -> compute_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.halt_ids\b", r"\1.compute_manager.halt_ids", "halt_ids -> compute_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.settings_dirty\b", r"\1.init_manager.settings_dirty", "dirty -> init_manager"),
    (r"\b((?:self\.)?(?:host|mw|window|main_window|mock_host|dummy_host|mock_parser_host|mock_window|MockWindow))(?:\.state_manager)?\.(undo_stack|redo_stack|push_undo_state|undo|redo)\b", r"\1.edit_actions_manager.\2", "undo/redo -> edit_actions"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.save_to_file\b", r"\1.io_manager.save_to_file", "save -> io_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.load_from_file\b", r"\1.io_manager.load_from_file", "load -> io_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.set_mode\b", r"\1.ui_manager.set_mode", "set_mode -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.refresh_2d_view\b", r"\1.ui_manager.refresh_2d_view", "refresh_2d -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_select_mode\b", r"\1.ui_manager.activate_select_mode", "select_mode -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_draw_mode\b", r"\1.ui_manager.activate_draw_mode", "draw_mode -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_erase_mode\b", r"\1.ui_manager.activate_erase_mode", "erase_mode -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_charge_plus_mode\b", r"\1.ui_manager.activate_charge_plus_mode", "charge_p -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_charge_minus_mode\b", r"\1.ui_manager.activate_charge_minus_mode", "charge_m -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.activate_radical_mode\b", r"\1.ui_manager.activate_radical_mode", "radical -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.draw_molecule_3d\b", r"\1.view_3d_manager.draw_molecule_3d", "draw_3d -> view_3d_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.update_chiral_labels\b", r"\1.view_3d_manager.update_chiral_labels", "chiral -> view_3d_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\._enable_3d_edit_actions\b", r"\1.ui_manager._enable_3d_edit_actions", "enable_3d_actions -> ui_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\._enable_3d_features\b", r"\1.ui_manager._enable_3d_features", "enable_3d_features -> ui_manager"),

    # UI Components
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.convert_button\b", r"\1.init_manager.convert_button", "convert_btn -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.optimize_3d_button\b", r"\1.init_manager.optimize_3d_button", "opt3d_btn -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.cleanup_button\b", r"\1.init_manager.cleanup_button", "cleanup_btn -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.export_button\b", r"\1.init_manager.export_button", "export_btn -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.toolbar\b", r"\1.init_manager.toolbar", "toolbar -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.toolbar_bottom\b", r"\1.init_manager.toolbar_bottom", "toolbar_bottom -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.measurement_action\b", r"\1.init_manager.measurement_action", "measurement -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.edit_3d_action\b", r"\1.init_manager.edit_3d_action", "edit3d -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.splitter\b", r"\1.init_manager.splitter", "splitter -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.plugin_toolbar\b", r"\1.init_manager.plugin_toolbar", "plugin_tb -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.formula_label\b", r"\1.init_manager.formula_label", "formula -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.left_panel\b", r"\1.init_manager.left_panel", "left_pan -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.right_panel\b", r"\1.init_manager.right_panel", "right_pan -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.undo_action\b", r"\1.init_manager.undo_action", "undo_act -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.redo_action\b", r"\1.init_manager.redo_action", "redo_act -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.copy_action\b", r"\1.init_manager.copy_action", "copy_act -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.cut_action\b", r"\1.init_manager.cut_action", "cut_act -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.paste_action\b", r"\1.init_manager.paste_action", "paste_act -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.analysis_action\b", r"\1.init_manager.analysis_action", "analysis -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.style_button\b", r"\1.init_manager.style_button", "style_btn -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.tool_group\b", r"\1.init_manager.tool_group", "tool_group -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.mode_actions\b", r"\1.init_manager.mode_actions", "mode_actions -> init_manager"),
    (r"\b((?:self\.)?(?:host|main_window|window|mw))\.other_atom_action\b", r"\1.init_manager.other_atom_action", "other_atom -> init_manager"),

    # 3. String-based Lookups (hasattr, getattr, setattr)
    # This addresses the "logic flips" where hasattr(self.host, 'plotter') was failing.
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]data['\"]\)", r"hasattr(\1.state_manager, 'data')", "hasattr(data)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]scene['\"]\)", r"hasattr(\1.init_manager, 'scene')", "hasattr(scene)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]plotter['\"]\)", r"hasattr(\1.view_3d_manager, 'plotter')", "hasattr(plotter)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]view_2d['\"]\)", r"hasattr(\1.init_manager, 'view_2d')", "hasattr(view_2d)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]optimize_3d_button['\"]\)", r"hasattr(\1.init_manager, 'optimize_3d_button')", "hasattr(opt3d_btn)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]cleanup_button['\"]\)", r"hasattr(\1.init_manager, 'cleanup_button')", "hasattr(cleanup_btn)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]export_button['\"]\)", r"hasattr(\1.init_manager, 'export_button')", "hasattr(export_btn)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]convert_button['\"]\)", r"hasattr(\1.init_manager, 'convert_button')", "hasattr(convert_btn)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]undo_action['\"]\)", r"hasattr(\1.init_manager, 'undo_action')", "hasattr(undo)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]redo_action['\"]\)", r"hasattr(\1.init_manager, 'redo_action')", "hasattr(redo)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]copy_action['\"]\)", r"hasattr(\1.init_manager, 'copy_action')", "hasattr(copy)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]cut_action['\"]\)", r"hasattr(\1.init_manager, 'cut_action')", "hasattr(cut)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]paste_action['\"]\)", r"hasattr(\1.init_manager, 'paste_action')", "hasattr(paste)"),
    (r"\bhasattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]splitter['\"]\)", r"hasattr(\1.init_manager, 'splitter')", "hasattr(splitter)"),
    (r"\bgetattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]data['\"]", r"getattr(\1.state_manager, 'data'", "getattr(data)"),
    (r"\bgetattr\(((?:self\.)?(?:host|main_window|window|mw)),\s*['\"]plotter['\"]", r"getattr(\1.view_3d_manager, 'plotter'", "getattr(plotter)"),
    (r"\bgetattr\((host|main_window|window|mw),\s*['\"]_calculating_text_actor['\"]", r"getattr(\1.compute_manager, '_calculating_text_actor'", "getattr(calc_text)"),
]

def refactor_file(fpath: Path, apply: bool):
    try:
        text = fpath.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = fpath.read_text(encoding="shift-jis", errors="replace")
    
    original_text = text
    changes = []

    # 1. Apply Module/Class Renames
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
    parser = argparse.ArgumentParser(description="Refactor for manager composition.")
    parser.add_argument("--apply", action="store_true", help="Apply changes to files")
    args = parser.parse_args()

    files_to_check = []
    files_to_check.extend(sorted((ROOT / "tests" / "unit").glob("*.py")))
    files_to_check.extend(sorted((ROOT / "tests" / "integration").glob("*.py")))
    files_to_check.extend(sorted(TESTS_GUI.glob("*.py")))
    files_to_check.extend(sorted(SRC_UI.glob("*.py")))

    total_files = 0
    total_changes = 0

    for fpath in files_to_check:
        if fpath.name in ["main_window.py", "sitecustomize.py"]: continue
        
        changes = refactor_file(fpath, args.apply)
        if changes:
            total_files += 1
            total_changes += len(changes)
            rel = fpath.relative_to(ROOT)
            status = "[FIXED]" if args.apply else "[DRY-RUN]"
            print(f"\n{rel}:")
            # Only print first 5 changes per file to avoid spam
            for c in changes[:5]:
                print(f"  {status} {c}")
            if len(changes) > 5:
                print(f"  ... and {len(changes)-5} more changes.")

    print(f"\nRefactored {total_files} files with {total_changes} changes.")
    if not args.apply and total_changes > 0:
        print("Run with --apply to save changes.")

if __name__ == "__main__":
    main()
