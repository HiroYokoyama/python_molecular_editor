#!/usr/bin/env python3
"""
Quick codemod for tmp/plugins composition-architecture connection fixes.

Focus:
- replace removed `plugin_manager.get_context(...)`
- route former bare `mw.*` methods/attrs through managers
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path("tmp/plugins")
PLUGIN_CONTEXT_IMPORT = "from moleditpy.plugins.plugin_interface import PluginContext"

REPLACEMENTS = [
    # Removed API: plugin_manager.get_context(...)
    (
        r"\b(\w+)\s*=\s*(\w+)\.plugin_manager\.get_context\(([^)]+)\)",
        r"\1 = PluginContext(\2.plugin_manager, \3)",
    ),
    # Bare draw call -> composition manager
    (r"\bself\.mw\.draw_molecule_3d\(", "self.mw.view_3d_manager.draw_molecule_3d("),
    (r"\bself\.main_window\.draw_molecule_3d\(", "self.main_window.view_3d_manager.draw_molecule_3d("),
    (r"\bmw\.draw_molecule_3d\(", "mw.view_3d_manager.draw_molecule_3d("),
    (r"\bmain_window\.draw_molecule_3d\(", "main_window.view_3d_manager.draw_molecule_3d("),
    # Former direct state/edit methods
    (r"\bself\.mw\.push_undo_state\(", "self.mw.edit_actions_manager.push_undo_state("),
    (r"\bself\.main_window\.push_undo_state\(", "self.main_window.edit_actions_manager.push_undo_state("),
    (r"\bmw\.push_undo_state\(", "mw.edit_actions_manager.push_undo_state("),
    (r"\bmain_window\.push_undo_state\(", "main_window.edit_actions_manager.push_undo_state("),
    (r"\bself\.mw\.update_window_title\(", "self.mw.state_manager.update_window_title("),
    (r"\bself\.main_window\.update_window_title\(", "self.main_window.state_manager.update_window_title("),
    (r"\bmw\.update_window_title\(", "mw.state_manager.update_window_title("),
    (r"\bmain_window\.update_window_title\(", "main_window.state_manager.update_window_title("),
    (r"\bself\.mw\.update_realtime_info\(", "self.mw.state_manager.update_realtime_info("),
    (r"\bself\.main_window\.update_realtime_info\(", "self.main_window.state_manager.update_realtime_info("),
    (r"\bmw\.update_realtime_info\(", "mw.state_manager.update_realtime_info("),
    (r"\bmain_window\.update_realtime_info\(", "main_window.state_manager.update_realtime_info("),
    # current_file_path moved to init_manager
    (r"\bself\.mw\.current_file_path\b", "self.mw.init_manager.current_file_path"),
    (r"\bself\.main_window\.current_file_path\b", "self.main_window.init_manager.current_file_path"),
    (r"\bmw\.current_file_path\b", "mw.init_manager.current_file_path"),
    (r"\bmain_window\.current_file_path\b", "main_window.init_manager.current_file_path"),
    # Additional scanner-driven routes
    (r"\bself\.mw\.restore_ui_for_editing\(", "self.mw.ui_manager.restore_ui_for_editing("),
    (r"\bmw\.restore_ui_for_editing\(", "mw.ui_manager.restore_ui_for_editing("),
    (r"\bself\.mw\.clear_2d_measurement_labels\(", "self.mw.edit_3d_manager.clear_2d_measurement_labels("),
    (r"\bmw\.clear_2d_measurement_labels\(", "mw.edit_3d_manager.clear_2d_measurement_labels("),
    (r"\bself\.mw\.data\b", "self.mw.state_manager.data"),
    (r"\bmw\.data\b", "mw.state_manager.data"),
    (r"\bself\.mw\.save_project\(", "self.mw.io_manager.save_project("),
    (r"\bmw\.save_project\(", "mw.io_manager.save_project("),
    (r"\bself\.mw\.selected_atoms_3d\b", "self.mw.edit_3d_manager.selected_atoms_3d"),
    (r"\bself\.mw\.selected_atoms_for_measurement\b", "self.mw.edit_3d_manager.selected_atoms_for_measurement"),
    (r"\bself\.mw\.toggle_measurement_mode\(", "self.mw.edit_3d_manager.toggle_measurement_mode("),
    (r"\bself\.mw\.estimate_bonds_from_distances\(", "self.mw.io_manager.estimate_bonds_from_distances("),
    (r"\bmw\.estimate_bonds_from_distances\(", "mw.io_manager.estimate_bonds_from_distances("),
    (r"\bself\.mw\.minimize_2d_panel\(", "self.mw.ui_manager.minimize_2d_panel("),
    (r"\bmw\.minimize_2d_panel\(", "mw.ui_manager.minimize_2d_panel("),
    # bad import route found by scanner
    (
        r"from moleditpy\.modules import main_window_main_init",
        "from moleditpy.ui import main_window_init as main_window_main_init",
    ),
]


def ensure_plugin_context_import(text: str) -> str:
    if "PluginContext(" not in text or PLUGIN_CONTEXT_IMPORT in text:
        return text

    lines = text.splitlines()
    insert_at = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("#!") or s.startswith("# -*-"):
            insert_at = i + 1
            continue
        if s.startswith("import ") or s.startswith("from "):
            insert_at = i + 1
            continue
        if s == "":
            continue
        break

    lines.insert(insert_at, PLUGIN_CONTEXT_IMPORT)
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def process_file(path: Path) -> tuple[bool, int]:
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="cp932")

    original = text
    count = 0

    for pattern, repl in REPLACEMENTS:
        text, n = re.subn(pattern, repl, text)
        count += n

    text = ensure_plugin_context_import(text)
    # scanner-specific correction: initial_settings belongs on MainWindow, not init_manager
    text = text.replace(
        "mw.init_manager.initial_settings = mw.init_manager.settings.copy()",
        "mw.initial_settings = mw.init_manager.settings.copy()",
    )
    if text != original:
        path.write_text(text, encoding="utf-8")
        return True, count
    return False, 0


def ensure_compat_run_stubs() -> int:
    touched = 0
    targets = [
        Path("tmp/plugins/pyscf_calculator/__init__.py"),
        Path("tmp/plugins/Utility/reaction_sketcher/__init__.py"),
    ]
    stub = (
        "\n\ndef run(mw):\n"
        "    if hasattr(mw, 'host'):\n"
        "        mw = mw.host\n"
        "    return\n"
    )
    for path in targets:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        if "def run(" in text:
            continue
        path.write_text(text + stub, encoding="utf-8")
        touched += 1
        print(f"updated: {path} (run stub)")
    return touched


def main() -> int:
    files = 0
    edits = 0
    for py in ROOT.rglob("*.py"):
        changed, n = process_file(py)
        if changed:
            files += 1
            edits += n
            print(f"updated: {py} ({n} replacement(s))")
    files += ensure_compat_run_stubs()
    print(f"done: files={files}, replacements={edits}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
