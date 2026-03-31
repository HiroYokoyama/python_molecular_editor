"""
fix_plugin_api.py
Fixes all tmp/plugins files to use the correct manager sub-objects
instead of calling deprecated direct attributes on MainWindow.

Run with:  python fix_plugin_api.py
"""
import re
import os
import ast
import glob

PLUGIN_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp", "plugins")

# ---------------------------------------------------------------------------
# 1. Token replacements applied to EVERY .py file
# ---------------------------------------------------------------------------
GLOBAL_REPLACEMENTS = [
    # has_unsaved_changes
    (r'\bself\.mw\.has_unsaved_changes\b',  'self.mw.state_manager.has_unsaved_changes'),
    (r'(?<!\w)mw\.has_unsaved_changes\b',   'mw.state_manager.has_unsaved_changes'),

    # update_realtime_info
    (r'\bself\.mw\.update_realtime_info\(\)',  'self.mw.state_manager.update_realtime_info()'),
    (r'(?<!\w)mw\.update_realtime_info\(\)',   'mw.state_manager.update_realtime_info()'),

    # update_window_title
    (r'\bself\.mw\.update_window_title\(\)',   'self.mw.state_manager.update_window_title()'),
    (r'(?<!\w)mw\.update_window_title\(\)',    'mw.state_manager.update_window_title()'),

    # update_undo_redo_actions
    (r'\bself\.mw\.update_undo_redo_actions\(\)',  'self.mw.edit_actions_manager.update_undo_redo_actions()'),
    (r'(?<!\w)mw\.update_undo_redo_actions\(\)',   'mw.edit_actions_manager.update_undo_redo_actions()'),

    # push_undo_state
    (r'\bself\.mw\.push_undo_state\(\)',  'self.mw.edit_actions_manager.push_undo_state()'),
    (r'(?<!\w)mw\.push_undo_state\(\)',   'mw.edit_actions_manager.push_undo_state()'),

    # check_chemistry_problems_fallback
    (r'\bself\.mw\.check_chemistry_problems_fallback\(\)',  'self.mw.compute_manager.check_chemistry_problems_fallback()'),
    (r'(?<!\w)mw\.check_chemistry_problems_fallback\(\)',   'mw.compute_manager.check_chemistry_problems_fallback()'),

    # trigger_conversion (method reference AND call)
    (r'\bself\.mw\.trigger_conversion\b',  'self.mw.compute_manager.trigger_conversion'),
    (r'(?<!\w)mw\.trigger_conversion\b',   'mw.compute_manager.trigger_conversion'),

    # on_calculation_finished
    (r'\bself\.mw\.on_calculation_finished\b',  'self.mw.compute_manager.on_calculation_finished'),
    (r'(?<!\w)mw\.on_calculation_finished\b',   'mw.compute_manager.on_calculation_finished'),

    # convert_button
    (r'\bself\.mw\.convert_button\b',  'self.mw.init_manager.convert_button'),
    (r'(?<!\w)mw\.convert_button\b',   'mw.init_manager.convert_button'),

    # splitter (reaction_sketcher and any other plugin)
    (r'\bself\.main_window\.splitter\b',  'self.main_window.init_manager.splitter'),
    (r'(?<!\w)mw\.splitter\b',            'mw.init_manager.splitter'),

    # view_2d
    (r'\bself\.main_window\.view_2d\b',   'self.main_window.init_manager.view_2d'),
    (r'(?<!\w)mw\.view_2d\b',             'mw.init_manager.view_2d'),

    # push_undo_state (self.main_window. prefix variant)
    (r'\bself\.main_window\.push_undo_state\(\)',  'self.main_window.edit_actions_manager.push_undo_state()'),
    # hasattr guard variant — keep the attribute name consistent
    (r"hasattr\(self\.main_window,\s*'push_undo_state'\)",
     "hasattr(self.main_window.edit_actions_manager, 'push_undo_state')"),

    # moleditpy.modules.* → new module paths
    (r'from moleditpy\.modules\.atom_item import',         'from moleditpy.ui.atom_item import'),
    (r'from moleditpy\.modules\.bond_item import',         'from moleditpy.ui.bond_item import'),
    (r'from moleditpy\.modules import atom_item',          'from moleditpy.ui import atom_item'),
    (r'from moleditpy\.modules import bond_item',          'from moleditpy.ui import bond_item'),
    (r'from moleditpy\.modules\.molecular_data import',    'from moleditpy.core.molecular_data import'),
    (r'from moleditpy\.modules\.main_window_compute import MainWindowCompute',
     'from moleditpy.ui.compute_logic import ComputeManager as MainWindowCompute'),
    (r'from moleditpy\.modules\.main_window_molecular_parsers import MainWindowMolecularParsers',
     'from moleditpy.ui.io_logic import IOManager as MainWindowMolecularParsers'),
    (r'from moleditpy\.modules import main_window_main_init',
     'from moleditpy.ui import main_window_init as main_window_main_init'),
    (r'from moleditpy\.modules\.view_2d import View2D',
     'from moleditpy.ui.zoomable_view import ZoomableView as View2D'),
    (r'from moleditpy\.modules\.constants import',        'from moleditpy.utils.constants import'),
    (r'from moleditpy\.modules\.main_window_ui_manager import MainWindowUiManager',
     'from moleditpy.ui.ui_manager import UIManager as MainWindowUiManager'),
    (r'from moleditpy\.modules\.main_window_app_state import MainWindowAppState',
     'from moleditpy.ui.app_state import StateManager as MainWindowAppState'),
    (r'from moleditpy\.modules\.main_window_export import MainWindowExport',
     'from moleditpy.ui.export_logic import ExportManager as MainWindowExport'),
    (r'from moleditpy\.modules\.main_window_edit_actions import MainWindowEditActions',
     'from moleditpy.ui.edit_actions_logic import EditActionsManager as MainWindowEditActions'),
    (r'from moleditpy\.modules\.molecule_scene import MoleculeScene',
     'from moleditpy.ui.molecule_scene import MoleculeScene'),

    # compute_manager
    (r'\bself\.mw\.create_atom_id_mapping\b',  'self.mw.compute_manager.create_atom_id_mapping'),
    (r'(?<!\w)mw\.create_atom_id_mapping\b',   'mw.compute_manager.create_atom_id_mapping'),

    # string_importer_manager (was main_window_string_importers)
    (r'\bmain_window_string_importers\b', 'string_importer_manager'),

    # --- mw.METHOD → mw.MANAGER.METHOD (remaining patterns) ---
    # state_manager
    (r'\bself\.mw\.data\b',                           'self.mw.state_manager.data'),
    (r'(?<!\w)mw\.data\b',                            'mw.state_manager.data'),
    (r'\bself\.mw\.reset_undo_stack\b',               'self.mw.state_manager.reset_undo_stack'),
    (r'(?<!\w)mw\.reset_undo_stack\b',                'mw.state_manager.reset_undo_stack'),
    (r'\bself\.mw\.create_json_data\b',               'self.mw.state_manager.create_json_data'),
    (r'(?<!\w)mw\.create_json_data\b',                'mw.state_manager.create_json_data'),
    (r'\bself\.mw\.load_from_json_data\b',            'self.mw.state_manager.load_from_json_data'),
    (r'(?<!\w)mw\.load_from_json_data\b',             'mw.state_manager.load_from_json_data'),
    # edit_actions_manager
    (r'\bself\.mw\.clear_2d_editor\b',                'self.mw.edit_actions_manager.clear_2d_editor'),
    (r'(?<!\w)mw\.clear_2d_editor\b',                 'mw.edit_actions_manager.clear_2d_editor'),
    (r'\bself\.mw\.clear_all\b',                      'self.mw.edit_actions_manager.clear_all'),
    (r'(?<!\w)mw\.clear_all\b',                       'mw.edit_actions_manager.clear_all'),
    (r'\bself\.mw\.undo\b',                           'self.mw.edit_actions_manager.undo'),
    (r'(?<!\w)mw\.undo\b',                            'mw.edit_actions_manager.undo'),
    (r'\bself\.mw\.redo\b',                           'self.mw.edit_actions_manager.redo'),
    (r'(?<!\w)mw\.redo\b',                            'mw.edit_actions_manager.redo'),
    # view_3d_manager
    (r'\bself\.mw\.draw_molecule_3d\b',               'self.mw.view_3d_manager.draw_molecule_3d'),
    (r'(?<!\w)mw\.draw_molecule_3d\b',                'mw.view_3d_manager.draw_molecule_3d'),
    (r'\bself\.mw\.fit_to_view\b',                    'self.mw.view_3d_manager.fit_to_view'),
    (r'(?<!\w)mw\.fit_to_view\b',                     'mw.view_3d_manager.fit_to_view'),
    (r'\bself\.mw\.glyph_source\b',                   'self.mw.view_3d_manager.glyph_source'),
    (r'(?<!\w)mw\.glyph_source\b',                    'mw.view_3d_manager.glyph_source'),
    (r'\bself\.mw\.atom_positions_3d\b',              'self.mw.view_3d_manager.atom_positions_3d'),
    (r'(?<!\w)mw\.atom_positions_3d\b',               'mw.view_3d_manager.atom_positions_3d'),
    (r'\bself\.mw\.apply_3d_settings\b',              'self.mw.view_3d_manager.apply_3d_settings'),
    (r'(?<!\w)mw\.apply_3d_settings\b',               'mw.view_3d_manager.apply_3d_settings'),
    (r'\bself\.mw\.update_chiral_labels\b',           'self.mw.view_3d_manager.update_chiral_labels'),
    (r'(?<!\w)mw\.update_chiral_labels\b',            'mw.view_3d_manager.update_chiral_labels'),
    (r'\bself\.mw\.setup_3d_hover\b',                 'self.mw.view_3d_manager.setup_3d_hover'),
    (r'(?<!\w)mw\.setup_3d_hover\b',                  'mw.view_3d_manager.setup_3d_hover'),
    (r'\bself\.mw\.current_3d_style\b',               'self.mw.view_3d_manager.current_3d_style'),
    (r'(?<!\w)mw\.current_3d_style\b',                'mw.view_3d_manager.current_3d_style'),
    # ui_manager
    (r'\bself\.mw\.set_mode\b',                       'self.mw.ui_manager.set_mode'),
    (r'(?<!\w)mw\.set_mode\b',                        'mw.ui_manager.set_mode'),
    (r'\bself\.mw\.restore_ui_for_editing\b',         'self.mw.ui_manager.restore_ui_for_editing'),
    (r'(?<!\w)mw\.restore_ui_for_editing\b',          'mw.ui_manager.restore_ui_for_editing'),
    (r'\bself\.mw\.toggle_3d_edit_mode\b',            'self.mw.ui_manager.toggle_3d_edit_mode'),
    (r'(?<!\w)mw\.toggle_3d_edit_mode\b',             'mw.ui_manager.toggle_3d_edit_mode'),
    (r'\bself\.mw\.activate_select_mode\b',           'self.mw.ui_manager.activate_select_mode'),
    (r'(?<!\w)mw\.activate_select_mode\b',            'mw.ui_manager.activate_select_mode'),
    # edit_3d_manager
    (r'\bself\.mw\.selected_atoms_3d\b',              'self.mw.edit_3d_manager.selected_atoms_3d'),
    (r'(?<!\w)mw\.selected_atoms_3d\b',               'mw.edit_3d_manager.selected_atoms_3d'),
    (r'\bself\.mw\.update_3d_selection_display\b',    'self.mw.edit_3d_manager.update_3d_selection_display'),
    (r'(?<!\w)mw\.update_3d_selection_display\b',     'mw.edit_3d_manager.update_3d_selection_display'),
    (r'\bself\.mw\.selected_atoms_for_measurement\b', 'self.mw.edit_3d_manager.selected_atoms_for_measurement'),
    (r'(?<!\w)mw\.selected_atoms_for_measurement\b',  'mw.edit_3d_manager.selected_atoms_for_measurement'),
    (r'\bself\.mw\.clear_2d_measurement_labels\b',    'self.mw.edit_3d_manager.clear_2d_measurement_labels'),
    (r'(?<!\w)mw\.clear_2d_measurement_labels\b',     'mw.edit_3d_manager.clear_2d_measurement_labels'),
    # io_manager
    (r'\bself\.mw\.estimate_bonds_from_distances\b',  'self.mw.io_manager.estimate_bonds_from_distances'),
    (r'(?<!\w)mw\.estimate_bonds_from_distances\b',   'mw.io_manager.estimate_bonds_from_distances'),
    # init_manager
    (r'\bself\.mw\.current_file_path\b',              'self.mw.init_manager.current_file_path'),
    (r'(?<!\w)mw\.current_file_path\b',               'mw.init_manager.current_file_path'),
    (r'\bself\.mw\.save_settings\b',                  'self.mw.init_manager.save_settings'),
    (r'(?<!\w)mw\.save_settings\b',                   'mw.init_manager.save_settings'),
    (r'\bself\.mw\.settings_dirty\b',                 'self.mw.init_manager.settings_dirty'),
    (r'(?<!\w)mw\.settings_dirty\b',                  'mw.init_manager.settings_dirty'),
    (r'\bself\.mw\.mode_actions\b',                   'self.mw.init_manager.mode_actions'),
    (r'(?<!\w)mw\.mode_actions\b',                    'mw.init_manager.mode_actions'),
    (r'\bself\.mw\.update_cpk_colors_from_settings\b', 'self.mw.init_manager.update_cpk_colors_from_settings'),
    (r'(?<!\w)mw\.update_cpk_colors_from_settings\b',  'mw.init_manager.update_cpk_colors_from_settings'),
]

# ---------------------------------------------------------------------------
# 2. Guard to inject at the TOP of every legacy run(mw) function
# ---------------------------------------------------------------------------
RUN_GUARD = "    if hasattr({p}, 'host'):\n        {p} = {p}.host\n"
RUN_GUARD_MARKER = "if hasattr({p}, 'host'):"


def inject_run_guard(text):
    """
    Find `def run(mw):` (or `def run(main_window):`) at module level and
    inject the MainInitManager→MainWindow guard as the first statement.
    """
    # Match def run(mw): or def run(main_window): with any body indented with 4 spaces
    pattern = re.compile(
        r'^(def run\((mw|main_window)\):[ \t]*\n)'
        r'((?:[ \t]+.*\n)*)',
        re.MULTILINE,
    )

    def replacer(m):
        sig = m.group(1)
        param = m.group(2)
        body = m.group(3)
        # Only inject if guard not already present
        marker = RUN_GUARD_MARKER.format(p=param)
        if marker in body:
            return m.group(0)
        guard = RUN_GUARD.format(p=param)
        return sig + guard + body

    return pattern.sub(replacer, text)


# ---------------------------------------------------------------------------
# 3. Main
# ---------------------------------------------------------------------------
def fix_file(path):
    with open(path, 'r', encoding='utf-8') as f:
        original = f.read()

    text = original
    counts = {}

    # Token replacements
    for pattern, replacement in GLOBAL_REPLACEMENTS:
        new_text, n = re.subn(pattern, replacement, text)
        if n:
            counts[pattern[:50]] = n
            text = new_text

    # run(mw) guard (for legacy plugins that have run())
    new_text = inject_run_guard(text)
    if new_text != text:
        counts['[run guard injected]'] = 1
        text = new_text

    # V3 plugin: has initialize() but no run() → add a no-op stub to silence
    # the "Missing attribute 'run'" error from the legacy plugin loader
    is_entry_point = (
        os.path.basename(path) == '__init__.py' or
        not any(sub in path for sub in [os.sep + sub + os.sep for sub in
                                        ['reaction_sketcher', 'orca_result_analyzer',
                                         'gaussian_fchk_mo_analyzer', 'orca_input_generator_pro',
                                         'nmr_predicator_nmrshiftdb2']])
    )
    has_initialize = bool(re.search(r'^def initialize\s*\(', text, re.MULTILINE))
    has_run = bool(re.search(r'^def run\s*\(', text, re.MULTILINE))
    if has_initialize and not has_run and is_entry_point:
        stub = ("\n\ndef run(mw):\n"
                "    \"\"\"Legacy compatibility stub — plugin uses initialize(context) API.\"\"\"\n"
                "    if hasattr(mw, 'host'):\n"
                "        mw = mw.host\n"
                "    pass\n")
        text = text + stub
        counts['[run stub added]'] = 1

    if text != original:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(text)
        total = sum(counts.values())
        rel = os.path.relpath(path, PLUGIN_BASE)
        print(f"  [{total:2d} changes] {rel}")
        for desc, n in counts.items():
            print(f"             {n}x  {desc}")
    return text != original


def main():
    py_files = sorted(glob.glob(os.path.join(PLUGIN_BASE, "**", "*.py"), recursive=True))
    changed = 0
    for path in py_files:
        if fix_file(path):
            changed += 1
    print(f"\nDone: {changed} file(s) modified out of {len(py_files)} scanned.")


if __name__ == "__main__":
    main()
