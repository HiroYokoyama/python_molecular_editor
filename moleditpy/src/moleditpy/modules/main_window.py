#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import traceback


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QMainWindow
)



from PyQt6.QtCore import (
    pyqtSignal, pyqtSlot
)


# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None
    
# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .main_window_app_state import MainWindowAppState
    from .main_window_compute import MainWindowCompute
    from .main_window_dialog_manager import MainWindowDialogManager
    from .main_window_edit_3d import MainWindowEdit3d
    from .main_window_edit_actions import MainWindowEditActions
    from .main_window_export import MainWindowExport
    from .main_window_main_init import MainWindowMainInit
    from .main_window_molecular_parsers import MainWindowMolecularParsers
    from .main_window_project_io import MainWindowProjectIo
    from .main_window_string_importers import MainWindowStringImporters
    from .main_window_ui_manager import MainWindowUiManager
    from .main_window_view_3d import MainWindowView3d
    from .main_window_view_loaders import MainWindowViewLoaders
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.main_window_app_state import MainWindowAppState
    from modules.main_window_compute import MainWindowCompute
    from modules.main_window_dialog_manager import MainWindowDialogManager
    from modules.main_window_edit_3d import MainWindowEdit3d
    from modules.main_window_edit_actions import MainWindowEditActions
    from modules.main_window_export import MainWindowExport
    from modules.main_window_main_init import MainWindowMainInit
    from modules.main_window_molecular_parsers import MainWindowMolecularParsers
    from modules.main_window_project_io import MainWindowProjectIo
    from modules.main_window_string_importers import MainWindowStringImporters
    from modules.main_window_ui_manager import MainWindowUiManager
    from modules.main_window_view_3d import MainWindowView3d
    from modules.main_window_view_loaders import MainWindowViewLoaders

class MainWindow(QMainWindow):

    # start_calculation carries the MOL block and an options object (second arg)
    start_calculation = pyqtSignal(str, object)
    def __init__(self, initial_file=None):
        # Initialize QMainWindow to ensure underlying Qt object is prepared
        # before helper classes interact with the window (fixes crash when
        # QMainWindow.__init__ not called).
        super().__init__()
        # --- MOVED TO main_window_main_init.py ---
        # The window functionality has been split across several helper
        # modules (main_window_*). Each helper is implemented as a class
        # whose methods were originally written to operate on the
        # MainWindow instance as 'self'. To maintain that behaviour we
        # create a small proxy (BoundFeature) that will forward call to
        # the helper class with the MainWindow instance as the first
        # argument.
        # Undo/Redo操作中に状態復元中であることを示すフラグ
        # 他のモジュールが呼び出される前に初期化する
        self._is_restoring_state = False

        class BoundFeature:
            """Bind a feature-class method calls to the MainWindow.

            Usage:
              self.main_window_view_3d = BoundFeature(MainWindowView3d, self)
              self.main_window_view_3d.set_3d_style('cpk')  # calls
              MainWindowView3d.set_3d_style(self, 'cpk')
            """

            def __init__(self, cls, host):
                self._cls = cls
                self._host = host

            def __getattr__(self, name):
                # Return a function which calls the class method with the
                # host (MainWindow instance) bound as first argument.
                attr = getattr(self._cls, name)
                if callable(attr):
                    return lambda *a, **k: attr(self._host, *a, **k)
                return attr

            def init(self, *a, **k):
                """Explicit initializer delegator for helper classes.

                Some helper classes define an __init__ method; when calling
                it via the BoundFeature instance, Python would resolve to the
                BoundFeature.__init__ method instead of the helper's.
                This helper method allows callers to forward initialization
                to the helper class implementation without touching _cls.
                """
                init_method = getattr(self._cls, '__init__', None)
                if callable(init_method):
                    return init_method(self._host, *a, **k)
                return None

        # Attach bound helpers to the main window; this keeps the
        # original call sites (e.g. self.main_window_compute.trigger_conversion)
        # working without changing them.
        self.main_window_main_init = BoundFeature(MainWindowMainInit, self)
        self.main_window_ui_manager = BoundFeature(MainWindowUiManager, self)
        self.main_window_view_3d = BoundFeature(MainWindowView3d, self)
        self.main_window_compute = BoundFeature(MainWindowCompute, self)
        self.main_window_edit_actions = BoundFeature(MainWindowEditActions, self)
        self.main_window_string_importers = BoundFeature(MainWindowStringImporters, self)
        self.main_window_molecular_parsers = BoundFeature(MainWindowMolecularParsers, self)
        self.main_window_view_loaders = BoundFeature(MainWindowViewLoaders, self)
        self.main_window_project_io = BoundFeature(MainWindowProjectIo, self)
        self.main_window_app_state = BoundFeature(MainWindowAppState, self)
        self.main_window_export = BoundFeature(MainWindowExport, self)
        self.main_window_dialog_manager = BoundFeature(MainWindowDialogManager, self)
        self.main_window_edit_3d = BoundFeature(MainWindowEdit3d, self)

        # Call the initialization method from main_window_main_init which
        # sets up the UI and the initial app state. Other helpers may
        # expect the UI to exist so initialize them after this call.
        try:
            # The helper's __init__ usually accepts the initial_file
            # argument and is expected to be invoked with the MainWindow
            # instance as the host. Because BoundFeature itself defines
            # __init__, attribute lookup would resolve to BoundFeature.__init__
            # instead of the helper class method. Call the helper class
            # __init__ explicitly with the main window instance as the
            # first argument.
            # Call the helper's __init__ via helper proxy to ensure the
            # initial host and arguments are forwarded correctly.
            self.main_window_main_init.init(initial_file)
        except Exception:
            # If main init fails, still continue so we can attempt
            # to catch / repair other issues via error messages.
            traceback.print_exc()

        # Initialize other helper modules which implement their own
        # __init__(self, main_window) if present.
        other_inits = [
            'main_window_view_3d', 'main_window_ui_manager', 'main_window_compute',
            'main_window_edit_actions', 'main_window_string_importers',
            'main_window_molecular_parsers', 'main_window_view_loaders',
            'main_window_project_io', 'main_window_app_state', 'main_window_export',
            'main_window_dialog_manager', 'main_window_edit_3d'
        ]
        for name in other_inits:
            try:
                # Call the helper's __init__ through the proxy helper to
                # avoid attribute-resolution issues.
                getattr(self, name).init()
            except Exception:
                # Ignore; many helpers only define __init__ when they
                # actually need it.
                pass

    def init_ui(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.init_ui()

    def init_menu_bar(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.init_menu_bar()

    def update_plugin_menu(self, plugin_menu):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.update_plugin_menu(plugin_menu)

    def init_worker_thread(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.init_worker_thread()

    def update_status_bar(self, message):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.update_status_bar(message)

    def set_mode(self, mode_str):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.set_mode(mode_str)

    def set_mode_and_update_toolbar(self, mode_str):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.set_mode_and_update_toolbar(mode_str)

    def set_3d_style(self, style_name):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.set_3d_style(style_name)

    def set_optimization_method(self, method_name):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.set_optimization_method(method_name)

    def copy_selection(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.copy_selection()

    def cut_selection(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.cut_selection()

    def paste_from_clipboard(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.paste_from_clipboard()

    def remove_hydrogen_atoms(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.remove_hydrogen_atoms()

    def add_hydrogen_atoms(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.add_hydrogen_atoms()

    def update_edit_menu_actions(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.update_edit_menu_actions()

    def show_convert_menu(self, pos):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.show_convert_menu(pos)

    def activate_select_mode(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.activate_select_mode()

    def _trigger_conversion_with_temp_mode(self, mode_key):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute._trigger_conversion_with_temp_mode(mode_key)

    def show_optimize_menu(self, pos):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.show_optimize_menu(pos)

    def _trigger_optimize_with_temp_method(self, method_key):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute._trigger_optimize_with_temp_method(method_key)

    def trigger_conversion(self):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.trigger_conversion()

    def halt_conversion(self):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.halt_conversion()

    def check_chemistry_problems_fallback(self):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.check_chemistry_problems_fallback()

    def optimize_3d_structure(self):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.optimize_3d_structure()

    def on_calculation_finished(self, result):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.on_calculation_finished(result)

    def create_atom_id_mapping(self):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.create_atom_id_mapping()

    @pyqtSlot(object)
    def on_calculation_error(self, result):
        # --- MOVED TO main_window_compute.py ---
        return self.main_window_compute.on_calculation_error(result)

    def eventFilter(self, obj, event):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.eventFilter(obj, event)

    def get_current_state(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.get_current_state()

    def set_state_from_data(self, state_data):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.set_state_from_data(state_data)

    def push_undo_state(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.push_undo_state()

    def update_window_title(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.update_window_title()

    def check_unsaved_changes(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.check_unsaved_changes()

    def reset_undo_stack(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.reset_undo_stack()

    def undo(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.undo()

    def redo(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.redo()

    def update_undo_redo_actions(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.update_undo_redo_actions()

    def update_realtime_info(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.update_realtime_info()

    def select_all(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.select_all()

    def show_about_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.show_about_dialog()

    def clear_all(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.clear_all()

    def clear_2d_editor(self, push_to_undo=True):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.clear_2d_editor(push_to_undo=True)

    def update_implicit_hydrogens(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.update_implicit_hydrogens()

    def import_smiles_dialog(self):
        # --- MOVED TO main_window_string_importers.py ---
        return self.main_window_string_importers.import_smiles_dialog()

    def import_inchi_dialog(self):
        # --- MOVED TO main_window_string_importers.py ---
        return self.main_window_string_importers.import_inchi_dialog()

    def load_from_smiles(self, smiles_string):
        # --- MOVED TO main_window_string_importers.py ---
        return self.main_window_string_importers.load_from_smiles(smiles_string)

    def load_from_inchi(self, inchi_string):
        # --- MOVED TO main_window_string_importers.py ---
        return self.main_window_string_importers.load_from_inchi(inchi_string)

    def fix_mol_counts_line(self, line: str) -> str:
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.fix_mol_counts_line(line)

    def fix_mol_block(self, mol_block: str) -> str:
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.fix_mol_block(mol_block)

    def load_mol_file(self, file_path=None):
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.load_mol_file(file_path)

    def load_mol_for_3d_viewing(self):
        # --- MOVED TO main_window_view_loaders.py ---
        return self.main_window_view_loaders.load_mol_for_3d_viewing()

    def load_xyz_for_3d_viewing(self, file_path=None):
        # --- MOVED TO main_window_view_loaders.py ---
        return self.main_window_view_loaders.load_xyz_for_3d_viewing(file_path)

    def load_xyz_file(self, file_path=None):
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.load_xyz_file(file_path)

    def estimate_bonds_from_distances(self, mol):
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.estimate_bonds_from_distances(mol)

    def save_project(self):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.save_project()

    def save_project_as(self):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.save_project_as()

    def save_raw_data(self):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.save_raw_data()

    def load_raw_data(self, file_path=None):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.load_raw_data(file_path)

    def save_as_json(self):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.save_as_json()

    def create_json_data(self):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.create_json_data()

    def load_json_data(self, file_path=None):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.load_json_data(file_path)

    def open_project_file(self, file_path=None):
        # --- MOVED TO main_window_project_io.py ---
        return self.main_window_project_io.open_project_file(file_path)

    def load_from_json_data(self, json_data):
        # --- MOVED TO main_window_app_state.py ---
        return self.main_window_app_state.load_from_json_data(json_data)

    def save_as_mol(self):
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.save_as_mol()

    def save_3d_as_mol(self):
        # --- MOVED TO main_window_view_loaders.py ---
        return self.main_window_view_loaders.save_3d_as_mol()

    def save_as_xyz(self):
        # --- MOVED TO main_window_molecular_parsers.py ---
        return self.main_window_molecular_parsers.save_as_xyz()

    def export_stl(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_stl()

    def export_obj_mtl(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_obj_mtl()

    def create_multi_material_obj(self, meshes_with_colors, obj_path, mtl_path):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.create_multi_material_obj(meshes_with_colors, obj_path, mtl_path)

    def export_color_stl(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_color_stl()

    def export_from_3d_view(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_from_3d_view()

    def export_from_3d_view_no_color(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_from_3d_view_no_color()

    def export_from_3d_view_with_colors(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_from_3d_view_with_colors()

    def export_2d_png(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_2d_png()

    def export_3d_png(self):
        # --- MOVED TO main_window_export.py ---
        return self.main_window_export.export_3d_png()

    def open_periodic_table_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_periodic_table_dialog()

    def set_atom_from_periodic_table(self, symbol): 
        self.set_mode(f'atom_{symbol}')

   
    def clean_up_2d_structure(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.clean_up_2d_structure()

    def resolve_overlapping_groups(self):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.resolve_overlapping_groups()

    def adjust_molecule_positions_to_avoid_collisions(self, mol, frags):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions.adjust_molecule_positions_to_avoid_collisions(mol, frags)

    def draw_molecule_3d(self, mol):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.draw_molecule_3d(mol)

    def _calculate_double_bond_offset(self, mol, bond, conf):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d._calculate_double_bond_offset(mol, bond, conf)

    def show_ez_labels_3d(self, mol):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.show_ez_labels_3d(mol)

    def toggle_chiral_labels_display(self, checked):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.toggle_chiral_labels_display(checked)

    def update_chiral_labels(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.update_chiral_labels()

    def toggle_atom_info_display(self, mode):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.toggle_atom_info_display(mode)

    def is_xyz_derived_molecule(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.is_xyz_derived_molecule()

    def has_original_atom_ids(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.has_original_atom_ids()

    def update_atom_id_menu_text(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.update_atom_id_menu_text()

    def update_atom_id_menu_state(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.update_atom_id_menu_state()

    def show_all_atom_info(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.show_all_atom_info()

    def clear_all_atom_info_labels(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.clear_all_atom_info_labels()

    def setup_3d_hover(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.setup_3d_hover()

    def open_analysis_window(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_analysis_window()

    def closeEvent(self, event):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.closeEvent(event)

    def zoom_in(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.zoom_in()

    def zoom_out(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.zoom_out()

    def reset_zoom(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.reset_zoom()

    def fit_to_view(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.fit_to_view()

    def draw_standard_3d_style(self, mol, style_override=None):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.draw_standard_3d_style(mol, style_override)

    def clear_measurement_selection(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.clear_measurement_selection()

    def toggle_3d_edit_mode(self, checked):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.toggle_3d_edit_mode(checked)

    def _setup_3d_picker(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager._setup_3d_picker()

    def _apply_chem_check_and_set_flags(self, mol, source_desc=None):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions._apply_chem_check_and_set_flags(mol, source_desc=None)

    def _clear_xyz_flags(self, mol=None):
        # --- MOVED TO main_window_edit_actions.py ---
        return self.main_window_edit_actions._clear_xyz_flags(mol=None)

    def load_mol_file_for_3d_viewing(self, file_path=None):
        # --- MOVED TO main_window_view_loaders.py ---
        return self.main_window_view_loaders.load_mol_file_for_3d_viewing(file_path)

    def load_command_line_file(self, file_path):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.load_command_line_file(file_path)

    def dragEnterEvent(self, event):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.dragEnterEvent(event)

    def dropEvent(self, event):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.dropEvent(event)

    def _enable_3d_edit_actions(self, enabled=True):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager._enable_3d_edit_actions(enabled=True)

    def _enable_3d_features(self, enabled=True):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager._enable_3d_features(enabled=True)

    def _enter_3d_viewer_ui_mode(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager._enter_3d_viewer_ui_mode()

    def restore_ui_for_editing(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.restore_ui_for_editing()

    def minimize_2d_panel(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.minimize_2d_panel()

    def restore_2d_panel(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.restore_2d_panel()

    def set_panel_layout(self, left_percent, right_percent):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.set_panel_layout(left_percent, right_percent)

    def toggle_2d_panel(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.toggle_2d_panel()

    def on_splitter_moved(self, pos, index):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.on_splitter_moved(pos, index)

    def open_template_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_template_dialog()

    def open_template_dialog_and_activate(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_template_dialog_and_activate()

    def save_2d_as_template(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.save_2d_as_template()

    def setup_splitter_tooltip(self):
        # --- MOVED TO main_window_ui_manager.py ---
        return self.main_window_ui_manager.setup_splitter_tooltip()

    def apply_initial_settings(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.apply_initial_settings()

    def update_cpk_colors_from_settings(self):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.update_cpk_colors_from_settings()

    def apply_3d_settings(self, redraw=True):
        # --- MOVED TO main_window_view_3d.py ---
        return self.main_window_view_3d.apply_3d_settings(redraw=True)

    def open_settings_dialog(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.open_settings_dialog()

    def reset_all_settings_menu(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.reset_all_settings_menu()

    def load_settings(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.load_settings()

    def save_settings(self):
        # --- MOVED TO main_window_main_init.py ---
        return self.main_window_main_init.save_settings()

    def toggle_measurement_mode(self, checked):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.toggle_measurement_mode(checked)

    def close_all_3d_edit_dialogs(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.close_all_3d_edit_dialogs()

    def handle_measurement_atom_selection(self, atom_idx):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.handle_measurement_atom_selection(atom_idx)

    def add_measurement_label(self, atom_idx, label_number):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.add_measurement_label(atom_idx, label_number)

    def update_measurement_labels_display(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.update_measurement_labels_display()

    def clear_measurement_selection(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.clear_measurement_selection()

    def update_2d_measurement_labels(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.update_2d_measurement_labels()

    def add_2d_measurement_label(self, atom_item, label_text):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.add_2d_measurement_label(atom_item, label_text)

    def clear_2d_measurement_labels(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.clear_2d_measurement_labels()

    def find_rdkit_atom_index(self, atom_item):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.find_rdkit_atom_index(atom_item)

    def calculate_and_display_measurements(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.calculate_and_display_measurements()

    def calculate_distance(self, atom1_idx, atom2_idx):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.calculate_distance(atom1_idx, atom2_idx)

    def calculate_angle(self, atom1_idx, atom2_idx, atom3_idx):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.calculate_angle(atom1_idx, atom2_idx, atom3_idx)

    def calculate_dihedral(self, atom1_idx, atom2_idx, atom3_idx, atom4_idx):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.calculate_dihedral(atom1_idx, atom2_idx, atom3_idx, atom4_idx)

    def display_measurement_text(self, measurement_lines):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.display_measurement_text(measurement_lines)

    def toggle_atom_selection_3d(self, atom_idx):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.toggle_atom_selection_3d(atom_idx)

    def clear_3d_selection(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.clear_3d_selection()

    def update_3d_selection_display(self):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.update_3d_selection_display()

    def planarize_selection(self, plane):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.planarize_selection(plane)

    def open_translation_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_translation_dialog()

    def open_move_group_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_move_group_dialog()

    def open_align_plane_dialog(self, plane):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_align_plane_dialog(plane)

    def open_planarize_dialog(self, plane=None):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_planarize_dialog(plane=None)

    def open_alignment_dialog(self, axis):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_alignment_dialog(axis)

    def open_bond_length_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_bond_length_dialog()

    def open_angle_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_angle_dialog()

    def open_dihedral_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_dihedral_dialog()

    def open_mirror_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_mirror_dialog()

    def open_constrained_optimization_dialog(self):
        # --- MOVED TO main_window_dialog_manager.py ---
        return self.main_window_dialog_manager.open_constrained_optimization_dialog()

    def remove_dialog_from_list(self, dialog):
        # --- MOVED TO main_window_edit_3d.py ---
        return self.main_window_edit_3d.remove_dialog_from_list(dialog)

