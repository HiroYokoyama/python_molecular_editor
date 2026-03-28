#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
Backward compatibility layer for moleditpy.modules
"""

import sys
import importlib
import types

# Define backward compatibility mapping
# Format: { 'old_module_name': 'new_full_import_path' }
_MAPPINGS = {
    # UI Components
    "about_dialog": "moleditpy.ui.about_dialog",
    "align_plane_dialog": "moleditpy.ui.align_plane_dialog",
    "alignment_dialog": "moleditpy.ui.alignment_dialog",
    "analysis_window": "moleditpy.ui.analysis_window",
    "angle_dialog": "moleditpy.ui.angle_dialog",
    "atom_item": "moleditpy.ui.atom_item",
    "bond_item": "moleditpy.ui.bond_item",
    "bond_length_dialog": "moleditpy.ui.bond_length_dialog",
    "color_settings_dialog": "moleditpy.ui.color_settings_dialog",
    "constrained_optimization_dialog": "moleditpy.ui.constrained_optimization_dialog",
    "custom_interactor_style": "moleditpy.ui.custom_interactor_style",
    "custom_qt_interactor": "moleditpy.ui.custom_qt_interactor",
    "dialog_3d_picking_mixin": "moleditpy.ui.dialog_3d_picking_mixin",
    "dialog_logic": "moleditpy.ui.dialog_logic",
    "dialog_manager": "moleditpy.ui.dialog_manager",
    "dihedral_dialog": "moleditpy.ui.dihedral_dialog",
    "main_window": "moleditpy.ui.main_window",
    "mirror_dialog": "moleditpy.ui.mirror_dialog",
    "molecular_scene_handler": "moleditpy.ui.molecular_scene_handler",
    "molecule_scene": "moleditpy.ui.molecule_scene",
    "move_group_dialog": "moleditpy.ui.move_group_dialog",
    "periodic_table_dialog": "moleditpy.ui.periodic_table_dialog",
    "planarize_dialog": "moleditpy.ui.planarize_dialog",
    "settings_dialog": "moleditpy.ui.settings_dialog",
    "template_preview_item": "moleditpy.ui.template_preview_item",
    "template_preview_view": "moleditpy.ui.template_preview_view",
    "translation_dialog": "moleditpy.ui.translation_dialog",
    "user_template_dialog": "moleditpy.ui.user_template_dialog",
    "zoomable_view": "moleditpy.ui.zoomable_view",
    "view_2d": "moleditpy.ui.zoomable_view",  # Alias for ZoomableView
    # Renamed UI Mixins
    "main_window_dialog_manager": "moleditpy.ui.dialog_manager",
    "main_window_edit_3d": "moleditpy.ui.edit_3d",
    "main_window_edit_actions": "moleditpy.ui.edit_actions",
    "main_window_export": "moleditpy.ui.export_logic",
    "main_window_main_init": "moleditpy.ui.main_window_init",
    "main_window_ui_manager": "moleditpy.ui.ui_manager",
    "main_window_view_3d": "moleditpy.ui.view_3d",
    "main_window_view_loaders": "moleditpy.ui.view_loaders",
    # Core Logic
    "calculation_worker": "moleditpy.ui.calculation_worker",
    "mol_geometry": "moleditpy.core.mol_geometry",
    "molecular_data": "moleditpy.core.molecular_data",
    "main_window_compute": "moleditpy.ui.compute_engine",
    "main_window_app_state": "moleditpy.ui.app_state",
    "main_window_molecular_parsers": "moleditpy.ui.molecular_parsers",
    "main_window_project_io": "moleditpy.ui.project_io",
    "main_window_string_importers": "moleditpy.ui.string_importers",
    # Utils
    "constants": "moleditpy.utils.constants",
    "sip_isdeleted_safe": "moleditpy.utils.sip_isdeleted_safe",
    "system_utils": "moleditpy.utils.system_utils",
    # Plugins
    "plugin_interface": "moleditpy.plugins.plugin_interface",
    "plugin_manager": "moleditpy.plugins.plugin_manager",
    "plugin_manager_window": "moleditpy.plugins.plugin_manager_window",
}


class _CompatibilityProxy(types.ModuleType):
    """Proxy module to handle imports on demand."""

    def __init__(self, name, target):
        super().__init__(name)
        self._target = target
        self._mod = None

    def _ensure_loaded(self):
        if self._mod is None:
            self._mod = importlib.import_module(self._target)
            # Handle class aliases (e.g. View2D)
            if self.__name__.endswith(".view_2d"):
                if not hasattr(self._mod, "View2D"):
                    setattr(
                        self._mod, "View2D", getattr(self._mod, "ZoomableView", None)
                    )
        return self._mod

    def __getattr__(self, name):
        return getattr(self._ensure_loaded(), name)

    def __dir__(self):
        return dir(self._ensure_loaded())


# Populate sys.modules and this package's namespace
for name, target in _MAPPINGS.items():
    proxy = _CompatibilityProxy(f"{__name__}.{name}", target)
    sys.modules[f"{__name__}.{name}"] = proxy
    setattr(sys.modules[__name__], name, proxy)

# Re-export key utilities directly in the package root
try:
    from moleditpy.utils.sip_isdeleted_safe import sip_isdeleted_safe
except ImportError:

    def sip_isdeleted_safe(obj):
        return False
