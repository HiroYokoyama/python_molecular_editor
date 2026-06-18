#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit import Chem

# PyQt6 Modules
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QMainWindow


try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .app_state import StateManager
    from .compute_logic import ComputeManager
    from .dialog_logic import DialogManager
    from .edit_3d_logic import Edit3DManager
    from .edit_actions_logic import EditActionsManager
    from .io_logic import IOManager
    from .export_logic import ExportManager
    from .main_window_init import MainInitManager
    from .string_importers import StringImporterManager
    from .ui_manager import UIManager
    from .view_3d_logic import View3DManager
    from .molecule_scene import MoleculeScene
    from ..core.molecular_data import MolecularData
    from .custom_qt_interactor import CustomQtInteractor
except (AttributeError, RuntimeError, TypeError, ImportError):
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.app_state import StateManager
    from moleditpy.ui.compute_logic import ComputeManager
    from moleditpy.ui.dialog_logic import DialogManager
    from moleditpy.ui.edit_3d_logic import Edit3DManager
    from moleditpy.ui.edit_actions_logic import EditActionsManager
    from moleditpy.ui.io_logic import IOManager
    from moleditpy.ui.export_logic import ExportManager
    from moleditpy.ui.main_window_init import MainInitManager
    from moleditpy.ui.string_importers import StringImporterManager
    from moleditpy.ui.ui_manager import UIManager
    from moleditpy.ui.view_3d_logic import View3DManager
    from moleditpy.ui.molecule_scene import MoleculeScene
    from moleditpy.core.molecular_data import MolecularData
    from moleditpy.ui.custom_qt_interactor import CustomQtInteractor


class MainWindow(QMainWindow):
    # start_calculation carries the MOL block and an options object (second arg)
    start_calculation = pyqtSignal(str, object)

    def __init__(
        self, initial_file: Optional[str] = None, safe_mode: bool = False
    ) -> None:
        QMainWindow.__init__(self)

        # Initialize properties
        self.is_restoring_state = False

        # Type annotations for static analysis / mypy
        self.initial_settings: Any
        self.is_xyz_derived: bool
        self.chem_check_tried: bool
        self.chem_check_failed: bool
        self.template_dialog: Any
        self.picking_consumed: bool
        self.initialization_complete: bool
        self.ih_update_counter: int
        self.plugin_manager: Any
        self.redraw_menu_action: Any
        self.translation_action: Any
        self.move_selected_atoms_action: Any
        self.move_group_action: Any
        self.align_menu: Any
        self.align_x_action: Any
        self.align_y_action: Any
        self.align_z_action: Any
        self.alignplane_xy_action: Any
        self.alignplane_xz_action: Any
        self.alignplane_yz_action: Any
        self.mirror_action: Any
        self.planarize_action: Any
        self.bond_length_action: Any
        self.angle_action: Any
        self.dihedral_action: Any
        self.constrained_opt_action: Any
        self.intermolecular_rdkit_action: Any
        self.atom_id_to_rdkit_idx_map: dict[int, int]

        # Initialize ad-hoc attributes at startup.
        self.initial_settings = None
        self.is_xyz_derived = False
        self.chem_check_tried = False
        self.chem_check_failed = False
        self.template_dialog = None
        self.picking_consumed = False
        self.initialization_complete = False
        self.ih_update_counter = 0
        self.plugin_manager = None
        self.redraw_menu_action = None
        self.translation_action = None
        self.move_selected_atoms_action = None
        self.move_group_action = None
        self.align_menu = None
        self.align_x_action = None
        self.align_y_action = None
        self.align_z_action = None
        self.alignplane_xy_action = None
        self.alignplane_xz_action = None
        self.alignplane_yz_action = None
        self.mirror_action = None
        self.planarize_action = None
        self.bond_length_action = None
        self.angle_action = None
        self.dihedral_action = None
        self.constrained_opt_action = None
        self.intermolecular_rdkit_action = None
        self.atom_id_to_rdkit_idx_map = {}

        self.export_manager = ExportManager(self)
        self.view_3d_manager = View3DManager(self)
        self.edit_3d_manager = Edit3DManager(self)
        self.edit_actions_manager = EditActionsManager(self)
        self.compute_manager = ComputeManager(self)
        self.dialog_manager = DialogManager(self)
        self.io_manager = IOManager(self)
        self.state_manager = StateManager(self)
        self.string_importer_manager = StringImporterManager(self)
        self.ui_manager = UIManager(self)

        self.init_manager = MainInitManager(
            self, initial_file=initial_file, safe_mode=safe_mode
        )

    # --- Mediator Methods ---

    def set_current_molecule(self, mol: Optional[Chem.Mol]) -> None:
        """Set the current 3D molecule and trigger rendering updates."""
        self.view_3d_manager.current_mol = mol

    def set_3d_atom_positions(self, positions: Any) -> None:
        """Update 3D atom positions."""
        self.view_3d_manager.atom_positions_3d = positions

    def clear_3d_view(self) -> None:
        """Clear 3D molecule state and plotter rendering."""
        self.view_3d_manager.current_mol = None
        if self.view_3d_manager.plotter:
            self.view_3d_manager.plotter.clear()

    def set_plotter_camera_position(self, camera_position: Any) -> None:
        """Set 3D view camera position."""
        if self.view_3d_manager.plotter:
            self.view_3d_manager.plotter.camera_position = camera_position

    def set_plotter_picker(self, picker: Any) -> None:
        """Set picker object on 3D view plotter."""
        if self.view_3d_manager.plotter:
            self.view_3d_manager.plotter.picker = picker

    def set_constraints_3d(self, constraints: List[Any]) -> None:
        """Update list of 3D geometry constraints."""
        self.edit_3d_manager.constraints_3d = constraints

    def get_constraints_3d(self) -> List[Any]:
        """Get the current list of 3D geometry constraints."""
        return self.edit_3d_manager.constraints_3d

    def set_has_unsaved_changes(self, value: bool) -> None:
        """Update project's unsaved changes status."""
        self.state_manager.has_unsaved_changes = value

    def set_current_file_path(self, path: Optional[str]) -> None:
        """Update the path of the currently open project file."""
        self.init_manager.current_file_path = path

    def get_current_file_path(self) -> Optional[str]:
        """Get the path of the currently open project file."""
        return self.init_manager.current_file_path

    def get_molecule_data(self) -> MolecularData:
        """Get current 2D/3D molecular data representation."""
        return self.state_manager.data

    def set_molecule_data(self, data: MolecularData) -> None:
        """Set molecular data and synchronize it with the 2D scene."""
        self.state_manager.data = data
        if self.init_manager.scene:
            self.init_manager.scene.data = data

    def set_settings_dirty(self, value: bool) -> None:
        """Mark settings as dirty."""
        self.init_manager.settings_dirty = value

    def set_is_2d_editable(self, value: bool) -> None:
        """Set whether 2D editor is editable."""
        self.ui_manager.is_2d_editable = value

    def set_optimization_method(self, method: str) -> None:
        """Set the active 3D optimization method."""
        self.init_manager.optimization_method = method

    def set_3d_edit_mode(self, enabled: bool) -> None:
        """Enable or disable 3D editing mode."""
        self.edit_3d_manager.is_3d_edit_mode = enabled

    def is_3d_measurement_mode(self) -> bool:
        """Check if 3D measurement mode is active."""
        return bool(self.edit_3d_manager.measurement_mode)

    def set_scene_mode(self, mode: str) -> None:
        """Set mode on 2D editor scene."""
        if self.init_manager.scene:
            self.init_manager.scene.mode = mode

    def set_scene_atom_symbol(self, symbol: str) -> None:
        """Set active atom symbol on 2D editor scene."""
        if self.init_manager.scene:
            self.init_manager.scene.current_atom_symbol = symbol

    def set_scene_bond_properties(self, order: int, stereo: int = 0) -> None:
        """Set active bond order and stereochemistry on 2D editor scene."""
        if self.init_manager.scene:
            self.init_manager.scene.bond_order = order
            self.init_manager.scene.bond_stereo = stereo

    def set_scene_user_template_data(self, data: Any) -> None:
        """Set template data on 2D editor scene."""
        if self.init_manager.scene:
            self.init_manager.scene.user_template_data = data

    def update_status_message(self, message: str, timeout: int = 0) -> None:
        """Show a message in the main window's status bar."""
        if self.statusBar():
            self.statusBar().showMessage(message, timeout)

    def set_last_successful_optimization_method(self, method: Optional[str]) -> None:
        """Set the last successful 3D optimization method."""
        self.compute_manager.last_successful_optimization_method = method

    def get_settings(self) -> Dict[str, Any]:
        """Get the settings dictionary."""
        return self.init_manager.settings

    def set_atom_id_to_rdkit_idx_map(self, mapping: Dict[int, int]) -> None:
        """Set 2D-to-3D atom ID mapping."""
        self.view_3d_manager.atom_id_to_rdkit_idx_map = mapping

    def save_state_snapshot(self) -> None:
        """Create a deep copy snapshot of the current state for undo/redo comparison."""
        import copy

        try:
            self.state_manager.saved_state = copy.deepcopy(
                self.state_manager.get_current_state()
            )
        except Exception:
            # Safe defensive fallback catching Exception
            pass

    def update_window_title(self) -> None:
        """Update main window title to reflect file path and save status."""
        self.state_manager.update_window_title()

    def set_show_chiral_labels(self, value: bool) -> None:
        """Set show chiral labels on 3D viewer."""
        self.view_3d_manager.show_chiral_labels = value

    def set_plotter(self, plotter: Optional[CustomQtInteractor]) -> None:
        """Set the 3D viewer plotter."""
        self.view_3d_manager.plotter = plotter

    def reset_active_calc_threads(self) -> None:
        """Reset active calculation threads list."""
        self.compute_manager.reset_active_threads()

    # --- Core Proxy Properties (Legacy Plugin Support Only. Bypassed by Core Logics) ---
    @property
    def current_mol(self) -> Optional[Chem.Mol]:
        """Proxy for current molecule. Not for core logic use."""
        return self.view_3d_manager.current_mol

    @current_mol.setter
    def current_mol(self, value: Any) -> None:
        """Proxy for current molecule setter. Not for core logic use."""
        self.view_3d_manager.current_mol = value

    @property
    def plotter(self) -> Optional[CustomQtInteractor]:
        """Proxy for 3D plotter. Not for core logic use."""
        return self.view_3d_manager.plotter

    @property
    def data(self) -> MolecularData:
        """Proxy for state data. Not for core logic use."""
        return self.state_manager.data  # type: ignore[return-value, no-any-return]

    @property
    def scene(self) -> Optional[MoleculeScene]:
        """Proxy for 2D scene. Not for core logic use."""
        return self.init_manager.scene

    def draw_molecule_3d(self, mol: Any) -> None:
        """Proxy for 3D rendering. Not for core logic use."""
        self.current_mol = mol
        self.view_3d_manager.draw_molecule_3d(mol)
