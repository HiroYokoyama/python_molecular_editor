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

from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit import Chem

# PyQt6 Modules
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QMainWindow

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None  # type: ignore[assignment]
    _sip_isdeleted = None

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
        self._is_restoring_state = False

        # Type annotations for static analysis / mypy
        self.initial_settings: Any
        self.is_xyz_derived: bool
        self.chem_check_tried: bool
        self.chem_check_failed: bool
        self._template_dialog: Any
        self._picking_consumed: bool
        self.initialization_complete: bool
        self._ih_update_counter: int
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
        self._temp_conv_mode: Any
        self._temp_optimization_method: Any

        # Initialize ad-hoc attributes directly in __dict__ at startup.
        # This bypasses class-level property descriptors mocked by testing fixtures
        # (before managers exist) and prevents overwriting values set during manager instantiation.
        self.__dict__["initial_settings"] = None
        self.__dict__["is_xyz_derived"] = False
        self.__dict__["chem_check_tried"] = False
        self.__dict__["chem_check_failed"] = False
        self.__dict__["_template_dialog"] = None
        self.__dict__["_picking_consumed"] = False
        self.__dict__["initialization_complete"] = False
        self.__dict__["_ih_update_counter"] = 0
        self.__dict__["plugin_manager"] = None
        self.__dict__["redraw_menu_action"] = None
        self.__dict__["translation_action"] = None
        self.__dict__["move_selected_atoms_action"] = None
        self.__dict__["move_group_action"] = None
        self.__dict__["align_menu"] = None
        self.__dict__["align_x_action"] = None
        self.__dict__["align_y_action"] = None
        self.__dict__["align_z_action"] = None
        self.__dict__["alignplane_xy_action"] = None
        self.__dict__["alignplane_xz_action"] = None
        self.__dict__["alignplane_yz_action"] = None
        self.__dict__["mirror_action"] = None
        self.__dict__["planarize_action"] = None
        self.__dict__["bond_length_action"] = None
        self.__dict__["angle_action"] = None
        self.__dict__["dihedral_action"] = None
        self.__dict__["constrained_opt_action"] = None
        self.__dict__["intermolecular_rdkit_action"] = None
        self.__dict__["atom_id_to_rdkit_idx_map"] = {}
        self.__dict__["_temp_conv_mode"] = None
        self.__dict__["_temp_optimization_method"] = None

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
