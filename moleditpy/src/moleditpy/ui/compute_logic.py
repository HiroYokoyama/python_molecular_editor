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
import logging
import contextlib
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from PyQt6.QtCore import QThread, QTimer, QPoint
from PyQt6.QtGui import QAction, QColor
from PyQt6.QtWidgets import QMenu, QMessageBox, QWidget
from rdkit import Chem

try:
    from . import OBABEL_AVAILABLE
except ImportError:
    from moleditpy.ui import OBABEL_AVAILABLE

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .calculation_worker import CalculationWorker
    from .mol_geometry import (
        identify_valence_problems,
        inject_ez_stereo_to_mol_block,
    )
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.calculation_worker import CalculationWorker
    from moleditpy.core.mol_geometry import (
        identify_valence_problems,
        inject_ez_stereo_to_mol_block,
    )


class ComputeManager:
    """Independent manager for molecular computations, ported from MainWindowCompute mixin."""

    def __init__(self, host: Any) -> None:
        self.host = host
        self.last_successful_optimization_method: Optional[str] = None
        self._active_calc_threads: List[QThread] = []
        self.halt_ids: Set[int] = set()
        self.next_conversion_id: int = 1
        self.active_worker_ids: Set[int] = set()

    def _safe_disconnect(self, signal: Any) -> None:
        """Safely disconnect a signal, silently ignoring RuntimeError."""
        try:
            signal.disconnect()
        except RuntimeError:
            pass

    def _remove_calculating_text(self) -> None:
        """Safely remove the 'Calculating...' text actor from the plotter."""
        actor = getattr(self, "_calculating_text_actor", None)
        if (
            actor
            and hasattr(self.host.view_3d_manager.plotter, "renderer")
            and self.host.view_3d_manager.plotter.renderer
        ):
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.host.view_3d_manager.plotter.renderer.RemoveActor(actor)
        with contextlib.suppress(AttributeError):
            if "_calculating_text_actor" in self.host.__dict__:
                delattr(self.host, "_calculating_text_actor")

    def _restore_button_ui(self) -> None:
        """Restore the Convert and Optimize buttons to their default state."""
        self._safe_disconnect(self.host.init_manager.convert_button.clicked)
        self.host.init_manager.convert_button.setText("Convert 2D to 3D")
        self.host.init_manager.convert_button.clicked.connect(self.trigger_conversion)
        self.host.init_manager.convert_button.setEnabled(True)

        if hasattr(self.host.init_manager, "optimize_3d_button"):
            self._safe_disconnect(self.host.init_manager.optimize_3d_button.clicked)
            self.host.init_manager.optimize_3d_button.setText("Optimize 3D")
            self.host.init_manager.optimize_3d_button.clicked.connect(
                self.optimize_3d_structure
            )
            self.host.init_manager.optimize_3d_button.setEnabled(True)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'optimize_3d_button' on object"
            )

    def _refresh_ui_state(self) -> None:
        """Consolidate UI state updates."""
        try:
            has_mol = self.host.view_3d_manager.current_mol is not None

            if hasattr(self.host.init_manager, "cleanup_button"):
                self.host.init_manager.cleanup_button.setEnabled(True)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'cleanup_button' on object"
                )

            self._restore_button_ui()

            if hasattr(self.host.init_manager, "optimize_3d_button"):
                self.host.init_manager.optimize_3d_button.setEnabled(has_mol)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'optimize_3d_button' on object"
                )
            if hasattr(self.host.init_manager, "export_button"):
                self.host.init_manager.export_button.setEnabled(has_mol)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'export_button' on object"
                )

            # ui_manager and its methods are guaranteed on the host
            self.host.ui_manager._enable_3d_features(has_mol)
            self.host.ui_manager._enable_3d_edit_actions(has_mol)

            if hasattr(self.host.init_manager, "analysis_action"):
                self.host.init_manager.analysis_action.setEnabled(has_mol)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'analysis_action' on object"
                )
            if hasattr(self.host.init_manager, "edit_3d_action"):
                self.host.init_manager.edit_3d_action.setEnabled(has_mol)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    "REPORT ERROR: Missing attribute 'edit_3d_action' on object"
                )

            # plotter and view_2d are fundamental host components
            if self.host.view_3d_manager.plotter:
                self.host.view_3d_manager.plotter.render()
            if self.host.init_manager.view_2d:
                self.host.init_manager.view_2d.setFocus()
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(f"Non-critical UI refresh error: {e}")

    def set_optimization_method(self, method_name: str) -> None:
        """Set preferred 3D optimization method."""
        if not method_name:
            return
        method = str(method_name).strip().upper()
        self.host.init_manager.optimization_method = method
        self.host.init_manager.settings["optimization_method"] = method
        self.host.init_manager.settings_dirty = True

        if (
            hasattr(self.host.init_manager, "opt3d_actions")
            and self.host.init_manager.opt3d_actions
        ):
            for k, act in self.host.init_manager.opt3d_actions.items():
                act.setChecked(k.upper() == method)

        label = self.host.init_manager.opt3d_method_labels.get(method, method)
        self.host.statusBar().showMessage(f"3D optimization method set to: {label}")

    def toggle_intermolecular_interaction_rdkit(self, checked: bool) -> None:
        """Toggle intermolecular interactions for RDKit optimization."""
        self.host.init_manager.settings["optimize_intermolecular_interaction_rdkit"] = (
            checked
        )
        self.host.init_manager.settings_dirty = True
        state_str = "Enabled" if checked else "Disabled"
        self.host.statusBar().showMessage(
            f"Intermolecular interaction for RDKit: {state_str}"
        )

    def show_convert_menu(self, pos: QPoint) -> None:
        """Temporary 3D conversion menu (right-click)."""
        if not self.host.init_manager.convert_button.isEnabled():
            return

        menu = QMenu(self.host)
        conv_options = [
            ("Fallback", "fallback"),
            ("RDKit only", "rdkit"),
            ("Open Babel only", "obabel"),
            ("Direct", "direct"),
        ]
        for label, key in conv_options:
            a = QAction(label, self.host)
            if key == "obabel" and not OBABEL_AVAILABLE:
                a.setEnabled(False)
            a.triggered.connect(
                lambda checked=False, k=key: self._trigger_conversion_with_temp_mode(k)
            )
            menu.addAction(a)
        menu.exec(self.host.init_manager.convert_button.mapToGlobal(pos))

    def _trigger_conversion_with_temp_mode(self, mode_key: str) -> None:
        self.host._temp_conv_mode = mode_key
        QTimer.singleShot(0, self.trigger_conversion)

    def show_optimize_menu(self, pos: QPoint) -> None:
        """Temporary 3D optimization menu (right-click)."""
        if not self.host.init_manager.optimize_3d_button.isEnabled():
            return

        menu = QMenu(self.host)
        opt_list = [
            ("MMFF94s (RDKit)", "MMFF_RDKIT"),
            ("MMFF94 (RDKit)", "MMFF94_RDKIT"),
            ("UFF (RDKit)", "UFF_RDKIT"),
            ("MMFF94s (Open Babel)", "MMFF94s_OBABEL"),
            ("MMFF94 (Open Babel)", "MMFF94_OBABEL"),
            ("UFF (Open Babel)", "UFF_OBABEL"),
            ("GAFF (Open Babel)", "GAFF_OBABEL"),
            ("Ghemical (Open Babel)", "GHEMICAL_OBABEL"),
        ]
        for label, key in opt_list:
            a = QAction(label, self.host)
            if (
                hasattr(self.host.init_manager, "opt3d_actions")
                and key in self.host.init_manager.opt3d_actions
            ):
                a.setEnabled(self.host.init_manager.opt3d_actions[key].isEnabled())
            a.triggered.connect(
                lambda checked=False, k=key: self._trigger_optimize_with_temp_method(k)
            )
            menu.addAction(a)
        menu.exec(self.host.init_manager.optimize_3d_button.mapToGlobal(pos))

    def _trigger_optimize_with_temp_method(self, method_key: str) -> None:
        self.host._temp_optimization_method = method_key
        QTimer.singleShot(0, self.optimize_3d_structure)

    def trigger_conversion(self) -> None:
        """Main entry point for 2D to 3D conversion."""
        self.last_successful_optimization_method = None
        self.host.edit_3d_manager.constraints_3d = []

        if not self.host.state_manager.data.atoms:
            self.host.view_3d_manager.plotter.clear()
            self.host.view_3d_manager.current_mol = None
            self.host.init_manager.analysis_action.setEnabled(False)
            self.host.statusBar().showMessage("3D view cleared.")
            self.host.init_manager.view_2d.setFocus()
            return

        # Reset modes
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)
        if self.host.edit_3d_manager.is_3d_edit_mode:
            self.host.init_manager.edit_3d_action.setChecked(False)
            self.host.ui_manager.toggle_3d_edit_mode(False)

        mol = self._prepare_rdkit_mol_for_conversion()
        if not mol:
            self._restore_button_ui()
            self.host.init_manager.cleanup_button.setEnabled(True)
            self._remove_calculating_text()
            return

        num_frags = len(Chem.GetMolFrags(mol))
        if num_frags > 1:
            self.host.statusBar().showMessage(
                f"Converting {num_frags} molecules to 3D with collision detection..."
            )
        else:
            self.host.statusBar().showMessage("Calculating 3D structure...")

        mol_block = self._setup_mol_block_for_worker(mol)

        run_id = self.next_conversion_id
        self.next_conversion_id = run_id + 1
        self.active_worker_ids.add(run_id)

        # UI Updates
        self.host.init_manager.convert_button.setText("Halt conversion")
        self._safe_disconnect(self.host.init_manager.convert_button.clicked)
        self.host.init_manager.convert_button.clicked.connect(self.halt_conversion)
        self.host.init_manager.cleanup_button.setEnabled(False)
        self.host.ui_manager._enable_3d_features(False)
        self.host.view_3d_manager.plotter.clear()
        self.host.view_3d_manager.current_mol = None

        # Add 'Calculating...' overlay
        bg_qcolor = QColor(
            self.host.init_manager.settings.get("background_color", "#919191")
        )
        text_color = (
            "black"
            if (bg_qcolor.isValid() and bg_qcolor.toHsl().lightness() > 128)
            else "white"
        )
        self.host.compute_manager._calculating_text_actor = (
            self.host.view_3d_manager.plotter.add_text(
                "Calculating...",
                position="lower_right",
                font_size=15,
                color=text_color,
                name="calculating_text",
            )
        )
        self.host.view_3d_manager.plotter.render()

        options = {
            "conversion_mode": self.host.__dict__.pop(
                "_temp_conv_mode",
                self.host.init_manager.settings.get("3d_conversion_mode", "fallback"),
            ),
            "optimization_method": self.host.__dict__.pop(
                "_temp_optimization_method", None
            )
            or getattr(self.host.init_manager, "optimization_method", "MMFF_RDKIT"),
            "optimize_intermolecular_interaction_rdkit": self.host.init_manager.settings.get(
                "optimize_intermolecular_interaction_rdkit", True
            ),
            "worker_id": run_id,
        }

        self._start_calculation_worker(mol_block, options, run_id)
        self.host.view_3d_manager.update_chiral_labels()
        self.host.init_manager.view_2d.setFocus()

    def halt_conversion(self) -> None:
        """Halt the in-progress conversion."""
        wids_to_halt = set(self.active_worker_ids)
        if wids_to_halt:
            self.halt_ids.update(wids_to_halt)
        self.active_worker_ids.clear()

        self._restore_button_ui()
        self.host.init_manager.cleanup_button.setEnabled(True)
        self._remove_calculating_text()
        self.host.statusBar().showMessage("Halted")

    def optimize_3d_structure(self) -> None:
        """Optimize 3D structure."""
        if not self.host.view_3d_manager.current_mol:
            self.host.statusBar().showMessage("No 3D molecule to optimize.")
            return

        if self.host.view_3d_manager.current_mol.GetNumConformers() == 0:
            self.host.statusBar().showMessage(
                "No conformer found. Generate 3D structure first."
            )
            return

        method = self.host.__dict__.pop("_temp_optimization_method", None) or getattr(
            self.host.init_manager, "optimization_method", "MMFF_RDKIT"
        )
        method = method.upper() if method else "MMFF_RDKIT"

        # Validate method against known labels (from init_manager) and registered plugins
        _init_mgr = getattr(self.host, "init_manager", None)
        _init_methods = set(getattr(_init_mgr, "opt3d_method_labels", None) or {})
        _plugin_mgr = getattr(self.host, "plugin_manager", None)
        _plugin_methods = set(getattr(_plugin_mgr, "optimization_methods", {}) or {})
        _all_known = _init_methods | _plugin_methods | {"OPTIMIZE_ONLY"}
        if _all_known and method not in _all_known:
            self.host.statusBar().showMessage(
                f"Selected optimization method '{method}' is not available."
            )
            return

        self.host.statusBar().showMessage(f"Optimizing 3D structure ({method})...")

        mol_block = Chem.MolToMolBlock(
            self.host.view_3d_manager.current_mol, includeStereo=True
        )
        options = {
            "conversion_mode": "optimize_only",
            "optimization_method": method,
            "optimize_intermolecular_interaction_rdkit": self.host.init_manager.settings.get(
                "optimize_intermolecular_interaction_rdkit", True
            ),
        }

        run_id = int(getattr(self, "next_conversion_id", 1))
        self.next_conversion_id = run_id + 1
        options["worker_id"] = run_id

        if hasattr(self.host, "active_worker_ids"):
            self.host.compute_manager.active_worker_ids.add(run_id)
        else:
            self.host.compute_manager.active_worker_ids = {run_id}

        if hasattr(self.host.init_manager, "optimize_3d_button"):
            self.host.init_manager.optimize_3d_button.setText("Halt optimize")
            self._safe_disconnect(self.host.init_manager.optimize_3d_button.clicked)
            self.host.init_manager.optimize_3d_button.clicked.connect(
                self.halt_conversion
            )
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'optimize_3d_button' on object"
            )

        self.host.ui_manager._enable_3d_features(False)
        # Re-enable the button so it can be clicked to Halt
        if hasattr(self.host.init_manager, "optimize_3d_button"):
            self.host.init_manager.optimize_3d_button.setEnabled(True)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'optimize_3d_button' on object"
            )

        self._start_calculation_worker(mol_block, options, run_id)

    def _prepare_rdkit_mol_for_conversion(self) -> Optional[Any]:
        mol = self.host.state_manager.data.to_rdkit_mol(use_2d_stereo=False)
        if not mol or mol.GetNumAtoms() == 0:
            self.check_chemistry_problems_fallback()
            return None

        self.original_atom_properties = {}
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            with contextlib.suppress(KeyError):
                self.original_atom_properties[i] = atom.GetIntProp("_original_atom_id")

        problems = Chem.DetectChemistryProblems(mol)
        if problems:
            self._handle_chemistry_problems(mol, problems)
            return None

        self.host.init_manager.scene.clear_all_problem_flags()
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            logging.error(f"Sanitization failed: {e}")
            self.host.statusBar().showMessage("Error: Invalid chemical structure.")
            return None
        return mol

    def _handle_chemistry_problems(self, mol: Any, problems: List[Any]) -> None:
        self.host.init_manager.scene.clear_all_problem_flags()
        msg = f"Error: {len(problems)} chemistry problem(s) found (e.g., hypervalency). Fix the 2D layout before converting."
        self.host.statusBar().showMessage(msg)

        with contextlib.suppress(Exception):
            QMessageBox.critical(self.host, "Chemistry Problem", msg)

        for prob in problems:
            with contextlib.suppress(Exception):
                atom_idx = prob.GetAtomIdx()
                rd_atom = mol.GetAtomWithIdx(atom_idx)
                orig_id = rd_atom.GetIntProp("_original_atom_id")
                item = self.host.state_manager.data.atoms[orig_id]["item"]
                item.has_problem = True
                item.update()
        self.host.init_manager.view_2d.setFocus()

    def _setup_mol_block_for_worker(self, mol: Any) -> str:
        mol_block = self.host.state_manager.data.to_mol_block()
        if not mol_block:
            mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
        return inject_ez_stereo_to_mol_block(
            mol_block, mol, self.host.state_manager.data.bonds
        )

    def _start_calculation_worker(
        self, mol_block: str, options: Dict[str, Any], run_id: int
    ) -> None:
        thread = QThread()
        worker = CalculationWorker()
        worker.halt_ids = self.halt_ids
        worker.moveToThread(thread)
        worker.status_update.connect(self.host.ui_manager.update_status_bar)

        def _cleanup() -> None:
            thread.quit()
            thread.finished.connect(thread.deleteLater)
            worker.deleteLater()
            with contextlib.suppress(ValueError, AttributeError):
                self._active_calc_threads.remove(thread)

        def _on_finished(result: Any) -> None:
            try:
                self.on_calculation_finished(result)
            finally:
                _cleanup()

        def _on_error(msg: Any) -> None:
            try:
                self.on_calculation_error(msg)
            finally:
                _cleanup()

        worker.finished.connect(_on_finished)
        worker.error.connect(_on_error)
        thread.start()
        QTimer.singleShot(10, lambda: worker.start_work.emit(mol_block, options))
        self._active_calc_threads.append(thread)

    def on_calculation_finished(
        self, result: Union[Chem.Mol, Tuple[int, Chem.Mol]]
    ) -> None:
        worker_id, mol = result if isinstance(result, tuple) else (None, result)
        if worker_id is not None:
            if worker_id not in self.active_worker_ids:
                # Still show something if the user wants to see 'staled' result
                self.host.statusBar().showMessage(
                    f"Ignored halted worker result (ID = {worker_id})"
                )
                return  # stale worker, ignore
            self.active_worker_ids.discard(worker_id)
            self.halt_ids.discard(worker_id)

        self.host.view_3d_manager.current_mol = mol
        self.host.is_xyz_derived = False

        # Restore properties
        if hasattr(self, "original_atom_properties") and mol:
            for i, orig_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    mol.GetAtomWithIdx(i).SetIntProp("_original_atom_id", orig_id)

        self.create_atom_id_mapping()
        self.host.view_3d_manager.update_chiral_labels()
        self.host.view_3d_manager.draw_molecule_3d(mol)
        self._remove_calculating_text()
        self._refresh_ui_state()
        self.host.edit_actions_manager.push_undo_state()
        self.host.view_3d_manager.plotter.reset_camera()

        # Record the successful optimization method from mol property or current setting
        try:
            method_key = None
            if mol and mol.HasProp("_pme_optimization_method"):
                method_key = mol.GetProp("_pme_optimization_method")
            if not method_key:
                method_key = getattr(self, "optimization_method", None) or getattr(
                    self.host.init_manager, "optimization_method", None
                )
            if method_key:
                labels = getattr(self.host.init_manager, "opt3d_method_labels", {})
                self.last_successful_optimization_method = labels.get(
                    method_key, method_key
                )
        except (AttributeError, TypeError):
            pass

        # Update 3D interactive features (hovers, menus)
        self.host.view_3d_manager.setup_3d_hover()
        self.host.view_3d_manager.update_atom_id_menu_text()
        self.host.view_3d_manager.update_atom_id_menu_state()

    def on_calculation_error(self, message: Union[str, Tuple[int, str]]) -> None:
        # Accept either a string or (worker_id, message) tuple from the worker signal
        if isinstance(message, tuple) and len(message) == 2:
            worker_id, msg = message
            if worker_id not in self.active_worker_ids:
                # Still cleanup overlay/buttons even if stale
                self._remove_calculating_text()
                self._restore_button_ui()
                # If it was a halt or error from a stale worker, show with ID
                if msg == "Halt" or msg == "Halted":
                    self.host.statusBar().showMessage(
                        f"Ignored halted worker (ID = {worker_id})"
                    )
                else:
                    self.host.statusBar().showMessage(
                        f"Ignored stale worker error: {msg} (ID = {worker_id})"
                    )
                return  # stale worker, ignore
            self.active_worker_ids.discard(worker_id)
        else:
            msg = str(message)

        self._remove_calculating_text()
        self._restore_button_ui()

        if msg == "Halt" or msg == "Halted":
            self.host.statusBar().showMessage("Halted")
        else:
            self.host.statusBar().showMessage(f"Calculation Error: {msg}")

        # Offer UFF fallback when MMFF fails
        if "MMFF" in msg.upper() and "fail" in msg.lower():
            try:
                reply = QMessageBox.question(
                    self.host if isinstance(self.host, QWidget) else None,
                    "Retry with UFF?",
                    f"{msg}\n\nWould you like to retry using UFF (RDKit) instead?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if reply == QMessageBox.StandardButton.Yes:
                    self._temp_optimization_method = "UFF_RDKIT"
                    self.optimize_3d_structure()
                    return
            except (TypeError, RuntimeError):
                pass

        with contextlib.suppress(TypeError, RuntimeError):
            QMessageBox.critical(self.host, "Calculation Error", msg)

    def create_atom_id_mapping(self) -> None:
        """Map 2D atom IDs to 3D RDKit indices."""
        if not self.host.view_3d_manager.current_mol:
            return

        self.host.atom_id_to_rdkit_idx_map = {}

        # Create mapping from RDKit properties
        for i in range(self.host.view_3d_manager.current_mol.GetNumAtoms()):
            rdkit_atom = self.host.view_3d_manager.current_mol.GetAtomWithIdx(i)
            try:
                original_atom_id = rdkit_atom.GetIntProp("_original_atom_id")
                self.host.atom_id_to_rdkit_idx_map[original_atom_id] = i
            except KeyError:
                # Skip if property missing
                continue

    def check_chemistry_problems_fallback(self) -> None:
        problem_atom_ids = identify_valence_problems(
            self.host.state_manager.data.atoms, self.host.state_manager.data.bonds
        )
        if problem_atom_ids:
            for aid in problem_atom_ids:
                item = self.host.state_manager.data.atoms[aid].get("item")
                if item:
                    item.has_problem = True
                    item.update()
            self.host.statusBar().showMessage(
                f"Error: {len(problem_atom_ids)} chemistry problems found."
            )
        self.host.init_manager.view_2d.setFocus()

    def update_aromatic_rings(self) -> None:
        """Update aromatic ring visualization."""
        try:
            mol = self.host.state_manager.data.to_rdkit_mol()
            if not mol:
                return
            Chem.SanitizeMol(mol)
            ri = mol.GetRingInfo()
            aromatic_rings = []
            for ring in ri.AtomRings():
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    aromatic_rings.append(ring)
            self.host.init_manager.scene.set_aromatic_rings(aromatic_rings)
        except Exception as e:
            logging.debug(f"Aromatic ring update failed: {e}")

    def select_connected_atoms(self) -> None:
        """Select all atoms connected to the current selection."""
        selected_items = self.host.init_manager.scene.selectedItems()
        if not selected_items:
            return

        atom_ids = set()
        for item in selected_items:
            if hasattr(item, "atom_id"):
                atom_ids.add(item.atom_id)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error("REPORT ERROR: Missing attribute 'atom_id' on item")

        if not atom_ids:
            return

        connected = set(atom_ids)
        stack = list(atom_ids)
        while stack:
            curr = stack.pop()
            for neighbor in self.host.state_manager.data.adjacency_list.get(curr, []):
                if neighbor not in connected:
                    connected.add(neighbor)
                    stack.append(neighbor)

        for aid in connected:
            item = self.host.state_manager.data.atoms[aid]["item"]
            item.setSelected(True)
