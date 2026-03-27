#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import contextlib
from PyQt6.QtCore import QThread, QTimer
from PyQt6.QtGui import QAction, QColor
from PyQt6.QtWidgets import QMenu, QMessageBox
from rdkit import Chem

try:
    from . import OBABEL_AVAILABLE
except ImportError:
    from moleditpy.ui import OBABEL_AVAILABLE

try:
    # package relative imports
    from .calculation_worker import CalculationWorker
    from .mol_geometry import (
        identify_valence_problems,
        inject_ez_stereo_to_mol_block,
    )
except ImportError:
    from moleditpy.core.calculation_worker import CalculationWorker
    from moleditpy.core.mol_geometry import (
        identify_valence_problems,
        inject_ez_stereo_to_mol_block,
    )


class ComputeManager:
    """Independent manager for molecular computations, ported from MainWindowCompute mixin."""

    def __init__(self, host):
        self.host = host
        self.last_successful_optimization_method = None
        self._active_calc_threads = []
        self.halt_ids = set()

    def __getattr__(self, name):
        """Delegate back to host for attributes not found on this manager."""
        if name == "host":
            raise AttributeError(name)
        return getattr(self.host, name)

    def _safe_disconnect(self, signal):
        """Safely disconnect a signal, silently ignoring RuntimeError."""
        try:
            signal.disconnect()
        except RuntimeError:
            pass

    def _remove_calculating_text(self):
        """Safely remove the 'Calculating...' text actor from the plotter."""
        actor = getattr(self.host, "_calculating_text_actor", None)
        if (
            actor
            and hasattr(self.host.plotter, "renderer")
            and self.host.plotter.renderer
        ):
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.host.plotter.renderer.RemoveActor(actor)
        if hasattr(self.host, "_calculating_text_actor"):
            delattr(self.host, "_calculating_text_actor")

    def _restore_button_ui(self):
        """Restore the Convert and Optimize buttons to their default state."""
        self._safe_disconnect(self.host.convert_button.clicked)
        self.host.convert_button.setText("Convert 2D to 3D")
        self.host.convert_button.clicked.connect(self.host.trigger_conversion)
        self.host.convert_button.setEnabled(True)

        if hasattr(self.host, "optimize_3d_button"):
            self._safe_disconnect(self.host.optimize_3d_button.clicked)
            self.host.optimize_3d_button.setText("Optimize 3D")
            self.host.optimize_3d_button.clicked.connect(
                self.host.optimize_3d_structure
            )
            self.host.optimize_3d_button.setEnabled(True)

    def _refresh_ui_state(self):
        """Consolidate UI state updates."""
        try:
            has_mol = self.host.current_mol is not None

            if hasattr(self.host, "cleanup_button"):
                self.host.cleanup_button.setEnabled(True)

            self._restore_button_ui()

            if hasattr(self.host, "optimize_3d_button"):
                self.host.optimize_3d_button.setEnabled(has_mol)
            if hasattr(self.host, "export_button"):
                self.host.export_button.setEnabled(has_mol)

            if hasattr(self.host, "_enable_3d_features"):
                self.host._enable_3d_features(has_mol)
            if hasattr(self.host, "_enable_3d_edit_actions"):
                self.host._enable_3d_edit_actions(has_mol)

            if hasattr(self.host, "analysis_action"):
                self.host.analysis_action.setEnabled(has_mol)
            if hasattr(self.host, "edit_3d_action"):
                self.host.edit_3d_action.setEnabled(has_mol)

            self.host.plotter.render()
            if self.host.view_2d:
                self.host.view_2d.setFocus()
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(f"Non-critical UI refresh error: {e}")

    def set_optimization_method(self, method_name):
        """Set preferred 3D optimization method."""
        if not method_name:
            return
        method = str(method_name).strip().upper()
        self.host.optimization_method = method
        self.host.settings["optimization_method"] = method
        self.host.settings_dirty = True

        if hasattr(self.host, "opt3d_actions") and self.host.opt3d_actions:
            for k, act in self.host.opt3d_actions.items():
                act.setChecked(k.upper() == method)

        label = getattr(self.host, "opt3d_method_labels", {}).get(method, method)
        self.host.statusBar().showMessage(f"3D optimization method set to: {label}")

    def toggle_intermolecular_interaction_rdkit(self, checked):
        """Toggle intermolecular interactions for RDKit optimization."""
        self.host.settings["optimize_intermolecular_interaction_rdkit"] = checked
        self.host.settings_dirty = True
        state_str = "Enabled" if checked else "Disabled"
        self.host.statusBar().showMessage(
            f"Intermolecular interaction for RDKit: {state_str}"
        )

    def show_convert_menu(self, pos):
        """Temporary 3D conversion menu (right-click)."""
        if not self.host.convert_button.isEnabled():
            return

        menu = QMenu(self.host)
        conv_options = [
            ("RDKit -> Open Babel -> Direct (fallback)", "fallback"),
            ("RDKit only", "rdkit"),
            ("Open Babel only", "obabel"),
            ("Direct (use 2D coords + add H)", "direct"),
        ]
        for label, key in conv_options:
            a = QAction(label, self.host)
            if key == "obabel" and not OBABEL_AVAILABLE:
                a.setEnabled(False)
            a.triggered.connect(
                lambda checked=False, k=key: self._trigger_conversion_with_temp_mode(k)
            )
            menu.addAction(a)
        menu.exec(self.host.convert_button.mapToGlobal(pos))

    def _trigger_conversion_with_temp_mode(self, mode_key):
        self.host._temp_conv_mode = mode_key
        QTimer.singleShot(0, self.host.trigger_conversion)

    def show_optimize_menu(self, pos):
        """Temporary 3D optimization menu (right-click)."""
        if not self.host.optimize_3d_button.isEnabled():
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
            if hasattr(self.host, "opt3d_actions") and key in self.host.opt3d_actions:
                a.setEnabled(self.host.opt3d_actions[key].isEnabled())
            a.triggered.connect(
                lambda checked=False, k=key: self._trigger_optimize_with_temp_method(k)
            )
            menu.addAction(a)
        menu.exec(self.host.optimize_3d_button.mapToGlobal(pos))

    def _trigger_optimize_with_temp_method(self, method_key):
        self.host._temp_optimization_method = method_key
        QTimer.singleShot(0, self.host.optimize_3d_structure)

    def trigger_conversion(self):
        """Main entry point for 2D to 3D conversion."""
        self.last_successful_optimization_method = None
        self.host.constraints_3d = []

        if not self.host.data.atoms:
            self.host.plotter.clear()
            self.host.current_mol = None
            self.host.analysis_action.setEnabled(False)
            self.host.statusBar().showMessage("3D view cleared.")
            self.host.view_2d.setFocus()
            return

        # Reset modes
        if self.host.measurement_mode:
            self.host.measurement_action.setChecked(False)
            self.host.toggle_measurement_mode(False)
        if self.host.is_3d_edit_mode:
            self.host.edit_3d_action.setChecked(False)
            self.host.toggle_3d_edit_mode(False)

        mol = self._prepare_rdkit_mol_for_conversion()
        if not mol:
            return

        num_frags = len(Chem.GetMolFrags(mol))
        if num_frags > 1:
            self.host.statusBar().showMessage(
                f"Converting {num_frags} molecules to 3D with collision detection..."
            )
        else:
            self.host.statusBar().showMessage("Calculating 3D structure...")

        mol_block = self._setup_mol_block_for_worker(mol)

        run_id = int(getattr(self.host, "next_conversion_id", 1))
        self.host.next_conversion_id = run_id + 1
        if not hasattr(self.host, "active_worker_ids"):
            self.host.active_worker_ids = set()
        self.host.active_worker_ids.add(run_id)

        # UI Updates
        self.host.convert_button.setText("Halt conversion")
        self._safe_disconnect(self.host.convert_button.clicked)
        self.host.convert_button.clicked.connect(self.halt_conversion)
        self.host.cleanup_button.setEnabled(False)
        self.host._enable_3d_features(False)
        self.host.plotter.clear()
        self.host.current_mol = None

        # Add 'Calculating...' overlay
        bg_qcolor = QColor(self.host.settings.get("background_color", "#919191"))
        text_color = (
            "black"
            if (bg_qcolor.isValid() and bg_qcolor.toHsl().lightness() > 128)
            else "white"
        )
        self.host._calculating_text_actor = self.host.plotter.add_text(
            "Calculating...",
            position="lower_right",
            font_size=15,
            color=text_color,
            name="calculating_text",
        )
        self.host.plotter.render()

        options = {
            "conversion_mode": self.host.__dict__.pop(
                "_temp_conv_mode",
                self.host.settings.get("3d_conversion_mode", "fallback"),
            ),
            "optimization_method": self.host.__dict__.pop(
                "_temp_optimization_method", None
            )
            or getattr(self.host, "optimization_method", "MMFF_RDKIT"),
            "optimize_intermolecular_interaction_rdkit": self.host.settings.get(
                "optimize_intermolecular_interaction_rdkit", True
            ),
            "worker_id": run_id,
        }

        self._start_calculation_worker(mol_block, options, run_id)
        self.host.push_undo_state()
        self.host.update_chiral_labels()
        self.host.view_2d.setFocus()

    def halt_conversion(self):
        """Halt the in-progress conversion."""
        wids_to_halt = set(getattr(self.host, "active_worker_ids", set()))
        if wids_to_halt:
            self.halt_ids.update(wids_to_halt)
        if hasattr(self.host, "active_worker_ids"):
            self.host.active_worker_ids.clear()

        self._restore_button_ui()
        self.host.cleanup_button.setEnabled(True)
        self._remove_calculating_text()
        self.host.statusBar().showMessage("Halted")

    def optimize_3d_structure(self):
        """Optimize 3D structure."""
        if not self.host.current_mol:
            self.host.statusBar().showMessage("No 3D molecule to optimize.")
            return

        method = self.host.__dict__.pop("_temp_optimization_method", None) or getattr(
            self.host, "optimization_method", "MMFF_RDKIT"
        )
        method = method.upper() if method else "MMFF_RDKIT"

        self.host.statusBar().showMessage(f"Optimizing 3D structure ({method})...")

        mol_block = Chem.MolToMolBlock(self.host.current_mol, includeStereo=True)
        options = {
            "conversion_mode": "optimize_only",
            "optimization_method": method,
            "optimize_intermolecular_interaction_rdkit": self.host.settings.get(
                "optimize_intermolecular_interaction_rdkit", True
            ),
        }

        run_id = int(getattr(self.host, "next_conversion_id", 1))
        self.host.next_conversion_id = run_id + 1
        options["worker_id"] = run_id

        if hasattr(self.host, "active_worker_ids"):
            self.host.active_worker_ids.add(run_id)
        else:
            self.host.active_worker_ids = {run_id}

        if hasattr(self.host, "optimize_3d_button"):
            self.host.optimize_3d_button.setText("Halt optimize")
            self._safe_disconnect(self.host.optimize_3d_button.clicked)
            self.host.optimize_3d_button.clicked.connect(self.halt_conversion)

        self.host._enable_3d_features(False)
        self._start_calculation_worker(mol_block, options, run_id)

    def _prepare_rdkit_mol_for_conversion(self):
        mol = self.host.data.to_rdkit_mol(use_2d_stereo=False)
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

        self.host.scene.clear_all_problem_flags()
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            logging.error(f"Sanitization failed: {e}")
            self.host.statusBar().showMessage("Error: Invalid chemical structure.")
            return None
        return mol

    def _handle_chemistry_problems(self, mol, problems):
        self.host.scene.clear_all_problem_flags()
        self.host.statusBar().showMessage(
            f"Error: {len(problems)} chemistry problem(s) found."
        )
        for prob in problems:
            with contextlib.suppress(Exception):
                atom_idx = prob.GetAtomIdx()
                rd_atom = mol.GetAtomWithIdx(atom_idx)
                orig_id = rd_atom.GetIntProp("_original_atom_id")
                item = self.host.data.atoms[orig_id]["item"]
                item.has_problem = True
                item.update()
        self.host.view_2d.setFocus()

    def _setup_mol_block_for_worker(self, mol):
        mol_block = self.host.data.to_mol_block()
        if not mol_block:
            mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
        return inject_ez_stereo_to_mol_block(mol_block, mol, self.host.data.bonds)

    def _start_calculation_worker(self, mol_block, options, run_id):
        thread = QThread()
        worker = CalculationWorker()
        worker.halt_ids = self.halt_ids
        worker.moveToThread(thread)
        worker.status_update.connect(self.host.update_status_bar)

        def _cleanup():
            thread.quit()
            thread.finished.connect(thread.deleteLater)
            worker.deleteLater()
            with contextlib.suppress(ValueError, AttributeError):
                self._active_calc_threads.remove(thread)

        def _on_finished(result):
            try:
                self.on_calculation_finished(result)
            finally:
                _cleanup()

        def _on_error(msg):
            try:
                self.on_calculation_error(msg)
            finally:
                _cleanup()

        worker.finished.connect(_on_finished)
        worker.error.connect(_on_error)
        thread.start()
        QTimer.singleShot(10, lambda: worker.start_work.emit(mol_block, options))
        self._active_calc_threads.append(thread)

    def on_calculation_finished(self, result):
        worker_id, mol = result if isinstance(result, tuple) else (None, result)
        if worker_id is not None:
            active_ids = getattr(self.host, "active_worker_ids", set())
            if worker_id not in active_ids:
                self._remove_calculating_text()
                self._restore_button_ui()
                return
            active_ids.discard(worker_id)
            self.halt_ids.discard(worker_id)

        self.host.current_mol = mol
        self.host.is_xyz_derived = False

        # Restore properties
        if hasattr(self, "original_atom_properties") and mol:
            for i, orig_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    mol.GetAtomWithIdx(i).SetIntProp("_original_atom_id", orig_id)

        self.host.create_atom_id_mapping()
        self.host.update_chiral_labels()
        self.host.draw_molecule_3d(mol)
        self._remove_calculating_text()
        self._refresh_ui_state()
        self.host.push_undo_state()
        self.host.plotter.reset_camera()

    def on_calculation_error(self, message):
        self._remove_calculating_text()
        self._restore_button_ui()
        self.host.statusBar().showMessage(f"Calculation Error: {message}")
        QMessageBox.critical(self.host, "Calculation Error", message)

    def create_atom_id_mapping(self):
        """Map 2D atom IDs to 3D RDKit indices."""
        if not self.host.current_mol:
            return

        self.host.atom_id_to_rdkit_idx_map = {}

        # Create mapping from RDKit properties
        for i in range(self.host.current_mol.GetNumAtoms()):
            rdkit_atom = self.host.current_mol.GetAtomWithIdx(i)
            try:
                original_atom_id = rdkit_atom.GetIntProp("_original_atom_id")
                self.host.atom_id_to_rdkit_idx_map[original_atom_id] = i
            except KeyError:
                # Skip if property missing
                continue

    def check_chemistry_problems_fallback(self):
        problem_atom_ids = identify_valence_problems(
            self.host.data.atoms, self.host.data.bonds
        )
        if problem_atom_ids:
            for aid in problem_atom_ids:
                item = self.host.data.atoms[aid].get("item")
                if item:
                    item.has_problem = True
                    item.update()
            self.host.statusBar().showMessage(
                f"Error: {len(problem_atom_ids)} chemistry problems found."
            )
        self.host.view_2d.setFocus()

    def update_aromatic_rings(self):
        """Update aromatic ring visualization."""
        try:
            mol = self.host.data.to_rdkit_mol()
            if not mol:
                return
            Chem.SanitizeMol(mol)
            ri = mol.GetRingInfo()
            aromatic_rings = []
            for ring in ri.AtomRings():
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    aromatic_rings.append(ring)
            self.host.scene.set_aromatic_rings(aromatic_rings)
        except Exception as e:
            logging.debug(f"Aromatic ring update failed: {e}")

    def select_connected_atoms(self):
        """Select all atoms connected to the current selection."""
        selected_items = self.host.scene.selectedItems()
        if not selected_items:
            return

        atom_ids = set()
        for item in selected_items:
            if hasattr(item, "atom_id"):
                atom_ids.add(item.atom_id)

        if not atom_ids:
            return

        connected = set(atom_ids)
        stack = list(atom_ids)
        while stack:
            curr = stack.pop()
            for neighbor in self.host.data.adjacency_list.get(curr, []):
                if neighbor not in connected:
                    connected.add(neighbor)
                    stack.append(neighbor)

        for aid in connected:
            item = self.host.data.atoms[aid]["item"]
            item.setSelected(True)
