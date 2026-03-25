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
# Module separated from MainWindow (main_window.py)
# Functional class: MainWindowCompute
"""

# RDKit imports (explicit to satisfy flake8 and used features)
import logging
from PyQt6.QtCore import QThread, QTimer
from PyQt6.QtGui import QAction, QColor

# PyQt6 Modules
from PyQt6.QtWidgets import QApplication, QMenu, QMessageBox
import contextlib
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from . import OBABEL_AVAILABLE
except ImportError:
    from modules import OBABEL_AVAILABLE


try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except (AttributeError, RuntimeError, TypeError):
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .calculation_worker import CalculationWorker
    from .mol_geometry import is_problematic_valence
except ImportError:
    # Fallback to absolute imports for script-style execution
    from modules.calculation_worker import CalculationWorker
    from modules.mol_geometry import is_problematic_valence


class MainWindowCompute(object):
    """Functional class separated from main_window.py"""
    
    # Default initial state 
    last_successful_optimization_method = None

    # ------------------------------------------------------------------
    # Private helpers for safe state management
    # ------------------------------------------------------------------

    def _safe_disconnect(self, signal):
        """Safely disconnect a signal, silently ignoring RuntimeError if nothing is connected."""
        try:
            signal.disconnect()
        except RuntimeError:
            pass  # nothing was connected, ignore

    def _remove_calculating_text(self):
        """Safely remove the 'Calculating...' text actor from the plotter."""
        actor = getattr(self, "_calculating_text_actor", None)
        if actor and hasattr(self.plotter, "renderer") and self.plotter.renderer:
            # Suppress potential errors if the actor or renderer is already destroyed during teardown
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.plotter.renderer.RemoveActor(actor)
        # Remove attribute safely without risk of AttributeError
        self.__dict__.pop("_calculating_text_actor", None)

    def _restore_button_ui(self):
        """Restore the Convert and Optimize buttons to their default state."""
        self._safe_disconnect(self.convert_button.clicked)
        self.convert_button.setText("Convert 2D to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        self.convert_button.setEnabled(True)

        if hasattr(self, "optimize_3d_button"):
            self._safe_disconnect(self.optimize_3d_button.clicked)
            self.optimize_3d_button.setText("Optimize 3D")
            self.optimize_3d_button.clicked.connect(self.optimize_3d_structure)
            self.optimize_3d_button.setEnabled(True)

    # ------------------------------------------------------------------

    def set_optimization_method(self, method_name):  
        """Set preferred 3D optimization method and persist to settings.

        Supported values: 'GAFF', 'MMFF'
        """
        # Normalize input and validate
        if not method_name:
            return
        method = str(method_name).strip().upper()
        valid_methods = (
            "MMFF_RDKIT",
            "MMFF94_RDKIT",
            "UFF_RDKIT",
            "UFF_OBABEL",
            "GAFF_OBABEL",
            "MMFF94s_OBABEL",
            "MMFF94_OBABEL",
            "GHEMICAL_OBABEL",
        )
        if method not in valid_methods:
            # Unknown method: ignore but notify
            self.statusBar().showMessage(
                f"Unknown 3D optimization method: {method_name}"
            )
            return

        # Update internal state (store canonical uppercase key)
        self.optimization_method = method

        # Persist to settings
        self.settings["optimization_method"] = self.optimization_method
        self.settings_dirty = True

        # Update menu checked state if actions mapping exists
        if hasattr(self, "opt3d_actions") and self.opt3d_actions:
            for k, act in self.opt3d_actions.items():
                act.setChecked(k.upper() == method)

        # Also show user-friendly label if available
        label = getattr(self, "opt3d_method_labels", {}).get(
            self.optimization_method, self.optimization_method
        )
        self.statusBar().showMessage(f"3D optimization method set to: {label}")

    def toggle_intermolecular_interaction_rdkit(self, checked):  
        """Toggle whether intermolecular interactions are considered for RDKit optimization."""
        self.settings["optimize_intermolecular_interaction_rdkit"] = checked
        self.settings_dirty = True
        state_str = "Enabled" if checked else "Disabled"
        self.statusBar().showMessage(f"Intermolecular interaction for RDKit: {state_str}")
    
    def show_convert_menu(self, pos):  
        """Temporary 3D conversion menu (right-click). Not persisted."""
        # If button is disabled (during calculation), do not show menu
        if not self.convert_button.isEnabled():
            return            
        
        menu = QMenu(self)
        conv_options = [
            ("RDKit -> Open Babel -> Direct (fallback)", "fallback"),
            ("RDKit only", "rdkit"),
            ("Open Babel only", "obabel"),
            ("Direct (use 2D coords + add H)", "direct"),
        ]
        for label, key in conv_options:
            a = QAction(label, self)
            if key == "obabel" and not OBABEL_AVAILABLE:
                a.setEnabled(False)
            a.triggered.connect(
                lambda checked=False, k=key: self._trigger_conversion_with_temp_mode(k)
            )
            menu.addAction(a)

        # Show menu at button position
        if hasattr(self, "convert_button") and _sip_isdeleted and not _sip_isdeleted(self.convert_button):
            menu.exec_(self.convert_button.mapToGlobal(pos))

    def _trigger_conversion_with_temp_mode(self, mode_key):  
        # store temporary override and invoke conversion
        self._temp_conv_mode = mode_key
        # Call the normal conversion entry point (it will consume the temp)
        QTimer.singleShot(0, self.trigger_conversion)

    def show_optimize_menu(self, pos):  
        """Temporary 3D optimization menu (right-click). Not persisted."""
        # If button is disabled (during calculation), do not show menu
        if not self.optimize_3d_button.isEnabled():
            return

        menu = QMenu(self)
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
            a = QAction(label, self)
            # If opt3d_actions exist, reflect their enabled state
            if hasattr(self, "opt3d_actions") and key in self.opt3d_actions:
                a.setEnabled(self.opt3d_actions[key].isEnabled())
            a.triggered.connect(
                lambda checked=False,
                k=key: self._trigger_optimize_with_temp_method(k)
            )
            menu.addAction(a)

        # Add Plugin Optimization Methods
        if (
            getattr(self, "plugin_manager", None)
            and getattr(self.plugin_manager, "optimization_methods", None)
        ):
            methods = self.plugin_manager.optimization_methods
            if methods:
                menu.addSeparator()
                for method_name, info in methods.items():
                    a = QAction(info.get("label", method_name), self)
                    a.triggered.connect(
                        lambda checked=False,
                        k=method_name: self._trigger_optimize_with_temp_method(k)
                    )
                    menu.addAction(a)

        # Show menu at button position
        if hasattr(self, "optimize_3d_button") and _sip_isdeleted and not _sip_isdeleted(self.optimize_3d_button):
            menu.exec_(self.optimize_3d_button.mapToGlobal(pos))

    def _trigger_optimize_with_temp_method(self, method_key):  
        # store temporary override and invoke optimization
        self._temp_optimization_method = method_key
        # Run optimize on next event loop turn so UI updates first
        QTimer.singleShot(0, self.optimize_3d_structure)

    def _prepare_rdkit_mol_for_conversion(self):
        """Convert current data to RDKit mol and check for chemistry problems."""
        mol = self.data.to_rdkit_mol(use_2d_stereo=False)

        # Check chemistry if mol object creation fails
        if not mol or mol.GetNumAtoms() == 0:
            # Run fallback chemistry check
            self.check_chemistry_problems_fallback()
            return None

        # Save atom properties (lost in worker)
        self.original_atom_properties = {}
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            try:
                original_id = atom.GetIntProp("_original_atom_id")
                self.original_atom_properties[i] = original_id
            except KeyError:
                continue

        problems = Chem.DetectChemistryProblems(mol)
        if problems:
            self.main_window_compute._handle_chemistry_problems(mol, problems)
            return None

        # Clear flags and run 3D conversion if no problems
        self.scene.clear_all_problem_flags()

        try:
            Chem.SanitizeMol(mol)
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.error(f"Sanitization failed: {e}")
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus()
            return None

        return mol

    def _handle_chemistry_problems(self, mol, problems):
        """Mark chemistry problems in the 2D scene."""
        self.scene.clear_all_problem_flags()
        self.statusBar().showMessage(f"Error: {len(problems)} chemistry problem(s) found.")
        self.scene.clearSelection()

        for prob in problems:
            try:
                atom_idx = prob.GetAtomIdx()
                rdkit_atom = mol.GetAtomWithIdx(atom_idx)
                if rdkit_atom.HasProp("_original_atom_id"):
                    original_id = rdkit_atom.GetIntProp("_original_atom_id")
                    if original_id in self.data.atoms:
                        item = self.data.atoms[original_id]["item"]
                        if item:
                            item.has_problem = True
                            item.update()
            except (AttributeError, RuntimeError, ValueError, TypeError):
                continue
        self.view_2d.setFocus()

    def _setup_mol_block_for_worker(self, mol):
        """Generate high-quality MOL block with E/Z stereo for the worker."""
        mol_block = self.data.to_mol_block()
        if not mol_block:
            mol_block = Chem.MolToMolBlock(mol, includeStereo=True)

        mol_lines = mol_block.split("\n")
        ez_bond_info = {}
        for (id1, id2), bond_data in self.data.bonds.items():
            if bond_data.get("stereo") in [3, 4]:  # E/Z labels
                rdkit_idx1, rdkit_idx2 = None, None
                for atom in mol.GetAtoms():
                    if atom.HasProp("_original_atom_id"):
                        orig_id = atom.GetIntProp("_original_atom_id")
                        if orig_id == id1: rdkit_idx1 = atom.GetIdx()
                        elif orig_id == id2: rdkit_idx2 = atom.GetIdx()

                if rdkit_idx1 is not None and rdkit_idx2 is not None:
                    rdkit_bond = mol.GetBondBetweenAtoms(rdkit_idx1, rdkit_idx2)
                    if rdkit_bond and rdkit_bond.GetBondType() == Chem.BondType.DOUBLE:
                        ez_bond_info[rdkit_bond.GetIdx()] = bond_data["stereo"]

        if ez_bond_info:
            insert_idx = len(mol_lines) - 1  # Before M  END
            for bond_idx, stereo_type in ez_bond_info.items():
                cfg_value = 1 if stereo_type == 3 else 2  # 1=Z, 2=E in MOL format
                cfg_line = f"M  CFG  1 {bond_idx + 1:3d}   {cfg_value}"
                mol_lines.insert(insert_idx, cfg_line)
                insert_idx += 1
            mol_block = "\n".join(mol_lines)
        return mol_block

    def _start_calculation_worker(self, mol_block, options, run_id):
        """Initialize and start the background calculation worker."""
        try:
            thread = QThread()
            worker = CalculationWorker()
            worker.halt_ids = self.halt_ids
            worker.moveToThread(thread)

            worker.status_update.connect(self.update_status_bar)

            def _cleanup(t=thread, w=worker):
                t.quit()
                t.finished.connect(t.deleteLater)
                w.deleteLater()
                with contextlib.suppress(ValueError, AttributeError):
                    self._active_calc_threads.remove(t)

            def _on_finished(result):
                try: self.on_calculation_finished(result)
                finally: _cleanup()

            def _on_error(error_msg):
                try: self.on_calculation_error(error_msg)
                finally: _cleanup()

            worker.error.connect(_on_error)
            worker.finished.connect(_on_finished)
            thread.start()
            QTimer.singleShot(10, lambda: worker.start_work.emit(mol_block, options))
            self._active_calc_threads.append(thread)
        except (RuntimeError, TypeError, AttributeError) as e:
            logging.exception(f"Failed to start calculation worker: {e}")
            self.on_calculation_error(str(e))

    def trigger_conversion(self):
        self.last_successful_optimization_method = None
        self.constraints_3d = []

        if not self.data.atoms:
            self.plotter.clear()
            self.current_mol = None
            self.analysis_action.setEnabled(False)
            self.statusBar().showMessage("3D view cleared.")
            self.view_2d.setFocus()
            return

        # Reset states
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)

        mol = self.main_window_compute._prepare_rdkit_mol_for_conversion()
        if not mol:
            return

        num_frags = len(Chem.GetMolFrags(mol))
        if num_frags > 1:
            self.statusBar().showMessage(
                f"Converting {num_frags} molecules to 3D with collision detection..."
            )
        else:
            self.statusBar().showMessage("Calculating 3D structure...")

        mol_block = self.main_window_compute._setup_mol_block_for_worker(mol)
        
        run_id = int(getattr(self, "next_conversion_id", 1))
        self.next_conversion_id = run_id + 1
        if not hasattr(self, "active_worker_ids"):
            self.active_worker_ids = set()
        self.active_worker_ids.add(run_id)

        # UI Updates
        self.convert_button.setText("Halt conversion")
        self._safe_disconnect(self.convert_button.clicked)
        self.convert_button.clicked.connect(self.halt_conversion)
        self.cleanup_button.setEnabled(False)
        self._enable_3d_features(False)
        self.plotter.clear()

        # Add 'Calculating...' overlay
        bg_qcolor = QColor(self.settings.get("background_color", "#919191"))
        text_color = "black" if (bg_qcolor.isValid() and bg_qcolor.toHsl().lightness() > 128) else "white"
        self._calculating_text_actor = self.plotter.add_text(
            "Calculating...", position="lower_right", font_size=15, color=text_color, name="calculating_text"
        )
        self.plotter.render()

        # Build options
        conv_mode = self.__dict__.pop("_temp_conv_mode", self.settings.get("3d_conversion_mode", "fallback"))
        opt_method = self.__dict__.pop("_temp_optimization_method", None) or getattr(self, "optimization_method", "MMFF_RDKIT")
        
        options = {
            "conversion_mode": conv_mode,
            "optimization_method": opt_method,
            "optimize_intermolecular_interaction_rdkit": self.settings.get("optimize_intermolecular_interaction_rdkit", True),
            "worker_id": run_id
        }

        self.main_window_compute._start_calculation_worker(mol_block, options, run_id)
        self.push_undo_state()
        self.update_chiral_labels()
        self.view_2d.setFocus()

    def halt_conversion(self):  
        """Halt the in-progress conversion."""
        # Add active worker IDs to halt_ids
        wids_to_halt = set(getattr(self, "active_worker_ids", set()))
        if wids_to_halt:
            self.halt_ids.update(wids_to_halt)
        if hasattr(self, "active_worker_ids"):
            self.active_worker_ids.clear()

        # Restore UI
        self._restore_button_ui()
        self.cleanup_button.setEnabled(True)

        # Remove 'Calculating...' text
        self._remove_calculating_text()

        self.statusBar().showMessage("Halted")

    def check_chemistry_problems_fallback(self):
        """Fallback chemistry check when RDKit fails."""
        # Suppress non-critical errors during fallback chemistry check to avoid crashing the UI
        with contextlib.suppress(AttributeError, RuntimeError, TypeError, ValueError):
            self.scene.clear_all_problem_flags()
            problem_atoms = []

            for atom_id, atom_data in self.data.atoms.items():
                atom_item = atom_data.get("item")
                if atom_item:
                    symbol = atom_data["symbol"]
                    charge = atom_data.get("charge", 0)
                    bond_count = sum(b.get("order", 1) for k, b in self.data.bonds.items() if atom_id in k)

                    if is_problematic_valence(symbol, bond_count, charge):
                        problem_atoms.append(atom_item)

            if problem_atoms:
                for atom_item in problem_atoms:
                    atom_item.has_problem = True
                    atom_item.update()
                self.statusBar().showMessage(f"Error: {len(problem_atoms)} chemistry problem(s) found (valence issues).")
            else:
                self.statusBar().showMessage("Error: Invalid chemical structure (RDKit conversion failed).")

            self.scene.clearSelection()
            self.view_2d.setFocus()

    def optimize_3d_structure(self):
        """Optimize 3D structure using force fields."""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to optimize.")
            return

        # If a prior chemical/sanitization check was attempted and failed, do not run optimization
        if getattr(self, "chem_check_tried", False) and getattr(
            self, "chem_check_failed", False
        ):
            self.statusBar().showMessage(
                "3D optimization disabled: molecule failed chemical sanitization."
            )
            # Ensure the Optimize 3D button is disabled to reflect this
            if hasattr(self, "optimize_3d_button"):
                self.optimize_3d_button.setEnabled(False)
            return

        try:
            # Allow a temporary optimization method override (right-click menu)
            method = self.__dict__.pop("_temp_optimization_method", None) or getattr(
                self, "optimization_method", "MMFF_RDKIT"
            )
            method = method.upper() if method else "MMFF_RDKIT"

            # Check for conformer
            if self.current_mol.GetNumConformers() == 0:
                self.statusBar().showMessage(
                    "No conformer found: cannot optimize. Embed molecule first."
                )
                return

            # If it's a plugin method, we still run it synchronously for now as we don't know
            # if the plugin callback is thread-safe.
            if (
                hasattr(self, "plugin_manager")
                and hasattr(self.plugin_manager, "optimization_methods")
                and method in self.plugin_manager.optimization_methods
            ):
                self.statusBar().showMessage(f"Optimizing with plugin method: {method}...")
                QApplication.processEvents()
                info = self.plugin_manager.optimization_methods[method]
                callback = info["callback"]
                try:
                    success = callback(self.current_mol)
                    if not success:
                        self.statusBar().showMessage(
                            f"Optimization method '{method}' returned failure."
                        )
                    else:
                        self.draw_molecule_3d(self.current_mol)
                        self.statusBar().showMessage(f"Optimization ({method}) successful.")
                        self.push_undo_state()
                        self.view_2d.setFocus()
                    return
                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    self.statusBar().showMessage(
                        f"Plugin optimization error ({method}): {e}"
                    )
                    return

            # For RDKit and OBabel methods, use CalculationWorker to avoid UI freeze and allow Halt
            if "_RDKIT" in method or "_OBABEL" in method:
                self.statusBar().showMessage(f"Optimizing 3D structure ({method})...")
                
                # Use current_mol as the source MOL block for optimization
                mol_block = Chem.MolToMolBlock(self.current_mol, includeStereo=True)
                
                options = {
                    "conversion_mode": "optimize_only",
                    "optimization_method": method,
                    "optimize_intermolecular_interaction_rdkit": self.settings.get("optimize_intermolecular_interaction_rdkit", True)
                }
                
                # Assign unique run ID
                run_id = int(getattr(self, "next_conversion_id", 1))
                self.next_conversion_id = run_id + 1
                options["worker_id"] = run_id

                if hasattr(self, "active_worker_ids"):
                    self.active_worker_ids.add(run_id)
                else:
                    self.active_worker_ids = {run_id}

                # Change button to 'Halt'
                if hasattr(self, "optimize_3d_button"):
                    self.optimize_3d_button.setText("Halt optimize")
                    self._safe_disconnect(self.optimize_3d_button.clicked)
                    self.optimize_3d_button.clicked.connect(self.halt_conversion)

                # Disable features
                self._enable_3d_features(False)

                # Setup worker thread
                try:
                    thread = QThread()
                    worker = CalculationWorker()
                    worker.halt_ids = self.halt_ids
                    worker.moveToThread(thread)
                    
                    worker.status_update.connect(self.update_status_bar)
                    
                    def _on_opt_worker_finished(result, w=worker, t=thread):
                        # Re-enable optimize button UI
                        if hasattr(self, "optimize_3d_button"):
                            self._safe_disconnect(self.optimize_3d_button.clicked)
                            self.optimize_3d_button.setText("Optimize 3D")
                            self.optimize_3d_button.clicked.connect(self.optimize_3d_structure)
                        
                        self.on_calculation_finished(result)

                        t.quit()
                        t.finished.connect(t.deleteLater)
                        w.deleteLater()
                        with contextlib.suppress(ValueError, AttributeError):
                            self._active_calc_threads.remove(t)

                    def _on_opt_worker_error(error_msg, w=worker, t=thread):
                        # Re-enable optimize button UI
                        if hasattr(self, "optimize_3d_button"):
                            self._safe_disconnect(self.optimize_3d_button.clicked)
                            self.optimize_3d_button.setText("Optimize 3D")
                            self.optimize_3d_button.clicked.connect(self.optimize_3d_structure)

                        self.on_calculation_error(error_msg)

                        t.quit()
                        t.finished.connect(t.deleteLater)
                        w.deleteLater()
                        with contextlib.suppress(ValueError, AttributeError):
                            self._active_calc_threads.remove(t)

                    worker.error.connect(_on_opt_worker_error)
                    worker.finished.connect(_on_opt_worker_finished)
                    thread.start()
                    
                    # Start work
                    QTimer.singleShot(10, lambda w=worker, m=mol_block, o=options: w.start_work.emit(m, o))
                    
                    self._active_calc_threads.append(thread)
                except (RuntimeError, TypeError, AttributeError) as e:
                    # Catch specific worker/thread creation failures.
                    self.on_calculation_error(str(e))
                
                return
            else:
                self.statusBar().showMessage(
                    "Selected optimization method is not available. Use MMFF94 (RDKit) or UFF (RDKit)."
                )
                return
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            self.statusBar().showMessage(f"3D optimization error: {e}")
            self.view_2d.setFocus()

    def on_calculation_finished(self, result):
        # Handle result tuple or single mol
        worker_id = None
        mol = None
        if isinstance(result, tuple) and len(result) == 2:
            worker_id, mol = result
        else:
            mol = result

        # Discard stale results
        if worker_id is not None:
            active_ids = getattr(self, "active_worker_ids", set())
            if worker_id not in active_ids:
                # Cleanup 'Calculating...' UI for stale result
                self._remove_calculating_text()
                self._restore_button_ui()
                self.cleanup_button.setEnabled(True)
                self.statusBar().showMessage(
                    "Ignored result from stale conversion."
                )
                return

        # Cleanup worker IDs
        if worker_id is not None:
            getattr(self, "active_worker_ids", set()).discard(worker_id)
            getattr(self, "halt_ids", set()).discard(worker_id)

        self.dragged_atom_info = None
        self.current_mol = mol
        self.is_xyz_derived = False  # Not XYZ-derived

        # Record the optimization method used for this conversion if available.
        opt_method = None
        if mol is not None and hasattr(mol, "HasProp") and mol.HasProp("_pme_optimization_method"):
            opt_method = mol.GetProp("_pme_optimization_method")

        # Store the optimization method using its user-friendly UI label if available.
        if opt_method:
            method_key = str(opt_method).upper()
            label = getattr(self, "opt3d_method_labels", {}).get(method_key, opt_method)
            self.last_successful_optimization_method = label
        else:
            self.last_successful_optimization_method = None

        # Restore atom properties (lost during worker process)
        if hasattr(self, "original_atom_properties"):
            for i, original_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(i)
                    atom.SetIntProp("_original_atom_id", original_id)

        # Create atom ID mapping
        self.create_atom_id_mapping()

        # Set chiral centers from 2D stereo
        skip_chem_property = False
        if mol is not None and hasattr(mol, "HasProp") and mol.HasProp("_xyz_skip_checks"):
            # Suppress potential errors if the property is malformed 
            with contextlib.suppress(RuntimeError, TypeError, ValueError, KeyError):
                skip_chem_property = bool(mol.GetIntProp("_xyz_skip_checks"))

        if not skip_chem_property:
            try:
                if mol.GetNumConformers() > 0:
                    # Save 2D stereo as property before 3D calculation
                    for bond in mol.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            bond.SetIntProp("_original_2d_stereo", bond.GetStereo())

                    # Assign stereochemistry respecting 2D info
                    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)

                self.update_chiral_labels()
            except (AttributeError, RuntimeError, TypeError):
                # Suppress non-critical stereochemistry assignment errors for malformed fragments.
                pass
        self.draw_molecule_3d(mol)

        # Collision avoidance handled by worker
        try:
            frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
            if len(frags) > 1:
                self.statusBar().showMessage(
                    f"{len(frags)} molecules converted with collision avoidance (background)."
                )
        except (AttributeError, RuntimeError, TypeError):
            # Suppress errors if fragment count retrieval fails for an invalid molecule state.
            pass
        # Remove 'Calculating...' text and refresh
        self._remove_calculating_text()
        try:
            self.plotter.render()
        except (AttributeError, RuntimeError, TypeError):
            # Suppress transient plotter render errors if the 3D scene is not yet ready.
            pass
        if self.last_successful_optimization_method:
            self.statusBar().showMessage(f"3D calculation ({self.last_successful_optimization_method}) successful.")
        else:
            self.statusBar().showMessage("3D calculation successful.")

        # Restore button UI
        self._restore_button_ui()

        self.push_undo_state()
        self.view_2d.setFocus()
        self.cleanup_button.setEnabled(True)

        # Enable 3D features
        self._enable_3d_features(True)

        self.plotter.reset_camera()

        # Setup 3D hover
        self.setup_3d_hover()

        # Update menu items
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

    def create_atom_id_mapping(self):
        """Map 2D atom IDs to 3D RDKit indices."""
        if not self.current_mol:
            return

        self.atom_id_to_rdkit_idx_map = {}

        # Create mapping from RDKit properties
        for i in range(self.current_mol.GetNumAtoms()):
            rdkit_atom = self.current_mol.GetAtomWithIdx(i)
            try:
                original_atom_id = rdkit_atom.GetIntProp("_original_atom_id")
                self.atom_id_to_rdkit_idx_map[original_atom_id] = i
            except KeyError:
                # Skip if property missing
                continue

    def on_calculation_error(self, result):
        """Handle worker error or halt."""
        worker_id = None
        error_message = ""
        if isinstance(result, tuple) and len(result) == 2:
            worker_id, error_message = result
        else:
            error_message = str(result)

        # If this error is from a stale/previous worker (not in active set), ignore it.
        if worker_id is not None and worker_id not in getattr(
            self, "active_worker_ids", set()
        ):
            # Stale/late error from a previously-halted worker; ignore to avoid clobbering newer runs
            print(f"Ignored stale error from worker {worker_id}: {error_message}")
            return

        # Cleanup plotter and 'Calculating...' text
        try:
            if self.current_mol is None:
                self.plotter.clear()
            else:
                # Re-render if molecule exists
                self.plotter.render()
        except (AttributeError, RuntimeError, TypeError):
            # Suppress non-critical error
            pass
        # Remove 'Calculating...' text
        self._remove_calculating_text()

        self.dragged_atom_info = None

        # Cleanup worker ID
        if worker_id is not None:
            getattr(self, "active_worker_ids", set()).discard(worker_id)

        # If a halt message and there are no active workers left, the user
        # already saw the halt message — suppress duplicate noise.
        low = (error_message or "").lower()
        if "halt" in low and not getattr(self, "active_worker_ids", set()):
            return

        # Interactive Fallback to UFF if optimization explicitly failed
        if "Optimization with" in error_message and "failed" in error_message:
            try:
                self.statusBar().showMessage(f"Optimization failed: {error_message}")
                reply = QMessageBox.question(
                    self,
                    "Optimization Failed",
                    f"{error_message}\n\nRDKit's UFF algorithm is often more robust for unusual elements.\nWould you like to try optimizing with UFF (RDKit) instead?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.Yes,
                )
                if reply == QMessageBox.StandardButton.Yes:
                    self._temp_optimization_method = "UFF_RDKIT"
                    self.optimize_3d_structure()
                    # Return immediately so the rest of the error cleanup doesn't disable buttons permanently
                    return
            except (AttributeError, RuntimeError, TypeError):
                # Suppress non-critical error
                pass
        if error_message == "Halted":
            self.statusBar().showMessage("Halted")
        else:
            self.statusBar().showMessage(f"Error: {error_message}")

        # Restore button UI
        self.cleanup_button.setEnabled(True)
        self._restore_button_ui()

        # Enable 3D features if valid molecule exists
        if self.current_mol is None:
            if hasattr(self, "optimize_3d_button"):
                self.optimize_3d_button.setEnabled(False)
            if hasattr(self, "export_button"):
                self.export_button.setEnabled(False)
            # Disable 3D features if no molecule
            self._enable_3d_features(False)
            # Disable 3D edit actions
            self._enable_3d_edit_actions(False)
        else:
            # We HAVE a molecule, ensure features are enabled
            self._enable_3d_features(True)

        # Handle menu item states
        if self.current_mol is None:
            if hasattr(self, "analysis_action"):
                self.analysis_action.setEnabled(False)
            if hasattr(self, "edit_3d_action"):
                self.edit_3d_action.setEnabled(False)
        else:
            # Ensure they are enabled if we have a molecule
            if hasattr(self, "analysis_action"):
                self.analysis_action.setEnabled(True)
            if hasattr(self, "edit_3d_action"):
                self.edit_3d_action.setEnabled(True)

        # Refresh UI
        try:
            self.plotter.render()
        except (AttributeError, RuntimeError, TypeError):
            # Suppress non-critical error
            pass
        # Return focus to 2D editor
        self.view_2d.setFocus()
