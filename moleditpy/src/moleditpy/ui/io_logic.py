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
import os
import json
import pickle
from typing import Any, Optional

from PyQt6.QtWidgets import QFileDialog, QMessageBox

class IOManager:
    """Independent manager for IO actions (Load/Save), ported from Mixins."""

    def __init__(self, host: Any) -> None:
        self.host = host

    def fix_mol_counts_line(self, line: str) -> str:
        if "V3000" in line or "V2000" in line: return line
        prefix = line.rstrip().ljust(33)[0:33]
        return prefix + " V2000"

    def fix_mol_block(self, mol_block: str) -> str:
        lines = mol_block.splitlines()
        if len(lines) < 4: return mol_block
        lines[3] = self.fix_mol_counts_line(lines[3])
        return "\n".join(lines)

    def load_xyz_file(self, file_path: str) -> Optional[Any]:
        # Implementation from molecular_parsers.py
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, rdGeometry
            with open(file_path, "r", encoding="utf-8") as f:
                raw_lines = f.readlines()
            lines = [ln.strip() for ln in raw_lines if not ln.strip().startswith("#")]
            while lines and not lines[0]: lines.pop(0)
            if len(lines) < 2: return None
            num_atoms = int(lines[0])
            atoms_data = []
            for i, line in enumerate(lines[2 : 2 + num_atoms]):
                parts = line.split()
                if len(parts) < 4: continue
                symbol = parts[0].capitalize()
                atoms_data.append((symbol, float(parts[1]), float(parts[2]), float(parts[3])))
            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms_data))
            for i, (symbol, x, y, z) in enumerate(atoms_data):
                atom = Chem.Atom(symbol)
                mol.AddAtom(atom)
                conf.SetAtomPosition(i, rdGeometry.Point3D(x, y, z))
            mol.AddConformer(conf)
            # Simplified bond estimation for XYV Load
            from rdkit.Chem import rdDetermineBonds
            try: rdDetermineBonds.DetermineBonds(mol, charge=0)
            except: pass
            return mol.GetMol()
        except: return None

    def save_project(self) -> None:
        """Save (Ctrl+S) - Defaults to PMEPRJ format."""
        if not self.host.data.atoms and not self.host.current_mol:
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return
            
        native_exts = [".pmeprj", ".pmeraw"]
        if self.host.current_file_path and any(
            self.host.current_file_path.lower().endswith(ext) for ext in native_exts
        ):
            try:
                if self.host.current_file_path.lower().endswith(".pmeraw"):
                    save_data = self.host.state_manager.get_current_state()
                    with open(self.host.current_file_path, "wb") as f:
                        pickle.dump(save_data, f)
                else:
                    json_data = self.host.state_manager.create_json_data()
                    with open(self.host.current_file_path, "w", encoding="utf-8") as f:
                        json.dump(json_data, f, indent=2, ensure_ascii=False)

                self.host.has_unsaved_changes = False
                self.host.state_manager.update_window_title()
                self.host.statusBar().showMessage(f"Project saved to {self.host.current_file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving project: {e}")
        else:
            self.save_project_as()

    def save_project_as(self) -> None:
        """Save As (Ctrl+Shift+S)"""
        if not self.host.data.atoms and not self.host.current_mol:
            self.host.statusBar().showMessage("Error: Nothing to save.")
            return

        default_dir = os.path.dirname(self.host.current_file_path) if self.host.current_file_path else ""
        file_path, _ = QFileDialog.getSaveFileName(
            self.host, "Save Project As", default_dir, "PME Project (*.pmeprj);;PME Raw (*.pmeraw);;All Files (*)"
        )

        if file_path:
            self.host.current_file_path = file_path
            self.save_project()

    def open_project(self) -> None:
        """Open an existing project file."""
        if not self.host.state_manager.check_unsaved_changes():
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self.host, "Open Project", "", "Project Files (*.pmeprj *.pmeraw);;All Files (*)"
        )
        if not file_path:
            return

        try:
            if file_path.lower().endswith(".pmeraw"):
                with open(file_path, "rb") as f:
                    state_data = pickle.load(f)
                self.host.state_manager.set_state_from_data(state_data)
            else:
                with open(file_path, "r", encoding="utf-8") as f:
                    json_data = json.load(f)
                self.host.state_manager.load_from_json_data(json_data)

            self.host.current_file_path = file_path
            self.host.has_unsaved_changes = False
            self.host.state_manager.update_window_title()
            self.host.edit_actions_manager.reset_undo_stack()
            self.host.statusBar().showMessage(f"Loaded {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.critical(self.host, "Load Error", f"Could not load project: {e}")

    def save_raw_data(self) -> None:
        """Export as PME Raw pickle format."""
        file_path, _ = QFileDialog.getSaveFileName(
            self.host, "Export Raw Data", "", "PME Raw (*.pmeraw);;All Files (*)"
        )
        if file_path:
            try:
                state = self.host.state_manager.get_current_state()
                with open(file_path, "wb") as f:
                    pickle.dump(state, f)
                self.host.statusBar().showMessage(f"Raw data exported to {file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Export error: {e}")

    def load_mol_file(self, file_path: Optional[str] = None) -> None:
        """Regular 2D MOL file loading logic."""
        if not self.host.state_manager.check_unsaved_changes():
            return

        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self.host, "Import MOL File", "", "Chemical Files (*.mol *.sdf);;All Files (*)"
            )
            if not file_path:
                return

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from PyQt6.QtCore import QPointF, QTimer

            if file_path.lower().endswith(".mol"):
                with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
            else:
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)

            if mol is None:
                raise ValueError("Failed to read molecule from file.")

            Chem.Kekulize(mol)
            self.host.ui_manager.restore_ui_for_editing()
            self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)
            self.host.current_mol = None
            self.host.plotter.clear()
            self.host.analysis_action.setEnabled(False)

            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)

            SCALE_FACTOR = 50.0
            view_center = self.host.view_2d.mapToScene(self.host.view_2d.viewport().rect().center())
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = sum(p.x for p in positions) / len(positions) if positions else 0.0
            mol_center_y = sum(p.y for p in positions) / len(positions) if positions else 0.0

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                scene_x = ((pos.x - mol_center_x) * SCALE_FACTOR) + view_center.x()
                scene_y = (-(pos.y - mol_center_y) * SCALE_FACTOR) + view_center.y()
                atom_id = self.host.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y), charge=atom.GetFormalCharge())
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                stereo = 0
                b_dir = bond.GetBondDir()
                if b_dir == Chem.BondDir.BEGINWEDGE: stereo = 1
                elif b_dir == Chem.BondDir.BEGINDASH: stereo = 2
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ: stereo = 3
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE: stereo = 4
                self.host.scene.create_bond(self.host.data.atoms[rdkit_idx_to_my_id[bond.GetBeginAtomIdx()]]["item"],
                                          self.host.data.atoms[rdkit_idx_to_my_id[bond.GetEndAtomIdx()]]["item"],
                                          bond_order=int(bond.GetBondTypeAsDouble()), bond_stereo=stereo)

            self.host.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.host.edit_actions_manager.reset_undo_stack()
            self.host.current_file_path = file_path
            self.host.has_unsaved_changes = False
            self.host.state_manager.update_window_title()
            self.host.scene.update_all_items()
            QTimer.singleShot(0, self.host.view_3d_manager.fit_to_view)
        except Exception as e:
            self.host.statusBar().showMessage(f"Error loading file: {e}")

    def save_as_mol(self) -> None:
        """Save current 2D structure as MOL file."""
        try:
            from moleditpy.utils.constants import VERSION
            mol_block = self.host.data.to_mol_block()
            if not mol_block:
                self.host.statusBar().showMessage("Error: No 2D data to save.")
                return
            lines = mol_block.split("\n")
            if len(lines) > 1 and "RDKit" in lines[1]:
                lines[1] = f"  MoleditPy Ver. {VERSION}  2D"
            file_path, _ = QFileDialog.getSaveFileName(self.host, "Save 2D MOL File", "", "MOL Files (*.mol);;All Files (*)")
            if file_path:
                if not file_path.lower().endswith(".mol"): file_path += ".mol"
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
                self.host.statusBar().showMessage(f"2D data saved to {file_path}")
        except Exception as e:
            self.host.statusBar().showMessage(f"Error saving MOL: {e}")

    def load_xyz_for_3d_viewing(self, file_path: Optional[str] = None) -> None:
        """Load XYZ file and display in 3D viewer."""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self.host, "Load 3D XYZ (View Only)", "", "XYZ Files (*.xyz);;All Files (*)")
            if not file_path: return

        try:
            mol = self.load_xyz_file(file_path)
            if mol is None: raise ValueError("Failed to create molecule from XYZ file.")
            
            self.host.edit_actions_manager.clear_2d_editor(push_to_undo=False)
            self.host.current_mol = mol
            
            if hasattr(self.host, "draw_molecule_3d"):
                self.host.view_3d_manager.draw_molecule_3d(self.host.current_mol)
            if hasattr(self.host, "_enter_3d_viewer_ui_mode"):
                self.host.ui_manager._enter_3d_viewer_ui_mode()
            self.host.statusBar().showMessage(f"3D Viewer Mode: Loaded {os.path.basename(file_path)}")
            self.host.current_file_path = file_path
            self.host.state_manager.update_window_title()
        except Exception as e:
            self.host.statusBar().showMessage(f"XYZ Load failed: {e}")

    def load_mol_file_for_3d_viewing(self, file_path: Optional[str] = None) -> None:
        """Open MOL/SDF file in 3D viewer."""
        if not self.host.state_manager.check_unsaved_changes(): return
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self.host, "Open MOL/SDF File", "", "MOL/SDF Files (*.mol *.sdf);;All Files (*)")
            if not file_path: return
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            if file_path.lower().endswith(".sdf"):
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                mol = next(suppl, None)
            else:
                with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                fixed_block = self.fix_mol_block(raw)
                mol = Chem.MolFromMolBlock(fixed_block, sanitize=True, removeHs=False)
            
            if mol is None: raise ValueError("Failed to load molecule.")
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol)
            
            self.host.current_mol = mol
            self.host.view_3d_manager.draw_molecule_3d(mol)
            self.host.ui_manager._enter_3d_viewer_ui_mode()
            self.host.statusBar().showMessage(f"Loaded {file_path} in 3D viewer")
            self.host.current_file_path = file_path
            self.host.state_manager.update_window_title()
        except Exception as e:
            self.host.statusBar().showMessage(f"3D MOL Load failed: {e}")

    def save_3d_as_mol(self) -> None:
        """Save current 3D structure as MOL file."""
        if not self.host.current_mol:
            self.host.statusBar().showMessage("Error: No 3D structure to save.")
            return
        from moleditpy.utils.constants import VERSION
        from rdkit import Chem
        file_path, _ = QFileDialog.getSaveFileName(self.host, "Save 3D MOL File", "", "MOL Files (*.mol);;All Files (*)")
        if file_path:
            try:
                mol_block = Chem.MolToMolBlock(self.host.current_mol, includeStereo=True)
                lines = mol_block.split("\n")
                if len(lines) > 1 and "RDKit" in lines[1]:
                    lines[1] = f"  MoleditPy Ver. {VERSION}  3D"
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
                self.host.statusBar().showMessage(f"3D data saved to {file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving 3D MOL: {e}")

    def save_as_xyz(self) -> None:
        """Save current 3D structure as XYZ file."""
        if not self.host.current_mol:
            self.host.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(self.host, "Save 3D XYZ File", "", "XYZ Files (*.xyz);;All Files (*)")
        if file_path:
            if not file_path.lower().endswith(".xyz"):
                file_path += ".xyz"
            try:
                from moleditpy.utils.constants import VERSION
                from rdkit.Chem import Descriptors, Chem
                conf = self.host.current_mol.GetConformer()
                num_atoms = self.host.current_mol.GetNumAtoms()
                xyz_lines = [str(num_atoms)]
                
                charge = 0
                if self.host.current_mol.HasProp("_xyz_charge"):
                    charge = self.host.current_mol.GetIntProp("_xyz_charge")
                else:
                    try: charge = Chem.GetFormalCharge(self.host.current_mol)
                    except: pass
                
                multiplicity = 1
                try: multiplicity = Descriptors.NumRadicalElectrons(self.host.current_mol) + 1
                except: pass

                xyz_lines.append(f"chrg = {charge}  mult = {multiplicity} | Generated by MoleditPy Ver. {VERSION}")
                for i in range(num_atoms):
                    pos = conf.GetAtomPosition(i)
                    symbol = self.host.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

                with open(file_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(xyz_lines) + "\n")
                self.host.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e:
                self.host.statusBar().showMessage(f"Error saving XYZ: {e}")

# Backward-compat aliases
MainWindowProjectIo = IOManager
MainWindowMolecularParsers = IOManager
MainWindowViewLoaders = IOManager
