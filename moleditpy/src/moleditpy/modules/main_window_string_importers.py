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
main_window_string_importers.py
Mixin class separated from main_window.py
"""
# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem

# PyQt6 Modules
from PyQt6.QtCore import QPointF, QTimer
from PyQt6.QtWidgets import QInputDialog

try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None


# --- Classes ---
class MainWindowStringImporters:
    """Mixin class separated from main_window.py."""

    def import_smiles_dialog(self):
        """Dialog for SMILES input."""
        smiles, ok = QInputDialog.getText(self, "Import SMILES", "Enter SMILES string:")
        if ok and smiles:
            self.load_from_smiles(smiles)

    def import_inchi_dialog(self):
        """Dialog for InChI input."""
        inchi, ok = QInputDialog.getText(self, "Import InChI", "Enter InChI string:")
        if ok and inchi:
            self.load_from_inchi(inchi)

    def load_from_smiles(self, smiles_string):
        """Load molecule from SMILES string to 2D editor."""
        if not self.check_unsaved_changes():
            return  # User cancelled

        cleaned_smiles = smiles_string.strip()

        try:
            mol = Chem.MolFromSmiles(cleaned_smiles)
            if mol is None:
                if not cleaned_smiles:
                    raise ValueError("SMILES string was empty.")
                raise ValueError("Invalid SMILES string.")

            AllChem.Compute2DCoords(mol)
            Chem.Kekulize(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid SMILES: {e}")
            return
        except (RuntimeError, TypeError, AttributeError) as e:
            self.statusBar().showMessage(f"Error parsing SMILES: {e}")
            return


        try:
            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)

            conf = mol.GetConformer()
            SCALE_FACTOR = 50.0

            view_center = self.view_2d.mapToScene(
                self.view_2d.viewport().rect().center()
            )
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = (
                sum(p.x for p in positions) / len(positions) if positions else 0.0
            )
            mol_center_y = (
                sum(p.y for p in positions) / len(positions) if positions else 0.0
            )

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()

                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y

                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()

                atom_id = self.scene.create_atom(
                    atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge
                )
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                # Single bond stereo
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1  # Wedge
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2  # Dash
                # Double bond E/Z
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3  # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4  # E

                if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                    a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                    a1_item = self.data.atoms[a1_id]["item"]
                    a2_item = self.data.atoms[a2_id]["item"]
                    self.scene.create_bond(
                        a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo
                    )

            self.statusBar().showMessage("Successfully loaded from SMILES.")
            self.scene.update_all_items()
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            self.statusBar().showMessage(f"Error loading from SMILES: {e}")



    def load_from_inchi(self, inchi_string):
        """Load molecule from InChI string to 2D editor."""
        if not self.check_unsaved_changes():
            return  # User cancelled

        cleaned_inchi = inchi_string.strip()

        try:
            mol = Chem.MolFromInchi(cleaned_inchi)
            if mol is None:
                if not cleaned_inchi:
                    raise ValueError("InChI string was empty.")
                raise ValueError("Invalid InChI string.")

            AllChem.Compute2DCoords(mol)
            Chem.Kekulize(mol)

            AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)
            conf = mol.GetConformer()
            AllChem.WedgeMolBonds(mol, conf)
        except ValueError as e:
            self.statusBar().showMessage(f"Invalid InChI: {e}")
            return
        except (RuntimeError, TypeError, AttributeError) as e:
            self.statusBar().showMessage(f"Error parsing InChI: {e}")
            return


        try:
            self.restore_ui_for_editing()
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            self.analysis_action.setEnabled(False)

            conf = mol.GetConformer()
            SCALE_FACTOR = 50.0

            view_center = self.view_2d.mapToScene(
                self.view_2d.viewport().rect().center()
            )
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            mol_center_x = (
                sum(p.x for p in positions) / len(positions) if positions else 0.0
            )
            mol_center_y = (
                sum(p.y for p in positions) / len(positions) if positions else 0.0
            )

            rdkit_idx_to_my_id = {}
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                charge = atom.GetFormalCharge()

                relative_x = pos.x - mol_center_x
                relative_y = pos.y - mol_center_y

                scene_x = (relative_x * SCALE_FACTOR) + view_center.x()
                scene_y = (-relative_y * SCALE_FACTOR) + view_center.y()

                atom_id = self.scene.create_atom(
                    atom.GetSymbol(), QPointF(scene_x, scene_y), charge=charge
                )
                rdkit_idx_to_my_id[i] = atom_id

            for bond in mol.GetBonds():
                b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                b_type = bond.GetBondTypeAsDouble()
                b_dir = bond.GetBondDir()
                stereo = 0
                # Single bond stereo
                if b_dir == Chem.BondDir.BEGINWEDGE:
                    stereo = 1  # Wedge
                elif b_dir == Chem.BondDir.BEGINDASH:
                    stereo = 2  # Dash
                # Double bond E/Z
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        stereo = 3  # Z
                    elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                        stereo = 4  # E

                if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                    a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                    a1_item = self.data.atoms[a1_id]["item"]
                    a2_item = self.data.atoms[a2_id]["item"]
                    self.scene.create_bond(
                        a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo
                    )

            self.statusBar().showMessage("Successfully loaded from InChI.")
            self.scene.update_all_items()
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.update_window_title()
            QTimer.singleShot(0, self.fit_to_view)

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            self.statusBar().showMessage(f"Error loading from InChI: {e}")


