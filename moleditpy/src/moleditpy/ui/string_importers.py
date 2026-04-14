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
class StringImporterManager:
    """Mixin for string-based molecular input (SMILES, InChI)."""

    def __init__(self, host):
        self.host = host

    def import_smiles_dialog(self) -> None:
        """Dialog for SMILES input."""
        smiles, ok = QInputDialog.getText(
            self.host, "Import SMILES", "Enter SMILES string:"
        )
        if ok and smiles:
            self.load_from_smiles(smiles)

    def import_inchi_dialog(self) -> None:
        """Dialog for InChI input."""
        inchi, ok = QInputDialog.getText(
            self.host, "Import InChI", "Enter InChI string:"
        )
        if ok and inchi:
            self.load_from_inchi(inchi)

    def _placement_center(self) -> QPointF:
        """Return the scene point where the next imported molecule should be centered.

        If the editor already contains atoms, offset to the right of the
        rightmost atom so the new fragment does not overlap.  Otherwise use
        the current viewport center.
        """
        existing = self.host.state_manager.data.atoms
        if existing:

            def _x(v):
                return v["pos"].x() if hasattr(v["pos"], "x") else v["pos"][0]

            def _y(v):
                return v["pos"].y() if hasattr(v["pos"], "y") else v["pos"][1]

            max_x = max(_x(v) for v in existing.values())
            avg_y = sum(_y(v) for v in existing.values()) / len(existing)
            return QPointF(max_x + 80.0, avg_y)

        return self.host.init_manager.view_2d.mapToScene(
            self.host.init_manager.view_2d.viewport().rect().center()
        )

    def _place_mol_bonds(self, mol, rdkit_idx_to_my_id: dict) -> None:
        """Create bonds in the 2D scene from an RDKit molecule."""
        for bond in mol.GetBonds():
            b_idx, e_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            b_type = bond.GetBondTypeAsDouble()
            b_dir = bond.GetBondDir()
            stereo = 0
            if b_dir == Chem.BondDir.BEGINWEDGE:
                stereo = 1
            elif b_dir == Chem.BondDir.BEGINDASH:
                stereo = 2
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                    stereo = 3
                elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                    stereo = 4

            if b_idx in rdkit_idx_to_my_id and e_idx in rdkit_idx_to_my_id:
                a1_id = rdkit_idx_to_my_id[b_idx]
                a2_id = rdkit_idx_to_my_id[e_idx]
                a1_item = self.host.state_manager.data.atoms[a1_id]["item"]
                a2_item = self.host.state_manager.data.atoms[a2_id]["item"]
                self.host.init_manager.scene.create_bond(
                    a1_item, a2_item, bond_order=int(b_type), bond_stereo=stereo
                )

    def _place_mol_atoms(self, mol, conf, place_center: QPointF) -> dict:
        """Create atoms in the 2D scene and return rdkit_idx → atom_id mapping."""
        SCALE_FACTOR = 50.0
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
            scene_x = ((pos.x - mol_center_x) * SCALE_FACTOR) + place_center.x()
            scene_y = (-(pos.y - mol_center_y) * SCALE_FACTOR) + place_center.y()
            atom_id = self.host.init_manager.scene.create_atom(
                atom.GetSymbol(),
                QPointF(scene_x, scene_y),
                charge=atom.GetFormalCharge(),
            )
            rdkit_idx_to_my_id[i] = atom_id
        return rdkit_idx_to_my_id

    def load_from_smiles(self, smiles_string: str) -> None:
        """Add a molecule from a SMILES string to the 2D editor."""
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
            self.host.statusBar().showMessage(f"Invalid SMILES: {e}")
            return
        except (RuntimeError, TypeError, AttributeError) as e:
            self.host.statusBar().showMessage(f"Error parsing SMILES: {e}")
            return

        try:
            place_center = self._placement_center()
            rdkit_idx_to_my_id = self._place_mol_atoms(
                mol, mol.GetConformer(), place_center
            )
            self._place_mol_bonds(mol, rdkit_idx_to_my_id)

            self.host.statusBar().showMessage("Successfully loaded from SMILES.")
            self.host.init_manager.scene.update_all_items()
            self.host.edit_actions_manager.push_undo_state()
            QTimer.singleShot(0, self.host.view_3d_manager.fit_to_view)

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            self.host.statusBar().showMessage(f"Error loading from SMILES: {e}")

    def load_from_inchi(self, inchi_string: str) -> None:
        """Add a molecule from an InChI string to the 2D editor."""
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
            self.host.statusBar().showMessage(f"Invalid InChI: {e}")
            return
        except (RuntimeError, TypeError, AttributeError) as e:
            self.host.statusBar().showMessage(f"Error parsing InChI: {e}")
            return

        try:
            place_center = self._placement_center()
            rdkit_idx_to_my_id = self._place_mol_atoms(
                mol, mol.GetConformer(), place_center
            )
            self._place_mol_bonds(mol, rdkit_idx_to_my_id)

            self.host.statusBar().showMessage("Successfully loaded from InChI.")
            self.host.init_manager.scene.update_all_items()
            self.host.edit_actions_manager.push_undo_state()
            QTimer.singleShot(0, self.host.view_3d_manager.fit_to_view)

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            self.host.statusBar().showMessage(f"Error loading from InChI: {e}")


# Backward-compat aliases
MainWindowStringImporters = StringImporterManager
