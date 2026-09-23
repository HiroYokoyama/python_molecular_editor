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
from dataclasses import dataclass
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from .molecular_data import MolecularData


@dataclass(frozen=True)
class ChiralityMismatch:
    """A stereocenter whose 3D configuration differs from the 2D drawing."""

    atom_id: int
    symbol: str
    rdkit_index: Optional[int]
    drawn: str
    actual: Optional[str]


def _cip_labels_by_atom_id(mol: Chem.Mol) -> Dict[int, str]:
    """Map editor atom IDs to CIP labels ('R'/'S') for every labelled center."""
    rdCIPLabeler.AssignCIPLabels(mol)
    labels: Dict[int, str] = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode") and atom.HasProp("_original_atom_id"):
            labels[atom.GetIntProp("_original_atom_id")] = atom.GetProp("_CIPCode")
    return labels


def drawn_chirality(data: MolecularData) -> Dict[int, str]:
    """CIP labels of the stereocenters the 2D drawing specifies with wedges."""
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    if mol is None:
        return {}
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return _cip_labels_by_atom_id(mol)


def actual_chirality(mol_3d: Chem.Mol) -> Dict[int, str]:
    """CIP labels of the stereocenters in a 3D structure, read from coordinates."""
    if mol_3d.GetNumConformers() == 0:
        return {}
    mol = Chem.Mol(mol_3d)
    Chem.AssignStereochemistryFrom3D(mol)
    return _cip_labels_by_atom_id(mol)


def find_chirality_mismatches(
    data: MolecularData, mol_3d: Chem.Mol
) -> List[ChiralityMismatch]:
    """Stereocenters drawn in 2D whose configuration the 3D structure does not keep.

    Only centers the drawing specifies (wedge or hash) are compared; an
    unspecified center can come out either way. A center that is no longer a
    stereocenter in 3D is reported with ``actual`` None.
    """
    try:
        drawn = drawn_chirality(data)
        if not drawn:
            return []
        actual = actual_chirality(mol_3d)
    except (RuntimeError, ValueError) as e:
        logging.warning("Chirality check failed: %s", e)
        return []

    index_of = {
        atom.GetIntProp("_original_atom_id"): atom.GetIdx()
        for atom in mol_3d.GetAtoms()
        if atom.HasProp("_original_atom_id")
    }
    return [
        ChiralityMismatch(
            atom_id=atom_id,
            symbol=str(data.atoms[atom_id]["symbol"]),
            rdkit_index=index_of.get(atom_id),
            drawn=label,
            actual=actual.get(atom_id),
        )
        for atom_id, label in sorted(drawn.items())
        if actual.get(atom_id) != label and atom_id in data.atoms
    ]
