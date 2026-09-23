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
from typing import Dict, FrozenSet, List, Optional, Tuple

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


@dataclass(frozen=True)
class EZMismatch:
    """A labelled double bond whose 3D configuration differs from the label."""

    atom_ids: Tuple[int, int]
    symbols: Tuple[str, str]
    rdkit_bond_index: Optional[int]
    rdkit_atom_indices: Optional[Tuple[int, int]]
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
    # Nothing is specified without a wedge or hash; skip CIP labelling.
    if not any(bond.get("stereo") in (1, 2) for bond in data.bonds.values()):
        return {}
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


def _bond_cip_by_atom_ids(
    mol: Chem.Mol,
) -> Dict[FrozenSet[int], Tuple[int, str]]:
    """Map editor atom-ID pairs to (bond index, CIP E/Z) for every labelled bond.

    Expects CIP labels already assigned on *mol*.
    """
    out: Dict[FrozenSet[int], Tuple[int, str]] = {}
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        if not (
            bond.HasProp("_CIPCode")
            and begin.HasProp("_original_atom_id")
            and end.HasProp("_original_atom_id")
        ):
            continue
        key = frozenset(
            (begin.GetIntProp("_original_atom_id"), end.GetIntProp("_original_atom_id"))
        )
        out[key] = (bond.GetIdx(), bond.GetProp("_CIPCode"))
    return out


def drawn_ez(data: MolecularData) -> Dict[Tuple[int, int], str]:
    """E/Z labels ("E"/"Z") of the drawn double bonds that have a CIP E/Z.

    A label on a bond with two identical groups at one end is dropped: that
    bond has no E/Z, so the label cannot be checked.
    """
    labelled = {
        (id1, id2): ("Z" if bond["stereo"] == 3 else "E")
        for (id1, id2), bond in data.bonds.items()
        if bond.get("stereo") in (3, 4)
    }
    if not labelled:
        return {}
    drawn_mol = data.to_rdkit_mol(use_2d_stereo=False)
    if drawn_mol is None:
        return {}
    rdCIPLabeler.AssignCIPLabels(drawn_mol)
    comparable = _bond_cip_by_atom_ids(drawn_mol)
    return {
        pair: label for pair, label in labelled.items() if frozenset(pair) in comparable
    }


def actual_ez(mol_3d: Chem.Mol) -> Dict[int, str]:
    """CIP E/Z ("E"/"Z") of every stereo double bond in a 3D structure.

    Keyed by RDKit bond index and read from the coordinates, with the same
    labeller as find_ez_mismatches, so a label and the check always agree.
    """
    if mol_3d.GetNumConformers() == 0:
        return {}
    probe = Chem.Mol(mol_3d)
    Chem.AssignStereochemistryFrom3D(probe)
    rdCIPLabeler.AssignCIPLabels(probe)
    return {
        bond.GetIdx(): bond.GetProp("_CIPCode")
        for bond in probe.GetBonds()
        if bond.HasProp("_CIPCode") and bond.GetProp("_CIPCode") in ("E", "Z")
    }


def find_ez_mismatches(data: MolecularData, mol_3d: Chem.Mol) -> List[EZMismatch]:
    """Double bonds labelled E or Z in 2D whose 3D configuration differs.

    Only bonds with a CIP E/Z are compared (see drawn_ez).
    """
    if mol_3d.GetNumConformers() == 0:
        return []
    try:
        drawn = drawn_ez(data)
        if not drawn:
            return []
        probe = Chem.Mol(mol_3d)
        Chem.AssignStereochemistryFrom3D(probe)
        rdCIPLabeler.AssignCIPLabels(probe)
        actual = _bond_cip_by_atom_ids(probe)
    except (RuntimeError, ValueError) as e:
        logging.warning("E/Z check failed: %s", e)
        return []

    bond_of: Dict[FrozenSet[int], Chem.Bond] = {
        frozenset(
            (
                b.GetBeginAtom().GetIntProp("_original_atom_id"),
                b.GetEndAtom().GetIntProp("_original_atom_id"),
            )
        ): b
        for b in mol_3d.GetBonds()
        if b.GetBeginAtom().HasProp("_original_atom_id")
        and b.GetEndAtom().HasProp("_original_atom_id")
    }
    mismatches: List[EZMismatch] = []
    for (id1, id2), label in sorted(drawn.items()):
        if id1 not in data.atoms or id2 not in data.atoms:
            continue
        key = frozenset((id1, id2))
        got = actual.get(key)
        got_label = got[1] if got else None
        if got_label == label:
            continue
        bond = bond_of.get(key)
        mismatches.append(
            EZMismatch(
                atom_ids=(id1, id2),
                symbols=(
                    str(data.atoms[id1]["symbol"]),
                    str(data.atoms[id2]["symbol"]),
                ),
                rdkit_bond_index=bond.GetIdx() if bond is not None else None,
                rdkit_atom_indices=(
                    (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    if bond is not None
                    else None
                ),
                drawn=label,
                actual=got_label,
            )
        )
    return mismatches
