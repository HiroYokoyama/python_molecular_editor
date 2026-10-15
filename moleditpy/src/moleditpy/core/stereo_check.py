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


@dataclass(frozen=True)
class WedgeForLabel:
    """How to draw a bond so that the stereocenter at its narrow end reads as R/S."""

    key: Tuple[int, int]  # (stereocenter, other end): data.bonds key, narrow end first
    stereo: int  # 1 = wedge, 2 = hash
    cleared: Tuple[Tuple[int, int], ...]  # other wedges at the center to make plain


def wedge_for_label(
    data: MolecularData,
    bond_key: Tuple[int, int],
    label: str,
    centers: Optional[List[int]] = None,
) -> Optional[WedgeForLabel]:
    """Wedge or hash on *bond_key* that makes one of its ends CIP *label*.

    *centers* is the order in which the bond's ends are tried (default: the
    key order). Each end is tried first with the other wedges kept; if they
    make *label* unreachable, then with this center's other wedges and hashes
    made plain, so the requested label wins. Returns None when neither end
    can be *label* (not a stereocenter). *data* is left unchanged.
    """
    if bond_key not in data.bonds or label not in ("R", "S"):
        return None
    ends = list(centers) if centers else list(bond_key)
    original = list(data.bonds.items())
    template = dict(data.bonds[bond_key], order=1)
    trial_bonds = {k: dict(v) for k, v in original if k != bond_key}
    try:
        for center in ends:
            if center not in bond_key:
                continue
            other = bond_key[1] if center == bond_key[0] else bond_key[0]
            own = tuple(
                k
                for k, v in trial_bonds.items()
                if k[0] == center and v.get("stereo") in (1, 2)
            )
            for cleared in ((), own) if own else ((),):
                for stereo in (1, 2):
                    data.bonds.clear()
                    for k, v in trial_bonds.items():
                        data.bonds[k] = dict(v, stereo=0) if k in cleared else dict(v)
                    data.bonds[(center, other)] = dict(template, stereo=stereo)
                    if drawn_chirality(data).get(center) == label:
                        return WedgeForLabel((center, other), stereo, cleared)
        return None
    except (RuntimeError, ValueError, TypeError) as e:
        logging.warning("R/S wedge search failed: %s", e)
        return None
    finally:
        data.bonds.clear()
        data.bonds.update(original)


def bond_label(data: MolecularData, bond_key: Tuple[int, int]) -> Optional[str]:
    """CIP label the wedge/hash *bond_key* gives the center at its narrow end."""
    bond = data.bonds.get(bond_key)
    if bond is None or bond.get("stereo") not in (1, 2):
        return None
    return drawn_chirality(data).get(bond_key[0])


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


def ambiguous_wedge_centers(data: MolecularData) -> List[int]:
    """Atom IDs of stereocenters whose wedges/hashes specify no configuration.

    These are drawings RDKit rejects as ambiguous (for example a wedge and a
    hash on opposite bonds), so the 3D structure gets an arbitrary
    configuration there. Atoms that are not stereocenters are not reported.
    """
    if not any(bond.get("stereo") in (1, 2) for bond in data.bonds.values()):
        return []
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    if mol is None:
        return []
    probe = Chem.Mol(mol)
    Chem.AssignChiralTypesFromBondDirs(probe)
    centers = {
        bond.GetBeginAtomIdx()
        for bond in mol.GetBonds()
        if bond.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    }
    found: List[int] = []
    for idx in sorted(centers):
        atom = probe.GetAtomWithIdx(idx)
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        # A stereocenter keeps an arbitrary tag through perception.
        trial = Chem.Mol(mol)
        trial.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        Chem.AssignStereochemistry(trial, cleanIt=True, force=True)
        trial_atom = trial.GetAtomWithIdx(idx)
        if trial_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and (
            trial_atom.HasProp("_original_atom_id")
        ):
            found.append(trial_atom.GetIntProp("_original_atom_id"))
    return found


def find_ambiguous_centers(
    data: MolecularData, mol_3d: Chem.Mol
) -> List[ChiralityMismatch]:
    """Ambiguously wedged stereocenters, with ``drawn`` set to "ambiguous".

    ``actual`` is the configuration the 3D structure happened to get.
    """
    try:
        atom_ids = ambiguous_wedge_centers(data)
        if not atom_ids:
            return []
        actual = actual_chirality(mol_3d)
    except (RuntimeError, ValueError) as e:
        logging.warning("Ambiguous wedge check failed: %s", e)
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
            drawn="ambiguous",
            actual=actual.get(atom_id),
        )
        for atom_id in atom_ids
        if atom_id in data.atoms
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
