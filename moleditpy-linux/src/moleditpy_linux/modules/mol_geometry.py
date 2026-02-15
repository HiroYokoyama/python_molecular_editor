#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""Pure-logic molecular-geometry helpers.

This module is intentionally free of any GUI imports so that it can
be unit-tested in isolation and shared across dialogs and main-window
submodules without circular dependencies.
"""

import numpy as np

# ------------------------------------------------------------------
# Graph traversal
# ------------------------------------------------------------------


def get_connected_group(mol, start_atom, exclude=None):
    """Return the set of atom indices reachable from *start_atom*
    without passing through *exclude*.

    Uses breadth-first search on the RDKit ``mol`` bond graph.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The molecule whose connectivity is queried.
    start_atom : int
        Index of the atom to start the traversal from.
    exclude : int or None
        If given, this atom is treated as a wall — the traversal
        will not enter or pass through it.

    Returns
    -------
    set[int]
        Indices of all atoms in the connected component containing
        *start_atom* (excluding *exclude* if it was encountered).
    """
    visited = set()
    to_visit = [start_atom]

    while to_visit:
        current = to_visit.pop()
        if current in visited or current == exclude:
            continue

        visited.add(current)

        atom = mol.GetAtomWithIdx(current)
        for bond in atom.GetBonds():
            other_idx = bond.GetOtherAtomIdx(current)
            if other_idx not in visited and other_idx != exclude:
                to_visit.append(other_idx)

    return visited


# ------------------------------------------------------------------
# Dihedral (torsion) angle
# ------------------------------------------------------------------


def calculate_dihedral(positions, i1, i2, i3, i4):
    """Compute the dihedral angle defined by four atom indices.

    Parameters
    ----------
    positions : array-like, shape (N, 3)
        Atom coordinates indexed by atom index.  Accepts any object
        that supports ``positions[i]`` → 3-element sequence
        (numpy array, list of lists, RDKit conformer positions, …).
    i1, i2, i3, i4 : int
        Atom indices.  The torsion is measured around the *i2–i3*
        central bond.

    Returns
    -------
    float
        Dihedral angle in **degrees**, in the range (−180, 180].
    """
    pos1 = np.asarray(positions[i1], dtype=float)
    pos2 = np.asarray(positions[i2], dtype=float)
    pos3 = np.asarray(positions[i3], dtype=float)
    pos4 = np.asarray(positions[i4], dtype=float)

    # Vectors between consecutive atoms
    v1 = pos2 - pos1  # 1->2
    v2 = pos3 - pos2  # 2->3 (central bond)
    v3 = pos4 - pos3  # 3->4

    # Normalize the central bond vector
    v2_norm = v2 / np.linalg.norm(v2)

    # Plane normal vectors
    n1 = np.cross(v1, v2)  # Normal to plane 1-2-3
    n2 = np.cross(v2, v3)  # Normal to plane 2-3-4

    n1_len = np.linalg.norm(n1)
    n2_len = np.linalg.norm(n2)

    if n1_len == 0 or n2_len == 0:
        return 0.0  # Atoms are collinear

    n1 = n1 / n1_len
    n2 = n2 / n2_len

    cos_angle = np.clip(np.dot(n1, n2), -1.0, 1.0)
    sin_angle = np.dot(np.cross(n1, n2), v2_norm)

    angle_rad = np.arctan2(sin_angle, cos_angle)
    return float(np.degrees(angle_rad))


# ------------------------------------------------------------------
# Valence sanity check
# ------------------------------------------------------------------

# Maximum *total bond order* before the atom is flagged, keyed by
# (symbol, requires_neutral).  Entries with ``requires_neutral=True``
# are only checked when the formal charge is zero.
_VALENCE_LIMITS = {
    "C": (4, False),
    "N": (3, True),
    "O": (2, True),
    "H": (1, False),
    "F": (1, True),
    "Cl": (1, True),
    "Br": (1, True),
    "I": (1, True),
}


def is_problematic_valence(symbol, bond_count, charge=0):
    """Return ``True`` if the atom's total bond order exceeds its
    typical maximum valence.

    Parameters
    ----------
    symbol : str
        Atomic symbol (e.g. ``'C'``, ``'N'``).
    bond_count : int or float
        Sum of bond orders attached to this atom.
    charge : int
        Formal charge on the atom.

    Returns
    -------
    bool
    """
    entry = _VALENCE_LIMITS.get(symbol)
    if entry is None:
        return False
    max_bonds, neutral_only = entry
    if neutral_only and charge != 0:
        return False
    return bond_count > max_bonds
