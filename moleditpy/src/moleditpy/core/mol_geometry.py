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
# Primitive geometry helpers
# ------------------------------------------------------------------


def calc_distance(pos1, pos2) -> float:
    """Return the Euclidean distance between two 3-D positions.

    Parameters
    ----------
    pos1, pos2 : array-like, shape (3,)
        Cartesian coordinates.

    Returns
    -------
    float
        Distance in the same units as the input coordinates.
    """
    return float(
        np.linalg.norm(np.asarray(pos2, dtype=float) - np.asarray(pos1, dtype=float))
    )


def calc_angle_deg(pos1, pos2_vertex, pos3) -> float:
    """Return the angle pos1–pos2_vertex–pos3 in degrees.

    The angle is measured at *pos2_vertex* and is always in [0, 180].
    Numerical clipping is applied before ``arccos`` to guard against
    floating-point values slightly outside [−1, 1].

    Parameters
    ----------
    pos1, pos2_vertex, pos3 : array-like, shape (3,)
        Cartesian coordinates of the three atoms.

    Returns
    -------
    float
        Angle in degrees.
    """
    a = np.asarray(pos1, dtype=float)
    b = np.asarray(pos2_vertex, dtype=float)
    c = np.asarray(pos3, dtype=float)
    vec1 = a - b
    vec2 = c - b
    denom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    if denom == 0.0:
        return 0.0
    cos_a = np.clip(np.dot(vec1, vec2) / denom, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_a)))


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
# Bond-angle adjustment (difference-based rotation)
# ------------------------------------------------------------------


def rodrigues_rotate(v, axis, angle):
    """Rotate vector *v* around a unit *axis* by *angle* radians.

    Implements Rodrigues' rotation formula:

        v_rot = v·cos θ + (axis × v)·sin θ + axis·(axis · v)·(1 − cos θ)

    Parameters
    ----------
    v : ndarray, shape (3,)
        The vector to rotate.
    axis : ndarray, shape (3,)
        Unit rotation axis.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    ndarray, shape (3,)
        The rotated vector.
    """
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)


def adjust_bond_angle(
    positions, idx_a, idx_b, idx_c, target_angle_deg, atom_indices_to_move
):
    """Adjust the A–B–C bond angle to *target_angle_deg* using a
    difference-based rotation.

    The algorithm avoids 3D rotational ambiguity by computing the
    difference Δθ between the current angle and the target, then
    rotating only by that delta around the normal to the A-B-C plane.

    Steps
    -----
    1. Compute BA = A − B and BC = C − B.
    2. Calculate the current angle θ_current from BA and BC.
    3. Obtain the plane normal  n = BA × BC  (normalised).
       If |n| ≈ 0 the atoms are collinear and no unique rotation
       plane exists; the function returns without modifying anything.
    4. Δθ = target_angle (rad) − θ_current.
    5. Translate all movable atoms so that B is at the origin.
    6. Rotate each movable atom by Δθ around n using Rodrigues'
       rotation formula.
    7. Translate back so that B returns to its original position.

    Parameters
    ----------
    positions : ndarray, shape (N, 3)
        Atom coordinates.  **Modified in-place.**
    idx_a, idx_b, idx_c : int
        Atom indices defining the angle.  *idx_b* is the vertex
        (central) atom.
    target_angle_deg : float
        Desired bond angle in degrees.
    atom_indices_to_move : iterable of int
        Indices of atoms whose coordinates will be rotated (typically
        atom C and its connected sub-structure).

    Returns
    -------
    float
        The applied rotation Δθ in **radians**.  Returns 0.0 when the
        atoms are collinear and no rotation was applied.
    """
    pos_a = positions[idx_a]
    pos_b = positions[idx_b]  # vertex
    pos_c = positions[idx_c]

    # Vectors from vertex to the two arms
    ba = pos_a - pos_b
    bc = pos_c - pos_b

    # Current angle
    norm_ba = np.linalg.norm(ba)
    norm_bc = np.linalg.norm(bc)
    if norm_ba == 0 or norm_bc == 0:
        return 0.0

    cos_current = np.clip(np.dot(ba, bc) / (norm_ba * norm_bc), -1.0, 1.0)
    current_angle_rad = np.arccos(cos_current)

    # Normal vector of the A-B-C plane
    rotation_axis = np.cross(ba, bc)
    axis_norm = np.linalg.norm(rotation_axis)

    if axis_norm < 1e-12:
        # Atoms are collinear — pick an arbitrary axis perpendicular
        # to BA so the rotation can still proceed.
        ba_unit = ba / norm_ba
        # Choose a candidate not parallel to BA
        candidate = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(ba_unit, candidate)) > 0.9:
            candidate = np.array([0.0, 1.0, 0.0])
        rotation_axis = np.cross(ba_unit, candidate)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    else:
        rotation_axis = rotation_axis / axis_norm

    # Delta rotation
    target_angle_rad = np.radians(target_angle_deg)
    delta_theta = target_angle_rad - current_angle_rad

    # Translate so B is at origin, rotate, then translate back
    for idx in atom_indices_to_move:
        rel = positions[idx] - pos_b  # step 5: translate
        rotated = rodrigues_rotate(rel, rotation_axis, delta_theta)  # step 6
        positions[idx] = pos_b + rotated  # step 7: translate back

    return float(delta_theta)


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


# ------------------------------------------------------------------
# MOL Block & Chemical Logic (Pure)
# ------------------------------------------------------------------


def inject_ez_stereo_to_mol_block(mol_block, rdkit_mol, bonds_data):
    """Generate a modified MOL block with 'M CFG' lines for E/Z stereochemistry.

    Parameters
    ----------
    mol_block : str
        The original MOL block string.
    rdkit_mol : rdkit.Chem.Mol
        The RDKit molecule object (must have '_original_atom_id' properties).
    bonds_data : dict
        The bonds dictionary from MolecularData (key is tuple of atom IDs).

    Returns
    -------
    str
        The modified MOL block string.
    """
    mol_lines = mol_block.split("\n")
    ez_bond_info = {}

    # Map atom_id -> rdkit_idx once
    atom_id_to_rdkit_idx = {}
    for atom in rdkit_mol.GetAtoms():
        if atom.HasProp("_original_atom_id"):
            atom_id_to_rdkit_idx[atom.GetIntProp("_original_atom_id")] = atom.GetIdx()

    for (id1, id2), bond_info in bonds_data.items():
        stereo_type = bond_info.get("stereo", 0)
        if stereo_type in [3, 4]:  # E/Z labels
            idx1 = atom_id_to_rdkit_idx.get(id1)
            idx2 = atom_id_to_rdkit_idx.get(id2)

            if idx1 is not None and idx2 is not None:
                from rdkit import Chem

                rdkit_bond = rdkit_mol.GetBondBetweenAtoms(idx1, idx2)
                if rdkit_bond and rdkit_bond.GetBondType() == Chem.BondType.DOUBLE:
                    ez_bond_info[rdkit_bond.GetIdx()] = stereo_type

    if not ez_bond_info:
        return mol_block

    # Insert CFG lines before M  END
    insert_idx = -1
    for i, line in enumerate(mol_lines):
        if "M  END" in line:
            insert_idx = i
            break

    if insert_idx == -1:
        insert_idx = len(mol_lines)

    for bond_idx, stereo_type in sorted(ez_bond_info.items()):
        cfg_value = 1 if stereo_type == 3 else 2  # 1=Z, 2=E in MOL format
        cfg_line = f"M  CFG  1 {bond_idx + 1:3d}   {cfg_value}"
        mol_lines.insert(insert_idx, cfg_line)
        insert_idx += 1

    return "\n".join(mol_lines)


def identify_valence_problems(atoms_data, bonds_data):
    """Identify atoms with problematic valence.

    Parameters
    ----------
    atoms_data : dict
        The atoms dictionary from MolecularData.
    bonds_data : dict
        The bonds dictionary from MolecularData.

    Returns
    -------
    list
        IDs of atoms with problematic valence.
    """
    problem_atom_ids = []

    # Pre-calculate bond orders per atom
    bond_orders = {}
    for (id1, id2), bond in bonds_data.items():
        order = bond.get("order", 1)
        bond_orders[id1] = bond_orders.get(id1, 0) + order
        bond_orders[id2] = bond_orders.get(id2, 0) + order

    for atom_id, data in atoms_data.items():
        symbol = data["symbol"]
        charge = data.get("charge", 0)
        total_order = bond_orders.get(atom_id, 0)

        if is_problematic_valence(symbol, total_order, charge):
            problem_atom_ids.append(atom_id)

    return problem_atom_ids
