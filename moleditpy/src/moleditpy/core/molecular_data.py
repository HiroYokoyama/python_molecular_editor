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
import math
from typing import Any, Dict, List, Optional, Tuple, Union
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ..utils.constants import ANGSTROM_PER_PIXEL


class PointTuple(Tuple[float, ...]):
    """Backward-compatible tuple that allows .x() and .y() access like QPointF."""

    def x(self) -> float:
        """Return the x coordinate as a float."""
        return float(self[0])

    def y(self) -> float:
        """Return the y coordinate as a float."""
        return float(self[1])


def _mdl_radical_code(radical_electrons: int) -> int:
    """Map a count of unpaired electrons to an MDL radical code."""
    # MDL encodes 1 = singlet, 2 = doublet, 3 = triplet.
    if radical_electrons <= 0:
        return 0
    return 2 if radical_electrons == 1 else 3


def _mdl_property_lines(atom_records: List[Dict[str, Any]]) -> List[str]:
    """Build the M CHG / M RAD property lines for a V2000 block."""
    lines = []
    for tag, key, encode in (
        ("CHG", "charge", int),
        ("RAD", "radical", _mdl_radical_code),
    ):
        entries = [
            (idx, encode(rec[key]))
            for idx, rec in enumerate(atom_records, start=1)
            if rec[key]
        ]
        # The MDL property block allows at most 8 entries per line.
        for start in range(0, len(entries), 8):
            chunk = entries[start : start + 8]
            lines.append(
                f"M  {tag}{len(chunk):3d}"
                + "".join(f"{idx:4d}{val:4d}" for idx, val in chunk)
            )
    return lines


def _to_v2000_block(
    atom_records: List[Dict[str, Any]], bond_records: List[Dict[str, Any]]
) -> str:
    """Render collected atom/bond records as an MDL V2000 MOL block."""
    lines = ["", "  MoleditPy", ""]
    lines.append(
        f"{len(atom_records):3d}{len(bond_records):3d}  0  0  0  0  0  0  0  0999 V2000"
    )
    for rec in atom_records:
        # Charge and radical travel in the M CHG / M RAD blocks below, so every
        # per-atom field here stays zero.
        lines.append(
            f"{rec['x']:10.4f}{rec['y']:10.4f}{0.0:10.4f} {rec['symbol']:<3}"
            " 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for rec in bond_records:
        lines.append(
            f"{rec['a1']:3d}{rec['a2']:3d}{rec['order']:3d}{rec['stereo']:3d}  0  0  0"
        )
    lines.extend(_mdl_property_lines(atom_records))
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def _to_v3000_block(
    atom_records: List[Dict[str, Any]], bond_records: List[Dict[str, Any]]
) -> str:
    """Render collected atom/bond records as an MDL V3000 MOL block."""
    lines = ["", "  MoleditPy", "", "  0  0  0  0  0  0  0  0  0  0999 V3000"]
    lines.append("M  V30 BEGIN CTAB")
    lines.append(f"M  V30 COUNTS {len(atom_records)} {len(bond_records)} 0 0 0")
    lines.append("M  V30 BEGIN ATOM")
    for idx, rec in enumerate(atom_records, start=1):
        extra = ""
        if rec["charge"]:
            extra += f" CHG={int(rec['charge'])}"
        if rec["radical"]:
            extra += f" RAD={_mdl_radical_code(rec['radical'])}"
        lines.append(
            f"M  V30 {idx} {rec['symbol']} "
            f"{rec['x']:.4f} {rec['y']:.4f} 0.0000 0{extra}"
        )
    lines.append("M  V30 END ATOM")
    lines.append("M  V30 BEGIN BOND")
    for idx, rec in enumerate(bond_records, start=1):
        cfg = ""
        if rec["stereo"] == 1:
            cfg = " CFG=1"
        elif rec["stereo"] == 6:
            cfg = " CFG=3"
        lines.append(f"M  V30 {idx} {rec['order']} {rec['a1']} {rec['a2']}{cfg}")
    lines.append("M  V30 END BOND")
    lines.append("M  V30 END CTAB")
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def align_ez_to_cip(mol: Chem.Mol, wanted: Dict[int, str]) -> None:
    """Make labelled double bonds read as the wanted CIP E/Z, flipping if needed.

    *wanted* maps bond index to "E" or "Z". The bonds' STEREOE/STEREOZ is
    encoded against whichever neighbors were picked as stereo atoms, which
    need not be the CIP-higher ones (on CH3-C(Cl)=CH-CH3 the methyl is picked
    over Cl), so a "Z" label could encode E. Ask the CIP labeler what each
    bond currently means and flip the ones that read the other way. Bonds
    with no CIP E/Z (two identical substituents on one end) are left alone.
    """
    if not wanted:
        return
    probe = Chem.Mol(mol)
    try:
        rdCIPLabeler.AssignCIPLabels(probe)
    except (RuntimeError, ValueError) as e:
        logging.warning("CIP labelling for E/Z labels failed: %s", e)
        return
    flipped = {
        Chem.BondStereo.STEREOZ: Chem.BondStereo.STEREOE,
        Chem.BondStereo.STEREOE: Chem.BondStereo.STEREOZ,
    }
    for bond_idx, label in wanted.items():
        probe_bond = probe.GetBondWithIdx(bond_idx)
        if not probe_bond.HasProp("_CIPCode"):
            continue
        if probe_bond.GetProp("_CIPCode") != label:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetStereo() in flipped:
                bond.SetStereo(flipped[bond.GetStereo()])


_WEDGE_Z = {Chem.BondDir.BEGINWEDGE: 1.0, Chem.BondDir.BEGINDASH: -1.0}

# How far past a straight line the two bonds without a wedge may open and
# still count as a T-shape rather than a fan the wedge sits inside.
_T_SHAPE_TOLERANCE_DEG = 5.0


def _tag_from_points(
    mol: Chem.Mol,
    center: int,
    points: Dict[int, Tuple[float, float, float]],
    virtual: Optional[Tuple[float, float, float]] = None,
) -> Chem.ChiralType:
    """Chiral tag of *center* with its neighbors at *points* (center-relative).

    *virtual* is the position of an implicit H or lone pair. The center is
    placed at the centroid of the four positions, inside the tetrahedron, so
    any three neighbor vectors are independent.
    """
    corners = list(points.values()) + ([virtual] if virtual is not None else [])
    centroid = tuple(sum(p[k] for p in corners) / len(corners) for k in range(3))
    probe = Chem.Mol(mol)
    conf = Chem.Conformer(probe.GetNumAtoms())
    conf.Set3D(True)
    for idx, pos in points.items():
        conf.SetAtomPosition(idx, pos)
    conf.SetAtomPosition(center, centroid)
    probe.RemoveAllConformers()
    conf_id = probe.AddConformer(conf, assignId=True)
    Chem.AssignAtomChiralTagsFromStructure(probe, confId=conf_id)
    return probe.GetAtomWithIdx(center).GetChiralTag()


def _largest_gap_deg(directions: List[Tuple[float, float]]) -> float:
    """Largest angle between angularly adjacent 2D directions, in degrees."""
    angles = sorted(math.degrees(math.atan2(dy, dx)) % 360.0 for dx, dy in directions)
    gaps = [b - a for a, b in zip(angles, angles[1:])]
    gaps.append(360.0 - angles[-1] + angles[0])
    return max(gaps)


def wedge_reading(mol: Chem.Mol, center: int) -> Optional[Chem.ChiralType]:
    """Chiral tag the wedges/hashes starting at *center* specify, or None.

    Reads the drawing the way a chemist does: plain bonds in the paper plane,
    a wedge toward the viewer, a hash away. With three drawn neighbors the
    implicit H (or lone pair) points opposite the wedges. With four, each
    stereo bond is read in turn as the one out of the plane; the reading
    counts only when they all agree.

    None means the drawing does not decide: no stereo bond starts here, a
    wedge and a hash cancel on a three-neighbor center, the wedge sits inside
    a fan of under 180 degrees, or the stereo bonds of a four-neighbor center
    contradict each other (a wedge and a hash on opposite bonds).
    """
    atom = mol.GetAtomWithIdx(center)
    if atom.GetDegree() not in (3, 4) or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    origin = conf.GetAtomPosition(center)
    directions: Dict[int, Tuple[float, float]] = {}
    lift: Dict[int, float] = {}
    for bond in atom.GetBonds():
        nbr = bond.GetOtherAtomIdx(center)
        pos = conf.GetAtomPosition(nbr)
        dx, dy = pos.x - origin.x, pos.y - origin.y
        length = math.hypot(dx, dy)
        if length < 1e-6:
            return None
        directions[nbr] = (dx / length, dy / length)
        lift[nbr] = (
            _WEDGE_Z.get(bond.GetBondDir(), 0.0)
            if bond.GetBeginAtomIdx() == center
            else 0.0
        )
    stereo_nbrs = [nbr for nbr, z in lift.items() if z]
    if not stereo_nbrs:
        return None
    flat = {nbr: (dx, dy, 0.0) for nbr, (dx, dy) in directions.items()}

    if len(directions) == 3:
        net = sum(lift.values())
        if not net:
            return None
        if _largest_gap_deg(list(directions.values())) > 180.0 + _T_SHAPE_TOLERANCE_DEG:
            return None
        return _tag_from_points(mol, center, flat, (0.0, 0.0, -math.copysign(1.0, net)))

    readings = set()
    for nbr in stereo_nbrs:
        points = dict(flat)
        points[nbr] = (0.0, 0.0, lift[nbr])
        readings.add(_tag_from_points(mol, center, points))
    return readings.pop() if len(readings) == 1 else None


def _bond_states(mol: Chem.Mol, indices: List[int]) -> List[Tuple[int, int, Any]]:
    """Begin, end and direction of the bonds at *indices*."""
    return [
        (
            mol.GetBondWithIdx(i).GetBeginAtomIdx(),
            mol.GetBondWithIdx(i).GetEndAtomIdx(),
            mol.GetBondWithIdx(i).GetBondDir(),
        )
        for i in indices
    ]


def _rdkit_tags(mol: Chem.Mol) -> List[Chem.ChiralType]:
    """Chiral tag RDKit's wedge perception gives each atom of *mol*."""
    probe = Chem.Mol(mol)
    for atom in probe.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignChiralTypesFromBondDirs(probe)
    return [atom.GetChiralTag() for atom in probe.GetAtoms()]


def _clear_own_wedges(mol: Chem.Mol, centers: Dict[int, Chem.ChiralType]) -> None:
    """Remove the wedges/hashes that start at *centers*."""
    for center in centers:
        for bond in mol.GetAtomWithIdx(center).GetBonds():
            if bond.GetBeginAtomIdx() == center and bond.GetBondDir() in _WEDGE_Z:
                bond.SetBondDir(Chem.BondDir.NONE)


def _rdkit_rewedge(mol: Chem.Mol, fixes: Dict[int, Chem.ChiralType]) -> None:
    """Let RDKit wedge each fixed center for its wanted tag (may reverse bonds)."""
    _clear_own_wedges(mol, fixes)
    for center, tag in fixes.items():
        mol.GetAtomWithIdx(center).SetChiralTag(tag)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    for center in fixes:
        mol.GetAtomWithIdx(center).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)


def _pick_wedges(mol: Chem.Mol, fixes: Dict[int, Chem.ChiralType]) -> None:
    """Wedge or hash one free bond per fixed center so it reads its wanted tag.

    A candidate counts only when RDKit and wedge_reading() agree on it.
    """
    for center, tag in fixes.items():
        for bond in mol.GetAtomWithIdx(center).GetBonds():
            if (
                bond.GetBeginAtomIdx() != center
                or bond.GetBondDir() != Chem.BondDir.NONE
                or bond.GetBondType() != Chem.BondType.SINGLE
            ):
                continue
            for direction in _WEDGE_Z:
                bond.SetBondDir(direction)
                if (
                    _rdkit_tags(mol)[center] == tag
                    and wedge_reading(mol, center) == tag
                ):
                    break
                bond.SetBondDir(Chem.BondDir.NONE)
            else:
                continue
            break


def resolve_wedge_chirality(mol: Chem.Mol) -> None:
    """Re-wedge stereocenters whose drawn wedges RDKit reads differently.

    RDKit's wedge perception is unstable on some unambiguous drawings: with
    the two plain bonds of a three-neighbor center drawn in a straight line
    (a T-shape, which the editor's angle snapping produces), it reports
    "ambiguous stereochemistry" or inverts the center, while the same drawing
    bent by a few degrees reads correctly. Where wedge_reading() gives a
    definite answer that RDKit's perception does not reproduce, the wedges at
    that center are replaced by ones RDKit picks for the intended
    configuration on the same coordinates (or, where RDKit declines, a wedge
    or hash on a free bond of the center that both readings agree on), so
    every later reader (MOL block, 3D embedding, stereo check) gets the drawn
    isomer. Nothing is changed unless every other stereocenter and every
    bond away from the fixed centers stays exactly as drawn.
    """
    if mol.GetNumConformers() == 0:
        return
    centers = sorted(
        {b.GetBeginAtomIdx() for b in mol.GetBonds() if b.GetBondDir() in _WEDGE_Z}
    )
    if not centers:
        return
    try:
        read = _rdkit_tags(mol)
        fixes: Dict[int, Chem.ChiralType] = {}
        for center in centers:
            wanted = wedge_reading(mol, center)
            if wanted is not None and wanted != read[center]:
                fixes[center] = wanted
        if not fixes:
            return

        # Bonds the re-wedge must not touch: everything away from the fixed
        # centers, and any wedge/hash an adjacent center draws onto a bond it
        # shares with a fixed center.
        protected = [
            b.GetIdx()
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() not in fixes
            and (b.GetEndAtomIdx() not in fixes or b.GetBondDir() in _WEDGE_Z)
        ]
        before = _bond_states(mol, protected)
        expected = [fixes.get(i, tag) for i, tag in enumerate(read)]

        def holds(trial: Chem.Mol) -> bool:
            """Fixed centers read right; every other center and bond is as drawn."""
            return (
                _bond_states(trial, protected) == before
                and _rdkit_tags(trial) == expected
            )

        # WedgeMolBonds may reverse bonds, which cannot be undone, so try it
        # on a copy and touch the real molecule only once it holds.
        trial = Chem.Mol(mol)
        _rdkit_rewedge(trial, fixes)
        if holds(trial):
            _rdkit_rewedge(mol, fixes)
            return
        # WedgeMolBonds declines some centers, e.g. a T-shape whose one
        # unambiguous bond already carries a neighbor's hash. Try each free
        # bond of the center instead; this never reverses a bond.
        trial = Chem.Mol(mol)
        _clear_own_wedges(trial, fixes)
        _pick_wedges(trial, fixes)
        if holds(trial):
            for center in fixes:
                for bond in mol.GetAtomWithIdx(center).GetBonds():
                    if bond.GetBeginAtomIdx() == center:
                        bond.SetBondDir(
                            trial.GetBondWithIdx(bond.GetIdx()).GetBondDir()
                        )
            return
        logging.debug("Re-wedging stereocenters %s did not hold", sorted(fixes))
    except (RuntimeError, ValueError) as e:
        logging.warning("Resolving wedge stereochemistry failed: %s", e)


class MolecularData:
    """In-memory graph of atoms and bonds, independent of any UI framework."""

    atoms: Dict[int, Dict[str, Any]]
    bonds: Dict[Tuple[int, int], Dict[str, Any]]
    adjacency_list: Dict[int, List[int]]
    next_atom_id: int

    def __init__(self) -> None:
        """Initialize MolecularData structure with atoms, bonds, and coordinates."""
        self.atoms = {}
        self.bonds = {}
        self.next_atom_id = 0
        self.adjacency_list = {}

    def add_atom(
        self,
        symbol: str,
        pos: Union[Any, Tuple[float, float]],
        charge: int = 0,
        radical: int = 0,
    ) -> int:
        """Add a new atom and return its auto-assigned integer ID."""
        atom_id = self.next_atom_id
        # Internalize position as raw floats to decouple from UI types (QPointF)
        if hasattr(pos, "x") and hasattr(pos, "y"):
            raw_pos = PointTuple((float(pos.x()), float(pos.y())))
        else:
            raw_pos = PointTuple((float(pos[0]), float(pos[1])))

        self.atoms[atom_id] = {
            "symbol": symbol,
            "pos": raw_pos,
            "charge": charge,
            "radical": radical,
        }
        self.adjacency_list[atom_id] = []
        self.next_atom_id += 1
        return atom_id

    def set_atom_pos(self, atom_id: int, pos: Union[Any, Tuple[float, float]]) -> None:
        """Update atom position using raw floats or QPointF."""
        if atom_id in self.atoms:
            if hasattr(pos, "x") and hasattr(pos, "y"):
                self.atoms[atom_id]["pos"] = PointTuple(
                    (float(pos.x()), float(pos.y()))
                )
            else:
                self.atoms[atom_id]["pos"] = PointTuple((float(pos[0]), float(pos[1])))

    def add_bond(
        self, id1: int, id2: int, order: Union[int, float] = 1, stereo: int = 0
    ) -> Tuple[Tuple[int, int], str]:
        """Add or update a bond between two atoms; return its canonical key and 'created'/'updated'."""
        # For stereo bonds, do not sort because ID order determines direction.
        # For non-stereo bonds, sort to normalize the key.
        if stereo == 0:
            if id1 > id2:
                id1, id2 = id2, id1

        bond_data = {"order": order, "stereo": stereo}

        # Check if it's a new bond, considering reverse direction keys.
        is_new_bond = (id1, id2) not in self.bonds and (id2, id1) not in self.bonds
        if is_new_bond:
            if id1 in self.adjacency_list and id2 in self.adjacency_list:
                self.adjacency_list[id1].append(id2)
                self.adjacency_list[id2].append(id1)

        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)].update(bond_data)
            return (id1, id2), "updated"
        if (id2, id1) in self.bonds:
            # Re-key reversed entry (stereo direction changed) to avoid duplicates
            existing = self.bonds.pop((id2, id1))
            existing.update(bond_data)
            self.bonds[(id1, id2)] = existing
            return (id1, id2), "updated"
        self.bonds[(id1, id2)] = bond_data
        return (id1, id2), "created"

    def remove_atom(self, atom_id: int) -> None:
        """Remove an atom and all bonds involving it from the data model."""
        if atom_id in self.atoms:
            # Safely get neighbors before deleting the atom's own entry
            neighbors = self.adjacency_list.get(atom_id, [])
            for neighbor_id in neighbors:
                if (
                    neighbor_id in self.adjacency_list
                    and atom_id in self.adjacency_list[neighbor_id]
                ):
                    self.adjacency_list[neighbor_id].remove(atom_id)

            # Now, safely delete the atom's own entry from the adjacency list
            if atom_id in self.adjacency_list:
                del self.adjacency_list[atom_id]

            if atom_id in self.atoms:
                del self.atoms[atom_id]

            # Remove bonds involving this atom
            bonds_to_remove = [key for key in self.bonds if atom_id in key]
            for key in bonds_to_remove:
                self.bonds.pop(key, None)

    def remove_bond(self, id1: int, id2: int) -> None:
        """Remove the bond between two atoms, handling stereo and non-stereo key variants."""
        # Look for directional stereo bonds (forward/reverse) and normalized non-stereo bond keys.
        key_to_remove = None
        if (id1, id2) in self.bonds:
            key_to_remove = (id1, id2)
        elif (id2, id1) in self.bonds:
            key_to_remove = (id2, id1)

        if key_to_remove:
            if id1 in self.adjacency_list and id2 in self.adjacency_list[id1]:
                self.adjacency_list[id1].remove(id2)
            if id2 in self.adjacency_list and id1 in self.adjacency_list[id2]:
                self.adjacency_list[id2].remove(id1)
            self.bonds.pop(key_to_remove, None)

    def to_rdkit_mol(self, use_2d_stereo: bool = True) -> Optional[Chem.Mol]:
        """
        use_2d_stereo: True estimates E/Z from 2D coordinates (as before). False prioritizes E/Z labels.
        Call with use_2d_stereo=False for 3D conversion.
        """
        if not self.atoms:
            return None
        mol = Chem.RWMol()

        # atoms
        atom_id_to_idx_map = {}
        for atom_id, data in self.atoms.items():
            try:
                atom = Chem.Atom(data["symbol"])
            except (RuntimeError, ValueError):
                # RDKit doesn't support this symbol. Return None to trigger
                # manual MoleditPy fallback (with 'MoleditPy' header).
                return None
            atom.SetFormalCharge(data.get("charge", 0))
            atom.SetNumRadicalElectrons(data.get("radical", 0))
            atom.SetIntProp("_original_atom_id", atom_id)
            atom.SetNoImplicit(False)  # Allow RDKit to perceive implicit valence
            idx = mol.AddAtom(atom)
            atom_id_to_idx_map[atom_id] = idx

        # save bonds & stereo info (label info is kept here)
        bond_stereo_info: Dict[
            int, Dict[str, Any]
        ] = {}  # bond_idx -> {'type', 'atom_ids', 'bond_data'}
        for (id1, id2), bond_data in self.bonds.items():
            if id1 not in atom_id_to_idx_map or id2 not in atom_id_to_idx_map:
                continue
            idx1, idx2 = atom_id_to_idx_map[id1], atom_id_to_idx_map[id2]

            order_val = float(bond_data["order"])
            order = {
                1.0: Chem.BondType.SINGLE,
                1.5: Chem.BondType.AROMATIC,
                2.0: Chem.BondType.DOUBLE,
                3.0: Chem.BondType.TRIPLE,
            }.get(order_val, Chem.BondType.SINGLE)

            bond_idx = mol.AddBond(idx1, idx2, order) - 1

            # If stereo label exists, keep details for bond_idx (used later)
            if "stereo" in bond_data and bond_data["stereo"] in [1, 2, 3, 4]:
                bond_stereo_info[bond_idx] = {
                    "type": int(bond_data["stereo"]),
                    "atom_ids": (id1, id2),
                    "bond_data": bond_data,
                }

        # sanitize
        final_mol = mol.GetMol()
        try:
            Chem.SanitizeMol(final_mol)
        except (RuntimeError, ValueError, TypeError):
            # Sanitization failure: return None to trigger manual MOL block fallback
            return None

        # add 2D conformer
        # Convert from scene pixels to angstroms when creating RDKit conformer.
        conf = Chem.Conformer(final_mol.GetNumAtoms())
        conf.Set3D(False)
        for atom_id, data in self.atoms.items():
            if atom_id in atom_id_to_idx_map:
                idx = atom_id_to_idx_map[atom_id]
                pos = data.get("pos")
                if pos:
                    # pos may be a tuple (x, y) or a QPointF (from old state deserialization)
                    if hasattr(pos, "x") and hasattr(pos, "y"):
                        ax = pos.x() * ANGSTROM_PER_PIXEL
                        ay = -pos.y() * ANGSTROM_PER_PIXEL
                    elif isinstance(pos, (list, tuple)) and len(pos) >= 2:
                        ax = pos[0] * ANGSTROM_PER_PIXEL
                        ay = -pos[1] * ANGSTROM_PER_PIXEL
                    else:
                        continue
                    conf.SetAtomPosition(idx, (ax, ay, 0.0))
        final_mol.AddConformer(conf)

        # Stereochemistry setting prioritizing E/Z labels ---
        # First, record bonds with E/Z labels
        ez_labeled_bonds = set()
        for bond_idx, info in bond_stereo_info.items():
            if info["type"] in [3, 4]:
                ez_labeled_bonds.add(bond_idx)

        # Estimate E/Z from 2D coordinates only if use_2d_stereo=True and no E/Z label exists
        if use_2d_stereo:
            Chem.SetDoubleBondNeighborDirections(final_mol, final_mol.GetConformer(0))
        else:
            # 3D conversion: Disable coordinate-based estimation completely if E/Z labels exist
            if ez_labeled_bonds:
                # If E/Z labels exist, clear BondDir for all bonds to disable coordinate-based estimation
                for b in final_mol.GetBonds():
                    b.SetBondDir(Chem.BondDir.NONE)
            else:
                # Perform coordinate-based estimation only if no E/Z labels exist
                Chem.SetDoubleBondNeighborDirections(
                    final_mol, final_mol.GetConformer(0)
                )

        # Helper: Pick neighbors prioritizing heavy atoms
        def pick_preferred_neighbor(atom: Chem.Atom, exclude_idx: int) -> Optional[int]:
            """Select preferred neighbor atom for geometry alignment or reference."""
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                if nbr.GetAtomicNum() > 1:
                    return int(nbr.GetIdx())
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() != exclude_idx:
                    return int(nbr.GetIdx())
            return None

        # Overwrite based on labels (E/Z has highest priority) ---
        ez_wanted: Dict[int, str] = {}
        for bond_idx, info in bond_stereo_info.items():
            stereo_type = info["type"]
            bond = final_mol.GetBondWithIdx(bond_idx)

            # Case with single bond wedge/dash labels (1/2)
            if stereo_type in [1, 2]:
                if stereo_type == 1:
                    bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif stereo_type == 2:
                    bond.SetBondDir(Chem.BondDir.BEGINDASH)
                continue

            # Double bond E/Z labels (3/4)
            if stereo_type in [3, 4]:
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue

                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()

                label_bond_data = info.get("bond_data") or {}
                stereo_atoms_specified = label_bond_data.get("stereo_atoms")

                if stereo_atoms_specified:
                    try:
                        a1_id, a2_id = stereo_atoms_specified
                        neigh1_idx = atom_id_to_idx_map.get(a1_id)
                        neigh2_idx = atom_id_to_idx_map.get(a2_id)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        neigh1_idx = None
                        neigh2_idx = None
                else:
                    neigh1_idx = pick_preferred_neighbor(
                        final_mol.GetAtomWithIdx(begin_atom_idx), end_atom_idx
                    )
                    neigh2_idx = pick_preferred_neighbor(
                        final_mol.GetAtomWithIdx(end_atom_idx), begin_atom_idx
                    )

                if neigh1_idx is None or neigh2_idx is None:
                    continue

                bond.SetStereoAtoms(neigh1_idx, neigh2_idx)
                if stereo_type == 3:
                    bond.SetStereo(Chem.BondStereo.STEREOZ)
                elif stereo_type == 4:
                    bond.SetStereo(Chem.BondStereo.STEREOE)
                ez_wanted[bond_idx] = "Z" if stereo_type == 3 else "E"

                # Preserve explicit wedge/dash labels, including those applied later.
                for center_idx, neighbor_idx in (
                    (begin_atom_idx, neigh1_idx),
                    (end_atom_idx, neigh2_idx),
                ):
                    adjacent = final_mol.GetBondBetweenAtoms(center_idx, neighbor_idx)
                    if adjacent is not None and bond_stereo_info.get(
                        adjacent.GetIdx(), {}
                    ).get("type") not in (1, 2):
                        adjacent.SetBondDir(Chem.BondDir.NONE)

        resolve_wedge_chirality(final_mol)

        # The labels are CIP E/Z; the neighbors picked above are not always
        # the CIP-higher ones.
        align_ez_to_cip(final_mol, ez_wanted)

        # Finalization (cache update + stereochemistry reassignment)
        final_mol.UpdatePropertyCache(strict=False)

        # During 3D conversion (use_2d_stereo=False), apply force=True if E/Z labels exist
        if not use_2d_stereo and ez_labeled_bonds:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=True)
        else:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=False)
        return final_mol

    def to_mol_block(self) -> Optional[str]:
        """Serialize the molecule to an MDL MOL block string, or None on failure."""
        mol = self.to_rdkit_mol()
        if mol:
            try:
                return Chem.MolToMolBlock(mol, includeStereo=True)  # type: ignore[no-any-return]
            except (RuntimeError, ValueError, TypeError) as e:
                logging.warning(
                    f"RDKit MolBlock generation failed: {e}"
                )  # Fallback to manual MolBlock generation if RDKit fails

        if not self.atoms:
            return None

        # Counts line and bond indices must only reflect atoms actually written
        atom_map: Dict[int, int] = {}
        atom_records: List[Dict[str, Any]] = []
        for old_id, atom in self.atoms.items():
            # Convert scene pixel coordinates to angstroms when emitting MOL block
            pos = atom.get("pos")
            if not pos:
                continue

            if hasattr(pos, "x") and hasattr(pos, "y"):
                x_px, y_px = pos.x(), -pos.y()
            elif isinstance(pos, (list, tuple)) and len(pos) >= 2:
                x_px, y_px = pos[0], -pos[1]
            else:
                continue

            atom_map[old_id] = len(atom_records)
            atom_records.append(
                {
                    "x": x_px * ANGSTROM_PER_PIXEL,
                    "y": y_px * ANGSTROM_PER_PIXEL,
                    "symbol": atom["symbol"],
                    "charge": int(atom.get("charge", 0) or 0),
                    "radical": int(atom.get("radical", 0) or 0),
                }
            )

        bond_records: List[Dict[str, Any]] = []
        for (id1, id2), bond in self.bonds.items():
            if id1 not in atom_map or id2 not in atom_map:
                continue
            # Bond order may be a float (1.5 = aromatic); MDL uses code 4.
            order_val = float(bond["order"])
            stereo_code = 0
            bond_stereo = bond.get("stereo", 0)
            if bond_stereo == 1:
                stereo_code = 1
            elif bond_stereo == 2:
                stereo_code = 6

            bond_records.append(
                {
                    "a1": atom_map[id1] + 1,
                    "a2": atom_map[id2] + 1,
                    "order": 4 if order_val == 1.5 else int(order_val),
                    "stereo": stereo_code,
                }
            )

        # V2000 packs the counts into three-character fields, so anything past
        # 999 would run the atom count into the bond count and corrupt the file.
        if len(atom_records) > 999 or len(bond_records) > 999:
            return _to_v3000_block(atom_records, bond_records)
        return _to_v2000_block(atom_records, bond_records)

    def to_template_dict(
        self, name: str, version: str = "1.0", application_version: str = ""
    ) -> Dict[str, Any]:
        """Convert current structure to a dictionary for template storage."""
        import datetime

        atoms_data = []
        for atom_id, atom_info in self.atoms.items():
            pos = atom_info["pos"]
            # pos is guaranteed to be a tuple (x, y) due to internal caching in add_atom/set_atom_pos
            atoms_data.append(
                {
                    "id": atom_id,
                    "symbol": atom_info["symbol"],
                    "x": pos[0],
                    "y": pos[1],
                    "charge": atom_info.get("charge", 0),
                    "radical": atom_info.get("radical", 0),
                }
            )

        bonds_data = []
        for (id1, id2), bond_info in self.bonds.items():
            bonds_data.append(
                {
                    "atom1": id1,
                    "atom2": id2,
                    "order": bond_info["order"],
                    "stereo": bond_info.get("stereo", 0),
                }
            )

        return {
            "format": "PME Template",
            "version": version,
            "application": "MoleditPy",
            "application_version": application_version,
            "name": name,
            "created": str(datetime.datetime.now()),
            "atoms": atoms_data,
            "bonds": bonds_data,
        }
