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
from typing import Any, Dict, List, Optional, Tuple, Union
from rdkit import Chem

from ..utils.constants import ANGSTROM_PER_PIXEL


class PointTuple(tuple):
    """Backward-compatible tuple that allows .x() and .y() access like QPointF."""

    def x(self) -> float:
        """Return the x coordinate as a float."""
        return float(self[0])

    def y(self) -> float:
        """Return the y coordinate as a float."""
        return float(self[1])


class MolecularData:
    """In-memory graph of atoms and bonds, independent of any UI framework."""

    atoms: Dict[int, Dict[str, Any]]
    bonds: Dict[Tuple[int, int], Dict[str, Any]]
    adjacency_list: Dict[int, List[int]]
    next_atom_id: int

    def __init__(self) -> None:
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
            # Same bond stored under the reversed key (e.g. the stereo
            # direction changed): re-key it instead of creating a duplicate.
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
        bond_stereo_info = {}  # bond_idx -> {'type': int, 'atom_ids': (id1,id2), 'bond_data': bond_data}
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

                bond_data: dict[str, Any] = info.get("bond_data") or {}  # type: ignore[assignment]
                stereo_atoms_specified = bond_data.get("stereo_atoms")

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

                # Clear BondDir (wedge/dash) of adjacent single bonds assigned via coordinates to avoid label conflicts
                b1 = final_mol.GetBondBetweenAtoms(begin_atom_idx, neigh1_idx)
                b2 = final_mol.GetBondBetweenAtoms(end_atom_idx, neigh2_idx)
                if b1 is not None:
                    b1.SetBondDir(Chem.BondDir.NONE)
                if b2 is not None:
                    b2.SetBondDir(Chem.BondDir.NONE)

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

        # Emit atom lines first so the counts line and bond indices only
        # reflect atoms that actually have a usable position.
        atom_map = {}
        atom_lines = []
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

            x, y = x_px * ANGSTROM_PER_PIXEL, y_px * ANGSTROM_PER_PIXEL
            z, symbol = 0.0, atom["symbol"]
            charge = atom.get("charge", 0)

            chg_code = 0
            if charge == 3:
                chg_code = 1
            elif charge == 2:
                chg_code = 2
            elif charge == 1:
                chg_code = 3
            elif charge == -1:
                chg_code = 5
            elif charge == -2:
                chg_code = 6
            elif charge == -3:
                chg_code = 7

            atom_map[old_id] = len(atom_lines)
            atom_lines.append(
                f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0{chg_code:3d}  0  0  0  0  0  0  0\n"
            )

        bond_lines = []
        for (id1, id2), bond in self.bonds.items():
            if id1 not in atom_map or id2 not in atom_map:
                continue
            idx1, idx2 = atom_map[id1] + 1, atom_map[id2] + 1
            # Bond order may be a float (1.5 = aromatic); V2000 uses code 4.
            order_val = float(bond["order"])
            order = 4 if order_val == 1.5 else int(order_val)
            stereo_code = 0
            bond_stereo = bond.get("stereo", 0)
            if bond_stereo == 1:
                stereo_code = 1
            elif bond_stereo == 2:
                stereo_code = 6

            bond_lines.append(f"{idx1:3d}{idx2:3d}{order:3d}{stereo_code:3d}  0  0  0\n")

        mol_block = "\n  MoleditPy\n\n"
        mol_block += (
            f"{len(atom_lines):3d}{len(bond_lines):3d}  0  0  0  0  0  0  0  0999 V2000\n"
        )
        mol_block += "".join(atom_lines)
        mol_block += "".join(bond_lines)
        mol_block += "M  END\n"
        return mol_block

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
