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
from rdkit import Chem

try:
    from ..utils.constants import ANGSTROM_PER_PIXEL
except ImportError:
    from moleditpy.utils.constants import ANGSTROM_PER_PIXEL


class PointTuple(tuple):
    """Backward-compatible tuple that allows .x() and .y() access like QPointF."""
    def x(self):
        return self[0]
    def y(self):
        return self[1]

class MolecularData:
    def __init__(self):
        self.atoms = {}
        self.bonds = {}
        self._next_atom_id = 0
        self.adjacency_list = {}

    def add_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self._next_atom_id
        # Internalize position as raw floats to decouple from UI types (QPointF)
        if hasattr(pos, "x") and hasattr(pos, "y"):
            raw_pos = PointTuple((float(pos.x()), float(pos.y())))
        else:
            raw_pos = PointTuple((float(pos[0]), float(pos[1])))

        self.atoms[atom_id] = {
            "symbol": symbol,
            "pos": raw_pos,
            "item": None,
            "charge": charge,
            "radical": radical,
        }
        self.adjacency_list[atom_id] = []
        self._next_atom_id += 1
        return atom_id

    def set_atom_pos(self, atom_id, pos):
        """Update atom position using raw floats or QPointF."""
        if atom_id in self.atoms:
            if hasattr(pos, "x") and hasattr(pos, "y"):
                self.atoms[atom_id]["pos"] = PointTuple((float(pos.x()), float(pos.y())))
            else:
                self.atoms[atom_id]["pos"] = PointTuple((float(pos[0]), float(pos[1])))

    def add_bond(self, id1, id2, order=1, stereo=0):
        # For stereo bonds, do not sort because ID order determines direction.
        # For non-stereo bonds, sort to normalize the key.
        if stereo == 0:
            if id1 > id2:
                id1, id2 = id2, id1

        bond_data = {"order": order, "stereo": stereo, "item": None}

        # Check if it's a new bond, considering reverse direction keys.
        is_new_bond = (id1, id2) not in self.bonds and (id2, id1) not in self.bonds
        if is_new_bond:
            if id1 in self.adjacency_list and id2 in self.adjacency_list:
                self.adjacency_list[id1].append(id2)
                self.adjacency_list[id2].append(id1)

        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)].update(bond_data)
            return (id1, id2), "updated"
        else:
            self.bonds[(id1, id2)] = bond_data
            return (id1, id2), "created"

    def remove_atom(self, atom_id):
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

    def remove_bond(self, id1, id2):
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

    def to_rdkit_mol(self, use_2d_stereo=True):
        """
        use_2d_stereo: True estimates E/Z from 2D coordinates (as before). False prioritizes E/Z labels.
        Call with use_2d_stereo=False for 3D conversion.
        """
        if not self.atoms:
            return None
        mol = Chem.RWMol()

        # atoms ---
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
            idx = mol.AddAtom(atom)
            atom_id_to_idx_map[atom_id] = idx

        # save bonds & stereo info (label info is kept here) ---
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

        # sanitize ---
        final_mol = mol.GetMol()
        try:
            Chem.SanitizeMol(final_mol)
        except (RuntimeError, ValueError, TypeError):
            # Sanitization failure: return None to trigger manual MOL block fallback
            return None

        # add 2D conformer ---
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
        def pick_preferred_neighbor(atom, exclude_idx):
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                if nbr.GetAtomicNum() > 1:
                    return nbr.GetIdx()
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() != exclude_idx:
                    return nbr.GetIdx()
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

                bond_data = info.get("bond_data", {}) or {}
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

    def update_ring_info_2d(self):
        """Update is_in_ring and ring_center for all BondItems based on 2D topology."""
        if not self.atoms or not self.bonds:
            return

        # 1. Generate RDKit molecule for topology analysis
        mol = self.to_rdkit_mol(use_2d_stereo=False)
        if not mol:
            # Fallback: reset all ring info if molecule generation fails
            for bond_data in self.bonds.values():
                bond_item = bond_data.get("item")
                if bond_item:
                    bond_item.is_in_ring = False
                    bond_item.ring_center = None
            return

        # 2. Extract ring information
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        bond_rings = ring_info.BondRings()

        # 3. Create mapping from RDKit atom index to editor atom item
        rdkit_idx_to_item = {}
        for atom in mol.GetAtoms():
            if atom.HasProp("_original_atom_id"):
                orig_id = atom.GetIntProp("_original_atom_id")
                if orig_id in self.atoms:
                    rdkit_idx_to_item[atom.GetIdx()] = self.atoms[orig_id]["item"]

        # 4. Map RDKit bond index to editor bond item
        rdkit_bond_idx_to_item = {}
        for bidx, rdkit_bond in enumerate(mol.GetBonds()):
            a1_idx = rdkit_bond.GetBeginAtomIdx()
            a2_idx = rdkit_bond.GetEndAtomIdx()
            if a1_idx in rdkit_idx_to_item and a2_idx in rdkit_idx_to_item:
                # Find corresponding editor bond item
                # This is slightly expensive but done once per update
                item1 = rdkit_idx_to_item[a1_idx]
                item2 = rdkit_idx_to_item[a2_idx]
                id1, id2 = item1.atom_id, item2.atom_id
                key = (id1, id2) if (id1, id2) in self.bonds else (id2, id1)
                if key in self.bonds:
                    rdkit_bond_idx_to_item[bidx] = self.bonds[key].get("item")

        # 5. Initialize/Reset all bond items
        for bond_data in self.bonds.values():
            bond_item = bond_data.get("item")
            if bond_item:
                bond_item.is_in_ring = False
                bond_item.ring_center = None

        # 6. Apply ring information
        for a_ring, b_ring in zip(atom_rings, bond_rings):
            # Calculate ring center (geometric mean of atom positions)
            positions = []
            for aidx in a_ring:
                item = rdkit_idx_to_item.get(aidx)
                if item and hasattr(item, "pos"):
                    pos = item.pos()
                    if pos is not None:
                        positions.append(pos)

            if not positions:
                continue

            center_x = sum(p.x() for p in positions) / len(positions)
            center_y = sum(p.y() for p in positions) / len(positions)
            ring_center = (center_x, center_y)  # Use tuple (x, y) instead of QPointF

            # Update all bonds in this ring
            for bidx in b_ring:
                bond_item = rdkit_bond_idx_to_item.get(bidx)
                if bond_item:
                    bond_item.is_in_ring = True
                    # Note: Simplified; a bond might be part of multiple rings.
                    # The inner-bond logic usually picks the smallest ring.
                    # Since we iterate through all rings, the last one wins.
                    # RDKit's AtomRings returns SSSR (Smallest Set of Smallest Rings),
                    # so this is usually correct for 2D drawing.
                    bond_item.ring_center = ring_center

    def to_mol_block(self):
        mol = self.to_rdkit_mol()
        if mol:
            try:
                return Chem.MolToMolBlock(mol, includeStereo=True)
            except (RuntimeError, ValueError, TypeError) as e:
                logging.warning(
                    f"RDKit MolBlock generation failed: {e}"
                )  # Fallback to manual MolBlock generation if RDKit fails

        if not self.atoms:
            return None
        atom_map = {old_id: new_id for new_id, old_id in enumerate(self.atoms.keys())}
        num_atoms, num_bonds = len(self.atoms), len(self.bonds)
        mol_block = "\n  MoleditPy\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
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

            mol_block += f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0{chg_code:3d}  0  0  0  0  0  0  0\n"

        for (id1, id2), bond in self.bonds.items():
            idx1, idx2, order = atom_map[id1] + 1, atom_map[id2] + 1, bond["order"]
            stereo_code = 0
            bond_stereo = bond.get("stereo", 0)
            if bond_stereo == 1:
                stereo_code = 1
            elif bond_stereo == 2:
                stereo_code = 6

            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}{stereo_code:3d}  0  0  0\n"

        mol_block += "M  END\n"
        return mol_block

    def to_template_dict(self, name, version="1.0", application_version=""):
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
