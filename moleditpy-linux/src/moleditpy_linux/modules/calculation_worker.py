#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import math
import re
import numpy as np

from PyQt6.QtCore import QObject, pyqtSignal, pyqtSlot

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem, rdGeometry
from rdkit.DistanceGeometry import DoTriangleSmoothing

try: # pragma: no cover
    from . import OBABEL_AVAILABLE
except Exception: # pragma: no cover
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE: # pragma: no cover
    try:
        import os
        import glob
        import openbabel
        from openbabel import pybel
        
        # The python wheel often misses setting up the data directory.
        # Check multiple potential locations for BABEL_DATADIR.
        ob_base = os.path.dirname(openbabel.__file__)
        data_candidates = [
            os.path.join(ob_base, "bin", "data"), # Typical for Windows wheel
            os.path.join(ob_base, "share", "openbabel", "*") # Typical for macOS/Linux
        ]
        
        found_datadir = None
        for pattern in data_candidates:
            matches = glob.glob(pattern)
            for m in matches:
                if os.path.isdir(m):
                    found_datadir = m
                    break
            if found_datadir:
                break
        
        if found_datadir:
            os.environ["BABEL_DATADIR"] = found_datadir
        
        os.environ["BABEL_LIBDIR"] = ob_base
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print(
            "Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled."
        )
else: # pragma: no cover
    pybel = None


class WorkerHaltError(Exception):
    """Custom exception raised when a calculation worker is requested to halt."""
    pass


# Human-readable labels for optimization method keys (mirrors opt3d_method_labels in the UI layer)
_OPT_METHOD_LABELS = {
    "MMFF_RDKIT":     "MMFF94s (RDKit)",
    "MMFF94_RDKIT":   "MMFF94 (RDKit)",
    "UFF_RDKIT":      "UFF (RDKit)",
    "MMFF94S_OBABEL": "MMFF94s (Open Babel)",
    "MMFF94_OBABEL":  "MMFF94 (Open Babel)",
    "UFF_OBABEL":     "UFF (Open Babel)",
    "GAFF_OBABEL":    "GAFF (Open Babel)",
    "GHEMICAL_OBABEL": "Ghemical (Open Babel)",
}


def _adjust_collision_avoidance(rd_mol, check_halted_cb, safe_status_cb):
    """
    Optimized collision avoidance using spatial partitioning (grid-based).
    This avoids the O(F^2 * N^2) complexity of the previous implementation.
    """
    try:
        frags = Chem.GetMolFrags(rd_mol, asMols=False, sanitizeFrags=False)
        if len(frags) <= 1:
            return

        safe_status_cb(f"Resolving potential collisions among {len(frags)} fragments...")

        conf = rd_mol.GetConformer()
        pt = Chem.GetPeriodicTable()

        # 1. Precalculate fragment data
        frag_data = []
        for f_indices in frags:
            pos_list = []
            radii_list = []
            for idx in f_indices:
                p = conf.GetAtomPosition(idx)
                pos_list.append([p.x, p.y, p.z])
                try:
                    radii_list.append(pt.GetRvdw(rd_mol.GetAtomWithIdx(idx).GetAtomicNum()))
                except Exception:
                    radii_list.append(1.5)

            pos_np = np.array(pos_list)
            radii_np = np.array(radii_list)
            frag_data.append({
                "indices": f_indices,
                "positions": pos_np,
                "radii": radii_np,
                "max_radius": np.max(radii_np) if len(radii_np) > 0 else 1.5
            })

        # Parameters
        scale = 1.1  # Collision threshold scale
        max_iters = 50
        grid_size = 5.0

        for iteration in range(max_iters):
            if check_halted_cb():
                raise WorkerHaltError("Halted")

            moved = False
            grid = {}
            for i, fd in enumerate(frag_data):
                fd_min = np.min(fd["positions"], axis=0)
                fd_max = np.max(fd["positions"], axis=0)
                margin = fd["max_radius"] * scale
                fd["bbox_min"] = fd_min - margin
                fd["bbox_max"] = fd_max + margin
                g_min = (fd["bbox_min"] / grid_size).astype(int)
                g_max = (fd["bbox_max"] / grid_size).astype(int)
                for gx in range(g_min[0], g_max[0] + 1):
                    for gy in range(g_min[1], g_max[1] + 1):
                        for gz in range(g_min[2], g_max[2] + 1):
                            cell = (gx, gy, gz)
                            if cell not in grid:
                                grid[cell] = []
                            grid[cell].append(i)

            all_push_vectors = [np.zeros(3) for _ in range(len(frag_data))]
            processed_pairs = set()

            for cell_indices in grid.values():
                if len(cell_indices) < 2:
                    continue
                for idx_idx_i, i in enumerate(cell_indices):
                    for j in cell_indices[idx_idx_i + 1:]:
                        pair = tuple(sorted((i, j)))
                        if pair in processed_pairs:
                            continue
                        processed_pairs.add(pair)
                        if np.any(frag_data[i]["bbox_min"] > frag_data[j]["bbox_max"]) or \
                           np.any(frag_data[j]["bbox_min"] > frag_data[i]["bbox_max"]):
                            continue
                        fd_i = frag_data[i]
                        fd_j = frag_data[j]
                        push_i = np.zeros(3)
                        push_j = np.zeros(3)
                        collision_count = 0
                        for p_idx_i, p_i in enumerate(fd_i["positions"]):
                            r_i = fd_i["radii"][p_idx_i]
                            for p_idx_j, p_j in enumerate(fd_j["positions"]):
                                r_j = fd_j["radii"][p_idx_j]
                                diff = p_i - p_j
                                dist_sq = np.dot(diff, diff)
                                min_dist = (r_i + r_j) * scale
                                if dist_sq < 0.0001:
                                    diff = np.random.uniform(-0.1, 0.1, 3)
                                    dist_sq = np.dot(diff, diff)
                                if dist_sq < min_dist * min_dist:
                                    dist = np.sqrt(dist_sq)
                                    if dist < 0.0001:
                                        dist = 0.0001
                                    push_mag = (min_dist - dist) / 2.0
                                    vec = (diff / dist) * push_mag
                                    push_i += vec
                                    push_j -= vec
                                    collision_count += 1
                        if collision_count > 0:
                            all_push_vectors[i] += push_i / collision_count
                            all_push_vectors[j] += push_j / collision_count
                            moved = True

            if not moved:
                break

            for i, push in enumerate(all_push_vectors):
                if np.any(push != 0):
                    frag_data[i]["positions"] += push
                    for local_idx, global_idx in enumerate(frag_data[i]["indices"]):
                        conf.SetAtomPosition(global_idx, frag_data[i]["positions"][local_idx].tolist())

        safe_status_cb("Collision avoidance completed.")
    except WorkerHaltError:
        raise
    except Exception as e: # pragma: no cover
        import traceback
        traceback.print_exc()
        safe_status_cb(f"Collision avoidance warning: {e}")


def _iterative_optimize(mol, method, check_halted_cb, safe_status_cb, max_iters=4000, chunk_size=100, options=None):
    """Perform force field optimization in small chunks to avoid UI freezing and allow halts."""
    try:
        # Respect the "Consider Intermolecular Interaction for RDKit" setting
        ignore_interfrag = not options.get("optimize_intermolecular_interaction_rdkit", True) if options else False

        if method in ("MMFF", "MMFF94", "MMFF94S", "MMFF94s"):
            mmff_variant = "MMFF94" if method == "MMFF94" else "MMFF94s"
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmff_variant)
            if props is None:
                raise ValueError(f"Failed to generate MMFF properties for variant {mmff_variant}.")
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0, ignoreInterfragInteractions=ignore_interfrag)
        elif method == "UFF":
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0, ignoreInterfragInteractions=ignore_interfrag)
        else:
            raise ValueError(f"Unknown optimization method: {method}")

        if ff is None:
            raise ValueError(f"Failed to setup force field for method: {method}")

        ff.Initialize()
        iters_done = 0
        while iters_done < max_iters:
            if check_halted_cb():
                raise WorkerHaltError("Halted")
            res = ff.Minimize(maxIts=chunk_size)
            iters_done += chunk_size
            if res == 0:
                break
            import time
            time.sleep(0.001)

        return True
    except WorkerHaltError:
        raise
    except Exception as e:
        safe_status_cb(f"Iterative optimization ({method}) error: {e}")
        return False


def _iterative_optimize_obabel(mol, method, check_halted_cb, safe_status_cb, max_iters=4000, chunk_size=100):
    """Perform force field optimization using OpenBabel in chunks to avoid UI freezing."""
    try:
        if not OBABEL_AVAILABLE or pybel is None:
            raise RuntimeError("OpenBabel is not available.")
            
        ff_name = method.lower()
        
        # Convert RDKit mol to OpenBabel pybel mol
        mol_block = Chem.MolToMolBlock(mol)
        ob_mol = pybel.readstring("mol", mol_block)
        
        # Set up forcefield
        ff = pybel.ob.OBForceField.FindForceField(ff_name)
        if not ff:
            raise ValueError(f"Forcefield '{ff_name}' not found in OpenBabel.")
            
        if not ff.Setup(ob_mol.OBMol):
            error_msg = f"Failed to load parameters for '{ff_name}'."
            if ff_name.lower() == "gaff":
                error_msg += " (Check if OpenBabel 'gaff.dat' and 'gaff.prm' parameter files are missing in your environment.)"
            safe_status_cb(error_msg)
            raise RuntimeError(error_msg)
            
        ff.SteepestDescent(100) # Initial stabilization
        
        iters_done = 0
        while iters_done < max_iters:
            if check_halted_cb():
                raise WorkerHaltError("Halted")
            
            ff.ConjugateGradients(chunk_size)
            iters_done += chunk_size
            
            import time
            time.sleep(0.001)
            
        ff.GetCoordinates(ob_mol.OBMol)
        
        # Copy coordinates back to RDKit mol
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            atom = ob_mol.OBMol.GetAtom(i + 1)
            conf.SetAtomPosition(i, [atom.GetX(), atom.GetY(), atom.GetZ()])
            
        return True
    except WorkerHaltError:
        raise
    except Exception as e:
        safe_status_cb(f"Iterative optimization ({method}) error (OpenBabel): {e}")
        return False


def _parse_explicit_stereo(mol_block):
    """Parse explicit stereochemistry from the MOL block."""
    explicit_stereo = {}
    mol_lines = mol_block.split("\n")
    for line in mol_lines:
        if line.startswith("M  CFG"):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    bond_idx = int(parts[3]) - 1  # MOL format is 1-indexed
                    cfg_value = int(parts[4])
                    # cfg_value: 1=Z, 2=E in MOL format
                    if cfg_value == 1:
                        explicit_stereo[bond_idx] = Chem.BondStereo.STEREOZ
                    elif cfg_value == 2:
                        explicit_stereo[bond_idx] = Chem.BondStereo.STEREOE
                except (ValueError, IndexError):
                    continue
    return explicit_stereo


def _apply_explicit_stereo(mol, explicit_stereo):
    """Apply explicit stereochemistry to the molecule."""
    for bond_idx, stereo_type in explicit_stereo.items():
        if bond_idx < mol.GetNumBonds():
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Find suitable stereo atoms
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()

                # Pick heavy atom neighbors preferentially
                begin_neighbors = [
                    nbr
                    for nbr in begin_atom.GetNeighbors()
                    if nbr.GetIdx() != end_atom.GetIdx()
                ]
                end_neighbors = [
                    nbr
                    for nbr in end_atom.GetNeighbors()
                    if nbr.GetIdx() != begin_atom.GetIdx()
                ]

                if begin_neighbors and end_neighbors:
                    # Prefer heavy atoms
                    begin_heavy = [
                        n for n in begin_neighbors if n.GetAtomicNum() > 1
                    ]
                    end_heavy = [
                        n for n in end_neighbors if n.GetAtomicNum() > 1
                    ]

                    stereo_atom1 = (
                        begin_heavy[0] if begin_heavy else begin_neighbors[0]
                    ).GetIdx()
                    stereo_atom2 = (
                        end_heavy[0] if end_heavy else end_neighbors[0]
                    ).GetIdx()

                    bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                    bond.SetStereo(stereo_type)


def _perform_direct_conversion(mol_block, mol, options, _check_halted, _safe_status):
    """
    Perform direct 3D conversion using 2D coordinates and adding missing H without embedding.
    Returns the processed mol object if successful, or raises Exception on failure.
    """
    # 1) Parse MOL block *with* existing hydrogens (removeHs=False)
    #    to get coordinates for *all existing* atoms.
    parsed_coords = []  # all-atom coordinates (x, y, z)
    stereo_dirs = []  # list of (begin_idx, end_idx, stereo_flag)

    base2d_all = None
    try:
        # Parse including hydrogen atoms
        base2d_all = Chem.MolFromMolBlock(
            mol_block, removeHs=False, sanitize=True
        )
    except Exception:
        try:
            base2d_all = Chem.MolFromMolBlock(
                mol_block, removeHs=False, sanitize=False
            )
        except Exception:
            base2d_all = None

    if base2d_all is not None and base2d_all.GetNumConformers() > 0:
        oconf = base2d_all.GetConformer()
        for i in range(base2d_all.GetNumAtoms()):
            p = oconf.GetAtomPosition(i)
            parsed_coords.append((float(p.x), float(p.y), 0.0))

    # 2) Parse wedge/dash bond information (using all atoms)
    try:
        lines = mol_block.splitlines()
        counts_idx = None

        for i, ln in enumerate(lines[:40]):
            if re.match(r"^\s*\d+\s+\d+", ln):
                counts_idx = i
                break

        if counts_idx is not None:
            parts = lines[counts_idx].split()
            try:
                natoms = int(parts[0])
                nbonds = int(parts[1])
            except Exception:
                natoms = nbonds = 0

            # All-atom map (MOL 1-based index -> 0-based index)
            atom_map = {i + 1: i for i in range(natoms)}

            bond_start = counts_idx + 1 + natoms
            for j in range(
                min(nbonds, max(0, len(lines) - bond_start))
            ):
                bond_line = lines[bond_start + j]
                try:
                    m = re.match(
                        r"^\s*(\d+)\s+(\d+)\s+(\d+)(?:\s+(-?\d+))?",
                        bond_line,
                    )
                    if m:
                        try:
                            atom1_mol = int(
                                m.group(1)
                            )  # 1-based MOL index
                            atom2_mol = int(
                                m.group(2)
                            )  # 1-based MOL index
                        except Exception:
                            continue
                        try:
                            stereo_raw = (
                                int(m.group(4))
                                if m.group(4) is not None
                                else 0
                            )
                        except Exception:
                            stereo_raw = 0
                    else:
                        fields = bond_line.split()
                        if len(fields) >= 4:
                            try:
                                atom1_mol = int(
                                    fields[0]
                                )  # 1-based MOL index
                                atom2_mol = int(
                                    fields[1]
                                )  # 1-based MOL index
                            except Exception:
                                continue
                            try:
                                stereo_raw = (
                                    int(fields[3])
                                    if len(fields) > 3
                                    else 0
                                )
                            except Exception:
                                stereo_raw = 0
                        else:
                            continue

                    # Normalize V2000 stereo notation
                    if stereo_raw == 1:
                        stereo_flag = 1  # Wedge
                    elif stereo_raw == 2:
                        stereo_flag = 6  # Dash (6 is Dash in V2000)
                    else:
                        stereo_flag = stereo_raw

                    # Check using all-atom map
                    if atom1_mol in atom_map and atom2_mol in atom_map:
                        idx1 = atom_map[atom1_mol]
                        idx2 = atom_map[atom2_mol]
                        if stereo_flag in (
                            1,
                            6,
                        ):  # Wedge (1) or Dash (6)
                            stereo_dirs.append(
                                (idx1, idx2, stereo_flag)
                            )
                except Exception:
                    continue
    except Exception:
        stereo_dirs = []

    # Fallback for parsed_coords (if RDKit parse failed)
    if not parsed_coords:
        try:
            lines = mol_block.splitlines()
            counts_idx = None
            for i, ln in enumerate(lines[:40]):
                if re.match(r"^\s*\d+\s+\d+", ln):
                    counts_idx = i
                    break
            if counts_idx is not None:
                parts = lines[counts_idx].split()
                try:
                    natoms = int(parts[0])
                except Exception:
                    natoms = 0
                atom_start = counts_idx + 1
                for j in range(
                    min(natoms, max(0, len(lines) - atom_start))
                ):
                    atom_line = lines[atom_start + j]
                    try:
                        x = float(atom_line[0:10].strip())
                        y = float(atom_line[10:20].strip())
                        z = float(atom_line[20:30].strip())
                    except Exception:
                        fields = atom_line.split()
                        if len(fields) >= 4:
                            try:
                                x = float(fields[0])
                                y = float(fields[1])
                                z = float(fields[2])
                            except Exception:
                                continue
                        else:
                            continue
                    # Do not skip hydrogen atoms
                    parsed_coords.append((x, y, z))
        except Exception:
            parsed_coords = []

    if not parsed_coords:
        raise ValueError(
            "Failed to parse coordinates from MOL block for direct conversion."
        )

    # 3) mol is already in AddHs state
    #    Get original atom count (including H) from parsed_coords length
    num_existing_atoms = len(parsed_coords)

    # 4) Create conformer
    conf = Chem.Conformer(mol.GetNumAtoms())

    for i in range(mol.GetNumAtoms()):
        if i < num_existing_atoms:
            # Existing atoms (including H): set 2D coordinates (z=0)
            x, y, z_ignored = parsed_coords[i]
            try:
                conf.SetAtomPosition(
                    i, rdGeometry.Point3D(float(x), float(y), 0.0)
                )
            except Exception:  # pragma: no cover
                import traceback
                traceback.print_exc()
        else:
            # Newly added H atoms: place near the parent atom
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() == 1:
                neighs = [
                    n
                    for n in atom.GetNeighbors()
                    if n.GetIdx() < num_existing_atoms
                ]
                heavy_pos_found = False
                for nb in neighs:  # Parent atom (heavy atom or existing H)
                    try:
                        nb_idx = nb.GetIdx()
                        # Already filtered in neighs
                        nbpos = conf.GetAtomPosition(nb_idx)
                        # Geometry-based placement:
                        # Compute an "empty" direction around the parent atom by
                        # summing existing bond unit vectors and taking the
                        # opposite. If degenerate, pick a perpendicular or
                        # fallback vector. Rotate slightly if multiple Hs already
                        # attached to avoid overlap.
                        parent_idx = nb_idx
                        try:
                            parent_pos = conf.GetAtomPosition(
                                parent_idx
                            )
                            parent_atom = mol.GetAtomWithIdx(parent_idx)
                            # collect unit vectors to already-placed neighbors (idx < i)
                            vecs = []
                            for nbr in parent_atom.GetNeighbors():
                                nidx = nbr.GetIdx()
                                if nidx == i:
                                    continue
                                # only consider neighbors whose positions are already set
                                if nidx < i:
                                    try:
                                        p = conf.GetAtomPosition(nidx)
                                        vx = float(p.x) - float(
                                            parent_pos.x
                                        )
                                        vy = float(p.y) - float(
                                            parent_pos.y
                                        )
                                        nrm = math.hypot(vx, vy)
                                        if nrm > 1e-6:
                                            vecs.append(
                                                (vx / nrm, vy / nrm)
                                            )
                                    except Exception:
                                        continue

                            if vecs:
                                sx = sum(v[0] for v in vecs)
                                sy = sum(v[1] for v in vecs)
                                fx = -sx
                                fy = -sy
                                fn = math.hypot(fx, fy)
                                if fn < 1e-6:
                                    # degenerate: pick a perpendicular to first bond
                                    fx = -vecs[0][1]
                                    fy = vecs[0][0]
                                    fn = math.hypot(fx, fy)
                                fx /= fn
                                fy /= fn

                                # Avoid placing multiple Hs at identical directions
                                existing_h_count = sum(
                                    1
                                    for nbr in parent_atom.GetNeighbors()
                                    if nbr.GetIdx() < i
                                    and nbr.GetAtomicNum() == 1
                                )
                                angle = existing_h_count * (
                                    math.pi / 6.0
                                )  # 30deg steps
                                cos_a = math.cos(angle)
                                sin_a = math.sin(angle)
                                rx = fx * cos_a - fy * sin_a
                                ry = fx * sin_a + fy * cos_a

                                bond_length = 1.0
                                conf.SetAtomPosition(
                                    i,
                                    rdGeometry.Point3D(
                                        float(parent_pos.x)
                                        + rx * bond_length,
                                        float(parent_pos.y)
                                        + ry * bond_length,
                                        0.3,
                                    ),
                                )
                            else:
                                # No existing placed neighbors: fallback to small offset
                                conf.SetAtomPosition(
                                    i,
                                    rdGeometry.Point3D(
                                        float(parent_pos.x) + 0.5,
                                        float(parent_pos.y) + 0.5,
                                        0.3,
                                    ),
                                )

                            heavy_pos_found = True
                            break
                        except Exception:
                            # fall back to trying the next neighbor if any
                            continue
                    except Exception:
                        continue
                if not heavy_pos_found:
                    # Fallback (near origin)
                    try:
                        conf.SetAtomPosition(
                            i, rdGeometry.Point3D(0.0, 0.0, 0.10)
                        )
                    except Exception:  # pragma: no cover
                        import traceback
                        traceback.print_exc()
    # 5) Apply Z-offset for Wedge/Dash constraints
    try:
        stereo_z_offset = 1.5  # wedge -> +1.5, dash -> -1.5
        for begin_idx, end_idx, stereo_flag in stereo_dirs:
            try:
                # Indices should be within existing atoms
                if (
                    begin_idx >= num_existing_atoms
                    or end_idx >= num_existing_atoms
                ):
                    continue

                if stereo_flag not in (1, 6):
                    continue

                sign = 1.0 if stereo_flag == 1 else -1.0
                
                # Apply Z-offset to end_idx (the atom at the end of the stereo notation)
                pos = conf.GetAtomPosition(end_idx)
                newz = float(pos.z) + (
                    stereo_z_offset * sign
                )  # Add offset to the existing Z-coordinate (usually 0)
                conf.SetAtomPosition(
                    end_idx,
                    rdGeometry.Point3D(
                        float(pos.x), float(pos.y), float(newz)
                    ),
                )
            except Exception:
                continue
    except Exception:  # pragma: no cover
        import traceback
        traceback.print_exc()
    # Replace conformer and finish
    try:
        mol.RemoveAllConformers()
    except Exception:  # pragma: no cover
        import traceback
        traceback.print_exc()
    mol.AddConformer(conf, assignId=True)

    # Optimization (respects do_optimize flag)
    do_optimize = options.get("do_optimize", True) if options else True
    if do_optimize:
        # Collision avoidance before optimization
        _adjust_collision_avoidance(mol, _check_halted, _safe_status)
        _safe_status("Direct conversion: optimizing geometry...")
        if _check_halted():
            raise WorkerHaltError("Halted")

        # Determine backend and actual method
        opt_method = options.get("optimization_method", "MMFF94s_RDKIT") if options else "MMFF94s_RDKIT"
        backend = "OBABEL" if "OBABEL" in opt_method.upper() else "RDKIT"
        
        # Default to MMFF94s
        method_key = "MMFF94s"
        if "MMFF94" in opt_method.upper():
            method_key = "MMFF94" if "MMFF94S" not in opt_method.upper() else "MMFF94s"
        elif "UFF" in opt_method.upper():
            method_key = "UFF"
        elif "GAFF" in opt_method.upper():
            method_key = "GAFF"
        elif "GHEMICAL" in opt_method.upper():
            method_key = "GHEMICAL"

        # Apply Force Field Optimization
        if backend == "OBABEL":
            try:
                # Add property for UI feedback
                mol.SetProp("_pme_optimization_method", opt_method)
            except Exception:
                pass
            _safe_status(f"Applying force field optimization ({method_key} / OpenBabel)...")
            opt_success = _iterative_optimize_obabel(mol, method_key, _check_halted, _safe_status)
            if not opt_success:
                _safe_status(f"Warning: Optimization with {opt_method} failed. Using unoptimized structure.")
                try:
                    mol.ClearProp("_pme_optimization_method")
                except Exception:
                    pass
        else: # RDKit backend
            try:
                # Add property for UI feedback
                mol.SetProp("_pme_optimization_method", opt_method)
            except Exception:
                pass
            _safe_status(f"Applying force field optimization ({method_key} / RDKit)...")
            opt_success = _iterative_optimize(mol, method_key, _check_halted, _safe_status, options=options)
            if not opt_success:
                _safe_status(f"Warning: Optimization with {opt_method} failed. Using unoptimized structure.")
                try:
                    mol.ClearProp("_pme_optimization_method")
                except Exception:
                    pass

    if _check_halted():
        raise WorkerHaltError("Halted")
    return mol


def _perform_optimize_only(mol, options, worker_id, _check_halted, _safe_status, _safe_finished):
    """Perform optimization on an existing 3D structure and emit finished signal."""
    _safe_status("Optimizing existing 3D structure...")
    opt_method = str(options.get("optimization_method", "MMFF_RDKIT")).upper()
    
    # Determine backend and actual method
    backend = "OBABEL" if "OBABEL" in opt_method else "RDKIT"
    
    # Default to MMFF94s
    method_key = "MMFF94s"
    if "MMFF94" in opt_method:
        method_key = "MMFF94" if "MMFF94S" not in opt_method else "MMFF94s"
    elif "UFF" in opt_method:
        method_key = "UFF"
    elif "GAFF" in opt_method:
        method_key = "GAFF"
    elif "GHEMICAL" in opt_method:
        method_key = "GHEMICAL"

    if backend == "OBABEL":
        try:
            mol.SetProp("_pme_optimization_method", opt_method)
        except Exception:
            pass
        opt_success = _iterative_optimize_obabel(mol, method_key, _check_halted, _safe_status)
    else:
        try:
            mol.SetProp("_pme_optimization_method", opt_method)
        except Exception:
            pass
        opt_success = _iterative_optimize(mol, method_key, _check_halted, _safe_status, options=options)
    
    if not opt_success:
        raise Exception(f"Optimization with {opt_method} failed. You can change the method in the settings.")
    
    if _check_halted():
        raise WorkerHaltError("Halted")
    try:
        _safe_finished((worker_id, mol))
    except Exception:
        _safe_finished(mol)
    friendly = _OPT_METHOD_LABELS.get(opt_method.upper(), opt_method)
    _safe_status(f"Optimization completed ({friendly}).")


def _perform_obabel_conversion(mol_block, conversion_mode, opt_method, worker_id, options, _check_halted, _safe_status, _safe_finished):
    """
    Perform Open Babel 3D conversion and optimization.
    Returns True if successful, False if it should fall back.
    Raises exception on terminal failure.
    """
    _safe_status(
        "RDKit embedding failed or disabled. Attempting Open Babel..."
    )
    try:
        if not OBABEL_AVAILABLE:
            raise RuntimeError(
                "Open Babel (pybel) is not available in this Python environment."
            )
        ob_mol = pybel.readstring("mol", mol_block)
        try:
            ob_mol.addh()
        except Exception:  # pragma: no cover
            import traceback
            traceback.print_exc()
        ob_mol.make3D()
        # Skip OpenBabel's redundant manual localopt here.
        # We will optimize rd_mol immediately after conversion using the selected user method.
        molblock_ob = ob_mol.write("mol")
        rd_mol = Chem.MolFromMolBlock(molblock_ob, removeHs=False)
        if rd_mol is None:
            raise ValueError("Open Babel produced invalid MOL block.")
        rd_mol = Chem.AddHs(rd_mol)
        # Collision avoidance before optimization
        _adjust_collision_avoidance(rd_mol, _check_halted, _safe_status)
        try:
            opt_method_raw = opt_method or "MMFF94s_RDKIT"
            backend = "OBABEL" if "OBABEL" in opt_method_raw.upper() else "RDKIT"

            method_key = "MMFF94s"
            if "MMFF94" in opt_method_raw.upper():
                method_key = "MMFF94" if "MMFF94S" not in opt_method_raw.upper() else "MMFF94s"
            elif "UFF" in opt_method_raw.upper():
                method_key = "UFF"
            elif "GAFF" in opt_method_raw.upper():
                method_key = "GAFF"
            elif "GHEMICAL" in opt_method_raw.upper():
                method_key = "GHEMICAL"

            if backend == "OBABEL":
                try:
                    rd_mol.SetProp("_pme_optimization_method", opt_method_raw)
                except Exception:
                    pass
                _safe_status(f"Optimizing with OpenBabel ({method_key})...")
                opt_success = _iterative_optimize_obabel(rd_mol, method_key, _check_halted, _safe_status)
            else:
                try:
                    rd_mol.SetProp("_pme_optimization_method", opt_method_raw)
                except Exception:
                    pass
                _safe_status(f"Optimizing with RDKit ({method_key})...")
                opt_success = _iterative_optimize(rd_mol, method_key, _check_halted, _safe_status, options=options)

            if not opt_success:
                _safe_status(f"Warning: Optimization with {opt_method_raw} failed. Using unoptimized structure.")
                try:
                    rd_mol.ClearProp("_pme_optimization_method")
                except Exception:
                    pass
        except WorkerHaltError:
            raise
        except Exception as opt_err:
            _safe_status(f"Warning: Optimization failed: {opt_err}. Using unoptimized structure.")
            try:
                rd_mol.ClearProp("_pme_optimization_method")
            except Exception:
                pass

        _safe_status(
            "Open Babel embedding succeeded. Warning: Conformation accuracy may be limited."
        )
        # CRITICAL: Check for halt *before* emitting finished signal
        if _check_halted():
            raise WorkerHaltError("Halted")

        try:
            _safe_finished((worker_id, rd_mol))
        except Exception:
            _safe_finished(rd_mol)
        return True
    except WorkerHaltError:
        raise
    except Exception as ob_err:
        if conversion_mode == "obabel":
            # obabel-only mode: no further fallback
            raise RuntimeError(f"Open Babel 3D conversion failed: {ob_err}")
        # fallback mode: continue to direct conversion below
        _safe_status(
            f"Open Babel unavailable or failed ({ob_err}). "
            "Falling back to direct conversion..."
        )
    return False


class CalculationWorker(QObject):
    status_update = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(object)
    start_work = pyqtSignal(str, object)

    def __init__(self, parent=None): # pragma: no cover
        super().__init__(parent)
        try:
            self.start_work.connect(self.run_calculation)
        except Exception:
            import traceback
            traceback.print_exc()

    @pyqtSlot(str, object)
    def run_calculation(self, mol_block, options=None):
        worker_id = None

        # The worker may be asked to halt via a shared set `halt_ids` and
        # identifies its own run by options['worker_id'] (int).
        def _check_halted():
            try:
                halt_ids = getattr(self, "halt_ids", None)
                if worker_id is None:
                    if getattr(self, "halt_all", False):
                        return True
                    if halt_ids is None:
                        return False
                    # Support both None-in-set and string sentinel for compatibility
                    return (None in halt_ids) or ("ALL" in halt_ids)

                if halt_ids is None:
                    return False
                return worker_id in halt_ids
            except Exception:
                return False

        # Safe-emission helpers: do nothing if this worker has been halted.
        def _safe_status(msg): # pragma: no cover
            try:
                if _check_halted():
                    raise WorkerHaltError("Halted")
                self.status_update.emit(msg)
            except WorkerHaltError:
                raise
            except Exception:
                # Swallow any signal-emission errors to avoid crashing the worker
                pass

        def _safe_finished(payload):  # pragma: no cover
            try:
                if _check_halted():
                    raise WorkerHaltError("Halted")
                try:
                    self.finished.emit(payload)
                except WorkerHaltError:
                    raise
                except TypeError:
                    try:
                        if isinstance(payload, (list, tuple)) and len(payload) >= 2:
                            self.finished.emit(payload[1])
                        else:
                            self.finished.emit(payload)
                    except WorkerHaltError:
                        raise
                    except Exception:  # pragma: no cover
                        import traceback
                        traceback.print_exc()
            except WorkerHaltError:
                raise
            except Exception:  # pragma: no cover
                import traceback
                traceback.print_exc()

        def _safe_error(msg):  # pragma: no cover
            try:
                if msg != "Halted" and _check_halted():
                    raise WorkerHaltError("Halted")
                try:
                    self.error.emit((worker_id, msg))
                except WorkerHaltError:
                    raise
                except Exception:
                    try:
                        self.error.emit(msg)
                    except WorkerHaltError:
                        raise
                    except Exception:  # pragma: no cover
                        import traceback
                        traceback.print_exc()
            except WorkerHaltError:
                raise
            except Exception:  # pragma: no cover
                import traceback
                traceback.print_exc()

        try:
            worker_id = None
            try:
                worker_id = options.get("worker_id") if options else None
            except Exception:
                worker_id = None

            _warned_no_worker_id = False
            if worker_id is None: # pragma: no cover
                try:
                    # best-effort, swallow any errors (signals may not be connected)
                    self.status_update.emit(
                        "Warning: worker started without 'worker_id'; will listen for global halt signals."
                    )
                except Exception:
                    import traceback
                    traceback.print_exc()
                _warned_no_worker_id = True

            # ── 1. Validate input ──────────────────────────────────────────
            if options is None:
                options = {}
            conversion_mode = options.get("conversion_mode", "fallback")
            params = None
            if not mol_block:
                raise ValueError("No atoms to convert.")

            _safe_status("Creating 3D structure...")

            # ── 2. Parse molecule ──────────────────────────────────────────

            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None:
                raise ValueError("Failed to create molecule from MOL block.")

            # Check early whether this run has been requested to halt
            if _check_halted():
                raise WorkerHaltError("Halted")

            # ── 3. Parse and apply explicit stereochemistry ────────────────
            explicit_stereo = _parse_explicit_stereo(mol_block)
            _apply_explicit_stereo(mol, explicit_stereo)
            # Do NOT call AssignStereochemistry here — it overrides our explicit labels

            if _check_halted():
                raise WorkerHaltError("Halted")

            # ── 4. Mode: optimize existing 3D structure ────────────────────
            if conversion_mode == "optimize_only":
                _perform_optimize_only(mol, options, worker_id, _check_halted, _safe_status, _safe_finished)
                return

            # ── 5. Add hydrogens + re-apply stereo ────────────────────────
            mol = Chem.AddHs(mol)

            if _check_halted():
                raise WorkerHaltError("Halted")

            # Re-apply after AddHs which may renumber atoms
            _apply_explicit_stereo(mol, explicit_stereo)

            # ── 6. Mode: direct 2D→3D without embedding ────────────────────
            if conversion_mode == "direct":
                _safe_status(
                    "Direct conversion: using 2D coordinates + adding missing H (no embedding)."
                )
                try:
                    mol = _perform_direct_conversion(mol_block, mol, options, _check_halted, _safe_status)

                    try:
                        _safe_finished((worker_id, mol))
                    except WorkerHaltError:
                        raise
                    except Exception:
                        _safe_finished(mol)
                    _safe_status("Direct conversion completed.")
                    return
                except WorkerHaltError:
                    raise
                except Exception as e:
                    _safe_status(f"Direct conversion failed: {e}")

            # ── 7. RDKit ETKDG embedding ───────────────────────────────────
            params = AllChem.ETKDGv2()
            params.randomSeed = 42
            params.useExpTorsionAnglePrefs = True
            params.useBasicKnowledge = True
            params.enforceChirality = True  # critical for stereo preservation

            # Snapshot stereo before embedding so we can restore it afterwards
            original_stereo_info = []
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        stereo_atoms = bond.GetStereoAtoms()
                        original_stereo_info.append(
                            (bond.GetIdx(), stereo_type, stereo_atoms)
                        )

            # Also snapshot stereo bonds not in explicit_stereo
            for bond in mol.GetBonds():
                if (
                    bond.GetBondType() == Chem.BondType.DOUBLE
                    and bond.GetStereo() != Chem.BondStereo.STEREONONE
                    and bond.GetIdx() not in explicit_stereo
                ):
                    stereo_atoms = bond.GetStereoAtoms()
                    original_stereo_info.append(
                        (bond.GetIdx(), bond.GetStereo(), stereo_atoms)
                    )

            if conversion_mode in ("fallback", "rdkit"):
                _safe_status("RDKit: Embedding 3D coordinates...")

            if _check_halted():
                raise WorkerHaltError("Halted")

            conf_id = -1

            # First attempt: Standard ETKDG with stereo enforcement
            try:
                # Only attempt RDKit embedding if mode allows
                if conversion_mode in ("fallback", "rdkit"):
                    conf_id = AllChem.EmbedMolecule(mol, params)
                else:
                    conf_id = -1
                # Final check before returning success
                if _check_halted():
                    raise WorkerHaltError("Halted")
            except WorkerHaltError:
                raise
            except Exception as e:
                # Standard embedding failed; report and continue to fallback attempts
                _safe_status(f"Standard embedding failed: {e}")

                # Second attempt: Use constraint embedding if available (only when RDKit is allowed)
                if conf_id == -1 and conversion_mode in ("fallback", "rdkit"):
                    try:
                        # Create distance constraints for double bonds to enforce E/Z geometry
                        bounds_matrix = AllChem.GetMoleculeBoundsMatrix(mol)

                        # Add constraints for E/Z bonds
                        for bond_idx, stereo, stereo_atoms in original_stereo_info:
                            bond = mol.GetBondWithIdx(bond_idx)
                            if len(stereo_atoms) == 2:
                                atom1_idx = bond.GetBeginAtomIdx()
                                atom2_idx = bond.GetEndAtomIdx()
                                neighbor1_idx = stereo_atoms[0]
                                neighbor2_idx = stereo_atoms[1]

                                # For Z (cis): neighbors should be closer
                                # For E (trans): neighbors should be farther
                                if stereo == Chem.BondStereo.STEREOZ:
                                    # Z configuration: set shorter distance constraint
                                    target_dist = 3.0  # Angstroms
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = min(
                                        bounds_matrix[neighbor1_idx][neighbor2_idx],
                                        target_dist,
                                    )
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = min(
                                        bounds_matrix[neighbor2_idx][neighbor1_idx],
                                        target_dist,
                                    )
                                elif stereo == Chem.BondStereo.STEREOE:
                                    # E configuration: set longer distance constraint
                                    target_dist = 5.0  # Angstroms
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = max(
                                        bounds_matrix[neighbor1_idx][neighbor2_idx],
                                        target_dist,
                                    )
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = max(
                                        bounds_matrix[neighbor2_idx][neighbor1_idx],
                                        target_dist,
                                    )

                        DoTriangleSmoothing(bounds_matrix)
                        conf_id = AllChem.EmbedMolecule(mol, bounds_matrix, params)
                        _safe_status("Constraint-based embedding succeeded")
                    except Exception:
                        # Constraint embedding failed: only raise error if mode is 'rdkit', otherwise allow fallback
                        _safe_status("RDKit: Constraint embedding failed")
                        if conversion_mode == "rdkit":
                            raise RuntimeError("RDKit: Constraint embedding failed")
                        conf_id = -1

            # Attempt 3: basic ETKDGv2 (no stereo constraints)
            if conf_id == -1:
                try:
                    if conversion_mode in ("fallback", "rdkit"):
                        basic_params = AllChem.ETKDGv2()
                        basic_params.randomSeed = 42
                        conf_id = AllChem.EmbedMolecule(mol, basic_params)
                    else:
                        conf_id = -1
                except Exception:  # pragma: no cover
                    import traceback
                    traceback.print_exc()

            # ── 8. RDKit succeeded: optimize + emit ────────────────────────
            opt_method = None
            try:
                opt_method = options.get("optimization_method") if options else None
            except Exception:
                opt_method = None

            if conf_id != -1:
                # Success with RDKit: optimize and finish
                # CRITICAL: Restore original stereochemistry after embedding (explicit labels first)
                for bond_idx, stereo, stereo_atoms in original_stereo_info:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if len(stereo_atoms) == 2:
                        bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                    bond.SetStereo(stereo)

                # Collision avoidance before optimization
                _adjust_collision_avoidance(mol, _check_halted, _safe_status)

                try:
                    opt_method_raw = opt_method or "MMFF94s_RDKIT"
                    backend = "OBABEL" if "OBABEL" in opt_method_raw.upper() else "RDKIT"
                    
                    method_key = "MMFF94s"
                    if "MMFF94" in opt_method_raw.upper():
                        method_key = "MMFF94" if "MMFF94S" not in opt_method_raw.upper() else "MMFF94s"
                    elif "UFF" in opt_method_raw.upper():
                        method_key = "UFF"
                    elif "GAFF" in opt_method_raw.upper():
                        method_key = "GAFF"
                    elif "GHEMICAL" in opt_method_raw.upper():
                        method_key = "GHEMICAL"

                    if backend == "OBABEL":
                        try:
                            mol.SetProp("_pme_optimization_method", opt_method_raw)
                        except Exception:  # pragma: no cover
                            pass
                        _safe_status(f"Optimizing with OpenBabel ({method_key})...")
                        opt_success = _iterative_optimize_obabel(mol, method_key, _check_halted, _safe_status)
                    else:
                        try:
                            mol.SetProp("_pme_optimization_method", opt_method_raw)
                        except Exception:  # pragma: no cover
                            pass
                        _safe_status(f"Optimizing with RDKit ({method_key})...")
                        opt_success = _iterative_optimize(mol, method_key, _check_halted, _safe_status, options=options)

                    if not opt_success:
                        raise Exception(f"Optimization with {opt_method_raw} failed. You can change the method in the settings.")
                except WorkerHaltError:
                    raise
                except Exception as opt_err:
                    if conversion_mode == "rdkit":
                        raise Exception(f"Optimization with {opt_method_raw} failed: {opt_err}")
                    _safe_status(f"Optimization failed: {opt_err}. Falling back...")  # pragma: no cover
                    try:  # pragma: no cover
                        mol.ClearProp("_pme_optimization_method")
                    except Exception:  # pragma: no cover
                        pass
                    # Allow fallback to proceed instead of crashing
                    conf_id = -1  # pragma: no cover
                
                # CRITICAL: Restore stereochemistry again after optimization (explicit labels priority)
                if conf_id != -1:
                    for bond_idx, stereo, stereo_atoms in original_stereo_info:
                        bond = mol.GetBondWithIdx(bond_idx)
                        if len(stereo_atoms) == 2:
                            bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                        bond.SetStereo(stereo)

                    # CRITICAL: Check for halt *before* emitting finished signal
                    if _check_halted():
                        raise WorkerHaltError("Halted")

                    try:
                        _safe_finished((worker_id, mol))
                    except WorkerHaltError:
                        raise
                    except Exception:
                        _safe_finished(mol)
                    _safe_status("RDKit 3D conversion succeeded.")
                    return

            # ── 9. RDKit failed: try Open Babel ───────────────────────────
            if conf_id == -1 and conversion_mode in ("fallback", "obabel"):
                success = _perform_obabel_conversion(mol_block, conversion_mode, opt_method, worker_id, options, _check_halted, _safe_status, _safe_finished)
                if success:
                    return

            if conf_id == -1 and conversion_mode == "rdkit":
                raise RuntimeError("RDKit 3D conversion failed (rdkit-only mode)")

            # ── 10. Last resort: direct conversion ─────────────────────────
            if conf_id == -1 and conversion_mode == "fallback":
                _safe_status(
                    "All embedding methods failed. Using direct conversion as last resort..."
                )
                direct_opts = dict(options) if options else {}
                direct_opts["conversion_mode"] = "direct"
                self.run_calculation(mol_block, direct_opts)
                return

        except WorkerHaltError:
            _safe_error("Halted")
            return
        except Exception as e:
            _safe_error(str(e))
