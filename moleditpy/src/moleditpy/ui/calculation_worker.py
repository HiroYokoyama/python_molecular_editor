#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import contextlib
import math
import re
import numpy as np
import sys
import subprocess
import time
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

from PyQt6.QtCore import QObject, pyqtSignal, pyqtSlot

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem, rdGeometry
from rdkit.DistanceGeometry import DoTriangleSmoothing

try:
    from .. import OBABEL_AVAILABLE
except ImportError:
    from moleditpy import OBABEL_AVAILABLE

# Only import pybel on demand
if OBABEL_AVAILABLE:
    # Suppress potential import errors if Open Babel is not correctly installed or configured
    with contextlib.suppress(ImportError, OSError, RuntimeError):
        import os
        import glob
        import openbabel
        from openbabel import pybel

        ob_base = os.path.dirname(openbabel.__file__)
        data_candidates = [
            os.path.join(ob_base, "bin", "data"),
            os.path.join(ob_base, "share", "openbabel", "*"),
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
    if not pybel:
        OBABEL_AVAILABLE = False
else:
    pybel = None


class WorkerHaltError(Exception):
    """Custom exception raised when a calculation worker is requested to halt."""
    pass


# Human-readable labels for optimization method keys (mirrors opt3d_method_labels in the UI layer)
_OPT_METHOD_LABELS = {
    "MMFF_RDKIT": "MMFF94s (RDKit)",
    "MMFF94_RDKIT": "MMFF94 (RDKit)",
    "UFF_RDKIT": "UFF (RDKit)",
    "MMFF94S_OBABEL": "MMFF94s (Open Babel)",
    "MMFF94_OBABEL": "MMFF94 (Open Babel)",
    "UFF_OBABEL": "UFF (Open Babel)",
    "GAFF_OBABEL": "GAFF (Open Babel)",
    "GHEMICAL_OBABEL": "Ghemical (Open Babel)",
}


def _adjust_collision_avoidance(
    rd_mol: Any,
    check_halted_cb: Callable[[], bool],
    safe_status_cb: Callable[[str], None],
) -> None:
    """Optimized collision avoidance using spatial partitioning (grid-based)."""
    try:
        frags = Chem.GetMolFrags(rd_mol, asMols=False, sanitizeFrags=False)
        if len(frags) <= 1:
            return

        safe_status_cb(
            f"Resolving potential collisions among {len(frags)} fragments..."
        )
        conf = rd_mol.GetConformer()
        pt = Chem.GetPeriodicTable()

        frag_data = []
        for f_indices in frags:
            pos_list, radii_list = [], []
            for idx in f_indices:
                p = conf.GetAtomPosition(idx)
                pos_list.append([p.x, p.y, p.z])
                atomic_num = rd_mol.GetAtomWithIdx(idx).GetAtomicNum()
                radii_list.append(pt.GetRvdw(atomic_num) if atomic_num > 0 else 1.5)

            pos_np = np.array(pos_list)
            radii_np = np.array(radii_list)
            frag_data.append(
                {
                    "indices": f_indices,
                    "positions": pos_np,
                    "radii": radii_np,
                    "max_radius": np.max(radii_np) if len(radii_np) > 0 else 1.5,
                }
            )

        scale, max_iters, grid_size = 1.1, 50, 5.0
        for _ in range(max_iters):
            if check_halted_cb():
                raise WorkerHaltError("Halted")

            moved: bool = False
            grid: Dict[Any, Any] = {}
            for i, fd in enumerate(frag_data):
                fd_min, fd_max = (
                    np.min(fd["positions"], axis=0),
                    np.max(fd["positions"], axis=0),
                )
                margin = fd["max_radius"] * scale
                fd["bbox_min"], fd["bbox_max"] = fd_min - margin, fd_max + margin
                g_min, g_max = (
                    (fd["bbox_min"] / grid_size).astype(int),
                    (fd["bbox_max"] / grid_size).astype(int),
                )
                for gx in range(g_min[0], g_max[0] + 1):
                    for gy in range(g_min[1], g_max[1] + 1):
                        for gz in range(g_min[2], g_max[2] + 1):
                            grid.setdefault((gx, gy, gz), []).append(i)

            all_push_vectors = [np.zeros(3) for _ in range(len(frag_data))]
            processed_pairs = set()

            for cell_indices in grid.values():
                if len(cell_indices) < 2:
                    continue
                for idx_idx_i, i in enumerate(cell_indices):
                    for j in cell_indices[idx_idx_i + 1 :]:
                        pair = tuple(sorted((i, j)))
                        if pair in processed_pairs:
                            continue
                        processed_pairs.add(pair)
                        if np.any(
                            frag_data[i]["bbox_min"] > frag_data[j]["bbox_max"]
                        ) or np.any(
                            frag_data[j]["bbox_min"] > frag_data[i]["bbox_max"]
                        ):
                            continue
                        fd_i, fd_j = frag_data[i], frag_data[j]
                        push_i, push_j, collision_count = np.zeros(3), np.zeros(3), 0
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
                                    dist = np.sqrt(dist_sq) or 0.0001
                                    vec = (diff / dist) * ((min_dist - dist) / 2.0)
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
                if np.any(push):
                    frag_data[i]["positions"] += push
                    for p_idx, global_idx in enumerate(frag_data[i]["indices"]):
                        conf.SetAtomPosition(
                            global_idx, frag_data[i]["positions"][p_idx]
                        )

        safe_status_cb("Collision avoidance completed.")
    except WorkerHaltError:
        raise
    except (AttributeError, RuntimeError, TypeError, ValueError):
        # Suppress non-critical errors during fragment collision resolution if state is inconsistent
        pass


def _iterative_optimize(
    mol: Any,
    method: str,
    check_halted_cb: Callable[[], bool],
    safe_status_cb: Callable[[str], None],
    max_iters: int = 4000,
    chunk_size: int = 100,
    options: Optional[Dict[str, Any]] = None,
) -> bool:
    """Iteratively optimize the molecule in chunks to allow for responsive halting."""
    safe_status_cb(f"Optimizing with {method}...")
    try:
        ignore_interfrag = not (options or {}).get(
            "optimize_intermolecular_interaction_rdkit", True
        )
        ff = None
        if method == "UFF":
            ff = AllChem.UFFGetMoleculeForceField(
                mol, confId=0, ignoreInterfragInteractions=ignore_interfrag
            )
        elif "MMFF94" in method:
            mmff_variant = "MMFF94" if method == "MMFF94" else "MMFF94s"
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmff_variant)
            if props:
                ff = AllChem.MMFFGetMoleculeForceField(
                    mol, props, confId=0, ignoreInterfragInteractions=ignore_interfrag
                )
            else:
                safe_status_cb(
                    "MMFF94 skipped (unsupported atoms). Reverting to UFF..."
                )
                ff = AllChem.UFFGetMoleculeForceField(
                    mol, confId=0, ignoreInterfragInteractions=ignore_interfrag
                )
                method = "UFF (Fallback)"

        if not ff:
            safe_status_cb(f"Failed to setup {method} Force Field.")
            return False

        iters_done = 0
        while iters_done < max_iters:
            if check_halted_cb():
                raise WorkerHaltError("Halted")
            chunk = min(chunk_size, max_iters - iters_done)
            res = ff.Minimize(maxIts=chunk)
            iters_done += chunk
            time.sleep(0.001)

            if res == 0:
                break
        return True
    except WorkerHaltError:
        raise
    except (RuntimeError, ValueError, TypeError) as e:
        safe_status_cb(f"Iterative optimization ({method}) error: {e}")
        return False


def _iterative_optimize_obabel(
    mol: Any,
    method: str,
    check_halted_cb: Callable[[], bool],
    safe_status_cb: Callable[[str], None],
    max_iters: int = 4000,
    chunk_size: int = 100,
    options: Optional[Dict[str, Any]] = None,
) -> bool:
    """Perform force field optimization using OpenBabel in chunks."""
    try:
        if not OBABEL_AVAILABLE or not pybel:
            raise RuntimeError("OpenBabel is not available.")
        ff = pybel.ob.OBForceField.FindForceField(method.lower())
        if not ff:
            raise ValueError(f"Forcefield '{method}' not found.")

        ob_mol = pybel.readstring("mol", Chem.MolToMolBlock(mol))
        if not ff.Setup(ob_mol.OBMol):
            raise RuntimeError(f"Failed to load parameters for '{method}'.")

        ff.SteepestDescent(100)
        for _ in range(0, max_iters, chunk_size):
            if check_halted_cb():
                raise WorkerHaltError("Halted")
            ff.ConjugateGradients(chunk_size)
            time.sleep(0.001)

        ff.GetCoordinates(ob_mol.OBMol)
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            at = ob_mol.OBMol.GetAtom(i + 1)
            conf.SetAtomPosition(i, [at.GetX(), at.GetY(), at.GetZ()])
        return True
    except WorkerHaltError:
        raise
    except (RuntimeError, ValueError, TypeError, ImportError) as e:
        safe_status_cb(f"Iterative optimization ({method}) error (OpenBabel): {e}")
        return False


def _parse_explicit_stereo(mol_block: str) -> Dict[int, Any]:
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


def _apply_explicit_stereo(mol: Any, explicit_stereo: Dict[int, Any]) -> None:
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
                    begin_heavy = [n for n in begin_neighbors if n.GetAtomicNum() > 1]
                    end_heavy = [n for n in end_neighbors if n.GetAtomicNum() > 1]

                    stereo_atom1 = (
                        begin_heavy[0] if begin_heavy else begin_neighbors[0]
                    ).GetIdx()
                    stereo_atom2 = (
                        end_heavy[0] if end_heavy else end_neighbors[0]
                    ).GetIdx()

                    bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                    bond.SetStereo(stereo_type)


def _perform_direct_conversion(
    mol_block: str,
    mol: Any,
    options: Optional[Dict[str, Any]],
    _check_halted: Callable[[], bool],
    _safe_status: Callable[[str], None],
) -> Optional[Any]:
    """Direct 3D conversion using 2D coordinates and adding missing H without embedding."""
    parsed_coords, stereo_dirs = [], []
    # Best-effort attempt to parse 2D coordinates from MOL block for direct conversion.
    # We suppress AttributeError/RuntimeError if the block is malformed or RDKit fails internally.
    with contextlib.suppress(AttributeError, RuntimeError, ValueError, TypeError):
        base2d_all = Chem.MolFromMolBlock(
            mol_block, removeHs=False, sanitize=True
        ) or Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
        if base2d_all and base2d_all.GetNumConformers() > 0:
            oconf = base2d_all.GetConformer()
            for i in range(base2d_all.GetNumAtoms()):
                p = oconf.GetAtomPosition(i)
                parsed_coords.append((float(p.x), float(p.y), 0.0))

    # Manual parse of bond stereo flags from MOL block text.
    # Suppress IndexError/ValueError if the block format is unexpected or lines are truncated.
    with contextlib.suppress(IndexError, ValueError, AttributeError):
        lines = mol_block.splitlines()
        counts_idx = next(
            (i for i, ln in enumerate(lines[:40]) if re.match(r"^\s*\d+\s+\d+", ln)),
            None,
        )
        if counts_idx is not None:
            parts = lines[counts_idx].split()
            natoms, nbonds = int(parts[0]), int(parts[1])
            atom_map = {i + 1: i for i in range(natoms)}
            bond_start = counts_idx + 1 + natoms
            for j in range(min(nbonds, max(0, len(lines) - bond_start))):
                bond_line = lines[bond_start + j]
                m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)(?:\s+(-?\d+))?", bond_line)
                if m:
                    a1, a2 = int(m.group(1)), int(m.group(2))
                    stereo_raw = int(m.group(4)) if m.group(4) else 0
                else:
                    f = bond_line.split()
                    if len(f) < 4:
                        continue
                    a1, a2, stereo_raw = int(f[0]), int(f[1]), int(f[3])

                stereo_flag = (
                    1 if stereo_raw == 1 else (6 if stereo_raw == 2 else stereo_raw)
                )
                if a1 in atom_map and a2 in atom_map and stereo_flag in (1, 6):
                    stereo_dirs.append((atom_map[a1], atom_map[a2], stereo_flag))

    if not parsed_coords:
        # Fallback manual parse of atom coordinates from MOL block text.
        # Suppress IndexError/ValueError if the block format is unexpected.
        with contextlib.suppress(IndexError, ValueError, AttributeError):
            lines = mol_block.splitlines()
            idx = next(
                (
                    i
                    for i, ln in enumerate(lines[:40])
                    if re.match(r"^\s*\d+\s+\d+", ln)
                ),
                None,
            )
            if idx is not None:
                n = int(lines[idx].split()[0])
                for j in range(min(n, len(lines) - (idx + 1))):
                    ln = lines[idx + 1 + j]
                    try:
                        x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
                    except ValueError:
                        f = ln.split()
                        x, y, z = float(f[0]), float(f[1]), float(f[2])
                    parsed_coords.append((x, y, z))

    if not parsed_coords:
        raise ValueError("Failed to parse coordinates for direct conversion.")

    num_existing = len(parsed_coords)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        if i < num_existing:
            x, y, _ = parsed_coords[i]
            conf.SetAtomPosition(i, rdGeometry.Point3D(float(x), float(y), 0.0))
        else:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() == 1:
                neighs = [n for n in atom.GetNeighbors() if n.GetIdx() < num_existing]
                if neighs:
                    parent_idx = neighs[0].GetIdx()
                    parent_pos = conf.GetAtomPosition(parent_idx)
                    vecs = []
                    for nbr in mol.GetAtomWithIdx(parent_idx).GetNeighbors():
                        if nbr.GetIdx() < i:
                            p = conf.GetAtomPosition(nbr.GetIdx())
                            vx, vy = (
                                float(p.x) - float(parent_pos.x),
                                float(p.y) - float(parent_pos.y),
                            )
                            nrm = math.hypot(vx, vy)
                            if nrm > 1e-6:
                                vecs.append((vx / nrm, vy / nrm))

                    if vecs:
                        fx, fy = -sum(v[0] for v in vecs), -sum(v[1] for v in vecs)
                        fn = math.hypot(fx, fy)
                        if fn < 1e-6:
                            fx, fy, fn = (
                                -vecs[0][1],
                                vecs[0][0],
                                math.hypot(-vecs[0][1], vecs[0][0]),
                            )
                        angle = sum(
                            1
                            for n in mol.GetAtomWithIdx(parent_idx).GetNeighbors()
                            if n.GetIdx() < i and n.GetAtomicNum() == 1
                        ) * (math.pi / 6.0)
                        rx, ry = (
                            (fx / fn) * math.cos(angle) - (fy / fn) * math.sin(angle),
                            (fx / fn) * math.sin(angle) + (fy / fn) * math.cos(angle),
                        )
                        conf.SetAtomPosition(
                            i,
                            rdGeometry.Point3D(
                                float(parent_pos.x) + rx, float(parent_pos.y) + ry, 0.3
                            ),
                        )
                    else:
                        conf.SetAtomPosition(
                            i,
                            rdGeometry.Point3D(
                                float(parent_pos.x) + 0.5,
                                float(parent_pos.y) + 0.5,
                                0.3,
                            ),
                        )
                else:
                    conf.SetAtomPosition(i, rdGeometry.Point3D(0.0, 0.0, 0.1))

    # Best-effort application of stereo directions to Z-coordinates for direct conversion.
    with contextlib.suppress(AttributeError, RuntimeError, TypeError, IndexError):
        for b, e, flag in stereo_dirs:
            if b < num_existing and e < num_existing:
                pos = conf.GetAtomPosition(e)
                conf.SetAtomPosition(
                    e,
                    rdGeometry.Point3D(
                        float(pos.x),
                        float(pos.y),
                        float(pos.z) + (1.5 if flag == 1 else -1.5),
                    ),
                )

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    if (options or {}).get("do_optimize", True):
        _adjust_collision_avoidance(mol, _check_halted, _safe_status)
        if _check_halted():
            raise WorkerHaltError("Halted")

        opt_method = (options or {}).get("optimization_method", "MMFF94s_RDKIT")
        backend = "OBABEL" if "OBABEL" in opt_method.upper() else "RDKIT"
        method_key = (
            "UFF"
            if "UFF" in opt_method.upper()
            else (
                "GAFF"
                if "GAFF" in opt_method.upper()
                else ("GHEMICAL" if "GHEMICAL" in opt_method.upper() else "MMFF94s")
            )
        )
        if "MMFF94" in opt_method.upper() and "MMFF94S" not in opt_method.upper():
            method_key = "MMFF94"

        # Best-effort property assignment for UI feedback
        with contextlib.suppress(AttributeError, RuntimeError, TypeError):
            mol.SetProp("_pme_optimization_method", opt_method)
            mol.SetProp("_pme_conversion_backend", "Direct")
        _safe_status(f"Optimizing ({method_key} / {backend})...")

        opt_func = (
            _iterative_optimize_obabel if backend == "OBABEL" else _iterative_optimize
        )
        if not opt_func(
            mol,
            method_key,
            _check_halted,
            _safe_status,
            options=options if backend == "RDKIT" else None,
        ):
            raise RuntimeError(f"Optimization with {opt_method} failed.")

    if _check_halted():
        raise WorkerHaltError("Halted")
    return mol


def _perform_optimize_only(
    mol: Any,
    options: Optional[Dict[str, Any]],
    worker_id: Any,
    _check_halted: Callable[[], bool],
    _safe_status: Callable[[str], None],
    _safe_finished: Callable[[Any], None],
) -> None:
    """Perform optimization on an existing 3D structure."""
    _safe_status("Optimizing existing 3D structure...")
    opt_method = str((options or {}).get("optimization_method", "MMFF_RDKIT")).upper()
    backend = "OBABEL" if "OBABEL" in opt_method else "RDKIT"
    method_key = (
        "UFF"
        if "UFF" in opt_method
        else (
            "GAFF"
            if "GAFF" in opt_method
            else ("GHEMICAL" if "GHEMICAL" in opt_method else "MMFF94s")
        )
    )
    if "MMFF94" in opt_method and "MMFF94S" not in opt_method:
        method_key = "MMFF94"

    # Best-effort property assignment for UI feedback
    with contextlib.suppress(AttributeError, RuntimeError, TypeError):
        mol.SetProp("_pme_optimization_method", opt_method)

    opt_func = (
        _iterative_optimize_obabel if backend == "OBABEL" else _iterative_optimize
    )
    if not opt_func(
        mol,
        method_key,
        _check_halted,
        _safe_status,
        options=options if backend == "RDKIT" else None,
    ):
        _safe_status(
            f"Warning: Optimization with {opt_method} failed. Structure preserved."
        )
        with contextlib.suppress(Exception):
            mol.ClearProp("_pme_optimization_method")

    # Final status message before finishing (to ensure it doesn't overwrite error/halt messages)
    opt_label = _OPT_METHOD_LABELS.get(opt_method, opt_method)
    _safe_status(f"Process completed (Existing 3D Structure / {opt_label}).")
    _safe_finished((worker_id, mol))


def _perform_obabel_conversion(
    mol_block: str,
    mode: str,
    opt_method: Optional[str],
    worker_id: Any,
    options: Optional[Dict[str, Any]],
    _check_halted: Callable[[], bool],
    _safe_status: Callable[[str], None],
    _safe_finished: Callable[[Any], None],
) -> bool:
    """Perform Open Babel 3D conversion and optimization."""
    _safe_status("Attempting Open Babel conversion...")
    try:
        if not OBABEL_AVAILABLE or not pybel:
            raise RuntimeError("Open Babel not available.")

        # Isolate make3D in a subprocess to prevent hard C++ aborts from crashing the main GUI
        script = """
import sys, contextlib
from openbabel import pybel
import os
# Inherit babel variables if any
if "BABEL_DATADIR" in os.environ: os.environ["BABEL_DATADIR"] = os.environ["BABEL_DATADIR"]
mol_block = sys.stdin.read()
ob_mol = pybel.readstring("mol", mol_block)
with contextlib.suppress(Exception): ob_mol.addh()
ob_mol.make3D()
print(ob_mol.write("mol"))
        """

        try:
            result = subprocess.run(
                [sys.executable, "-c", script],
                input=mol_block,
                text=True,
                capture_output=True,
                timeout=20,  # Prevent infinite hangs in OBabel
            )
            if result.returncode != 0 or not result.stdout.strip():
                raise RuntimeError(
                    f"Subprocess crashed or failed. Error: {result.stderr.strip()}"
                )
            out_mol_block = result.stdout
        except subprocess.TimeoutExpired:
            raise RuntimeError("Open Babel make3D() timed out.")
        except Exception as e:
            raise RuntimeError(f"Open Babel isolated execution failed: {e}")

        rd_mol = Chem.AddHs(Chem.MolFromMolBlock(out_mol_block, removeHs=False))
        if not rd_mol:
            raise ValueError("Open Babel produced invalid MOL.")

        _adjust_collision_avoidance(rd_mol, _check_halted, _safe_status)
        opt_method = opt_method or "MMFF94s_RDKIT"
        backend = "OBABEL" if "OBABEL" in opt_method.upper() else "RDKIT"
        method_key = (
            "UFF"
            if "UFF" in opt_method.upper()
            else (
                "GAFF"
                if "GAFF" in opt_method.upper()
                else ("GHEMICAL" if "GHEMICAL" in opt_method.upper() else "MMFF94s")
            )
        )
        if "MMFF94" in opt_method.upper() and "MMFF94S" not in opt_method.upper():
            method_key = "MMFF94"

        # Best-effort property assignment for UI feedback
        with contextlib.suppress(AttributeError, RuntimeError, TypeError):
            rd_mol.SetProp("_pme_optimization_method", opt_method)
            rd_mol.SetProp("_pme_conversion_backend", "Open Babel")
        _safe_status(f"Optimizing ({method_key} / {backend})...")

        opt_func = (
            _iterative_optimize_obabel if backend == "OBABEL" else _iterative_optimize
        )
        if not opt_func(
            rd_mol,
            method_key,
            _check_halted,
            _safe_status,
            options=options if backend == "RDKIT" else None,
        ):
            _safe_status("Warning: Optimization failed. Using unoptimized structure.")
            # Best-effort property cleanup if optimization was skipped
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                rd_mol.ClearProp("_pme_optimization_method")

        if _check_halted():
            raise WorkerHaltError("Halted")
        # Final status message before finishing (to ensure it doesn't overwrite error/halt messages)
        opt_label = _OPT_METHOD_LABELS.get(opt_method, opt_method)
        _safe_status(f"Process completed (Open Babel Conversion / {opt_label}).")
        _safe_finished((worker_id, rd_mol))
        return True
    except WorkerHaltError:
        raise
    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
        if mode == "obabel":
            raise RuntimeError(f"Open Babel conversion failed: {e}")
        _safe_status(f"Open Babel failed: {e}. Falling back...")
    return False


class CalculationWorker(QObject):
    status_update = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(object)
    start_work = pyqtSignal(str, object)

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.halt_ids: Optional[Set[int]] = None
        self.start_work.connect(self.run_calculation)

    @pyqtSlot(str, object)
    def run_calculation(
        self, mol_block: str, options: Optional[Dict[str, Any]] = None
    ) -> None:
        """Main entry point for 3D coordinate generation and optimization."""
        options = options or {}
        w_id = options.get("worker_id")

        def _check_halted() -> bool:
            h_ids = getattr(self, "halt_ids", None)
            if getattr(self, "halt_all", False):
                return True  # type: ignore[return-value]
            if h_ids is None:
                return False
            return bool(
                ("ALL" in h_ids)
                or (None in h_ids)
                or (w_id is not None and w_id in h_ids)
            )

        def _safe_status(msg: str) -> None:
            if _check_halted():
                raise WorkerHaltError("Halted")
            with contextlib.suppress(AttributeError, RuntimeError):
                self.status_update.emit(msg)

        def _safe_finished(payload: Any) -> None:
            if _check_halted():
                raise WorkerHaltError("Halted")
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.finished.emit(payload)

        def _safe_error(msg: str) -> None:
            # If we're already halting, don't raise another error
            if msg == "Halted":
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    self.error.emit((w_id, "Halted"))
                return

            # If a new halt was requested during an error, raise it to the local handler
            if _check_halted():
                raise WorkerHaltError("Halted")

            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.error.emit((w_id, msg))

        helpers = {
            "check_halted": _check_halted,
            "status": _safe_status,
            "finished": _safe_finished,
        }

        try:
            if not mol_block:
                raise ValueError("No atoms to convert.")
            _safe_status("Creating 3D structure...")

            # 1. Prepare Molecule (Parsing & Stereo)
            mol, ex_stereo = self._prepare_molecule_for_calc(mol_block, helpers)
            mode = options.get("conversion_mode", "fallback")

            # 2. Optimization Only Mode
            if mode == "optimize_only":
                _perform_optimize_only(
                    mol, options, w_id, _check_halted, _safe_status, _safe_finished
                )
                return

            # 3. Add Hydrogens for full conversion modes
            mol = Chem.AddHs(mol)
            if _check_halted():
                raise WorkerHaltError("Halted")
            _apply_explicit_stereo(mol, ex_stereo)

            # 4. Mode-specific Workflows
            if mode == "direct":
                self._run_direct_workflow(mol_block, mol, options, helpers)
                return

            # 5. RDKit Workflow (Primary)
            if mode in ("fallback", "rdkit"):
                success = self._run_rdkit_workflow(mol, ex_stereo, options, helpers)
                if success:
                    return
                if mode == "rdkit":
                    raise RuntimeError(
                        "RDKit 3D conversion failed: Could not generate (embed) 3D coordinates."
                    )

            # 6. Open Babel Workflow (Fallback)
            if mode in ("fallback", "obabel"):
                success = self._run_obabel_workflow(mol_block, options, helpers)
                if success:
                    return

            # 7. Final Fallback to Direct
            if mode == "fallback":
                _safe_status("Falling back to direct conversion...")
                self._run_direct_workflow(mol_block, mol, options, helpers)

        except WorkerHaltError:
            # Swallow here; the loop has already been notified
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.error.emit((w_id, "Halted"))
        except (
            Exception,
            RuntimeError,
            ValueError,
            TypeError,
            AttributeError,
            ImportError,
            OSError,
            UnicodeDecodeError,
        ) as e:
            try:
                _safe_error(str(e))
            except WorkerHaltError:
                # Swallowed if safe_error itself detected a halt
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    self.error.emit((w_id, "Halted"))

    def _prepare_molecule_for_calc(
        self, mol_block: str, helpers: Dict[str, Any]
    ) -> Tuple[Any, Dict[int, Any]]:
        """Parse MOL block and extract explicit stereochemistry info."""
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if not mol:
            raise ValueError("Failed to parse MOL block.")
        if helpers["check_halted"]():
            raise WorkerHaltError("Halted")

        ex_stereo = _parse_explicit_stereo(mol_block)
        _apply_explicit_stereo(mol, ex_stereo)
        if helpers["check_halted"]():
            raise WorkerHaltError("Halted")
        return mol, ex_stereo

    def _run_rdkit_workflow(
        self,
        mol: Any,
        ex_stereo: Dict[int, Any],
        options: Dict[str, Any],
        helpers: Dict[str, Any],
    ) -> bool:
        """Execute RDKit ETKDG embedding and force-field optimization."""
        _check_halted = helpers["check_halted"]
        _safe_status = helpers["status"]
        _safe_finished = helpers["finished"]
        w_id = options.get("worker_id")

        params = AllChem.ETKDGv2()
        params.randomSeed = 42
        params.useExpTorsionAnglePrefs = params.useBasicKnowledge = (
            params.enforceChirality
        ) = True

        # Track existing/explicit stereo to restore after embedding
        orig_stereo = []
        for b_idx, s_type in ex_stereo.items():
            if b_idx < mol.GetNumBonds():
                b = mol.GetBondWithIdx(b_idx)
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    orig_stereo.append((b.GetIdx(), s_type, b.GetStereoAtoms()))
        for b in mol.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and b.GetStereo() != Chem.BondStereo.STEREONONE
                and b.GetIdx() not in ex_stereo
            ):
                orig_stereo.append((b.GetIdx(), b.GetStereo(), b.GetStereoAtoms()))

        _safe_status("RDKit: Embedding 3D coordinates...")
        if _check_halted():
            raise WorkerHaltError("Halted")

        conf_id = -1
        with contextlib.suppress(RuntimeError, ValueError):
            conf_id = AllChem.EmbedMolecule(mol, params)

        # Fallback embedding strategies
        if conf_id == -1:
            with contextlib.suppress(
                AttributeError, RuntimeError, ValueError, TypeError
            ):
                bm = AllChem.GetMoleculeBoundsMatrix(mol)
                for b_idx, s, satoms in orig_stereo:
                    if len(satoms) == 2:
                        t = 3.0 if s == Chem.BondStereo.STEREOZ else 5.0
                        bm[satoms[0]][satoms[1]] = bm[satoms[1]][satoms[0]] = t
                DoTriangleSmoothing(bm)
                conf_id = AllChem.EmbedMolecule(mol, bm, params)

        if conf_id == -1:
            with contextlib.suppress(
                AttributeError, RuntimeError, ValueError, TypeError
            ):
                conf_id = AllChem.EmbedMolecule(mol, AllChem.ETKDGv2(randomSeed=42))

        if conf_id == -1:
            return False

        # Restoration and Optimization
        for b_idx, s, satoms in orig_stereo:
            b = mol.GetBondWithIdx(b_idx)
            if len(satoms) == 2:
                b.SetStereoAtoms(satoms[0], satoms[1])
            b.SetStereo(s)

        _adjust_collision_avoidance(mol, _check_halted, _safe_status)
        opt_method = options.get("optimization_method", "MMFF94s_RDKIT")
        backend = "OBABEL" if "OBABEL" in opt_method.upper() else "RDKIT"
        method_key = (
            "UFF"
            if "UFF" in opt_method.upper()
            else (
                "GAFF"
                if "GAFF" in opt_method.upper()
                else ("GHEMICAL" if "GHEMICAL" in opt_method.upper() else "MMFF94s")
            )
        )
        if "MMFF94" in opt_method.upper() and "MMFF94S" not in opt_method.upper():
            method_key = "MMFF94"

        with contextlib.suppress(AttributeError, RuntimeError, TypeError):
            mol.SetProp("_pme_optimization_method", opt_method)
            mol.SetProp("_pme_conversion_backend", "RDKit")
        _safe_status(f"Optimizing ({method_key} / {backend})...")

        opt_func = (
            _iterative_optimize_obabel if backend == "OBABEL" else _iterative_optimize
        )
        if not opt_func(
            mol,
            method_key,
            _check_halted,
            _safe_status,
            options=options if backend == "RDKIT" else None,
        ):
            _safe_status("Warning: Optimization failed. Using unoptimized structure.")
            with contextlib.suppress(Exception):
                mol.ClearProp("_pme_optimization_method")

        # Final stereo restoration check
        for b_idx, s, satoms in orig_stereo:
            b = mol.GetBondWithIdx(b_idx)
            if len(satoms) == 2:
                b.SetStereoAtoms(satoms[0], satoms[1])
            b.SetStereo(s)

        if _check_halted():
            raise WorkerHaltError("Halted")
        # Final status message before finishing (to ensure it doesn't overwrite error/halt messages)
        opt_label = _OPT_METHOD_LABELS.get(opt_method, opt_method)
        _safe_status(f"Process completed (RDKit Conversion / {opt_label}).")
        _safe_finished((w_id, mol))
        return True

    def _run_obabel_workflow(
        self, mol_block: str, options: Dict[str, Any], helpers: Dict[str, Any]
    ) -> bool:
        """Execute Open Babel 3D conversion."""
        return _perform_obabel_conversion(
            mol_block,
            options.get("conversion_mode", "fallback"),
            options.get("optimization_method"),
            options.get("worker_id"),
            options,
            helpers["check_halted"],
            helpers["status"],
            helpers["finished"],
        )

    def _run_direct_workflow(
        self, mol_block: str, mol: Any, options: Dict[str, Any], helpers: Dict[str, Any]
    ) -> None:
        """Execute direct 2D->3D conversion."""
        mol = _perform_direct_conversion(
            mol_block, mol, options, helpers["check_halted"], helpers["status"]
        )
        # Final status message before finishing (to ensure it doesn't overwrite error/halt messages)
        _safe_status = helpers["status"]
        opt_method = (options or {}).get("optimization_method", "MMFF94s_RDKIT")
        if (options or {}).get("do_optimize", True):
            opt_label = _OPT_METHOD_LABELS.get(opt_method, opt_method)
            _safe_status(f"Process completed (Direct 2D->3D Conversion / {opt_label}).")
        else:
            _safe_status("Process completed (Direct 2D->3D Conversion).")
        helpers["finished"]((options.get("worker_id"), mol))
