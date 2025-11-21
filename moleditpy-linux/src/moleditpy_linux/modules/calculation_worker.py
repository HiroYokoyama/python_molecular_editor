#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import QObject

from PyQt6.QtCore import pyqtSignal, pyqtSlot

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdGeometry
from rdkit.DistanceGeometry import DoTriangleSmoothing
import math
import re


# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
# Prefer package-relative import when running as `python -m moleditpy` and
# fall back to a top-level import when running as a script. This mirrors the
# import style used in other modules and keeps the package robust.
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None

class CalculationWorker(QObject):
    status_update = pyqtSignal(str)
    finished = pyqtSignal(object)
    error = pyqtSignal(object)  # emit (worker_id, msg) tuples for robustness
    # Per-worker start signal to avoid sharing a single MainWindow signal
    # among many worker instances (which causes race conditions and stale
    # workers being started on a single emission).
    start_work = pyqtSignal(str, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Connect the worker's own start signal to its run slot. This
        # guarantees that only this worker will respond when start_work
        # is emitted (prevents cross-talk between workers).
        try:
            self.start_work.connect(self.run_calculation)
        except Exception:
            # Be defensive: if connection fails, continue; the caller may
            # fallback to emitting directly.
            pass

    @pyqtSlot(str, object)
    def run_calculation(self, mol_block, options=None):
        try:
            # The worker may be asked to halt via a shared set `halt_ids` and
            # identifies its own run by options['worker_id'] (int).
            worker_id = None
            try:
                worker_id = options.get('worker_id') if options else None
            except Exception:
                worker_id = None

            # If a caller starts a worker without providing a worker_id, treat
            # it as a "global" worker that can still be halted via a global
            # halt flag. Emit a single status warning so callers know that
            # the worker was started without an identifier.
            _warned_no_worker_id = False
            if worker_id is None:
                try:
                    # best-effort, swallow any errors (signals may not be connected)
                    self.status_update.emit("Warning: worker started without 'worker_id'; will listen for global halt signals.")
                except Exception:
                    pass
                _warned_no_worker_id = True

            def _check_halted():
                try:
                    halt_ids = getattr(self, 'halt_ids', None)
                    # If worker_id is None, allow halting via a global mechanism:
                    #  - an explicit attribute `halt_all` set to True on the worker
                    #  - the shared `halt_ids` set containing None or the sentinel 'ALL'
                    if worker_id is None:
                        if getattr(self, 'halt_all', False):
                            return True
                        if halt_ids is None:
                            return False
                        # Support both None-in-set and string sentinel for compatibility
                        return (None in halt_ids) or ('ALL' in halt_ids)

                    if halt_ids is None:
                        return False
                    return (worker_id in halt_ids)
                except Exception:
                    return False

            # Safe-emission helpers: do nothing if this worker has been halted.
            def _safe_status(msg):
                try:
                    if _check_halted():
                        return
                    self.status_update.emit(msg)
                except Exception:
                    # Swallow any signal-emission errors to avoid crashing the worker
                    pass

            def _safe_finished(payload):
                try:
                    # Attempt to emit the payload; preserve existing fallback behavior
                    try:
                        self.finished.emit(payload)
                    except TypeError:
                        # Some slots/old code may expect a single-molecule arg; try that too
                        try:
                            # If payload was a tuple like (worker_id, mol), try sending the second element
                            if isinstance(payload, (list, tuple)) and len(payload) >= 2:
                                self.finished.emit(payload[1])
                            else:
                                self.finished.emit(payload)
                        except Exception:
                            pass
                except Exception:
                    pass

            def _safe_error(msg):
                try:
                    # Emit a tuple containing the worker_id (may be None) and the message
                    try:
                        self.error.emit((worker_id, msg))
                    except Exception:
                        # Fallback to emitting the raw message if tuple emission fails for any reason
                        try:
                            self.error.emit(msg)
                        except Exception:
                            pass
                except Exception:
                    pass

            # options: dict-like with keys: 'conversion_mode' -> 'fallback'|'rdkit'|'obabel'|'direct'
            if options is None:
                options = {}
            conversion_mode = options.get('conversion_mode', 'fallback')
            # Ensure params exists in all code paths (some RDKit calls below
            # reference `params` and earlier editing introduced a path where
            # it might not be defined). Initialize to None here and assign
            # a proper ETKDG params object later where needed.
            params = None
            if not mol_block:
                raise ValueError("No atoms to convert.")
            
            _safe_status("Creating 3D structure...")

            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None:
                raise ValueError("Failed to create molecule from MOL block.")

            # Check early whether this run has been requested to halt
            if _check_halted():
                raise RuntimeError("Halted")

            # CRITICAL FIX: Extract and restore explicit E/Z labels from MOL block
            # Parse M CFG lines to get explicit stereo labels
            explicit_stereo = {}
            mol_lines = mol_block.split('\n')
            for line in mol_lines:
                if line.startswith('M  CFG'):
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

            # Force explicit stereo labels regardless of coordinates
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Find suitable stereo atoms
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        
                        # Pick heavy atom neighbors preferentially
                        begin_neighbors = [nbr for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_atom.GetIdx()]
                        end_neighbors = [nbr for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_atom.GetIdx()]
                        
                        if begin_neighbors and end_neighbors:
                            # Prefer heavy atoms
                            begin_heavy = [n for n in begin_neighbors if n.GetAtomicNum() > 1]
                            end_heavy = [n for n in end_neighbors if n.GetAtomicNum() > 1]
                            
                            stereo_atom1 = (begin_heavy[0] if begin_heavy else begin_neighbors[0]).GetIdx()
                            stereo_atom2 = (end_heavy[0] if end_heavy else end_neighbors[0]).GetIdx()
                            
                            bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                            bond.SetStereo(stereo_type)

            # Do NOT call AssignStereochemistry here as it overrides our explicit labels

            mol = Chem.AddHs(mol)

            # Check after adding Hs (may be a long operation)
            if _check_halted():
                raise RuntimeError("Halted")

            # CRITICAL: Re-apply explicit stereo after AddHs which may renumber atoms
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Re-find suitable stereo atoms after hydrogen addition
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        
                        # Pick heavy atom neighbors preferentially
                        begin_neighbors = [nbr for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_atom.GetIdx()]
                        end_neighbors = [nbr for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_atom.GetIdx()]
                        
                        if begin_neighbors and end_neighbors:
                            # Prefer heavy atoms
                            begin_heavy = [n for n in begin_neighbors if n.GetAtomicNum() > 1]
                            end_heavy = [n for n in end_neighbors if n.GetAtomicNum() > 1]
                            
                            stereo_atom1 = (begin_heavy[0] if begin_heavy else begin_neighbors[0]).GetIdx()
                            stereo_atom2 = (end_heavy[0] if end_heavy else end_neighbors[0]).GetIdx()
                            
                            bond.SetStereoAtoms(stereo_atom1, stereo_atom2)
                            bond.SetStereo(stereo_type)

            # Direct mode: construct a 3D conformer without embedding by using
            # the 2D coordinates from the MOL block (z=0) and placing added
            # hydrogens close to their parent heavy atoms with a small z offset.
            # This avoids 3D embedding entirely and is useful for quick viewing
            # when stereochemistry/geometry refinement is not desired.
            if conversion_mode == 'direct':
                _safe_status("Direct conversion: using 2D coordinates + adding missing H (no embedding).")
                try:
                    # 1) Parse MOL block *with* existing hydrogens (removeHs=False)
                    #    to get coordinates for *all existing* atoms.
                    parsed_coords = []  # all-atom coordinates (x, y, z)
                    stereo_dirs = []    # list of (begin_idx, end_idx, stereo_flag)
                    
                    base2d_all = None
                    try:
                        # H原子を含めてパース
                        base2d_all = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
                    except Exception:
                        try:
                            base2d_all = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
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
                            
                            # 全原子マップ (MOL 1-based index -> 0-based index)
                            atom_map = {i + 1: i for i in range(natoms)}
                            
                            bond_start = counts_idx + 1 + natoms
                            for j in range(min(nbonds, max(0, len(lines) - bond_start))):
                                bond_line = lines[bond_start + j]
                                try:
                                    m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)(?:\s+(-?\d+))?", bond_line)
                                    if m:
                                        try:
                                            atom1_mol = int(m.group(1))  # 1-based MOL index
                                            atom2_mol = int(m.group(2))  # 1-based MOL index
                                        except Exception:
                                            continue
                                        try:
                                            stereo_raw = int(m.group(4)) if m.group(4) is not None else 0
                                        except Exception:
                                            stereo_raw = 0
                                    else:
                                        fields = bond_line.split()
                                        if len(fields) >= 4:
                                            try:
                                                atom1_mol = int(fields[0])  # 1-based MOL index
                                                atom2_mol = int(fields[1])  # 1-based MOL index
                                            except Exception:
                                                continue
                                            try:
                                                stereo_raw = int(fields[3]) if len(fields) > 3 else 0
                                            except Exception:
                                                stereo_raw = 0
                                        else:
                                            continue

                                    # V2000の立体表記を正規化
                                    if stereo_raw == 1:
                                        stereo_flag = 1 # Wedge
                                    elif stereo_raw == 2:
                                        stereo_flag = 6 # Dash (V2000では 6 がDash)
                                    else:
                                        stereo_flag = stereo_raw

                                    # 全原子マップでチェック
                                    if atom1_mol in atom_map and atom2_mol in atom_map:
                                        idx1 = atom_map[atom1_mol]
                                        idx2 = atom_map[atom2_mol]
                                        if stereo_flag in (1, 6): # Wedge (1) or Dash (6)
                                            stereo_dirs.append((idx1, idx2, stereo_flag))
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
                                for j in range(min(natoms, max(0, len(lines) - atom_start))):
                                    atom_line = lines[atom_start + j]
                                    try:
                                        x = float(atom_line[0:10].strip()); y = float(atom_line[10:20].strip()); z = float(atom_line[20:30].strip())
                                    except Exception:
                                        fields = atom_line.split()
                                        if len(fields) >= 4:
                                            try:
                                                x = float(fields[0]); y = float(fields[1]); z = float(fields[2])
                                            except Exception:
                                                continue
                                        else:
                                            continue
                                    # H原子もスキップしない
                                    parsed_coords.append((x, y, z))
                        except Exception:
                            parsed_coords = []
                    
                    if not parsed_coords:
                        raise ValueError("Failed to parse coordinates from MOL block for direct conversion.")

                    # 3) `mol` は既に AddHs された状態
                    #    元の原子数 (H含む) を parsed_coords の長さから取得
                    num_existing_atoms = len(parsed_coords)

                    # 4) コンフォーマを作成
                    conf = Chem.Conformer(mol.GetNumAtoms())

                    for i in range(mol.GetNumAtoms()):
                        if i < num_existing_atoms:
                            # 既存原子 (H含む): 2D座標 (z=0) を設定
                            x, y, z_ignored = parsed_coords[i]
                            try:
                                conf.SetAtomPosition(i, rdGeometry.Point3D(float(x), float(y), 0.0))
                            except Exception:
                                pass
                        else:
                            # 新規追加されたH原子: 親原子の近くに配置
                            atom = mol.GetAtomWithIdx(i)
                            if atom.GetAtomicNum() == 1:
                                neighs = [n for n in atom.GetNeighbors() if n.GetIdx() < num_existing_atoms]
                                heavy_pos_found = False
                                for nb in neighs: # 親原子 (重原子または既存H)
                                    try:
                                        nb_idx = nb.GetIdx()
                                        # if nb_idx < num_existing_atoms: # チェックは不要 (neighs で既にフィルタ済み)
                                        nbpos = conf.GetAtomPosition(nb_idx)
                                        # Geometry-based placement:
                                        # Compute an "empty" direction around the parent atom by
                                        # summing existing bond unit vectors and taking the
                                        # opposite. If degenerate, pick a perpendicular or
                                        # fallback vector. Rotate slightly if multiple Hs already
                                        # attached to avoid overlap.
                                        parent_idx = nb_idx
                                        try:
                                            parent_pos = conf.GetAtomPosition(parent_idx)
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
                                                        vx = float(p.x) - float(parent_pos.x)
                                                        vy = float(p.y) - float(parent_pos.y)
                                                        nrm = math.hypot(vx, vy)
                                                        if nrm > 1e-6:
                                                            vecs.append((vx / nrm, vy / nrm))
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
                                                fx /= fn; fy /= fn

                                                # Avoid placing multiple Hs at identical directions
                                                existing_h_count = sum(1 for nbr in parent_atom.GetNeighbors()
                                                                       if nbr.GetIdx() < i and nbr.GetAtomicNum() == 1)
                                                angle = existing_h_count * (math.pi / 6.0)  # 30deg steps
                                                cos_a = math.cos(angle); sin_a = math.sin(angle)
                                                rx = fx * cos_a - fy * sin_a
                                                ry = fx * sin_a + fy * cos_a

                                                bond_length = 1.0
                                                conf.SetAtomPosition(i, rdGeometry.Point3D(
                                                    float(parent_pos.x) + rx * bond_length,
                                                    float(parent_pos.y) + ry * bond_length,
                                                    0.3
                                                ))
                                            else:
                                                # No existing placed neighbors: fallback to small offset
                                                conf.SetAtomPosition(i, rdGeometry.Point3D(
                                                    float(parent_pos.x) + 0.5,
                                                    float(parent_pos.y) + 0.5,
                                                    0.3
                                                ))

                                            heavy_pos_found = True
                                            break
                                        except Exception:
                                            # fall back to trying the next neighbor if any
                                            continue
                                    except Exception:
                                        continue
                                if not heavy_pos_found:
                                    # フォールバック (原点近く)
                                    try:
                                        conf.SetAtomPosition(i, rdGeometry.Point3D(0.0, 0.0, 0.10))
                                    except Exception:
                                        pass
                    
                    # 5) Wedge/Dash の Zオフセットを適用
                    try:
                        stereo_z_offset = 1.5  # wedge -> +1.5, dash -> -1.5
                        for begin_idx, end_idx, stereo_flag in stereo_dirs:
                            try:
                                # インデックスは既存原子内のはず
                                if begin_idx >= num_existing_atoms or end_idx >= num_existing_atoms:
                                    continue
                                    
                                if stereo_flag not in (1, 6):
                                    continue
                                
                                sign = 1.0 if stereo_flag == 1 else -1.0
                                
                                # end_idx (立体表記の終点側原子) にZオフセットを適用
                                pos = conf.GetAtomPosition(end_idx)
                                newz = float(pos.z) + (stereo_z_offset * sign) # 既存のZ=0にオフセットを加算
                                conf.SetAtomPosition(end_idx, rdGeometry.Point3D(float(pos.x), float(pos.y), float(newz)))
                            except Exception:
                                continue
                    except Exception:
                        pass
                    
                    # コンフォーマを入れ替えて終了
                    try:
                        mol.RemoveAllConformers()
                    except Exception:
                        pass
                    mol.AddConformer(conf, assignId=True)
                    
                    if _check_halted():
                        raise RuntimeError("Halted (after optimization)")
                    try:
                        _safe_finished((worker_id, mol))
                    except Exception:
                        _safe_finished(mol)
                    _safe_status("Direct conversion completed.")
                    return
                except Exception as e:
                    _safe_status(f"Direct conversion failed: {e}")

            params = AllChem.ETKDGv2()
            params.randomSeed = 42
            # CRITICAL: Force ETKDG to respect the existing stereochemistry
            params.useExpTorsionAnglePrefs = True
            params.useBasicKnowledge = True
            params.enforceChirality = True  # This is critical for stereo preservation
            
            # Store original stereochemistry before embedding (prioritizing explicit labels)
            original_stereo_info = []
            for bond_idx, stereo_type in explicit_stereo.items():
                if bond_idx < mol.GetNumBonds():
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        stereo_atoms = bond.GetStereoAtoms()
                        original_stereo_info.append((bond.GetIdx(), stereo_type, stereo_atoms))
            
            # Also store any other stereo bonds not in explicit_stereo
            for bond in mol.GetBonds():
                if (bond.GetBondType() == Chem.BondType.DOUBLE and 
                    bond.GetStereo() != Chem.BondStereo.STEREONONE and
                    bond.GetIdx() not in explicit_stereo):
                    stereo_atoms = bond.GetStereoAtoms()
                    original_stereo_info.append((bond.GetIdx(), bond.GetStereo(), stereo_atoms))
            
            # Only report RDKit-specific messages when RDKit embedding will be
            # attempted. For other conversion modes, emit clearer, non-misleading
            # status messages so the UI doesn't show "RDKit" when e.g. direct
            # coordinates or Open Babel will be used.
            if conversion_mode in ('fallback', 'rdkit'):
                _safe_status("RDKit: Embedding 3D coordinates...")
            elif conversion_mode == 'obabel':
                pass
            else:
                # direct mode (or any other explicit non-RDKit mode)
                pass
            if _check_halted():
                raise RuntimeError("Halted")
            
            # Try multiple times with different approaches if needed
            conf_id = -1
            
            # First attempt: Standard ETKDG with stereo enforcement
            try:
                # Only attempt RDKit embedding if mode allows
                if conversion_mode in ('fallback', 'rdkit'):
                    conf_id = AllChem.EmbedMolecule(mol, params)
                else:
                    conf_id = -1
                # Final check before returning success
                if _check_halted():
                    raise RuntimeError("Halted")
            except Exception as e:
                # Standard embedding failed; report and continue to fallback attempts
                _safe_status(f"Standard embedding failed: {e}")

                # Second attempt: Use constraint embedding if available (only when RDKit is allowed)
                if conf_id == -1 and conversion_mode in ('fallback', 'rdkit'):
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
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = min(bounds_matrix[neighbor1_idx][neighbor2_idx], target_dist)
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = min(bounds_matrix[neighbor2_idx][neighbor1_idx], target_dist)
                                elif stereo == Chem.BondStereo.STEREOE:
                                    # E configuration: set longer distance constraint  
                                    target_dist = 5.0  # Angstroms
                                    bounds_matrix[neighbor1_idx][neighbor2_idx] = max(bounds_matrix[neighbor1_idx][neighbor2_idx], target_dist)
                                    bounds_matrix[neighbor2_idx][neighbor1_idx] = max(bounds_matrix[neighbor2_idx][neighbor1_idx], target_dist)

                        DoTriangleSmoothing(bounds_matrix)
                        conf_id = AllChem.EmbedMolecule(mol, bounds_matrix, params)
                        _safe_status("Constraint-based embedding succeeded")
                    except Exception:
                        # Constraint embedding failed: only raise error if mode is 'rdkit', otherwise allow fallback
                        _safe_status("RDKit: Constraint embedding failed")
                        if conversion_mode == 'rdkit':
                            raise RuntimeError("RDKit: Constraint embedding failed")
                        conf_id = -1
                    
            # Fallback: Try basic embedding
            if conf_id == -1:
                try:
                    if conversion_mode in ('fallback', 'rdkit'):
                        basic_params = AllChem.ETKDGv2()
                        basic_params.randomSeed = 42
                        conf_id = AllChem.EmbedMolecule(mol, basic_params)
                    else:
                        conf_id = -1
                except Exception:
                    pass
            '''
            if conf_id == -1:
                        _safe_status("Initial embedding failed, retrying with ignoreSmoothingFailures=True...")
                # Try again with ignoreSmoothingFailures instead of random-seed retries
                params.ignoreSmoothingFailures = True
                # Use a deterministic seed to avoid random-coordinate behavior here
                params.randomSeed = 0
                conf_id = AllChem.EmbedMolecule(mol, params)

            if conf_id == -1:
                self.status_update.emit("Random-seed retry failed, attempting with random coordinates...")
                try:
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True, ignoreSmoothingFailures=True)
                except TypeError:
                    # Some RDKit versions expect useRandomCoords in params
                    params.useRandomCoords = True
                    conf_id = AllChem.EmbedMolecule(mol, params)
            '''

            # Determine requested MMFF variant from options (fall back to MMFF94s)
            opt_method = None
            try:
                opt_method = options.get('optimization_method') if options else None
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
                
                try:
                    mmff_variant = "MMFF94s"
                    if opt_method and str(opt_method).upper() == 'MMFF94_RDKIT':
                        mmff_variant = "MMFF94"
                    if _check_halted():
                        raise RuntimeError("Halted")
                    AllChem.MMFFOptimizeMolecule(mol, mmffVariant=mmff_variant)
                except Exception:
                    # fallback to UFF if MMFF fails
                    try:
                        if _check_halted():
                            raise RuntimeError("Halted")
                        AllChem.UFFOptimizeMolecule(mol)
                    except Exception:
                        pass
                
                # CRITICAL: Restore stereochemistry again after optimization (explicit labels priority)
                for bond_idx, stereo, stereo_atoms in original_stereo_info:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if len(stereo_atoms) == 2:
                        bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
                    bond.SetStereo(stereo)
                
                # Do NOT call AssignStereochemistry here as it would override our explicit labels
                # Include worker_id so the main thread can ignore stale results
                # CRITICAL: Check for halt *before* emitting finished signal
                if _check_halted():
                    raise RuntimeError("Halted (after optimization)")
                try:
                    _safe_finished((worker_id, mol))
                except Exception:
                    # Fallback to legacy single-arg emit
                    _safe_finished(mol)
                _safe_status("RDKit 3D conversion succeeded.")
                return

            # If RDKit did not produce a conf and OBabel is allowed, try Open Babel
            if conf_id == -1 and conversion_mode in ('fallback', 'obabel'):
                _safe_status("RDKit embedding failed or disabled. Attempting Open Babel...")
                try:
                    if not OBABEL_AVAILABLE:
                        raise RuntimeError("Open Babel (pybel) is not available in this Python environment.")
                    ob_mol = pybel.readstring("mol", mol_block)
                    try:
                        ob_mol.addh()
                    except Exception:
                        pass
                    ob_mol.make3D()
                    try:
                        _safe_status("Optimizing with Open Babel (MMFF94)...")
                        if _check_halted():
                            raise RuntimeError("Halted")
                        ob_mol.localopt(forcefield='mmff94', steps=500)
                    except Exception:
                        try:
                            _safe_status("MMFF94 failed, falling back to UFF...")
                            if _check_halted():
                                raise RuntimeError("Halted")
                            ob_mol.localopt(forcefield='uff', steps=500)
                        except Exception:
                            _safe_status("UFF optimization also failed.")
                    molblock_ob = ob_mol.write("mol")
                    rd_mol = Chem.MolFromMolBlock(molblock_ob, removeHs=False)
                    if rd_mol is None:
                        raise ValueError("Open Babel produced invalid MOL block.")
                    rd_mol = Chem.AddHs(rd_mol)
                    try:
                        mmff_variant = "MMFF94s"
                        if opt_method and str(opt_method).upper() == 'MMFF94_RDKIT':
                            mmff_variant = "MMFF94"
                        if _check_halted():
                            raise RuntimeError("Halted")
                        AllChem.MMFFOptimizeMolecule(rd_mol, mmffVariant=mmff_variant)
                    except Exception:
                        try:
                            if _check_halted():
                                raise RuntimeError("Halted")
                            AllChem.UFFOptimizeMolecule(rd_mol)
                        except Exception:
                            pass
                    _safe_status("Open Babel embedding succeeded. Warning: Conformation accuracy may be limited.")
                    # CRITICAL: Check for halt *before* emitting finished signal
                    if _check_halted():
                        raise RuntimeError("Halted (after optimization)")
                    try:
                        _safe_finished((worker_id, rd_mol))
                    except Exception:
                        _safe_finished(rd_mol)
                    return
                except Exception as ob_err:
                    raise RuntimeError(f"Open Babel 3D conversion failed: {ob_err}")

            if conf_id == -1 and conversion_mode == 'rdkit':
                raise RuntimeError("RDKit 3D conversion failed (rdkit-only mode)")

        except Exception as e:
            _safe_error(str(e))
