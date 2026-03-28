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
import math
import logging
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from PyQt6.QtCore import Qt, QPointF, QLineF, QRectF
from PyQt6.QtGui import QCursor
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsLineItem, QApplication

try:
    from .atom_item import AtomItem
    from .bond_item import BondItem
except ImportError:
    from moleditpy_linux.ui.atom_item import AtomItem
    from moleditpy_linux.ui.bond_item import BondItem

try:
    from ..utils.sip_isdeleted_safe import sip_isdeleted_safe
except ImportError:
    from moleditpy_linux.utils.sip_isdeleted_safe import sip_isdeleted_safe

try:
    from ..utils.constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE
except ImportError:
    from moleditpy_linux.utils.constants import (
        DEFAULT_BOND_LENGTH,
        SNAP_DISTANCE,
        SUM_TOLERANCE,
    )

class TemplateMixin:
    """
    Mixin class that handles all template and fragment insertion logic for MoleculeScene.
    Because this is a Mixin, `self` refers directly to the MoleculeScene instance.
    """

    def clear_template_preview(self) -> None:
        """Remove all ghost lines for template preview."""
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(
                item, "_is_template_preview", False
            ):
                if sip_isdeleted_safe(item):
                    continue
                sc = getattr(item, "scene", lambda: None)()
                if sc is None:
                    continue
                try:
                    self.removeItem(item)
                except (RuntimeError, ValueError, TypeError) as e:
                    # Best-effort: ignore removal errors during teardown if underlying C++ object is already gone
                    logging.debug(f"Could not remove template preview item: {e}")
                    pass
        self.template_context = {}
        if hasattr(self, "template_preview"):
            self.template_preview.hide()

    def _calculate_6ring_rotation(self, num_points: int, bonds_info: List[Tuple[int, int, int]], atom_items: List[Optional[AtomItem]]) -> int:
        """
        Calculate the best rotation for a 6-ring template (like benzene)
        to match existing bond orders and ensure chemical safety.
        """
        existing_orders = {}
        for k, (i_idx, j_idx, _) in enumerate(bonds_info):
            if i_idx < len(atom_items) and j_idx < len(atom_items):
                a, b = atom_items[i_idx], atom_items[j_idx]
                if a and b:
                    eb = self.find_bond_between(a, b)
                    if eb:
                        existing_orders[k] = getattr(eb, "order", 1)

        if not existing_orders:
            return 0

        orig_orders = [o for (_, _, o) in bonds_info]
        best_rot = 0
        max_score = -999

        # Case A: Fused (>= 2 edges)
        if len(existing_orders) >= 2:
            for rot in range(num_points):
                match_double_count = 0
                match_bonus = 0
                mismatch_penalty = 0
                safe_connection_score = 0

                # Identify template-side indices corresponding to fusing k
                used_template_indices = set(
                    (k + rot) % num_points for k in existing_orders
                )

                for t_idx in used_template_indices:
                    # Neighbor indices in the template
                    adj_l = (t_idx - 1) % num_points
                    adj_r = (t_idx + 1) % num_points

                    # Connection safety: Unused template neighbors that are single bonds
                    if adj_l not in used_template_indices and orig_orders[adj_l] == 1:
                        safe_connection_score += 5000
                    if adj_r not in used_template_indices and orig_orders[adj_r] == 1:
                        safe_connection_score += 5000

                for k, exist_order in existing_orders.items():
                    template_ord = orig_orders[(k + rot) % num_points]
                    if template_ord == exist_order:
                        match_bonus += 100
                        if exist_order == 2:
                            match_double_count += 1
                    else:
                        mismatch_penalty += 50

                current_score = (
                    safe_connection_score
                    + (match_double_count * 1000)
                    + match_bonus
                    - mismatch_penalty
                )
                if current_score > max_score:
                    max_score = current_score
                    best_rot = rot

        # Case B: 1-edge fuse
        elif len(existing_orders) == 1:
            k_fuse = next(iter(existing_orders.keys()))
            exist_order = existing_orders[k_fuse]
            for rot in range(num_points):
                current_score = 0
                template_ord = orig_orders[(k_fuse + rot) % num_points]

                # Leg order matching (prefer alternating or double-double for conjugation)
                if (
                    (exist_order == 1 and template_ord == 2)
                    or (exist_order == 2 and template_ord == 1)
                    or (exist_order == 2 and template_ord == 2)
                ):
                    current_score += 100

                # Check adjacent template edges for alternating arrangement
                m_adj1 = (k_fuse - 1 + rot) % num_points
                m_adj2 = (k_fuse + 1 + rot) % num_points
                if exist_order == 1:
                    if orig_orders[m_adj1] == 2:
                        current_score += 50
                    if orig_orders[m_adj2] == 2:
                        current_score += 50
                elif exist_order == 2:
                    if orig_orders[m_adj1] == 1:
                        current_score += 50
                    if orig_orders[m_adj2] == 1:
                        current_score += 50

                if current_score > max_score:
                    max_score = current_score
                    best_rot = rot

        return best_rot

    def _should_overwrite_benzene_bond(self, exist_b):
        """
        Enforce policy for benzene template insertion.
        Overwrite existing single bonds only if they can participate in the template's aromatic system
        without violating valence (no other double bonds on the atoms).
        """
        if exist_b.order != 1:
            return False

        atom1, atom2 = exist_b.atom1, exist_b.atom2

        # Atoms at both ends must not have other double bonds
        atom1_has_other_double = any(
            b.order == 2 for b in atom1.bonds if b is not exist_b
        )
        atom2_has_other_double = any(
            b.order == 2 for b in atom2.bonds if b is not exist_b
        )

        return not atom1_has_other_double and not atom2_has_other_double

    def add_molecule_fragment(
        self,
        points: List[Union[QPointF, Tuple[float, float]]],
        bonds_info: List[Tuple[int, int, int]],
        existing_items: Optional[List[AtomItem]] = None,
        symbol: str = "C",
    ) -> List[AtomItem]:
        """Add a molecular fragment (e.g., benzene template) to the scene.

        - Enforce policy of not changing existing bond orders.
        - Benzene template rotation is decided based on fused bond orders
          so that exactly two new double bonds are created.
        """

        num_points = len(points)
        atom_items = [None] * num_points

        is_benzene_template = num_points == 6 and any(o == 2 for _, _, o in bonds_info)

        def coords(p):
            if hasattr(p, "x") and hasattr(p, "y"):
                return (p.x(), p.y())
            try:
                return (p[0], p[1])
            except (AttributeError, RuntimeError, ValueError, TypeError):
                raise ValueError("point has no x/y")

        def dist_pts(a, b):
            ax, ay = coords(a)
            bx, by = coords(b)
            return math.hypot(ax - bx, ay - by)

        # --- 1) Map already clicked existing_items to template vertices ---
        existing_items = existing_items or []
        used_indices = set()
        ref_lengths = [
            dist_pts(points[i], points[j])
            for i, j, _ in bonds_info
            if i < num_points and j < num_points
        ]
        avg_len = (sum(ref_lengths) / len(ref_lengths)) if ref_lengths else 20.0
        map_threshold = max(0.5 * avg_len, 8.0)

        for ex_item in existing_items:
            try:
                ex_pos = ex_item.pos()
                best_idx, best_d = -1, float("inf")
                for i, p in enumerate(points):
                    if i in used_indices:
                        continue
                    d = dist_pts(p, ex_pos)
                    if best_d is None or d < best_d:
                        best_d, best_idx = d, i
                if best_idx != -1 and best_d <= max(map_threshold, 1.5 * avg_len):
                    atom_items[best_idx] = ex_item
                    used_indices.add(best_idx)
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress point distance mapping errors during fragment addition

        # --- 2) Enumerate existing atoms in the scene from self.data.atoms and map them ---
        mapped_atoms = {it for it in atom_items if it is not None}
        for i, p in enumerate(points):
            if atom_items[i] is not None:
                continue

            nearby = None
            best_d = float("inf")

            for atom_data in self.data.atoms.values():
                a_item = atom_data.get("item")
                if not a_item or a_item in mapped_atoms:
                    continue
                try:
                    d = dist_pts(p, a_item.pos())
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    continue
                if d < best_d:
                    best_d, nearby = d, a_item

            if nearby and best_d <= map_threshold:
                atom_items[i] = nearby
                mapped_atoms.add(nearby)

        # --- 3) Create missing vertices ---
        for i, p in enumerate(points):
            if atom_items[i] is None:
                atom_id = self.create_atom(symbol, p)
                atom_items[i] = self.data.atoms[atom_id]["item"]

        # --- 4) Determine bond array for the template (benzene rotation alignment processing) ---
        template_bonds_to_use = list(bonds_info)
        is_6ring = num_points == 6 and len(bonds_info) == 6
        template_has_double = any(o == 2 for (_, _, o) in bonds_info)

        if is_6ring and template_has_double:
            best_rot = self._calculate_6ring_rotation(
                num_points, bonds_info, atom_items
            )

            # Reflect final rotation
            orig_orders = [o for (_, _, o) in bonds_info]
            new_tb = []
            for m in range(num_points):
                i_idx, j_idx, _ = bonds_info[m]
                new_order = orig_orders[(m + best_rot) % num_points]
                new_tb.append((i_idx, j_idx, new_order))
            template_bonds_to_use = new_tb

        # --- 5) Bond Creation/Update ---
        for id1_idx, id2_idx, order in template_bonds_to_use:
            if id1_idx < len(atom_items) and id2_idx < len(atom_items):
                a_item, b_item = atom_items[id1_idx], atom_items[id2_idx]
                if not a_item or not b_item or a_item is b_item:
                    continue

                id1, id2 = a_item.atom_id, b_item.atom_id
                if id1 > id2:
                    id1, id2 = id2, id1

                exist_b = self.find_bond_between(a_item, b_item)

                if exist_b:
                    # Keep existing bonds by default unless it's a safe benzene overwrite
                    should_overwrite = (
                        is_benzene_template
                        and self._should_overwrite_benzene_bond(exist_b)
                    )

                    if should_overwrite:
                        # Update bond order if conditions met
                        exist_b.order = order
                        exist_b.stereo = 0
                        self.data.bonds[(id1, id2)]["order"] = order
                        self.data.bonds[(id1, id2)]["stereo"] = 0
                        exist_b.update()
                    else:
                        # Keep existing bond if conditions not met
                        continue
                else:
                    # Create new bond
                    self.create_bond(a_item, b_item, bond_order=order, bond_stereo=0)

        # --- 6) Redraw ---
        for at in atom_items:
            try:
                if at:
                    at.update_style()
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress style update errors during fragment addition

        return atom_items

    def update_template_preview(self, pos: QPointF) -> None:
        mode_parts = self.mode.split("_")

        # Check if this is a user template
        if len(mode_parts) >= 3 and mode_parts[1] == "user":
            self.update_user_template_preview(pos)
            return

        is_aromatic = False
        if mode_parts[1] == "benzene":
            n = 6
            is_aromatic = True
        else:
            try:
                n = int(mode_parts[1])
            except ValueError:
                return

        items_under = self.items(pos)  # top-most first
        item = None
        for it in items_under:
            if isinstance(it, (AtomItem, BondItem)):
                item = it
                break

        points, bonds_info = [], []
        l = DEFAULT_BOND_LENGTH
        self.template_context = {}

        if isinstance(item, AtomItem):
            p0 = item.pos()
            continuous_angle = math.atan2(pos.y() - p0.y(), pos.x() - p0.x())
            snap_angle_rad = math.radians(15)
            snapped_angle = round(continuous_angle / snap_angle_rad) * snap_angle_rad
            p1 = p0 + QPointF(l * math.cos(snapped_angle), l * math.sin(snapped_angle))
            points = self._calculate_polygon_from_edge(p0, p1, n)
            self.template_context["items"] = [item]

        elif isinstance(item, BondItem):
            # Snap to bond
            p0, p1 = item.atom1.pos(), item.atom2.pos()
            points = self._calculate_polygon_from_edge(
                p0, p1, n, cursor_pos=pos, use_existing_length=True
            )
            self.template_context["items"] = [item.atom1, item.atom2]

        else:
            angle_step = 2 * math.pi / n
            start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2
            points = [
                pos
                + QPointF(
                    l * math.cos(start_angle + i * angle_step),
                    l * math.sin(start_angle + i * angle_step),
                )
                for i in range(n)
            ]

        if points:
            if is_aromatic:
                bonds_info = [
                    (i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)
                ]
            else:
                bonds_info = [(i, (i + 1) % n, 1) for i in range(n)]

            self.template_context["points"] = points
            self.template_context["bonds_info"] = bonds_info

            self.template_preview.set_geometry(points, is_aromatic)

            self.template_preview.show()
            if self.views():
                self.views()[0].viewport().update()
        else:
            self.template_preview.hide()
            if self.views():
                self.views()[0].viewport().update()

    def _calculate_polygon_from_edge(
        self,
        p0: QPointF,
        p1: QPointF,
        n: int,
        cursor_pos: Optional[QPointF] = None,
        use_existing_length: bool = False,
    ) -> List[QPointF]:
        if n < 3:
            return []
        v_edge = p1 - p0
        edge_length = math.sqrt(v_edge.x() ** 2 + v_edge.y() ** 2)
        if edge_length == 0:
            return []

        target_length = edge_length if use_existing_length else DEFAULT_BOND_LENGTH

        v_edge = (v_edge / edge_length) * target_length

        if not use_existing_length:
            p1 = p0 + v_edge

        points = [p0, p1]

        interior_angle = (n - 2) * math.pi / n
        rotation_angle = math.pi - interior_angle

        if cursor_pos:
            # Note: v_edge is normalized, direction is same
            v_cursor = cursor_pos - p0
            cross_product_z = (p1 - p0).x() * v_cursor.y() - (
                p1 - p0
            ).y() * v_cursor.x()
            if cross_product_z < 0:
                rotation_angle = -rotation_angle

        cos_a, sin_a = math.cos(rotation_angle), math.sin(rotation_angle)

        current_p, current_v = p1, v_edge
        for _ in range(n - 2):
            new_vx = current_v.x() * cos_a - current_v.y() * sin_a
            new_vy = current_v.x() * sin_a + current_v.y() * cos_a
            current_v = QPointF(new_vx, new_vy)
            current_p = current_p + current_v
            points.append(current_p)
        return points

    def add_user_template_fragment(self, context):
        """Place user template fragment"""
        points = context.get("points", [])
        bonds_info = context.get("bonds_info", [])
        atoms_data = context.get("atoms_data", [])
        attachment_atom = context.get("attachment_atom")

        if not points or not atoms_data:
            return

        # Create atoms
        atom_id_map = {}  # template id -> scene atom id

        for i, (pos, atom_data) in enumerate(zip(points, atoms_data)):
            # Skip first atom if attaching to existing atom
            if i == 0 and attachment_atom:
                atom_id_map[atom_data["id"]] = attachment_atom.atom_id
                continue

            symbol = atom_data.get("symbol", "C")
            charge = atom_data.get("charge", 0)
            radical = atom_data.get("radical", 0)

            atom_id = self.data.add_atom(symbol, pos, charge, radical)
            atom_id_map[atom_data["id"]] = atom_id

            # Create visual atom item
            atom_item = AtomItem(atom_id, symbol, pos, charge, radical)
            self.data.atoms[atom_id]["item"] = atom_item
            self.addItem(atom_item)

        # Create bonds (bonds_info is always id-based)
        # Create index-to-id conversion table first
        index_to_id = [atom_data.get("id", i) for i, atom_data in enumerate(atoms_data)]
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                # Convert 1st/2nd elements of bonds_info to id if they are indices
                atom1_idx = bond_info[0]
                atom2_idx = bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                stereo = bond_info[3] if len(bond_info) > 3 else 0

                # Index-to-id conversion (leave as-is if already id)
                if isinstance(atom1_idx, int) and atom1_idx < len(index_to_id):
                    template_atom1_id = index_to_id[atom1_idx]
                else:
                    template_atom1_id = atom1_idx
                if isinstance(atom2_idx, int) and atom2_idx < len(index_to_id):
                    template_atom2_id = index_to_id[atom2_idx]
                else:
                    template_atom2_id = atom2_idx

                atom1_id = atom_id_map.get(template_atom1_id)
                atom2_id = atom_id_map.get(template_atom2_id)

                if atom1_id is not None and atom2_id is not None:
                    # Skip if bond already exists
                    existing_bond = None
                    if (atom1_id, atom2_id) in self.data.bonds:
                        existing_bond = (atom1_id, atom2_id)
                    elif (atom2_id, atom1_id) in self.data.bonds:
                        existing_bond = (atom2_id, atom1_id)

                    if not existing_bond:
                        bond_key, _ = self.data.add_bond(
                            atom1_id, atom2_id, order, stereo
                        )
                        # Create visual bond item
                        atom1_item = self.data.atoms[atom1_id]["item"]
                        atom2_item = self.data.atoms[atom2_id]["item"]
                        if atom1_item and atom2_item:
                            bond_item = BondItem(atom1_item, atom2_item, order, stereo)
                            self.data.bonds[bond_key]["item"] = bond_item
                            self.addItem(bond_item)
                            atom1_item.bonds.append(bond_item)
                            atom2_item.bonds.append(bond_item)

        # Update atom visuals
        for atom_id in atom_id_map.values():
            if atom_id in self.data.atoms and self.data.atoms[atom_id]["item"]:
                self.data.atoms[atom_id]["item"].update_style()

    def update_user_template_preview(self, pos):
        """Update user template preview"""
        # Robust preview: avoid self.data.atoms for preview-only atoms
        if not hasattr(self, "user_template_data") or not self.user_template_data:
            return

        template_data = self.user_template_data
        atoms = template_data.get("atoms", [])
        bonds = template_data.get("bonds", [])

        if not atoms:
            return

        # Find attachment point (first atom or clicked item)
        items_under = self.items(pos)
        attachment_atom = None
        for item in items_under:
            if isinstance(item, AtomItem):
                attachment_atom = item
                break

        # Calculate template positions
        points = []
        # Find template bounds for centering
        center_x = 0.0
        center_y = 0.0
        if atoms:
            min_x = min(atom["x"] for atom in atoms)
            max_x = max(atom["x"] for atom in atoms)
            min_y = min(atom["y"] for atom in atoms)
            max_y = max(atom["y"] for atom in atoms)
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
        # Position template
        if attachment_atom:
            # Attach to existing atom
            attach_pos = attachment_atom.pos()
            offset_x = attach_pos.x() - atoms[0]["x"]
            offset_y = attach_pos.y() - atoms[0]["y"]
        else:
            # Center at cursor position
            offset_x = pos.x() - center_x
            offset_y = pos.y() - center_y
        # Calculate atom positions
        for atom in atoms:
            new_pos = QPointF(atom["x"] + offset_x, atom["y"] + offset_y)
            points.append(new_pos)
        # Create atom ID to index mapping (for preview only)
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get("id", i)
            atom_id_to_index[atom_id] = i
        # Generate bonds_info from template bonds
        bonds_info = []
        for bond in bonds:
            atom1_idx = atom_id_to_index.get(bond["atom1"])
            atom2_idx = atom_id_to_index.get(bond["atom2"])
            if atom1_idx is not None and atom2_idx is not None:
                order = bond.get("order", 1)
                stereo = bond.get("stereo", 0)
                bonds_info.append((atom1_idx, atom2_idx, order, stereo))
        # Preview: draw lines from points and bonds_info
        # Save placement context
        self.template_context = {
            "points": points,
            "bonds_info": bonds_info,
            "atoms_data": atoms,
            "attachment_atom": attachment_atom,
        }
        # Clear legacy preview items
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(
                item, "_is_template_preview", False
            ):
                self.removeItem(item)

        # Draw preview using TemplatePreviewItem
        self.template_preview.set_user_template_geometry(points, bonds_info, atoms)
        self.template_preview.show()
        if self.views():
            self.views()[0].viewport().update()


class KeyboardMixin:
    """
    Mixin class that handles all keyboard events for MoleculeScene.
    Because this is a Mixin, `self` refers directly to the MoleculeScene instance.
    """

    def _calculate_new_atom_position(self, start_atom, bond_length):
        """
        Calculate the position for a new atom based on the surroundings of start_atom.
        Returns the offset QPointF.
        """
        start_pos = start_atom.pos()
        l = bond_length
        new_pos_offset = QPointF(0, -l)  # Default offset (up)

        # Get non-H neighbors
        neighbor_positions = []
        for bond in start_atom.bonds:
            other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
            if hasattr(other_atom, "symbol") and other_atom.symbol != "H":  # Ignore H
                neighbor_positions.append(other_atom.pos())

        num_non_H_neighbors = len(neighbor_positions)

        if num_non_H_neighbors == 0:
            # Zero bonds: default direction (up)
            new_pos_offset = QPointF(0, -l)

        elif num_non_H_neighbors == 1:
            # One bond: ~120/60 degree angle
            bond = start_atom.bonds[0]
            other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
            existing_bond_vector = start_pos - other_atom.pos()

            # Rotate 60° clockwise from existing bond
            angle_rad = math.radians(60)
            cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
            vx, vy = existing_bond_vector.x(), existing_bond_vector.y()
            new_vx, new_vy = vx * cos_a - vy * sin_a, vx * sin_a + vy * cos_a
            rotated_vector = QPointF(new_vx, new_vy)
            line = QLineF(QPointF(0, 0), rotated_vector)
            line.setLength(l)
            new_pos_offset = line.p2()

        elif num_non_H_neighbors == 3:
            bond_vectors_sum = QPointF(0, 0)
            for pos in neighbor_positions:
                # Vector from start_pos to neighbor_pos
                vec = pos - start_pos
                # Convert to unit vector
                line_to_other = QLineF(QPointF(0, 0), vec)
                if line_to_other.length() > 0:
                    line_to_other.setLength(1.0)
                    bond_vectors_sum += line_to_other.p2()

            # SUM_TOLERANCE is now a module-level constant
            if bond_vectors_sum.manhattanLength() > SUM_TOLERANCE:
                new_direction_line = QLineF(QPointF(0, 0), -bond_vectors_sum)
                new_direction_line.setLength(l)
                new_pos_offset = new_direction_line.p2()
            else:
                new_pos_offset = QPointF(l * 0.7071, -l * 0.7071)

        else:  # 2, 4+ bonds: skeleton continuation or over-bonding
            bond_vectors_sum = QPointF(0, 0)
            for bond in start_atom.bonds:
                other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                line_to_other = QLineF(start_pos, other_atom.pos())
                if line_to_other.length() > 0:
                    line_to_other.setLength(1.0)
                    bond_vectors_sum += line_to_other.p2() - line_to_other.p1()

            if bond_vectors_sum.manhattanLength() > 0.01:
                new_direction_line = QLineF(QPointF(0, 0), -bond_vectors_sum)
                new_direction_line.setLength(l)
                new_pos_offset = new_direction_line.p2()
            else:
                # Default (up) if sum is zero
                new_pos_offset = QPointF(0, -l)

        return new_pos_offset

    def keyPressEvent(self, event: Any) -> None:
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())
        key = event.key()
        modifiers = event.modifiers()

        if not self.window.is_2d_editable:
            return

        try:
            if key == Qt.Key.Key_4:
                # Case 1: Cursor over atom/bond (one-shot placement)
                if isinstance(item_at_cursor, (AtomItem, BondItem)):
                    # Set benzene template parameters
                    n, is_aromatic = 6, True
                    points, bonds_info, existing_items = [], [], []

                    # Calculate placement like update_template_preview
                    if isinstance(item_at_cursor, AtomItem):
                        p0 = item_at_cursor.pos()
                        l = DEFAULT_BOND_LENGTH
                        direction = QLineF(p0, cursor_pos).unitVector()
                        p1 = (
                            p0 + direction.p2() * l
                            if direction.length() > 0
                            else p0 + QPointF(l, 0)
                        )
                        points = self._calculate_polygon_from_edge(p0, p1, n)
                        existing_items = [item_at_cursor]

                    elif isinstance(item_at_cursor, BondItem):
                        p0, p1 = item_at_cursor.atom1.pos(), item_at_cursor.atom2.pos()
                        points = self._calculate_polygon_from_edge(
                            p0, p1, n, cursor_pos=cursor_pos, use_existing_length=True
                        )
                        existing_items = [item_at_cursor.atom1, item_at_cursor.atom2]

                    if points:
                        bonds_info = [
                            (i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)
                        ]

                        # Add fragment at calculated position
                        self.add_molecule_fragment(
                            points, bonds_info, existing_items=existing_items
                        )
                        self.update_all_items()
                        self.window.push_undo_state()

                # Case 2: Cursor over empty space (mode switch)
                else:
                    self.window.set_mode_and_update_toolbar("template_benzene")

                event.accept()
                return

            # --- 0a. Change radical (.) ---
            if key == Qt.Key.Key_Period:
                target_atoms = []
                selected = self.selectedItems()
                if selected:
                    target_atoms = [
                        item for item in selected if isinstance(item, AtomItem)
                    ]
                elif isinstance(item_at_cursor, AtomItem):
                    target_atoms = [item_at_cursor]

                if target_atoms:
                    for atom in target_atoms:
                        # Toggle radical state (0 -> 1 -> 2 -> 0)
                        atom.prepareGeometryChange()
                        atom.radical = (atom.radical + 1) % 3
                        self.data.atoms[atom.atom_id]["radical"] = atom.radical
                        atom.update_style()
                    self.update_all_items()
                    self.window.push_undo_state()
                    event.accept()
                    return

            # --- 0b. Change charge (+/-) ---
            if key == Qt.Key.Key_Plus or key == Qt.Key.Key_Minus:
                target_atoms = []
                selected = self.selectedItems()
                if selected:
                    target_atoms = [
                        item for item in selected if isinstance(item, AtomItem)
                    ]
                elif isinstance(item_at_cursor, AtomItem):
                    target_atoms = [item_at_cursor]

                if target_atoms:
                    delta = 1 if key == Qt.Key.Key_Plus else -1
                    for atom in target_atoms:
                        atom.prepareGeometryChange()
                        atom.charge += delta
                        self.data.atoms[atom.atom_id]["charge"] = atom.charge
                        atom.update_style()
                    self.update_all_items()
                    self.window.push_undo_state()
                    event.accept()
                    return

            # --- 1. Atom operations (change symbol) ---
            if isinstance(item_at_cursor, AtomItem):
                new_symbol = None
                if (
                    modifiers == Qt.KeyboardModifier.NoModifier
                    and key in self.key_to_symbol_map
                ):
                    new_symbol = self.key_to_symbol_map[key]
                elif (
                    modifiers == Qt.KeyboardModifier.ShiftModifier
                    and key in self.key_to_symbol_map_shift
                ):
                    new_symbol = self.key_to_symbol_map_shift[key]

                if new_symbol and item_at_cursor.symbol != new_symbol:
                    try:
                        item_at_cursor.prepareGeometryChange()

                        item_at_cursor.symbol = new_symbol
                        self.data.atoms[item_at_cursor.atom_id]["symbol"] = new_symbol
                        item_at_cursor.update_style()

                        atoms_to_update = {item_at_cursor}
                        for bond in item_at_cursor.bonds:
                            bond.update()
                            other_atom = (
                                bond.atom1
                                if bond.atom2 is item_at_cursor
                                else bond.atom2
                            )
                            atoms_to_update.add(other_atom)

                        for atom in atoms_to_update:
                            atom.update_style()

                        self.update_all_items()
                        self.window.push_undo_state()
                        event.accept()
                        return
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.error(
                            f"Error changing atom symbol via key {key}: {e}",
                            exc_info=True,
                        )
                        pass

            # --- 2. Bond operations (change order/stereo) ---
            target_bonds = []
            if isinstance(item_at_cursor, BondItem):
                target_bonds = [item_at_cursor]
            else:
                target_bonds = [
                    it for it in self.selectedItems() if isinstance(it, BondItem)
                ]

            if target_bonds:
                any_bond_changed = False
                for bond in target_bonds:
                    # 1. Identify current bond key in data model
                    id1, id2 = bond.atom1.atom_id, bond.atom2.atom_id
                    current_key = None
                    if (id1, id2) in self.data.bonds:
                        current_key = (id1, id2)
                    elif (id2, id1) in self.data.bonds:
                        current_key = (id2, id1)

                    if not current_key:
                        continue

                    # 2. Save previous state
                    old_order, old_stereo = bond.order, bond.stereo

                    # 3. Update BondItem properties based on key
                    if key == Qt.Key.Key_W:
                        if bond.stereo == 1:
                            bond_data = self.data.bonds.pop(current_key)
                            new_key = (current_key[1], current_key[0])
                            self.data.bonds[new_key] = bond_data
                            bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                            bond.update_position()
                            was_reversed = True
                        else:
                            bond.order = 1
                            bond.stereo = 1

                    elif key == Qt.Key.Key_D:
                        if bond.stereo == 2:
                            bond_data = self.data.bonds.pop(current_key)
                            new_key = (current_key[1], current_key[0])
                            self.data.bonds[new_key] = bond_data
                            bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                            bond.update_position()
                            was_reversed = True
                        else:
                            bond.order = 1
                            bond.stereo = 2

                    elif key == Qt.Key.Key_1 and (bond.order != 1 or bond.stereo != 0):
                        bond.order = 1
                        bond.stereo = 0
                    elif key == Qt.Key.Key_2 and (bond.order != 2 or bond.stereo != 0):
                        bond.order = 2
                        bond.stereo = 0
                    elif key == Qt.Key.Key_3 and bond.order != 3:
                        bond.order = 3
                        bond.stereo = 0

                    # 4. Update data model if changed
                    if old_order != bond.order or old_stereo != bond.stereo:
                        any_bond_changed = True

                        # 5. Remove data with old key
                        bond_data = self.data.bonds.pop(current_key)
                        bond_data["order"] = bond.order
                        bond_data["stereo"] = bond.stereo

                        # 6. Determine new key and re-register
                        new_key_id1, new_key_id2 = (
                            bond.atom1.atom_id,
                            bond.atom2.atom_id,
                        )
                        if bond.stereo == 0:
                            if new_key_id1 > new_key_id2:
                                new_key_id1, new_key_id2 = new_key_id2, new_key_id1

                        new_key = (new_key_id1, new_key_id2)
                        self.data.bonds[new_key] = bond_data

                        bond.update()

                if any_bond_changed:
                    self.update_all_items()
                    self.window.push_undo_state()

                if key in [
                    Qt.Key.Key_1,
                    Qt.Key.Key_2,
                    Qt.Key.Key_3,
                    Qt.Key.Key_W,
                    Qt.Key.Key_D,
                ]:
                    event.accept()
                    return

            if isinstance(self.hovered_item, BondItem) and self.hovered_item.order == 2:
                if event.key() == Qt.Key.Key_Z:
                    self.update_bond_stereo(self.hovered_item, 3)  # Z-isomer
                    self.update_all_items()
                    self.window.push_undo_state()
                    event.accept()
                    return
                elif event.key() == Qt.Key.Key_E:
                    self.update_bond_stereo(self.hovered_item, 4)  # E-isomer
                    self.update_all_items()
                    self.window.push_undo_state()
                    event.accept()
                    return

            # --- 3. Atom operations (add atom) ---
            if key in [Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3]:
                target_order = 1
                if key == Qt.Key.Key_2:
                    target_order = 2
                elif key == Qt.Key.Key_3:
                    target_order = 3

                start_atom = None
                if isinstance(item_at_cursor, AtomItem):
                    start_atom = item_at_cursor
                else:
                    selected_atoms = [
                        item
                        for item in self.selectedItems()
                        if isinstance(item, AtomItem)
                    ]
                    if len(selected_atoms) == 1:
                        start_atom = selected_atoms[0]

                if start_atom:
                    start_pos = start_atom.pos()
                    l = DEFAULT_BOND_LENGTH
                    new_pos_offset = self._calculate_new_atom_position(start_atom, l)

                    # SNAP_DISTANCE is a module-level constant
                    target_pos = start_pos + new_pos_offset

                    # Find nearby atom
                    near_atom = self.find_atom_near(target_pos, tol=SNAP_DISTANCE)

                    if near_atom and near_atom is not start_atom:
                        # Bond if exists
                        self.create_bond(
                            start_atom,
                            near_atom,
                            bond_order=target_order,
                            bond_stereo=0,
                        )
                    else:
                        # Create new atom and bond
                        new_atom_id = self.create_atom("C", target_pos)
                        new_atom_item = self.data.atoms[new_atom_id]["item"]
                        self.create_bond(
                            start_atom,
                            new_atom_item,
                            bond_order=target_order,
                            bond_stereo=0,
                        )

                    self.clearSelection()
                    self.update_all_items()
                    self.window.push_undo_state()
                    event.accept()
                    return

            # --- 4. Global operations (delete, mode switch) ---
            if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
                if self.temp_line:
                    try:
                        if not sip_isdeleted_safe(self.temp_line):
                            try:
                                if (
                                    getattr(self.temp_line, "scene", None)
                                    and self.temp_line.scene()
                                ):
                                    self.removeItem(self.temp_line)
                            except (
                                AttributeError,
                                RuntimeError,
                                ValueError,
                                TypeError,
                            ) as e:
                                logging.debug(
                                    f"Suppressed exception: {e}"
                                )  # Suppress bond visual state sync errors
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        try:
                            self.removeItem(self.temp_line)
                        except (
                            AttributeError,
                            RuntimeError,
                            ValueError,
                            TypeError,
                        ) as e:
                            logging.debug(
                                f"Suppressed exception: {e}"
                            )  # Suppress atom visual state sync errors

                    self.temp_line = None
                    self.start_atom = None
                    self.start_pos = None
                    self.initial_positions_in_event = {}
                    event.accept()
                    return

                items_to_process = set(self.selectedItems())
                # Include item under cursor in deletion
                if item_at_cursor and isinstance(item_at_cursor, (AtomItem, BondItem)):
                    items_to_process.add(item_at_cursor)

                if self.delete_items(items_to_process):
                    self.update_all_items()
                    self.window.push_undo_state()
                    self.window.statusBar().showMessage("Deleted selected items.")

                # Clear scene if no atoms left
                if not self.data.atoms:
                    # 1. Remove all graphics items
                    self.clear()

                    # 2. Re-initialize required items
                    self.reinitialize_items()

                    # 3. Reset temporary states
                    self.temp_line = None
                    self.start_atom = None
                    self.start_pos = None
                    self.initial_positions_in_event = {}

                    # Event handled
                    event.accept()
                    return

                # Force redraw
                if self.views():
                    self.views()[0].viewport().update()
                    QApplication.processEvents()

                    event.accept()
                    return

            if key == Qt.Key.Key_Space:
                if self.mode != "select":
                    self.window.activate_select_mode()
                else:
                    self.window.select_all()
                event.accept()
                return

            # Global drawing mode switch
            mode_to_set = None

            # 1. Switch to atom mode
            symbol_for_mode_change = None
            if (
                modifiers == Qt.KeyboardModifier.NoModifier
                and key in self.key_to_symbol_map
            ):
                symbol_for_mode_change = self.key_to_symbol_map[key]
            elif (
                modifiers == Qt.KeyboardModifier.ShiftModifier
                and key in self.key_to_symbol_map_shift
            ):
                symbol_for_mode_change = self.key_to_symbol_map_shift[key]

            if symbol_for_mode_change:
                mode_to_set = f"atom_{symbol_for_mode_change}"

            # 2. Switch to bond mode
            elif (
                modifiers == Qt.KeyboardModifier.NoModifier
                and key in self.key_to_bond_mode_map
            ):
                bond_data = self.key_to_bond_mode_map[key]
                mode_to_set = f"bond_{bond_data[0]}_{bond_data[1]}"

            # Execute mode change
            if mode_to_set:
                if hasattr(self.window, "set_mode_and_update_toolbar"):
                    self.window.set_mode_and_update_toolbar(mode_to_set)
                    event.accept()
                    return

            # Fallback
            super().keyPressEvent(event)

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.error(
                f"Unexpected error in MoleculeScene.keyPressEvent: {e}", exc_info=True
            )
            event.ignore()


class SceneQueryMixin:
    """
    Mixin class for spatial queries and basic item lifecycle.
    """

    # -------------------------------------------------------------------------
    # CUT AND PASTE THE FOLLOWING METHODS FROM molecule_scene.py HERE:
    # -------------------------------------------------------------------------
    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self.data.add_atom(symbol, pos, charge=charge, radical=radical)
        atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
        self.data.atoms[atom_id]["item"] = atom_item
        self.addItem(atom_item)
        return atom_id

    def create_bond(self, start_atom, end_atom, bond_order=None, bond_stereo=None):
        if start_atom is None or end_atom is None:
            logging.error("Error: Cannot create bond with None atoms")
            return

        try:
            exist_b = self.find_bond_between(start_atom, end_atom)
            if exist_b:
                return

            # Use the order specified in the argument if provided, otherwise use current defaults
            order_to_use = self.bond_order if bond_order is None else bond_order
            stereo_to_use = self.bond_stereo if bond_stereo is None else bond_stereo

            key, status = self.data.add_bond(
                start_atom.atom_id, end_atom.atom_id, order_to_use, stereo_to_use
            )
            if status == "created":
                bond_item = BondItem(start_atom, end_atom, order_to_use, stereo_to_use)
                self.data.bonds[key]["item"] = bond_item
                if hasattr(start_atom, "bonds"):
                    start_atom.bonds.append(bond_item)
                if hasattr(end_atom, "bonds"):
                    end_atom.bonds.append(bond_item)
                self.addItem(bond_item)

            if hasattr(start_atom, "update_style"):
                start_atom.update_style()
            if hasattr(end_atom, "update_style"):
                end_atom.update_style()

        except (AttributeError, RuntimeError, ValueError) as e:
            logging.error(f"Error creating bond: {e}", exc_info=True)
            self.update_all_items()

    def delete_items(self, items_to_delete):
        """Safely delete specified items (atoms/bonds) in order"""
        if not items_to_delete:
            return False

        # 1. Sanitize input: only keep live QGraphicsItem wrappers
        sanitized = set()
        for it in items_to_delete:
            if it is not None and not sip_isdeleted_safe(it):
                if isinstance(it, (AtomItem, BondItem, QGraphicsItem)):
                    sanitized.add(it)
        items_to_delete = sanitized

        try:
            atoms_to_delete = {it for it in items_to_delete if isinstance(it, AtomItem)}
            bonds_to_delete = {it for it in items_to_delete if isinstance(it, BondItem)}

            # Include bonds attached to atoms being deleted
            for atom in list(atoms_to_delete):
                if hasattr(atom, "bonds") and atom.bonds:
                    for b in list(atom.bonds):
                        bonds_to_delete.add(b)

            # Determine surviving atoms whose bond lists need pruning
            atoms_to_update = set()
            for bond in list(bonds_to_delete):
                a1 = getattr(bond, "atom1", None)
                a2 = getattr(bond, "atom2", None)
                if a1 and a1 not in atoms_to_delete:
                    atoms_to_update.add(a1)
                if a2 and a2 not in atoms_to_delete:
                    atoms_to_update.add(a2)

            # 2. Update surviving atoms' bond lists
            for atom in list(atoms_to_update):
                if sip_isdeleted_safe(atom):
                    continue
                if hasattr(atom, "bonds") and atom.bonds:
                    try:
                        # Filter out bonds being deleted and SIP-stale wrappers
                        atom.bonds[:] = [
                            b
                            for b in atom.bonds
                            if not sip_isdeleted_safe(b) and b not in bonds_to_delete
                        ]
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress SIP-stale bond removal errors

                if hasattr(atom, "update_style"):
                    atom.update_style()

            # 3. Remove from data model
            for bond in list(bonds_to_delete):
                a1 = getattr(bond, "atom1", None)
                a2 = getattr(bond, "atom2", None)
                if a1 and a2 and hasattr(self, "data"):
                    try:
                        # Try both directions
                        if not self.data.remove_bond(a1.atom_id, a2.atom_id):
                            self.data.remove_bond(a2.atom_id, a1.atom_id)
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress bond model data removal errors

            for atom in list(atoms_to_delete):
                if hasattr(atom, "atom_id") and hasattr(self, "data"):
                    try:
                        self.data.remove_atom(atom.atom_id)
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress data model/graphic removal errors

            try:
                self._ih_update_counter = getattr(self, "_ih_update_counter", 0) + 1
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Suppress non-critical internal update counter increment errors.
                # This counter is only used for UI throttling and is non-critical for data integrity.
                pass
            # 4. Remove graphic items from the scene
            current_scene_items = set(self.items())

            # Helper to safely remove and hide items
            def safe_remove_and_hide(item_set):
                for item in list(item_set):
                    if sip_isdeleted_safe(item):
                        continue
                    if item in current_scene_items:
                        try:
                            self.removeItem(item)
                        except (RuntimeError, ValueError, TypeError) as e:
                            logging.debug(
                                f"Suppressed exception: {e}"
                            )  # Suppress graphic removal errors (stale pointers)

                    try:
                        item.hide()
                        if (
                            not hasattr(self, "_deleted_items")
                            or self._deleted_items is None
                        ):
                            self._deleted_items = []
                        self._deleted_items.append(item)
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress data model/graphic removal errors

            safe_remove_and_hide(bonds_to_delete)
            safe_remove_and_hide(atoms_to_delete)

            for atom in list(atoms_to_update):
                if hasattr(atom, "update_style"):
                    atom.update_style()

            self.update_all_items()
            return True

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.error(f"Error during delete_items operation: {e}", exc_info=True)
            self.update_all_items()
            return False

    def find_atom_near(self, pos, tol=14.0):
        # Create a small search rectangle around the position
        search_rect = QRectF(pos.x() - tol, pos.y() - tol, 2 * tol, 2 * tol)
        nearby_items = self.items(search_rect)

        for it in nearby_items:
            if isinstance(it, AtomItem):
                # Check the precise distance only for candidate items
                if QLineF(it.pos(), pos).length() <= tol:
                    return it
        return None

    def find_bond_between(self, atom1, atom2):
        for b in atom1.bonds:
            if (b.atom1 is atom1 and b.atom2 is atom2) or (
                b.atom1 is atom2 and b.atom2 is atom1
            ):
                return b
        return None

    def update_bond_stereo(self, bond_item, new_stereo):
        """Update bond stereochemistry"""
        if bond_item is None:
            return

        try:
            if bond_item.order != 2 or bond_item.stereo == new_stereo:
                return

            a1 = getattr(bond_item, "atom1", None)
            a2 = getattr(bond_item, "atom2", None)
            if (
                not a1
                or not a2
                or not hasattr(a1, "atom_id")
                or not hasattr(a2, "atom_id")
            ):
                return

            id1, id2 = a1.atom_id, a2.atom_id
            key_to_update = (id1, id2)

            if key_to_update not in self.data.bonds:
                key_to_update = (id2, id1)
                if key_to_update not in self.data.bonds:
                    if hasattr(self.window, "statusBar"):
                        self.window.statusBar().showMessage(
                            f"Warning: Bond {id1}-{id2} not found in model.", 3000
                        )
                    return

            # Update data model and visual representation
            self.data.bonds[key_to_update]["stereo"] = new_stereo
            bond_item.set_stereo(new_stereo)
            self.data_changed_in_event = True

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.error(
                f"Error updating bond stereo for bond {bond_item}: {e}", exc_info=True
            )
            if hasattr(self.window, "statusBar"):
                self.window.statusBar().showMessage(f"Error: {e}", 5000)
            self.update_all_items()
