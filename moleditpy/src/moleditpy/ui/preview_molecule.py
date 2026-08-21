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
from typing import Any, Dict, List, Sequence, Tuple

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene

from .atom_item import AtomItem
from .bond_item import BondItem


class PreviewScene(QGraphicsScene):
    """Scene for preview-only atom/bond items.

    AtomItem and BondItem read their font, colour and spacing through
    scene().get_setting(), so any scene holding them has to answer it. ``host``
    is the editor scene, so previews follow the user's 2D settings.
    """

    def __init__(self, host: Any = None) -> None:
        super().__init__()
        self.host: Any = host

    def get_setting(self, key: str, default: Any = None) -> Any:
        """Forward a settings lookup to the editor scene."""
        host = self.host
        if host is None or not hasattr(host, "get_setting"):
            return default
        try:
            return host.get_setting(key, default)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return default

    def views(self) -> List[Any]:
        """Report the editor's views so zoom-dependent radii match the real items."""
        host = self.host
        if host is not None and hasattr(host, "views"):
            try:
                return list(host.views())
            except (AttributeError, RuntimeError, TypeError):
                return []
        return list(super().views())


def ensure_preview_settings(scene: Any, host: Any = None) -> None:
    """Give a plain QGraphicsScene the get_setting() the preview items call."""
    if hasattr(scene, "get_setting"):
        return

    def get_setting(key: str, default: Any = None) -> Any:
        if host is None or not hasattr(host, "get_setting"):
            return default
        try:
            # The host outlives the preview normally, but a closed editor would
            # otherwise raise out of paint()
            return host.get_setting(key, default)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return default

    scene.get_setting = get_setting


def build_preview_items(
    atoms: Sequence[Dict[str, Any]],
    bonds: Sequence[Sequence[int]],
) -> Tuple[List[AtomItem], List[BondItem], Dict[int, Tuple[int, ...]]]:
    """Build the atom/bond items that draw a template the way the editor draws it.

    ``atoms`` entries need ``pos`` plus the usual ``symbol``/``charge``/``radical``;
    ``bonds`` are ``(index1, index2[, order[, stereo]])``. The items are returned
    unparented — the caller decides which scene they belong to — together with the
    ring each bond belongs to, so a moved preview can recompute its ring centres.
    """
    atom_items: List[AtomItem] = []
    for i, atom_data in enumerate(atoms):
        item = AtomItem(
            i,
            atom_data.get("symbol", "C"),
            QPointF(atom_data.get("pos", QPointF())),
            charge=int(atom_data.get("charge", 0) or 0),
            radical=int(atom_data.get("radical", 0) or 0),
        )
        item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable, False)
        item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        item.setAcceptHoverEvents(False)
        atom_items.append(item)

    bond_items: List[BondItem] = []
    for bond in bonds:
        if len(bond) < 2:
            continue
        i, j = bond[0], bond[1]
        if not isinstance(i, int) or not isinstance(j, int):
            continue
        if not (0 <= i < len(atom_items) and 0 <= j < len(atom_items)) or i == j:
            continue
        order = int(bond[2]) if len(bond) > 2 else 1
        stereo = int(bond[3]) if len(bond) > 3 else 0
        atom1, atom2 = atom_items[i], atom_items[j]
        bond_item = BondItem(atom1, atom2, order, stereo)
        bond_item.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        bond_item.setAcceptHoverEvents(False)
        # create_bond() does this too; without it carbons never turn invisible
        atom1.bonds.append(bond_item)
        atom2.bonds.append(bond_item)
        bond_items.append(bond_item)

    ring_atoms = apply_preview_topology(atom_items, bond_items)
    for item in atom_items:
        item.update_style()
    return atom_items, bond_items, ring_atoms


def preview_content_rect(
    atom_items: Sequence[AtomItem], bond_items: Sequence[BondItem]
) -> QRectF:
    """Return the rect the preview actually draws into, without any hit padding."""
    rect = QRectF()
    for atom in atom_items:
        rect = rect.united(atom.visual_rect().translated(atom.pos()))
    for bond in bond_items:
        line = bond.get_line_in_local_coords()
        rect = rect.united(
            QRectF(line.p1(), line.p2()).normalized().translated(bond.pos())
        )
    return rect


def apply_preview_topology(
    atom_items: Sequence[AtomItem], bond_items: Sequence[BondItem]
) -> Dict[int, Tuple[int, ...]]:
    """Fill in implicit hydrogens and ring info for the preview.

    Deliberately no sanitization: a preview must draw whatever the template says,
    so the valence is never checked and an odd one still gets its H count.
    """
    try:
        from rdkit import Chem  # pylint: disable=import-outside-toplevel

        bond_types = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.AROMATIC,
        }
        mol = Chem.RWMol()
        for item in atom_items:
            atom = Chem.Atom(item.symbol)
            atom.SetFormalCharge(item.charge)
            atom.SetNumRadicalElectrons(item.radical)
            mol.AddAtom(atom)
        index_of = {id(item): i for i, item in enumerate(atom_items)}
        for bond in bond_items:
            mol.AddBond(
                index_of[id(bond.atom1)],
                index_of[id(bond.atom2)],
                bond_types.get(bond.order, Chem.BondType.SINGLE),
            )
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
    except (ValueError, RuntimeError, TypeError, AttributeError, ImportError):
        logging.debug("Preview topology unavailable", exc_info=True)
        return {}

    for item, atom in zip(atom_items, mol.GetAtoms()):
        item.implicit_h_count = atom.GetTotalNumHs()

    ring_info = mol.GetRingInfo()
    ring_atoms: Dict[int, Tuple[int, ...]] = {}
    best_ring_size: Dict[int, int] = {}
    for atom_ring, bond_ring in zip(ring_info.AtomRings(), ring_info.BondRings()):
        positions = [atom_items[idx].pos() for idx in atom_ring]
        center = (
            sum(p.x() for p in positions) / len(positions),
            sum(p.y() for p in positions) / len(positions),
        )
        for bond_idx in bond_ring:
            if bond_idx >= len(bond_items):
                continue
            bond = bond_items[bond_idx]
            bond.is_in_ring = True
            if best_ring_size.get(bond_idx, len(atom_items) + 1) > len(atom_ring):
                bond.ring_center = center
                best_ring_size[bond_idx] = len(atom_ring)
                ring_atoms[bond_idx] = tuple(atom_ring)
    return ring_atoms
