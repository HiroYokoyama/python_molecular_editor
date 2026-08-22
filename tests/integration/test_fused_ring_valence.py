"""Fusing a ring template must leave a molecule that is chemically possible."""

import math

import pytest
from PyQt6.QtCore import QPointF
from rdkit import Chem

BOND_LEN = 75.0  # the app's default 2D bond length
APOTHEM = BOND_LEN * math.sqrt(3) / 2.0
START_ANGLE = math.radians(-120.0)  # matches the ring template's start angle
STEP = math.radians(60.0)

RING_A = (0.0, 0.0)
RING_B = (112.5, -64.95)  # shares the a1-a2 edge with ring A

# Naphthalene as a valid Kekulé structure
NAPHTHALENE_BONDS = [
    ("a5", "a0", 2),
    ("a0", "a1", 1),
    ("a1", "a2", 2),
    ("a2", "a3", 1),
    ("a3", "a4", 2),
    ("a4", "a5", 1),
    ("a1", "b0", 1),
    ("b0", "b1", 2),
    ("b1", "b2", 1),
    ("b2", "b3", 2),
    ("b3", "a2", 1),
]
# Every bond on the rim, with the ring centre it belongs to
PERIPHERAL_BONDS = [
    ("a5", "a0", RING_A),
    ("a0", "a1", RING_A),
    ("a2", "a3", RING_A),
    ("a3", "a4", RING_A),
    ("a4", "a5", RING_A),
    ("a1", "b0", RING_B),
    ("b0", "b1", RING_B),
    ("b1", "b2", RING_B),
    ("b2", "b3", RING_B),
    ("b3", "a2", RING_B),
]


def _hex_points(centre):
    cx, cy = centre
    return [
        QPointF(
            cx + BOND_LEN * math.cos(START_ANGLE + i * STEP),
            cy + BOND_LEN * math.sin(START_ANGLE + i * STEP),
        )
        for i in range(6)
    ]


def _build_naphthalene(scene):
    ring_a, ring_b = _hex_points(RING_A), _hex_points(RING_B)
    coords = {
        "a0": ring_a[0],
        "a1": ring_a[1],
        "a2": ring_a[2],
        "a3": ring_a[3],
        "a4": ring_a[4],
        "a5": ring_a[5],
        "b0": ring_b[0],
        "b1": ring_b[1],
        "b2": ring_b[2],
        "b3": ring_b[3],
    }
    atoms = {n: scene.atom_items[scene.create_atom("C", p)] for n, p in coords.items()}
    for first, second, order in NAPHTHALENE_BONDS:
        scene.create_bond(atoms[first], atoms[second], bond_order=order)
    return atoms


def _fuse_benzene_onto(scene, atoms, first, second, centre):
    """Hover a benzene template over the given rim bond and place it."""
    mid = (atoms[first].pos() + atoms[second].pos()) / 2.0
    away = mid - QPointF(*centre)
    length = math.hypot(away.x(), away.y()) or 1.0
    scene.update_template_preview(
        mid + QPointF(away.x() / length * APOTHEM, away.y() / length * APOTHEM)
    )
    context = scene.template_context
    scene.add_molecule_fragment(
        list(context["points"]),
        list(context["bonds_info"]),
        context.get("items") or [],
    )


@pytest.mark.parametrize("first,second,centre", PERIPHERAL_BONDS)
def test_fusing_benzene_onto_naphthalene_stays_valid(window, first, second, centre):
    """No rim bond may give a carbon a fifth bond.

    A Kekulé ring always carries a double next to its fused bond, so fusing
    onto a single bond whose atoms each already hold a double used to hand one
    of them order 5. Those bonds are added single instead, which leaves the new
    carbons as methylene.
    """
    scene = window.init_manager.scene
    atoms = _build_naphthalene(scene)
    window.ui_manager.set_mode("template_benzene")

    _fuse_benzene_onto(scene, atoms, first, second, centre)

    overloaded = {
        atom_id: sum(b.order for b in item.bonds)
        for atom_id, item in scene.atom_items.items()
        if sum(b.order for b in item.bonds) > 4
    }
    assert not overloaded, f"hypervalent carbons after fusing on {first}-{second}"

    mol = window.state_manager.data.to_rdkit_mol()
    assert mol is not None, f"fusing on {first}-{second} produced no molecule"
    assert not Chem.DetectChemistryProblems(mol)


def test_fusing_keeps_aromatic_ring_when_it_fits(window):
    """The fallback must not fire where a full Kekulé ring is still valid."""
    scene = window.init_manager.scene
    atoms = _build_naphthalene(scene)
    window.ui_manager.set_mode("template_benzene")

    # a5-a0 is a rim double bond: the template can alternate around it
    _fuse_benzene_onto(scene, atoms, "a5", "a0", RING_A)

    added = [
        i for i in scene.atom_items if i not in {a.atom_id for a in atoms.values()}
    ]
    assert len(added) == 4, "fusing should add four new carbons"
    new_doubles = sum(
        1 for item in scene.atom_items.values() for b in item.bonds if b.order == 2
    )
    assert new_doubles > 0, "the fused ring lost all of its double bonds"


def test_fusing_onto_a_crowded_bond_stays_valid(window):
    """An upgrade must fit too, not merely avoid a neighbouring double.

    x-a-b-y with a-w and b-z leaves a and b holding three single bonds each.
    Re-ordering a-b to a double takes both to four, so the ring hanging off
    them has nowhere to go; the upgrade has to be declined.
    """
    scene = window.init_manager.scene
    ring = _hex_points(RING_A)
    first = scene.atom_items[scene.create_atom("C", ring[0])]
    second = scene.atom_items[scene.create_atom("C", ring[1])]
    subs = [
        scene.atom_items[scene.create_atom("C", ring[0] + QPointF(-BOND_LEN, 0))],
        scene.atom_items[
            scene.create_atom("C", ring[0] + QPointF(-BOND_LEN, BOND_LEN))
        ],
        scene.atom_items[scene.create_atom("C", ring[1] + QPointF(BOND_LEN, 0))],
        scene.atom_items[scene.create_atom("C", ring[1] + QPointF(BOND_LEN, BOND_LEN))],
    ]
    scene.create_bond(first, second, bond_order=1)
    scene.create_bond(first, subs[0], bond_order=1)
    scene.create_bond(first, subs[1], bond_order=1)
    scene.create_bond(second, subs[2], bond_order=1)
    scene.create_bond(second, subs[3], bond_order=1)

    window.ui_manager.set_mode("template_benzene")
    scene.update_template_preview(QPointF(*RING_A))
    context = scene.template_context
    scene.add_molecule_fragment(
        list(context["points"]),
        list(context["bonds_info"]),
        context.get("items") or [],
    )

    overloaded = {
        atom_id: sum(b.order for b in item.bonds)
        for atom_id, item in scene.atom_items.items()
        if sum(b.order for b in item.bonds) > 4
    }
    assert not overloaded, "fusing onto a fully substituted bond overfilled an atom"
    assert scene.find_bond_between(first, second).order == 1, (
        "the shared bond had no room to become a double"
    )
