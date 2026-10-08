"""Regression coverage for restored connectivity and scene reset ownership."""

import pytest
import copy
from PyQt6.QtCore import QPointF
from PyQt6 import sip

from moleditpy.core.molecular_data import MolecularData
from moleditpy.ui.molecule_scene import MoleculeScene


@pytest.mark.parametrize("json_format", [False, True])
def test_restored_connectivity_supports_new_bonds(app, json_format):
    """Both restoration formats retain directional bonds and isolated atoms."""
    scene = MoleculeScene(MolecularData(), None)
    if json_format:
        scene.restore_atoms_and_bonds_from_json(
            [{"id": i, "symbol": "C", "x": i * 75, "y": 0} for i in range(3)],
            [{"atom1": 1, "atom2": 0, "order": 1, "stereo": 1}],
        )
    else:
        scene.restore_atoms_and_bonds(
            {i: {"symbol": "C", "pos": (i * 75, 0)} for i in range(3)},
            {(1, 0): {"order": 1, "stereo": 1}},
        )
    assert scene.data.adjacency_list == {0: [1], 1: [0], 2: []}
    scene.data.add_bond(1, 2)
    assert scene.data.adjacency_list == {0: [1], 1: [0, 2], 2: [1]}


def test_clear_drops_deleted_items_and_drag_references(app):
    """Clearing flagged atoms, then restoring fewer items, leaves no Qt zombies."""
    scene = MoleculeScene(MolecularData(), None)
    ids = [scene.create_atom("C", QPointF(i * 75, 0)) for i in range(3)]
    atoms = [scene.atom_items[i] for i in ids]
    scene.create_bond(atoms[0], atoms[1])
    bond = next(iter(scene.bond_items.values()))
    atoms[2].has_problem = True
    scene.start_atom = scene.hovered_item = scene.chain_start_atom = atoms[2]
    scene.initial_positions_in_event = {atoms[2]: atoms[2].pos()}
    scene.template_context = {"items": atoms}
    scene.clear()
    assert all(sip.isdeleted(item) for item in [*atoms, bond])
    assert scene.atom_items == scene.bond_items == {}
    assert scene.start_atom is scene.hovered_item is scene.chain_start_atom is None
    assert not scene.initial_positions_in_event
    scene.data = MolecularData()
    scene.reinitialize_items()
    scene.restore_atoms_and_bonds({0: {"symbol": "N", "pos": (0, 0)}}, {})
    assert set(scene.atom_items) == {0}
    assert not scene.bond_items
    assert scene.clear_all_problem_flags() is False


@pytest.mark.parametrize("stereo", [1, 2])
def test_benzene_fuses_on_reverse_direction_stereo_bond(app, stereo):
    """Reverse wedge/dash model keys stay consistent with the fused bond item."""
    scene = MoleculeScene(MolecularData(), None)
    scene.create_atom("C", QPointF(0, 0))
    scene.create_atom("C", QPointF(75, 0))
    a, b = scene.atom_items.values()
    scene.create_bond(b, a, bond_order=1, bond_stereo=stereo)
    points = scene._calculate_polygon_from_edge(a.pos(), b.pos(), 6)
    bonds = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    scene.add_molecule_fragment(points, bonds, [a, b])
    assert len(scene.data.atoms) == len(scene.data.bonds) == 6
    assert (1, 0) in scene.data.bonds
    for key, item in scene.bond_items.items():
        assert scene.data.bonds[key] == {"order": item.order, "stereo": item.stereo}


def test_fragment_failure_rolls_back_partial_model_and_scene(app, monkeypatch):
    """A failure after creating atoms/bonds restores the complete pre-edit state."""
    scene = MoleculeScene(MolecularData(), None)
    scene.create_atom("N", QPointF(0, 0))
    scene.atom_items[0].setSelected(True)
    before = copy.deepcopy(vars(scene.data))
    create_bond = scene.create_bond
    calls = 0

    def fail_second_bond(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("injected scene write failure")
        return create_bond(*args, **kwargs)

    monkeypatch.setattr(scene, "create_bond", fail_second_bond)
    with pytest.raises(RuntimeError, match="injected scene write failure"):
        scene.add_molecule_fragment(
            [QPointF(0, 0), QPointF(75, 0), QPointF(150, 0)],
            [(0, 1, 1), (1, 2, 1)],
            [scene.atom_items[0]],
        )
    assert vars(scene.data) == before
    assert set(scene.atom_items) == {0}
    assert not scene.bond_items
    assert scene.atom_items[0].isSelected()
