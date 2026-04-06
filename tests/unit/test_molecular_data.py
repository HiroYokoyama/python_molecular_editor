"""Tests for MolecularData — CRUD operations, adjacency, RDKit round-trip, and value validation."""

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from moleditpy.core.molecular_data import MolecularData
from moleditpy.utils.constants import ANGSTROM_PER_PIXEL
from PyQt6.QtCore import QPointF


# =============================================================================
# Basic CRUD
# =============================================================================


def test_add_atom_returns_incrementing_ids():
    """Verify atom IDs are sequential integers starting from 0."""
    data = MolecularData()
    id0 = data.add_atom("C", QPointF(0, 0))
    id1 = data.add_atom("O", QPointF(10, 0))
    id2 = data.add_atom("N", QPointF(20, 0))
    assert (id0, id1, id2) == (0, 1, 2)
    assert len(data.atoms) == 3


def test_add_atom_stores_properties():
    """Verify all atom properties (symbol, position, charge, radical) are stored."""
    data = MolecularData()
    aid = data.add_atom("N", QPointF(5.5, -3.2), charge=-1, radical=1)
    atom = data.atoms[aid]
    assert atom["symbol"] == "N"
    assert atom["pos"][0] == pytest.approx(5.5)
    assert atom["pos"][1] == pytest.approx(-3.2)
    assert atom["charge"] == -1
    assert atom["radical"] == 1


def test_add_bond_creates_adjacency():
    """Verify bonds update the adjacency list bidirectionally."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    data.add_bond(a, b, order=1)
    assert b in data.adjacency_list[a]
    assert a in data.adjacency_list[b]


def test_add_bond_normalizes_non_stereo_key():
    """Non-stereo bonds should have normalized keys (smaller ID first)."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    key, status = data.add_bond(b, a, order=1, stereo=0)  # Reversed order
    assert key == (a, b)  # Should be normalized
    assert status == "created"


def test_add_bond_preserves_stereo_key_direction():
    """Stereo bonds should preserve the original key direction."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    key, status = data.add_bond(b, a, order=1, stereo=1)  # Wedge
    assert key == (b, a)  # NOT normalized — direction matters


def test_add_bond_update_existing():
    """Adding a bond with same atoms should update, not duplicate."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    data.add_bond(a, b, order=1)
    key, status = data.add_bond(a, b, order=2)
    assert status == "updated"
    assert data.bonds[key]["order"] == 2


def test_remove_atom_cleans_adjacency_and_bonds():
    """Removing an atom should also remove its bonds and adjacency entries."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("O", QPointF(10, 0))
    c = data.add_atom("C", QPointF(20, 0))
    data.add_bond(a, b, order=1)
    data.add_bond(b, c, order=2)

    data.remove_atom(b)

    assert b not in data.atoms
    assert b not in data.adjacency_list
    assert a not in data.adjacency_list.get(b, [])
    # Bonds involving b should be gone
    for key in data.bonds:
        assert b not in key


def test_remove_bond_cleans_adjacency():
    """Removing a bond should update adjacency but leave atoms intact."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    data.add_bond(a, b, order=1)

    data.remove_bond(a, b)

    assert len(data.bonds) == 0
    assert b not in data.adjacency_list[a]
    assert a not in data.adjacency_list[b]
    # Atoms still exist
    assert a in data.atoms and b in data.atoms


def test_remove_bond_reverse_lookup():
    """remove_bond should find bond regardless of key order."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))
    data.add_bond(a, b, order=1)

    data.remove_bond(b, a)  # Reversed order
    assert len(data.bonds) == 0


def test_remove_atom_non_existent():
    """remove_atom gracefully handles removing an atom that doesn't exist."""
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    # Should not raise exception
    data.remove_atom(999)
    assert len(data.atoms) == 1


def test_remove_bond_non_existent():
    """remove_bond gracefully handles removing a bond that doesn't exist."""
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    data.add_atom("O", QPointF(10, 0))
    # Should not raise exception
    data.remove_bond(0, 1)
    assert len(data.bonds) == 0


def test_to_mol_block_handles_sanitization_failure(monkeypatch):
    """to_mol_block falls back if RDKit molecule generation fails."""
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))

    # Mock to_rdkit_mol to simulate RDKit sanitization failure (returns None)
    monkeypatch.setattr(data, "to_rdkit_mol", lambda **kwargs: None)

    mol_block = data.to_mol_block()
    assert mol_block is not None
    assert "MoleditPy" in mol_block
    assert "V2000" in mol_block


# =============================================================================
# RDKit round-trip: compare against RDKit reference values
# =============================================================================


def test_to_rdkit_mol_ethanol():
    """Build ethanol in MolecularData and compare RDKit molecular formula."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    o = data.add_atom("O", QPointF(100, 0))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c2, o, order=1)

    mol = data.to_rdkit_mol()
    assert mol is not None

    # Compare formula against RDKit's own ethanol
    ref_mol = Chem.MolFromSmiles("CCO")
    assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(
        ref_mol
    )


def test_to_rdkit_mol_benzene():
    """Build benzene (alternating single/double) and compare against RDKit reference."""
    data = MolecularData()
    # Hexagonal positions
    import math

    positions = []
    for i in range(6):
        angle = math.pi / 3 * i
        positions.append(QPointF(50 * math.cos(angle), 50 * math.sin(angle)))

    ids = [data.add_atom("C", p) for p in positions]
    for i in range(6):
        order = 2 if i % 2 == 0 else 1
        data.add_bond(ids[i], ids[(i + 1) % 6], order=order)

    mol = data.to_rdkit_mol()
    assert mol is not None

    ref_mol = Chem.MolFromSmiles("c1ccccc1")
    # After sanitization, both should have aromatic carbons
    assert mol.GetNumAtoms() == 6
    assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(
        ref_mol
    )


def test_to_rdkit_mol_charged_atom():
    """Verify formal charges are correctly transferred to RDKit mol."""
    data = MolecularData()
    n = data.add_atom("N", QPointF(0, 0), charge=1)
    data.add_atom("H", QPointF(10, 0))
    data.add_atom("H", QPointF(0, 10))
    data.add_atom("H", QPointF(-10, 0))
    data.add_atom("H", QPointF(0, -10))
    for i in range(1, 5):
        data.add_bond(n, i, order=1)

    mol = data.to_rdkit_mol()
    assert mol is not None

    # Find N atom and check charge
    n_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == "N"]
    assert len(n_atoms) == 1
    assert n_atoms[0].GetFormalCharge() == 1


def test_to_rdkit_mol_radical():
    """Verify radical electrons are correctly transferred."""
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0), radical=1)

    mol = data.to_rdkit_mol()
    assert mol is not None
    assert mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1


def test_to_rdkit_mol_double_bond():
    """Verify double bond order is correctly transferred."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=2)

    mol = data.to_rdkit_mol()
    assert mol is not None
    assert mol.GetBondBetweenAtoms(0, 1).GetBondTypeAsDouble() == 2.0


def test_to_rdkit_mol_triple_bond():
    """Verify triple bond: acetylene C#C."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=3)

    mol = data.to_rdkit_mol()
    assert mol is not None
    assert mol.GetBondBetweenAtoms(0, 1).GetBondTypeAsDouble() == 3.0

    ref = Chem.MolFromSmiles("C#C")
    assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)


def test_to_rdkit_mol_coordinate_conversion():
    """Verify pixel-to-angstrom coordinate conversion is exact."""
    data = MolecularData()
    px, py = 200.0, -150.0
    data.add_atom("C", QPointF(px, py))

    mol = data.to_rdkit_mol()
    conf = mol.GetConformer()
    pos = conf.GetAtomPosition(0)

    assert pos.x == pytest.approx(px * ANGSTROM_PER_PIXEL, abs=1e-6)
    assert pos.y == pytest.approx(-py * ANGSTROM_PER_PIXEL, abs=1e-6)  # Y is inverted


def test_to_rdkit_mol_preserves_original_atom_id():
    """Verify _original_atom_id property is set on RDKit atoms."""
    data = MolecularData()
    id0 = data.add_atom("C", QPointF(0, 0))
    id1 = data.add_atom("O", QPointF(50, 0))
    data.add_bond(id0, id1, order=1)

    mol = data.to_rdkit_mol()
    for atom in mol.GetAtoms():
        assert atom.HasProp("_original_atom_id")
        orig_id = atom.GetIntProp("_original_atom_id")
        assert orig_id in [id0, id1]


def test_to_rdkit_mol_empty_returns_none():
    """Empty MolecularData should return None."""
    data = MolecularData()
    assert data.to_rdkit_mol() is None


def test_to_mol_block_contains_all_atoms():
    """Verify MOL block text contains all expected atoms."""
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    data.add_atom("N", QPointF(50, 0))
    data.add_atom("O", QPointF(100, 0))
    data.add_bond(0, 1, order=1)
    data.add_bond(1, 2, order=2)

    mol_block = data.to_mol_block()
    assert mol_block is not None
    assert "V2000" in mol_block or "MoleditPy" in mol_block
    # Must round-trip through RDKit and back
    ref_mol = Chem.MolFromMolBlock(mol_block, sanitize=False)
    if ref_mol:
        symbols = set(a.GetSymbol() for a in ref_mol.GetAtoms())
        assert symbols == {"C", "N", "O"}


def test_molecular_weight_matches_rdkit():
    """Build acetic acid and compare MW against RDKit reference."""
    data = MolecularData()
    # CH3COOH
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    o1 = data.add_atom("O", QPointF(100, 0))
    o2 = data.add_atom("O", QPointF(50, -50))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c2, o1, order=1)
    data.add_bond(c2, o2, order=2)

    mol = data.to_rdkit_mol()
    ref = Chem.MolFromSmiles("CC(=O)O")

    # Same heavy-atom molecular weight (no Hs)
    assert Descriptors.HeavyAtomMolWt(mol) == pytest.approx(
        Descriptors.HeavyAtomMolWt(ref), abs=0.01
    )


def test_to_template_dict():
    """Verify template serialization dictionary format and content."""
    data = MolecularData()
    c = data.add_atom("C", QPointF(1.0, 2.0), charge=1, radical=0)
    o = data.add_atom("O", QPointF(10.0, 20.0))
    data.add_bond(c, o, order=1, stereo=1)  # Wedge

    tmpl = data.to_template_dict(
        "Test Template", version="2.0", application_version="1.2.3"
    )

    assert tmpl["format"] == "PME Template"
    assert tmpl["version"] == "2.0"
    assert tmpl["application_version"] == "1.2.3"
    assert tmpl["name"] == "Test Template"
    assert "created" in tmpl

    assert len(tmpl["atoms"]) == 2
    assert tmpl["atoms"][0]["symbol"] == "C"
    assert tmpl["atoms"][0]["x"] == 1.0
    assert tmpl["atoms"][0]["y"] == 2.0
    assert tmpl["atoms"][0]["charge"] == 1

    assert len(tmpl["bonds"]) == 1
    assert tmpl["bonds"][0]["order"] == 1
    assert tmpl["bonds"][0]["stereo"] == 1
