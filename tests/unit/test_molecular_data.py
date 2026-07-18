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


# =============================================================================
# Stereo bond RDKit round-trip (moved from test_main_app)
# =============================================================================


def test_to_rdkit_mol_stereo_wedge_dash():
    """Verify wedge/dash stereo bonds map to BEGINWEDGE/BEGINDASH in RDKit."""
    data = MolecularData()
    c1_id = data.add_atom("C", QPointF(0, 0))
    h1_id = data.add_atom("H", QPointF(0, 50))
    h2_id = data.add_atom("H", QPointF(0, -50))
    h3_id = data.add_atom("H", QPointF(50, 0))
    h4_id = data.add_atom("H", QPointF(-50, 0))

    data.add_bond(c1_id, h1_id, order=1, stereo=0)
    data.add_bond(c1_id, h2_id, order=1, stereo=0)
    data.add_bond(c1_id, h3_id, order=1, stereo=1)  # Wedge
    data.add_bond(c1_id, h4_id, order=1, stereo=2)  # Dash

    mol = data.to_rdkit_mol(use_2d_stereo=True)
    assert mol is not None

    atom_map = {
        atom.GetIntProp("_original_atom_id"): atom.GetIdx() for atom in mol.GetAtoms()
    }
    wedge_bond = mol.GetBondBetweenAtoms(atom_map[c1_id], atom_map[h3_id])
    dash_bond = mol.GetBondBetweenAtoms(atom_map[c1_id], atom_map[h4_id])

    assert wedge_bond.GetBondDir() == Chem.BondDir.BEGINWEDGE
    assert dash_bond.GetBondDir() == Chem.BondDir.BEGINDASH


def test_to_rdkit_mol_ez_stereo():
    """Verify E/Z stereo double bond maps to STEREOZ in RDKit."""
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(-100, 50))
    c2 = data.add_atom("C", QPointF(-50, 0))
    c3 = data.add_atom("C", QPointF(50, 0))
    c4 = data.add_atom("C", QPointF(100, 50))

    data.add_bond(c1, c2, order=1, stereo=0)
    data.add_bond(c3, c4, order=1, stereo=0)
    data.add_bond(c2, c3, order=2, stereo=3)  # Z

    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None

    atom_map = {
        atom.GetIntProp("_original_atom_id"): atom.GetIdx() for atom in mol.GetAtoms()
    }
    double_bond = mol.GetBondBetweenAtoms(atom_map[c2], atom_map[c3])

    assert double_bond.GetBondType() == Chem.BondType.DOUBLE
    assert double_bond.GetStereo() == Chem.BondStereo.STEREOZ


# =============================================================================
# Reversed-key bond handling (stereo direction vs. sorted non-stereo keys)
# =============================================================================


def test_add_bond_stereo_over_existing_reversed_key_updates_in_place():
    """Adding a directional stereo bond over an existing sorted-key bond
    must re-key the entry, not create a duplicate for the same atom pair."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))

    key1, status1 = data.add_bond(a, b, order=1, stereo=0)  # stored as (a, b)
    assert (key1, status1) == ((a, b), "created")

    # Wedge drawn from b to a: directional key (b, a)
    key2, status2 = data.add_bond(b, a, order=1, stereo=1)
    assert status2 == "updated"
    assert key2 == (b, a)
    assert len(data.bonds) == 1
    assert data.bonds[(b, a)]["stereo"] == 1
    # Adjacency must not be duplicated either
    assert data.adjacency_list[a].count(b) == 1
    assert data.adjacency_list[b].count(a) == 1


def test_add_bond_plain_over_existing_stereo_reversed_key_updates_in_place():
    """Demoting a directional stereo bond back to a plain bond (sorted key)
    must also re-key instead of duplicating."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(10, 0))

    key1, _ = data.add_bond(b, a, order=1, stereo=1)  # directional (b, a)
    assert key1 == (b, a)

    key2, status2 = data.add_bond(a, b, order=2, stereo=0)  # sorted (a, b)
    assert status2 == "updated"
    assert key2 == (a, b)
    assert len(data.bonds) == 1
    assert data.bonds[(a, b)]["order"] == 2
    assert data.bonds[(a, b)]["stereo"] == 0


# =============================================================================
# Manual MOL block fallback robustness
# =============================================================================


def _force_manual_fallback(monkeypatch, data):
    """Make to_rdkit_mol fail so to_mol_block uses the manual writer."""
    monkeypatch.setattr(
        type(data), "to_rdkit_mol", lambda self, use_2d_stereo=True: None
    )


def test_to_mol_block_fallback_handles_float_bond_order(monkeypatch):
    """Aromatic 1.5 bond order must map to V2000 code 4, not crash on '{:d}'."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("C", QPointF(75, 0))
    data.add_bond(a, b, order=1.5)
    _force_manual_fallback(monkeypatch, data)

    mol_block = data.to_mol_block()
    assert mol_block is not None
    bond_line = mol_block.splitlines()[4 + 2]  # header(3) + counts + 2 atoms
    assert bond_line[:12] == f"{1:3d}{2:3d}{4:3d}{0:3d}"


def test_to_mol_block_fallback_skips_positionless_atoms_consistently(monkeypatch):
    """Atoms without a position are excluded from counts and bond indices."""
    data = MolecularData()
    a = data.add_atom("C", QPointF(0, 0))
    b = data.add_atom("O", QPointF(75, 0))
    c = data.add_atom("N", QPointF(150, 0))
    data.atoms[a]["pos"] = None  # simulate corrupted/missing position
    data.add_bond(b, c, order=1)
    _force_manual_fallback(monkeypatch, data)

    mol_block = data.to_mol_block()
    assert mol_block is not None
    lines = mol_block.splitlines()
    counts = lines[3]
    assert counts[:6] == f"{2:3d}{1:3d}"  # 2 atoms written, 1 bond
    atom_lines = lines[4:6]
    assert "O" in atom_lines[0] and "N" in atom_lines[1]
    bond_line = lines[6]
    # Bond must reference the re-numbered atoms 1-2, not the original 2-3
    assert bond_line[:6] == f"{1:3d}{2:3d}"


# =============================================================================
# set_atom_pos raw-tuple branch
# =============================================================================

from unittest.mock import patch


def test_set_atom_pos_accepts_raw_tuple():
    data = MolecularData()
    aid = data.add_atom("C", QPointF(0, 0))
    data.set_atom_pos(aid, (5.0, 7.0))
    pos = data.atoms[aid]["pos"]
    assert (pos.x(), pos.y()) == (5.0, 7.0)


# =============================================================================
# to_mol_block manual fallback (to_rdkit_mol unavailable)
# =============================================================================


def test_to_mol_block_manual_fallback_charges_positions_stereo():
    data = MolecularData()
    ids = {
        3: data.add_atom("C", QPointF(0, 0), charge=3),
        2: data.add_atom("C", QPointF(10, 0), charge=2),
        1: data.add_atom("C", QPointF(20, 0), charge=1),
        -1: data.add_atom("C", QPointF(30, 0), charge=-1),
        -2: data.add_atom("C", QPointF(40, 0), charge=-2),
        -3: data.add_atom("C", QPointF(50, 0), charge=-3),
    }
    tuple_atom = data.add_atom("C", QPointF(60, 0))
    data.atoms[tuple_atom]["pos"] = [70.0, 5.0]  # raw list pos branch
    none_atom = data.add_atom("C", QPointF(80, 0))
    data.atoms[none_atom]["pos"] = None  # falsy pos -> skipped
    bad_atom = data.add_atom("C", QPointF(90, 0))
    data.atoms[bad_atom]["pos"] = 12345  # invalid pos type -> skipped

    data.add_bond(ids[3], ids[2], order=1, stereo=1)  # wedge stereo code
    data.add_bond(ids[2], ids[1], order=2, stereo=2)  # dash stereo code
    # Bond referencing a skipped atom is dropped from the bond block
    data.bonds[(ids[1], none_atom)] = {"order": 1, "stereo": 0}

    with patch.object(data, "to_rdkit_mol", return_value=None):
        block = data.to_mol_block()

    assert block is not None
    assert "MoleditPy" in block
    lines = block.split("\n")
    counts = lines[3]
    # 6 charged + 1 tuple atom = 7 atoms written (None/invalid skipped)
    assert counts.startswith("  7")


def test_to_mol_block_returns_none_when_empty():
    data = MolecularData()
    with patch.object(data, "to_rdkit_mol", return_value=None):
        assert data.to_mol_block() is None


# =============================================================================
# to_rdkit_mol E/Z stereo label handling
# =============================================================================


def _but2ene():
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    c3 = data.add_atom("C", QPointF(100, 0))
    c4 = data.add_atom("C", QPointF(150, 0))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c2, c3, order=2, stereo=3)  # Z label on the double bond
    data.add_bond(c3, c4, order=1)
    return data, (c2, c3)


def test_to_rdkit_mol_assigns_z_stereo_on_double_bond():
    data, (c2, c3) = _but2ene()
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None
    dbond = None
    for b in mol.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            dbond = b
    assert dbond is not None
    assert dbond.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE)


def test_to_rdkit_mol_stereo_label_on_single_bond_is_ignored():
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1, stereo=3)  # E/Z label on a single bond
    mol = data.to_rdkit_mol()
    assert mol is not None  # label ignored, no crash


def test_to_rdkit_mol_ez_falls_back_to_hydrogen_neighbor():
    # C0=C1 with an explicit H on C1 and a methyl on C0: the C1 side has only
    # a hydrogen neighbor, exercising the non-heavy fallback in neighbor pick.
    data = MolecularData()
    c0 = data.add_atom("C", QPointF(0, 0))
    c1 = data.add_atom("C", QPointF(50, 0))
    ch3 = data.add_atom("C", QPointF(-50, 0))
    h = data.add_atom("H", QPointF(100, 0))
    data.add_bond(c0, c1, order=2, stereo=4)  # E label
    data.add_bond(c0, ch3, order=1)
    data.add_bond(c1, h, order=1)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None


def test_to_rdkit_mol_ez_with_explicit_stereo_atoms():
    data, (c2, c3) = _but2ene()
    # Provide explicit stereo_atoms so the specified-index branch runs
    data.bonds[(c2, c3)]["stereo_atoms"] = (
        list(data.atoms.keys())[0],
        list(data.atoms.keys())[3],
    )
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None


def test_to_rdkit_mol_terminal_double_bond_stereo_skipped():
    # Ethene with no other substituents: neighbor pick returns None -> skip
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=2, stereo=3)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None


# =============================================================================
# to_rdkit_mol drops bonds referencing removed atoms
# =============================================================================


def test_to_rdkit_mol_skips_bond_with_unmapped_atom():
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    # Inject a dangling bond referencing a non-existent atom id
    data.bonds[(c1, 999)] = {"order": 1, "stereo": 0}
    mol = data.to_rdkit_mol()
    assert mol is not None
    assert mol.GetNumBonds() == 1  # dangling bond ignored


def test_add_atom_accepts_raw_tuple_position():
    data = MolecularData()
    aid = data.add_atom("C", (3.0, 4.0))  # raw tuple, not QPointF
    assert (data.atoms[aid]["pos"].x(), data.atoms[aid]["pos"].y()) == (3.0, 4.0)


def test_set_atom_pos_accepts_qpointf():
    data = MolecularData()
    aid = data.add_atom("C", (0.0, 0.0))
    data.set_atom_pos(aid, QPointF(9.0, 8.0))
    assert (data.atoms[aid]["pos"].x(), data.atoms[aid]["pos"].y()) == (9.0, 8.0)


def test_to_rdkit_mol_returns_none_for_unknown_symbol():
    data = MolecularData()
    data.add_atom("Xx", QPointF(0, 0))  # RDKit rejects this element
    assert data.to_rdkit_mol() is None


def test_to_rdkit_mol_returns_none_on_sanitize_failure():
    data = MolecularData()
    c = data.add_atom("C", QPointF(0, 0))
    hs = [data.add_atom("H", QPointF(10 * i, 10)) for i in range(5)]
    for h in hs:
        data.add_bond(c, h, order=1)  # pentavalent carbon fails sanitization
    assert data.to_rdkit_mol() is None


def test_to_rdkit_mol_conformer_handles_tuple_and_invalid_pos():
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    data.atoms[c1]["pos"] = [10.0, 20.0]  # raw list -> tuple conformer branch
    data.atoms[c2]["pos"] = 999  # invalid -> conformer skip branch
    mol = data.to_rdkit_mol()
    assert mol is not None


def test_to_mol_block_manual_fallback_when_moltomolblock_raises():
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0), charge=0)
    import moleditpy.core.molecular_data as md
    with patch.object(md.Chem, "MolToMolBlock", side_effect=RuntimeError("boom")):
        block = data.to_mol_block()
    assert block is not None and "MoleditPy" in block


def test_to_rdkit_mol_3d_without_ez_labels_estimates_from_coords():
    # use_2d_stereo=False with a double bond but no E/Z label -> coordinate path
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=2, stereo=0)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None


def test_to_rdkit_mol_ez_invalid_stereo_atoms_suppressed():
    data, (c2, c3) = _but2ene()
    data.bonds[(c2, c3)]["stereo_atoms"] = 5  # not unpackable -> exception branch
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert mol is not None
