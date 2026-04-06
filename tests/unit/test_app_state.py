"""Tests for StateManager — get_current_state, JSON serialization round-trip."""

from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.app_state import StateManager
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock


def _make_state_manager(host):
    """Create a StateManager with host properly set."""
    sm = StateManager(host)
    sm.data = host.state_manager.data  # mirror conftest setup
    return sm


def _add_atom(host, symbol, x, y, charge=0, radical=0):
    """Add a real atom to host.state_manager.data and return its id."""
    aid = host.state_manager.data.add_atom(
        symbol, QPointF(x, y), charge=charge, radical=radical
    )
    # add a mock item so JSON serialization (which accesses item.pos()) works
    host.state_manager.data.atoms[aid]["item"] = MagicMock(pos=lambda: QPointF(x, y))
    return aid


def _add_bond(host, id1, id2, order=1, stereo=0):
    host.state_manager.data.add_bond(id1, id2, order=order, stereo=stereo)


# =============================================================================
# get_current_state
# =============================================================================


def test_get_current_state_captures_atoms(mock_parser_host):
    """get_current_state should capture all atoms with correct properties."""
    sm = _make_state_manager(mock_parser_host)
    _add_atom(mock_parser_host, "C", 10, 20)
    _add_atom(mock_parser_host, "O", 30, 40, charge=-1)

    state = sm.get_current_state()

    atoms = state["atoms"]
    assert len(atoms) == 2
    symbols = {d["symbol"] for d in atoms.values()}
    assert symbols == {"C", "O"}


def test_get_current_state_captures_bonds(mock_parser_host):
    """get_current_state should capture bonds with correct order."""
    sm = _make_state_manager(mock_parser_host)
    c1 = _add_atom(mock_parser_host, "C", 0, 0)
    c2 = _add_atom(mock_parser_host, "C", 50, 0)
    _add_bond(mock_parser_host, c1, c2, order=2)

    state = sm.get_current_state()

    assert len(state["bonds"]) == 1
    bond_data = next(iter(state["bonds"].values()))
    assert bond_data["order"] == 2


def test_get_current_state_with_3d_mol(mock_parser_host):
    """State should include 3D molecule binary when current_mol is set."""
    sm = _make_state_manager(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    for i in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(i).SetIntProp("_original_atom_id", i)
    mock_parser_host.view_3d_manager.current_mol = mol

    state = sm.get_current_state()

    mol_3d_binary = state["mol_3d"]
    assert isinstance(mol_3d_binary, bytes) and len(mol_3d_binary) > 0
    restored = Chem.Mol(mol_3d_binary)
    assert restored.GetNumAtoms() == mol.GetNumAtoms()
    assert "mol_3d_atom_ids" in state


def test_get_current_state_includes_version(mock_parser_host):
    """State should always include version string."""
    sm = _make_state_manager(mock_parser_host)
    state = sm.get_current_state()
    assert "version" in state


def test_get_current_state_captures_constraints(mock_parser_host):
    """3D constraints should be serialized in JSON-safe format."""
    sm = _make_state_manager(mock_parser_host)
    mock_parser_host.edit_3d_manager.constraints_3d = [
        ("Distance", (0, 1), 1.5, 1.0e5),
        ("Angle", (0, 1, 2), 109.5, 1.0e5),
    ]

    state = sm.get_current_state()

    assert len(state["constraints_3d"]) == 2
    for c in state["constraints_3d"]:
        assert isinstance(c, list)


def test_get_current_state_empty(mock_parser_host):
    """Empty state should produce valid structure with no atoms/bonds."""
    sm = _make_state_manager(mock_parser_host)
    state = sm.get_current_state()
    assert state["atoms"] == {}
    assert state["bonds"] == {}


# =============================================================================
# JSON round-trip: create_json_data → load_from_json_data
# =============================================================================


def test_json_roundtrip_preserves_atoms(mock_parser_host):
    """Atoms should survive JSON serialization round-trip."""
    sm = _make_state_manager(mock_parser_host)
    _add_atom(mock_parser_host, "C", 10, 20)
    _add_atom(mock_parser_host, "N", 30, 40, charge=1)

    json_data = sm.create_json_data()

    mock_parser_host.data.atoms.clear()
    mock_parser_host.data.bonds.clear()
    sm.load_from_json_data(json_data)

    assert len(mock_parser_host.data.atoms) == 2
    symbols = {d["symbol"] for d in mock_parser_host.data.atoms.values()}
    assert symbols == {"C", "N"}
    charged = [
        d for d in mock_parser_host.data.atoms.values() if d.get("charge", 0) != 0
    ]
    assert len(charged) == 1
    assert charged[0]["charge"] == 1


def test_json_roundtrip_preserves_bonds(mock_parser_host):
    """Bond order should survive JSON round-trip."""
    sm = _make_state_manager(mock_parser_host)
    c1 = _add_atom(mock_parser_host, "C", 0, 0)
    c2 = _add_atom(mock_parser_host, "C", 50, 0)
    _add_bond(mock_parser_host, c1, c2, order=2)

    json_data = sm.create_json_data()

    mock_parser_host.data.atoms.clear()
    mock_parser_host.data.bonds.clear()
    sm.load_from_json_data(json_data)

    assert len(mock_parser_host.data.bonds) == 1
    bond = next(iter(mock_parser_host.data.bonds.values()))
    assert bond["order"] == 2


def test_json_roundtrip_preserves_radical(mock_parser_host):
    """Radical electrons should survive JSON round-trip."""
    sm = _make_state_manager(mock_parser_host)
    _add_atom(mock_parser_host, "C", 0, 0, radical=2)

    json_data = sm.create_json_data()
    mock_parser_host.data.atoms.clear()
    mock_parser_host.data.bonds.clear()
    sm.load_from_json_data(json_data)

    assert len(mock_parser_host.data.atoms) == 1
    assert next(iter(mock_parser_host.data.atoms.values()))["radical"] == 2
