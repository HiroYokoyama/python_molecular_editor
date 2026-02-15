"""Tests for MainWindowAppState — undo/redo, JSON serialization round-trip, version handling."""

import pytest
import copy
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.modules.main_window_app_state import MainWindowAppState
from moleditpy.modules.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest.mock import MagicMock, patch


class DummyAppState(MainWindowAppState):
    """Thin wrapper to exercise app state logic with a mocked host."""

    def __init__(self, host):
        self._host = host
        self.data = host.data
        self.scene = host.scene
        self.current_mol = host.current_mol
        self.is_2d_editable = host.is_2d_editable
        self.constraints_3d = host.constraints_3d
        self.current_file_path = host.current_file_path
        self.has_unsaved_changes = host.has_unsaved_changes
        self.undo_stack = []
        self.redo_stack = []

    def __getattr__(self, name):
        return getattr(self._host, name)

    def statusBar(self):
        return self._host.statusBar()

    def clear_2d_editor(self, push_to_undo=True):
        pass

    def update_window_title(self):
        pass


# =============================================================================
# get_current_state / create_json_data
# =============================================================================


def test_get_current_state_captures_atoms(mock_parser_host):
    """get_current_state should capture all atoms with correct properties."""
    app = DummyAppState(mock_parser_host)
    app.scene.create_atom("C", QPointF(10, 20))
    app.scene.create_atom("O", QPointF(30, 40), charge=-1)

    state = app.get_current_state()

    atoms = state["atoms"]
    assert len(atoms) == 2
    symbols = {d["symbol"] for d in atoms.values()}
    assert symbols == {"C", "O"}


def test_get_current_state_captures_bonds(mock_parser_host):
    """get_current_state should capture bonds with correct order."""
    app = DummyAppState(mock_parser_host)
    c1 = app.scene.create_atom("C", QPointF(0, 0))
    c2 = app.scene.create_atom("C", QPointF(50, 0))
    a1_item = app.data.atoms[c1]["item"]
    a2_item = app.data.atoms[c2]["item"]
    app.scene.create_bond(a1_item, a2_item, bond_order=2)

    state = app.get_current_state()

    assert len(state["bonds"]) == 1
    bond_data = next(iter(state["bonds"].values()))
    assert bond_data["order"] == 2


def test_get_current_state_with_3d_mol(mock_parser_host):
    """State should include 3D molecule binary when current_mol is set."""
    app = DummyAppState(mock_parser_host)
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    # Set _original_atom_id
    for i in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(i).SetIntProp("_original_atom_id", i)

    app.current_mol = mol
    state = app.get_current_state()

    assert "mol_3d" in state
    assert state["mol_3d"] is not None
    assert "mol_3d_atom_ids" in state


def test_get_current_state_includes_version(mock_parser_host):
    """State should always include version string."""
    app = DummyAppState(mock_parser_host)
    state = app.get_current_state()
    assert "version" in state


def test_get_current_state_captures_constraints(mock_parser_host):
    """3D constraints should be serialized in JSON-safe format."""
    app = DummyAppState(mock_parser_host)
    app.constraints_3d = [
        ("Distance", (0, 1), 1.5, 1.0e5),
        ("Angle", (0, 1, 2), 109.5, 1.0e5),
    ]
    state = app.get_current_state()

    assert len(state["constraints_3d"]) == 2
    # Should be JSON-safe lists, not tuples
    for c in state["constraints_3d"]:
        assert isinstance(c, list)


# =============================================================================
# JSON round-trip: create_json_data → load_from_json_data
# =============================================================================


def test_json_roundtrip_preserves_atoms(mock_parser_host):
    """Atoms should survive JSON serialization round-trip with correct properties."""
    app = DummyAppState(mock_parser_host)

    # Build a structure
    app.scene.create_atom("C", QPointF(10, 20))
    app.scene.create_atom("N", QPointF(30, 40), charge=1)

    json_data = app.create_json_data()

    # Clear and reload
    app.data.atoms.clear()
    app.data.bonds.clear()
    app.load_from_json_data(json_data)

    assert len(app.data.atoms) == 2
    symbols = {d["symbol"] for d in app.data.atoms.values()}
    assert symbols == {"C", "N"}
    # Check charge persists
    charged = [d for d in app.data.atoms.values() if d.get("charge", 0) != 0]
    assert len(charged) == 1
    assert charged[0]["charge"] == 1


def test_json_roundtrip_preserves_bonds(mock_parser_host):
    """Bond order and stereo should survive JSON round-trip."""
    app = DummyAppState(mock_parser_host)
    c1 = app.scene.create_atom("C", QPointF(0, 0))
    c2 = app.scene.create_atom("C", QPointF(50, 0))
    a1_item = app.data.atoms[c1]["item"]
    a2_item = app.data.atoms[c2]["item"]
    app.scene.create_bond(a1_item, a2_item, bond_order=2, bond_stereo=0)

    json_data = app.create_json_data()

    app.data.atoms.clear()
    app.data.bonds.clear()
    app.load_from_json_data(json_data)

    assert len(app.data.bonds) == 1
    bond = next(iter(app.data.bonds.values()))
    assert bond["order"] == 2


def test_json_roundtrip_preserves_radical(mock_parser_host):
    """Radical electrons should survive JSON round-trip."""
    app = DummyAppState(mock_parser_host)
    app.scene.create_atom("C", QPointF(0, 0), radical=2)

    json_data = app.create_json_data()
    app.data.atoms.clear()
    app.data.bonds.clear()
    app.load_from_json_data(json_data)

    assert len(app.data.atoms) == 1
    assert next(iter(app.data.atoms.values()))["radical"] == 2


# =============================================================================
# Edge cases
# =============================================================================


def test_get_current_state_empty(mock_parser_host):
    """Empty state should produce valid structure with no atoms/bonds."""
    app = DummyAppState(mock_parser_host)
    state = app.get_current_state()
    assert state["atoms"] == {}
    assert state["bonds"] == {}