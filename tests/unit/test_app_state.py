"""Tests for StateManager — get_current_state, JSON serialization round-trip."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.app_state import StateManager
from moleditpy.core.molecular_data import MolecularData
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


# =============================================================================
# Persistence tests (merged from test_app_state_persistence.py)
# =============================================================================


class DummyMainWindow(StateManager):
    def __init__(self):
        self.host = self
        self.data = MolecularData()

        self.state_manager = self
        self.init_manager = MagicMock()
        self.ui_manager = MagicMock()
        self.edit_actions_manager = MagicMock()
        self.view_3d_manager = MagicMock()
        self.edit_3d_manager = MagicMock()
        self.compute_manager = MagicMock()
        self.io_manager = MagicMock()
        self.export_manager = MagicMock()
        self.dialog_manager = MagicMock()
        self.plugin_manager = MagicMock()
        self.plugin_manager.save_handlers = {}
        self.plugin_manager.load_handlers = {}

        self.init_manager.scene = MagicMock()
        self.init_manager.scene.atom_items = {}
        self.init_manager.scene.bond_items = {}

        def mock_restore_atoms_and_bonds(raw_atoms, raw_bonds):
            for atom_id, data in raw_atoms.items():
                self.data.atoms[atom_id] = {
                    "symbol": data["symbol"],
                    "pos": tuple(data["pos"]),
                    "charge": data.get("charge", 0),
                    "radical": data.get("radical", 0),
                }
            for key_tuple, data in raw_bonds.items():
                self.data.bonds[key_tuple] = {
                    "order": data.get("order", 1),
                    "stereo": data.get("stereo", 0),
                }

        def mock_restore_atoms_and_bonds_from_json(atoms_2d, bonds_2d):
            for atom_data in atoms_2d:
                atom_id = atom_data["id"]
                self.data.atoms[atom_id] = {
                    "symbol": atom_data["symbol"],
                    "pos": (float(atom_data["x"]), float(atom_data["y"])),
                    "charge": atom_data.get("charge", 0),
                    "radical": atom_data.get("radical", 0),
                }
            for bond_data in bonds_2d:
                atom1_id = bond_data["atom1"]
                atom2_id = bond_data["atom2"]
                self.data.bonds[(atom1_id, atom2_id)] = {
                    "order": bond_data["order"],
                    "stereo": bond_data.get("stereo", 0),
                }

        self.init_manager.scene.restore_atoms_and_bonds.side_effect = (
            mock_restore_atoms_and_bonds
        )
        self.init_manager.scene.restore_atoms_and_bonds_from_json.side_effect = (
            mock_restore_atoms_and_bonds_from_json
        )
        self.init_manager.view_2d = MagicMock()
        self.init_manager.settings = MagicMock()
        self.view_3d_manager.view_3d = MagicMock()
        self.view_3d_manager.plotter = MagicMock()
        self.view_3d_manager.current_mol = None
        self.view_3d_manager.atom_positions_3d = None

        self.host.init_manager.current_file_path = None
        self.host.state_manager.has_unsaved_changes = False
        self.edit_3d_manager.constraints_3d = []
        self.is_2d_editable = True
        self.edit_actions_manager.undo_stack = []
        self.edit_actions_manager.redo_stack = []
        self.is_restoring_state = False
        self.initialization_complete = True
        self._preserved_plugin_data = {}
        self.init_manager.formula_label = MagicMock()
        self.compute_manager.last_successful_optimization_method = None
        self.edit_actions_manager.clear_2d_editor.side_effect = self._do_clear_2d_editor
        self.edit_actions_manager.update_implicit_hydrogens = MagicMock()

    def _do_clear_2d_editor(self, push_to_undo=True):
        self.data.atoms.clear()
        self.data.bonds.clear()
        self.data._next_atom_id = 0

    def statusBar(self):
        return MagicMock()

    def setWindowTitle(self, title):
        pass

    def update_window_title(self):
        pass

    def update_implicit_hydrogens(self):
        pass

    def update_chiral_labels(self):
        pass

    def update_2d_measurement_labels(self):
        pass

    def update_realtime_info(self):
        pass

    def update_undo_redo_actions(self):
        pass

    def _enable_3d_features(self, val):
        pass

    def _enable_3d_edit_actions(self, val):
        pass

    def clear_2d_editor(self, push_to_undo=True):
        self._do_clear_2d_editor(push_to_undo)

    def restore_ui_for_editing(self):
        pass

    def create_atom_id_mapping(self):
        pass

    def draw_molecule_3d(self, mol):
        pass

    def create_json_data(self):
        return super().create_json_data()

    def load_from_json_data(self, data):
        return super().load_from_json_data(data)

    def get_current_state(self):
        return super().get_current_state()

    def set_state_from_data(self, data):
        return super().set_state_from_data(data)

    def set_current_molecule(self, mol):
        self.view_3d_manager.current_mol = mol

    def set_3d_atom_positions(self, positions):
        self.view_3d_manager.atom_positions_3d = positions

    def clear_3d_view(self):
        self.view_3d_manager.current_mol = None

    def set_constraints_3d(self, constraints):
        self.edit_3d_manager.constraints_3d = constraints

    def get_constraints_3d(self):
        return self.edit_3d_manager.constraints_3d

    def set_has_unsaved_changes(self, value):
        self.state_manager.has_unsaved_changes = value

    def set_current_file_path(self, path):
        self.init_manager.current_file_path = path

    def get_current_file_path(self):
        return self.init_manager.current_file_path

    def set_is_2d_editable(self, value):
        self.ui_manager.is_2d_editable = value

    def update_status_message(self, message, timeout=0):
        pass

    def set_last_successful_optimization_method(self, method):
        self.compute_manager.last_successful_optimization_method = method

    def get_settings(self):
        return self.init_manager.settings

    @property
    def scene(self):
        return self.init_manager.scene


@pytest.fixture
def dummy_window(app):
    return DummyMainWindow()


def test_pmeprj_serialization_roundtrip(dummy_window):
    """Test full project serialization/deserialization (PMEPRJ)."""
    mw = dummy_window

    aid1 = mw.state_manager.data.add_atom("C", QPointF(10, 20))
    aid2 = mw.state_manager.data.add_atom("O", QPointF(30, 40))
    mw.state_manager.data.add_bond(aid1, aid2, order=1)

    mol = Chem.MolFromSmiles("CO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", aid1)
    mol.GetAtomWithIdx(1).SetIntProp("_original_atom_id", aid2)
    mw.view_3d_manager.current_mol = mol

    mw.edit_3d_manager.constraints_3d = [("DISTANCE", (0, 1), 1.43, 8.0e5)]
    mw.compute_manager.last_successful_optimization_method = "MMFF94s"
    mw._preserved_plugin_data = {"TestPlugin": {"val": 42}}

    json_data = mw.create_json_data()

    assert json_data["format"] == "PME Project"
    assert "2d_structure" in json_data
    assert "3d_structure" in json_data
    assert json_data["3d_structure"]["num_conformers"] == 1

    mw.clear_2d_editor(push_to_undo=False)
    mw.view_3d_manager.current_mol = None
    mw.edit_3d_manager.constraints_3d = []

    mw.load_from_json_data(json_data)

    assert len(mw.state_manager.data.atoms) == 2
    assert mw.state_manager.data.atoms[aid1]["symbol"] == "C"
    assert mw.state_manager.data.atoms[aid2]["symbol"] == "O"
    assert mw.state_manager.data.atoms[aid1]["pos"][0] == 10

    assert mw.view_3d_manager.current_mol is not None
    assert mw.view_3d_manager.current_mol.GetNumAtoms() == mol.GetNumAtoms()
    assert (
        mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id")
        == aid1
    )
    assert (
        mw.view_3d_manager.current_mol.GetAtomWithIdx(1).GetIntProp("_original_atom_id")
        == aid2
    )

    assert mw.view_3d_manager.atom_positions_3d is not None
    assert mw.view_3d_manager.atom_positions_3d.shape == (mol.GetNumAtoms(), 3)

    assert len(mw.edit_3d_manager.constraints_3d) == 1
    assert mw.edit_3d_manager.constraints_3d[0][0] == "DISTANCE"
    assert mw.edit_3d_manager.constraints_3d[0][2] == 1.43
    assert mw.compute_manager.last_successful_optimization_method == "MMFF94s"
    assert mw._preserved_plugin_data["TestPlugin"]["val"] == 42


def test_undo_state_binary_roundtrip(dummy_window):
    """Test the internal binary state serialization used for Undo/Redo."""
    mw = dummy_window

    aid = mw.state_manager.data.add_atom("N", QPointF(5, 5))

    mol = Chem.MolFromSmiles("N")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", aid)
    mw.view_3d_manager.current_mol = mol

    state = mw.get_current_state()

    assert "mol_3d" in state
    assert isinstance(state["mol_3d"], bytes)
    assert state["mol_3d_atom_ids"][0] == aid

    mw.clear_2d_editor(push_to_undo=False)
    mw.set_state_from_data(state)

    assert mw.state_manager.data.atoms[aid]["symbol"] == "N"
    assert (
        mw.view_3d_manager.current_mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id")
        == aid
    )


def test_legacy_version_handling(dummy_window):
    """Verify that version mismatch warnings are triggered (but don't crash)."""
    mw = dummy_window

    json_data = {
        "format": "PME Project",
        "version": "99.0.0",
        "2d_structure": {"atoms": [], "bonds": []},
    }

    with MagicMock() as mock_msg:
        mw.warning_message_box = mock_msg
        mw.load_from_json_data(json_data)
