"""Unit tests for implicit hydrogen management."""

import pytest

from moleditpy.core.molecular_data import MolecularData
from moleditpy.ui.edit_actions_logic import EditActionsManager
from PyQt6.QtCore import QPointF
from unittest.mock import patch

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")


class DummyEditActions(EditActionsManager):
    def __init__(self, host):
        super().__init__(host)
        self._host = host
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene
        self.view_2d = host.init_manager.view_2d
        self.plotter = host.view_3d_manager.plotter
        self.settings = host.init_manager.settings
        self.current_mol = host.view_3d_manager.current_mol
        self.host.init_manager.current_file_path = host.init_manager.current_file_path
        self.host.state_manager.has_unsaved_changes = (
            host.state_manager.has_unsaved_changes
        )
        self.main_window_edit_actions = self

    def statusBar(self):
        return self._host.statusBar()

    def push_undo_state(self):
        self._host.push_undo_state()


def test_add_hydrogen_atoms_app_logic(mock_parser_host):
    """Verify add_hydrogen_atoms creates items in the scene based on RDKit logic."""
    actions = DummyEditActions(mock_parser_host)
    # Add a Carbon atom via scene to trigger item mock setup
    actions.scene.create_atom("C", QPointF(100, 100))

    # Run the app logic
    actions.add_hydrogen_atoms()

    # scene.create_atom should have been called 4 times for 'H'
    h_calls = [
        c
        for c in actions.scene.create_atom.call_args_list
        if c.args[0] == "H" or c.kwargs.get("symbol") == "H"
    ]
    assert len(h_calls) == 4
    # scene.create_bond should have 4 calls
    assert actions.scene.create_bond.call_count == 4


def test_remove_hydrogen_atoms_app_logic(mock_parser_host):
    """Verify remove_hydrogen_atoms finds and deletes H items using app logic."""
    actions = DummyEditActions(mock_parser_host)

    # Add C and H via scene
    c_id = actions.scene.create_atom("C", QPointF(100, 100))
    h_id = actions.scene.create_atom("H", QPointF(150, 100))

    # Create bond using items (mocked items should exist now)
    actions.scene.create_bond(
        actions.scene.atom_items[c_id], actions.scene.atom_items[h_id]
    )

    # We need to ensure sip_isdeleted_safe doesn't block deletion in test
    with patch(
        "moleditpy.ui.edit_actions_logic.sip_isdeleted_safe", return_value=False
    ):
        actions.remove_hydrogen_atoms()

    # Should call scene.delete_items with a set containing the H item
    actions.scene.delete_items.assert_called()
    deleted_set = actions.scene.delete_items.call_args[0][0]
    # In mock_parser_host, host.state_manager.data.atoms[h_id]['item'] is a MagicMock
    h_item = actions.scene.atom_items[h_id]
    assert h_item in deleted_set
    assert actions.scene.atom_items[c_id] not in deleted_set


def _build_2d_data(smiles):
    """Load *smiles* into a MolecularData the way string_importers does."""
    ref = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(ref)
    Chem.Kekulize(ref)
    data = MolecularData()
    idx_to_id = {}
    conf = ref.GetConformer()
    for i in range(ref.GetNumAtoms()):
        atom = ref.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        idx_to_id[i] = data.add_atom(
            atom.GetSymbol(),
            (pos.x * 50.0, -pos.y * 50.0),
            charge=atom.GetFormalCharge(),
        )
    for bond in ref.GetBonds():
        data.add_bond(
            idx_to_id[bond.GetBeginAtomIdx()],
            idx_to_id[bond.GetEndAtomIdx()],
            order=int(bond.GetBondTypeAsDouble()),
        )
    return ref, data


# Sanitization assigns an aromatic N-H to NumExplicitHs, where reading only
# GetNumImplicitHs() reported 0 and dropped the H from the 2D label.
@pytest.mark.parametrize(
    "smiles",
    [
        "NC1=CC=CC2=C1C(=O)NNC2=O",  # luminol
        "c1cc[nH]c1",  # pyrrole
        "c1cnc[nH]1",  # imidazole
        "c1ccc2[nH]ccc2c1",  # indole
        "O=c1cc[nH]c(=O)[nH]1",  # uracil
        "Nc1ncnc2[nH]cnc12",  # adenine
        "Cc1ccccc1",  # toluene, no aromatic N-H
        "CC(=O)O",  # acetic acid
        "[O-]C(=O)c1ccccc1",  # charged
    ],
)
def test_h_label_counts_match_molecule(mock_parser_host, smiles):
    """The 2D H label must equal the hydrogens the molecule actually carries."""
    ref, data = _build_2d_data(smiles)
    mol = data.to_rdkit_mol()
    assert mol is not None

    actions = DummyEditActions(mock_parser_host)
    counts = actions._compute_h_counts(mol)

    expected = {
        atom.GetIntProp("_original_atom_id"): ref.GetAtomWithIdx(
            atom.GetIntProp("_original_atom_id")
        ).GetTotalNumHs()
        for atom in mol.GetAtoms()
    }
    assert counts == expected


def test_h_label_ignores_hydrogens_drawn_as_atoms(mock_parser_host):
    """H drawn as its own atom must not also be counted in the label."""
    data = MolecularData()
    c_id = data.add_atom("C", (0.0, 0.0))
    for i in range(4):
        h_id = data.add_atom("H", (50.0 * (i + 1), 0.0))
        data.add_bond(c_id, h_id, order=1)

    mol = data.to_rdkit_mol()
    assert mol is not None
    actions = DummyEditActions(mock_parser_host)
    assert actions._compute_h_counts(mol)[c_id] == 0


def test_add_hydrogen_atoms_covers_aromatic_nh(mock_parser_host):
    """Add Hydrogens must not skip an aromatic N-H (luminol is C8H7N3O2)."""
    actions = DummyEditActions(mock_parser_host)
    ref, _ = _build_2d_data("NC1=CC=CC2=C1C(=O)NNC2=O")

    conf = ref.GetConformer()
    idx_to_id = {}
    for i in range(ref.GetNumAtoms()):
        atom = ref.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        idx_to_id[i] = actions.scene.create_atom(
            atom.GetSymbol(), QPointF(pos.x * 50.0, -pos.y * 50.0)
        )
    for bond in ref.GetBonds():
        actions.scene.create_bond(
            actions.scene.atom_items[idx_to_id[bond.GetBeginAtomIdx()]],
            actions.scene.atom_items[idx_to_id[bond.GetEndAtomIdx()]],
            bond_order=int(bond.GetBondTypeAsDouble()),
        )

    before = len(
        [c for c in actions.scene.create_atom.call_args_list if c.args[0] == "H"]
    )
    with patch(
        "moleditpy.ui.edit_actions_logic.sip_isdeleted_safe", return_value=False
    ):
        actions.add_hydrogen_atoms()
    after = [c for c in actions.scene.create_atom.call_args_list if c.args[0] == "H"]

    nitrogens = [idx_to_id[a.GetIdx()] for a in ref.GetAtoms() if a.GetSymbol() == "N"]
    added_to = [
        bond_call.args[0].atom_id
        for bond_call in actions.scene.create_bond.call_args_list
        if bond_call.args[0].atom_id in nitrogens
    ]
    assert len(after) - before == 7
    # The two ring N-H plus the two amine H: every N must receive hydrogen.
    assert set(nitrogens) == set(added_to)
