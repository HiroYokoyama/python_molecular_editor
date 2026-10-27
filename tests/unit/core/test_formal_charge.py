"""Formal-charge inference and editor action regression tests."""

import copy
from unittest.mock import MagicMock, patch

import pytest
from PySide6.QtWidgets import QMainWindow
from rdkit import Chem

from moleditpy.core.formal_charge import infer_charge_adjustments
from moleditpy.core.molecular_data import MolecularData
from moleditpy.ui.edit_actions_logic import EditActionsManager
from moleditpy.ui.main_window_init import MainInitManager


def star(symbol, bonds, charge=0, radical=0):
    data = MolecularData()
    center = data.add_atom(symbol, (0, 0), charge=charge, radical=radical)
    for i, order in enumerate(bonds):
        neighbor = data.add_atom("C", (i + 1, 1))
        data.add_bond(center, neighbor, order)
    return data, center


@pytest.mark.parametrize(
    "symbol,bonds,charge,expected",
    [
        ("N", [1, 1, 1, 1], 0, 1),
        ("N", [2, 1, 1], 0, 1),
        ("O", [1, 1, 1], 0, 1),
        ("B", [1, 1, 1, 1], 0, -1),
        ("C", [1, 1, 1, 1], 1, 0),
        ("N", [1, 1, 1, 1], 1, None),
        ("N", [1, 1, 1], 0, None),
        ("N", [1, 1, 1], 1, None),
        ("O", [1], 0, None),
        ("O", [1], -1, None),
        ("C", [1, 1, 1], 1, None),
        ("C", [1, 1, 1, 1, 1], 0, None),
        ("N", [1, 1, 1], 2, None),
        ("Fe", [1, 1, 1, 1], 0, None),
    ],
)
def test_charge_inference_preserves_graph(symbol, bonds, charge, expected):
    data, center = star(symbol, bonds, charge)
    before = copy.deepcopy((data.atoms, data.bonds))
    assert infer_charge_adjustments(data) == (
        {} if expected is None else {center: expected}
    )
    assert (data.atoms, data.bonds) == before
    if expected is not None:
        data.atoms[center]["charge"] = expected
        assert data.to_rdkit_mol() is not None
        assert infer_charge_adjustments(data) == {}


def test_aromatic_cases_without_unique_local_valence_are_preserved():
    mol = Chem.MolFromSmiles("C[n+]1ccccc1")
    data = MolecularData()
    ids = [data.add_atom(atom.GetSymbol(), (0, 0)) for atom in mol.GetAtoms()]
    for bond in mol.GetBonds():
        data.add_bond(
            ids[bond.GetBeginAtomIdx()],
            ids[bond.GetEndAtomIdx()],
            bond.GetBondTypeAsDouble(),
        )
    assert infer_charge_adjustments(data) == {}


def test_radicals_and_unsupported_neighbors_are_preserved():
    data, center = star("N", [1, 1, 1, 1], radical=1)
    assert infer_charge_adjustments(data) == {}
    data.atoms[center]["radical"] = 0
    data.atoms[next(i for i in data.atoms if i != center)]["symbol"] = "R"
    assert infer_charge_adjustments(data) == {}


def test_disconnected_invalid_fragment_does_not_block_repair():
    data, center = star("N", [1, 1, 1, 1])
    data.add_atom("R", (10, 10))
    assert data.to_rdkit_mol() is None
    assert infer_charge_adjustments(data) == {center: 1}


def test_empty_graph():
    assert infer_charge_adjustments(MolecularData()) == {}


def test_action_updates_model_and_item_once(mock_parser_host):
    host = mock_parser_host
    data, center = star("N", [1, 1, 1, 1])
    host.state_manager.data = data
    item = MagicMock()
    host.init_manager.scene.atom_items = {center: item}
    manager = EditActionsManager(host)
    with patch.object(manager, "push_undo_state") as push:
        manager.adjust_formal_charges()
        assert data.atoms[center]["charge"] == 1
        assert item.charge == 1
        item.prepareGeometryChange.assert_called_once()
        item.update_style.assert_called_once()
        push.assert_called_once()
        manager.adjust_formal_charges()
        push.assert_called_once()


def test_action_empty_graph_does_not_push_history(mock_parser_host):
    manager = EditActionsManager(mock_parser_host)
    with patch.object(manager, "push_undo_state") as push:
        manager.adjust_formal_charges()
        push.assert_not_called()
    assert "No molecule" in mock_parser_host.statusBar().showMessage.call_args.args[0]


def test_action_charge_is_recorded_for_undo_redo(mock_parser_host):
    host = mock_parser_host
    data, center = star("N", [1, 1, 1, 1])
    host.state_manager.data = data
    host.is_restoring_state = False
    host.state_manager.get_current_state.side_effect = lambda: {
        "atoms": copy.deepcopy(data.atoms),
        "bonds": copy.deepcopy(data.bonds),
        "_next_atom_id": data.next_atom_id,
    }

    def restore(snapshot):
        data.atoms = copy.deepcopy(snapshot["atoms"])
        data.bonds = copy.deepcopy(snapshot["bonds"])

    host.state_manager.set_state_from_data.side_effect = restore
    manager = EditActionsManager(host)
    with patch.object(manager, "update_implicit_hydrogens"):
        manager.reset_history()
        manager.adjust_formal_charges()
        assert len(manager.undo_stack) == 2
        assert data.atoms[center]["charge"] == 1
        manager.undo()
        assert data.atoms[center]["charge"] == 0
        manager.redo()
        assert data.atoms[center]["charge"] == 1


def test_edit_menu_placement_and_trigger(app):
    window = QMainWindow()
    window.edit_actions_manager = MagicMock()
    window.compute_manager = MagicMock()
    manager = object.__new__(MainInitManager)
    manager.host = window
    try:
        manager._init_edit_menu(window.menuBar())
        menu = window.menuBar().actions()[0].menu()
        actions = menu.actions()
        labels = [action.text() for action in actions]
        idx = labels.index("Add Hydrogens")
        assert labels[idx : idx + 3] == [
            "Add Hydrogens",
            "Remove Hydrogens",
            "Adjust Charge",
        ]
        actions[idx + 2].trigger()
        assert actions[idx + 2].shortcut().toString() == "Ctrl+E"
        window.edit_actions_manager.adjust_formal_charges.assert_called_once()
    finally:
        window.close()
        window.deleteLater()
