"""
Shared contract tests for any BasePickingDialog subclass.

Consuming test class must provide a make_dialog fixture that returns (dlg, mol, mw).
The factory must accept keyword args: mol=None, preselected_atoms=None.
"""

import numpy as np
from unittest.mock import MagicMock, patch
from rdkit import Chem
from rdkit.Chem import AllChem


def make_ethane():
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def make_mock_mw(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    return mw


class PickingDialogContractTests:
    def test_pick_adds_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.on_atom_picked(0)
        assert 0 in dlg.selected_atoms

    def test_repick_removes_atom(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.on_atom_picked(0)
        dlg.on_atom_picked(0)
        assert 0 not in dlg.selected_atoms

    def test_multiple_picks_accumulate(self, make_dialog):
        dlg, _, _ = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        assert dlg.selected_atoms == {0, 1, 2, 3}

    def test_preselected_atoms_loaded(self, make_dialog):
        dlg, _, _ = make_dialog(preselected_atoms=[0, 1, 2])
        assert dlg.selected_atoms == {0, 1, 2}

    def test_clear_empties_selection(self, make_dialog):
        dlg, _, _ = make_dialog()
        for i in range(4):
            dlg.on_atom_picked(i)
        with patch.object(type(dlg), "clear_atom_labels"):
            dlg.clear_selection()
        assert len(dlg.selected_atoms) == 0
        assert not dlg.apply_button.isEnabled()

    def test_select_all_selects_every_atom(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.select_all_atoms()
        assert dlg.selected_atoms == set(range(mol.GetNumAtoms()))

    def test_select_all_enables_apply(self, make_dialog):
        dlg, mol, _ = make_dialog()
        with patch.object(type(dlg), "show_atom_labels"):
            dlg.select_all_atoms()
        assert dlg.apply_button.isEnabled()

    def test_zero_atoms_disables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms.clear()
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_two_atoms_disables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1}
        dlg.update_display()
        assert not dlg.apply_button.isEnabled()

    def test_three_atoms_enables_apply(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1, 2}
        dlg.update_display()
        assert dlg.apply_button.isEnabled()

    def test_count_shown_in_label(self, make_dialog):
        dlg, _, _ = make_dialog()
        dlg.selected_atoms = {0, 1, 2, 3}
        dlg.update_display()
        assert "4" in dlg.selection_label.text()
