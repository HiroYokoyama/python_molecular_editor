import pytest
import sys
import numpy as np
from rdkit import Chem
from PyQt6.QtCore import Qt
from moleditpy.modules.mirror_dialog import MirrorDialog


@pytest.mark.skipif(sys.platform == "win32", reason="RPC fatal crash on Windows CI")
def test_chiral_labels_toggle_3d(window, qtbot):
    """
    Test that toggling 'Show Chiral Labels' correctly displays/hides labels in 3D.
    """
    window.load_from_smiles("C[C@H](O)CC")
    window.trigger_conversion()
    qtbot.waitUntil(lambda: window.current_mol is not None, timeout=15000)

    # Reset mock to clear any initialization calls
    window.plotter.add_point_labels.reset_mock()

    assert window.show_chiral_labels is False

    def labels_drawn():
        return any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.plotter.add_point_labels.call_args_list
        )

    assert not labels_drawn()

    window.toggle_chiral_action.setChecked(True)
    window.toggle_chiral_action.triggered.emit(True)

    assert window.show_chiral_labels is True

    mol = window.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assert len(chiral_centers) == 1
    initial_label = chiral_centers[0][1]
    assert initial_label in ["R", "S"]

    qtbot.waitUntil(labels_drawn, timeout=5000)

    chiral_call = next(
        call
        for call in window.plotter.add_point_labels.call_args_list
        if call.kwargs.get("name") == "chiral_labels"
    )
    labels = chiral_call.args[1]
    assert initial_label in labels

    # Commented out to avoid fatal exception in headless tests
    # atom_idx, label = chiral_centers[0]
    # atom_id = mol.GetAtomWithIdx(atom_idx).GetIntProp("_original_atom_id")
    # assert window.data.atoms[atom_id]['item'].chiral_label == initial_label


@pytest.mark.skipif(sys.platform == "win32", reason="RPC fatal crash on Windows CI")
def test_chiral_labels_mirror_inversion_3d(window, qtbot):
    """
    Test that mirror transformation inverts the chiral label in 3D.
    """
    window.load_from_smiles("C[C@H](O)CC")
    window.trigger_conversion()
    qtbot.waitUntil(lambda: window.current_mol is not None, timeout=15000)

    window.toggle_chiral_action.setChecked(True)
    window.toggle_chiral_action.triggered.emit(True)

    qtbot.waitUntil(
        lambda: any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.plotter.add_point_labels.call_args_list
        ),
        timeout=5000,
    )
    mol = window.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    initial_label = chiral_centers[0][1]
    assert initial_label in ["R", "S"]

    dialog = MirrorDialog(window.current_mol, window)
    dialog.xy_radio.setChecked(True)

    window.plotter.add_point_labels.reset_mock()
    dialog.apply_mirror()
    dialog.close()  # Close dialog

    mol = window.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assert len(chiral_centers) == 1
    new_label = chiral_centers[0][1]
    assert new_label != initial_label
    assert new_label in ["R", "S"]

    qtbot.waitUntil(
        lambda: any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.plotter.add_point_labels.call_args_list
        ),
        timeout=5000,
    )

    chiral_call = next(
        call
        for call in window.plotter.add_point_labels.call_args_list
        if call.kwargs.get("name") == "chiral_labels"
    )
    labels = chiral_call.args[1]
    assert new_label in labels
    assert initial_label not in labels

    # Commented out to avoid fatal exception in headless tests
    # atom_idx, label = chiral_centers[0]
    # atom_id = mol.GetAtomWithIdx(atom_idx).GetIntProp("_original_atom_id")
    # assert window.data.atoms[atom_id]['item'].chiral_label == new_label
