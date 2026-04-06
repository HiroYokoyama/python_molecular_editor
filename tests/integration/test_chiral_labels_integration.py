from rdkit import Chem
from moleditpy.ui.mirror_dialog import MirrorDialog


def test_chiral_labels_toggle_3d(window, qtbot):
    """
    Test that toggling 'Show Chiral Labels' correctly displays/hides labels in 3D.
    """
    window.string_importer_manager.load_from_smiles("C[C@H](O)CC")
    window.compute_manager.trigger_conversion()
    qtbot.waitUntil(lambda: window.view_3d_manager.current_mol is not None, timeout=15000)

    # Reset mock to clear any initialization calls
    window.view_3d_manager.plotter.add_point_labels.reset_mock()

    assert window.view_3d_manager.show_chiral_labels is False

    def labels_drawn():
        return any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.view_3d_manager.plotter.add_point_labels.call_args_list
        )

    assert not labels_drawn()

    window.init_manager.toggle_chiral_action.setChecked(True)
    window.init_manager.toggle_chiral_action.triggered.emit(True)

    assert window.view_3d_manager.show_chiral_labels is True

    mol = window.view_3d_manager.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assert len(chiral_centers) == 1
    initial_label = chiral_centers[0][1]
    assert initial_label in ["R", "S"]

    qtbot.waitUntil(labels_drawn, timeout=5000)

    chiral_call = next(
        call
        for call in window.view_3d_manager.plotter.add_point_labels.call_args_list
        if call.kwargs.get("name") == "chiral_labels"
    )
    assert len(chiral_call.args) > 1, "add_point_labels was not called with a labels argument"
    labels = chiral_call.args[1]
    assert initial_label in labels

    # Commented out to avoid fatal exception in headless tests
    # atom_idx, label = chiral_centers[0]
    # atom_id = mol.GetAtomWithIdx(atom_idx).GetIntProp("_original_atom_id")
    # assert window.state_manager.data.atoms[atom_id]['item'].chiral_label == initial_label


def test_chiral_labels_mirror_inversion_3d(window, qtbot):
    """
    Test that mirror transformation inverts the chiral label in 3D.
    """
    window.string_importer_manager.load_from_smiles("C[C@H](O)CC")
    window.compute_manager.trigger_conversion()
    qtbot.waitUntil(lambda: window.view_3d_manager.current_mol is not None, timeout=15000)

    window.init_manager.toggle_chiral_action.setChecked(True)
    window.init_manager.toggle_chiral_action.triggered.emit(True)

    qtbot.waitUntil(
        lambda: any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.view_3d_manager.plotter.add_point_labels.call_args_list
        ),
        timeout=5000,
    )
    mol = window.view_3d_manager.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    initial_label = chiral_centers[0][1]
    assert initial_label in ["R", "S"]

    dialog = MirrorDialog(window.view_3d_manager.current_mol, window)
    dialog.xy_radio.setChecked(True)

    window.view_3d_manager.plotter.add_point_labels.reset_mock()
    dialog.apply_mirror()
    dialog.close()  # Close dialog

    mol = window.view_3d_manager.current_mol
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    assert len(chiral_centers) == 1
    new_label = chiral_centers[0][1]
    assert new_label != initial_label
    assert new_label in ["R", "S"]

    qtbot.waitUntil(
        lambda: any(
            call.kwargs.get("name") == "chiral_labels"
            for call in window.view_3d_manager.plotter.add_point_labels.call_args_list
        ),
        timeout=5000,
    )

    chiral_call = next(
        call
        for call in window.view_3d_manager.plotter.add_point_labels.call_args_list
        if call.kwargs.get("name") == "chiral_labels"
    )
    assert len(chiral_call.args) > 1, "add_point_labels was not called with a labels argument"
    labels = chiral_call.args[1]
    assert new_label in labels
    assert initial_label not in labels

    # Commented out to avoid fatal exception in headless tests
    # atom_idx, label = chiral_centers[0]
    # atom_id = mol.GetAtomWithIdx(atom_idx).GetIntProp("_original_atom_id")
    # assert window.state_manager.data.atoms[atom_id]['item'].chiral_label == new_label
