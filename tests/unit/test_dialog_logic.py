import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from moleditpy.ui.bond_length_dialog import BondLengthDialog
from moleditpy.ui.alignment_dialog import AlignmentDialog
from moleditpy.ui.angle_dialog import AngleDialog
from moleditpy.ui.dihedral_dialog import DihedralDialog
from moleditpy.ui.translation_dialog import TranslationDialog
from moleditpy.ui.move_group_dialog import MoveGroupDialog
from moleditpy.core.mol_geometry import (
    calc_distance,
    calculate_dihedral,
    calc_angle_deg,
)


@pytest.fixture
def mol():
    m = Chem.MolFromSmiles("CCCCCO")  # Hexanol-like
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    return m


def test_bond_length_adjustment_logic(mock_parser_host, mol):
    """Test the geometric logic of bond length adjustment directly."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.bond_length_dialog.BondLengthDialog.init_ui"):
        dialog = BondLengthDialog(mol, window)
        dialog.atom1_idx = 0
        dialog.atom2_idx = 1

        dialog.atom1_fix_group_radio = MagicMock()
        dialog.atom1_fix_group_radio.isChecked.return_value = True
        dialog.atom1_fix_radio = MagicMock()
        dialog.atom1_fix_radio.isChecked.return_value = False
        dialog.both_groups_radio = MagicMock()
        dialog.both_groups_radio.isChecked.return_value = False

        dialog.apply_geometry_update(2.0)

        final_dist = calc_distance(
            mol.GetConformer().GetAtomPosition(0), mol.GetConformer().GetAtomPosition(1)
        )
        assert pytest.approx(final_dist, abs=1e-3) == 2.0


def test_alignment_logic(mock_parser_host, mol):
    """Test the geometry logic for aligning a bond to a specific axis."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.alignment_dialog.AlignmentDialog.init_ui"):
        dialog = AlignmentDialog(mol, window, axis="x")
        dialog.selected_atoms = {0, 1}

        with patch("moleditpy.ui.alignment_dialog.QMessageBox"):
            dialog.apply_alignment()

        new_conf = mol.GetConformer()
        p0 = np.array(new_conf.GetAtomPosition(0))
        p1 = np.array(new_conf.GetAtomPosition(1))

        assert np.allclose(p0, [0, 0, 0], atol=1e-7)
        assert p1[0] > 0
        assert pytest.approx(p1[1], abs=1e-7) == 0
        assert pytest.approx(p1[2], abs=1e-7) == 0


def test_angle_adjustment_logic(mock_parser_host):
    """Test the geometric logic of bond angle adjustment directly."""
    window = mock_parser_host
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.angle_dialog.AngleDialog.init_ui"):
        dialog = AngleDialog(mol, window)
        dialog.atom1_idx = 1
        dialog.atom2_idx = 0
        dialog.atom3_idx = 2
        dialog.both_groups_radio = MagicMock()
        dialog.both_groups_radio.isChecked.return_value = False
        dialog.rotate_atom_radio = MagicMock()
        dialog.rotate_atom_radio.isChecked.return_value = False

        initial_angle = calc_angle_deg(
            conf.GetAtomPosition(1), conf.GetAtomPosition(0), conf.GetAtomPosition(2)
        )
        assert initial_angle != pytest.approx(120.0)

        dialog.apply_geometry_update(120.0)

        final_angle = calc_angle_deg(
            mol.GetConformer().GetAtomPosition(1),
            mol.GetConformer().GetAtomPosition(0),
            mol.GetConformer().GetAtomPosition(2),
        )
        assert pytest.approx(final_angle, abs=1e-2) == 120.0


def test_dihedral_adjustment_logic(mock_parser_host):
    """Test the geometric logic of dihedral angle adjustment."""
    window = mock_parser_host
    mol = Chem.MolFromSmiles("OO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.dihedral_dialog.DihedralDialog.init_ui"):
        dialog = DihedralDialog(mol, window)
        dialog.atom1_idx = 2
        dialog.atom2_idx = 0
        dialog.atom3_idx = 1
        dialog.atom4_idx = 3
        dialog.both_groups_radio = MagicMock()
        dialog.both_groups_radio.isChecked.return_value = False
        dialog.rotate_atom_radio = MagicMock()
        dialog.rotate_atom_radio.isChecked.return_value = False
        dialog.rotate_group_radio = MagicMock()
        dialog.rotate_group_radio.isChecked.return_value = True
        dialog.both_groups_radio = MagicMock()
        dialog.both_groups_radio.isChecked.return_value = False

        dialog.apply_geometry_update(180.0)

        final_dihedral = calculate_dihedral(
            mol.GetConformer().GetPositions(), 2, 0, 1, 3
        )
        assert pytest.approx(abs(final_dihedral), abs=1e-2) == 180.0


def test_translation_logic(mock_parser_host, mol):
    """Test the geometric logic of centroid-based translation."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )
    window.view_3d_manager.current_mol = mol

    with patch("moleditpy.ui.translation_dialog.TranslationDialog.init_ui"):
        dialog = TranslationDialog(mol, window)
        dialog.selected_atoms = {0, 1}
        dialog.dx_input = MagicMock()
        dialog.dy_input = MagicMock()
        dialog.dz_input = MagicMock()
        dialog.dx_input.text.return_value = "10.0"
        dialog.dy_input.text.return_value = "10.0"
        dialog.dz_input.text.return_value = "10.0"
        p0_initial = np.array(conf.GetAtomPosition(0))
        p5_initial = np.array(conf.GetAtomPosition(5))
        translation_vec = np.array([10.0, 10.0, 10.0])
        expected_pos0 = p0_initial + translation_vec

        with patch("moleditpy.ui.translation_dialog.QMessageBox"):
            dialog.apply_translation()

        new_pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
        new_pos5 = np.array(mol.GetConformer().GetAtomPosition(5))

        # Atom 0 should have moved by delta
        assert np.allclose(new_pos0, expected_pos0, atol=1e-7)
        # Atom 5 (not selected) should NOT have moved
        assert np.allclose(new_pos5, p5_initial, atol=1e-7)


def test_move_group_logic(mock_parser_host, mol):
    """Test the translation and rotation logic in MoveGroupDialog."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with (
        patch("moleditpy.ui.move_group_dialog.MoveGroupDialog.init_ui"),
        patch("moleditpy.ui.move_group_dialog.MoveGroupDialog.show_atom_labels"),
    ):
        dialog = MoveGroupDialog(mol, window)
        dialog.group_atoms = {0, 1, 2}

        dialog.x_trans_input = MagicMock()
        dialog.y_trans_input = MagicMock()
        dialog.z_trans_input = MagicMock()
        dialog.x_trans_input.text.return_value = "5.0"
        dialog.y_trans_input.text.return_value = "0.0"
        dialog.z_trans_input.text.return_value = "0.0"

        initial_pos0 = np.array(conf.GetAtomPosition(0))
        initial_pos5 = np.array(conf.GetAtomPosition(5))

        with patch("moleditpy.ui.move_group_dialog.QMessageBox"):
            dialog.apply_translation()

        assert np.allclose(
            np.array(mol.GetConformer().GetAtomPosition(0)), initial_pos0 + [5, 0, 0]
        )
        assert np.allclose(
            np.array(mol.GetConformer().GetAtomPosition(5)), initial_pos5
        )

        dialog.x_rot_input = MagicMock()
        dialog.y_rot_input = MagicMock()
        dialog.z_rot_input = MagicMock()
        dialog.x_rot_input.text.return_value = "0.0"
        dialog.y_rot_input.text.return_value = "0.0"
        dialog.z_rot_input.text.return_value = "90.0"

        pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
        pos1 = np.array(mol.GetConformer().GetAtomPosition(1))
        pos2 = np.array(mol.GetConformer().GetAtomPosition(2))
        centroid = (pos0 + pos1 + pos2) / 3.0

        with patch("moleditpy.ui.move_group_dialog.QMessageBox"):
            dialog.apply_rotation()

        new_pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
        relative_pos = pos0 - centroid
        expected_rotated = (
            np.array([-relative_pos[1], relative_pos[0], relative_pos[2]]) + centroid
        )
        assert np.allclose(new_pos0, expected_rotated, atol=1e-7)


# ---------------------------------------------------------------------------
# TranslationDialog — Absolute tab
# ---------------------------------------------------------------------------

def _make_abs_dialog(mol, window):
    """Create a TranslationDialog with init_ui patched, tabs mocked for Absolute tab."""
    window.view_3d_manager.current_mol = mol
    with patch("moleditpy.ui.translation_dialog.TranslationDialog.init_ui"):
        dialog = TranslationDialog(mol, window)
    # Simulate Absolute tab active
    dialog.tabs = MagicMock()
    dialog.tabs.currentIndex.return_value = 0  # _TAB_ABSOLUTE
    dialog.abs_selection_label = MagicMock()
    dialog.abs_apply_btn = MagicMock()
    dialog.abs_x_input = MagicMock()
    dialog.abs_y_input = MagicMock()
    dialog.abs_z_input = MagicMock()
    dialog.move_mol_checkbox = MagicMock()
    return dialog


def test_abs_single_atom_enforcement(mock_parser_host, mol):
    """Clicking a second atom in Absolute tab replaces the first selection."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)

    dialog.on_atom_picked(0)
    assert dialog.selected_atoms == {0}

    dialog.on_atom_picked(3)
    assert dialog.selected_atoms == {3}, "Second pick should replace first"


def test_abs_coordinate_inputs_populated_on_pick(mock_parser_host, mol):
    """Selecting an atom auto-fills X/Y/Z inputs with its current position."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)
    dialog.on_atom_picked(1)

    pos = conf.GetAtomPosition(1)
    dialog.abs_x_input.setText.assert_called_once_with(f"{pos.x:.4f}")
    dialog.abs_y_input.setText.assert_called_once_with(f"{pos.y:.4f}")
    dialog.abs_z_input.setText.assert_called_once_with(f"{pos.z:.4f}")


def test_abs_apply_moves_entire_molecule_by_default(mock_parser_host, mol):
    """apply_absolute with move_mol=True shifts every atom by the same delta."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)
    dialog.selected_atoms = {0}
    dialog.move_mol_checkbox.isChecked.return_value = True

    atom0_pos = np.array(conf.GetAtomPosition(0))
    atom5_pos = np.array(conf.GetAtomPosition(5))
    target = atom0_pos + np.array([5.0, 0.0, 0.0])

    dialog.abs_x_input.text.return_value = f"{target[0]:.4f}"
    dialog.abs_y_input.text.return_value = f"{target[1]:.4f}"
    dialog.abs_z_input.text.return_value = f"{target[2]:.4f}"

    with patch("moleditpy.ui.translation_dialog.QMessageBox"):
        dialog.apply_absolute()

    new_pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
    new_pos5 = np.array(mol.GetConformer().GetAtomPosition(5))

    assert np.allclose(new_pos0, atom0_pos + [5, 0, 0], atol=1e-4)
    assert np.allclose(new_pos5, atom5_pos + [5, 0, 0], atol=1e-4)


def test_abs_apply_moves_only_selected_atom(mock_parser_host, mol):
    """apply_absolute with move_mol=False moves only the selected atom."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)
    dialog.selected_atoms = {0}
    dialog.move_mol_checkbox.isChecked.return_value = False

    atom0_pos = np.array(conf.GetAtomPosition(0))
    atom5_pos = np.array(conf.GetAtomPosition(5))
    target = atom0_pos + np.array([0.0, 3.0, 0.0])

    dialog.abs_x_input.text.return_value = f"{target[0]:.4f}"
    dialog.abs_y_input.text.return_value = f"{target[1]:.4f}"
    dialog.abs_z_input.text.return_value = f"{target[2]:.4f}"

    with patch("moleditpy.ui.translation_dialog.QMessageBox"):
        dialog.apply_absolute()

    new_pos0 = np.array(mol.GetConformer().GetAtomPosition(0))
    new_pos5 = np.array(mol.GetConformer().GetAtomPosition(5))

    assert np.allclose(new_pos0, atom0_pos + [0, 3, 0], atol=1e-4)
    assert np.allclose(new_pos5, atom5_pos, atol=1e-4), "Unselected atom must not move"


def test_abs_apply_noop_when_no_delta(mock_parser_host, mol):
    """apply_absolute does nothing when target equals current position."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)
    dialog.selected_atoms = {0}
    dialog.move_mol_checkbox.isChecked.return_value = True

    pos = conf.GetAtomPosition(0)
    # Use full precision so the parsed delta is within np.allclose(delta, 0) tolerance
    dialog.abs_x_input.text.return_value = repr(pos.x)
    dialog.abs_y_input.text.return_value = repr(pos.y)
    dialog.abs_z_input.text.return_value = repr(pos.z)

    original_positions = mol.GetConformer().GetPositions().copy()

    with patch("moleditpy.ui.translation_dialog.QMessageBox"):
        dialog.apply_absolute()

    assert np.allclose(mol.GetConformer().GetPositions(), original_positions, atol=1e-7)


def test_abs_apply_requires_one_atom(mock_parser_host, mol):
    """apply_absolute shows a warning and does nothing when no atom is selected."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    dialog = _make_abs_dialog(mol, window)
    dialog.selected_atoms = set()
    dialog.abs_x_input.text.return_value = "1.0"
    dialog.abs_y_input.text.return_value = "1.0"
    dialog.abs_z_input.text.return_value = "1.0"

    original_positions = mol.GetConformer().GetPositions().copy()

    with patch("moleditpy.ui.translation_dialog.QMessageBox") as mock_mb:
        dialog.apply_absolute()
        mock_mb.warning.assert_called_once()

    assert np.allclose(mol.GetConformer().GetPositions(), original_positions, atol=1e-7)


def test_tab_switch_clears_selection(mock_parser_host, mol):
    """Switching tabs clears selected_atoms and atom labels."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.translation_dialog.TranslationDialog.init_ui"):
        dialog = TranslationDialog(mol, window)

    dialog.tabs = MagicMock()
    dialog.tabs.currentIndex.return_value = 1  # Delta tab after switch
    dialog.abs_selection_label = MagicMock()
    dialog.abs_apply_btn = MagicMock()
    dialog.delta_selection_label = MagicMock()
    dialog.apply_button = MagicMock()
    dialog.selected_atoms = {0, 1, 2}

    with patch.object(dialog, "clear_atom_labels") as mock_clear:
        dialog._on_tab_changed(1)
        assert dialog.selected_atoms == set()
        mock_clear.assert_called_once()


def test_delta_tab_atom_pick_toggles(mock_parser_host, mol):
    """In Delta tab, picking an already-selected atom deselects it."""
    window = mock_parser_host
    conf = mol.GetConformer()
    window.view_3d_manager.atom_positions_3d = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    )

    with patch("moleditpy.ui.translation_dialog.TranslationDialog.init_ui"):
        dialog = TranslationDialog(mol, window)

    dialog.tabs = MagicMock()
    dialog.tabs.currentIndex.return_value = 1  # _TAB_DELTA
    dialog.delta_selection_label = MagicMock()
    dialog.apply_button = MagicMock()
    dialog.selected_atoms = {0}

    with patch.object(dialog, "show_atom_labels"), patch.object(dialog, "update_display"):
        dialog.on_atom_picked(0)

    assert 0 not in dialog.selected_atoms, "Re-picking selected atom should deselect it"
