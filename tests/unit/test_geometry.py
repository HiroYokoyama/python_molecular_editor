"""Unit tests for mol_geometry coordinate math and 3D calculations."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.ui.calculation_worker import CalculationWorker
from moleditpy.core.molecular_data import MolecularData
from PyQt6.QtCore import QPointF
from unittest import mock as _mock


def test_3d_bond_lengths(qtbot):
    """Verify optimized 3D coordinates yield physical bond lengths."""
    # Create Ethane C-C bond length test
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(10, 0))
    data.add_bond(c1, c2, order=1)

    mol_block = data.to_mol_block()
    print(f"DEBUG: mol_block='{mol_block}'")

    worker = CalculationWorker()

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, {"conversion_mode": "rdkit"})

    result = blocker.args[0]
    if isinstance(result, tuple):
        mol = result[1]
    else:
        mol = result

    assert mol is not None
    conf = mol.GetConformer()

    # RDKit indices might change, but here it should be 0 and 1 (heavy atoms)
    # plus hydrogens. We find C-C bond.
    cc_bond = None
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtom().GetSymbol() == "C"
            and bond.GetEndAtom().GetSymbol() == "C"
        ):
            cc_bond = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            break

    assert cc_bond is not None
    dist = rdMolTransforms.GetBondLength(conf, cc_bond[0], cc_bond[1])

    # Ethane C-C length is approx 1.54 A. Tightening tolerance to roughly 0.03 A
    assert 1.51 < dist < 1.57

    # Add sp3 angle verification (e.g., H-C-H or H-C-C should be ~109.5)
    # We find a C atom and its neighbors
    c_idx = cc_bond[0]
    neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()]
    if len(neighbors) >= 2:
        angle = rdMolTransforms.GetAngleDeg(conf, neighbors[0], c_idx, neighbors[1])
        # Typical range for sp3 in RDKit ETKDG is roughly 107-112
        assert 107 < angle < 112


from moleditpy.ui.mirror_dialog import MirrorDialog
from unittest.mock import MagicMock
from PyQt6.QtWidgets import QWidget


def test_mirror_dialog_logic(qtbot):
    """Verify MirrorDialog correctly manipulates coordinates and UI (software logic)."""
    # (S)-2-butanol
    mol = Chem.MolFromSmiles("C[C@H](O)CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Mock main_window using QWidget to satisfy PyQt6
    main_window = QWidget()
    main_window.statusBar = MagicMock()
    main_window.statusBar.return_value = MagicMock()

    # Initialize all managers to avoid AttributeErrors
    main_window.view_3d_manager = MagicMock()
    main_window.view_3d_manager.draw_molecule_3d = MagicMock()
    main_window.view_3d_manager.update_chiral_labels = MagicMock()

    main_window.state_manager = MagicMock()
    main_window.edit_actions_manager = MagicMock()
    main_window.init_manager = MagicMock()
    main_window.ui_manager = MagicMock()
    main_window.compute_manager = MagicMock()

    # Initial state
    conf = mol.GetConformer()
    x_before = conf.GetAtomPosition(0).x

    # Instantiate dialog
    dialog = MirrorDialog(mol, main_window)
    qtbot.addWidget(dialog)

    # select YZ plane (id=2)
    dialog.yz_radio.setChecked(True)

    # Apply mirror
    dialog.apply_mirror()

    # Verify coordinates inverted on X axis
    x_after = conf.GetAtomPosition(0).x
    assert x_after == pytest.approx(-x_before, abs=1e-3)

    # Verify MoleditPy specific recovery calls
    assert main_window.view_3d_manager.draw_molecule_3d.called
    assert main_window.view_3d_manager.update_chiral_labels.called
    assert main_window.edit_actions_manager.push_undo_state.called


from moleditpy.ui.planarize_dialog import PlanarizeDialog
import numpy as np


def test_planarize_logic(qtbot):
    """Verify planarize functionality using the actual PlanarizeDialog logic."""
    # (R)-2-butanol
    mol = Chem.MolFromSmiles("C[C@H](O)CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()

    # Mock main_window
    main_window = QWidget()
    positions = mol.GetConformer().GetPositions()

    # Initialize all managers
    main_window.view_3d_manager = MagicMock()
    main_window.state_manager = MagicMock()
    main_window.edit_actions_manager = MagicMock()
    main_window.init_manager = MagicMock()
    main_window.ui_manager = MagicMock()
    main_window.compute_manager = MagicMock()

    main_window.view_3d_manager.atom_positions_3d = positions
    main_window.view_3d_manager.plotter = MagicMock()
    main_window.view_3d_manager.draw_molecule_3d = MagicMock()
    main_window.view_3d_manager.update_chiral_labels = MagicMock()

    # Instantiate dialog with all atoms selected
    dialog = PlanarizeDialog(
        mol, main_window, preselected_atoms=range(mol.GetNumAtoms())
    )
    qtbot.addWidget(dialog)

    # Apply planarize
    # We mock QMessageBox to prevent blocking
    with _mock.patch("PyQt6.QtWidgets.QMessageBox.information"):
        dialog.apply_planarize()

    new_positions = conf.GetPositions()
    centroid = np.mean(new_positions, axis=0)
    centered = new_positions - centroid

    u, s, vh = np.linalg.svd(centered)
    # The smallest singular value s[-1] should be near 0 if they are planar
    assert s[-1] < 1e-10

    # Also verify that main_window methods were called
    assert main_window.view_3d_manager.draw_molecule_3d.called
    assert main_window.view_3d_manager.update_chiral_labels.called
    assert main_window.edit_actions_manager.push_undo_state.called


# ------------------------------------------------------------------
# rodrigues_rotate & adjust_bond_angle tests
# ------------------------------------------------------------------
from moleditpy.core.mol_geometry import adjust_bond_angle, rodrigues_rotate


def test_rodrigues_rotate_90_deg():
    """Rotate [1,0,0] by 90° around [0,0,1] → expect [0,1,0]."""
    v = np.array([1.0, 0.0, 0.0])
    axis = np.array([0.0, 0.0, 1.0])
    result = rodrigues_rotate(v, axis, np.pi / 2)
    np.testing.assert_allclose(result, [0.0, 1.0, 0.0], atol=1e-12)


def test_rodrigues_rotate_identity():
    """Rotate by 0° → vector unchanged."""
    v = np.array([3.0, -1.0, 2.0])
    axis = np.array([0.0, 1.0, 0.0])
    result = rodrigues_rotate(v, axis, 0.0)
    np.testing.assert_allclose(result, v, atol=1e-12)


def test_adjust_bond_angle_simple():
    """Set a 90° angle to 120° and verify."""
    # A at (1,0,0), B at origin, C at (0,1,0) → 90°
    positions = np.array(
        [
            [1.0, 0.0, 0.0],  # A (idx 0)
            [0.0, 0.0, 0.0],  # B (idx 1)
            [0.0, 1.0, 0.0],  # C (idx 2)
        ]
    )
    target = 120.0
    adjust_bond_angle(positions, 0, 1, 2, target, {2})

    # Verify resulting angle
    ba = positions[0] - positions[1]
    bc = positions[2] - positions[1]
    cos_a = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle_deg = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
    assert angle_deg == pytest.approx(target, abs=1e-8)

    # B should not have moved
    np.testing.assert_allclose(positions[1], [0, 0, 0], atol=1e-12)

    # Bond length of C from B should be preserved (was 1.0)
    assert np.linalg.norm(positions[2] - positions[1]) == pytest.approx(1.0, abs=1e-12)


def test_adjust_bond_angle_with_group():
    """Move multiple atoms; verify relative geometry is preserved."""
    # A(0), B(1) at origin, C(2), D(3) attached to C
    positions = np.array(
        [
            [1.0, 0.0, 0.0],  # A
            [0.0, 0.0, 0.0],  # B (vertex)
            [0.0, 1.0, 0.0],  # C
            [0.0, 2.0, 0.0],  # D (attached to C)
        ]
    )
    cd_before = np.linalg.norm(positions[3] - positions[2])

    adjust_bond_angle(positions, 0, 1, 2, 60.0, {2, 3})

    # Angle should be 60°
    ba = positions[0] - positions[1]
    bc = positions[2] - positions[1]
    cos_a = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle_deg = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
    assert angle_deg == pytest.approx(60.0, abs=1e-8)

    # C-D distance should be preserved
    cd_after = np.linalg.norm(positions[3] - positions[2])
    assert cd_after == pytest.approx(cd_before, abs=1e-12)


def test_adjust_bond_angle_collinear():
    """Collinear atoms → fallback axis is used, rotation still succeeds."""
    positions = np.array(
        [
            [2.0, 0.0, 0.0],  # A
            [0.0, 0.0, 0.0],  # B
            [-1.0, 0.0, 0.0],  # C  (collinear with A-B, angle = 180°)
        ]
    )
    target = 90.0
    delta = adjust_bond_angle(positions, 0, 1, 2, target, {2})

    # Should have rotated (non-zero delta)
    assert delta != 0.0

    # Verify resulting angle is 90°
    ba = positions[0] - positions[1]
    bc = positions[2] - positions[1]
    cos_a = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle_deg = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
    assert angle_deg == pytest.approx(target, abs=1e-8)

    # Bond length of C from B should be preserved (was 1.0)
    assert np.linalg.norm(positions[2] - positions[1]) == pytest.approx(1.0, abs=1e-12)


def test_optimize_2d_coords():
    """Verify 2D coordinate optimization generating coordinates for a simple molecule."""
    from rdkit import Chem
    from moleditpy.core.mol_geometry import optimize_2d_coords

    mol = Chem.MolFromSmiles("C1=CC=CC=C1")
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetIntProp("_original_atom_id", i)

    new_pos = optimize_2d_coords(mol)
    assert len(new_pos) == 6
    for pos in new_pos.values():
        assert len(pos) == 2


def test_calculate_best_fit_plane_projection():
    """Verify orthogonal projection onto a best-fit plane."""
    import numpy as np
    from moleditpy.core.mol_geometry import calculate_best_fit_plane_projection

    points = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3], [1, 2, 3]])
    centroid = np.mean(points, axis=0)
    centered = points - centroid
    u, s, vh = np.linalg.svd(centered)
    normal = vh[-1]
    projected = calculate_best_fit_plane_projection(centered, normal, centroid)

    # Check planarity: SVD of centered projected points should have s[-1] approx 0
    p_centered = projected - np.mean(projected, axis=0)
    pu, ps, pvh = np.linalg.svd(p_centered)
    assert ps[-1] < 1e-10


def test_rotate_2d_points():
    """Verify 2D rotation of point maps."""
    import numpy as np
    from moleditpy.core.mol_geometry import rotate_2d_points

    points = {1: (1, 0), 2: (0, 1)}
    rotated = rotate_2d_points(points, 0, 0, 90)
    np.testing.assert_allclose(rotated[1], [0, 1], atol=1e-12)
    np.testing.assert_allclose(rotated[2], [-1, 0], atol=1e-12)


def test_resolve_2d_overlaps():
    """Verify 2D overlap resolution logic handles collisions correctly."""
    from moleditpy.core.mol_geometry import resolve_2d_overlaps

    atom_ids = {1, 2}
    positions = {1: (0, 0), 2: (0, 0.1)}  # Overlapping if threshold > 0.1
    adj = {1: [], 2: []}  # Not connected
    # If they are NOT bonded, they should move
    moves = resolve_2d_overlaps(
        atom_ids,
        positions,
        adj,
        overlap_threshold=0.5,
        move_distance=1.0,
        has_bond_check_func=lambda i, j: False,
    )
    assert len(moves) > 0
    # One fragment should be planned to move
    assert len(moves[0][0]) == 1


# =============================================================================
# identify_valence_problems — RDKit-backed detection (red-box flagging)
# =============================================================================


def _data_with_center(sym, nbonds, charge=0, radical=0, order=1):
    from moleditpy.core.molecular_data import MolecularData

    d = MolecularData()
    center = d.add_atom(sym, (0, 0), charge=charge, radical=radical)
    for i in range(nbonds):
        c = d.add_atom("C", (50 * (i + 1), 0))
        d.add_bond(center, c, order=order)
    return d, center


@pytest.mark.parametrize(
    "sym,nbonds,charge,radical",
    [
        ("N", 3, 0, 1),  # trimethylamine radical (reported bug)
        ("N", 5, 1, 0),  # charged hypervalent N+
        ("S", 7, 0, 0),  # beyond expanded octet
        ("B", 5, 0, 0),
        ("O", 2, 0, 1),  # oxygen radical with full valence
        ("C", 4, 0, 1),  # carbon radical with full valence
        ("C", 5, 0, 0),
        ("O", 3, 0, 0),
    ],
)
def test_identify_valence_problems_flags_invalid(sym, nbonds, charge, radical):
    """Hypervalent/radical-overloaded atoms are flagged for any element."""
    from moleditpy.core.mol_geometry import identify_valence_problems

    d, center = _data_with_center(sym, nbonds, charge=charge, radical=radical)
    assert center in identify_valence_problems(d.atoms, d.bonds)


@pytest.mark.parametrize(
    "sym,nbonds,charge,radical",
    [
        ("N", 3, 0, 0),  # trimethylamine
        ("N", 4, 1, 0),  # ammonium-like
        ("S", 6, 0, 0),  # valid expanded octet
        ("C", 3, 0, 1),  # methyl-like radical
        ("O", 2, 0, 0),
    ],
)
def test_identify_valence_problems_accepts_valid(sym, nbonds, charge, radical):
    """Chemically valid atoms (incl. charges/radicals) are not flagged."""
    from moleditpy.core.mol_geometry import identify_valence_problems

    d, center = _data_with_center(sym, nbonds, charge=charge, radical=radical)
    assert center not in identify_valence_problems(d.atoms, d.bonds)


# ---------------------------------------------------------------------------
# Degenerate-input guards (pure geometry)
# ---------------------------------------------------------------------------

import numpy as _np


def test_calc_angle_deg_zero_length_returns_zero():
    from moleditpy.core.mol_geometry import calc_angle_deg

    # pos1 coincides with the vertex -> zero-length arm -> 0.0
    assert calc_angle_deg([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]) == 0.0


def test_get_connected_group_respects_exclude_wall():
    from moleditpy.core.mol_geometry import get_connected_group

    # Linear chain C0-C1-C2-C3; excluding C1 isolates C0
    mol = Chem.MolFromSmiles("CCCC")
    assert get_connected_group(mol, 0, exclude=1) == {0}
    # Without the wall the whole chain is reachable
    assert get_connected_group(mol, 0) == {0, 1, 2, 3}


def test_adjust_bond_angle_zero_length_arm_returns_zero():
    from moleditpy.core.mol_geometry import adjust_bond_angle

    positions = _np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float
    )  # A coincides with vertex B
    result = adjust_bond_angle(positions, 0, 1, 2, 120.0, [2])
    assert result == 0.0


def test_calculate_dihedral_collinear_returns_zero():
    from moleditpy.core.mol_geometry import calculate_dihedral

    positions = _np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        dtype=float,
    )  # all four points collinear -> zero normals
    assert calculate_dihedral(positions, 0, 1, 2, 3) == 0.0


def test_adjust_dihedral_wraps_positive_delta():
    from moleditpy.core.mol_geometry import adjust_dihedral, calculate_dihedral

    # Non-planar butane-like arrangement so the dihedral is well-defined
    positions = _np.array(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
        dtype=float,
    )
    before = calculate_dihedral(positions, 0, 1, 2, 3)
    target = before + 350.0
    adjust_dihedral(positions, 0, 1, 2, 3, target, [3])
    after = calculate_dihedral(positions, 0, 1, 2, 3)
    # +350 wraps to -10 (short way round); after should match target mod 360
    assert abs(((after - target + 180) % 360) - 180) < 1.0


def test_adjust_dihedral_no_op_when_already_at_target():
    from moleditpy.core.mol_geometry import adjust_dihedral, calculate_dihedral

    positions = _np.array(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
        dtype=float,
    )
    current = calculate_dihedral(positions, 0, 1, 2, 3)
    assert adjust_dihedral(positions, 0, 1, 2, 3, current, [3]) == 0.0


def test_adjust_dihedral_degenerate_axis_returns_zero():
    from moleditpy.core.mol_geometry import adjust_dihedral

    positions = _np.array(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
        dtype=float,
    )  # i2 == i3 -> zero-length rotation axis
    assert adjust_dihedral(positions, 0, 1, 2, 3, 90.0, [3]) == 0.0


# ---------------------------------------------------------------------------
# is_problematic_valence
# ---------------------------------------------------------------------------


def test_is_problematic_valence_unknown_symbol_is_never_flagged():
    from moleditpy.core.mol_geometry import is_problematic_valence

    assert is_problematic_valence("Xx", 99) is False


def test_is_problematic_valence_charged_atom_skips_neutral_only_limit():
    from moleditpy.core.mol_geometry import is_problematic_valence

    # Nitrogen's limit is neutral-only; a charged nitrogen bypasses the check
    assert is_problematic_valence("N", 4, charge=1) is False


def test_is_problematic_valence_over_limit_neutral_flagged():
    from moleditpy.core.mol_geometry import is_problematic_valence

    assert is_problematic_valence("C", 5, charge=0) is True


# ---------------------------------------------------------------------------
# identify_valence_problems fallback heuristic (RDKit path unavailable)
# ---------------------------------------------------------------------------


def test_identify_valence_problems_heuristic_fallback():
    from moleditpy.core import mol_geometry

    # Carbon 0 bonded to five hydrogens: only C overloads (each H stays at 1)
    atoms = {0: {"symbol": "C", "charge": 0}}
    atoms.update({i: {"symbol": "H"} for i in range(1, 6)})
    bonds = {(0, i): {"order": 1} for i in range(1, 6)}
    with _mock.patch.object(
        mol_geometry, "_identify_problems_rdkit", return_value=None
    ):
        result = mol_geometry.identify_valence_problems(atoms, bonds)
    assert result == [0]


def test_identify_problems_rdkit_returns_none_for_unknown_symbol():
    from moleditpy.core.mol_geometry import _identify_problems_rdkit

    atoms = {1: {"symbol": "Zz"}}  # not a real element -> construction fails
    assert _identify_problems_rdkit(atoms, {}) is None


def test_identify_problems_rdkit_detects_overbonded_carbon():
    from moleditpy.core.mol_geometry import _identify_problems_rdkit

    atoms = {i: {"symbol": "C"} for i in range(6)}
    bonds = {(0, i): {"order": 1} for i in range(1, 6)}  # 5 bonds on atom 0
    result = _identify_problems_rdkit(atoms, bonds)
    assert 0 in result


# ---------------------------------------------------------------------------
# inject_ez_stereo_to_mol_block — no M END anchor
# ---------------------------------------------------------------------------


def test_inject_ez_stereo_no_ez_bonds_returns_original():
    from moleditpy.core.mol_geometry import inject_ez_stereo_to_mol_block

    mol = Chem.MolFromSmiles("CC")
    for atom in mol.GetAtoms():
        atom.SetIntProp("_original_atom_id", atom.GetIdx())
    block = "orig block\nM  END\n"
    # No stereo bonds -> unchanged
    assert inject_ez_stereo_to_mol_block(block, mol, {(0, 1): {"stereo": 0}}) == block


def test_inject_ez_stereo_appends_when_no_m_end():
    from moleditpy.core.mol_geometry import inject_ez_stereo_to_mol_block

    mol = Chem.MolFromSmiles("C/C=C/C")
    for atom in mol.GetAtoms():
        atom.SetIntProp("_original_atom_id", atom.GetIdx())
    # Double bond is atoms 1-2 in this SMILES
    block = "just a header without a terminator"
    out = inject_ez_stereo_to_mol_block(block, mol, {(1, 2): {"stereo": 3}})
    assert "M  CFG" in out


# ---------------------------------------------------------------------------
# resolve_2d_overlaps — bonded-skip and multi-fragment split
# ---------------------------------------------------------------------------


def test_resolve_2d_overlaps_skips_bonded_pairs():
    from moleditpy.core.mol_geometry import resolve_2d_overlaps

    positions = {1: (0.0, 0.0), 2: (0.1, 0.0)}  # overlapping
    adjacency = {1: [2], 2: [1]}
    # has_bond_check reports them bonded -> pair skipped -> no moves
    result = resolve_2d_overlaps(
        [1, 2], positions, adjacency, has_bond_check_func=lambda a, b: True
    )
    assert result == []


def test_resolve_2d_overlaps_moves_unbonded_fragment():
    from moleditpy.core.mol_geometry import resolve_2d_overlaps

    # Two disconnected fragments (1-2) and (3-4) whose atoms 2 and 3 overlap
    positions = {
        1: (0.0, 0.0),
        2: (10.0, 0.0),
        3: (10.1, 0.0),
        4: (20.0, 0.0),
    }
    adjacency = {1: [2], 2: [1], 3: [4], 4: [3]}
    result = resolve_2d_overlaps(
        [1, 2, 3, 4], positions, adjacency, has_bond_check_func=lambda a, b: False
    )
    assert len(result) == 1
    moved_ids, vector = result[0]
    # Only the overlapping atoms form the group; the higher-id side moves
    assert moved_ids == {3}
    assert vector == (-20, 20)


def test_get_connected_group_revisits_ring_atoms():
    from moleditpy.core.mol_geometry import get_connected_group

    mol = Chem.MolFromSmiles("C1CCCCC1")  # cyclohexane: BFS re-pops visited atoms
    assert get_connected_group(mol, 0) == {0, 1, 2, 3, 4, 5}


def test_adjust_dihedral_wraps_negative_delta():
    from moleditpy.core.mol_geometry import adjust_dihedral, calculate_dihedral

    positions = _np.array(
        [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
        dtype=float,
    )
    before = calculate_dihedral(positions, 0, 1, 2, 3)
    target = before - 350.0  # delta < -180 -> +360 branch
    adjust_dihedral(positions, 0, 1, 2, 3, target, [3])
    after = calculate_dihedral(positions, 0, 1, 2, 3)
    assert abs(((after - target + 180) % 360) - 180) < 1.0


def test_inject_ez_stereo_inserts_cfg_before_m_end():
    from moleditpy.core.mol_geometry import inject_ez_stereo_to_mol_block

    mol = Chem.MolFromSmiles("C/C=C/C")
    for atom in mol.GetAtoms():
        atom.SetIntProp("_original_atom_id", atom.GetIdx())
    block = "header\nbody\nM  END\n"
    out = inject_ez_stereo_to_mol_block(block, mol, {(1, 2): {"stereo": 4}})
    lines = out.split("\n")
    cfg_idx = next(i for i, l in enumerate(lines) if "M  CFG" in l)
    end_idx = next(i for i, l in enumerate(lines) if "M  END" in l)
    assert cfg_idx < end_idx  # CFG inserted just before the terminator


def test_identify_problems_rdkit_skips_bond_with_unknown_atom():
    from moleditpy.core.mol_geometry import _identify_problems_rdkit

    atoms = {1: {"symbol": "C"}}
    bonds = {(1, 99): {"order": 1}}  # atom 99 absent -> bond skipped, no crash
    assert _identify_problems_rdkit(atoms, bonds) == []


def test_resolve_2d_overlaps_splits_group_with_internal_bond():
    from moleditpy.core.mol_geometry import resolve_2d_overlaps

    # Atoms 1,2,3 all mutually overlap; 1-2 are bonded (one fragment), 3 is separate
    positions = {1: (0.0, 0.0), 2: (0.1, 0.0), 3: (0.2, 0.0)}
    adjacency = {1: [2], 2: [1], 3: []}
    result = resolve_2d_overlaps(
        [1, 2, 3], positions, adjacency, has_bond_check_func=lambda a, b: False
    )
    # The representative overlapping pair (1,2) lies within one bonded fragment,
    # so no fragment is moved.
    assert result == []


def test_calc_distance_basic():
    from moleditpy.core.mol_geometry import calc_distance

    assert calc_distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]) == pytest.approx(5.0)


def test_calc_angle_deg_right_angle():
    from moleditpy.core.mol_geometry import calc_angle_deg

    angle = calc_angle_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    assert angle == pytest.approx(90.0)


def test_resolve_2d_overlaps_single_fragment_group_not_moved():
    from moleditpy.core.mol_geometry import resolve_2d_overlaps

    # Two overlapping atoms that are bonded to each other form one fragment,
    # so the group cannot be split and nothing moves.
    positions = {1: (0.0, 0.0), 2: (0.1, 0.0)}
    adjacency = {1: [2], 2: [1]}
    result = resolve_2d_overlaps(
        [1, 2], positions, adjacency, has_bond_check_func=lambda a, b: False
    )
    assert result == []
