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
    positions = np.array([
        [1.0, 0.0, 0.0],  # A (idx 0)
        [0.0, 0.0, 0.0],  # B (idx 1)
        [0.0, 1.0, 0.0],  # C (idx 2)
    ])
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
    positions = np.array([
        [1.0, 0.0, 0.0],   # A
        [0.0, 0.0, 0.0],   # B (vertex)
        [0.0, 1.0, 0.0],   # C
        [0.0, 2.0, 0.0],   # D (attached to C)
    ])
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
    positions = np.array([
        [2.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0],  # B
        [-1.0, 0.0, 0.0], # C  (collinear with A-B, angle = 180°)
    ])
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
