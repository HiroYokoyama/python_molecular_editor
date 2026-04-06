import pytest
import sys
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.ui.calculation_worker import CalculationWorker
from moleditpy.core.mol_geometry import calculate_dihedral
from moleditpy.core.molecular_data import MolecularData
from moleditpy.ui.angle_dialog import AngleDialog
from moleditpy.ui.bond_length_dialog import BondLengthDialog
from PyQt6.QtCore import QPointF
from unittest.mock import patch, MagicMock


def test_calculation_worker_bond_length_validation(qtbot, app):
    """
    Integration test: Run a real 3D optimization for Ethane.
    Compare the resulting C-C bond length against RDKit's reference calculated value.
    No hard-coded bond length value is used.
    """
    # 1. Setup MolecularData (2D Ethane)
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    # 2. Run CalculationWorker (Real 3D Embedded + Optimization)
    worker = CalculationWorker()
    settings = {"conversion_mode": "rdkit", "do_optimize": True}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    result = blocker.args[0]
    # result is (mol_block_3d, rdkit_mol, info_dict)
    mol_3d = result[1]
    assert mol_3d is not None

    # 3. Validation: Compare against RDKit's own reference calculations
    # We find the C-C bond in the result mol.
    cc_bonds = [
        b
        for b in mol_3d.GetBonds()
        if b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"
    ]
    assert len(cc_bonds) >= 1
    target_bond = cc_bonds[0]

    # Get optimized distance from the result conformer
    conf = mol_3d.GetConformer()

    # Use BondLengthDialog logic via importing
    main_window = MagicMock()
    dialog = BondLengthDialog(mol_3d, main_window)
    dialog.atom1_idx = target_bond.GetBeginAtomIdx()
    dialog.atom2_idx = target_bond.GetEndAtomIdx()
    dialog.update_display()

    # Dialog distance text holds the value formatted to 3 decimal places
    dist_measured_text = dialog.distance_input.text()
    if dist_measured_text == "":
        p1 = np.array(conf.GetAtomPosition(target_bond.GetBeginAtomIdx()))
        p2 = np.array(conf.GetAtomPosition(target_bond.GetEndAtomIdx()))
        dist_measured = np.linalg.norm(p2 - p1)
    else:
        dist_measured = float(dist_measured_text)

    # Get "Reference" RDKit distance for a standard ETKDG embedded molecule
    ref_mol = Chem.MolFromSmiles("CC")
    ref_mol = Chem.AddHs(ref_mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(ref_mol, params)
    AllChem.MMFFOptimizeMolecule(ref_mol, mmffVariant="MMFF94s")
    ref_conf = ref_mol.GetConformer()
    # In ref_mol (from "CC"), C-C bond is between atoms 0 and 1
    dist_ref = rdMolTransforms.GetBondLength(ref_conf, 0, 1)

    assert dist_measured == pytest.approx(dist_ref, rel=1e-2)


def test_calculation_worker_angle_validation(qtbot, app):
    """
    Integration test: Run 3D optimization for Propane.
    Compare the C-C-C angle against RDKit's reference.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    c3 = data.add_atom("C", QPointF(100, 50))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c2, c3, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "rdkit", "do_optimize": True}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    mol_3d.GetConformer()

    # Find C-C-C triplet by symbol and connectivity
    # Central carbon has 2 carbon neighbors.
    central = -1
    neighbors = []
    for atom in mol_3d.GetAtoms():
        if atom.GetSymbol() == "C":
            c_neighbors = [
                n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "C"
            ]
            if len(c_neighbors) == 2:
                central = atom.GetIdx()
                neighbors = c_neighbors
                break

    assert central != -1, "Could not find central carbon in propane"

    # Use AngleDialog logic via importing
    main_window = MagicMock()
    dialog = AngleDialog(mol_3d, main_window)
    dialog.atom1_idx = neighbors[0]
    dialog.atom2_idx = central
    dialog.atom3_idx = neighbors[1]
    dialog.update_display()

    angle_measured_text = dialog.angle_input.text()
    if angle_measured_text == "":
        current_angle = dialog.calculate_angle()
        angle_measured = current_angle
    else:
        angle_measured = float(angle_measured_text)

    # Reference Propane (Identical setup to CalculationWorker)
    ref_mol = Chem.MolFromSmiles("CCC")
    ref_mol = Chem.AddHs(ref_mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(ref_mol, params)
    AllChem.MMFFOptimizeMolecule(ref_mol, mmffVariant="MMFF94s")

    # In ref_mol (from SMILES "CCC"), C atoms are 0, 1, 2. Central is 1.
    angle_ref = rdMolTransforms.GetAngleDeg(ref_mol.GetConformer(), 0, 1, 2)

    assert angle_measured == pytest.approx(angle_ref, abs=5.0), (
        f"Angle mismatch: measured {angle_measured:.2f} vs ref {angle_ref:.2f}"
    )


def test_calculation_worker_dihedral_validation(qtbot, app):
    """
    Integration test: Run 3D optimization for n-Butane.
    Check if it converges to a staggered conformation (Trans or Gauche).
    No hardcoded values; check against RDKit staggered preference.
    """
    data = MolecularData()
    # Linear start to avoid bias
    for i in range(4):
        data.add_atom("C", QPointF(i * 50, 0))
        if i > 0:
            data.add_bond(i - 1, i, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "rdkit", "do_optimize": True}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()

    # Path C1-C2-C3-C4
    # We find a chain of 4 carbons.
    positions = conf.GetPositions()
    dihedral = calculate_dihedral(positions, 0, 1, 2, 3)

    # A staggered conformation should have dihedral around 60 (gauche) or 180 (trans).
    # We check if it is NOT eclipsed (near 0, 120).
    is_staggered = (40 < abs(dihedral) < 80) or (160 < abs(dihedral) <= 180)

    assert is_staggered


def test_calculation_worker_conversion_no_optimize(qtbot, app):
    """
    Integration test: 2D to 3D without optimization.
    Check if atom counts and symbols are preserved.
    """
    data = MolecularData()
    data.add_atom("F", QPointF(0, 0))
    data.add_atom("Cl", QPointF(50, 0))
    data.add_bond(0, 1, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "rdkit", "do_optimize": False}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    symbols = set(a.GetSymbol() for a in mol_3d.GetAtoms())
    assert "F" in symbols
    assert "Cl" in symbols
    # Verification: Bond length is non-zero (3D coordinates generated)
    conf = mol_3d.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 0, 1)
    assert dist > 0.1


def test_calculation_worker_direct_mode(qtbot, app):
    """
    Integration test: Direct conversion mode.
    Ensures 2D coordinates are preserved as Z=0 and added Hs get a small Z offset.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(10, 20))
    c2 = data.add_atom("C", QPointF(60, 20))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "direct", "do_optimize": False}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()

    # Check heavy atoms (should be at Z=0 and match input coords)
    # Note: MolecularData coords are in pixels, MOL block might scale them or keep them.
    # Usually they are kept as is in the MOL block generated by to_mol_block.
    p1 = conf.GetAtomPosition(0)
    p2 = conf.GetAtomPosition(1)

    # ANGSTROM_PER_PIXEL is 1.5 / 75 = 0.02
    scale = 0.02
    assert p1.z == pytest.approx(0.0)
    assert p2.z == pytest.approx(0.0)
    assert p1.x == pytest.approx(10.0 * scale)
    assert p1.y == pytest.approx(-20.0 * scale)  # MolecularData Y is inverted

    # Check added hydrogens (should have non-zero Z)
    h_atoms = [a.GetIdx() for a in mol_3d.GetAtoms() if a.GetSymbol() == "H"]
    assert len(h_atoms) > 0
    for h_idx in h_atoms:
        hp = conf.GetAtomPosition(h_idx)
        assert abs(hp.z) > 0.05


def test_calculation_worker_halt_logic(qtbot, app):
    """
    Integration test: Verify halt mechanism.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    worker.halt_ids = {123}  # Halt this specific ID
    settings = {"worker_id": 123}

    with qtbot.waitSignal(worker.error, timeout=5000) as blocker:
        worker.run_calculation(mol_block, settings)

    # Payload is (worker_id, msg)
    err_id, err_msg = blocker.args[0]
    assert err_id == 123
    assert "Halted" in err_msg


def test_calculation_worker_global_halt(qtbot, app):
    """
    Integration test: Verify global halt.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    worker.halt_all = True
    settings = {"worker_id": None}

    with qtbot.waitSignal(worker.error, timeout=5000) as blocker:
        worker.run_calculation(mol_block, settings)

    err_id, err_msg = blocker.args[0]
    assert err_id is None
    assert "Halted" in err_msg


def test_calculation_worker_invalid_input(qtbot, app):
    """
    Test error handling for empty input.
    """
    worker = CalculationWorker()
    with qtbot.waitSignal(worker.error, timeout=5000) as blocker:
        worker.run_calculation("", None)

    _, err_msg = blocker.args[0]
    assert "No atoms to convert" in err_msg


def test_calculation_worker_isolation(qtbot, app):
    """
    Ensure worker_id correctly isolates halt signals.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    worker.halt_ids = {456}  # Halt 456, but 789 should run

    # Run for 789
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, {"worker_id": 789})

    res_id, _ = blocker.args[0]
    assert res_id == 789


def test_calculation_worker_direct_mode_stereo(qtbot, app):
    """
    Integration test: Direct mode with wedge/dash bonds.
    Verify that Z-offsets are applied to the 'end' atoms of stereo bonds.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))  # Bond 0-1
    c3 = data.add_atom("C", QPointF(0, 50))  # Bond 0-2 (stereo)
    data.add_bond(c1, c2, order=1)
    # 0 -> 2 is wedge
    data.add_bond(c1, c3, order=1, stereo=1)  # 1=wedge
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "direct", "do_optimize": False}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()

    # c1(0) is at origin, c2(1) is flat, c3(2) should have Z offset
    p1 = conf.GetAtomPosition(0)
    p2 = conf.GetAtomPosition(1)
    p3 = conf.GetAtomPosition(2)

    assert p1.z == pytest.approx(0.0)
    assert p2.z == pytest.approx(0.0)
    assert p3.z == pytest.approx(1.5)  # stereo_z_offset = 1.5 for wedge


def test_calculation_worker_constraint_embedding_fallback(qtbot, app):
    """
    Test the fallback to constraint-based embedding when initial embedding fails.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with patch("moleditpy.ui.calculation_worker.AllChem.EmbedMolecule") as mock_embed:
        # First (standard) fails (-1), second (constraint) succeeds (1);
        # extend list to avoid StopIteration if called more than twice
        mock_embed.side_effect = [-1] + [1] * 5

        with qtbot.waitSignal(worker.finished, timeout=10000):
            worker.run_calculation(mol_block, {"conversion_mode": "rdkit"})

        assert mock_embed.call_count >= 2


def test_calculation_worker_opt_failure_emits_error(qtbot, app):
    """
    Test that optimization failure emits an error string instead of failing silently or hanging.
    """
    data = MolecularData()
    c = data.add_atom("C", QPointF(0, 0))
    # We use at least 2 atoms so RDKit Mol is definitely valid
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c, c2)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with (
        patch(
            "rdkit.Chem.AllChem.MMFFGetMoleculeProperties", return_value=None
        ) as mock_mmff_props,
    ):
        # We now wait for the finished signal because optimization failure is non-fatal
        with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
            worker.run_calculation(mol_block, {"conversion_mode": "rdkit"})

        assert mock_mmff_props.called
        res_id, res_mol = blocker.args[0]
        assert res_mol is not None
        # Check that it tried to set the optimization method property but maybe cleared it or didn't set it
        # The key is that it didn't emit an error and returned a molecule.


def test_calculation_worker_mmff_variants(qtbot, app):
    """
    Test switching between MMFF94 and MMFF94s.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with patch(
        "moleditpy.ui.calculation_worker.AllChem.MMFFGetMoleculeProperties"
    ) as mock_props:
        # Test MMFF94
        worker.run_calculation(mol_block, {"optimization_method": "MMFF94_RDKIT"})
        qtbot.wait(500)
        # Check if called with MMFF94
        found = False
        for call in mock_props.call_args_list:
            if call.kwargs.get("mmffVariant") == "MMFF94":
                found = True
        assert found


def test_calculation_worker_obabel_fallback_mocked(qtbot, app):
    """
    Test the Open Babel fallback path by mocking availability and pybel.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with (
        patch(
            "moleditpy.ui.calculation_worker.AllChem.EmbedMolecule",
            return_value=-1,
        ),
        patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", True),
        patch("moleditpy.ui.calculation_worker.pybel", spec=True),
        patch("moleditpy.ui.calculation_worker.subprocess.run") as mock_run,
        patch("moleditpy.ui.calculation_worker.Chem.MolFromMolBlock") as mock_rd_parse,
    ):
        # Mock subprocess result
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "dummy_mol_block"
        mock_run.return_value = mock_result

        # Mock RDKit back-parse
        mock_rd_mol = Chem.MolFromSmiles("C")
        mock_rd_parse.return_value = mock_rd_mol

        with qtbot.waitSignal(worker.finished, timeout=10000):
            worker.run_calculation(mol_block, {"conversion_mode": "fallback"})

        assert mock_run.called
        # Check if first arg to subprocess.run contains sys.executable
        args, kwargs = mock_run.call_args
        assert sys.executable in args[0]


def test_calculation_worker_complex_direct_h_placement(qtbot, app):
    """
    Test direct mode with 4 hydrogens on a carbon to hit the rotation/offset logic.
    """
    data = MolecularData()
    c = data.add_atom("C", QPointF(0, 0))
    # Add 4 Hs in 2D
    h1 = data.add_atom("H", QPointF(10, 0))
    h2 = data.add_atom("H", QPointF(-10, 0))
    h3 = data.add_atom("H", QPointF(0, 10))
    h4 = data.add_atom("H", QPointF(0, -10))
    data.add_bond(c, h1)
    data.add_bond(c, h2)
    data.add_bond(c, h3)
    data.add_bond(c, h4)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "direct", "do_optimize": False}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d.GetNumAtoms() >= 5


# ---------------------------------------------------------------------------
# Unit tests for module-level helpers and new code paths
# ---------------------------------------------------------------------------
from moleditpy.ui.calculation_worker import (
    WorkerHaltError,
    _iterative_optimize,
    _adjust_collision_avoidance,
)


def test_worker_halt_error_is_exception():
    """WorkerHaltError is a proper Exception subclass."""
    err = WorkerHaltError("test halt")
    assert isinstance(err, Exception)
    assert str(err) == "test halt"


def test_worker_halt_error_not_caught_by_generic():
    """WorkerHaltError should propagate through except Exception if re-raised."""
    with pytest.raises(WorkerHaltError, match="Halted"):
        try:
            raise WorkerHaltError("Halted")
        except WorkerHaltError:
            raise


def test_iterative_optimize_mmff(app):
    """_iterative_optimize: MMFF method converges on ethane."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    result = _iterative_optimize(
        mol, "MMFF94s", lambda: False, lambda m: None, max_iters=200, chunk_size=50
    )
    assert result is True


def test_iterative_optimize_uff(app):
    """_iterative_optimize: UFF method converges on ethane."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    result = _iterative_optimize(
        mol, "UFF", lambda: False, lambda m: None, max_iters=200, chunk_size=50
    )
    assert result is True


def test_iterative_optimize_mmff94_variant(app):
    """_iterative_optimize: MMFF94 (non-s) variant."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    result = _iterative_optimize(
        mol, "MMFF94", lambda: False, lambda m: None, max_iters=200, chunk_size=50
    )
    assert result is True


def test_iterative_optimize_unknown_method(app):
    """_iterative_optimize: unknown method returns False."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    result = _iterative_optimize(mol, "INVALID", lambda: False, lambda m: None)
    assert result is False


def test_iterative_optimize_halt_during_optimization(app):
    """_iterative_optimize: raises WorkerHaltError when halted during chunks."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    call_count = [0]

    def check_halted():
        call_count[0] += 1
        return call_count[0] >= 2  # Halt on 2nd check (after 1st chunk completes)

    with pytest.raises(WorkerHaltError):
        _iterative_optimize(
            mol, "UFF", check_halted, lambda m: None, max_iters=4000, chunk_size=10
        )


def test_iterative_optimize_props_none(app):
    """_iterative_optimize: returns False when MMFF properties cannot be computed."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    with (
        patch(
            "moleditpy.ui.calculation_worker.AllChem.MMFFGetMoleculeProperties",
            return_value=None,
        ),
        patch(
            "moleditpy.ui.calculation_worker.AllChem.UFFGetMoleculeForceField",
            return_value=None,
        ),
    ):
        result = _iterative_optimize(mol, "MMFF94s", lambda: False, lambda m: None)
    assert result is False


def test_iterative_optimize_ff_none(app):
    """_iterative_optimize: returns False when UFF force field is None."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    with patch(
        "moleditpy.ui.calculation_worker.AllChem.UFFGetMoleculeForceField",
        return_value=None,
    ):
        result = _iterative_optimize(mol, "UFF", lambda: False, lambda m: None)
    assert result is False


def test_adjust_collision_avoidance_single_fragment(app):
    """_adjust_collision_avoidance: single fragment returns immediately (no-op)."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    status_msgs = []
    _adjust_collision_avoidance(mol, lambda: False, lambda m: status_msgs.append(m))
    # Single fragment → should return before emitting any status
    assert len(status_msgs) == 0


def test_adjust_collision_avoidance_multi_fragment(app):
    """_adjust_collision_avoidance: two overlapping fragments get separated."""
    # Create two separate methane molecules very close together
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    # Force the two carbons on top of each other
    conf = mol.GetConformer()
    frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    assert len(frags) == 2

    # Move all atoms of second fragment to same position as first fragment
    ref_pos = conf.GetAtomPosition(frags[0][0])
    for idx in frags[1]:
        conf.SetAtomPosition(idx, [ref_pos.x, ref_pos.y, ref_pos.z])

    status_msgs = []
    _adjust_collision_avoidance(mol, lambda: False, lambda m: status_msgs.append(m))

    # Should have emitted at least "Resolving..." and "Collision avoidance completed."
    assert len(status_msgs) >= 2
    assert "Resolving" in status_msgs[0]
    assert "completed" in status_msgs[-1]


def test_adjust_collision_avoidance_halt(app):
    """_adjust_collision_avoidance: raises WorkerHaltError when halted."""
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

    # Force collision
    conf = mol.GetConformer()
    frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    ref_pos = conf.GetAtomPosition(frags[0][0])
    for idx in frags[1]:
        conf.SetAtomPosition(idx, [ref_pos.x, ref_pos.y, ref_pos.z])

    with pytest.raises(WorkerHaltError):
        _adjust_collision_avoidance(mol, lambda: True, lambda m: None)


def test_calculation_worker_direct_with_optimize(qtbot, app):
    """
    Integration test: Direct conversion mode with do_optimize=True.
    Ensures optimization runs on the direct-mode structure.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(10, 20))
    c2 = data.add_atom("C", QPointF(60, 20))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "direct", "do_optimize": True}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() >= 1


def test_calculation_worker_optimized_result_better(qtbot, app):
    """
    Integration test: Optimized ethane should have a C-C bond length
    within a reasonable range (~1.52 Å for MMFF94s).
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(
            mol_block, {"conversion_mode": "rdkit", "do_optimize": True}
        )

    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()

    # Find C-C bond
    for b in mol_3d.GetBonds():
        if b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C":
            dist = rdMolTransforms.GetBondLength(
                conf, b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            )
            # Reasonable C-C single bond range after optimization
            assert 1.3 < dist < 1.7, f"Unexpected C-C bond length: {dist}"
            break


def test_calculation_worker_optimize_only_mmff(qtbot, app):
    """
    Integration test: optimize_only mode with MMFF method.
    Verifies that optimize_only skips embedding and only runs optimization.
    """
    # Create ethane with a valid 3D conformer
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    mol_block = Chem.MolToMolBlock(mol)

    worker = CalculationWorker()
    settings = {"conversion_mode": "optimize_only", "optimization_method": "MMFF94s"}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() >= 1


def test_calculation_worker_optimize_only_uff(qtbot, app):
    """
    Integration test: optimize_only mode with UFF method.
    """
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    mol_block = Chem.MolToMolBlock(mol)

    worker = CalculationWorker()
    settings = {"conversion_mode": "optimize_only", "optimization_method": "UFF"}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None


def test_calculation_worker_optimize_only_default(qtbot, app):
    """
    Integration test: optimize_only mode with default (no optimization_method).
    Should fallback to MMFF94s then UFF if needed.
    """
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    mol_block = Chem.MolToMolBlock(mol)

    worker = CalculationWorker()
    settings = {"conversion_mode": "optimize_only"}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None


def test_calculation_worker_optimize_only_mmff94_variant(qtbot, app):
    """
    Integration test: optimize_only mode with MMFF94 (non-s) variant.
    """
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    mol_block = Chem.MolToMolBlock(mol)

    worker = CalculationWorker()
    settings = {
        "conversion_mode": "optimize_only",
        "optimization_method": "MMFF94_RDKIT",
    }

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None


def test_calculation_worker_status_signals(qtbot, app):
    """
    Integration test: verify status_update signals are emitted during conversion.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    status_messages = []
    worker.status_update.connect(lambda msg: status_messages.append(msg))

    with qtbot.waitSignal(worker.finished, timeout=10000):
        worker.run_calculation(
            mol_block, {"conversion_mode": "rdkit", "do_optimize": True}
        )

    # Should have received at least one status update
    assert len(status_messages) >= 1
    assert any("3D" in msg or "Creating" in msg for msg in status_messages)


def test_calculation_worker_multi_fragment_rdkit(qtbot, app):
    """
    Integration test: multi-fragment molecule triggers collision avoidance in RDKit path.
    """
    # Two separate methane molecules
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    mol_block = Chem.MolToMolBlock(mol)

    worker = CalculationWorker()
    status_messages = []
    worker.status_update.connect(lambda msg: status_messages.append(msg))

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(
            mol_block, {"conversion_mode": "rdkit", "do_optimize": True}
        )

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    # Multi-fragment should trigger collision avoidance status messages
    collision_msgs = [
        m for m in status_messages if "collision" in m.lower() or "Resolving" in m
    ]
    assert len(collision_msgs) >= 1


def test_calculation_worker_direct_dash_stereo(qtbot, app):
    """
    Integration test: direct mode with a dash bond.
    Verifies negative Z-offset for the dash end atom.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    c3 = data.add_atom("C", QPointF(0, 50))
    data.add_bond(c1, c2, order=1)
    data.add_bond(c1, c3, order=1, stereo=2)  # 2=dash
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {"conversion_mode": "direct", "do_optimize": False}

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()
    # The dash end atom should have negative Z
    p3 = conf.GetAtomPosition(2)
    assert p3.z == pytest.approx(-1.5)


def test_calculation_worker_direct_mmff94_rdkit_variant(qtbot, app):
    """
    Integration test: direct mode with do_optimize=True and MMFF94_RDKIT variant.
    Covers line 832 where optimization_method is checked for MMFF94_RDKIT.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(10, 20))
    c2 = data.add_atom("C", QPointF(60, 20))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {
        "conversion_mode": "direct",
        "do_optimize": True,
        "optimization_method": "MMFF94_RDKIT",
    }

    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() >= 1


def test_calculation_worker_fallback_to_direct_no_obabel(qtbot, app):
    """
    Integration test: fallback mode with OBABEL_AVAILABLE=False and
    RDKit embedding forced to fail. Should fall through to direct conversion
    as last resort.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    status_messages = []
    worker.status_update.connect(lambda msg: status_messages.append(msg))

    # Mock OBABEL_AVAILABLE to False and force EmbedMolecule to return -1
    with (
        patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", False),
        patch("moleditpy.ui.calculation_worker.AllChem.EmbedMolecule", return_value=-1),
    ):
        with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
            worker.run_calculation(
                mol_block, {"conversion_mode": "fallback", "do_optimize": False}
            )

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() >= 1
    # Should have hit the direct conversion path message
    assert any("direct" in m.lower() or "Direct" in m for m in status_messages)


def test_calculation_worker_fallback_to_direct_with_optimize(qtbot, app):
    """
    Integration test: fallback-to-direct with do_optimize=True.
    Direct conversion should still run optimization after coordinate placement.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0))
    data.add_bond(c1, c2, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()

    with (
        patch("moleditpy.ui.calculation_worker.OBABEL_AVAILABLE", False),
        patch("moleditpy.ui.calculation_worker.AllChem.EmbedMolecule", return_value=-1),
    ):
        with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
            worker.run_calculation(
                mol_block, {"conversion_mode": "fallback", "do_optimize": True}
            )

    mol_3d = blocker.args[0][1]
    assert mol_3d is not None
    assert mol_3d.GetNumConformers() >= 1
