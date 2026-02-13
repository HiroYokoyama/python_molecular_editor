import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from moleditpy.modules.calculation_worker import CalculationWorker
from moleditpy.modules.molecular_data import MolecularData
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
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    result = blocker.args[0]
    # result is (mol_block_3d, rdkit_mol, info_dict)
    mol_3d = result[1]
    assert mol_3d is not None
    
    # 3. Validation: Compare against RDKit's own reference calculations
    # We find the C-C bond in the result mol.
    cc_bonds = [b for b in mol_3d.GetBonds() if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C']
    assert len(cc_bonds) >= 1
    target_bond = cc_bonds[0]
    
    # Get optimized distance from the result conformer
    conf = mol_3d.GetConformer()
    dist_measured = rdMolTransforms.GetBondLength(conf, target_bond.GetBeginAtomIdx(), target_bond.GetEndAtomIdx())
    
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
    
    print(f"DEBUG: Measured C-C = {dist_measured:.4f}, Reference C-C = {dist_ref:.4f}")
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
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()
    
    # Find C-C-C triplet by symbol and connectivity
    # Central carbon has 2 carbon neighbors.
    central = -1
    neighbors = []
    for atom in mol_3d.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(c_neighbors) == 2:
                central = atom.GetIdx()
                neighbors = c_neighbors
                break
    
    assert central != -1, "Could not find central carbon in propane"
    angle_measured = rdMolTransforms.GetAngleDeg(conf, neighbors[0], central, neighbors[1])
    
    # Reference Propane (Identical setup to CalculationWorker)
    ref_mol = Chem.MolFromSmiles("CCC")
    ref_mol = Chem.AddHs(ref_mol)
    params = AllChem.ETKDGv2()
    params.randomSeed = 42
    AllChem.EmbedMolecule(ref_mol, params)
    AllChem.MMFFOptimizeMolecule(ref_mol, mmffVariant="MMFF94s")
    
    # In ref_mol (from SMILES "CCC"), C atoms are 0, 1, 2. Central is 1.
    angle_ref = rdMolTransforms.GetAngleDeg(ref_mol.GetConformer(), 0, 1, 2)
    
    print(f"DEBUG: Measured C-C-C = {angle_measured:.4f}, Reference = {angle_ref:.4f}")
    assert angle_measured == pytest.approx(angle_ref, abs=5.0), \
        f"Angle mismatch: measured {angle_measured:.2f} vs ref {angle_ref:.2f}"

def test_calculation_worker_dihedral_validation(qtbot, app):
    """
    Integration test: Run 3D optimization for n-Butane.
    Check if it converges to a staggered conformation (Trans or Gauche).
    No hardcoded values; check against RDKit staggered preference.
    """
    data = MolecularData()
    # Linear start to avoid bias
    for i in range(4):
        data.add_atom("C", QPointF(i*50, 0))
        if i > 0:
            data.add_bond(i-1, i, order=1)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'rdkit', 'do_optimize': True}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    
    mol_3d = blocker.args[0][1]
    conf = mol_3d.GetConformer()
    
    # Path C1-C2-C3-C4
    # We find a chain of 4 carbons.
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3)
    
    # A staggered conformation should have dihedral around 60 (gauche) or 180 (trans).
    # We check if it is NOT eclipsed (near 0, 120).
    is_staggered = (40 < abs(dihedral) < 80) or (160 < abs(dihedral) <= 180)
    
    print(f"DEBUG: n-Butane Dihedral = {dihedral:.4f}")
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
    settings = {'conversion_mode': 'rdkit', 'do_optimize': False}
    
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
    settings = {'conversion_mode': 'direct'}
    
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
    assert p1.y == pytest.approx(-20.0 * scale) # MolecularData Y is inverted
    
    # Check added hydrogens (should have non-zero Z)
    h_atoms = [a.GetIdx() for a in mol_3d.GetAtoms() if a.GetSymbol() == 'H']
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
    worker.halt_ids = {123} # Halt this specific ID
    settings = {'worker_id': 123}
    
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
    settings = {'worker_id': None}
    
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
    worker.halt_ids = {456} # Halt 456, but 789 should run
    
    # Run for 789
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, {'worker_id': 789})
    
    res_id, _ = blocker.args[0]
    assert res_id == 789

def test_calculation_worker_direct_mode_stereo(qtbot, app):
    """
    Integration test: Direct mode with wedge/dash bonds.
    Verify that Z-offsets are applied to the 'end' atoms of stereo bonds.
    """
    data = MolecularData()
    c1 = data.add_atom("C", QPointF(0, 0))
    c2 = data.add_atom("C", QPointF(50, 0)) # Bond 0-1
    c3 = data.add_atom("C", QPointF(0, 50)) # Bond 0-2 (stereo)
    data.add_bond(c1, c2, order=1)
    # 0 -> 2 is wedge
    data.add_bond(c1, c3, order=1, stereo=1) # 1=wedge
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'direct'}
    
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
    assert p3.z == pytest.approx(1.5) # stereo_z_offset = 1.5 for wedge

def test_calculation_worker_constraint_embedding_fallback(qtbot, app):
    """
    Test the fallback to constraint-based embedding when initial embedding fails.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    
    with patch("moleditpy.modules.calculation_worker.AllChem.EmbedMolecule") as mock_embed:
        # First (standard) fails (-1), second (constraint) succeeds (1)
        mock_embed.side_effect = [-1, 1]
        
        with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
            worker.run_calculation(mol_block, {'conversion_mode': 'rdkit'})
        
        assert mock_embed.call_count >= 2

def test_calculation_worker_uff_fallback(qtbot, app):
    """
    Test fallback to UFF when MMFF optimization fails.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    
    with patch("moleditpy.modules.calculation_worker.AllChem.MMFFOptimizeMolecule") as mock_mmff, \
         patch("moleditpy.modules.calculation_worker.AllChem.UFFOptimizeMolecule") as mock_uff:
        
        mock_mmff.side_effect = Exception("MMFF Failed")
        
        with qtbot.waitSignal(worker.finished, timeout=10000):
            worker.run_calculation(mol_block, {'conversion_mode': 'rdkit'})
        
        assert mock_mmff.called
        assert mock_uff.called

def test_calculation_worker_mmff_variants(qtbot, app):
    """
    Test switching between MMFF94 and MMFF94s.
    """
    data = MolecularData()
    data.add_atom("C", QPointF(0, 0))
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    
    with patch("moleditpy.modules.calculation_worker.AllChem.MMFFOptimizeMolecule") as mock_mmff:
        # Test MMFF94
        worker.run_calculation(mol_block, {'optimization_method': 'MMFF94_RDKIT'})
        qtbot.wait(500)
        # Check if called with MMFF94
        found = False
        for call in mock_mmff.call_args_list:
            if call.kwargs.get('mmffVariant') == 'MMFF94':
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
    
    with patch("moleditpy.modules.calculation_worker.AllChem.EmbedMolecule", return_value=-1), \
         patch("moleditpy.modules.calculation_worker.OBABEL_AVAILABLE", True), \
         patch("moleditpy.modules.calculation_worker.pybel") as mock_pybel, \
         patch("moleditpy.modules.calculation_worker.Chem.MolFromMolBlock") as mock_rd_parse:
        
        mock_ob_mol = MagicMock()
        mock_pybel.readstring.return_value = mock_ob_mol
        mock_ob_mol.write.return_value = "dummy_mol_block"
        
        # Mock RDKit back-parse
        mock_rd_mol = Chem.MolFromSmiles("C")
        mock_rd_parse.return_value = mock_rd_mol
        
        with qtbot.waitSignal(worker.finished, timeout=10000):
            worker.run_calculation(mol_block, {'conversion_mode': 'fallback'})
        
        assert mock_pybel.readstring.called
        assert mock_ob_mol.make3D.called

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
    data.add_bond(c, h1); data.add_bond(c, h2); data.add_bond(c, h3); data.add_bond(c, h4)
    mol_block = data.to_mol_block()

    worker = CalculationWorker()
    settings = {'conversion_mode': 'direct'}
    
    with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
        worker.run_calculation(mol_block, settings)
    
    mol_3d = blocker.args[0][1]
    assert mol_3d.GetNumAtoms() >= 5

# def test_calculation_worker_signal_type_error_fallback(qtbot, app):
#     """
#     Test the fallback in _safe_finished when a TypeError occurs (simulating old signal signature).
#     """
#     data = MolecularData()
#     data.add_atom("C", QPointF(0, 0))
#     mol_block = data.to_mol_block()
# 
#     worker = CalculationWorker()
#     
#     with patch.object(worker.finished, 'emit') as mock_emit:
#         # First call fails with TypeError, second call (fallback) should pass
#         mock_emit.side_effect = [TypeError("Simulated signature mismatch"), None]
#         
#         with qtbot.waitSignal(worker.finished, timeout=10000) as blocker:
#             worker.run_calculation(mol_block, {'conversion_mode': 'rdkit'})
#         
#         assert mock_emit.call_count >= 2
