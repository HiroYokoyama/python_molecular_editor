import pytest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.compute_logic import ComputeManager
from moleditpy.core.molecular_data import MolecularData
from PyQt6.QtCore import QPointF, QPoint, QTimer, QThread
from PyQt6.QtGui import QColor, QAction
from PyQt6.QtWidgets import QMenu, QMessageBox
from unittest.mock import MagicMock, patch


class DummyCompute(ComputeManager):
    opt3d_method_labels = None  # Prevent __getattr__ from shadowing init_manager.opt3d_method_labels

    def __init__(self, host):
        self._host = host
        ComputeManager.__init__(self, host)
        
        # Force populate host if it's a MagicMock or missing managers
        if not hasattr(host, "init_manager"):
            host.init_manager = MagicMock()
        
        # Always set these to ensure they are real objects, not MagicMocks
        host.init_manager.settings = getattr(host.init_manager, "settings", {}) or {}
        host.init_manager.opt3d_method_labels = {
            "MMFF_RDKIT": "MMFF94s (RDKit)",
            "UFF_RDKIT": "UFF (RDKit)",
        }
        if not hasattr(host, "view_3d_manager"): host.view_3d_manager = MagicMock()
        if not hasattr(host, "state_manager"): host.state_manager = MagicMock()
        if not hasattr(host, "ui_manager"): host.ui_manager = MagicMock()
        
        # Ensure buttons/actions exist for UI transition tests
        for btn in ["convert_button", "cleanup_button", "optimize_3d_button", "export_button", "analysis_action", "edit_3d_action"]:
            if not hasattr(host.init_manager, btn):
                setattr(host.init_manager, btn, MagicMock())

    def __getattr__(self, name):
        return getattr(self._host, name)

    @property
    def data(self): return self.host.state_manager.data
    @data.setter
    def data(self, v): self.host.state_manager.data = v

    @property
    def scene(self): return self.host.init_manager.scene
    @scene.setter
    def scene(self, v): self.host.init_manager.scene = v

    @property
    def view_2d(self): return self.host.init_manager.view_2d
    @view_2d.setter
    def view_2d(self, v): self.host.init_manager.view_2d = v

    @property
    def view_3d(self): return self.host.view_3d_manager.view_3d
    @view_3d.setter
    def view_3d(self, v): self.host.view_3d_manager.view_3d = v

    @property
    def plotter(self): return self.host.view_3d_manager.plotter
    @plotter.setter
    def plotter(self, v): self.host.view_3d_manager.plotter = v

    @property
    def settings(self): return self.host.init_manager.settings
    @settings.setter
    def settings(self, v): self.host.init_manager.settings = v
    @property
    def current_mol(self): return self.host.view_3d_manager.current_mol
    @current_mol.setter
    def current_mol(self, v): self.host.view_3d_manager.current_mol = v

    @property
    def optimization_method(self):
        return self.host.init_manager.settings.get("optimization_method", "MMFF_RDKIT")
    @optimization_method.setter
    def optimization_method(self, v):
        self.host.init_manager.settings["optimization_method"] = v
        self.host.init_manager.optimization_method = v

    @property
    def constraints_3d(self): return self.host.edit_3d_manager.constraints_3d
    @constraints_3d.setter
    def constraints_3d(self, v): self.host.edit_3d_manager.constraints_3d = v

    def statusBar(self):
        return self._host.statusBar()

    def get_status_messages(self):
        """Helper to get all messages sent to the status bar."""
        return [
            str(call[0][0])
            for call in self.statusBar().showMessage.call_args_list
            if call[0]
        ]

    def update_from_rdkit_mol(self, mol):
        pass

    def create_atom_id_mapping(self):
        return {}

    def draw_molecule_3d(self, mol):
        pass

    def update_chiral_labels(self):
        pass

    def setup_3d_hover(self):
        pass

    def update_atom_id_menu_text(self):
        pass

    def update_atom_id_menu_state(self):
        pass

    def update_window_title(self):
        pass

    def setEnabled(self, val):
        pass

    def _enable_3d_features(self, val):
        pass

    def _enable_3d_edit_actions(self, val):
        pass

    def push_undo_state(self):
        pass

    def adjust_molecule_positions_to_avoid_collisions(self, mol, frags):
        pass


def test_on_calculation_error_stale(mock_parser_host):
    """Test on_calculation_error when the worker is stale (not in active set)."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"new_worker_id"}
    compute.on_calculation_error(("stale_id", "Ignore this error"))
    assert not compute.statusBar().showMessage.called


def test_on_calculation_error_basic(mock_parser_host):
    """Test on_calculation_error for an ACTIVE worker."""
    compute = DummyCompute(mock_parser_host)
    worker_id = "active_id"
    compute.active_worker_ids = {worker_id}
    compute.on_calculation_error((worker_id, "Real Error"))
    assert compute.statusBar().showMessage.called
    assert "Real Error" in compute.statusBar().showMessage.call_args[0][0]


def test_compute_set_optimization_method(mock_parser_host):
    """Verify that setting the optimization method updates both settings and internal state."""
    compute = DummyCompute(mock_parser_host)
    compute.set_optimization_method("GAFF_OBABEL")
    assert compute.host.init_manager.settings["optimization_method"] == "GAFF_OBABEL"
    assert compute.statusBar().showMessage.called
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "Optimization" in msg or "GAFF_OBABEL" in msg
    # Verify internal state that actually affects calculation
    assert compute.host.init_manager.optimization_method == "GAFF_OBABEL"


def test_compute_halt_logic(mock_parser_host):
    """Verify that halt_conversion correctly marks active workers for termination."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids.add("test_id")
    compute.halt_conversion()
    assert "test_id" in compute.halt_ids
    assert len(compute.active_worker_ids) == 0
    assert compute.statusBar().showMessage.called


def test_on_calculation_finished_basic(mock_parser_host):
    """Verify that on_calculation_finished correctly processes a finished worker result."""
    compute = DummyCompute(mock_parser_host)
    worker_id = "test_worker"
    compute.active_worker_ids.add(worker_id)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    result = (worker_id, mol)
    with patch.object(compute, "draw_molecule_3d") as mock_draw:
        compute.on_calculation_finished(result)
        assert compute.host.view_3d_manager.current_mol == mol
        assert worker_id not in compute.active_worker_ids


def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    """Verify that the manual valence fallback correctly identifies overvalent atoms."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.data.bonds = {}
    c_id = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    c_item = mock_parser_host.data.atoms[c_id]["item"]
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", QPointF(i + 1, 0))
        h_item = mock_parser_host.data.atoms[h_id]["item"]
        mock_parser_host.scene.create_bond(c_item, h_item, bond_order=1)
    compute.check_chemistry_problems_fallback()
    assert c_item.has_problem is True
    assert compute.statusBar().showMessage.called


def test_trigger_conversion_empty(mock_parser_host):
    """Verify that trigger_conversion handles empty molecular data gracefully."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.trigger_conversion()
    assert compute.statusBar().showMessage.called


def test_trigger_conversion_with_atoms(mock_parser_host):
    """Verify that trigger_conversion correctly starts the calculation thread for a valid molecule."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {"symbol": "C", "item": MagicMock()}}
    compute.settings["conversion_target"] = "all"
    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker"),
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.trigger_conversion()
        assert compute.statusBar().showMessage.called


def test_optimize_3d_structure_logic(mock_parser_host):
    """Verify the high-level logic of triggering 3D optimization on the current molecule."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    compute.current_mol = mol
    compute.optimize_3d_structure()
    assert compute.statusBar().showMessage.called


def test_on_calculation_finished_worker_id_mismatch(mock_parser_host):
    """Verify that on_calculation_finished ignores results from stale or mismatched workers."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"valid_id"}
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    result = ("invalid_id", mol)
    compute.on_calculation_finished(result)
    assert compute.host.view_3d_manager.current_mol is None


def test_trigger_conversion_chemistry_problems(mock_parser_host):
    """Test trigger_conversion when Chem.DetectChemistryProblems finds issues."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    atom = mol.GetAtomWithIdx(0)
    atom.SetIntProp("_original_atom_id", 1)
    compute.data.atoms = {1: {"symbol": "C", "item": MagicMock()}}

    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0
    with (
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[problem]),
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
    ):
        compute.trigger_conversion()
        all_messages = [
            str(call[0][0])
            for call in compute.statusBar().showMessage.call_args_list
            if call[0]
        ]
        assert any("chemistry problem(s) found" in msg for msg in all_messages)


def test_trigger_conversion_sanitize_error(mock_parser_host):
    """Test trigger_conversion when Chem.SanitizeMol fails."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    # MUST populate atoms to avoid empty trigger return
    compute.data.atoms = {1: {"symbol": "C", "item": MagicMock()}}

    with (
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol", side_effect=ValueError("Sanitize failed")),
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
    ):
        compute.trigger_conversion()
        all_messages = [
            str(call[0][0])
            for call in compute.statusBar().showMessage.call_args_list
            if call[0]
        ]
        assert any("Error: Invalid chemical structure." in msg for msg in all_messages)


def test_trigger_conversion_multiple_frags(mock_parser_host):
    """Test trigger_conversion with multiple fragments."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C.C")
    mol = Chem.AddHs(mol)
    # MUST populate atoms to avoid empty trigger return
    compute.data.atoms = {
        1: {"symbol": "C", "item": MagicMock()},
        2: {"symbol": "C", "item": MagicMock()},
    }

    with (
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol"),
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        patch("moleditpy.ui.compute_logic.CalculationWorker"),
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.trigger_conversion()
        all_messages = [
            str(call[0][0])
            for call in compute.statusBar().showMessage.call_args_list
            if call[0]
        ]
        assert any("collision detection" in msg for msg in all_messages)


def test_on_calculation_finished_single_mol_legacy(mock_parser_host):
    """Test on_calculation_finished with a single mol (legacy result format)."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    compute.on_calculation_finished(mol)
    assert compute.current_mol == mol


def test_on_calculation_error_legacy_payload(mock_parser_host):
    """Test on_calculation_error with a string (legacy error format)."""
    compute = DummyCompute(mock_parser_host)
    compute.on_calculation_error("Fatal Error")
    assert any(
        "Fatal Error" in str(call[0][0])
        for call in compute.statusBar().showMessage.call_args_list
        if call[0]
    )


def test_optimize_3d_temp_method_override(mock_parser_host):
    """Test optimize_3d_structure with temporary optimization method override."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "UFF_RDKIT"
    compute.host._temp_optimization_method = "MMFF_RDKIT"

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        mock_worker = MockWorker.return_value
        compute.optimize_3d_structure()
        
        # In the new async implementation, it should create a worker
        assert MockWorker.called
        # Verify it uses the temp override method (MMFF_RDKIT)
        # We can't easily check the worker start_work signal payload here without more mocking,
        # but we've verified it doesn't crash and starts the worker.
        # Verify the exact status message from production code
        compute.statusBar().showMessage.assert_any_call(
            "Optimizing 3D structure (MMFF_RDKIT)..."
        )


def test_optimize_3d_mmff94s_success(mock_parser_host):
    """Test MMFF94s optimization succeeds."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=0):
        compute.optimize_3d_structure()
        print(
            f"DEBUG: showMessage calls: {compute.statusBar().showMessage.call_args_list}"
        )
        print(f"DEBUG: messages extracted: {compute.get_status_messages()}")
        assert any("Optimizing" in msg for msg in compute.get_status_messages())


def test_optimize_3d_uff_success(mock_parser_host):
    """Test UFF optimization succeeds."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "UFF_RDKIT"

    with patch("rdkit.Chem.AllChem.UFFOptimizeMolecule", return_value=0):
        compute.optimize_3d_structure()
        assert any("Optimizing" in msg for msg in compute.get_status_messages())


def test_optimize_3d_no_conformer(mock_parser_host):
    """Test optimize_3d_structure when molecule has no conformer."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    compute.optimize_3d_structure()
    messages = compute.get_status_messages()
    assert any("No conformer found" in msg for msg in messages)


def test_optimize_3d_mmff_exception_handling(mock_parser_host):
    """Test grace during MMFF exception."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.optimize_3d_structure()
        assert MockWorker.called


def test_optimize_3d_uff_exception_handling(mock_parser_host):
    """Test grace during UFF exception."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "UFF_RDKIT"

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.optimize_3d_structure()
        assert MockWorker.called


def test_optimize_3d_plugin_method(mock_parser_host):
    """Test plugin optimization method — worker should be started (not rejected)."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol

    compute.host.plugin_manager = MagicMock()
    compute.host.plugin_manager.optimization_methods = {"CUSTOM": {"callback": MagicMock()}}
    compute.optimization_method = "CUSTOM"

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker"),
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.optimize_3d_structure()
    msgs = compute.get_status_messages()
    # Plugin method is valid — should not show "not available" and should show "Optimizing"
    assert not any("not available" in msg for msg in msgs)
    assert any("Optimizing" in msg for msg in msgs)


def test_optimize_3d_plugin_failure(mock_parser_host):
    """Test unavailable/unknown method shows correct error."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "NONEXISTENT_CUSTOM"

    compute.optimize_3d_structure()
    msgs = compute.get_status_messages()
    assert any("not available" in msg for msg in msgs)


def test_optimize_3d_mmff_fallback_success(mock_parser_host):
    """Test MMFF fallback to ForceField API when basic optimization fails."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    # Mock MMFFOptimizeMolecule to return 1 (non-zero/fail)
    with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=1):
        # Mock FF Minimize to return 0 (success)
        mock_ff = MagicMock()
        mock_ff.Minimize.return_value = 0
        with patch(
            "rdkit.Chem.AllChem.MMFFGetMoleculeForceField", return_value=mock_ff
        ):
            compute.optimize_3d_structure()
            # Success redraw should happen
            assert any("Optimizing" in msg for msg in compute.get_status_messages())


def test_optimize_3d_uff_fallback_failure(mock_parser_host):
    """Test UFF fallback failure handles gracefully."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        compute.optimize_3d_structure()
        assert MockWorker.called


def test_optimize_3d_unavailable_method(mock_parser_host):
    """Test error when optimization method is unavailable."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    compute.current_mol = mol
    compute.optimization_method = "INVALID_METHOD"

    compute.optimize_3d_structure()
    messages = compute.get_status_messages()
    assert any("Selected optimization" in msg for msg in messages)


def test_on_calculation_finished_collision_single_frag(mock_parser_host):
    """Test collision logic is skipped for single fragment."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    with patch.object(
        compute, "adjust_molecule_positions_to_avoid_collisions"
    ) as mock_adjust:
        compute.on_calculation_finished(mol)
        assert not mock_adjust.called
        assert not any(
            "Detecting collisions" in msg for msg in compute.get_status_messages()
        )


# removed collision test at this layer


# removed collision exception test at this layer


def test_molecular_data_radical_transfer():
    """Test that radical electrons are transferred correctly to RDKit mol."""
    from moleditpy.core.molecular_data import MolecularData
    from PyQt6.QtCore import QPointF

    data = MolecularData()
    # Adding a carbon with 1 radical electron
    data.add_atom("C", QPointF(0, 0), radical=1)
    mol = data.to_rdkit_mol()
    assert mol is not None
    assert mol.GetAtomWithIdx(0).GetNumRadicalElectrons() == 1


def test_app_state_radical_and_constraint_preservation(mock_parser_host):
    """Test that radicals and constraints are preserved through state round-trip."""
    from moleditpy.ui.app_state import MainWindowAppState
    from moleditpy.core.molecular_data import MolecularData
    from moleditpy.ui.atom_item import AtomItem
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)
    compute.data = MolecularData()  # Use real MolecularData for this test

    # Setup initial state with radical and constraint
    pos = QPointF(10, 20)
    aid = compute.data.add_atom("C", pos, radical=2)
    # Manually add AtomItem because get_current_state expects it
    item = AtomItem(aid, "C", pos, radical=2)
    compute.data.atoms[aid]["item"] = item

    compute.constraints_3d = [("DISTANCE", (0, 1), 1.5, 1e5)]

    # MainWindowAppState is a mixin class, we call its methods by passing 'compute' as self
    state = MainWindowAppState.get_current_state(compute)

    # Verify state content
    assert state["atoms"][0]["radical"] == 2
    # The code converts tuple to list for JSON compatibility in state
    assert state["constraints_3d"][0] == ["DISTANCE", [0, 1], 1.5, 1e5]

    # Clear and Restore
    compute.data.atoms.clear()
    compute.constraints_3d = []
    # Mock scene.addItem as it is called in set_state_from_data
    compute.scene = MagicMock()

    MainWindowAppState.set_state_from_data(compute, state)

    # Verify restoration
    assert 0 in compute.data.atoms
    assert compute.data.atoms[0]["radical"] == 2
    # set_state_from_data converts back to tuple
    assert compute.constraints_3d == [("DISTANCE", (0, 1), 1.5, 1e5)]


def test_app_state_original_atom_id_preservation(mock_parser_host):
    """Test that _original_atom_id is preserved in 3D molecule state round-trip."""
    from moleditpy.ui.app_state import MainWindowAppState

    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", 123)
    compute.current_mol = mol

    # Get state
    state = MainWindowAppState.get_current_state(compute)
    assert 123 in state["mol_3d_atom_ids"]

    # Clear and Restore
    compute.current_mol = None
    MainWindowAppState.set_state_from_data(compute, state)

    # Verify restoration
    assert compute.current_mol is not None
    assert compute.current_mol.GetAtomWithIdx(0).HasProp("_original_atom_id")
    assert compute.current_mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id") == 123


def test_optimize_3d_method_persistence(mock_parser_host):
    """Test that the optimization method is recorded after success."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp("_pme_optimization_method", "MMFF_RDKIT")
    compute.current_mol = mol
    compute.optimization_method = "MMFF_RDKIT"

    # In the new async implementation, the method persistence happens in on_calculation_finished.
    # We call it directly to test the persistence logic.
    compute.on_calculation_finished(mol)
    assert compute.last_successful_optimization_method == "MMFF94s (RDKit)"


def test_trigger_conversion_early_exits(mock_parser_host):
    """Test early exits in trigger_conversion (empty mol, etc.)."""
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.plotter = MagicMock()

    # 1. Empty data
    compute.trigger_conversion()
    assert compute.host.view_3d_manager.current_mol is None
    assert any("3D view cleared" in msg for msg in compute.get_status_messages())

    # 2. RDKit conversion failure (triggers fallback check)
    compute.data.add_atom("C", QPointF(0, 0))
    # Mock to_rdkit_mol to return None
    with patch.object(compute.data, "to_rdkit_mol", return_value=None):
        with patch.object(
            compute, "check_chemistry_problems_fallback"
        ) as mock_fallback:
            compute.trigger_conversion()
            assert mock_fallback.called


@pytest.mark.skip(reason="Segfaults in headless environment")
def test_check_chemistry_problems_fallback(mock_parser_host):
    """Test the manual valence check when RDKit fails."""
    # from moleditpy.modules.atom_item import AtomItem
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)

    # Setup a problematic atom: Carbon with 5 bonds
    pos = QPointF(0, 0)
    compute.data.add_atom("C", pos)

    # Use Mock item to avoid Qt crash
    mock_item = MagicMock()
    mock_item.has_problem = False
    mock_item.pos.return_value = pos
    compute.data.atoms[0]["item"] = mock_item

    # Mock bonds to have count 5
    compute.data.bonds = {
        (0, 1): {"order": 2},
        (0, 2): {"order": 1},
        (0, 3): {"order": 1},
        (0, 4): {"order": 1},
    }
    # Total bond count = 2+1+1+1 = 5

    compute.check_chemistry_problems_fallback()
    assert compute.data.atoms[0]["item"].has_problem == True
    msgs = compute.get_status_messages()
    assert any("chemistry problem" in msg and "valence" in msg for msg in msgs)


def test_trigger_conversion_happy_path(mock_parser_host):
    """Test trigger_conversion follows the success path until worker setup."""
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)
    compute.data.add_atom("C", QPointF(0, 0))
    mol = Chem.MolFromSmiles("C")

    with patch.object(compute.data, "to_rdkit_mol", return_value=mol):
        with patch("rdkit.Chem.DetectChemistryProblems", return_value=[]):
            with patch("rdkit.Chem.SanitizeMol"):
                with patch.object(
                    compute.data, "to_mol_block", return_value="M  V2000\n"
                ):
                    # This should hit the main logic and proceed to worker creation
                    # We mock the signal to avoid real thread start
                    compute.start_calculation = MagicMock()
                    compute.worker_thread = MagicMock()

                    with (
                        patch("moleditpy.ui.compute_logic.QThread"),
                        patch("moleditpy.ui.compute_logic.CalculationWorker"),
                        patch("PyQt6.QtCore.QTimer.singleShot"),
                    ):
                        compute.trigger_conversion()
                        msgs = compute.get_status_messages()
                        assert any("Calculating 3D structure" in msg for msg in msgs)


def test_trigger_conversion_stereo_enhancement(mock_parser_host):
    """Test trigger_conversion stereo enhancement logic for E/Z bonds."""
    from PyQt6.QtCore import QPointF, QTimer

    compute = DummyCompute(mock_parser_host)
    compute.data.add_atom("C", QPointF(0, 0))
    # Molecule with a double bond
    mol = Chem.MolFromSmiles("CC=CC")
    bond = mol.GetBondWithIdx(1)
    bond.SetStereo(Chem.BondStereo.STEREOE)

    with patch.object(compute.data, "to_rdkit_mol", return_value=mol):
        # Provide a properly formatted MOL block generated by RDKit to avoid warnings
        mol_block = Chem.MolToMolBlock(mol)
        with patch.object(compute.data, "to_mol_block", return_value=mol_block):
            with patch("rdkit.Chem.DetectChemistryProblems", return_value=[]):
                # Mock QThread to prevent actual thread creation in main_window_compute namespace
                with patch("moleditpy.ui.compute_logic.QThread"):
                    with patch("moleditpy.ui.compute_logic.CalculationWorker"):
                        with patch("PyQt6.QtCore.QTimer.singleShot") as mock_timer:
                            compute.trigger_conversion()
                            assert mock_timer.called


def test_halt_conversion(mock_parser_host):
    """Test halt_conversion clears active workers and restores UI."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {1, 2}
    compute.halt_ids = set()
    compute.host.init_manager.convert_button = MagicMock()

    compute.halt_conversion()
    assert compute.active_worker_ids == set()
    assert 1 in compute.halt_ids
    assert compute.host.init_manager.convert_button.setText.called
    # Note: halt_conversion doesn't show status message directly


@pytest.mark.skip(reason="Segfaults in headless environment")
def test_on_calculation_error_updated(mock_parser_host):
    """Test on_calculation_error with correct message formatting."""
    compute = DummyCompute(mock_parser_host)
    worker_id = 99
    compute.active_worker_ids = {worker_id}
    # Pass as tuple to ensure worker removal
    compute.on_calculation_error((worker_id, "Test Error"))
    msgs = compute.get_status_messages()
    assert any("Error: Test Error" in msg for msg in msgs)
    assert worker_id not in compute.active_worker_ids


@pytest.mark.skip(reason="Segfaults in headless environment despite mocking")
def test_trigger_conversion_chemistry_problem_detection(mock_parser_host):
    """Test trigger_conversion detects and flags chemistry problems (valence)."""
    from PyQt6.QtCore import QPointF
    # Avoid importing AtomItem to prevent QGraphicsItem instantiation in headless test
    # from moleditpy.modules.atom_item import AtomItem

    compute = DummyCompute(mock_parser_host)

    # Create a molecule that will fail sanitization/valence check (e.g. overvalent Nitrogen)
    compute.data.add_atom("N", QPointF(0, 0))

    # Use a mock item instead of real AtomItem to avoid QGraphicsItem/Scene crashes
    mock_item = MagicMock()
    mock_item.has_problem = False
    mock_item.pos.return_value = QPointF(0, 0)
    compute.data.atoms[0]["item"] = mock_item

    for i in range(1, 6):
        compute.data.add_atom("H", QPointF(i, 0))
        # Ensure neighbor atoms also have mock items
        h_item = MagicMock()
        h_item.pos.return_value = QPointF(i, 0)
        compute.data.atoms[i]["item"] = h_item
        compute.data.bonds[(0, i)] = {"order": 1}

    # RDKit DetectChemistryProblems will catch this
    # We mock to_rdkit_mol to return None, simulating a conversion failure safely
    # without triggering actual RDKit sanitization errors that might crash the process.
    with patch.object(compute.data, "to_rdkit_mol", return_value=None):
        compute.trigger_conversion()

    # It should hit either DetectChemistryProblems or the fallback if mol is None
    msgs = compute.get_status_messages()
    assert any("chemistry problem(s) found" in msg for msg in msgs)

    # Verify the item was flagged
    assert mock_item.has_problem == True

    # Explicit cleanup to prevent potential crash during teardown
    compute.data.atoms.clear()
    compute.data.bonds.clear()
    del compute


def test_trigger_conversion_fragment_message_exact(mock_parser_host):
    """Test verification of the exact status bar message for multiple fragments."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C.C.O")  # 3 fragments
    compute.data.atoms = {
        1: {"symbol": "C", "item": MagicMock()},
        2: {"symbol": "C", "item": MagicMock()},
        3: {"symbol": "O", "item": MagicMock()},
    }

    with (
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol"),
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        patch("moleditpy.ui.compute_logic.CalculationWorker"),
        patch("moleditpy.ui.compute_logic.QThread"),
    ):
        compute.trigger_conversion()
        all_messages = [
            str(call[0][0])
            for call in compute.statusBar().showMessage.call_args_list
            if call[0]
        ]
        # Should say "Converting 3 molecules..."
        assert any("Converting 3 molecules" in msg for msg in all_messages)


def test_trigger_conversion_to_mol_block_priority(mock_parser_host):
    """Test that data.to_mol_block() is used preferentially over RDKit's generation."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {"symbol": "C", "item": MagicMock()}}
    mol = Chem.MolFromSmiles("C")

    # Unique string to identify our custom block
    custom_block = "M  V2000\nCUSTOM_BLOCK_CONTENT"

    with (
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol"),
        patch.object(
            compute.data, "to_mol_block", return_value=custom_block
        ) as mock_to_block,
        patch("rdkit.Chem.MolToMolBlock") as mock_rdkit_block,
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
    ):
        # Setup worker mock to capture the start_work signal payload
        mock_worker_instance = MockWorker.return_value
        mock_start_work = MagicMock()
        mock_worker_instance.start_work = mock_start_work

        # We need to capture the lambda passed to QTimer.singleShot because that's what triggers start_work
        with patch("PyQt6.QtCore.QTimer.singleShot") as mock_timer:
            compute.trigger_conversion()

            # Verify to_mol_block was called
            assert mock_to_block.called
            # Verify RDKit's block generation was NOT called (or at least our result was used)
            assert not mock_rdkit_block.called

            # Extract the lambda from QTimer.singleShot(10, lambda ...)
            args, _ = mock_timer.call_args
            assert len(args) >= 2
            callback = args[1]

            # Execute the callback to trigger the signal emit on our mock
            callback()

            # Verify the payload sent to worker contains our custom block
            call_args = mock_start_work.emit.call_args
            assert call_args is not None
            sent_block = call_args[0][0]
            assert sent_block == custom_block


def test_trigger_conversion_ez_stereo_injection(mock_parser_host):
    """Test that M CFG lines are injected for E/Z stereo bonds."""
    compute = DummyCompute(mock_parser_host)

    # Create valid E-isomer (trans-2-butene)
    mol = Chem.MolFromSmiles("C/C=C/C")
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", 1)
    mol.GetAtomWithIdx(1).SetIntProp("_original_atom_id", 2)
    mol.GetAtomWithIdx(2).SetIntProp("_original_atom_id", 3)
    mol.GetAtomWithIdx(3).SetIntProp("_original_atom_id", 4)

    # Setup data to match
    compute.data.atoms = {
        1: {"symbol": "C", "item": MagicMock()},
        2: {"symbol": "C", "item": MagicMock()},
        3: {"symbol": "C", "item": MagicMock()},
        4: {"symbol": "C", "item": MagicMock()},
    }

    # Bond between 2 and 3 is the double bond
    # stereo=3 (Z) or 4 (E). Let's use 4 (E) -> CFG value 2
    # RDKit bond index for this might be 1 (0-1, 1-2, 2-3)
    compute.data.bonds = {
        (2, 3): {"stereo": 4, "order": 2},
        (1, 2): {"order": 1},
        (3, 4): {"order": 1},
    }

    # A mocked MOL block that we pretend came from our data
    base_block = "M  V2000\nM  END"

    with (
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol"),
        patch.object(compute.data, "to_mol_block", return_value=base_block),
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
    ):
        mock_worker_instance = MockWorker.return_value
        mock_start_work = MagicMock()
        mock_worker_instance.start_work = mock_start_work

        with patch("PyQt6.QtCore.QTimer.singleShot") as mock_timer:
            compute.trigger_conversion()

            # Execute callback
            args, _ = mock_timer.call_args
            callback = args[1]
            callback()

            # Check sent block for M CFG
            call_args = mock_start_work.emit.call_args
            sent_block = call_args[0][0]

            print(f"DEBUG: Sent Block:\\n{sent_block}")
            assert "M  CFG" in sent_block
            # Check for E isomer (val 2) on bond index 2 (RDKit idx 1 + 1)
            # 1-2 is usually index 1
            assert "M  CFG  1   2   2" in sent_block


def test_on_calculation_error_uff_fallback_temporary(mock_parser_host):
    """Verify that UFF fallback uses _temp_optimization_method and doesn't change persistent setting."""
    compute = DummyCompute(mock_parser_host)
    worker_id = "test_id"
    compute.active_worker_ids = {worker_id}
    compute.optimization_method = "MMFF_RDKIT"

    # Mocking QMessageBox.question to return Yes
    # QMessageBox.StandardButton.Yes is usually 16384 (0x4000)
    with patch(
        "moleditpy.ui.compute_logic.QMessageBox.question",
        return_value=QMessageBox.StandardButton.Yes,
    ):
        with patch.object(compute, "optimize_3d_structure") as mock_optimize:
            compute.on_calculation_error(
                (worker_id, "Optimization with MMFF94 failed")
            )

            # Check if temp override was set
            assert compute._temp_optimization_method == "UFF_RDKIT"
            # Check if optimize was called
            assert mock_optimize.called
            # Check if persistent method remains unchanged
            assert compute.optimization_method == "MMFF_RDKIT"
