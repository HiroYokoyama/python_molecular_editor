"""Unit tests for ComputeManager optimization trigger logic."""

from rdkit import Chem
from rdkit.Chem import AllChem
from moleditpy.ui.compute_logic import ComputeManager
from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QMessageBox
from unittest.mock import MagicMock, patch


class DummyCompute(ComputeManager):
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
        if not hasattr(host, "view_3d_manager"):
            host.view_3d_manager = MagicMock()
        if not hasattr(host, "state_manager"):
            host.state_manager = MagicMock()
        if not hasattr(host, "ui_manager"):
            host.ui_manager = MagicMock()

        # Ensure buttons/actions exist for UI transition tests
        for btn in [
            "convert_button",
            "cleanup_button",
            "optimize_3d_button",
            "export_button",
            "analysis_action",
            "edit_3d_action",
        ]:
            if not hasattr(host.init_manager, btn):
                setattr(host.init_manager, btn, MagicMock())

    def __getattr__(self, name):
        return getattr(self._host, name)

    def _start_calculation_worker(self, mol_block, options, run_id) -> None:
        # Call the real _start_calculation_worker but with QThread.start mocked
        # to prevent background threads from starting, avoiding random segfaults/aborts.
        from unittest.mock import patch

        with patch("PyQt6.QtCore.QThread.start"):
            ComputeManager._start_calculation_worker(self, mol_block, options, run_id)

    @property
    def data(self):
        return self.host.state_manager.data

    @data.setter
    def data(self, v):
        self.host.state_manager.data = v

    @property
    def scene(self):
        return self.host.init_manager.scene

    @scene.setter
    def scene(self, v):
        self.host.init_manager.scene = v

    @property
    def view_2d(self):
        return self.host.init_manager.view_2d

    @view_2d.setter
    def view_2d(self, v):
        self.host.init_manager.view_2d = v

    @property
    def view_3d(self):
        return self.host.view_3d_manager.view_3d

    @view_3d.setter
    def view_3d(self, v):
        self.host.view_3d_manager.view_3d = v

    @property
    def plotter(self):
        return self.host.view_3d_manager.plotter

    @plotter.setter
    def plotter(self, v):
        self.host.view_3d_manager.plotter = v

    @property
    def settings(self):
        return self.host.init_manager.settings

    @settings.setter
    def settings(self, v):
        self.host.init_manager.settings = v

    @property
    def current_mol(self):
        return self.host.view_3d_manager.current_mol

    @current_mol.setter
    def current_mol(self, v):
        self.host.view_3d_manager.current_mol = v

    @property
    def optimization_method(self):
        return self.host.init_manager.settings.get("optimization_method", "MMFF_RDKIT")

    @optimization_method.setter
    def optimization_method(self, v):
        self.host.init_manager.settings["optimization_method"] = v
        self.host.init_manager.optimization_method = v

    @property
    def constraints_3d(self):
        return self.host.edit_3d_manager.constraints_3d

    @constraints_3d.setter
    def constraints_3d(self, v):
        self.host.edit_3d_manager.constraints_3d = v

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
    # Stale workers still show a status message (informational, not an error)
    compute.statusBar().showMessage.assert_called()
    msg = compute.statusBar().showMessage.call_args[0][0]
    assert "stale" in msg.lower() or "Ignored" in msg


def test_on_calculation_error_basic(mock_parser_host):
    """Test on_calculation_error for an ACTIVE worker."""
    compute = DummyCompute(mock_parser_host)
    worker_id = "active_id"
    compute.active_worker_ids = {worker_id}
    compute.on_calculation_error((worker_id, "Real Error"))
    assert compute.statusBar().showMessage.called
    assert "Real Error" in compute.statusBar().showMessage.call_args[0][0]


def test_on_calculation_error_reenables_export(mock_parser_host):
    """A failed optimization runs the full UI restore (which re-enables Export),
    not just the convert/optimize button restore."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"active_id"}
    with patch.object(compute, "_refresh_ui_state") as refresh:
        compute.on_calculation_error(("active_id", "Optimization failed"))
    refresh.assert_called_once()


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
    with patch.object(compute, "draw_molecule_3d"):
        compute.on_calculation_finished(result)
        assert compute.host.view_3d_manager.current_mol == mol
        assert worker_id not in compute.active_worker_ids


def test_check_chemistry_problems_fallback_detects(mock_parser_host):
    """Verify that the manual valence fallback correctly identifies overvalent atoms."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {}
    compute.data.bonds = {}
    c_id = mock_parser_host.scene.create_atom("C", QPointF(0, 0))
    c_item = mock_parser_host.scene.atom_items[c_id]
    for i in range(5):
        h_id = mock_parser_host.scene.create_atom("H", QPointF(i + 1, 0))
        h_item = mock_parser_host.scene.atom_items[h_id]
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
    compute.data.atoms = {1: {"symbol": "C"}}
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
    compute.data.atoms = {1: {"symbol": "C"}}

    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0
    with (
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[problem]),
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        # Patch the modal dialog: the MagicMock host is not a valid QWidget
        # parent, so the real QMessageBox.critical would raise TypeError.
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
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
    compute.data.atoms = {1: {"symbol": "C"}}

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
        1: {"symbol": "C"},
        2: {"symbol": "C"},
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

    with (
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("moleditpy.ui.compute_logic.QThread"),
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        MockWorker.return_value
        compute.optimize_3d_structure("MMFF_RDKIT")

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
    compute.host.plugin_manager.optimization_methods = {
        "CUSTOM": {"callback": MagicMock()}
    }
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


def test_on_calculation_error_mmff_fail_shows_uff_dialog(mock_parser_host):
    """When MMFF fails, on_calculation_error shows a UFF retry dialog instead of auto-fallback."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}

    with patch(
        "moleditpy.ui.compute_logic.QMessageBox.question",
        return_value=QMessageBox.StandardButton.No,
    ) as mock_dialog:
        compute.on_calculation_error(
            ("w1", "Optimization with MMFF94s (RDKit) failed.")
        )

    mock_dialog.assert_called_once()
    call_args = mock_dialog.call_args[0]
    assert "UFF" in call_args[2]


def test_on_calculation_error_mmff_fail_user_accepts_uff(mock_parser_host):
    """When user accepts UFF retry after MMFF failure, optimize_3d_structure is called with UFF_RDKIT."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}

    with (
        patch(
            "moleditpy.ui.compute_logic.QMessageBox.question",
            return_value=QMessageBox.StandardButton.Yes,
        ),
        patch.object(compute, "optimize_3d_structure") as mock_optimize,
    ):
        compute.on_calculation_error(
            ("w1", "Optimization with MMFF94s (RDKit) failed.")
        )

    mock_optimize.assert_called_once_with("UFF_RDKIT")


def test_on_calculation_error_mmff_fail_user_declines_uff(mock_parser_host):
    """When user declines UFF retry, a critical error dialog is shown instead."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}

    with (
        patch(
            "moleditpy.ui.compute_logic.QMessageBox.question",
            return_value=QMessageBox.StandardButton.No,
        ),
        patch("moleditpy.ui.compute_logic.QMessageBox.critical") as mock_critical,
    ):
        compute.on_calculation_error(
            ("w1", "Optimization with MMFF94s (RDKit) failed.")
        )

    mock_critical.assert_called_once()


def test_on_calculation_error_non_mmff_no_uff_dialog(mock_parser_host):
    """Non-MMFF errors do not trigger the UFF retry dialog."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}

    with patch("moleditpy.ui.compute_logic.QMessageBox.question") as mock_dialog:
        compute.on_calculation_error(("w1", "Some other error occurred."))

    mock_dialog.assert_not_called()


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
    from moleditpy.ui.app_state import StateManager
    from moleditpy.core.molecular_data import MolecularData
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)
    compute.data = MolecularData()  # Use real MolecularData for this test

    # Setup initial state with radical and constraint
    pos = QPointF(10, 20)
    aid = compute.data.add_atom("C", pos, radical=2)

    compute.constraints_3d = [("DISTANCE", (0, 1), 1.5, 1e5)]

    # StateManager is the canonical class; call its methods by passing 'compute' as self
    state = StateManager.get_current_state(compute)

    # Verify state content
    assert state["atoms"][0]["radical"] == 2
    # The code converts tuple to list for JSON compatibility in state
    assert state["constraints_3d"][0] == ["DISTANCE", [0, 1], 1.5, 1e5]

    # Clear and Restore
    compute.data.atoms.clear()
    compute.constraints_3d = []
    compute.scene = MagicMock()

    def mock_restore(raw_atoms, raw_bonds):
        for atom_id, data in raw_atoms.items():
            compute.data.atoms[atom_id] = {
                "symbol": data["symbol"],
                "pos": tuple(data["pos"]),
                "charge": data.get("charge", 0),
                "radical": data.get("radical", 0),
            }
        for key, data in raw_bonds.items():
            compute.data.bonds[key] = {
                "order": data.get("order", 1),
                "stereo": data.get("stereo", 0),
            }

    compute.scene.restore_atoms_and_bonds.side_effect = mock_restore

    StateManager.set_state_from_data(compute, state)

    # Verify restoration
    assert 0 in compute.data.atoms
    assert compute.data.atoms[0]["radical"] == 2
    # set_state_from_data converts back to tuple
    assert compute.constraints_3d == [("DISTANCE", (0, 1), 1.5, 1e5)]


def test_app_state_original_atom_id_preservation(mock_parser_host):
    """Test that _original_atom_id is preserved in 3D molecule state round-trip."""
    from moleditpy.ui.app_state import StateManager

    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    AllChem.EmbedMolecule(mol)
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", 123)
    compute.current_mol = mol

    # Get state
    state = StateManager.get_current_state(compute)
    assert 123 in state["mol_3d_atom_ids"]

    # Clear and Restore
    compute.current_mol = None
    StateManager.set_state_from_data(compute, state)

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


def test_check_chemistry_problems_fallback(mock_parser_host):
    """Test the manual valence check when RDKit fails."""
    # from moleditpy.modules.atom_item import AtomItem
    from PyQt6.QtCore import QPointF

    compute = DummyCompute(mock_parser_host)

    # Setup a problematic atom: Carbon with total bond order 5
    pos = QPointF(0, 0)
    compute.data.add_atom("C", pos)
    for i in range(4):
        compute.data.add_atom("C", QPointF(50.0 * (i + 1), 0))

    # Use Mock item to avoid Qt crash
    mock_item = MagicMock()
    mock_item.has_problem = False
    mock_item.pos.return_value = pos
    compute.host.init_manager.scene.atom_items[0] = mock_item

    # Bonds sum to order 5 on atom 0 (2+1+1+1)
    compute.data.bonds = {
        (0, 1): {"order": 2},
        (0, 2): {"order": 1},
        (0, 3): {"order": 1},
        (0, 4): {"order": 1},
    }

    compute.check_chemistry_problems_fallback()
    assert mock_item.has_problem == True
    msgs = compute.get_status_messages()
    assert any("chemistry problems found" in msg for msg in msgs)


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
    from PyQt6.QtCore import QPointF

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
    compute.host.init_manager.scene.atom_items[0] = mock_item

    for i in range(1, 6):
        compute.data.add_atom("H", QPointF(i, 0))
        compute.data.bonds[(0, i)] = {"order": 1}

    # RDKit DetectChemistryProblems will catch this
    # We mock to_rdkit_mol to return None, simulating a conversion failure safely
    # without triggering actual RDKit sanitization errors that might crash the process.
    with patch.object(compute.data, "to_rdkit_mol", return_value=None):
        compute.trigger_conversion()

    # It should hit either DetectChemistryProblems or the fallback if mol is None
    msgs = compute.get_status_messages()
    assert any("chemistry problems found" in msg for msg in msgs)

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
        1: {"symbol": "C"},
        2: {"symbol": "C"},
        3: {"symbol": "O"},
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
    compute.data.atoms = {1: {"symbol": "C"}}
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


def test_trigger_conversion_uses_default_when_temp_mode_is_none(mock_parser_host):
    """Preinitialized temp conversion mode should not suppress the saved setting."""
    compute = DummyCompute(mock_parser_host)
    compute.data.atoms = {1: {"symbol": "C"}}
    compute.settings["3d_conversion_mode"] = "fallback"
    mol = Chem.MolFromSmiles("C")

    with (
        patch.object(compute.data, "to_rdkit_mol", return_value=mol),
        patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
        patch("rdkit.Chem.SanitizeMol"),
        patch.object(
            compute.data, "to_mol_block", return_value=Chem.MolToMolBlock(mol)
        ),
        patch.object(compute, "_start_calculation_worker") as mock_start,
    ):
        compute.trigger_conversion(None)

    options = mock_start.call_args[0][1]
    assert options["conversion_mode"] == "fallback"


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
        1: {"symbol": "C"},
        2: {"symbol": "C"},
        3: {"symbol": "C"},
        4: {"symbol": "C"},
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
            compute.on_calculation_error((worker_id, "Optimization with MMFF94 failed"))

            # Check if optimize was called with UFF_RDKIT
            mock_optimize.assert_called_once_with("UFF_RDKIT")
            # Check if persistent method remains unchanged
            assert compute.optimization_method == "MMFF_RDKIT"


# =============================================================================
# Plugin-registered optimization methods (register_optimization_method)
# =============================================================================


def _mol_with_conformer():
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol)
    return mol


def _compute_with_plugin_method(mock_parser_host, callback):
    compute = DummyCompute(mock_parser_host)
    mock_parser_host.view_3d_manager.current_mol = _mol_with_conformer()
    mock_parser_host.plugin_manager.optimization_methods = {
        "MYOPT": {"plugin": "P", "callback": callback, "label": "My Optimizer"}
    }
    return compute


def test_plugin_optimization_method_invoked(mock_parser_host):
    """A plugin-registered method is dispatched with the current mol."""
    calls = []
    compute = _compute_with_plugin_method(
        mock_parser_host, lambda mol: calls.append(mol) or True
    )

    compute.optimize_3d_structure("MYOPT")

    assert calls == [mock_parser_host.view_3d_manager.current_mol]
    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_called_once()
    mock_parser_host.edit_actions_manager.push_undo_state.assert_called_once()
    assert compute.last_successful_optimization_method == "My Optimizer"


def test_plugin_optimization_failure_is_isolated(mock_parser_host):
    """A raising plugin callback is caught, logged, and reported via status."""

    def _boom(mol):
        raise RuntimeError("plugin exploded")

    compute = _compute_with_plugin_method(mock_parser_host, _boom)
    compute.optimize_3d_structure("MYOPT")  # must not raise

    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_not_called()
    messages = [
        str(c[0][0])
        for c in mock_parser_host.update_status_message.call_args_list
        if c[0]
    ]
    assert any("failed" in m for m in messages)


def test_plugin_optimization_false_return_reports_failure(mock_parser_host):
    """A callback returning False reports failure without redrawing."""
    compute = _compute_with_plugin_method(mock_parser_host, lambda mol: False)
    compute.optimize_3d_structure("MYOPT")

    mock_parser_host.view_3d_manager.draw_molecule_3d.assert_not_called()
    assert compute.last_successful_optimization_method is None


# =============================================================================
# Additional coverage: small helpers, menus, worker plumbing, error branches
# =============================================================================


def test_reset_active_threads(mock_parser_host):
    """reset_active_threads empties the tracked thread list."""
    compute = DummyCompute(mock_parser_host)
    compute._active_calc_threads = [MagicMock(), MagicMock()]
    compute.reset_active_threads()
    assert compute._active_calc_threads == []


def test_safe_disconnect_swallows_runtime_error(mock_parser_host):
    """_safe_disconnect must not propagate a RuntimeError from disconnect()."""
    compute = DummyCompute(mock_parser_host)
    signal = MagicMock()
    signal.disconnect.side_effect = RuntimeError("not connected")
    compute._safe_disconnect(signal)  # must not raise
    signal.disconnect.assert_called_once()


def test_remove_calculating_text_removes_actor(mock_parser_host):
    """_remove_calculating_text removes the actor from the plotter renderer."""
    compute = DummyCompute(mock_parser_host)
    actor = MagicMock()
    compute._calculating_text_actor = actor
    plotter = compute.host.view_3d_manager.plotter
    compute._remove_calculating_text()
    plotter.renderer.RemoveActor.assert_called_once_with(actor)
    assert compute._calculating_text_actor is None


def test_set_optimization_method_ignores_empty(mock_parser_host):
    """An empty method name is a no-op (no status message, no host update)."""
    compute = DummyCompute(mock_parser_host)
    compute.host.set_optimization_method = MagicMock()
    compute.set_optimization_method("")
    compute.host.set_optimization_method.assert_not_called()


def test_set_optimization_method_checks_matching_action(mock_parser_host):
    """The action matching the chosen method gets checked, others unchecked."""
    compute = DummyCompute(mock_parser_host)
    act_mmff = MagicMock()
    act_uff = MagicMock()
    compute.host.init_manager.opt3d_actions = {
        "MMFF_RDKIT": act_mmff,
        "UFF_RDKIT": act_uff,
    }
    compute.set_optimization_method("mmff_rdkit")
    act_mmff.setChecked.assert_called_with(True)
    act_uff.setChecked.assert_called_with(False)


def test_toggle_intermolecular_interaction_rdkit(mock_parser_host):
    """Toggling intermolecular interaction stores the flag and reports status."""
    compute = DummyCompute(mock_parser_host)
    compute.toggle_intermolecular_interaction_rdkit(True)
    assert (
        compute.host.get_settings()["optimize_intermolecular_interaction_rdkit"] is True
    )
    assert any("Enabled" in m for m in compute.get_status_messages())


def test_show_convert_menu_disabled_button_returns(mock_parser_host):
    """show_convert_menu bails out when the convert button is disabled."""
    compute = DummyCompute(mock_parser_host)
    compute.host.init_manager.convert_button.isEnabled.return_value = False
    with patch("moleditpy.ui.compute_logic.QMenu") as MockMenu:
        compute.show_convert_menu(QPointF(0, 0).toPoint())
    MockMenu.assert_not_called()


def test_show_convert_menu_builds_actions(mock_parser_host):
    """show_convert_menu builds one action per conversion option and execs the menu."""
    compute = DummyCompute(mock_parser_host)
    compute.host.init_manager.convert_button.isEnabled.return_value = True
    with (
        patch("moleditpy.ui.compute_logic.QMenu") as MockMenu,
        patch("moleditpy.ui.compute_logic.QAction") as MockAction,
    ):
        menu = MockMenu.return_value
        compute.show_convert_menu(QPointF(1, 2).toPoint())
    assert MockAction.call_count == 4  # fallback, rdkit, obabel, direct
    menu.exec.assert_called_once()


def test_show_optimize_menu_builds_actions(mock_parser_host):
    """show_optimize_menu mirrors the registered opt3d actions into a temporary menu."""
    compute = DummyCompute(mock_parser_host)
    compute.host.init_manager.optimize_3d_button.isEnabled.return_value = True
    src = MagicMock()
    src.text.return_value = "&MMFF94s"
    src.isEnabled.return_value = True
    compute.host.init_manager.opt3d_actions = {"MMFF_RDKIT": src}
    with (
        patch("moleditpy.ui.compute_logic.QMenu") as MockMenu,
        patch("moleditpy.ui.compute_logic.QAction") as MockAction,
    ):
        menu = MockMenu.return_value
        compute.show_optimize_menu(QPointF(0, 0).toPoint())
    MockAction.assert_called_once()
    menu.exec.assert_called_once()


def test_show_optimize_menu_disabled_returns(mock_parser_host):
    """show_optimize_menu bails out when the optimize button is disabled."""
    compute = DummyCompute(mock_parser_host)
    compute.host.init_manager.optimize_3d_button.isEnabled.return_value = False
    with patch("moleditpy.ui.compute_logic.QMenu") as MockMenu:
        compute.show_optimize_menu(QPointF(0, 0).toPoint())
    MockMenu.assert_not_called()


def test_trigger_optimize_with_temp_method_schedules_call(mock_parser_host):
    """_trigger_optimize_with_temp_method schedules optimize on the event loop."""
    compute = DummyCompute(mock_parser_host)
    with patch("moleditpy.ui.compute_logic.QTimer.singleShot") as mock_timer:
        compute._trigger_optimize_with_temp_method("UFF_RDKIT")
    mock_timer.assert_called_once()


def test_optimize_3d_structure_no_current_mol(mock_parser_host):
    """optimize_3d_structure reports when there is no 3D molecule at all."""
    compute = DummyCompute(mock_parser_host)
    compute.host.view_3d_manager.current_mol = None
    compute.optimize_3d_structure()
    assert any("No 3D molecule to optimize" in m for m in compute.get_status_messages())


def test_handle_chemistry_problems_flags_matched_item(mock_parser_host):
    """_handle_chemistry_problems flags the 2D scene item for a problem atom."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", 7)
    prob = MagicMock()
    prob.GetAtomIdx.return_value = 0
    item = MagicMock()
    item.has_problem = False
    compute.host.init_manager.scene.atom_items = {7: item}
    with patch("moleditpy.ui.compute_logic.QMessageBox.critical"):
        ComputeManager._handle_chemistry_problems(compute, mol, [prob])
    assert item.has_problem is True
    item.update.assert_called_once()


def test_setup_mol_block_falls_back_to_rdkit(mock_parser_host):
    """_setup_mol_block_for_worker uses RDKit generation when data yields no block."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.Compute2DCoords(mol)
    with patch.object(compute.data, "to_mol_block", return_value=""):
        block = ComputeManager._setup_mol_block_for_worker(compute, mol)
    assert isinstance(block, str) and "V2000" in block


def test_start_calculation_worker_wires_callbacks(mock_parser_host):
    """The real worker plumbing wires finished/error callbacks and cleans up."""
    compute = DummyCompute(mock_parser_host)
    with (
        patch("moleditpy.ui.compute_logic.QThread") as MockThread,
        patch("moleditpy.ui.compute_logic.CalculationWorker") as MockWorker,
        patch("PyQt6.QtCore.QTimer.singleShot"),
    ):
        thread = MockThread.return_value
        worker = MockWorker.return_value
        ComputeManager._start_calculation_worker(compute, "block", {"worker_id": 1}, 1)
        on_finished = worker.finished.connect.call_args[0][0]
        on_error = worker.error.connect.call_args[0][0]

    assert thread in compute._active_calc_threads

    with patch.object(compute, "on_calculation_finished") as fin:
        on_finished("result")
        fin.assert_called_once_with("result")
    # _cleanup removed the thread and quit it
    assert thread not in compute._active_calc_threads
    thread.quit.assert_called()

    # Second callback path (error) — cleanup remove now raises ValueError, suppressed
    with patch.object(compute, "on_calculation_error") as err:
        on_error("boom")
        err.assert_called_once_with("boom")


def test_on_calculation_finished_restores_atom_props(mock_parser_host):
    """Original atom ids stored before conversion are re-applied to the result mol."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    compute.original_atom_properties = {0: 42}
    compute.on_calculation_finished(mol)
    assert mol.GetAtomWithIdx(0).GetIntProp("_original_atom_id") == 42


def test_on_calculation_finished_method_lookup_error_is_safe(mock_parser_host):
    """A mol whose property access raises is handled without propagating."""
    compute = DummyCompute(mock_parser_host)
    mol = MagicMock()
    mol.HasProp.side_effect = TypeError("bad prop")
    mol.GetNumAtoms.return_value = 0
    compute.on_calculation_finished(mol)  # must not raise
    # Falls through the except; nothing crashes.
    assert compute.host.view_3d_manager.current_mol is mol


def test_on_calculation_error_active_halt_shows_halted(mock_parser_host):
    """An active worker signalling 'Halt' produces a plain 'Halted' status."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}
    compute.on_calculation_error(("w1", "Halt"))
    assert any(m == "Halted" for m in compute.get_status_messages())


def test_on_calculation_error_stale_halt_shows_id(mock_parser_host):
    """A stale worker signalling 'Halt' reports it was ignored with its id."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"active"}
    compute.on_calculation_error(("stale", "Halted"))
    assert any("Ignored halted worker" in m for m in compute.get_status_messages())


def test_on_calculation_error_uff_dialog_exception_is_caught(mock_parser_host):
    """A TypeError from the UFF retry dialog is swallowed and a critical shown."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}
    with (
        patch(
            "moleditpy.ui.compute_logic.QMessageBox.question",
            side_effect=TypeError("bad parent"),
        ),
        patch("moleditpy.ui.compute_logic.QMessageBox.critical") as mock_critical,
    ):
        compute.on_calculation_error(("w1", "MMFF setup failed"))
    mock_critical.assert_called_once()


def test_create_atom_id_mapping_no_mol_returns(mock_parser_host):
    """create_atom_id_mapping is a no-op when there is no current molecule."""
    compute = DummyCompute(mock_parser_host)
    compute.host.view_3d_manager.current_mol = None
    compute.host.atom_id_to_rdkit_idx_map = {"stale": 1}
    ComputeManager.create_atom_id_mapping(compute)
    # Early return leaves the map untouched (not reset to {}).
    assert compute.host.atom_id_to_rdkit_idx_map == {"stale": 1}


def test_create_atom_id_mapping_builds_map(mock_parser_host):
    """create_atom_id_mapping maps original atom ids to RDKit indices."""
    compute = DummyCompute(mock_parser_host)
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    mol.GetAtomWithIdx(0).SetIntProp("_original_atom_id", 10)
    mol.GetAtomWithIdx(1).SetIntProp("_original_atom_id", 20)
    compute.host.view_3d_manager.current_mol = mol
    compute.host.atom_id_to_rdkit_idx_map = {}
    ComputeManager.create_atom_id_mapping(compute)
    assert compute.host.atom_id_to_rdkit_idx_map[10] == 0
    assert compute.host.atom_id_to_rdkit_idx_map[20] == 1


def test_update_aromatic_rings_sets_rings(mock_parser_host):
    """update_aromatic_rings detects a benzene ring and forwards it to the scene."""
    compute = DummyCompute(mock_parser_host)
    benzene = Chem.MolFromSmiles("c1ccccc1")
    with patch.object(compute.data, "to_rdkit_mol", return_value=benzene):
        compute.update_aromatic_rings()
    compute.host.init_manager.scene.set_aromatic_rings.assert_called_once()
    rings = compute.host.init_manager.scene.set_aromatic_rings.call_args[0][0]
    assert len(rings) == 1 and len(rings[0]) == 6


def test_update_aromatic_rings_none_mol_is_noop(mock_parser_host):
    """update_aromatic_rings returns quietly when no mol can be built."""
    compute = DummyCompute(mock_parser_host)
    with patch.object(compute.data, "to_rdkit_mol", return_value=None):
        compute.update_aromatic_rings()
    compute.host.init_manager.scene.set_aromatic_rings.assert_not_called()


def test_select_connected_atoms_expands_selection(mock_parser_host):
    """select_connected_atoms grows the selection over the bond graph."""
    compute = DummyCompute(mock_parser_host)
    a0 = compute.data.add_atom("C", QPointF(0, 0))
    a1 = compute.data.add_atom("C", QPointF(1, 0))
    a2 = compute.data.add_atom("C", QPointF(2, 0))
    compute.data.add_bond(a0, a1, order=1)
    compute.data.add_bond(a1, a2, order=1)

    scene = compute.host.init_manager.scene
    sel_item = MagicMock()
    sel_item.atom_id = a0
    scene.selectedItems.return_value = [sel_item]
    items = {aid: MagicMock() for aid in (a0, a1, a2)}
    scene.atom_items = items

    compute.select_connected_atoms()
    for aid in (a0, a1, a2):
        items[aid].setSelected.assert_called_with(True)


def test_select_connected_atoms_no_selection_returns(mock_parser_host):
    """select_connected_atoms is a no-op when nothing is selected."""
    compute = DummyCompute(mock_parser_host)
    compute.host.init_manager.scene.selectedItems.return_value = []
    compute.select_connected_atoms()  # must not raise


# ---------------------------------------------------------------------------
# Focus-restore tests (2D editor regains focus after conversion errors)
# ---------------------------------------------------------------------------


def test_handle_chemistry_problems_defers_set_focus(mock_parser_host):
    """_handle_chemistry_problems must schedule a deferred setFocus on view_2d
    via QTimer.singleShot so that Qt finishes routing the modal-dialog result
    before focus is stolen back from the Convert button.
    """
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    atom = mol.GetAtomWithIdx(0)
    atom.SetIntProp("_original_atom_id", 1)

    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0

    singleshot_calls = []

    def capture_singleshot(delay, callback):
        singleshot_calls.append((delay, callback))

    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot", side_effect=capture_singleshot),
    ):
        compute._handle_chemistry_problems(mol, [problem])

    # A deferred call must have been scheduled targeting view_2d.setFocus
    view_2d = compute.host.init_manager.view_2d
    assert any(
        cb == view_2d.setFocus for _delay, cb in singleshot_calls
    ), "setFocus on view_2d was not scheduled via QTimer.singleShot"


def test_handle_chemistry_problems_sets_focus_immediately_absent_timer(mock_parser_host):
    """When QTimer.singleShot fires synchronously (delay=0 executed right away),
    view_2d.setFocus is ultimately called — verify the callback is the right target.
    """
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    atom = mol.GetAtomWithIdx(0)
    atom.SetIntProp("_original_atom_id", 1)

    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0

    # Execute the callback immediately to simulate the timer firing
    def immediate_singleshot(delay, callback):
        callback()

    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot", side_effect=immediate_singleshot),
    ):
        compute._handle_chemistry_problems(mol, [problem])

    compute.host.init_manager.view_2d.setFocus.assert_called()


def test_on_calculation_error_defers_set_focus_to_view_2d(mock_parser_host):
    """on_calculation_error must schedule a deferred setFocus on view_2d after
    showing the error dialog, so keyboard shortcuts in the 2D editor are restored.
    """
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w1"}

    singleshot_calls = []

    def capture_singleshot(delay, callback):
        singleshot_calls.append((delay, callback))

    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot", side_effect=capture_singleshot),
    ):
        compute.on_calculation_error(("w1", "RDKit 3D conversion failed"))

    view_2d = compute.host.init_manager.view_2d
    assert any(
        cb == view_2d.setFocus for _delay, cb in singleshot_calls
    ), "setFocus on view_2d was not scheduled via QTimer.singleShot in on_calculation_error"


def test_on_calculation_error_focus_fires_when_timer_executes(mock_parser_host):
    """When the deferred QTimer callback fires, view_2d.setFocus is called."""
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w2"}

    def immediate_singleshot(delay, callback):
        callback()

    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot", side_effect=immediate_singleshot),
    ):
        compute.on_calculation_error(("w2", "Embed failed"))

    compute.host.init_manager.view_2d.setFocus.assert_called()


def test_on_calculation_error_halt_restores_focus(mock_parser_host):
    """A 'Halt' error still schedules a deferred setFocus so the 2D editor
    remains keyboard-accessible after the user halts a conversion.
    """
    compute = DummyCompute(mock_parser_host)
    compute.active_worker_ids = {"w3"}

    singleshot_calls = []

    def capture_singleshot(delay, callback):
        singleshot_calls.append((delay, callback))

    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot", side_effect=capture_singleshot),
    ):
        compute.on_calculation_error(("w3", "Halt"))

    view_2d = compute.host.init_manager.view_2d
    assert any(
        cb == view_2d.setFocus for _delay, cb in singleshot_calls
    ), "setFocus on view_2d must be scheduled even after a Halt"


def test_chemistry_problems_focus_not_called_directly(mock_parser_host):
    """_handle_chemistry_problems must NOT call view_2d.setFocus() directly
    (synchronously); only the deferred QTimer path is allowed.
    This guards against regression to the pre-fix synchronous call.
    """
    compute = DummyCompute(mock_parser_host)
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    atom = mol.GetAtomWithIdx(0)
    atom.SetIntProp("_original_atom_id", 1)

    problem = MagicMock()
    problem.GetAtomIdx.return_value = 0

    # Do NOT execute the timer callback — capture it without firing
    with (
        patch("moleditpy.ui.compute_logic.QMessageBox.critical"),
        patch("moleditpy.ui.compute_logic.QTimer.singleShot"),  # no-op — callback never fires
    ):
        compute._handle_chemistry_problems(mol, [problem])

    # If setFocus was called synchronously it would appear here; it must not.
    compute.host.init_manager.view_2d.setFocus.assert_not_called()

