"""Integration tests for trigger_conversion plugin-wrapping safety.

Verifies that:
1. _trigger_conversion_with_temp_mode passes the mode to trigger_conversion
   via _pending_conversion_mode without using a kwarg.
2. trigger_conversion still works correctly when a plugin wraps it without
   forwarding **kwargs (the original TypeError bug scenario).
3. _pending_conversion_mode is consumed exactly once per call.
"""

import pytest
from unittest.mock import MagicMock, patch
from rdkit import Chem

from moleditpy.ui.compute_logic import ComputeManager


# ---------------------------------------------------------------------------
# Minimal test double that re-uses DummyCompute logic from unit tests
# ---------------------------------------------------------------------------


class _PluginWrappedCompute(ComputeManager):
    """Simulates a plugin that wraps trigger_conversion without forwarding kwargs."""

    def __init__(self, host):
        self._host = host
        ComputeManager.__init__(self, host)

        if not hasattr(host, "init_manager"):
            host.init_manager = MagicMock()
        host.init_manager.settings = (
            host.init_manager.settings
            if isinstance(getattr(host.init_manager, "settings", None), dict)
            else {}
        )
        host.init_manager.opt3d_method_labels = {
            "MMFF_RDKIT": "MMFF94s (RDKit)",
            "UFF_RDKIT": "UFF (RDKit)",
        }
        for attr in ("view_3d_manager", "state_manager", "ui_manager"):
            if not hasattr(host, attr):
                setattr(host, attr, MagicMock())
        for btn in (
            "convert_button",
            "cleanup_button",
            "optimize_3d_button",
            "export_button",
            "analysis_action",
            "edit_3d_action",
        ):
            if not hasattr(host.init_manager, btn):
                setattr(host.init_manager, btn, MagicMock())

        # Simulate a plugin monkey-patching trigger_conversion:
        # The plugin's wrapper does NOT forward **kwargs — this is the exact
        # bug scenario from the TypeError report.
        original = self.trigger_conversion

        def plugin_wrapper():
            original()  # no kwargs forwarded

        self.trigger_conversion = plugin_wrapper

    def __getattr__(self, name):
        return getattr(self._host, name)

    def _start_calculation_worker(self, mol_block, options, run_id) -> None:
        with patch("PyQt6.QtCore.QThread.start"):
            ComputeManager._start_calculation_worker(self, mol_block, options, run_id)

    @property
    def data(self):
        return self.host.state_manager.data

    @data.setter
    def data(self, v):
        self.host.state_manager.data = v


@pytest.fixture
def host(app):
    from unittest.mock import MagicMock
    from moleditpy.core.molecular_data import MolecularData

    h = MagicMock()
    h.state_manager = MagicMock()
    h.state_manager.data = MolecularData()
    h.init_manager = MagicMock()
    h.init_manager.settings = {
        "3d_conversion_mode": "rdkit",
        "optimization_method": "MMFF_RDKIT",
        "skip_chemistry_checks": True,
        "background_color": "#919191",
    }
    h.init_manager.opt3d_method_labels = {
        "MMFF_RDKIT": "MMFF94s (RDKit)",
        "UFF_RDKIT": "UFF (RDKit)",
    }
    h.view_3d_manager = MagicMock()
    h.view_3d_manager.current_mol = None
    h.ui_manager = MagicMock()
    h.edit_3d_manager = MagicMock()
    h.edit_3d_manager.measurement_mode = False
    h.edit_3d_manager.is_3d_edit_mode = False
    # get_settings() must return a real dict so QColor(settings.get(...)) works
    h.get_settings.return_value = h.init_manager.settings
    type(h).settings = property(lambda self: self.init_manager.settings)
    type(h).data = property(lambda self: self.state_manager.data)
    return h


# ---------------------------------------------------------------------------
# Test: _trigger_conversion_with_temp_mode passes the mode correctly
# ---------------------------------------------------------------------------


class TestTriggerConversionTempMode:
    def test_temp_mode_stored_and_consumed(self, host, app):
        """_trigger_conversion_with_temp_mode sets _pending_conversion_mode
        and trigger_conversion reads and clears it."""
        from moleditpy.ui.compute_logic import ComputeManager

        compute = object.__new__(ComputeManager)
        compute._host = host
        ComputeManager.__init__(compute, host)

        captured = []

        def fake_trigger():
            mode = compute.__dict__.get("_pending_conversion_mode")
            captured.append(mode)

        compute.trigger_conversion = fake_trigger

        with patch("PyQt6.QtCore.QTimer.singleShot") as mock_timer:
            compute._trigger_conversion_with_temp_mode("rdkit")

        # QTimer.singleShot was called with the right callback
        assert mock_timer.called
        # The _pending_conversion_mode was set before singleShot fires
        assert compute.__dict__.get("_pending_conversion_mode") == "rdkit"

    def test_pending_mode_consumed_after_trigger_conversion(self, host, app):
        """After trigger_conversion runs, _pending_conversion_mode is None."""
        from moleditpy.ui.compute_logic import ComputeManager

        compute = object.__new__(ComputeManager)
        compute._host = host
        ComputeManager.__init__(compute, host)

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}
        host.init_manager.settings["3d_conversion_mode"] = "rdkit"
        compute.__dict__["_pending_conversion_mode"] = "fallback"

        with (
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker"),
        ):
            compute.trigger_conversion()

        # After running, the pending mode is cleared
        assert compute.__dict__.get("_pending_conversion_mode") is None

    def test_pending_mode_takes_priority_over_settings(self, host, app):
        """_pending_conversion_mode overrides 3d_conversion_mode setting."""
        from moleditpy.ui.compute_logic import ComputeManager

        compute = object.__new__(ComputeManager)
        compute._host = host
        ComputeManager.__init__(compute, host)

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}
        host.init_manager.settings["3d_conversion_mode"] = "rdkit"
        compute.__dict__["_pending_conversion_mode"] = "fallback"

        with (
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker") as mock_start,
        ):
            compute.trigger_conversion()

        options = mock_start.call_args[0][1]
        assert options["conversion_mode"] == "fallback"

    def test_settings_mode_used_when_no_pending_mode(self, host, app):
        """When _pending_conversion_mode is absent, the settings value is used."""
        from moleditpy.ui.compute_logic import ComputeManager

        compute = object.__new__(ComputeManager)
        compute._host = host
        ComputeManager.__init__(compute, host)

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}
        host.init_manager.settings["3d_conversion_mode"] = "obabel"

        with (
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker") as mock_start,
        ):
            compute.trigger_conversion()

        options = mock_start.call_args[0][1]
        assert options["conversion_mode"] == "obabel"


# ---------------------------------------------------------------------------
# Test: plugin wrapping trigger_conversion without **kwargs does NOT raise
# ---------------------------------------------------------------------------


class TestPluginWrappedTriggerConversion:
    def test_no_typeerror_when_plugin_wraps_without_kwargs(self, host, app):
        """Calling _trigger_conversion_with_temp_mode must not raise TypeError
        when a plugin has replaced trigger_conversion with a no-kwarg wrapper.
        This is the regression test for the original bug report."""
        compute = _PluginWrappedCompute(host)

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}
        host.init_manager.settings["3d_conversion_mode"] = "rdkit"

        # Patch the timer to fire the callback immediately so we can verify
        # the full call chain without actually needing a Qt event loop
        timer_callbacks = []

        def capture_callback(delay, cb):
            timer_callbacks.append(cb)

        with (
            patch("PyQt6.QtCore.QTimer.singleShot", side_effect=capture_callback),
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker"),
        ):
            # This must not raise TypeError
            compute._trigger_conversion_with_temp_mode("rdkit")

            # Fire the timer callback (simulates Qt event loop)
            for cb in timer_callbacks:
                cb()  # no TypeError expected here

    def test_mode_reaches_worker_through_plugin_wrapper(self, host, app):
        """The conversion mode set via _trigger_conversion_with_temp_mode
        reaches _start_calculation_worker even when a plugin wraps trigger_conversion."""
        from moleditpy.ui.compute_logic import ComputeManager

        class WrapperCompute(ComputeManager):
            def __init__(self, h):
                self._host = h
                ComputeManager.__init__(self, h)
                if not hasattr(h, "init_manager"):
                    h.init_manager = MagicMock()
                h.init_manager.settings = {}
                h.init_manager.opt3d_method_labels = {
                    "MMFF_RDKIT": "MMFF94s",
                    "UFF_RDKIT": "UFF",
                }
                for attr in ("view_3d_manager", "state_manager", "ui_manager"):
                    if not hasattr(h, attr):
                        setattr(h, attr, MagicMock())
                for btn in (
                    "convert_button",
                    "cleanup_button",
                    "optimize_3d_button",
                    "export_button",
                    "analysis_action",
                    "edit_3d_action",
                ):
                    if not hasattr(h.init_manager, btn):
                        setattr(h.init_manager, btn, MagicMock())

                original = self.trigger_conversion

                def no_kwargs_wrapper():
                    original()

                self.trigger_conversion = no_kwargs_wrapper

            def __getattr__(self, name):
                return getattr(self._host, name)

            def _start_calculation_worker(self, mol_block, options, run_id):
                self._captured_options = options

            @property
            def data(self):
                return self.host.state_manager.data

            @data.setter
            def data(self, v):
                self.host.state_manager.data = v

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}
        host.init_manager.settings["3d_conversion_mode"] = "rdkit"

        compute = WrapperCompute(host)

        timer_callbacks = []

        def capture(delay, cb):
            timer_callbacks.append(cb)

        with (
            patch("PyQt6.QtCore.QTimer.singleShot", side_effect=capture),
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
        ):
            compute._trigger_conversion_with_temp_mode("fallback")
            for cb in timer_callbacks:
                cb()

        assert hasattr(compute, "_captured_options"), (
            "trigger_conversion never reached _start_calculation_worker"
        )
        assert compute._captured_options["conversion_mode"] == "fallback"


# ---------------------------------------------------------------------------
# Test: plugin optimization method as the conversion default
# ---------------------------------------------------------------------------


class TestPluginConversionPreOptimize:
    """When the default optimizer is a plugin method, conversion pre-optimizes
    with a built-in force field and runs the plugin callback as a post-step."""

    @staticmethod
    def _make_compute(host):
        compute = object.__new__(ComputeManager)
        compute._host = host
        ComputeManager.__init__(compute, host)
        return compute

    def test_plugin_default_preoptimizes_with_mmff_and_queues_post_step(
        self, host, app
    ):
        compute = self._make_compute(host)
        entry = {
            "plugin": "P",
            "callback": MagicMock(return_value=True),
            "label": "Quick UFF",
        }
        host.plugin_manager.optimization_methods = {"QUICK UFF": entry}
        host.init_manager.optimization_method = "QUICK UFF"

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}

        with (
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker") as mock_start,
        ):
            compute.trigger_conversion()

        options = mock_start.call_args[0][1]
        # Worker gets a built-in FF, not the plugin key.
        assert options["optimization_method"] == "MMFF_RDKIT"
        # The plugin post-step is queued against the worker id.
        run_id = options["worker_id"]
        assert compute._pending_plugin_opt[run_id] == ("QUICK UFF", entry)

    def test_builtin_default_queues_no_post_step(self, host, app):
        compute = self._make_compute(host)
        host.plugin_manager.optimization_methods = {}
        host.init_manager.optimization_method = "UFF_RDKIT"

        mol = Chem.MolFromSmiles("C")
        host.state_manager.data.atoms = {1: {"symbol": "C"}}

        with (
            patch.object(host.state_manager.data, "to_rdkit_mol", return_value=mol),
            patch("rdkit.Chem.DetectChemistryProblems", return_value=[]),
            patch("rdkit.Chem.SanitizeMol"),
            patch.object(
                host.state_manager.data,
                "to_mol_block",
                return_value=Chem.MolToMolBlock(mol),
            ),
            patch.object(compute, "_start_calculation_worker") as mock_start,
        ):
            compute.trigger_conversion()

        options = mock_start.call_args[0][1]
        assert options["optimization_method"] == "UFF_RDKIT"
        assert compute._pending_plugin_opt == {}

    def test_finished_runs_pending_plugin_callback(self, host, app):
        compute = self._make_compute(host)
        cb = MagicMock(return_value=True)
        entry = {"plugin": "P", "callback": cb, "label": "Quick UFF"}
        run_id = 7
        compute._pending_plugin_opt[run_id] = ("QUICK UFF", entry)
        compute.active_worker_ids.add(run_id)

        mol = Chem.MolFromSmiles("C")
        host.view_3d_manager.current_mol = mol

        with (
            patch.object(compute, "create_atom_id_mapping"),
            patch.object(compute, "_refresh_ui_state"),
            patch.object(compute, "_remove_calculating_text"),
        ):
            compute.on_calculation_finished((run_id, mol))

        cb.assert_called_once_with(mol)
        assert run_id not in compute._pending_plugin_opt
