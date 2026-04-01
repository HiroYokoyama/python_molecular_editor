"""
Unit tests for DialogManager (dialog_logic.py).
Covers:
  - _get_preselected_atoms_3d
  - show_about_dialog
  - open_periodic_table_dialog
  - open_analysis_window
  - open_template_dialog
  - open_template_dialog_and_activate
  - save_2d_as_template
  - open_translation_dialog
  - open_move_group_dialog
  - open_align_plane_dialog
  - open_planarize_dialog
  - open_alignment_dialog
  - open_bond_length_dialog
  - open_angle_dialog
  - open_dihedral_dialog
  - open_mirror_dialog
  - open_settings_dialog
  - open_color_settings_dialog
  - open_constrained_optimization_dialog

Strategy: all dialog classes are patched at the moleditpy.ui.dialog_logic.*
namespace.  Signals and host sub-managers are fully mocked via DummyHost.
No local imports exist in dialog_logic.py so all classes are patchable at the
module level.
"""

import os
import sys
import json
import pytest
import tempfile
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication, QMessageBox

# Make the local moleditpy package discoverable
workspace_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(workspace_src_path) and workspace_src_path not in sys.path:
    sys.path.insert(0, workspace_src_path)

from moleditpy.ui.dialog_logic import DialogManager


# ---------------------------------------------------------------------------
# Shared QApplication (session-scoped)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


# ---------------------------------------------------------------------------
# DummyHost
# ---------------------------------------------------------------------------

class DummyHost:
    def __init__(self):
        self.statusBar_mock = MagicMock()
        self.is_xyz_derived = False
        self._template_dialog = None

        self.init_manager = MagicMock()
        self.init_manager.settings_dir = os.path.join(tempfile.gettempdir(), "moleditpy_test_settings")
        self.init_manager.settings = {}

        self.state_manager = MagicMock()
        self.state_manager.data = MagicMock()
        self.state_manager.data.atoms = {"a1": {}}     # non-empty by default
        self.state_manager.has_unsaved_changes = False

        self.view_3d_manager = MagicMock()
        self.view_3d_manager.current_mol = MagicMock()  # non-None by default

        self.edit_3d_manager = MagicMock()
        self.edit_3d_manager.measurement_mode = False
        self.edit_3d_manager.selected_atoms_for_measurement = []
        self.edit_3d_manager.active_3d_dialogs = []

        self.ui_manager = MagicMock()
        self.edit_actions_manager = MagicMock()

    def statusBar(self):
        return self.statusBar_mock


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def dm(qapp):
    return DialogManager(DummyHost())


# ===========================================================================
# _get_preselected_atoms_3d
# ===========================================================================

class TestGetPreselectedAtoms3D:

    def test_returns_empty_when_no_selection(self, dm):
        dm.host.edit_3d_manager.selected_atoms_for_measurement = []
        assert dm._get_preselected_atoms_3d() == []

    def test_returns_selected_atoms(self, dm):
        dm.host.edit_3d_manager.selected_atoms_for_measurement = [1, 2, 3]
        result = dm._get_preselected_atoms_3d()
        assert result == [1, 2, 3]

    def test_logs_error_when_edit_3d_manager_missing(self, qapp):
        host = DummyHost()
        del host.edit_3d_manager
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.logging.error") as mock_log:
            result = dlgm._get_preselected_atoms_3d()
        assert result == []
        mock_log.assert_called_once()


# ===========================================================================
# show_about_dialog
# ===========================================================================

class TestShowAboutDialog:

    def test_creates_and_execs_dialog(self, dm):
        with patch("moleditpy.ui.dialog_logic.AboutDialog") as MockAbout:
            instance = MagicMock()
            MockAbout.return_value = instance
            dm.show_about_dialog()
        MockAbout.assert_called_once_with(dm.host, dm.host)
        instance.exec.assert_called_once()


# ===========================================================================
# open_periodic_table_dialog
# ===========================================================================

class TestOpenPeriodicTableDialog:

    def test_creates_connects_and_execs(self, dm):
        with patch("moleditpy.ui.dialog_logic.PeriodicTableDialog") as MockPT:
            instance = MagicMock()
            MockPT.return_value = instance
            dm.open_periodic_table_dialog()
        MockPT.assert_called_once_with(dm.host)
        instance.element_selected.connect.assert_called_once_with(
            dm.host.ui_manager.set_atom_from_periodic_table
        )
        instance.exec.assert_called_once()

    def test_unchecks_tool_group_action(self, dm):
        checked = MagicMock()
        dm.host.init_manager.tool_group.checkedAction.return_value = checked
        with patch("moleditpy.ui.dialog_logic.PeriodicTableDialog"):
            dm.open_periodic_table_dialog()
        checked.setChecked.assert_called_with(False)


# ===========================================================================
# open_analysis_window
# ===========================================================================

class TestOpenAnalysisWindow:

    def test_opens_when_mol_exists(self, dm):
        with patch("moleditpy.ui.dialog_logic.AnalysisWindow") as MockAW:
            instance = MagicMock()
            MockAW.return_value = instance
            dm.open_analysis_window()
        MockAW.assert_called_once_with(
            dm.host.view_3d_manager.current_mol,
            dm.host,
            is_xyz_derived=dm.host.is_xyz_derived,
        )
        instance.exec.assert_called_once()

    def test_shows_error_when_no_mol(self, dm):
        dm.host.view_3d_manager.current_mol = None
        with patch("moleditpy.ui.dialog_logic.AnalysisWindow") as MockAW:
            dm.open_analysis_window()
        MockAW.assert_not_called()
        dm.host.statusBar_mock.showMessage.assert_called_once()
        msg = dm.host.statusBar_mock.showMessage.call_args[0][0]
        assert "3D" in msg or "generate" in msg.lower()


# ===========================================================================
# open_template_dialog
# ===========================================================================

class TestOpenTemplateDialog:

    def test_creates_and_execs(self, dm):
        with patch("moleditpy.ui.dialog_logic.UserTemplateDialog") as MockUT:
            instance = MagicMock()
            MockUT.return_value = instance
            dm.open_template_dialog()
        MockUT.assert_called_once_with(dm.host, dm.host)
        instance.exec.assert_called_once()


# ===========================================================================
# open_template_dialog_and_activate
# ===========================================================================

class TestOpenTemplateDialogAndActivate:

    def test_creates_new_dialog_when_none_exists(self, dm):
        dm.host._template_dialog = None
        with patch("moleditpy.ui.dialog_logic.UserTemplateDialog") as MockUT:
            instance = MagicMock()
            MockUT.return_value = instance
            dm.open_template_dialog_and_activate()
        MockUT.assert_called_once_with(dm.host, dm.host)
        instance.show.assert_called_once()
        instance.finished.connect.assert_called_once()

    def test_raises_existing_visible_dialog(self, dm):
        existing = MagicMock()
        existing.isHidden.return_value = False
        dm.host._template_dialog = existing
        with patch("moleditpy.ui.dialog_logic.UserTemplateDialog") as MockUT:
            dm.open_template_dialog_and_activate()
        MockUT.assert_not_called()
        existing.raise_.assert_called_once()
        existing.activateWindow.assert_called_once()

    def test_on_finished_sets_mode_when_template_selected(self, dm):
        dm.host._template_dialog = None
        captured_cb = []

        with patch("moleditpy.ui.dialog_logic.UserTemplateDialog") as MockUT:
            instance = MagicMock()
            MockUT.return_value = instance
            instance.finished.connect.side_effect = lambda cb: captured_cb.append(cb)
            dm.open_template_dialog_and_activate()

        # Simulate template selection on the newly stored dialog
        dm.host._template_dialog.selected_template = {"name": "benzene"}
        assert captured_cb, "finished.connect was never called"
        captured_cb[0]()

        dm.host.ui_manager.set_mode.assert_called_once_with("template_user_benzene")
        dm.host.statusBar_mock.showMessage.assert_called_once()

    def test_on_finished_noop_when_no_template_selected(self, dm):
        dm.host._template_dialog = None
        captured_cb = []

        with patch("moleditpy.ui.dialog_logic.UserTemplateDialog") as MockUT:
            instance = MagicMock()
            MockUT.return_value = instance
            instance.finished.connect.side_effect = lambda cb: captured_cb.append(cb)
            dm.open_template_dialog_and_activate()

        dm.host._template_dialog.selected_template = None
        captured_cb[0]()
        dm.host.ui_manager.set_mode.assert_not_called()


# ===========================================================================
# save_2d_as_template
# ===========================================================================

class TestSave2DAsTemplate:

    def test_warns_when_no_atoms(self, qapp):
        host = DummyHost()
        host.state_manager.data.atoms = {}
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.QMessageBox.warning") as mock_warn:
            dlgm.save_2d_as_template()
        mock_warn.assert_called_once()
        args = mock_warn.call_args[0]
        assert "No structure" in args[2] or "template" in args[2].lower()

    def test_noop_on_cancelled_input(self, dm):
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("", False)):
            dm.save_2d_as_template()
        dm.host.state_manager.data.to_template_dict.assert_not_called()

    def test_noop_on_blank_name(self, dm):
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("   ", True)):
            dm.save_2d_as_template()
        dm.host.state_manager.data.to_template_dict.assert_not_called()

    def test_saves_template_file(self, dm, tmp_path):
        dm.host.init_manager.settings_dir = str(tmp_path)
        dm.host.state_manager.data.to_template_dict.return_value = {"name": "mytemplate", "atoms": []}
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("mytemplate", True)), \
             patch("moleditpy.ui.dialog_logic.QMessageBox.information"):
            dm.save_2d_as_template()
        saved = tmp_path / "user-templates" / "mytemplate.pmetmplt"
        assert saved.exists()
        assert json.loads(saved.read_text())["name"] == "mytemplate"

    def test_overwrites_after_yes_confirmation(self, dm, tmp_path):
        dm.host.init_manager.settings_dir = str(tmp_path)
        tpl_dir = tmp_path / "user-templates"
        tpl_dir.mkdir()
        f = tpl_dir / "mytemplate.pmetmplt"
        f.write_text("{}")
        dm.host.state_manager.data.to_template_dict.return_value = {"name": "mytemplate", "atoms": []}
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("mytemplate", True)), \
             patch("moleditpy.ui.dialog_logic.QMessageBox.question",
                   return_value=QMessageBox.StandardButton.Yes), \
             patch("moleditpy.ui.dialog_logic.QMessageBox.information"):
            dm.save_2d_as_template()
        assert "atoms" in json.loads(f.read_text())

    def test_skips_overwrite_on_no_confirmation(self, dm, tmp_path):
        dm.host.init_manager.settings_dir = str(tmp_path)
        tpl_dir = tmp_path / "user-templates"
        tpl_dir.mkdir()
        f = tpl_dir / "mytemplate.pmetmplt"
        f.write_text('{"original": true}')
        dm.host.state_manager.data.to_template_dict.return_value = {"name": "mytemplate", "atoms": []}
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("mytemplate", True)), \
             patch("moleditpy.ui.dialog_logic.QMessageBox.question",
                   return_value=QMessageBox.StandardButton.No):
            dm.save_2d_as_template()
        assert json.loads(f.read_text()) == {"original": True}

    def test_shows_error_on_exception(self, dm):
        dm.host.state_manager.data.to_template_dict.side_effect = AttributeError("boom")
        with patch("moleditpy.ui.dialog_logic.QInputDialog.getText",
                   return_value=("mytemplate", True)), \
             patch("moleditpy.ui.dialog_logic.QMessageBox.critical") as mock_crit:
            dm.save_2d_as_template()
        mock_crit.assert_called_once()


# ===========================================================================
# Modeless geometry dialogs — shared helper
# ===========================================================================

def _assert_modeless(dm, method, cls_name, *args):
    """Assert dialog is appended, shown, and all 3 signals connected."""
    with patch(f"moleditpy.ui.dialog_logic.{cls_name}") as MockDlg:
        instance = MagicMock()
        MockDlg.return_value = instance
        getattr(dm, method)(*args)
    instance.show.assert_called_once()
    assert instance.accepted.connect.call_count == 2
    instance.finished.connect.assert_called_once()
    assert instance in dm.host.edit_3d_manager.active_3d_dialogs


class TestModelessGeometryDialogs:

    def test_open_translation_dialog(self, dm):
        _assert_modeless(dm, "open_translation_dialog", "TranslationDialog")

    def test_translation_disables_measurement_mode(self, qapp):
        host = DummyHost()
        host.edit_3d_manager.measurement_mode = True
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.TranslationDialog"):
            dlgm.open_translation_dialog()
        host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)

    def test_open_move_group_dialog(self, dm):
        _assert_modeless(dm, "open_move_group_dialog", "MoveGroupDialog")

    def test_open_align_plane_dialog(self, dm):
        _assert_modeless(dm, "open_align_plane_dialog", "AlignPlaneDialog", "xy")

    def test_align_plane_message_contains_plane(self, qapp):
        host = DummyHost()
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.AlignPlaneDialog") as MockDlg:
            instance = MagicMock()
            MockDlg.return_value = instance
            cbs = []
            instance.accepted.connect.side_effect = lambda cb: cbs.append(cb)
            dlgm.open_align_plane_dialog("xz")
        cbs[0]()
        assert "XZ" in host.statusBar_mock.showMessage.call_args[0][0]

    def test_open_planarize_dialog(self, dm):
        _assert_modeless(dm, "open_planarize_dialog", "PlanarizeDialog")

    def test_open_alignment_dialog(self, dm):
        _assert_modeless(dm, "open_alignment_dialog", "AlignmentDialog", "x")

    def test_alignment_message_contains_axis(self, qapp):
        host = DummyHost()
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.AlignmentDialog") as MockDlg:
            instance = MagicMock()
            MockDlg.return_value = instance
            cbs = []
            instance.accepted.connect.side_effect = lambda cb: cbs.append(cb)
            dlgm.open_alignment_dialog("y")
        cbs[0]()
        assert "Y" in host.statusBar_mock.showMessage.call_args[0][0]

    def test_open_bond_length_dialog(self, dm):
        _assert_modeless(dm, "open_bond_length_dialog", "BondLengthDialog")

    def test_open_angle_dialog(self, dm):
        _assert_modeless(dm, "open_angle_dialog", "AngleDialog")

    def test_open_dihedral_dialog(self, dm):
        _assert_modeless(dm, "open_dihedral_dialog", "DihedralDialog")

    def test_accepted_status_messages(self, qapp):
        """Each dialog's first accepted lambda posts the right status bar message."""
        cases = [
            ("open_translation_dialog",   "TranslationDialog", [],     "Translation applied."),
            ("open_move_group_dialog",     "MoveGroupDialog",   [],     "Group transformation applied."),
            ("open_bond_length_dialog",    "BondLengthDialog",  [],     "Bond length adjusted."),
            ("open_angle_dialog",          "AngleDialog",       [],     "Angle adjusted."),
            ("open_dihedral_dialog",       "DihedralDialog",    [],     "Dihedral angle adjusted."),
            ("open_planarize_dialog",      "PlanarizeDialog",   [],     "Selection planarized to best-fit plane."),
        ]
        for method, cls, args, expected in cases:
            host = DummyHost()
            dlgm = DialogManager(host)
            with patch(f"moleditpy.ui.dialog_logic.{cls}") as MockDlg:
                instance = MagicMock()
                MockDlg.return_value = instance
                cbs = []
                instance.accepted.connect.side_effect = lambda cb: cbs.append(cb)
                getattr(dlgm, method)(*args)
            cbs[0]()
            actual = host.statusBar_mock.showMessage.call_args[0][0]
            assert actual == expected, f"{method}: {actual!r} != {expected!r}"

    def test_accepted_pushes_undo_state(self, qapp):
        """Second accepted lambda calls push_undo_state."""
        host = DummyHost()
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.TranslationDialog") as MockDlg:
            instance = MagicMock()
            MockDlg.return_value = instance
            cbs = []
            instance.accepted.connect.side_effect = lambda cb: cbs.append(cb)
            dlgm.open_translation_dialog()
        cbs[1]()
        host.edit_actions_manager.push_undo_state.assert_called_once()

    def test_finished_removes_dialog_from_list(self, qapp):
        """finished lambda calls remove_dialog_from_list with this dialog."""
        host = DummyHost()
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.TranslationDialog") as MockDlg:
            instance = MagicMock()
            MockDlg.return_value = instance
            fin_cbs = []
            instance.finished.connect.side_effect = lambda cb: fin_cbs.append(cb)
            dlgm.open_translation_dialog()
        fin_cbs[0]()
        host.edit_3d_manager.remove_dialog_from_list.assert_called_once_with(instance)


# ===========================================================================
# open_mirror_dialog
# ===========================================================================

class TestOpenMirrorDialog:

    def test_opens_when_mol_exists(self, dm):
        with patch("moleditpy.ui.dialog_logic.MirrorDialog") as MockM:
            instance = MagicMock()
            MockM.return_value = instance
            dm.open_mirror_dialog()
        MockM.assert_called_once_with(dm.host.view_3d_manager.current_mol, dm.host)
        instance.exec.assert_called_once()

    def test_shows_error_when_no_mol(self, dm):
        dm.host.view_3d_manager.current_mol = None
        with patch("moleditpy.ui.dialog_logic.MirrorDialog") as MockM:
            dm.open_mirror_dialog()
        MockM.assert_not_called()
        dm.host.statusBar_mock.showMessage.assert_called_with("No 3D molecule loaded.")

    def test_disables_measurement_mode(self, qapp):
        host = DummyHost()
        host.edit_3d_manager.measurement_mode = True
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.MirrorDialog"):
            dlgm.open_mirror_dialog()
        host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)


# ===========================================================================
# open_settings_dialog
# ===========================================================================

class TestOpenSettingsDialog:

    def test_creates_and_execs(self, dm):
        with patch("moleditpy.ui.dialog_logic.SettingsDialog") as MockSD:
            instance = MagicMock()
            MockSD.return_value = instance
            dm.open_settings_dialog()
        MockSD.assert_called_once_with(dm.host.init_manager.settings, parent=dm.host)
        instance.exec.assert_called_once()


# ===========================================================================
# open_color_settings_dialog
# ===========================================================================

class TestOpenColorSettingsDialog:

    def test_creates_and_execs(self, dm):
        with patch("moleditpy.ui.dialog_logic.ColorSettingsDialog") as MockCD:
            instance = MagicMock()
            MockCD.return_value = instance
            dm.open_color_settings_dialog()
        MockCD.assert_called_once_with(dm.host.init_manager.settings, parent=dm.host)
        instance.exec.assert_called_once()


# ===========================================================================
# open_constrained_optimization_dialog
# ===========================================================================

class TestOpenConstrainedOptimizationDialog:

    def test_opens_when_mol_exists(self, dm):
        with patch("moleditpy.ui.dialog_logic.ConstrainedOptimizationDialog") as MockCO:
            instance = MagicMock()
            MockCO.return_value = instance
            dm.open_constrained_optimization_dialog()
        MockCO.assert_called_once_with(
            dm.host.view_3d_manager.current_mol, dm.host, parent=dm.host
        )
        instance.show.assert_called_once()
        instance.finished.connect.assert_called_once()
        assert instance in dm.host.edit_3d_manager.active_3d_dialogs

    def test_shows_error_when_no_mol(self, dm):
        dm.host.view_3d_manager.current_mol = None
        with patch("moleditpy.ui.dialog_logic.ConstrainedOptimizationDialog") as MockCO:
            dm.open_constrained_optimization_dialog()
        MockCO.assert_not_called()
        dm.host.statusBar_mock.showMessage.assert_called_with("No 3D molecule loaded.")

    def test_disables_measurement_mode(self, qapp):
        host = DummyHost()
        host.edit_3d_manager.measurement_mode = True
        dlgm = DialogManager(host)
        with patch("moleditpy.ui.dialog_logic.ConstrainedOptimizationDialog"):
            dlgm.open_constrained_optimization_dialog()
        host.edit_3d_manager.toggle_measurement_mode.assert_called_with(False)

    def test_finished_removes_from_active_dialogs(self, dm):
        with patch("moleditpy.ui.dialog_logic.ConstrainedOptimizationDialog") as MockCO:
            instance = MagicMock()
            MockCO.return_value = instance
            fin_cbs = []
            instance.finished.connect.side_effect = lambda cb: fin_cbs.append(cb)
            dm.open_constrained_optimization_dialog()
        fin_cbs[0]()
        dm.host.edit_3d_manager.remove_dialog_from_list.assert_called_once_with(instance)
