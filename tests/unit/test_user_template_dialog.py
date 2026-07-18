"""Unit tests for UserTemplateDialog template management."""

import json
import os

import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QMessageBox

from moleditpy.ui.user_template_dialog import UserTemplateDialog


def _template_payload(name="Benzene"):
    return {
        "format": "PME Template",
        "version": "1.0",
        "name": name,
        "atoms": [
            {"id": 1, "symbol": "C", "x": 0.0, "y": 0.0, "charge": 0, "radical": 0},
            {"id": 2, "symbol": "C", "x": 50.0, "y": 0.0, "charge": 0, "radical": 0},
        ],
        "bonds": [{"atom1": 1, "atom2": 2, "order": 1, "stereo": 0}],
    }


@pytest.fixture
def make_dialog(app, tmp_path):
    created = []

    def _factory():
        mw = MagicMock()
        mw.init_manager.settings_dir = str(tmp_path)
        dlg = UserTemplateDialog(mw, parent=None)
        created.append(dlg)
        return dlg, mw

    yield _factory

    for dlg in created:
        dlg.selected_template = None
        dlg.close()


def _write_template(tmp_path, name="Benzene", filename=None):
    template_dir = tmp_path / "user-templates"
    template_dir.mkdir(exist_ok=True)
    path = template_dir / (filename or f"{name}.pmetmplt")
    path.write_text(json.dumps(_template_payload(name)), encoding="utf-8")
    return path


class TestTemplateFiles:
    def test_get_template_directory_created_on_demand(self, make_dialog, tmp_path):
        dlg, _mw = make_dialog()
        result = dlg.get_template_directory()
        assert result == os.path.join(str(tmp_path), "user-templates")
        assert os.path.isdir(result)

    def test_load_template_file_valid_json(self, make_dialog, tmp_path):
        dlg, _mw = make_dialog()
        path = _write_template(tmp_path)
        data = dlg.load_template_file(str(path))
        assert data["name"] == "Benzene"

    def test_load_template_file_invalid_json_returns_none(self, make_dialog, tmp_path):
        dlg, _mw = make_dialog()
        bad = tmp_path / "bad.pmetmplt"
        bad.write_text("{not json", encoding="utf-8")
        assert dlg.load_template_file(str(bad)) is None

    def test_save_template_file_roundtrip(self, make_dialog, tmp_path):
        dlg, _mw = make_dialog()
        path = tmp_path / "out.pmetmplt"
        assert dlg.save_template_file(str(path), _template_payload()) is True
        assert json.loads(path.read_text(encoding="utf-8"))["name"] == "Benzene"

    def test_save_template_file_bad_path_returns_false(self, make_dialog, tmp_path):
        dlg, _mw = make_dialog()
        bad_path = tmp_path / "no_such_dir" / "out.pmetmplt"
        assert dlg.save_template_file(str(bad_path), {}) is False


class TestLoadUserTemplates:
    def test_loads_only_pmetmplt_files(self, make_dialog, tmp_path):
        _write_template(tmp_path, "One")
        _write_template(tmp_path, "Two")
        (tmp_path / "user-templates" / "notes.txt").write_text("skip me")
        dlg, _mw = make_dialog()

        names = sorted(t["name"] for t in dlg.user_templates)
        assert names == ["One", "Two"]
        for t in dlg.user_templates:
            assert t["filename"].endswith(".pmetmplt")
            assert os.path.exists(t["filepath"])

    def test_grid_gets_one_widget_per_template(self, make_dialog, tmp_path):
        _write_template(tmp_path, "One")
        _write_template(tmp_path, "Two")
        dlg, _mw = make_dialog()
        assert dlg.template_layout.count() == 2

    def test_corrupt_file_skipped(self, make_dialog, tmp_path):
        _write_template(tmp_path, "Good")
        (tmp_path / "user-templates" / "bad.pmetmplt").write_text("{oops")
        dlg, _mw = make_dialog()
        assert [t["name"] for t in dlg.user_templates] == ["Good"]


class TestConvertStructureToTemplate:
    def test_converts_atoms_and_bonds(self, make_dialog):
        dlg, mw = make_dialog()
        mw.state_manager.data.atoms = {
            1: {"symbol": "N", "pos": (10.0, 20.0), "charge": 1, "radical": 0},
            2: {"symbol": "C", "pos": (60.0, 20.0)},
        }
        mw.state_manager.data.bonds = {(1, 2): {"order": 2, "stereo": 1}}

        data = dlg.convert_structure_to_template("MyFrag")

        assert data["format"] == "PME Template"
        assert data["name"] == "MyFrag"
        atoms = {a["id"]: a for a in data["atoms"]}
        assert atoms[1]["symbol"] == "N"
        assert atoms[1]["charge"] == 1
        assert atoms[2]["charge"] == 0  # default
        assert atoms[1]["x"] == 10.0 and atoms[1]["y"] == 20.0
        assert data["bonds"] == [{"atom1": 1, "atom2": 2, "order": 2, "stereo": 1}]


class TestSaveCurrentAsTemplate:
    def test_warns_when_no_structure(self, make_dialog):
        dlg, mw = make_dialog()
        mw.state_manager.data.atoms = {}
        with patch("moleditpy.ui.user_template_dialog.QMessageBox.warning") as warn:
            dlg.save_current_as_template()
        warn.assert_called_once()

    def test_cancelled_name_dialog_saves_nothing(self, make_dialog, tmp_path):
        dlg, mw = make_dialog()
        mw.state_manager.data.atoms = {
            1: {"symbol": "C", "pos": (0.0, 0.0)},
        }
        mw.state_manager.data.bonds = {}
        with patch(
            "moleditpy.ui.user_template_dialog.QInputDialog.getText",
            return_value=("Ignored", False),
        ):
            dlg.save_current_as_template()
        assert list((tmp_path / "user-templates").glob("*.pmetmplt")) == []

    def test_saves_named_template_and_reloads(self, make_dialog, tmp_path):
        dlg, mw = make_dialog()
        mw.state_manager.data.atoms = {
            1: {"symbol": "C", "pos": (0.0, 0.0)},
        }
        mw.state_manager.data.bonds = {}
        with (
            patch(
                "moleditpy.ui.user_template_dialog.QInputDialog.getText",
                return_value=("My Frag", True),
            ),
            patch("moleditpy.ui.user_template_dialog.QMessageBox.information") as info,
        ):
            dlg.save_current_as_template()

        saved = tmp_path / "user-templates" / "My_Frag.pmetmplt"
        assert saved.exists()
        assert json.loads(saved.read_text(encoding="utf-8"))["name"] == "My Frag"
        info.assert_called_once()
        assert [t["name"] for t in dlg.user_templates] == ["My Frag"]

    def test_overwrite_declined_keeps_existing_file(self, make_dialog, tmp_path):
        existing = _write_template(tmp_path, "My Frag", filename="My_Frag.pmetmplt")
        original_content = existing.read_text(encoding="utf-8")
        dlg, mw = make_dialog()
        mw.state_manager.data.atoms = {
            1: {"symbol": "O", "pos": (5.0, 5.0)},
        }
        mw.state_manager.data.bonds = {}
        with (
            patch(
                "moleditpy.ui.user_template_dialog.QInputDialog.getText",
                return_value=("My Frag", True),
            ),
            patch(
                "moleditpy.ui.user_template_dialog.QMessageBox.question",
                return_value=QMessageBox.StandardButton.No,
            ),
        ):
            dlg.save_current_as_template()
        assert existing.read_text(encoding="utf-8") == original_content


class TestDeleteSelectedTemplate:
    def test_noop_without_selection(self, make_dialog):
        dlg, _mw = make_dialog()
        with patch(
            "moleditpy.ui.user_template_dialog.QMessageBox.question"
        ) as question:
            dlg.delete_selected_template()
        question.assert_not_called()

    def test_confirmed_delete_removes_file_and_resets(self, make_dialog, tmp_path):
        path = _write_template(tmp_path, "Doomed")
        dlg, _mw = make_dialog()
        dlg.selected_template = dlg.user_templates[0]
        dlg.delete_button.setEnabled(True)
        with (
            patch(
                "moleditpy.ui.user_template_dialog.QMessageBox.question",
                return_value=QMessageBox.StandardButton.Yes,
            ),
            patch("moleditpy.ui.user_template_dialog.QMessageBox.information"),
        ):
            dlg.delete_selected_template()
        assert not path.exists()
        assert dlg.selected_template is None
        assert not dlg.delete_button.isEnabled()

    def test_declined_delete_keeps_file(self, make_dialog, tmp_path):
        path = _write_template(tmp_path, "Kept")
        dlg, _mw = make_dialog()
        dlg.selected_template = dlg.user_templates[0]
        with patch(
            "moleditpy.ui.user_template_dialog.QMessageBox.question",
            return_value=QMessageBox.StandardButton.No,
        ):
            dlg.delete_selected_template()
        assert path.exists()


class TestTemplateModeLifecycle:
    def test_use_template_activates_and_selects(self, make_dialog):
        dlg, _mw = make_dialog()
        dlg._activate_template_mode = MagicMock()
        payload = _template_payload()
        dlg.use_template(payload)
        dlg._activate_template_mode.assert_called_once_with(payload)
        assert dlg.selected_template is payload

    def test_cleanup_template_mode_resets_scene_and_mode(self, make_dialog):
        dlg, mw = make_dialog()
        dlg.selected_template = _template_payload()
        dlg.delete_button.setEnabled(True)
        scene = mw.init_manager.scene
        scene.views.return_value = []

        dlg.cleanup_template_mode()

        assert dlg.selected_template is None
        assert not dlg.delete_button.isEnabled()
        mw.ui_manager.set_mode_and_update_toolbar.assert_called_with("atom_C")
        assert scene.mode == "atom_C"
        assert scene.current_atom_symbol == "C"
        assert scene.user_template_data is None
        assert scene.template_context == {}
        scene.clear_template_preview.assert_called_once()
        scene.template_preview.hide.assert_called_once()
        scene.update.assert_called()

    def test_close_event_triggers_cleanup(self, make_dialog):
        dlg, _mw = make_dialog()
        dlg.cleanup_template_mode = MagicMock()
        dlg.close()
        dlg.cleanup_template_mode.assert_called_once()
