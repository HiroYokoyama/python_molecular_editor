"""
Extended tests for ColorSettingsDialog (ui/color_settings_dialog.py).

Covers uncovered paths:
  - pick_bs_bond_color: valid / invalid QColorDialog result
  - reset_all: restores bs button color from parent_window defaults
  - reset_all: without parent_window defaults to #7F7F7F
  - apply_changes: _reset_all_flag deletes cpk_colors key
  - apply_changes: no parent → early return (no crash)
  - apply_changes: mol present → calls draw_molecule_3d
  - apply_changes: _reset_all_flag without changed_bs_color → resets bond color
  - apply_changes: updates 2D scene items via update_style
  - accept: calls apply_changes then super().accept
  - on_element_clicked: invalid color → no change to changed_cpk
  - init: override color from settings applied to button
"""

from unittest.mock import patch, MagicMock
from moleditpy.ui.color_settings_dialog import ColorSettingsDialog
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


def _make_parent():
    parent = MagicMock()
    parent.init_manager.settings = dict(DEFAULT_SETTINGS)
    parent.init_manager.settings["cpk_colors"] = {}
    parent.default_settings = dict(DEFAULT_SETTINGS)
    parent.view_3d_manager.current_mol = None
    parent.statusBar.return_value = MagicMock()
    return parent


# ---------------------------------------------------------------------------
# pick_bs_bond_color
# ---------------------------------------------------------------------------

@patch("moleditpy.ui.color_settings_dialog.QColorDialog.getColor")
def test_pick_bs_bond_color_valid_updates_changed_bs(mock_get, app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_color = MagicMock()
    mock_color.isValid.return_value = True
    mock_color.name.return_value = "#aabbcc"
    mock_get.return_value = mock_color

    dialog.pick_bs_bond_color()
    assert dialog.changed_bs_color == "#aabbcc"
    assert "background-color: #aabbcc" in dialog.bs_button.styleSheet()


@patch("moleditpy.ui.color_settings_dialog.QColorDialog.getColor")
def test_pick_bs_bond_color_invalid_no_change(mock_get, app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_color = MagicMock()
    mock_color.isValid.return_value = False
    mock_get.return_value = mock_color

    dialog.pick_bs_bond_color()
    assert dialog.changed_bs_color is None


# ---------------------------------------------------------------------------
# reset_all
# ---------------------------------------------------------------------------

def test_reset_all_with_parent_uses_default_bond_color(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = MagicMock()
    parent.default_settings = {"ball_stick_bond_color": "#123456"}
    dialog.parent_window = parent

    dialog.reset_all()

    assert dialog.changed_bs_color == "#123456"
    assert "background-color: #123456" in dialog.bs_button.styleSheet()


def test_reset_all_without_parent_defaults_to_gray(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    dialog.parent_window = None

    dialog.reset_all()

    assert dialog.changed_bs_color == "#7F7F7F"
    assert "background-color: #7f7f7f" in dialog.bs_button.styleSheet().lower()


def test_reset_all_restores_element_buttons(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    dialog.parent_window = None
    dialog.changed_cpk["C"] = "#ff0000"

    dialog.reset_all()

    assert dialog.changed_cpk == {}
    assert dialog._reset_all_flag is True


# ---------------------------------------------------------------------------
# apply_changes
# ---------------------------------------------------------------------------

def test_apply_changes_no_parent_returns_early(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    dialog.parent_window = None
    dialog.apply_changes()   # should not raise


def test_apply_changes_reset_flag_deletes_cpk_colors(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    parent.init_manager.settings["cpk_colors"] = {"C": "#ff0000"}
    dialog.parent_window = parent
    dialog._reset_all_flag = True

    dialog.apply_changes()

    assert "cpk_colors" not in parent.init_manager.settings


def test_apply_changes_with_mol_calls_draw_molecule_3d(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    mol = MagicMock()
    parent.view_3d_manager.current_mol = mol
    dialog.parent_window = parent
    dialog.changed_cpk["C"] = "#ff0000"

    dialog.apply_changes()

    parent.view_3d_manager.draw_molecule_3d.assert_called()


def test_apply_changes_reset_flag_resets_bond_color(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    dialog.parent_window = parent
    dialog._reset_all_flag = True
    dialog.changed_bs_color = None   # no explicit bs change

    dialog.apply_changes()

    assert parent.init_manager.settings.get("ball_stick_bond_color") == "#7F7F7F"


def test_apply_changes_updates_2d_scene_items(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    dialog.parent_window = parent

    item_with_style = MagicMock()
    item_with_style.update_style = MagicMock()
    del item_with_style.update   # force update_style path
    item_plain = MagicMock()
    del item_plain.update_style   # no update_style → falls back to update()

    scene = MagicMock()
    scene.items.return_value = [item_with_style, item_plain]
    parent.init_manager.scene = scene

    dialog.apply_changes()

    item_with_style.update_style.assert_called()


def test_apply_changes_calls_update_cpk_colors(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    dialog.parent_window = parent

    dialog.apply_changes()

    parent.init_manager.update_cpk_colors_from_settings.assert_called()


def test_apply_changes_calls_save_settings_when_cpk_changed(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    parent = _make_parent()
    dialog.parent_window = parent
    dialog.changed_cpk["N"] = "#0000ff"

    dialog.apply_changes()

    parent.init_manager.save_settings.assert_called()


# ---------------------------------------------------------------------------
# accept
# ---------------------------------------------------------------------------

def test_accept_calls_apply_then_super(app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    dialog.parent_window = None

    with patch.object(dialog, "apply_changes") as mock_apply, \
         patch("moleditpy.ui.color_settings_dialog.QDialog.accept") as mock_super:
        dialog.accept()

    mock_apply.assert_called_once()
    mock_super.assert_called_once()


# ---------------------------------------------------------------------------
# on_element_clicked — invalid color
# ---------------------------------------------------------------------------

@patch("moleditpy.ui.color_settings_dialog.QColorDialog.getColor")
def test_on_element_clicked_invalid_color_no_change(mock_get, app):
    dialog = ColorSettingsDialog(current_settings={}, parent=None)
    mock_color = MagicMock()
    mock_color.isValid.return_value = False
    mock_get.return_value = mock_color

    btn = dialog.element_buttons["C"]
    with patch.object(dialog, "sender", return_value=btn):
        dialog.on_element_clicked()

    assert "C" not in dialog.changed_cpk


# ---------------------------------------------------------------------------
# init — override color from settings
# ---------------------------------------------------------------------------

def test_init_cpk_override_applied_to_button(app):
    settings = {"cpk_colors": {"C": "#112233"}}
    dialog = ColorSettingsDialog(current_settings=settings, parent=None)
    style = dialog.element_buttons["C"].styleSheet()
    assert "#112233" in style
