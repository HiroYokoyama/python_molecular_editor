"""
Unit tests for MoveDialogMixin (ui/move_dialog_mixin.py).

The mixin holds the translation, rotation, highlighting and selection-label
code that MoveGroupDialog and MoveSelectedAtomsDialog used to carry a copy of
each. These tests exercise it through both dialogs, so a change that suits one
and breaks the other cannot pass.

Covers:
  - MRO: the mixin must not shadow BasePickingDialog's _update_molecule_geometry
    or _push_undo, and its clear_atom_labels must still reach the base class
  - Per-dialog constants: highlight actor name, warning text, label wording
  - update_display: empty, short, and the >5-atom ellipsis boundary
  - apply_translation / apply_rotation: geometry, guards, and undo ordering
  - clear_selection: both atom sets and the drag state
"""

import os
import sys

import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from PyQt6.QtWidgets import QApplication

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from moleditpy.ui.move_dialog_mixin import MoveDialogMixin  # noqa: E402

_MIXIN = "moleditpy.ui.move_dialog_mixin"


def _propane():
    """Propane with 3D coords: 3 heavy atoms plus hydrogens."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _heptane():
    """Seven heavy atoms, so update_display crosses its five-atom ellipsis boundary."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCCCCCC")
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


def _dialog_classes():
    from moleditpy.ui.move_group_dialog import MoveGroupDialog
    from moleditpy.ui.move_selected_atoms_dialog import MoveSelectedAtomsDialog

    return [MoveGroupDialog, MoveSelectedAtomsDialog]


@pytest.fixture
def make_dialog(qapp):
    """Build either dialog with the 3D view mocked out."""
    created = []

    def _factory(cls, mol=None):
        _mol = mol if mol is not None else _propane()
        mw = MagicMock()
        mw._picking_consumed = False
        mw.view_3d_manager.atom_positions_3d = np.array(
            _mol.GetConformer().GetPositions(), dtype=float
        )
        mw.view_3d_manager.plotter = MagicMock()
        with (
            patch.object(cls, "show_atom_labels"),
            patch.object(cls, "clear_atom_labels"),
        ):
            dlg = cls(_mol, mw)
        created.append(dlg)
        return dlg, _mol, mw

    yield _factory

    for dlg in created:
        try:
            dlg.picking_enabled = False
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Wiring: the parts that silently break if the MRO is wrong
# ---------------------------------------------------------------------------


class TestMixinWiring:
    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_both_dialogs_take_their_behaviour_from_the_mixin(self, cls):
        """Neither dialog may keep its own copy of a method the mixin owns."""
        shared = (
            "update_display",
            "show_atom_labels",
            "clear_atom_labels",
            "reset_translation_inputs",
            "reset_rotation_inputs",
            "apply_translation",
            "apply_rotation",
            "clear_selection",
        )
        for name in shared:
            assert name not in vars(cls), f"{cls.__name__} still defines {name}"
            assert getattr(cls, name) is getattr(MoveDialogMixin, name)

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_mixin_does_not_shadow_the_base_class_helpers(self, cls):
        """_update_molecule_geometry and _push_undo must stay the base class's.

        The mixin annotates them so it can call them. Assigning them instead --
        a stub, or a Callable default -- would put them ahead of
        BasePickingDialog in the MRO, and every translation would silently stop
        writing coordinates back.
        """
        from moleditpy.ui.base_picking_dialog import BasePickingDialog

        for name in ("_update_molecule_geometry", "_push_undo"):
            assert name not in vars(MoveDialogMixin)
            assert getattr(cls, name) is getattr(BasePickingDialog, name)

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_clear_atom_labels_still_reaches_the_base_class(self, cls, make_dialog):
        """The mixin overrides it, so the base class's own cleanup must be called."""
        dlg, _mol, _mw = make_dialog(cls)
        from moleditpy.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin

        with patch.object(Dialog3DPickingMixin, "clear_atom_labels") as base:
            MoveDialogMixin.clear_atom_labels(dlg)

        base.assert_called_once()

    def test_the_two_dialogs_do_not_share_a_highlight_actor_name(self):
        """Both highlights live in one plotter; a shared name would erase the other."""
        group, selected = _dialog_classes()
        assert group.HIGHLIGHT_ACTOR != selected.HIGHLIGHT_ACTOR


# ---------------------------------------------------------------------------
# update_display
# ---------------------------------------------------------------------------


class TestUpdateDisplay:
    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_empty_selection_uses_the_dialog_s_own_wording(self, cls, make_dialog):
        dlg, _mol, _mw = make_dialog(cls)
        dlg.group_atoms.clear()

        dlg.update_display()

        assert dlg.selection_label.text() == cls.NO_SELECTION_TEXT

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_lists_every_atom_up_to_five(self, cls, make_dialog):
        dlg, _mol, _mw = make_dialog(cls)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0, 1, 2})

        dlg.update_display()

        text = dlg.selection_label.text()
        assert text.startswith(f"{cls.SELECTION_PREFIX}: 3 atoms - ")
        assert "C(0), C(1), C(2)" in text
        assert "..." not in text

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_more_than_five_atoms_is_truncated_with_an_ellipsis(self, cls, make_dialog):
        """Seven atoms: the count stays honest while the list is cut to five."""
        dlg, _mol, _mw = make_dialog(cls, _heptane())
        dlg.group_atoms.clear()
        dlg.group_atoms.update(range(7))

        dlg.update_display()

        text = dlg.selection_label.text()
        assert f"{cls.SELECTION_PREFIX}: 7 atoms - " in text
        assert text.endswith(" ...")
        assert "C(5)" not in text


# ---------------------------------------------------------------------------
# Translation and rotation
# ---------------------------------------------------------------------------


class TestApplyTranslation:
    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_moves_only_the_selected_atoms(self, cls, make_dialog):
        dlg, mol, _mw = make_dialog(cls)
        before = np.array(mol.GetConformer().GetPositions(), dtype=float)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0, 1})
        dlg.x_trans_input.setText("1.0")
        dlg.y_trans_input.setText("-2.0")
        dlg.z_trans_input.setText("0.5")

        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_translation()

        after = np.array(mol.GetConformer().GetPositions(), dtype=float)
        assert np.allclose(after[[0, 1]] - before[[0, 1]], [1.0, -2.0, 0.5])
        assert np.allclose(after[2:], before[2:])

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_undo_is_pushed_after_the_move_not_before(self, cls, make_dialog):
        """A snapshot taken first would record the geometry the user is leaving."""
        dlg, _mol, _mw = make_dialog(cls)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0})
        dlg.x_trans_input.setText("1.0")
        order = []

        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(
                type(dlg),
                "_update_molecule_geometry",
                lambda *_a: order.append("write"),
            ),
            patch.object(type(dlg), "_push_undo", lambda *_a: order.append("undo")),
        ):
            dlg.apply_translation()

        assert order == ["write", "undo"]

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_empty_selection_warns_in_the_dialog_s_own_words(self, cls, make_dialog):
        dlg, _mol, _mw = make_dialog(cls)
        dlg.group_atoms.clear()

        with patch(f"{_MIXIN}.QMessageBox") as msg_box:
            dlg.apply_translation()

        assert msg_box.warning.call_args[0][2] == cls.EMPTY_SELECTION_WARNING

    @pytest.mark.parametrize("cls", _dialog_classes())
    @pytest.mark.parametrize("text", ["not a number", "nan", "inf", "-inf"])
    def test_unparsable_input_warns_and_leaves_the_geometry_alone(
        self, cls, make_dialog, text
    ):
        # float() accepts "nan"/"inf", which used to move atoms to NaN/inf.
        dlg, mol, _mw = make_dialog(cls)
        before = np.array(mol.GetConformer().GetPositions(), dtype=float)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0})
        dlg.x_trans_input.setText(text)

        with patch(f"{_MIXIN}.QMessageBox") as msg_box:
            dlg.apply_translation()

        assert "valid translation" in msg_box.warning.call_args[0][2]
        assert np.allclose(np.array(mol.GetConformer().GetPositions()), before)


class TestApplyRotation:
    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_ninety_degrees_about_z_swaps_x_and_y_around_the_centroid(
        self, cls, make_dialog
    ):
        dlg, mol, _mw = make_dialog(cls)
        before = np.array(mol.GetConformer().GetPositions(), dtype=float)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0, 1, 2})
        dlg.z_rot_input.setText("90.0")
        centroid = before[[0, 1, 2]].mean(axis=0)

        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_rotation()

        after = np.array(mol.GetConformer().GetPositions(), dtype=float)
        for idx in (0, 1, 2):
            rel_before = before[idx] - centroid
            rel_after = after[idx] - centroid
            assert np.allclose(rel_after[0], -rel_before[1], atol=1e-8)
            assert np.allclose(rel_after[1], rel_before[0], atol=1e-8)
            assert np.allclose(rel_after[2], rel_before[2], atol=1e-8)

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_rotation_is_rigid_and_keeps_the_centroid(self, cls, make_dialog):
        """Distances within the selection and its centre of mass must not move."""
        dlg, mol, _mw = make_dialog(cls)
        before = np.array(mol.GetConformer().GetPositions(), dtype=float)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0, 1, 2})
        for widget, value in (
            (dlg.x_rot_input, "31.0"),
            (dlg.y_rot_input, "-17.0"),
            (dlg.z_rot_input, "64.0"),
        ):
            widget.setText(value)

        with patch.object(type(dlg), "show_atom_labels"):
            dlg.apply_rotation()

        after = np.array(mol.GetConformer().GetPositions(), dtype=float)
        sel = [0, 1, 2]
        assert np.allclose(before[sel].mean(axis=0), after[sel].mean(axis=0))
        assert np.allclose(
            np.linalg.norm(before[0] - before[2]),
            np.linalg.norm(after[0] - after[2]),
        )

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_empty_selection_warns_in_the_dialog_s_own_words(self, cls, make_dialog):
        dlg, _mol, _mw = make_dialog(cls)
        dlg.group_atoms.clear()

        with patch(f"{_MIXIN}.QMessageBox") as msg_box:
            dlg.apply_rotation()

        assert msg_box.warning.call_args[0][2] == cls.EMPTY_SELECTION_WARNING

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_unparsable_input_warns_and_leaves_the_geometry_alone(
        self, cls, make_dialog
    ):
        dlg, mol, _mw = make_dialog(cls)
        before = np.array(mol.GetConformer().GetPositions(), dtype=float)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0})
        dlg.x_rot_input.setText("")

        with patch(f"{_MIXIN}.QMessageBox") as msg_box:
            dlg.apply_rotation()

        assert "valid rotation" in msg_box.warning.call_args[0][2]
        assert np.allclose(np.array(mol.GetConformer().GetPositions()), before)


# ---------------------------------------------------------------------------
# Highlighting and selection lifecycle
# ---------------------------------------------------------------------------


class TestHighlightAndSelection:
    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_highlight_is_added_under_the_dialog_s_own_actor_name(
        self, cls, make_dialog
    ):
        dlg, _mol, mw = make_dialog(cls)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0, 1})

        with patch(f"{_MIXIN}.pv", MagicMock()):
            MoveDialogMixin.show_atom_labels(dlg)

        assert mw.view_3d_manager.plotter.add_mesh.call_args.kwargs["name"] == (
            cls.HIGHLIGHT_ACTOR
        )

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_nothing_is_drawn_for_an_empty_selection(self, cls, make_dialog):
        dlg, _mol, mw = make_dialog(cls)
        dlg.group_atoms.clear()

        with patch(f"{_MIXIN}.pv", MagicMock()):
            MoveDialogMixin.show_atom_labels(dlg)

        mw.view_3d_manager.plotter.add_mesh.assert_not_called()

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_missing_3d_positions_bail_out_instead_of_raising(self, cls, make_dialog):
        """The 3D view can be torn down while the dialog is still open."""
        dlg, _mol, mw = make_dialog(cls)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0})
        mw.view_3d_manager.atom_positions_3d = None

        with patch(f"{_MIXIN}.pv", MagicMock()):
            MoveDialogMixin.show_atom_labels(dlg)

        mw.view_3d_manager.plotter.add_mesh.assert_not_called()

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_a_closed_plotter_is_not_drawn_into(self, cls, make_dialog):
        """The plotter can be gone by the time a queued redraw arrives."""
        dlg, _mol, mw = make_dialog(cls)
        dlg.group_atoms.clear()
        dlg.group_atoms.update({0})
        mw.view_3d_manager.plotter = None

        with patch(f"{_MIXIN}.pv", MagicMock()):
            MoveDialogMixin.show_atom_labels(dlg)

        assert dlg.highlight_actor is None

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_clear_atom_labels_drops_the_actor_by_name_and_by_handle(
        self, cls, make_dialog
    ):
        dlg, _mol, mw = make_dialog(cls)
        actor = MagicMock()
        dlg.highlight_actor = actor

        MoveDialogMixin.clear_atom_labels(dlg)

        removed = [
            call.args[0]
            for call in mw.view_3d_manager.plotter.remove_actor.call_args_list
        ]
        assert cls.HIGHLIGHT_ACTOR in removed
        assert actor in removed
        assert dlg.highlight_actor is None

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_clear_selection_empties_both_sets_and_the_drag_state(
        self, cls, make_dialog
    ):
        """group_atoms is an alias in one dialog and a wider BFS set in the other."""
        dlg, _mol, _mw = make_dialog(cls)
        dlg.selected_atoms.add(0)
        dlg.group_atoms.update({0, 1})
        dlg.is_dragging_group = True
        dlg.drag_start_pos = (10, 10)

        with (
            patch.object(type(dlg), "clear_atom_labels"),
            patch.object(type(dlg), "update_display"),
        ):
            dlg.clear_selection()

        assert not dlg.selected_atoms
        assert not dlg.group_atoms
        assert dlg.is_dragging_group is False
        assert dlg.drag_start_pos is None

    @pytest.mark.parametrize("cls", _dialog_classes())
    def test_reset_inputs_zero_every_field(self, cls, make_dialog):
        dlg, _mol, _mw = make_dialog(cls)
        for widget in (
            dlg.x_trans_input,
            dlg.y_trans_input,
            dlg.z_trans_input,
            dlg.x_rot_input,
            dlg.y_rot_input,
            dlg.z_rot_input,
        ):
            widget.setText("5.0")

        dlg.reset_translation_inputs()
        dlg.reset_rotation_inputs()

        assert dlg.x_trans_input.text() == "0.0"
        assert dlg.y_trans_input.text() == "0.0"
        assert dlg.z_trans_input.text() == "0.0"
        assert dlg.x_rot_input.text() == "0.0"
        assert dlg.y_rot_input.text() == "0.0"
        assert dlg.z_rot_input.text() == "0.0"
