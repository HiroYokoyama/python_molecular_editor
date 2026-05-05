"""
Extended tests for MoveGroupDialog — fills coverage gaps in:
  - __init__ preselected_atoms path
  - init_ui widget structure / window title
  - update_display with exactly 5 atoms (no ellipsis boundary)
  - apply_rotation with X and Y axes
  - show_atom_labels with mocked plotter
  - clear_atom_labels (with / without highlight actor)
  - eventFilter: early-return paths (None plotter, None mol, dbl-click)
"""

import os
import sys
import numpy as np
import pytest
from unittest.mock import MagicMock, patch, PropertyMock

from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtWidgets import QApplication

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ethane():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _five_atom_mol():
    """Linear C5 with 3D coords; 5 heavy atoms so ellipsis boundary is exactly at 5."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCCCC")
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def _make_main_window(mol):
    mw = MagicMock()
    mw._picking_consumed = False
    mw.view_3d_manager.atom_positions_3d = np.array(
        mol.GetConformer().GetPositions(), dtype=float
    )
    mw.view_3d_manager.plotter = MagicMock()
    return mw


@pytest.fixture(scope="session")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def make_dialog(qapp):
    created = []

    def _factory(mol=None):
        from moleditpy.ui.move_group_dialog import MoveGroupDialog

        _mol = mol if mol is not None else _ethane()
        mw = _make_main_window(_mol)
        with (
            patch.object(MoveGroupDialog, "show_atom_labels"),
            patch.object(MoveGroupDialog, "clear_atom_labels"),
        ):
            dlg = MoveGroupDialog(_mol, mw)
        entry = (dlg, _mol, mw)
        created.append(entry)
        return entry

    yield _factory

    for dlg, *_ in created:
        try:
            dlg.picking_enabled = False
            dlg.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# __init__ / preselected_atoms
# ---------------------------------------------------------------------------


class TestInit:
    def test_preselected_atoms_triggers_on_atom_picked(self, make_dialog, qapp):
        from moleditpy.ui.move_group_dialog import MoveGroupDialog

        mol = _ethane()
        mw = _make_main_window(mol)
        # update_display is called before init_ui creates selection_label, so patch it
        with (
            patch.object(MoveGroupDialog, "show_atom_labels"),
            patch.object(MoveGroupDialog, "clear_atom_labels"),
            patch.object(MoveGroupDialog, "update_display"),
        ):
            dlg = MoveGroupDialog(mol, mw, preselected_atoms=[0])
        # BFS from atom 0 should have selected all atoms in ethane
        assert len(dlg.group_atoms) == mol.GetNumAtoms()
        dlg.picking_enabled = False
        dlg.close()

    def test_no_preselected_atoms_leaves_group_empty(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert len(dlg.group_atoms) == 0

    def test_initial_drag_state(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert dlg.is_dragging_group is False
        assert dlg.drag_start_pos is None
        assert dlg.potential_drag is False

    def test_window_title(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert dlg.windowTitle() == "Move Group"


# ---------------------------------------------------------------------------
# init_ui widgets
# ---------------------------------------------------------------------------


class TestInitUI:
    def test_translation_inputs_default_zero(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert dlg.x_trans_input.text() == "0.0"
        assert dlg.y_trans_input.text() == "0.0"
        assert dlg.z_trans_input.text() == "0.0"

    def test_rotation_inputs_default_zero(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert dlg.x_rot_input.text() == "0.0"
        assert dlg.y_rot_input.text() == "0.0"
        assert dlg.z_rot_input.text() == "0.0"

    def test_selection_label_initial_text(self, make_dialog):
        dlg, _, _ = make_dialog()
        assert "No group" in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# update_display boundary
# ---------------------------------------------------------------------------


class TestUpdateDisplayBoundary:
    def test_exactly_5_atoms_no_ellipsis(self, make_dialog):
        """Exactly 5 selected atoms must NOT show '...'."""
        mol = _five_atom_mol()
        dlg, mol, _ = make_dialog(mol=mol)
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)  # selects all 5 heavy atoms

        dlg.update_display()
        text = dlg.selection_label.text()
        assert "5" in text
        assert "..." not in text

    def test_more_than_5_atoms_has_ellipsis(self, make_dialog):
        """8-atom ethane (with H) must show '...' after 5th."""
        dlg, _, _ = make_dialog()
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)
        dlg.update_display()
        assert "..." in dlg.selection_label.text()


# ---------------------------------------------------------------------------
# apply_rotation — X and Y axes
# ---------------------------------------------------------------------------


class TestApplyRotationAxes:
    def _pick_all(self, dlg):
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

    def test_90deg_x_rotation_around_centroid(self, make_dialog):
        """90° rotation around X maps [0,1,0] → [0,0,1] relative to centroid."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import Geometry

        mol = Chem.MolFromSmiles("[H][H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(0.0, 1.0, 0.0))
        conf.SetAtomPosition(1, Geometry.Point3D(0.0, -1.0, 0.0))

        dlg, mol, mw = make_dialog(mol=mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            [[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]
        )

        self._pick_all(dlg)
        dlg.x_rot_input.setText("90.0")
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()

        after = mol.GetConformer().GetPositions()
        # centroid=[0,0,0]; Rx 90°: [0,1,0]→[0,0,1], [0,-1,0]→[0,0,-1]
        assert after[0] == pytest.approx([0.0, 0.0, 1.0], abs=1e-5)
        assert after[1] == pytest.approx([0.0, 0.0, -1.0], abs=1e-5)

    def test_90deg_y_rotation_around_centroid(self, make_dialog):
        """90° rotation around Y maps [1,0,0] → [0,0,-1] relative to centroid."""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import Geometry

        mol = Chem.MolFromSmiles("[H][H]")
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        conf.SetAtomPosition(0, Geometry.Point3D(1.0, 0.0, 0.0))
        conf.SetAtomPosition(1, Geometry.Point3D(-1.0, 0.0, 0.0))

        dlg, mol, mw = make_dialog(mol=mol)
        mw.view_3d_manager.atom_positions_3d = np.array(
            [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
        )

        self._pick_all(dlg)
        dlg.y_rot_input.setText("90.0")
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()

        after = mol.GetConformer().GetPositions()
        # centroid=[0,0,0]; Ry 90°: [1,0,0]→[0,0,-1], [-1,0,0]→[0,0,1]
        assert after[0] == pytest.approx([0.0, 0.0, -1.0], abs=1e-5)
        assert after[1] == pytest.approx([0.0, 0.0, 1.0], abs=1e-5)

    def test_combined_rotation_pushes_undo(self, make_dialog):
        dlg, _, mw = make_dialog()
        self._pick_all(dlg)
        dlg.x_rot_input.setText("30.0")
        dlg.y_rot_input.setText("45.0")
        dlg.z_rot_input.setText("60.0")
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.apply_rotation()
        mw.edit_actions_manager.push_undo_state.assert_called()


# ---------------------------------------------------------------------------
# show_atom_labels / clear_atom_labels
# ---------------------------------------------------------------------------


class TestAtomLabels:
    def _pick_all(self, dlg):
        with (
            patch.object(type(dlg), "show_atom_labels"),
            patch.object(type(dlg), "clear_atom_labels"),
        ):
            dlg.on_atom_picked(0)

    def test_show_atom_labels_calls_plotter_add_mesh(self, make_dialog):
        import pyvista as pv

        dlg, mol, mw = make_dialog()
        self._pick_all(dlg)
        mw.view_3d_manager.atom_positions_3d = np.array(
            mol.GetConformer().GetPositions(), dtype=float
        )
        mock_plotter = MagicMock()
        mw.view_3d_manager.plotter = mock_plotter
        # pv.Sphere is not in the mock pyvista; create=True allows adding it
        with patch("pyvista.Sphere", create=True, return_value=MagicMock()):
            dlg.show_atom_labels()
        mock_plotter.add_mesh.assert_called()
        mock_plotter.render.assert_called()

    def test_show_atom_labels_no_group_does_nothing(self, make_dialog):
        dlg, _, mw = make_dialog()
        dlg.group_atoms.clear()
        mw.view_3d_manager.plotter = MagicMock()
        dlg.show_atom_labels()
        mw.view_3d_manager.plotter.add_mesh.assert_not_called()

    def test_clear_atom_labels_removes_highlight_actor(self, make_dialog):
        dlg, _, mw = make_dialog()
        mock_actor = MagicMock()
        dlg.highlight_actor = mock_actor
        mock_plotter = MagicMock()
        mw.view_3d_manager.plotter = mock_plotter
        # Patch base class clear_atom_labels to avoid label list complexity
        with patch(
            "moleditpy.ui.dialog_3d_picking_mixin.Dialog3DPickingMixin.clear_atom_labels"
        ):
            dlg.clear_atom_labels()
        mock_plotter.remove_actor.assert_called()
        assert dlg.highlight_actor is None

    def test_clear_atom_labels_none_plotter_does_not_raise(self, make_dialog):
        dlg, _, mw = make_dialog()
        mw.view_3d_manager.plotter = None
        with patch(
            "moleditpy.ui.dialog_3d_picking_mixin.Dialog3DPickingMixin.clear_atom_labels"
        ):
            dlg.clear_atom_labels()  # should not raise


# ---------------------------------------------------------------------------
# eventFilter — early-return paths
# ---------------------------------------------------------------------------


class TestEventFilter:
    def _make_event(self, event_type):
        ev = MagicMock()
        ev.type.return_value = event_type
        return ev

    def test_returns_false_when_plotter_is_none(self, make_dialog):
        dlg, _, mw = make_dialog()
        mw.view_3d_manager.plotter = None
        ev = self._make_event(QEvent.Type.MouseButtonPress)
        result = dlg.eventFilter(MagicMock(), ev)
        assert result is False

    def test_returns_false_when_mol_is_none(self, make_dialog):
        dlg, _, mw = make_dialog()
        dlg.mol = None
        ev = self._make_event(QEvent.Type.MouseButtonPress)
        result = dlg.eventFilter(MagicMock(), ev)
        assert result is False

    def test_double_click_resets_state_and_returns_false(self, make_dialog):
        dlg, _, mw = make_dialog()
        dlg.is_dragging_group = True
        dlg.potential_drag = True

        ev = MagicMock()
        ev.type.return_value = QEvent.Type.MouseButtonDblClick

        # obj must match plotter.interactor for the filter to engage
        plotter = MagicMock()
        mw.view_3d_manager.plotter = plotter
        result = dlg.eventFilter(plotter.interactor, ev)
        assert result is False
        assert dlg.is_dragging_group is False
        assert dlg.potential_drag is False

    def test_non_interactor_obj_delegates_to_super(self, make_dialog):
        """Events on objects other than the plotter interactor use base behaviour."""
        dlg, _, mw = make_dialog()
        ev = MagicMock()
        ev.type.return_value = QEvent.Type.MouseButtonPress
        # Pass an unrelated object — should NOT be the plotter.interactor
        result = dlg.eventFilter(MagicMock(), ev)
        # Super's eventFilter returns False for unrecognised objects
        assert result is False
