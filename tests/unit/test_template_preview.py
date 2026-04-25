"""
Tests for template preview components.

Covers:
  TemplatePreviewItem (ui/template_preview_item.py):
    - __init__ defaults
    - set_geometry / set_user_template_geometry
    - boundingRect: empty polygon, regular polygon, user-template points
    - paint: None painter early-return
    - paint_regular_template: non-aromatic, aromatic, empty polygon
    - paint_user_template: no points, single/double/triple bonds, non-C atoms

  TemplatePreviewView (ui/template_preview_view.py):
    - __init__ defaults
    - set_template_data
    - resizeEvent: with / without original_scene_rect
    - refit_view: valid rect, empty rect, no rect
    - showEvent: with / without rect
    - redraw_with_current_size: no data, with data + parent_dialog

  UserTemplateDialog (ui/user_template_dialog.py):
    - load_template_file: valid JSON, missing file, invalid JSON
    - save_template_file: success, OSError
    - draw_template_preview: no atoms, single/double/triple bonds, non-C atom
    - convert_structure_to_template: returns correct dict
    - cleanup_template_mode: resets state, handles missing ui_manager
    - select_template: sets selection, enables delete button
    - use_template: calls _activate_template_mode
    - fit_preview_view_safely: valid rect, empty rect
    - refit_all_previews: iterates children, handles errors
    - delete_selected_template: no selection → early return
"""

import json
import os
import pytest
from unittest.mock import MagicMock, patch, call
from PyQt6.QtCore import QPointF, QRectF, QSize
from PyQt6.QtGui import QImage, QPainter, QResizeEvent, QShowEvent
from PyQt6.QtWidgets import QGraphicsScene

from moleditpy.ui.template_preview_item import TemplatePreviewItem
from moleditpy.ui.template_preview_view import TemplatePreviewView
from moleditpy.ui.user_template_dialog import UserTemplateDialog


# ---------------------------------------------------------------------------
# Helpers — QPainter on an offscreen image
# ---------------------------------------------------------------------------

def _painter() -> tuple:
    """Return (QImage, QPainter) ready for drawing; caller must call painter.end()."""
    img = QImage(200, 200, QImage.Format.Format_ARGB32)
    p = QPainter(img)
    return img, p


def _points(n: int = 6) -> list:
    """Return a regular polygon of n QPointFs."""
    import math
    return [
        QPointF(50 + 30 * math.cos(2 * math.pi * i / n),
                50 + 30 * math.sin(2 * math.pi * i / n))
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# TemplatePreviewItem — __init__ / geometry setters
# ---------------------------------------------------------------------------

def test_preview_item_init_defaults(app):
    item = TemplatePreviewItem()
    assert item.is_aromatic is False
    assert item.is_user_template is False
    assert item.user_template_points == []
    assert item.user_template_bonds == []
    assert item.user_template_atoms == []


def test_set_geometry_updates_polygon(app):
    item = TemplatePreviewItem()
    pts = _points(6)
    item.set_geometry(pts, is_aromatic=False)
    assert not item.polygon.isEmpty()
    assert item.is_aromatic is False
    assert item.is_user_template is False


def test_set_geometry_aromatic_flag(app):
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    assert item.is_aromatic is True


def test_set_user_template_geometry(app):
    item = TemplatePreviewItem()
    pts = _points(3)
    bonds = [(0, 1, 1), (1, 2, 2)]
    atoms = [{"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    assert item.is_user_template is True
    assert item.user_template_points == pts
    assert item.user_template_bonds == bonds
    assert item.user_template_atoms == atoms
    assert item.is_aromatic is False


# ---------------------------------------------------------------------------
# TemplatePreviewItem — boundingRect
# ---------------------------------------------------------------------------

def test_bounding_rect_empty_polygon(app):
    item = TemplatePreviewItem()
    rect = item.boundingRect()
    assert isinstance(rect, QRectF)


def test_bounding_rect_regular_polygon(app):
    item = TemplatePreviewItem()
    item.set_geometry(_points(6))
    rect = item.boundingRect()
    assert rect.width() > 0


def test_bounding_rect_user_template_with_points(app):
    item = TemplatePreviewItem()
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(100, 0), QPointF(50, 80)],
        [], []
    )
    rect = item.boundingRect()
    assert rect.width() > 0
    assert rect.height() > 0


def test_bounding_rect_user_template_no_points(app):
    item = TemplatePreviewItem()
    item.is_user_template = True
    item.user_template_points = []
    rect = item.boundingRect()
    assert isinstance(rect, QRectF)


# ---------------------------------------------------------------------------
# TemplatePreviewItem — paint (None painter guard)
# ---------------------------------------------------------------------------

def test_paint_none_painter_returns_early(app):
    item = TemplatePreviewItem()
    item.set_geometry(_points(6))
    # Should not raise
    item.paint(None, None, None)


# ---------------------------------------------------------------------------
# TemplatePreviewItem — paint_regular_template
# ---------------------------------------------------------------------------

def test_paint_regular_template_non_aromatic(app):
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=False)
    _, p = _painter()
    item.paint_regular_template(p)
    p.end()


def test_paint_regular_template_aromatic(app):
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    _, p = _painter()
    item.paint_regular_template(p)
    p.end()


def test_paint_regular_template_empty_polygon(app):
    item = TemplatePreviewItem()
    # polygon is empty by default — should not raise
    _, p = _painter()
    item.paint_regular_template(p)
    p.end()


# ---------------------------------------------------------------------------
# TemplatePreviewItem — paint_user_template
# ---------------------------------------------------------------------------

def test_paint_user_template_no_points_returns_early(app):
    item = TemplatePreviewItem()
    item.user_template_points = []
    _, p = _painter()
    item.paint_user_template(p)   # should not raise
    p.end()


def test_paint_user_template_single_bond(app):
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 1)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


def test_paint_user_template_double_bond(app):
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 2)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


def test_paint_user_template_triple_bond(app):
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 3)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


def test_paint_user_template_non_carbon_atom(app):
    item = TemplatePreviewItem()
    pts = [QPointF(50, 50), QPointF(100, 50)]
    bonds = [(0, 1, 1)]
    atoms = [{"symbol": "N"}, {"symbol": "O"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


def test_paint_user_template_bond_info_2_elements(app):
    """bond_info with only 2 elements (no order) defaults to single bond."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1)]   # no order field
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


def test_paint_user_template_out_of_range_indices(app):
    """Bond indices beyond point list should be skipped without error."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10)]
    bonds = [(0, 5, 1)]   # index 5 out of range
    atoms = [{"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_user_template(p)
    p.end()


# ---------------------------------------------------------------------------
# TemplatePreviewView — __init__ / set_template_data
# ---------------------------------------------------------------------------

def test_preview_view_init_defaults(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    assert view.original_scene_rect is None
    assert view.template_data is None
    assert view.parent_dialog is None


def test_preview_view_set_template_data(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    data = {"name": "test"}
    parent = MagicMock()
    view.set_template_data(data, parent)
    assert view.template_data == data
    assert view.parent_dialog is parent


# ---------------------------------------------------------------------------
# TemplatePreviewView — resizeEvent
# ---------------------------------------------------------------------------

def test_preview_view_resize_no_scene_rect_no_timer(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    ev = QResizeEvent(QSize(200, 200), QSize(100, 100))
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.resizeEvent(ev)
    mock_timer.assert_not_called()


def test_preview_view_resize_with_scene_rect_schedules_timer(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF(0, 0, 100, 100)
    ev = QResizeEvent(QSize(200, 200), QSize(100, 100))
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.resizeEvent(ev)
    mock_timer.assert_called_once()


# ---------------------------------------------------------------------------
# TemplatePreviewView — refit_view
# ---------------------------------------------------------------------------

def test_refit_view_no_rect_is_noop(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    view.refit_view()   # should not raise


def test_refit_view_empty_rect_is_noop(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF()   # empty
    view.refit_view()   # should not raise


def test_refit_view_valid_rect(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF(0, 0, 100, 100)
    view.refit_view()   # should not raise


# ---------------------------------------------------------------------------
# TemplatePreviewView — showEvent
# ---------------------------------------------------------------------------

def test_show_event_no_rect_no_timer(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    ev = QShowEvent()
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.showEvent(ev)
    mock_timer.assert_not_called()


def test_show_event_with_rect_schedules_timer(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF(0, 0, 100, 100)
    ev = QShowEvent()
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.showEvent(ev)
    mock_timer.assert_called_once()


# ---------------------------------------------------------------------------
# TemplatePreviewView — redraw_with_current_size
# ---------------------------------------------------------------------------

def test_redraw_no_data_returns_early(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.template_data = None
    view.parent_dialog = None
    view.redraw_with_current_size()   # should not raise


def test_redraw_with_data_calls_draw_template_preview(app):
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    data = {"atoms": [], "bonds": []}
    parent = MagicMock()
    view.template_data = data
    view.parent_dialog = parent
    view.redraw_with_current_size()
    parent.draw_template_preview.assert_called_once()


# ---------------------------------------------------------------------------
# Helpers — UserTemplateDialog
# ---------------------------------------------------------------------------

def _make_mw(tmp_path):
    mw = MagicMock()
    mw.init_manager.settings_dir = str(tmp_path)
    return mw


def _make_dlg(app, tmp_path):
    mw = _make_mw(tmp_path)
    with patch.object(UserTemplateDialog, "load_user_templates"):
        dlg = UserTemplateDialog(mw)
    return dlg, mw


# ---------------------------------------------------------------------------
# UserTemplateDialog — load_template_file
# ---------------------------------------------------------------------------

def test_load_template_file_valid(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    data = {"name": "benzene", "atoms": [], "bonds": []}
    fp = tmp_path / "benzene.pmetmplt"
    fp.write_text(json.dumps(data), encoding="utf-8")
    result = dlg.load_template_file(str(fp))
    assert result == data


def test_load_template_file_missing(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    result = dlg.load_template_file(str(tmp_path / "nonexistent.pmetmplt"))
    assert result is None


def test_load_template_file_invalid_json(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    fp = tmp_path / "bad.pmetmplt"
    fp.write_text("NOT JSON", encoding="utf-8")
    result = dlg.load_template_file(str(fp))
    assert result is None


# ---------------------------------------------------------------------------
# UserTemplateDialog — save_template_file
# ---------------------------------------------------------------------------

def test_save_template_file_success(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    data = {"name": "ethane"}
    fp = tmp_path / "ethane.pmetmplt"
    result = dlg.save_template_file(str(fp), data)
    assert result is True
    assert json.loads(fp.read_text()) == data


def test_save_template_file_oserror(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    result = dlg.save_template_file("/nonexistent_dir/bad.pmetmplt", {})
    assert result is False


# ---------------------------------------------------------------------------
# UserTemplateDialog — draw_template_preview
# ---------------------------------------------------------------------------

def _scene():
    return QGraphicsScene()


def test_draw_template_preview_no_atoms_adds_placeholder(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    dlg.draw_template_preview(sc, {"atoms": [], "bonds": []})
    assert sc.items()   # placeholder text item


def test_draw_template_preview_single_bond(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    template = {
        "atoms": [
            {"id": 0, "symbol": "C", "x": 0.0, "y": 0.0},
            {"id": 1, "symbol": "C", "x": 50.0, "y": 0.0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 1}],
    }
    dlg.draw_template_preview(sc, template)
    assert sc.items()


def test_draw_template_preview_double_bond(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    template = {
        "atoms": [
            {"id": 0, "symbol": "C", "x": 0.0, "y": 0.0},
            {"id": 1, "symbol": "C", "x": 50.0, "y": 0.0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 2}],
    }
    dlg.draw_template_preview(sc, template)
    assert sc.items()


def test_draw_template_preview_triple_bond(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    template = {
        "atoms": [
            {"id": 0, "symbol": "C", "x": 0.0, "y": 0.0},
            {"id": 1, "symbol": "C", "x": 50.0, "y": 0.0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 3}],
    }
    dlg.draw_template_preview(sc, template)
    assert sc.items()


def test_draw_template_preview_non_carbon_atom(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    template = {
        "atoms": [
            {"id": 0, "symbol": "N", "x": 0.0, "y": 0.0},
            {"id": 1, "symbol": "O", "x": 50.0, "y": 0.0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 1}],
    }
    dlg.draw_template_preview(sc, template)
    assert sc.items()


def test_draw_template_preview_zero_size_mol(app, tmp_path):
    """Single atom (mol_size=0) falls back to scale_factor=1.0 without crash."""
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    template = {
        "atoms": [{"id": 0, "symbol": "C", "x": 0.0, "y": 0.0}],
        "bonds": [],
    }
    dlg.draw_template_preview(sc, template)


# ---------------------------------------------------------------------------
# UserTemplateDialog — convert_structure_to_template
# ---------------------------------------------------------------------------

def test_convert_structure_to_template(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    mw.state_manager.data.atoms = {
        1: {"symbol": "C", "pos": (10.0, 20.0), "charge": 0, "radical": 0},
        2: {"symbol": "N", "pos": (30.0, 20.0), "charge": 1, "radical": 0},
    }
    mw.state_manager.data.bonds = {
        (1, 2): {"order": 2, "stereo": 0},
    }
    result = dlg.convert_structure_to_template("TestMol")
    assert result["name"] == "TestMol"
    assert len(result["atoms"]) == 2
    assert len(result["bonds"]) == 1
    assert result["bonds"][0]["order"] == 2
    assert result["format"] == "PME Template"


# ---------------------------------------------------------------------------
# UserTemplateDialog — cleanup_template_mode
# ---------------------------------------------------------------------------

def test_cleanup_template_mode_resets_selected(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    dlg.selected_template = {"name": "foo"}
    dlg.cleanup_template_mode()
    assert dlg.selected_template is None
    assert dlg.delete_button.isEnabled() is False


def test_cleanup_template_mode_calls_set_mode(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    mw.ui_manager = MagicMock()
    dlg.cleanup_template_mode()
    mw.ui_manager.set_mode_and_update_toolbar.assert_called_once_with("atom_C")


def test_cleanup_template_mode_no_ui_manager_no_crash(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    del mw.ui_manager
    dlg.cleanup_template_mode()   # should not raise


def test_cleanup_template_mode_scene_reset(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    scene = MagicMock()
    scene.views.return_value = [MagicMock()]
    mw.init_manager.scene = scene
    dlg.cleanup_template_mode()
    assert scene.mode == "atom_C"
    assert scene.user_template_data is None


# ---------------------------------------------------------------------------
# UserTemplateDialog — select_template
# ---------------------------------------------------------------------------

def test_select_template_sets_selected_and_enables_delete(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    data = {"name": "benzene", "atoms": [], "bonds": []}
    widget = MagicMock()
    with patch.object(dlg, "_activate_template_mode"):
        dlg.select_template(data, widget)
    assert dlg.selected_template == data
    assert dlg.delete_button.isEnabled() is True


def test_select_template_activates_template_mode(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    data = {"name": "ethane", "atoms": [], "bonds": []}
    widget = MagicMock()
    with patch.object(dlg, "_activate_template_mode") as mock_act:
        dlg.select_template(data, widget)
    mock_act.assert_called_once_with(data)


# ---------------------------------------------------------------------------
# UserTemplateDialog — use_template
# ---------------------------------------------------------------------------

def test_use_template_calls_activate(app, tmp_path):
    dlg, mw = _make_dlg(app, tmp_path)
    data = {"name": "methane", "atoms": [], "bonds": []}
    with patch.object(dlg, "_activate_template_mode") as mock_act:
        dlg.use_template(data)
    mock_act.assert_called_once_with(data)
    assert dlg.selected_template == data


# ---------------------------------------------------------------------------
# UserTemplateDialog — fit_preview_view_safely
# ---------------------------------------------------------------------------

def test_fit_preview_view_safely_valid_rect(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    view = MagicMock()
    rect = QRectF(0, 0, 100, 100)
    dlg.fit_preview_view_safely(view, rect)
    view.fitInView.assert_called_once()


def test_fit_preview_view_safely_empty_rect(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    view = MagicMock()
    dlg.fit_preview_view_safely(view, QRectF())   # empty → no-op
    view.fitInView.assert_not_called()


def test_fit_preview_view_safely_none_view(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    # None view should not raise (guard: `if view and ...`)
    dlg.fit_preview_view_safely(None, QRectF(0, 0, 10, 10))


# ---------------------------------------------------------------------------
# UserTemplateDialog — refit_all_previews
# ---------------------------------------------------------------------------

def test_refit_all_previews_no_items_no_crash(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    dlg.refit_all_previews()   # empty grid — should not raise


def test_refit_all_previews_calls_redraw(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    preview_view = MagicMock(spec=TemplatePreviewView)
    preview_view.redraw_with_current_size = MagicMock()

    child_widget = MagicMock()
    child_widget.findChildren.return_value = [preview_view]

    layout_item = MagicMock()
    layout_item.widget.return_value = child_widget

    with patch.object(dlg.template_layout, "count", return_value=1), \
         patch.object(dlg.template_layout, "itemAt", return_value=layout_item):
        dlg.refit_all_previews()

    preview_view.redraw_with_current_size.assert_called_once()


# ---------------------------------------------------------------------------
# UserTemplateDialog — delete_selected_template (no selection)
# ---------------------------------------------------------------------------

def test_delete_selected_no_selection_returns_early(app, tmp_path):
    dlg, _ = _make_dlg(app, tmp_path)
    dlg.selected_template = None
    with patch("moleditpy.ui.user_template_dialog.QMessageBox.question") as mock_q:
        dlg.delete_selected_template()
    mock_q.assert_not_called()
