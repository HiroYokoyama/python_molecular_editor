"""
Tests for template preview components.

Covers:
  TemplatePreviewItem (ui/template_preview_item.py):
    - __init__ defaults
    - set_geometry / set_user_template_geometry
    - boundingRect: empty polygon, regular polygon, user-template points
    - paint: None painter early-return
    - paint_ghost: ring and user templates, single/double/triple bonds, non-C atoms
    - set_chain_geometry / paint_chain: alkyl chain ghost + repeat-count label

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
import pytest
from unittest.mock import MagicMock, patch
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
        QPointF(
            50 + 30 * math.cos(2 * math.pi * i / n),
            50 + 30 * math.sin(2 * math.pi * i / n),
        )
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# TemplatePreviewItem — __init__ / geometry setters
# ---------------------------------------------------------------------------


def test_preview_item_init_defaults(app):
    """TemplatePreviewItem initialises with no ghost molecule and nothing to draw."""
    item = TemplatePreviewItem()
    assert item.ghost_atoms == []
    assert item.ghost_bonds == []
    assert item.is_chain is False
    assert item.boundingRect().isNull()


def test_set_geometry_builds_a_ring(app):
    """set_geometry turns the polygon into single-bonded ghost atoms."""
    item = TemplatePreviewItem()
    pts = _points(6)
    item.set_geometry(pts, is_aromatic=False)
    assert [atom.pos() for atom in item.ghost_atoms] == pts
    assert [bond.order for bond in item.ghost_bonds] == [1] * 6


def test_set_geometry_aromatic_alternates_bonds(app):
    """An aromatic ring previews as the Kekulé structure that gets placed."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    assert [bond.order for bond in item.ghost_bonds] == [2, 1, 2, 1, 2, 1]


def test_set_geometry_takes_given_bond_orders(app):
    """The caller can pass the orders the placement will use."""
    item = TemplatePreviewItem()
    orders = [(i, (i + 1) % 6, 1 if i % 2 == 0 else 2) for i in range(6)]
    item.set_geometry(_points(6), is_aromatic=True, bonds_info=orders)
    assert [bond.order for bond in item.ghost_bonds] == [1, 2, 1, 2, 1, 2]


def test_set_user_template_geometry(app):
    """set_user_template_geometry builds the template's atoms and bonds."""
    item = TemplatePreviewItem()
    pts = _points(3)
    bonds = [(0, 1, 1), (1, 2, 2)]
    atoms = [{"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    assert [atom.symbol for atom in item.ghost_atoms] == ["C", "N", "O"]
    assert [bond.order for bond in item.ghost_bonds] == [1, 2]
    assert item.is_chain is False


# ---------------------------------------------------------------------------
# TemplatePreviewItem — boundingRect
# ---------------------------------------------------------------------------


def test_bounding_rect_empty_polygon(app):
    """boundingRect returns a QRectF when the polygon is empty."""
    item = TemplatePreviewItem()
    rect = item.boundingRect()
    assert isinstance(rect, QRectF)


def test_bounding_rect_regular_polygon(app):
    """boundingRect returns a positive-width rect for a non-empty polygon."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6))
    rect = item.boundingRect()
    assert rect.width() > 0


def test_bounding_rect_user_template_with_points(app):
    """boundingRect returns positive width and height for a user template with points."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(100, 0), QPointF(50, 80)], [], []
    )
    rect = item.boundingRect()
    assert rect.width() > 0
    assert rect.height() > 0


def test_bounding_rect_user_template_no_points(app):
    """boundingRect returns a QRectF for a user template with no points."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry([], [], [])
    rect = item.boundingRect()
    assert isinstance(rect, QRectF)


# ---------------------------------------------------------------------------
# TemplatePreviewItem — paint (None painter guard)
# ---------------------------------------------------------------------------


def test_paint_none_painter_returns_early(app):
    """paint with a None painter returns without raising."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6))
    # Should not raise
    item.paint(None, None, None)


# ---------------------------------------------------------------------------
# TemplatePreviewItem — painting a ring template
# ---------------------------------------------------------------------------


def test_paint_ring_template_non_aromatic(app):
    """A plain ring template paints without raising."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=False)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_ring_template_aromatic(app):
    """An aromatic ring template paints without raising."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_without_geometry_draws_nothing(app):
    """A preview with no ghost molecule paints nothing and does not raise."""
    item = TemplatePreviewItem()
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


# ---------------------------------------------------------------------------
# TemplatePreviewItem — alkyl chain preview
# ---------------------------------------------------------------------------


def _blank_painter() -> tuple:
    """Return (white QImage, QPainter); caller must call painter.end()."""
    img = QImage(200, 200, QImage.Format.Format_ARGB32)
    img.fill(0xFFFFFFFF)
    return img, QPainter(img)


def _painted_pixel_count(img: QImage) -> int:
    return sum(
        1
        for y in range(img.height())
        for x in range(img.width())
        if img.pixel(x, y) != 0xFFFFFFFF
    )


def test_set_chain_geometry_stores_points_and_label(app):
    """set_chain_geometry switches the item into chain mode and clears other modes."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry([QPointF(0, 0)], [], [{"symbol": "C"}])
    pts = [QPointF(0, 0), QPointF(50, 20)]
    item.set_chain_geometry(pts, "n = 1", QPointF(60, 10))

    assert item.is_chain is True
    assert item.ghost_atoms == []
    assert item.chain_points == pts
    assert item.chain_label == "n = 1"
    assert item.chain_label_pos == QPointF(60, 10)


def test_set_geometry_leaves_chain_mode(app):
    """A ring template preview after a chain preview no longer paints as a chain."""
    item = TemplatePreviewItem()
    item.set_chain_geometry([QPointF(0, 0), QPointF(50, 0)], "n = 1", QPointF(60, 0))
    item.set_geometry(_points(6))
    assert item.is_chain is False


def test_bounding_rect_covers_chain_and_label(app):
    """boundingRect spans both the zigzag and the label drawn beside the cursor."""
    item = TemplatePreviewItem()
    item.set_chain_geometry([QPointF(0, 0), QPointF(60, 30)], "n = 1", QPointF(76, 18))
    rect = item.boundingRect()
    assert rect.contains(QPointF(0, 0))
    assert rect.contains(QPointF(76, 18))


def test_paint_chain_draws_the_zigzag_and_label(app):
    """paint_chain marks the canvas for both the chain and its label."""
    item = TemplatePreviewItem()
    item.set_chain_geometry(
        [QPointF(10, 100), QPointF(70, 60), QPointF(130, 100)],
        "n = 2",
        QPointF(140, 40),
    )
    img, p = _blank_painter()
    item.paint(p, None, None)
    p.end()
    assert _painted_pixel_count(img) > 0


def test_paint_chain_needs_at_least_one_bond(app):
    """A chain preview with a single vertex paints nothing."""
    item = TemplatePreviewItem()
    item.set_chain_geometry([QPointF(10, 10)], "n = 0", QPointF(20, 20))
    img, p = _blank_painter()
    item.paint(p, None, None)
    p.end()
    assert _painted_pixel_count(img) == 0


# ---------------------------------------------------------------------------
# TemplatePreviewItem — painting a user template
# ---------------------------------------------------------------------------


def test_paint_user_template_no_points_returns_early(app):
    """A user template with no points paints nothing and does not raise."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry([], [], [])
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_single_bond(app):
    """paint_ghost renders a single bond without raising."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 1)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_double_bond(app):
    """paint_ghost renders a double bond without raising."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 2)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_triple_bond(app):
    """paint_ghost renders a triple bond without raising."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1, 3)]
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_non_carbon_atom(app):
    """paint_ghost renders non-carbon atom labels without raising."""
    item = TemplatePreviewItem()
    pts = [QPointF(50, 50), QPointF(100, 50)]
    bonds = [(0, 1, 1)]
    atoms = [{"symbol": "N"}, {"symbol": "O"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_bond_info_2_elements(app):
    """bond_info with only 2 elements (no order) defaults to single bond."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10), QPointF(60, 10)]
    bonds = [(0, 1)]  # no order field
    atoms = [{"symbol": "C"}, {"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_paint_user_template_out_of_range_indices(app):
    """Bond indices beyond point list should be skipped without error."""
    item = TemplatePreviewItem()
    pts = [QPointF(10, 10)]
    bonds = [(0, 5, 1)]  # index 5 out of range
    atoms = [{"symbol": "C"}]
    item.set_user_template_geometry(pts, bonds, atoms)
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


# ---------------------------------------------------------------------------
# TemplatePreviewItem — ghost molecule built from the real AtomItem/BondItem
# ---------------------------------------------------------------------------


def test_ghost_items_mirror_the_template(app):
    """The preview builds one AtomItem per point and one BondItem per bond."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(50, 0), QPointF(75, 43)],
        [(0, 1, 2, 0), (1, 2, 1, 0)],
        [{"symbol": "N"}, {"symbol": "C"}, {"symbol": "C"}],
    )
    assert [a.symbol for a in item.ghost_atoms] == ["N", "C", "C"]
    assert [b.order for b in item.ghost_bonds] == [2, 1]
    # The bonds have to be registered on the atoms or carbons stay visible
    assert len(item.ghost_atoms[1].bonds) == 2


def test_ghost_atoms_get_implicit_hydrogens(app):
    """Implicit H counts come from RDKit, so O of an alcohol previews as OH."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(50, 0)],
        [(0, 1, 1, 0)],
        [{"symbol": "C"}, {"symbol": "O"}],
    )
    assert item.ghost_atoms[0].implicit_h_count == 3
    assert item.ghost_atoms[1].implicit_h_count == 1


def test_ghost_ring_bonds_carry_ring_center(app):
    """Ring double bonds need is_in_ring/ring_center to draw the inner short line."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    assert all(b.is_in_ring for b in item.ghost_bonds)
    assert all(b.ring_center is not None for b in item.ghost_bonds)
    # Kekulé alternation, matching what add_molecule_fragment creates
    assert [b.order for b in item.ghost_bonds] == [2, 1, 2, 1, 2, 1]


def test_impossible_valence_still_previews(app):
    """The preview draws what the template says — the valence is never checked."""
    item = TemplatePreviewItem()
    pts = [QPointF(0, 0)] + _points(5)
    bonds = [(0, i, 1, 0) for i in range(1, 6)]
    atoms = [{"symbol": "C"}] * 6
    item.set_user_template_geometry(pts, bonds, atoms)
    assert len(item.ghost_bonds) == 5
    # The five-bonded carbon simply gets no hydrogens; the rest keep theirs
    assert item.ghost_atoms[0].implicit_h_count == 0
    assert item.ghost_atoms[1].implicit_h_count == 3
    _, p = _painter()
    item.paint_ghost(p)
    p.end()


def test_geometry_change_is_announced_after_the_ghost_exists(app):
    """Regression: the scene cached this item's empty boundingRect during the
    rebuild (find_atom_near queries it), so Qt culled the preview and no ghost
    was painted in the editor. The rebuild must announce the geometry again
    once the ghost items are in place."""
    item = TemplatePreviewItem()
    seen = []
    item.prepareGeometryChange = lambda: seen.append(len(item.ghost_atoms))

    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(50, 0)],
        [(0, 1, 1, 0)],
        [{"symbol": "N"}, {"symbol": "C"}],
    )

    assert seen[-1] == 2


def test_moving_the_cursor_reuses_the_ghost_items(app):
    """A cursor move must not rebuild the ghost; rebuilding every move made big
    templates lag."""
    item = TemplatePreviewItem()
    bonds = [(0, 1, 1, 0)]
    atoms = [{"symbol": "N"}, {"symbol": "C"}]
    item.set_user_template_geometry([QPointF(0, 0), QPointF(50, 0)], bonds, atoms)
    first_atoms = list(item.ghost_atoms)

    item.set_user_template_geometry([QPointF(90, 30), QPointF(140, 30)], bonds, atoms)

    assert item.ghost_atoms == first_atoms
    assert item.ghost_atoms[0].pos() == QPointF(90, 30)
    assert item.ghost_bonds[0].pos() == QPointF(90, 30)


def test_changed_template_rebuilds_the_ghost(app):
    """A different template must not reuse the previous ghost items."""
    item = TemplatePreviewItem()
    points = [QPointF(0, 0), QPointF(50, 0)]
    item.set_user_template_geometry(
        points, [(0, 1, 1, 0)], [{"symbol": "N"}, {"symbol": "C"}]
    )
    first_atoms = list(item.ghost_atoms)

    item.set_user_template_geometry(
        points, [(0, 1, 2, 0)], [{"symbol": "O"}, {"symbol": "C"}]
    )

    assert item.ghost_atoms != first_atoms
    assert item.ghost_atoms[0].symbol == "O"
    assert item.ghost_bonds[0].order == 2


def test_ring_centres_follow_a_moved_ghost(app):
    """The inner line of a ring double bond is placed from the ring centre, so the
    centre has to move with the preview."""
    item = TemplatePreviewItem()
    item.set_geometry(_points(6), is_aromatic=True)
    before = item.ghost_bonds[0].ring_center

    item.set_geometry([p + QPointF(200, 0) for p in _points(6)], is_aromatic=True)
    after = item.ghost_bonds[0].ring_center

    assert after[0] == pytest.approx(before[0] + 200)
    assert after[1] == pytest.approx(before[1])


def test_first_atom_is_marked_for_user_templates_only(app):
    """The attachment atom is highlighted for user templates, not for rings."""
    item = TemplatePreviewItem()
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(50, 0)],
        [(0, 1, 1, 0)],
        [{"symbol": "N"}, {"symbol": "C"}],
    )
    assert item.mark_first_atom is True
    item.set_geometry(_points(6))
    assert item.mark_first_atom is False


def test_ghosts_stay_out_of_the_editor_scene(app):
    """Ghost atoms must never enter the editor scene or hit tests would find them."""
    scene = QGraphicsScene()
    item = TemplatePreviewItem()
    scene.addItem(item)
    item.set_user_template_geometry(
        [QPointF(0, 0), QPointF(50, 0)],
        [(0, 1, 1, 0)],
        [{"symbol": "C"}, {"symbol": "C"}],
    )
    assert item.ghost_atoms
    assert all(a.scene() is item.ghost_scene for a in item.ghost_atoms)
    assert scene.items() == [item]


# ---------------------------------------------------------------------------
# TemplatePreviewView — __init__ / set_template_data
# ---------------------------------------------------------------------------


def test_preview_view_init_defaults(app):
    """TemplatePreviewView initialises with original_scene_rect, template_data, and parent_dialog all None."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    assert view.original_scene_rect is None
    assert view.template_data is None
    assert view.parent_dialog is None


def test_preview_view_set_template_data(app):
    """set_template_data stores template_data and parent_dialog on the view."""
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
    """resizeEvent does not schedule a timer when original_scene_rect is None."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    ev = QResizeEvent(QSize(200, 200), QSize(100, 100))
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.resizeEvent(ev)
    mock_timer.assert_not_called()


def test_preview_view_resize_with_scene_rect_schedules_timer(app):
    """resizeEvent schedules a single-shot timer when original_scene_rect is set."""
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
    """refit_view does not raise when original_scene_rect is None."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    view.refit_view()  # should not raise


def test_refit_view_empty_rect_is_noop(app):
    """refit_view does not raise when original_scene_rect is an empty QRectF."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF()  # empty
    view.refit_view()  # should not raise


def test_refit_view_valid_rect(app):
    """refit_view does not raise for a valid non-empty QRectF."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = QRectF(0, 0, 100, 100)
    view.refit_view()  # should not raise


# ---------------------------------------------------------------------------
# TemplatePreviewView — showEvent
# ---------------------------------------------------------------------------


def test_show_event_no_rect_no_timer(app):
    """showEvent does not schedule a timer when original_scene_rect is None."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.original_scene_rect = None
    ev = QShowEvent()
    with patch("moleditpy.ui.template_preview_view.QTimer.singleShot") as mock_timer:
        view.showEvent(ev)
    mock_timer.assert_not_called()


def test_show_event_with_rect_schedules_timer(app):
    """showEvent schedules a single-shot timer when original_scene_rect is set."""
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
    """redraw_with_current_size returns without raising when template_data is None."""
    scene = QGraphicsScene()
    view = TemplatePreviewView(scene)
    view.template_data = None
    view.parent_dialog = None
    view.redraw_with_current_size()  # should not raise


def test_redraw_with_data_calls_draw_template_preview(app):
    """redraw_with_current_size calls parent_dialog.draw_template_preview with the stored data."""
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
    """load_template_file returns the parsed dict for a valid JSON template file."""
    dlg, _ = _make_dlg(app, tmp_path)
    data = {"name": "benzene", "atoms": [], "bonds": []}
    fp = tmp_path / "benzene.pmetmplt"
    fp.write_text(json.dumps(data), encoding="utf-8")
    result = dlg.load_template_file(str(fp))
    assert result == data


def test_load_template_file_missing(app, tmp_path):
    """load_template_file returns None for a nonexistent file path."""
    dlg, _ = _make_dlg(app, tmp_path)
    result = dlg.load_template_file(str(tmp_path / "nonexistent.pmetmplt"))
    assert result is None


def test_load_template_file_invalid_json(app, tmp_path):
    """load_template_file returns None when the file contains invalid JSON."""
    dlg, _ = _make_dlg(app, tmp_path)
    fp = tmp_path / "bad.pmetmplt"
    fp.write_text("NOT JSON", encoding="utf-8")
    result = dlg.load_template_file(str(fp))
    assert result is None


# ---------------------------------------------------------------------------
# UserTemplateDialog — save_template_file
# ---------------------------------------------------------------------------


def test_save_template_file_success(app, tmp_path):
    """save_template_file writes JSON to disk and returns True on success."""
    dlg, _ = _make_dlg(app, tmp_path)
    data = {"name": "ethane"}
    fp = tmp_path / "ethane.pmetmplt"
    result = dlg.save_template_file(str(fp), data)
    assert result is True
    assert json.loads(fp.read_text()) == data


def test_save_template_file_oserror(app, tmp_path):
    """save_template_file returns False when the file cannot be written."""
    dlg, _ = _make_dlg(app, tmp_path)
    result = dlg.save_template_file("/nonexistent_dir/bad.pmetmplt", {})
    assert result is False


# ---------------------------------------------------------------------------
# UserTemplateDialog — draw_template_preview
# ---------------------------------------------------------------------------


def _scene():
    return QGraphicsScene()


def test_draw_template_preview_no_atoms_adds_placeholder(app, tmp_path):
    """draw_template_preview adds a placeholder text item when the template has no atoms."""
    dlg, _ = _make_dlg(app, tmp_path)
    sc = _scene()
    dlg.draw_template_preview(sc, {"atoms": [], "bonds": []})
    assert sc.items()  # placeholder text item


def test_draw_template_preview_single_bond(app, tmp_path):
    """draw_template_preview renders a single bond without raising."""
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
    """draw_template_preview renders a double bond without raising."""
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
    """draw_template_preview renders a triple bond without raising."""
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
    """draw_template_preview renders non-carbon atom labels without raising."""
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
    """convert_structure_to_template returns a dict with the correct name, atoms, bonds, and format."""
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
    """cleanup_template_mode clears selected_template and disables the delete button."""
    dlg, mw = _make_dlg(app, tmp_path)
    dlg.selected_template = {"name": "foo"}
    dlg.cleanup_template_mode()
    assert dlg.selected_template is None
    assert dlg.delete_button.isEnabled() is False


def test_cleanup_template_mode_calls_set_mode(app, tmp_path):
    """cleanup_template_mode calls ui_manager.set_mode_and_update_toolbar with 'atom_C'."""
    dlg, mw = _make_dlg(app, tmp_path)
    mw.ui_manager = MagicMock()
    dlg.cleanup_template_mode()
    mw.ui_manager.set_mode_and_update_toolbar.assert_called_once_with("atom_C")


def test_cleanup_template_mode_no_ui_manager_no_crash(app, tmp_path):
    """cleanup_template_mode does not raise when ui_manager is absent."""
    dlg, mw = _make_dlg(app, tmp_path)
    del mw.ui_manager
    dlg.cleanup_template_mode()  # should not raise


def test_cleanup_template_mode_scene_reset(app, tmp_path):
    """cleanup_template_mode resets scene.mode and scene.user_template_data."""
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
    """select_template stores the template data and enables the delete button."""
    dlg, mw = _make_dlg(app, tmp_path)
    data = {"name": "benzene", "atoms": [], "bonds": []}
    widget = MagicMock()
    with patch.object(dlg, "_activate_template_mode"):
        dlg.select_template(data, widget)
    assert dlg.selected_template == data
    assert dlg.delete_button.isEnabled() is True


def test_select_template_activates_template_mode(app, tmp_path):
    """select_template calls _activate_template_mode with the chosen template data."""
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
    """use_template calls _activate_template_mode and sets selected_template."""
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
    """fit_preview_view_safely calls fitInView for a non-empty rect."""
    dlg, _ = _make_dlg(app, tmp_path)
    view = MagicMock()
    rect = QRectF(0, 0, 100, 100)
    dlg.fit_preview_view_safely(view, rect)
    view.fitInView.assert_called_once()


def test_fit_preview_view_safely_empty_rect(app, tmp_path):
    """fit_preview_view_safely skips fitInView for an empty QRectF."""
    dlg, _ = _make_dlg(app, tmp_path)
    view = MagicMock()
    dlg.fit_preview_view_safely(view, QRectF())  # empty → no-op
    view.fitInView.assert_not_called()


def test_fit_preview_view_safely_none_view(app, tmp_path):
    """fit_preview_view_safely does not raise when view is None."""
    dlg, _ = _make_dlg(app, tmp_path)
    # None view should not raise (guard: `if view and ...`)
    dlg.fit_preview_view_safely(None, QRectF(0, 0, 10, 10))


# ---------------------------------------------------------------------------
# UserTemplateDialog — refit_all_previews
# ---------------------------------------------------------------------------


def test_refit_all_previews_no_items_no_crash(app, tmp_path):
    """refit_all_previews does not raise when the template grid is empty."""
    dlg, _ = _make_dlg(app, tmp_path)
    dlg.refit_all_previews()  # empty grid — should not raise


def test_refit_all_previews_calls_redraw(app, tmp_path):
    """refit_all_previews calls redraw_with_current_size on each TemplatePreviewView child."""
    dlg, _ = _make_dlg(app, tmp_path)
    preview_view = MagicMock(spec=TemplatePreviewView)
    preview_view.redraw_with_current_size = MagicMock()

    child_widget = MagicMock()
    child_widget.findChildren.return_value = [preview_view]

    layout_item = MagicMock()
    layout_item.widget.return_value = child_widget

    with (
        patch.object(dlg.template_layout, "count", return_value=1),
        patch.object(dlg.template_layout, "itemAt", return_value=layout_item),
    ):
        dlg.refit_all_previews()

    preview_view.redraw_with_current_size.assert_called_once()


# ---------------------------------------------------------------------------
# UserTemplateDialog — delete_selected_template (no selection)
# ---------------------------------------------------------------------------


def test_delete_selected_no_selection_returns_early(app, tmp_path):
    """delete_selected_template returns early without prompting when no template is selected."""
    dlg, _ = _make_dlg(app, tmp_path)
    dlg.selected_template = None
    with patch("moleditpy.ui.user_template_dialog.QMessageBox.question") as mock_q:
        dlg.delete_selected_template()
    mock_q.assert_not_called()
