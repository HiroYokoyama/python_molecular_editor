"""Tests verifying that hover highlight and snapping use the same distance.

The ``AtomItem.shape()`` hit radius and ``find_atom_near()`` snap radius
must agree so that an atom is keyboard-editable exactly when the user sees
the hover highlight.  Both radii are expressed in screen pixels and derived
from the ``bond_snapping_distance_2d`` setting.
"""

import pytest
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QTransform
from moleditpy.ui.molecule_scene import MoleculeScene
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

DEFAULT_SNAP_PX = 14.0


@pytest.fixture
def scene_with_atom(app):
    """A MoleculeScene containing a single carbon atom at the origin.

    ``init_manager.settings`` is a real dict so ``get_setting()`` works
    exactly as it does in production code.
    """
    settings = {"bond_snapping_distance_2d": DEFAULT_SNAP_PX}

    mock_init_manager = MagicMock()
    mock_init_manager.settings = settings

    mock_window = MagicMock()
    mock_window.is_2d_editable = True
    mock_window.init_manager = mock_init_manager

    data = MagicMock()
    data.atoms = {}
    data.bonds = {}

    scene = MoleculeScene(data, mock_window)

    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()  # identity (m11 = 1.0)
    scene.views = MagicMock(return_value=[mock_view])

    atom = AtomItem(0, "C", QPointF(0, 0))
    atom.atom_id = 0
    scene.addItem(atom)

    yield scene, atom, mock_view, settings


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _glyph_extent(atom):
    """Half-width / half-height of the atom's drawn label, in scene units."""
    rect = atom.get_bg_ellipse_path().boundingRect()
    return rect.width() / 2, rect.height() / 2


def test_hover_and_snap_use_same_default_radius(scene_with_atom):
    """At default settings the hit radius must be the setting (glyph aside)."""
    scene, atom, _, _ = scene_with_atom

    assert atom.hit_radius() == pytest.approx(DEFAULT_SNAP_PX)

    # Snap: find_atom_near at the boundary distance, probed along the axis
    # where the glyph does not extend past the snap circle.
    _, glyph_h = _glyph_extent(atom)
    assert glyph_h > DEFAULT_SNAP_PX  # the label is taller than the snap circle
    assert (
        scene.find_atom_near(QPointF(DEFAULT_SNAP_PX - 0.1, 0), tol=DEFAULT_SNAP_PX)
        is atom
    ), "atom should be found inside snapping distance"
    far = QPointF(0, glyph_h + DEFAULT_SNAP_PX + 10.0)
    assert scene.find_atom_near(far, tol=DEFAULT_SNAP_PX) is None, (
        "atom should NOT be found far outside snapping distance and glyph"
    )


def test_hover_and_snap_follow_custom_setting(scene_with_atom):
    """When the setting is changed, both hover and snap must change together."""
    scene, atom, _, settings = scene_with_atom
    custom_dist = 25.0
    settings["bond_snapping_distance_2d"] = custom_dist

    assert atom.hit_radius() == pytest.approx(custom_dist)

    _, glyph_h = _glyph_extent(atom)
    assert scene.find_atom_near(QPointF(custom_dist - 0.1, 0), tol=custom_dist) is atom
    far = QPointF(0, max(custom_dist, glyph_h) + 10.0)
    assert scene.find_atom_near(far, tol=custom_dist) is None


@pytest.mark.parametrize("scale", [0.5, 1.0, 2.0, 4.0])
def test_hover_and_snap_stay_consistent_at_different_zoom_levels(
    scene_with_atom, scale
):
    """Hover shape and snapping must agree at every point, at any zoom.

    This is the invariant that matters: an atom is keyboard-editable exactly
    when the user sees its hover highlight.
    """
    scene, atom, mock_view, _ = scene_with_atom
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    # The pixel-sized radius stays constant on screen.
    assert atom.hit_radius() * scale == pytest.approx(DEFAULT_SNAP_PX)

    glyph_w, glyph_h = _glyph_extent(atom)
    reach = max(glyph_w, glyph_h, DEFAULT_SNAP_PX / scale) + 12.0
    probes = []
    for i in range(1, 13):
        d = reach * i / 12.0
        probes += [QPointF(d, 0), QPointF(0, d), QPointF(d * 0.7, d * 0.7)]

    for p in probes:
        hovered = atom.contains(atom.mapFromScene(p))
        snapped = scene.find_atom_near(p, tol=DEFAULT_SNAP_PX) is atom
        assert hovered == snapped, (
            f"scale={scale}: hover ({hovered}) and snap ({snapped}) disagree at {p}"
        )


@pytest.mark.parametrize("scale", [0.5, 1.0, 2.0, 4.0, 8.0])
def test_hit_shape_stays_inside_bounding_rect(scene_with_atom, scale):
    """Qt prunes hits by boundingRect, so shape() must never escape it."""
    scene, atom, mock_view, _ = scene_with_atom
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    bounding = atom.boundingRect()
    shape_rect = atom.shape().boundingRect()
    assert bounding.contains(shape_rect), (
        f"scale={scale}: shape {shape_rect} escapes boundingRect {bounding}"
    )


@pytest.mark.parametrize("scale", [4.0, 8.0, 16.0])
def test_drawn_label_is_always_clickable_when_zoomed_in(scene_with_atom, scale):
    """Zoomed in, the pixel radius shrinks below the glyph — it must still hit."""
    scene, atom, mock_view, _ = scene_with_atom
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    glyph_w, _ = _glyph_extent(atom)
    assert glyph_w > atom.hit_radius()  # the regression precondition

    on_glyph = QPointF(glyph_w * 0.7, 0)
    assert atom.contains(atom.mapFromScene(on_glyph)), "glyph must be hoverable"
    assert scene.find_atom_near(on_glyph, tol=DEFAULT_SNAP_PX) is atom, (
        "glyph must be snappable"
    )


@pytest.fixture
def scene_with_bond(app):
    """A MoleculeScene containing a single bond between two atoms."""
    settings = {"bond_snapping_distance_2d": DEFAULT_SNAP_PX}

    mock_init_manager = MagicMock()
    mock_init_manager.settings = settings

    mock_window = MagicMock()
    mock_window.is_2d_editable = True
    mock_window.init_manager = mock_init_manager

    data = MagicMock()
    data.atoms = {}
    data.bonds = {}

    scene = MoleculeScene(data, mock_window)

    mock_view = MagicMock()
    mock_view.transform.return_value = QTransform()  # identity (m11 = 1.0)
    scene.views = MagicMock(return_value=[mock_view])

    # Bond between (-50, 0) and (50, 0) (Horizontal)
    atom1 = AtomItem(0, "C", QPointF(-50, 0))
    atom2 = AtomItem(1, "C", QPointF(50, 0))
    scene.addItem(atom1)
    scene.addItem(atom2)

    bond = BondItem(atom1, atom2, order=1, stereo=0)
    scene.addItem(bond)

    yield scene, bond, mock_view, settings


def test_bond_hover_uses_default_radius(scene_with_bond):
    """At default settings, bond shape must use bond_snapping_distance_2d."""
    scene, bond, _, _ = scene_with_bond

    # For a perfectly horizontal line at y=0, the height of the shape bounding rect
    # is exactly equal to the stroker width. Since width = 2 * radius,
    # radius is height / 2.
    hover_radius = bond.shape().boundingRect().height() / 2

    assert hover_radius == pytest.approx(DEFAULT_SNAP_PX, abs=0.5)


def test_bond_hover_follows_custom_setting(scene_with_bond):
    """When the setting is changed, bond hover must change."""
    scene, bond, _, settings = scene_with_bond
    custom_dist = 25.0
    settings["bond_snapping_distance_2d"] = custom_dist

    hover_radius = bond.shape().boundingRect().height() / 2
    assert hover_radius == pytest.approx(custom_dist, abs=0.5)


@pytest.mark.parametrize("scale", [0.05, 0.25, 1.0, 4.0])
def test_bond_hit_shape_stays_inside_bounding_rect(scene_with_bond, scale):
    """Zoomed far out the hit band grows in scene units; boundingRect must follow."""
    _scene, bond, mock_view, settings = scene_with_bond
    settings["bond_snapping_distance_2d"] = 50.0  # slider maximum
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    bounding = bond.boundingRect()
    shape_rect = bond.shape().boundingRect()
    assert bounding.contains(shape_rect), (
        f"scale={scale}: shape {shape_rect} escapes boundingRect {bounding}"
    )


@pytest.mark.parametrize("scale", [0.5, 1.0, 2.0, 4.0])
def test_bond_hover_stays_consistent_at_different_zoom_levels(scene_with_bond, scale):
    """Bond radius must stay identical in screen-pixel terms at any zoom."""
    scene, bond, mock_view, _ = scene_with_bond
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    # Hover: shape() returns scene-unit dimensions; multiply by scale → screen px
    hover_scene_radius = bond.shape().boundingRect().height() / 2
    hover_screen_px = hover_scene_radius * scale

    assert hover_screen_px == pytest.approx(DEFAULT_SNAP_PX, abs=0.5)
