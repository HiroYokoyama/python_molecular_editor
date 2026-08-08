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


def test_hover_and_snap_use_same_default_radius(scene_with_atom):
    """At default settings both radii must equal bond_snapping_distance_2d."""
    scene, atom, _, _ = scene_with_atom

    # Hover: shape() radius (in scene units; at scale=1 ≡ screen px)
    hover_radius = atom.shape().boundingRect().width() / 2

    # Snap: find_atom_near at the boundary distance
    assert (
        scene.find_atom_near(QPointF(DEFAULT_SNAP_PX - 0.1, 0), tol=DEFAULT_SNAP_PX)
        is atom
    ), "atom should be found inside snapping distance"
    assert (
        scene.find_atom_near(QPointF(DEFAULT_SNAP_PX + 0.1, 0), tol=DEFAULT_SNAP_PX)
        is None
    ), "atom should NOT be found outside snapping distance"

    assert hover_radius == pytest.approx(DEFAULT_SNAP_PX, abs=0.5), (
        f"hover radius ({hover_radius}) must match snap distance ({DEFAULT_SNAP_PX})"
    )


def test_hover_and_snap_follow_custom_setting(scene_with_atom):
    """When the setting is changed, both hover and snap must change together."""
    scene, atom, _, settings = scene_with_atom
    custom_dist = 25.0
    settings["bond_snapping_distance_2d"] = custom_dist

    hover_radius = atom.shape().boundingRect().width() / 2

    assert scene.find_atom_near(QPointF(custom_dist - 0.1, 0), tol=custom_dist) is atom
    assert scene.find_atom_near(QPointF(custom_dist + 0.1, 0), tol=custom_dist) is None
    assert hover_radius == pytest.approx(custom_dist, abs=0.5)


@pytest.mark.parametrize("scale", [0.5, 1.0, 2.0, 4.0])
def test_hover_and_snap_stay_consistent_at_different_zoom_levels(
    scene_with_atom, scale
):
    """Both radii must stay identical in screen-pixel terms at any zoom."""
    scene, atom, mock_view, _ = scene_with_atom
    mock_view.transform.return_value = QTransform.fromScale(scale, scale)

    # Hover: shape() returns scene-unit radius; multiply by scale → screen px
    hover_scene_radius = atom.shape().boundingRect().width() / 2
    hover_screen_px = hover_scene_radius * scale

    # Snap: find_atom_near converts tol (screen px) to scene units internally.
    # The boundary in scene units is DEFAULT_SNAP_PX / scale.
    boundary_scene = DEFAULT_SNAP_PX / scale
    assert (
        scene.find_atom_near(QPointF(boundary_scene - 0.01, 0), tol=DEFAULT_SNAP_PX)
        is atom
    ), f"scale={scale}: atom should be found inside snap zone"
    assert (
        scene.find_atom_near(QPointF(boundary_scene + 0.5, 0), tol=DEFAULT_SNAP_PX)
        is None
    ), f"scale={scale}: atom should NOT be found outside snap zone"

    assert hover_screen_px == pytest.approx(DEFAULT_SNAP_PX, abs=0.5), (
        f"scale={scale}: hover screen radius ({hover_screen_px}) != "
        f"snap distance ({DEFAULT_SNAP_PX})"
    )
