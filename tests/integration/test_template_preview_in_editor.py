"""The template ghost preview has to actually appear in the real 2D editor."""

import math
from unittest.mock import patch

from PyQt6.QtCore import QPointF, QRectF
from PyQt6.QtGui import QColor, QImage, QPainter

PYRIDINE = {
    "name": "pyridine",
    "atoms": [
        {"id": 0, "symbol": "N", "x": 0, "y": 0},
        {"id": 1, "symbol": "C", "x": 43, "y": 25},
        {"id": 2, "symbol": "C", "x": 43, "y": 75},
        {"id": 3, "symbol": "C", "x": 0, "y": 100},
        {"id": 4, "symbol": "C", "x": -43, "y": 75},
        {"id": 5, "symbol": "C", "x": -43, "y": 25},
    ],
    "bonds": [
        {"atom1": 0, "atom2": 1, "order": 2},
        {"atom1": 1, "atom2": 2, "order": 1},
        {"atom1": 2, "atom2": 3, "order": 2},
        {"atom1": 3, "atom2": 4, "order": 1},
        {"atom1": 4, "atom2": 5, "order": 2},
        {"atom1": 5, "atom2": 0, "order": 1},
    ],
}


def _painted_pixels(scene, region):
    """Render a scene region and count the pixels that are not background."""
    image = QImage(240, 240, QImage.Format.Format_ARGB32)
    image.fill(QColor("white"))
    painter = QPainter(image)
    scene.render(painter, target=QRectF(image.rect()), source=region)
    painter.end()
    return sum(
        1
        for y in range(image.height())
        for x in range(image.width())
        if image.pixel(x, y) != 0xFFFFFFFF
    )


def test_user_template_preview_appears(window):
    """Moving the cursor in user-template mode shows the ghost molecule."""
    scene = window.init_manager.scene
    scene.user_template_data = PYRIDINE
    window.ui_manager.set_mode("template_user_pyridine")

    scene.update_template_preview(QPointF(300, 250))

    preview = scene.template_preview
    assert preview.isVisible()
    assert len(preview.ghost_atoms) == 6
    assert len(preview.ghost_bonds) == 6
    assert not preview.boundingRect().isEmpty()


def test_preview_is_indexed_where_it_is_drawn(window):
    """Regression: the scene index kept an empty rect from before the ghost was
    built, so items(exposed) never returned the preview and Qt culled it — the
    ghost was invisible in the editor even though the item was visible."""
    scene = window.init_manager.scene
    first = scene.create_atom("C", QPointF(0, 0))
    second = scene.create_atom("O", QPointF(50, 0))
    scene.create_bond(scene.atom_items[first], scene.atom_items[second], 1, 0)
    scene.update_all_items()

    scene.user_template_data = PYRIDINE
    window.ui_manager.set_mode("template_user_pyridine")
    # Paint once first: that is what puts the preview into the scene index
    _painted_pixels(scene, QRectF(-100, -100, 600, 600))

    where = QPointF(300, 250)
    scene.update_template_preview(where)

    preview = scene.template_preview
    assert preview.sceneBoundingRect().contains(where)
    assert preview in scene.items(preview.sceneBoundingRect())
    # Rendering the scene is what the editor does; the ghost must show up in it
    assert _painted_pixels(scene, preview.sceneBoundingRect()) > 0


def test_fused_atom_keeps_its_own_element(window):
    """A ring vertex landing on an existing atom previews that atom, not a carbon,
    and the editor keeps drawing it."""
    scene = window.init_manager.scene
    nitrogen = scene.create_atom("N", QPointF(0, 0))
    carbon = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[nitrogen], scene.atom_items[carbon], 1, 0)
    scene.update_all_items()

    window.ui_manager.set_mode("template_benzene")
    scene.update_template_preview(QPointF(25, 2))

    preview = scene.template_preview
    fused = [i for i, a in enumerate(preview.ghost_atoms) if a.symbol == "N"]
    assert fused, "the ring should have picked up the existing nitrogen"
    assert fused[0] in preview.existing_indices
    assert preview.ghost_atoms[fused[0]].pos() == scene.atom_items[nitrogen].pos()


def test_fusing_onto_another_atom_relabels_the_ghost(window):
    """The fused atom's hydrogens come from the editor, so a ghost reused across
    cursor moves must not keep the label of the atom it was fused to before."""
    scene = window.init_manager.scene
    lone = scene.create_atom("O", QPointF(0, 0))
    bonded = scene.create_atom("O", QPointF(300, 0))
    partner = scene.create_atom("C", QPointF(350, 0))
    scene.create_bond(scene.atom_items[bonded], scene.atom_items[partner], 1, 0)
    scene.atom_items[lone].implicit_h_count = 2
    scene.atom_items[bonded].implicit_h_count = 1
    scene.update_all_items()

    window.ui_manager.set_mode("template_benzene")
    preview = scene.template_preview

    scene.update_template_preview(scene.atom_items[lone].pos())
    fused = [i for i in preview.existing_indices]
    assert fused and preview.ghost_atoms[fused[0]].implicit_h_count == 2

    scene.update_template_preview(scene.atom_items[bonded].pos())
    fused = [i for i in preview.existing_indices]
    assert fused and preview.ghost_atoms[fused[0]].implicit_h_count == 1


def test_user_template_previews_its_own_atoms_over_existing_ones(window):
    """Only the first atom of a user template reuses an existing one; a vertex that
    merely lands on another becomes a new atom, so the ghost must not fuse it."""
    scene = window.init_manager.scene
    anchor = scene.create_atom("C", QPointF(0, 0))
    scene.create_atom("N", QPointF(50, 0))
    scene.update_all_items()

    scene.user_template_data = {
        "name": "ethanol_stub",
        "atoms": [
            {"id": 0, "symbol": "C", "x": 0, "y": 0},
            {"id": 1, "symbol": "C", "x": 50, "y": 0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 1}],
    }
    window.ui_manager.set_mode("template_user_ethanol_stub")
    scene.update_template_preview(scene.atom_items[anchor].pos())

    preview = scene.template_preview
    assert preview.ghost_atoms[1].symbol == "C"
    assert not preview.existing_indices

    scene.add_user_template_fragment(scene.template_context)
    assert len(scene.data.atoms) == 3


def test_ghost_labels_use_the_editors_font_settings(window):
    """update_style() reads the font through scene(), so the ghost items have to be
    styled after they join a scene or they keep the built-in default size."""
    scene = window.init_manager.scene
    window.init_manager.settings["atom_font_size_2d"] = 13
    scene.user_template_data = PYRIDINE
    window.ui_manager.set_mode("template_user_pyridine")

    scene.update_template_preview(QPointF(300, 250))

    assert scene.template_preview.ghost_atoms[0].font.pointSize() == 13


def test_fused_ring_preview_uses_the_placement_rotation(window):
    """A fused aromatic ring is rotated to fit the bonds already there, so the
    preview has to show the double bonds where the click will put them — while
    the placement context keeps the unrotated orders it rotates itself."""
    scene = window.init_manager.scene
    left = scene.create_atom("C", QPointF(0, 0))
    right = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[left], scene.atom_items[right], 1, 0)
    scene.update_all_items()

    window.ui_manager.set_mode("template_benzene")
    with patch.object(
        type(scene), "_calculate_6ring_rotation", return_value=1, create=False
    ):
        scene.update_template_preview(QPointF(25, 2))

    preview = scene.template_preview
    assert [bond.order for bond in preview.ghost_bonds] == [1, 2, 1, 2, 1, 2]
    assert [order for (_, _, order) in scene.template_context["bonds_info"]] == [
        2,
        1,
        2,
        1,
        2,
        1,
    ]


def test_leaving_template_mode_clears_the_dialog_selection(window):
    """Picking another mode has to unhighlight the template dialog's selection."""
    from moleditpy.ui.user_template_dialog import UserTemplateDialog

    with patch.object(UserTemplateDialog, "load_user_templates"):
        dialog = UserTemplateDialog(window, window)
    window.template_dialog = dialog
    dialog.selected_template = PYRIDINE
    dialog.delete_button.setEnabled(True)

    window.ui_manager.set_mode("select")

    assert dialog.selected_template is None
    assert not dialog.delete_button.isEnabled()


def test_using_a_user_template_checks_the_user_toolbar_button(window):
    """The toolbar has to say which mode is live: picking a user template checks
    USER and drops the highlight from whatever was active before."""
    from moleditpy.ui.user_template_dialog import UserTemplateDialog

    with patch.object(UserTemplateDialog, "load_user_templates"):
        dialog = UserTemplateDialog(window, window)
    window.template_dialog = dialog
    actions = window.init_manager.mode_actions
    window.ui_manager.set_mode_and_update_toolbar("atom_C")

    dialog.use_template(PYRIDINE)

    assert window.init_manager.scene.mode == "template_user_pyridine"
    assert actions["template_user"].isChecked()
    assert not actions["atom_C"].isChecked()


def test_escape_leaves_template_mode_like_the_close_button(window):
    """Esc calls reject() and never reaches closeEvent, so it used to leave the
    editor armed with a template the dialog had already forgotten."""
    from moleditpy.ui.user_template_dialog import UserTemplateDialog

    with patch.object(UserTemplateDialog, "load_user_templates"):
        dialog = UserTemplateDialog(window, window)
    window.template_dialog = dialog
    dialog.use_template(PYRIDINE)
    assert window.init_manager.scene.mode == "template_user_pyridine"

    dialog.reject()

    assert window.init_manager.scene.mode == "atom_C"
    assert dialog.selected_template is None
    assert not window.init_manager.mode_actions["template_user"].isChecked()


def test_ring_template_preview_appears(window):
    """A built-in ring template previews as a ghost ring, Kekulé bonds included."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("template_benzene")

    scene.update_template_preview(QPointF(320, 260))

    preview = scene.template_preview
    assert preview.isVisible()
    assert [b.order for b in preview.ghost_bonds] == [2, 1, 2, 1, 2, 1]
    assert all(b.is_in_ring for b in preview.ghost_bonds)


# ---------------------------------------------------------------------------
# Fusing: the ghost must show the order placement will leave on a shared bond
# ---------------------------------------------------------------------------

_HEX_SIDE = 75.0  # the app's default 2D bond length
_HEX_START = math.radians(-120.0)  # matches the ring template's start angle
_HEX_STEP = math.radians(60.0)
_KEKULE = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]


def _hex_points(centre=QPointF(0.0, 0.0)):
    """The six vertices a benzene template lands on when centred at *centre*."""
    return [
        centre
        + QPointF(
            _HEX_SIDE * math.cos(_HEX_START + i * _HEX_STEP),
            _HEX_SIDE * math.sin(_HEX_START + i * _HEX_STEP),
        )
        for i in range(6)
    ]


def _shared_ghost_orders(preview):
    """Orders of ghost bonds whose two ends both sit on existing atoms."""
    return [
        bond.order
        for bond in preview.ghost_bonds
        if bond.atom1.atom_id in preview.existing_indices
        and bond.atom2.atom_id in preview.existing_indices
    ]


def test_fusing_ghost_keeps_an_untouched_bond_order(window):
    """Fusing onto a ring must not draw a double over a bond that stays single.

    add_molecule_fragment refuses to re-order an existing bond whose atoms
    already carry a double, so the alternation the rotation picks never reaches
    the shared bond. The ghost used to draw it anyway, leaving a second line
    along every single bond of the ring being fused to.
    """
    scene = window.init_manager.scene
    points = _hex_points()
    atoms = [scene.atom_items[scene.create_atom("C", p)] for p in points]
    for i, j, order in _KEKULE:
        scene.create_bond(atoms[i], atoms[j], bond_order=order)
    window.ui_manager.set_mode("template_benzene")

    apothem = _HEX_SIDE * math.sqrt(3) / 2.0
    for i, j, order in _KEKULE:
        mid = (points[i] + points[j]) / 2.0
        length = math.hypot(mid.x(), mid.y()) or 1.0
        scene.update_template_preview(
            mid + QPointF(mid.x() / length * apothem, mid.y() / length * apothem)
        )
        preview = scene.template_preview

        shared = _shared_ghost_orders(preview)
        assert shared, f"nothing fused when hovering over bond ({i}, {j})"
        assert shared == [order], (
            f"ghost promises {shared} on bond ({i}, {j}), which stays {order}"
        )


def test_fusing_ghost_still_previews_an_allowed_overwrite(window):
    """A lone single bond does get upgraded, and the ghost must say so."""
    scene = window.init_manager.scene
    points = _hex_points()
    first = scene.atom_items[scene.create_atom("C", points[0])]
    second = scene.atom_items[scene.create_atom("C", points[1])]
    scene.create_bond(first, second, bond_order=1)
    window.ui_manager.set_mode("template_benzene")

    scene.update_template_preview(QPointF(0.0, 0.0))
    preview = scene.template_preview
    assert _shared_ghost_orders(preview) == [2]

    context = scene.template_context
    scene.add_molecule_fragment(
        list(context["points"]),
        list(context["bonds_info"]),
        context.get("items") or [],
    )

    assert scene.find_bond_between(first, second).order == 2


def test_fusing_ghost_does_not_redraw_an_unchanged_bond(window):
    """The ghost leaves a bond it will not change to the editor.

    Painting a ghost copy on top of the real bond doubles every line: over an
    existing double bond the overlay reads as a triple.
    """
    scene = window.init_manager.scene
    points = _hex_points()
    atoms = [scene.atom_items[scene.create_atom("C", p)] for p in points]
    for i, j, order in _KEKULE:
        scene.create_bond(atoms[i], atoms[j], bond_order=order)
    window.ui_manager.set_mode("template_benzene")

    # Fuse onto bond (0, 1), which is a double and stays one
    apothem = _HEX_SIDE * math.sqrt(3) / 2.0
    mid = (points[0] + points[1]) / 2.0
    length = math.hypot(mid.x(), mid.y()) or 1.0
    scene.update_template_preview(
        mid + QPointF(mid.x() / length * apothem, mid.y() / length * apothem)
    )
    preview = scene.template_preview
    assert preview.isVisible()

    # Ghost items are indexed by template vertex, not by editor atom id
    shared_pairs = {
        frozenset((b.atom1.atom_id, b.atom2.atom_id))
        for b in preview.ghost_bonds
        if b.atom1.atom_id in preview.existing_indices
        and b.atom2.atom_id in preview.existing_indices
    }
    assert len(shared_pairs) == 1, "expected exactly one fused bond"
    assert shared_pairs <= preview.editor_drawn_bonds

    # A tight box on the shared bond must look the same with and without the ghost
    region = QRectF(mid.x() - 12, mid.y() - 12, 24, 24)
    with_ghost = _painted_pixels(scene, region)
    preview.hide()
    without_ghost = _painted_pixels(scene, region)

    assert with_ghost == without_ghost, (
        f"ghost added {with_ghost - without_ghost} px over an unchanged bond"
    )


def test_user_template_ghost_ignores_a_previous_fusion(window):
    """Switching template kinds must not carry the fused pairs across.

    The vertex pairs a fused ring records are indices, so they collide with a
    user template's own vertices and would silently hide one of its bonds.
    """
    scene = window.init_manager.scene
    points = _hex_points()
    first = scene.atom_items[scene.create_atom("C", points[0])]
    second = scene.atom_items[scene.create_atom("C", points[1])]
    scene.create_bond(first, second, bond_order=2)

    window.ui_manager.set_mode("template_benzene")
    scene.update_template_preview(QPointF(0.0, 0.0))
    assert scene.template_preview.editor_drawn_bonds, "expected a fused pair"

    scene.user_template_data = PYRIDINE
    window.ui_manager.set_mode("template_user_pyridine")
    scene.update_template_preview(QPointF(600.0, 600.0))

    assert scene.template_preview.editor_drawn_bonds == set()
    hidden = [
        (b.atom1.atom_id, b.atom2.atom_id)
        for b in scene.template_preview.ghost_bonds
        if frozenset((b.atom1.atom_id, b.atom2.atom_id))
        in scene.template_preview.editor_drawn_bonds
    ]
    assert not hidden, f"user-template ghost bonds suppressed: {hidden}"
