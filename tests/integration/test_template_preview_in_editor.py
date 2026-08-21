"""The template ghost preview has to actually appear in the real 2D editor."""

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


def test_ring_template_preview_appears(window):
    """A built-in ring template previews as a ghost ring, Kekulé bonds included."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("template_benzene")

    scene.update_template_preview(QPointF(320, 260))

    preview = scene.template_preview
    assert preview.isVisible()
    assert [b.order for b in preview.ghost_bonds] == [2, 1, 2, 1, 2, 1]
    assert all(b.is_in_ring for b in preview.ghost_bonds)
