"""The template ghost preview has to actually appear in the real 2D editor."""

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


def test_ring_template_preview_appears(window):
    """A built-in ring template previews as a ghost ring, Kekulé bonds included."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("template_benzene")

    scene.update_template_preview(QPointF(320, 260))

    preview = scene.template_preview
    assert preview.isVisible()
    assert [b.order for b in preview.ghost_bonds] == [2, 1, 2, 1, 2, 1]
    assert all(b.is_in_ring for b in preview.ghost_bonds)
