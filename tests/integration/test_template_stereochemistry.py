"""R/S annotations survive the complete user-template workflow."""

import pytest
from PySide6.QtCore import QPointF

from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.preview_molecule import PreviewScene
from moleditpy.ui.template_preview_item import TemplatePreviewItem
from moleditpy.ui.user_template_dialog import UserTemplateDialog


@pytest.mark.parametrize("label", ["R", "S", None])
@pytest.mark.parametrize("attach", [False, True])
def test_template_preserves_rs_through_save_preview_and_placement(
    window, qtbot, tmp_path, label, attach
):
    scene = window.init_manager.scene
    center = scene.create_atom("C", QPointF(0, 0))
    for symbol, point in zip(
        ("F", "Cl", "Br", "H"),
        (QPointF(50, 0), QPointF(0, 50), QPointF(-50, 0), QPointF(0, -50)),
    ):
        atom = scene.create_atom(symbol, point)
        scene.create_bond(scene.atom_items[center], scene.atom_items[atom])
    scene.set_atom_rs_label(scene.atom_items[center], label)
    core_template = window.state_manager.data.to_template_dict("Stereo")
    assert core_template["atoms"][0].get("rs") == label

    dialog = UserTemplateDialog(window, parent=window)
    qtbot.addWidget(dialog)
    template = dialog.convert_structure_to_template("Stereo")
    assert template["atoms"][0].get("rs") == label
    path = tmp_path / "stereo.pmetmplt"
    assert dialog.save_template_file(str(path), template)
    template = dialog.load_template_file(str(path))
    preview_scene = PreviewScene(scene)
    dialog.draw_template_preview(preview_scene, template)
    preview_atoms = [i for i in preview_scene.items() if isinstance(i, AtomItem)]
    assert next(i for i in preview_atoms if i.symbol == "C").rs_label == label

    ghost = TemplatePreviewItem()
    ghost.set_user_template_geometry([QPointF(0, 0)], [], [{"symbol": "C", "rs": "S"}])
    ghost.set_user_template_geometry(
        [QPointF(10, 0)], [], [{"symbol": "C", "rs": label}]
    )
    assert ghost.ghost_atoms[0].rs_label == label

    window.edit_actions_manager.clear_all(skip_check=True)
    if attach:
        existing_id = scene.create_atom("C", QPointF(200, 200))
        scene.set_atom_rs_label(scene.atom_items[existing_id], "R")
    scene.user_template_data = template
    scene.update_user_template_preview(QPointF(200, 200))
    assert (scene.template_context["attachment_atom"] is not None) == attach
    scene.add_user_template_fragment(scene.template_context)
    atoms = window.state_manager.data.atoms
    assert len(atoms) == 5
    placed_id = next(i for i, atom in atoms.items() if atom["symbol"] == "C")
    assert atoms[placed_id].get("rs") == label
    assert scene.atom_items[placed_id].rs_label == label
    assert len(window.state_manager.data.bonds) == 4
