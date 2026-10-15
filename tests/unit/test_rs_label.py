"""R/S atom labels: core priority over wedges, wedge alignment, scene tool."""

from unittest.mock import MagicMock, patch

import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QPainter, QTransform
from PyQt6.QtWidgets import QStyleOptionGraphicsItem
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from moleditpy.core import molecular_data as md
from moleditpy.core import stereo_check as sc
from moleditpy.core.molecular_data import (
    MolecularData,
    apply_rs_labels,
    rs_entry,
)
from moleditpy.core.stereo_check import (
    align_wedges_to_label,
    drawn_chirality,
    wedge_for_label,
)
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.molecular_scene_handler import next_rs_label
from moleditpy.ui.molecule_scene import MoleculeScene

# C0 with F1 (right), Cl2 (lower left), Br3 (upper left), O4 (down): one
# stereocenter, drawn in scene pixels (y down).
_CENTER_ATOMS = [
    ("C", (0, 0)),
    ("F", (50, 0)),
    ("Cl", (-25, 43)),
    ("Br", (-25, -43)),
    ("O", (0, 60)),
]


def _center_data():
    data = MolecularData()
    ids = [data.add_atom(sym, pos) for sym, pos in _CENTER_ATOMS]
    for other in ids[1:]:
        data.add_bond(ids[0], other)
    return data


def _cip_of_mol_block(block, idx=0):
    mol = Chem.MolFromMolBlock(block)
    rdCIPLabeler.AssignCIPLabels(mol)
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


# ---------------------------------------------------------------------------
# core: labels override wedges
# ---------------------------------------------------------------------------


def test_rs_entry_only_for_valid_labels():
    assert rs_entry({"rs": "R"}) == {"rs": "R"}
    assert rs_entry({"rs": "S"}) == {"rs": "S"}
    assert rs_entry({"rs": "x"}) == {}
    assert rs_entry({}) == {}


@pytest.mark.parametrize("label", ["R", "S"])
def test_label_alone_sets_configuration(label):
    data = _center_data()
    assert drawn_chirality(data) == {}
    data.atoms[0]["rs"] = label
    assert drawn_chirality(data) == {0: label}
    assert drawn_chirality(data, use_rs_labels=False) == {}
    # Written as wedges, so a MOL block (export, 3D conversion) carries it.
    assert _cip_of_mol_block(data.to_mol_block()) == label


def test_label_overrides_a_contradicting_wedge():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    drawn = drawn_chirality(data)[0]
    wanted = "R" if drawn == "S" else "S"
    data.atoms[0]["rs"] = wanted
    assert drawn_chirality(data) == {0: wanted}
    assert drawn_chirality(data, use_rs_labels=False) == {0: drawn}
    assert _cip_of_mol_block(data.to_mol_block()) == wanted
    mol = data.to_rdkit_mol(use_2d_stereo=False, use_rs_labels=False)
    assert mol is not None


def test_label_that_agrees_with_the_wedge_changes_nothing():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    drawn = drawn_chirality(data)[0]
    data.atoms[0]["rs"] = drawn
    with patch.object(md, "_rdkit_rewedge") as rewedge:
        assert drawn_chirality(data) == {0: drawn}
    rewedge.assert_not_called()


def test_label_on_a_non_stereocenter_is_ignored():
    data = _center_data()
    data.atoms[1]["rs"] = "R"  # fluorine
    assert drawn_chirality(data) == {}


def test_apply_rs_labels_skips_empty_and_conformerless():
    mol = Chem.MolFromSmiles("FC(Cl)(Br)O")
    apply_rs_labels(mol, {})
    apply_rs_labels(mol, {1: "R"})  # no conformer: nothing to wedge against
    assert all(b.GetBondDir() == Chem.BondDir.NONE for b in mol.GetBonds())


def test_tag_for_label_none_without_cip_code():
    mol = _center_data().to_rdkit_mol(use_rs_labels=False)
    tags = md._rdkit_tags(mol)
    with (
        patch.object(md.Chem, "AssignStereochemistry"),
        patch.object(md.rdCIPLabeler, "AssignCIPLabels"),
    ):
        assert md._tag_for_label(mol, 0, "R", tags) is None


def test_apply_rs_labels_falls_back_to_picking_wedges():
    data = _center_data()
    data.atoms[0]["rs"] = "R"
    with patch.object(md, "_rdkit_rewedge"):  # first strategy does nothing
        assert drawn_chirality(data) == {0: "R"}


def test_apply_rs_labels_warns_when_no_drawing_holds(caplog):
    data = _center_data()
    data.atoms[0]["rs"] = "R"
    with (
        patch.object(md, "_rdkit_rewedge"),
        patch.object(md, "_clear_and_pick_wedges"),
        caplog.at_level("WARNING"),
    ):
        data.to_rdkit_mol()
    assert "Could not draw R/S labels" in caplog.text


def test_apply_rs_labels_logs_rdkit_errors(caplog):
    data = _center_data()
    data.atoms[0]["rs"] = "R"
    with (
        patch.object(md, "_rdkit_tags", side_effect=RuntimeError("boom")),
        caplog.at_level("WARNING"),
    ):
        data.to_rdkit_mol()
    assert "Applying R/S labels failed" in caplog.text


# ---------------------------------------------------------------------------
# core: aligning existing wedges with a label
# ---------------------------------------------------------------------------


def test_align_without_wedges_keeps_plain_bonds():
    assert align_wedges_to_label(_center_data(), 0, "R") == {}


def test_align_rejects_unknown_label():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    assert align_wedges_to_label(data, 0, "X") == {}


def test_align_when_drawing_already_agrees():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    assert align_wedges_to_label(data, 0, drawn_chirality(data)[0]) == {}


def test_align_swaps_wedges_to_hashes():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    data.bonds[(0, 2)]["stereo"] = 1
    drawn = drawn_chirality(data)[0]
    wanted = "R" if drawn == "S" else "S"
    before = {k: dict(v) for k, v in data.bonds.items()}
    changes = align_wedges_to_label(data, 0, wanted)
    assert changes == {(0, 1): 2, (0, 2): 2}
    assert data.bonds == before  # data left unchanged
    for key, stereo in changes.items():
        data.bonds[key]["stereo"] = stereo
    assert drawn_chirality(data, use_rs_labels=False) == {0: wanted}


def test_align_redraws_one_bond_when_swapping_does_not_help():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    data.bonds[(0, 2)]["stereo"] = 2
    found = sc.WedgeForLabel((0, 1), 2, ((0, 2),))
    with (
        patch.object(sc, "drawn_chirality", return_value={}),
        patch.object(sc, "wedge_for_label", return_value=found) as search,
    ):
        changes = align_wedges_to_label(data, 0, "R")
    assert changes == {(0, 2): 0, (0, 1): 2}
    search.assert_called_once_with(data, (0, 1), "R", centers=[0])
    assert data.bonds[(0, 1)]["stereo"] == 1  # data restored after the swap


def test_align_gives_up_when_nothing_reads_as_label():
    data = _center_data()
    data.bonds[(0, 1)]["stereo"] = 1
    with (
        patch.object(sc, "drawn_chirality", return_value={}),
        patch.object(sc, "wedge_for_label", return_value=None),
    ):
        assert align_wedges_to_label(data, 0, "R") == {}


def test_wedge_for_label_guards():
    data = _center_data()
    assert wedge_for_label(data, (9, 9), "R") is None
    assert wedge_for_label(data, (0, 1), "X") is None
    # An end that is not part of the bond is skipped.
    assert wedge_for_label(data, (0, 1), "R", centers=[7]) is None


def test_wedge_for_label_clears_conflicting_hashes():
    data = _center_data()
    data.bonds[(0, 2)]["stereo"] = 2
    data.bonds[(0, 3)]["stereo"] = 2
    drawn = drawn_chirality(data, use_rs_labels=False)[0]
    wanted = "R" if drawn == "S" else "S"
    found = wedge_for_label(data, (0, 1), wanted)
    assert found is not None and found.key == (0, 1)
    assert found.cleared  # the hashes that fought the label


def test_wedge_for_label_logs_and_restores_on_error(caplog):
    data = _center_data()
    before = list(data.bonds.items())
    with (
        patch.object(sc, "drawn_chirality", side_effect=ValueError("bad")),
        caplog.at_level("WARNING"),
    ):
        assert wedge_for_label(data, (0, 1), "R") is None
    assert list(data.bonds.items()) == before
    assert "R/S wedge search failed" in caplog.text


# ---------------------------------------------------------------------------
# scene: the R/S tool
# ---------------------------------------------------------------------------


def _scene(host):
    scene = MoleculeScene(host.data, host)
    view = MagicMock()
    view.transform.return_value = QTransform()
    view.mapToScene.return_value = QPointF(0, 0)
    view.mapFromGlobal.return_value = QPointF(0, 0)
    scene.views = MagicMock(return_value=[view])
    return scene


def _center_scene(host):
    scene = _scene(host)
    items = [
        scene.atom_items[scene.create_atom(sym, QPointF(*pos))]
        for sym, pos in _CENTER_ATOMS
    ]
    for other in items[1:]:
        scene.create_bond(items[0], other, bond_order=1, bond_stereo=0)
    return scene, items


def _event(pos, button=Qt.MouseButton.LeftButton):
    event = MagicMock()
    event.scenePos.return_value = pos
    event.button.return_value = button
    return event


def test_next_rs_label_cycle():
    assert next_rs_label(None) == "R"
    assert next_rs_label("R") == "S"
    assert next_rs_label("S") is None


def test_apply_rs_label_sets_and_clears(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    assert scene.apply_rs_label(center, "R") is True
    assert scene.data.atoms[center.atom_id]["rs"] == "R"
    assert center.rs_label == "R"
    assert scene.apply_rs_label(center, "R") is False  # unchanged
    assert scene.apply_rs_label(center, "S") is True
    assert scene.data.atoms[center.atom_id]["rs"] == "S"
    assert scene.apply_rs_label(center, None) is True
    assert "rs" not in scene.data.atoms[center.atom_id]
    assert center.rs_label is None


def test_apply_rs_label_refuses_non_stereocenter(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    fluorine = items[1]
    assert scene.apply_rs_label(fluorine, "R") is False
    assert "rs" not in scene.data.atoms[fluorine.atom_id]
    mock_parser_host.statusBar().showMessage.assert_called()


def test_apply_rs_label_refusal_keeps_previous_label(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    scene.apply_rs_label(center, "R")
    with patch(
        "moleditpy.ui.molecular_scene_handler.drawn_chirality", return_value={}
    ):
        assert scene.apply_rs_label(center, "S") is False
    assert scene.data.atoms[center.atom_id]["rs"] == "R"


def test_apply_rs_label_ignores_unknown_atoms(mock_parser_host):
    scene = _scene(mock_parser_host)
    stray = MagicMock(atom_id=999)
    assert scene.apply_rs_label(stray, "R") is False
    assert scene.apply_rs_label(None, "R") is False
    scene.set_atom_rs_label(stray, "R")  # no record: nothing to do


def test_apply_rs_label_redraws_existing_wedges(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center, fluorine = items[0], items[1]
    key = (center.atom_id, fluorine.atom_id)
    scene.data.bonds[key]["stereo"] = 1
    scene.find_bond_between(center, fluorine).set_stereo(1)
    drawn = drawn_chirality(scene.data)[center.atom_id]
    wanted = "R" if drawn == "S" else "S"

    assert scene.apply_rs_label(center, wanted)

    assert scene.data.bonds[key]["stereo"] == 2
    assert scene.find_bond_between(center, fluorine).stereo == 2
    assert drawn_chirality(scene.data, use_rs_labels=False) == {
        center.atom_id: wanted
    }


def test_set_bond_stereo_by_key_rekeys_plain_bonds(mock_parser_host):
    scene = _scene(mock_parser_host)
    a = scene.atom_items[scene.create_atom("C", QPointF(0, 0))]
    b = scene.atom_items[scene.create_atom("C", QPointF(50, 0))]
    scene.create_bond(b, a, bond_order=1, bond_stereo=1)  # key (b, a), b > a
    bond = scene.find_bond_between(a, b)
    key = (b.atom_id, a.atom_id)
    assert key in scene.data.bonds

    scene._set_bond_stereo_by_key(key, 0)

    sorted_key = (a.atom_id, b.atom_id)
    assert sorted_key in scene.data.bonds and key not in scene.data.bonds
    assert scene.bond_items[sorted_key] is bond
    assert bond.atom1 is a and bond.stereo == 0


def test_set_bond_stereo_by_key_without_item_or_data(mock_parser_host):
    scene = _scene(mock_parser_host)
    scene._set_bond_stereo_by_key((1, 2), 1)  # no such bond
    scene.data.bonds[(5, 3)] = {"order": 1, "stereo": 1}
    scene._set_bond_stereo_by_key((5, 3), 0)
    assert scene.data.bonds[(3, 5)]["stereo"] == 0


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_click_cycles_label(mock_drag, mock_release, mock_press, mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    scene.mode = "rs_label"
    push = mock_parser_host.edit_actions_manager.push_undo_state
    push.reset_mock()
    seen = []
    with (
        patch.object(MoleculeScene, "itemAt", return_value=center),
        patch.object(MoleculeScene, "find_atom_near", return_value=center),
    ):
        for _ in range(3):
            event = _event(QPointF(0, 0))
            scene.mousePressEvent(event)
            scene.mouseReleaseEvent(event)
            seen.append(center.rs_label)
    assert seen == ["R", "S", None]
    assert push.call_count == 3


@patch("PyQt6.QtWidgets.QGraphicsScene.mousePressEvent")
@patch("PyQt6.QtWidgets.QGraphicsScene.mouseReleaseEvent")
@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_click_on_non_stereocenter_pushes_no_undo(
    mock_drag, mock_release, mock_press, mock_parser_host
):
    scene, items = _center_scene(mock_parser_host)
    scene.mode = "rs_label"
    push = mock_parser_host.edit_actions_manager.push_undo_state
    push.reset_mock()
    with (
        patch.object(MoleculeScene, "itemAt", return_value=items[1]),
        patch.object(MoleculeScene, "find_atom_near", return_value=items[1]),
    ):
        event = _event(QPointF(50, 0))
        scene.mousePressEvent(event)
        scene.mouseReleaseEvent(event)
    assert items[1].rs_label is None
    push.assert_not_called()


def test_double_click_advances_label(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    scene.mode = "rs_label"
    with patch.object(MoleculeScene, "itemAt", return_value=center):
        scene.mouseDoubleClickEvent(_event(QPointF(0, 0)))
        assert center.rs_label == "R"
        with patch.object(scene, "apply_rs_label", return_value=False):
            scene.mouseDoubleClickEvent(_event(QPointF(0, 0)))
    assert center.rs_label == "R"


@patch("moleditpy.ui.molecule_scene.QApplication.startDragDistance", return_value=1000)
def test_right_click_clears_label_without_deleting(mock_drag, mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    scene.apply_rs_label(center, "S")
    scene.mode = "rs_label"
    with patch.object(MoleculeScene, "itemAt", return_value=center):
        scene.mousePressEvent(_event(QPointF(0, 0), Qt.MouseButton.RightButton))
    assert center.rs_label is None
    assert center.atom_id in scene.data.atoms
    # Right-click on an unlabelled atom in R/S mode does not delete it either.
    with patch.object(MoleculeScene, "itemAt", return_value=items[1]):
        scene.mousePressEvent(_event(QPointF(50, 0), Qt.MouseButton.RightButton))
    assert items[1].atom_id in scene.data.atoms


def _key_event(key, modifiers=Qt.KeyboardModifier.NoModifier):
    event = MagicMock()
    event.key.return_value = key
    event.modifiers.return_value = modifiers
    return event


def test_r_key_cycles_hovered_atom(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    with (
        patch.object(MoleculeScene, "find_atom_near", return_value=center),
        patch.object(MoleculeScene, "itemAt", return_value=center),
    ):
        scene.keyPressEvent(_key_event(Qt.Key.Key_R))
        assert center.rs_label == "R"
        scene.keyPressEvent(_key_event(Qt.Key.Key_R))
        assert center.rs_label == "S"
        # Non-stereocenter: event consumed, nothing changes.
        with patch.object(scene, "apply_rs_label", return_value=False):
            event = _key_event(Qt.Key.Key_R)
            scene.keyPressEvent(event)
            event.accept.assert_called_once()
    assert center.symbol == "C"


# ---------------------------------------------------------------------------
# AtomItem drawing
# ---------------------------------------------------------------------------


def test_atom_item_label_geometry_and_paint(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    assert center._rs_label_rect() is None
    plain_rect = center.boundingRect()
    scene.apply_rs_label(center, "R")
    rect = center._rs_label_rect()
    assert rect is not None
    assert center.boundingRect().contains(rect)
    assert center.boundingRect() != plain_rect

    painter = MagicMock(spec=QPainter)
    center.paint(painter, QStyleOptionGraphicsItem(), None)
    assert painter.drawPath.call_count >= 2  # outline + text

    # Outline color falls back to white when the background is unreadable.
    fake_scene = MagicMock()
    fake_scene.get_setting.side_effect = lambda key, default=None: default
    fake_scene.backgroundBrush.side_effect = RuntimeError("gone")
    with patch.object(AtomItem, "scene", return_value=fake_scene):
        center._paint_rs_label(MagicMock(spec=QPainter))


def test_atom_item_label_without_bonds_sits_upper_right(mock_parser_host):
    scene = _scene(mock_parser_host)
    lone = scene.atom_items[scene.create_atom("C", QPointF(0, 0))]
    lone.rs_label = "S"
    rect = lone._rs_label_rect()
    assert rect.center().x() > 0 and rect.center().y() < 0


def test_atom_item_label_skips_bad_partners(mock_parser_host):
    scene, items = _center_scene(mock_parser_host)
    center = items[0]
    center.rs_label = "R"
    broken = MagicMock()
    broken.atom1 = center
    broken.atom2.pos.side_effect = RuntimeError("deleted")
    missing = MagicMock(atom1=center, atom2=None)
    stacked = MagicMock(atom1=center)
    stacked.atom2.pos.return_value = center.pos()
    center.bonds = [broken, missing, stacked]
    assert center._rs_label_rect() is not None
