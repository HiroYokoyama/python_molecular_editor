"""AtomItem label helpers shared by paint(), visual_rect() and the clip path."""

from types import SimpleNamespace

import pytest
from PyQt6.QtCore import QPointF

from moleditpy.ui.atom_item import AtomItem


def _atom(symbol="O", charge=0, radical=0, h=0, partner_dx=()):
    atom = AtomItem(0, symbol, QPointF(0, 0), charge=charge, radical=radical)
    atom.implicit_h_count = h
    for dx in partner_dx:
        atom.bonds.append(
            SimpleNamespace(atom1=atom, atom2=AtomItem(1, "C", QPointF(dx, 0)))
        )
    return atom


@pytest.mark.parametrize("charge, text", [(1, "+"), (-1, "-"), (2, "2+"), (-3, "3-")])
def test_charge_text(app, charge, text):
    assert _atom(charge=charge)._charge_text() == text


def test_hydrogen_part(app):
    assert _atom(h=0)._hydrogen_part() == ""
    assert _atom(h=1)._hydrogen_part() == "H"
    assert _atom(h=2)._hydrogen_part() == "H₂"
    # A bonded neutral carbon is skeletal: no label at all.
    assert _atom("C", h=3, partner_dx=[20])._hydrogen_part() == ""


def test_label_flips_away_from_bonds(app):
    assert _atom(h=1, partner_dx=[30])._label_flipped("H") is True
    assert _atom(h=1, partner_dx=[-30])._label_flipped("H") is False
    assert _atom(h=1, partner_dx=[30])._label_flipped("") is False


def test_label_flip_skips_missing_partner(app):
    atom = _atom(h=1, partner_dx=[-30])
    atom.bonds.append(SimpleNamespace(atom1=atom, atom2=None))
    assert atom._label_flipped("H") is False


def test_label_flip_skips_partner_without_position(app):
    """A partner whose pos() is None, or whose wrapper raises, is ignored."""
    atom = _atom(h=1, partner_dx=[-30])
    no_pos = SimpleNamespace(pos=lambda: None)
    atom.bonds.append(SimpleNamespace(atom1=atom, atom2=no_pos))

    def _gone():
        raise RuntimeError("wrapped C/C++ object has been deleted")

    atom.bonds.append(SimpleNamespace(atom1=atom, atom2=SimpleNamespace(pos=_gone)))
    # Only the real partner, on the left, counts: no flip.
    assert atom._label_flipped("H") is False
