"""The shared periodic-table layout used by the element picker and CPK Colors."""

from PyQt6.QtGui import QColor

from moleditpy.ui.periodic_table_dialog import element_button_style
from moleditpy.utils.constants import PERIODIC_TABLE_LAYOUT


def test_layout_has_every_element_once_in_its_own_cell():
    symbols = [s for s, _r, _c in PERIODIC_TABLE_LAYOUT]
    cells = [(r, c) for _s, r, c in PERIODIC_TABLE_LAYOUT]
    assert len(symbols) == 118
    assert len(set(symbols)) == 118
    assert len(set(cells)) == 118


def test_button_label_contrasts_with_fill():
    assert "color: white" in element_button_style(QColor("#000080"))
    assert "color: black" in element_button_style(QColor("#FFFF00"))
