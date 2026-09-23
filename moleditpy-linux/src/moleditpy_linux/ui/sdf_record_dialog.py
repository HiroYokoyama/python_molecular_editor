#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

from typing import List, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..utils.sdf_records import SdfRecord


class SdfRecordDialog(QDialog):
    """Pick one record of a multi-molecule SDF file.

    Mirrors the Open Babel Conversion Tool plugin's selector. Records RDKit
    could not read are listed, greyed out, so the numbering matches the file.
    """

    def __init__(
        self, records: List[SdfRecord], parent: Optional[QWidget] = None
    ) -> None:
        """Build the list for *records*, preselecting the first readable one."""
        super().__init__(parent)
        self.records = records
        self.setWindowTitle("Select Molecule")
        self.resize(400, 300)

        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(f"Found {len(records)} molecules. Please select one:", self)
        )

        self.list_widget = QListWidget(self)
        self.list_widget.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        for record in records:
            item = QListWidgetItem(self._item_text(record))
            if not record.readable:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.list_widget.addItem(item)
        self.list_widget.itemDoubleClicked.connect(lambda _item: self.accept())
        layout.addWidget(self.list_widget)

        first = next((r.index for r in records if r.readable), None)
        if first is not None:
            self.list_widget.setCurrentRow(first)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @staticmethod
    def _item_text(record: SdfRecord) -> str:
        """List text: number, name and formula, or a note that it is unreadable."""
        number = f"{record.index + 1}: {record.display_name()}"
        if not record.readable:
            return f"{number} (unreadable)"
        return f"{number} ({record.formula})"

    def selected_record(self) -> Optional[SdfRecord]:
        """The chosen readable record, or None when nothing usable is selected."""
        row = self.list_widget.currentRow()
        if 0 <= row < len(self.records) and self.records[row].readable:
            return self.records[row]
        return None
