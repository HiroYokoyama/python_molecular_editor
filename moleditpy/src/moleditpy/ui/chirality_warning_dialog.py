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

import logging
from typing import Any, List, Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..core.stereo_check import ChiralityMismatch, EZMismatch


class ChiralityWarningDialog(QDialog):
    """Non-modal, always-on-top warning that the 3D result has wrong stereo.

    Covers stereocenters (R/S) and labelled double bonds (E/Z). While it is
    open the 3D view shows chiral and E/Z labels, whatever the View menu says,
    with the wrong ones in red; closing it hands the labels back to the menu
    setting.
    """

    def __init__(
        self,
        main_window: Any,
        mismatches: List[ChiralityMismatch],
        total_centers: int,
        parent: Optional[QWidget] = None,
        ez_mismatches: Optional[List[EZMismatch]] = None,
        total_double_bonds: int = 0,
    ) -> None:
        """Build the warning for *mismatches* and *ez_mismatches*.

        *total_centers* and *total_double_bonds* count what the drawing
        specifies, for the "N of M" summary.
        """
        super().__init__(parent)
        self.main_window = main_window
        self.mismatches = mismatches
        self.ez_mismatches = ez_mismatches or []
        self.setWindowTitle("Stereochemistry Check")
        self.setModal(False)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._build_ui(total_centers + total_double_bonds)
        self.finished.connect(self._restore_chiral_labels)
        self._force_chiral_labels()

    def _rows(self) -> List[List[str]]:
        """Table rows: stereocenters first, then double bonds."""
        rows = [
            [
                f"{m.symbol} (ID {m.atom_id})",
                "" if m.rdkit_index is None else str(m.rdkit_index),
                m.drawn,
                m.actual or "none",
            ]
            for m in self.mismatches
        ]
        for b in self.ez_mismatches:
            (s1, s2), (i1, i2) = b.symbols, b.atom_ids
            index = (
                ""
                if b.rdkit_atom_indices is None
                else "=".join(str(i) for i in b.rdkit_atom_indices)
            )
            rows.append(
                [f"{s1}={s2} (ID {i1}={i2})", index, b.drawn, b.actual or "none"]
            )
        return rows

    def _build_ui(self, total: int) -> None:
        """Lay out the summary, the per-item table and the close button."""
        layout = QVBoxLayout(self)

        rows = self._rows()
        summary = QLabel(
            f"<b>The 3D structure is wrong.</b><br>"
            f"{len(rows)} of {total} stereo element(s) drawn in 2D (R/S centers, "
            f"E/Z double bonds) have a different configuration in 3D."
        )
        summary.setWordWrap(True)
        layout.addWidget(summary)

        table = QTableWidget(len(rows), 4, self)
        table.setHorizontalHeaderLabels(["Atom / Bond", "Index", "Drawn (2D)", "3D"])
        table.verticalHeader().setVisible(False)  # type: ignore[union-attr]
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        for row, cells in enumerate(rows):
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                table.setItem(row, col, item)
        header = table.horizontalHeader()
        if header is not None:
            header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(table)

        hint = QLabel(
            "Please flip these manually in the 3D structure. While this window "
            "is open, the 3D view shows chiral and E/Z labels: wrong ones in "
            "red, correct ones in the usual color."
        )
        hint.setWordWrap(True)
        layout.addWidget(hint)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close, self)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.resize(460, 320)

    def _view_3d(self) -> Any:
        """Return the main window's 3D view manager, or None."""
        return getattr(self.main_window, "view_3d_manager", None)

    def _redraw(self, view_3d: Any) -> None:
        """Redraw the current 3D molecule so the label setting takes effect."""
        mol = getattr(view_3d, "current_mol", None)
        if mol is not None:
            try:
                view_3d.draw_molecule_3d(mol)
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                logging.warning("Redraw for chiral labels failed: %s", e)

    def _force_chiral_labels(self) -> None:
        """Show stereo labels, wrong ones marked, for as long as the dialog is open."""
        view_3d = self._view_3d()
        if view_3d is None:
            return
        view_3d.chirality_mismatches = {
            m.rdkit_index: m.actual or "?"
            for m in self.mismatches
            if m.rdkit_index is not None
        }
        view_3d.ez_mismatches = {
            b.rdkit_bond_index: b.actual or "?"
            for b in self.ez_mismatches
            if b.rdkit_bond_index is not None
        }
        view_3d.show_chiral_labels = True
        self._redraw(view_3d)

    def _restore_chiral_labels(self) -> None:
        """Hand chiral label visibility back to the View menu setting."""
        view_3d = self._view_3d()
        if view_3d is None:
            return
        action = getattr(
            getattr(self.main_window, "init_manager", None),
            "toggle_chiral_action",
            None,
        )
        view_3d.chirality_mismatches = {}
        view_3d.ez_mismatches = {}
        view_3d.show_chiral_labels = (
            bool(action.isChecked()) if action is not None else False
        )
        # Always redraw: even with labels left on, the red marks must go.
        self._redraw(view_3d)
