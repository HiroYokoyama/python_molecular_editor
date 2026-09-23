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

from ..core.stereo_check import ChiralityMismatch


class ChiralityWarningDialog(QDialog):
    """Non-modal, always-on-top warning that the 3D result has wrong chirality.

    While it is open the 3D view shows chiral labels, whatever the View menu
    says, with the wrong centers in red; closing it hands the labels back to
    the menu setting.
    """

    def __init__(
        self,
        main_window: Any,
        mismatches: List[ChiralityMismatch],
        total_centers: int,
        parent: Optional[QWidget] = None,
    ) -> None:
        """Build the warning for *mismatches* out of *total_centers* drawn centers."""
        super().__init__(parent)
        self.main_window = main_window
        self.mismatches = mismatches
        self.setWindowTitle("Chirality Check")
        self.setModal(False)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._build_ui(total_centers)
        self.finished.connect(self._restore_chiral_labels)
        self._force_chiral_labels()

    def _build_ui(self, total_centers: int) -> None:
        """Lay out the summary, the per-atom table and the close button."""
        layout = QVBoxLayout(self)

        count = len(self.mismatches)
        summary = QLabel(
            f"<b>The 3D structure is wrong.</b><br>"
            f"{count} of {total_centers} stereocenter(s) drawn in 2D have a "
            f"different configuration in 3D."
        )
        summary.setWordWrap(True)
        layout.addWidget(summary)

        table = QTableWidget(count, 4, self)
        table.setHorizontalHeaderLabels(["Atom", "Index", "Drawn (2D)", "3D"])
        table.verticalHeader().setVisible(False)  # type: ignore[union-attr]
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        for row, m in enumerate(self.mismatches):
            index = "" if m.rdkit_index is None else str(m.rdkit_index)
            cells = [f"{m.symbol} (ID {m.atom_id})", index, m.drawn, m.actual or "none"]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                table.setItem(row, col, item)
        header = table.horizontalHeader()
        if header is not None:
            header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(table)

        hint = QLabel(
            "Please flip these stereocenters manually in the 3D structure. "
            "While this window is open, the 3D view shows chiral labels: "
            "wrong centers in red, correct ones in the usual color."
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
        """Show chiral labels, wrong centers marked, for as long as the dialog is open."""
        view_3d = self._view_3d()
        if view_3d is None:
            return
        view_3d.chirality_mismatches = {
            m.rdkit_index: m.actual or "?"
            for m in self.mismatches
            if m.rdkit_index is not None
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
        view_3d.show_chiral_labels = (
            bool(action.isChecked()) if action is not None else False
        )
        # Always redraw: even with labels left on, the red marks must go.
        self._redraw(view_3d)
