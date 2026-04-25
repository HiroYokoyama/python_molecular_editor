#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QCloseEvent, QKeyEvent
from PyQt6.QtWidgets import QDialog, QWidget
from rdkit import Chem, Geometry

try:
    from .dialog_3d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from moleditpy_linux.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin

if TYPE_CHECKING:
    from .main_window import MainWindow


class BasePickingDialog(Dialog3DPickingMixin, QDialog):
    """
    Base class for any dialog requiring 3D atom picking.
    Provides standard cleanup and event handling for picking filters and labels.
    """

    def __init__(
        self,
        mol: Chem.Mol,
        main_window: "MainWindow",
        parent: Optional[QWidget] = None,
    ) -> None:
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self._molecule_modified = (
            False  # Track if any modifications were made during this session
        )

    def keyPressEvent(self, event: Optional[QKeyEvent]) -> None:
        """Standard keyboard handler: Enter/Return triggers 'Apply'."""
        if event is None:
            return
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if hasattr(self, "apply_button") and self.apply_button.isEnabled():
                # Call the apply method (must be implemented by subclass or connected)
                self.apply_button.click()
            event.accept()
        else:
            QDialog.keyPressEvent(self, event)

    def closeEvent(self, event: Optional[QCloseEvent]) -> None:
        """Cleanup on window close."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self) -> None:
        """Cleanup on cancel."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()

    def accept(self) -> None:
        """Cleanup on OK."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()

    def _update_molecule_geometry(
        self, positions: Union[np.ndarray, dict[int, np.ndarray]]
    ) -> None:
        """
        Update the molecule's conformer and the 3D position cache, then redraw.
        :param positions: A numpy array or dictionary of all atom positions.
        """
        conf = self.mol.GetConformer()
        num_atoms = conf.GetNumAtoms()

        # 1. Update RDKit Conformer
        if isinstance(positions, dict):
            for i, p in positions.items():
                conf.SetAtomPosition(
                    i, Geometry.Point3D(float(p[0]), float(p[1]), float(p[2]))
                )
        else:
            for i in range(num_atoms):
                p = positions[i]
                conf.SetAtomPosition(
                    i, Geometry.Point3D(float(p[0]), float(p[1]), float(p[2]))
                )

        # 2. Update 3D Visualization cache
        try:
            cache = self.main_window.view_3d_manager.atom_positions_3d
            if cache is not None:
                if isinstance(positions, dict):
                    for i, p in positions.items():
                        cache[i] = p
                else:
                    cache[:] = positions[:]
        except (AttributeError, ValueError, TypeError, IndexError):
            # If for some reason the cache is incompatible, draw_molecule_3d below will rebuild it
            pass

        # 3. Redraw
        self.main_window.view_3d_manager.draw_molecule_3d(self.mol)
        self._molecule_modified = True

        # 4. Refresh chiral/cis-trans labels if applicable
        if hasattr(self.main_window.view_3d_manager, "update_chiral_labels"):
            self.main_window.view_3d_manager.update_chiral_labels()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'update_chiral_labels' on object"
            )

    def _push_undo(self) -> None:
        """Centralized undo logic to push current state to the undo stack."""
        if hasattr(self.main_window, "state_manager"):
            self.main_window.edit_actions_manager.push_undo_state()
            self._molecule_modified = False
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'state_manager' on self.main_window"
            )

    def done(self, result: int) -> None:
        """Override done to push a final undo state if the molecule was modified."""
        if self._molecule_modified:
            self._push_undo()
        super().done(result)
