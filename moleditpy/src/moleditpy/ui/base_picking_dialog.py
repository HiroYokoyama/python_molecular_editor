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
from PyQt6.QtWidgets import QDialog
from PyQt6.QtCore import Qt

try:
    from .dialog_3d_picking_mixin import Dialog3DPickingMixin
except ImportError:
    from moleditpy.ui.dialog_3d_picking_mixin import Dialog3DPickingMixin


class BasePickingDialog(Dialog3DPickingMixin, QDialog):
    """
    Base class for any dialog requiring 3D atom picking.
    Provides standard cleanup and event handling for picking filters and labels.
    """

    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self._molecule_modified = False  # Track if any modifications were made during this session

    def keyPressEvent(self, event):
        """Standard keyboard handler: Enter/Return triggers 'Apply'."""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if hasattr(self, "apply_button") and self.apply_button.isEnabled():
                # Call the apply method (must be implemented by subclass or connected)
                self.apply_button.click()
            event.accept()
        else:
            QDialog.keyPressEvent(self, event)

    def closeEvent(self, event):
        """Cleanup on window close."""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)

    def reject(self):
        """Cleanup on cancel."""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()

    def accept(self):
        """Cleanup on OK."""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()

    def _update_molecule_geometry(self, positions):
        """
        Update the molecule's conformer and the 3D position cache, then redraw.
        :param positions: A numpy array or dictionary of all atom positions.
        """
        from rdkit import Geometry
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
            if isinstance(positions, dict):
                for i, p in positions.items():
                    self.main_window.view_3d_manager.atom_positions_3d[i] = p
            else:
                self.main_window.view_3d_manager.atom_positions_3d[:] = positions[:]
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
            logging.error(f"REPORT ERROR: Missing attribute 'update_chiral_labels' on object")

    def _push_undo(self):
        """Centralized undo logic to push current state to the undo stack."""
        if hasattr(self.main_window, "state_manager"):
            self.main_window.edit_actions_manager.push_undo_state()
            self._molecule_modified = False
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(f"REPORT ERROR: Missing attribute 'state_manager' on self.main_window")

    def done(self, result):
        """Override done to push a final undo state if the molecule was modified."""
        if self._molecule_modified:
            self._push_undo()
        super().done(result)
