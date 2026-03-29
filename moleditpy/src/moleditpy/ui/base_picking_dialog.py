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
