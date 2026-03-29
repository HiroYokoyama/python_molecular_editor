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
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMessageBox

try:
    from .base_picking_dialog import BasePickingDialog
except ImportError:
    from moleditpy.ui.base_picking_dialog import BasePickingDialog


class GeometryBaseDialog(BasePickingDialog):
    """
    Base class for dialogs that adjust a single numerical value (Length, Angle, Dihedral).
    Handles synchronization between QLineEdit and QSlider and real-time 3D updates.
    """

    def __init__(self, mol, main_window, parent=None):
        super().__init__(mol, main_window, parent)
        self._slider_dragging = False
        self._snapshot_positions = None

    def _sync_input_to_slider(self, val, slider, scale=1.0, wrap=False):
        """Sync a numerical value from input text to the slider."""
        try:
            f_val = float(val)
            if wrap:
                f_val = (f_val + 180) % 360 - 180
            
            slider.blockSignals(True)
            slider.setValue(int(round(f_val * scale)))
            slider.blockSignals(False)
        except (ValueError, TypeError):
            pass

    def on_slider_pressed(self):
        """Prepare for a slider drag operation by saving a geometry snapshot."""
        if not self._is_selection_complete():
            return
            
        self._slider_dragging = True

        # Capture geometry snapshot to ensure stable axes during rotation/dragging
        self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()

    def on_slider_released(self):
        """Finalize a slider drag operation."""
        self._slider_dragging = False
        # Snapshot is usually kept until selection changes to preserve turn direction,
        # but subclasses can override this behavior if needed.
        self.main_window.view_3d_manager.draw_molecule_3d(self.mol)
        if hasattr(self.main_window.view_3d_manager, "update_chiral_labels"):
            self.main_window.view_3d_manager.update_chiral_labels()

    def on_slider_value_changed_click(self, value, input_box, scale=1.0):
        """
        Handle a discrete value change (e.g., clicking on the slider rail).
        Triggers a geometric update and saves undo state.
        """
        if self._slider_dragging:
            return  # Handled by on_slider_moved
            
        if not self._is_selection_complete():
            return

        # Ensure snapshot exists for consistency
        if self._snapshot_positions is None:
            self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()

        input_box.blockSignals(True)
        input_box.setText(f"{value / scale:.3f}")
        input_box.blockSignals(False)
        
        self.apply_geometry_update(float(value / scale))
        
        if hasattr(self.main_window.view_3d_manager, "update_chiral_labels"):
            self.main_window.view_3d_manager.update_chiral_labels()

    def on_slider_moved_realtime(self, value, input_box, scale=1.0):
        """Update geometry in real-time as the slider is dragged."""
        if not self._is_selection_complete():
            return

        input_box.blockSignals(True)
        input_box.setText(f"{value / scale:.3f}")
        input_box.blockSignals(False)

        self.apply_geometry_update(float(value / scale))

    def _is_selection_complete(self):
        """Must be implemented by subclass to check if enough atoms are picked."""
        raise NotImplementedError

    def apply_geometry_update(self, new_value):
        """Must be implemented by subclass to perform the actual RdKit/Numpy update."""
        raise NotImplementedError
