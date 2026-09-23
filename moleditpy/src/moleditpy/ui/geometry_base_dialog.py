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
from typing import TYPE_CHECKING, Any, Optional, Union

import numpy as np
from PyQt6.QtWidgets import QLineEdit, QSlider, QWidget
from rdkit import Chem

from .base_picking_dialog import BasePickingDialog

if TYPE_CHECKING:
    from .main_window import MainWindow


class GeometryBaseDialog(BasePickingDialog):
    """
    Base class for dialogs that adjust a single numerical value (Length, Angle, Dihedral).
    Handles synchronization between QLineEdit and QSlider and real-time 3D updates.
    """

    def __init__(
        self,
        mol: Chem.Mol,
        main_window: "MainWindow",
        parent: Optional[QWidget] = None,
    ) -> None:
        """Initialize base geometry manipulation dialog."""
        super().__init__(mol, main_window, parent)
        self._slider_dragging = False
        self._snapshot_positions = None
        # Positions this dialog last wrote, to recognise edits made elsewhere.
        self._last_written_positions: Any = None

    def _update_molecule_geometry(
        self, positions: Union[np.ndarray, dict[int, np.ndarray]]
    ) -> None:
        """Write *positions* and remember them as this dialog's own edit."""
        super()._update_molecule_geometry(positions)
        try:
            self._last_written_positions = np.asarray(
                self.mol.GetConformer().GetPositions(), dtype=float
            ).copy()
        except (AttributeError, RuntimeError, TypeError, ValueError):
            self._last_written_positions = None

    def _drop_stale_positions(self) -> None:
        """Re-capture saved positions when the molecule was moved elsewhere.

        The dialogs build each update from saved positions (the slider
        snapshot, and for angles and dihedrals the selection baseline) so the
        rotation axis stays stable. The dialog is modeless, though: after an
        optimization, an undo or another geometry dialog moves atoms, applying
        from the old copy would silently revert every one of those moves.
        """
        try:
            current = np.asarray(self.mol.GetConformer().GetPositions(), dtype=float)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            return
        if current.ndim != 2:
            return
        baseline = getattr(self, "_baseline_positions", None)
        for known in (self._last_written_positions, self._snapshot_positions, baseline):
            if (
                isinstance(known, np.ndarray)
                and known.shape == current.shape
                and np.allclose(known, current, atol=1e-6)
            ):
                return
        if self._snapshot_positions is not None:
            self._snapshot_positions = current.copy()
        if baseline is not None:
            self._baseline_positions = current.copy()
        self._last_written_positions = None

    def _sync_input_to_slider(
        self,
        val: Union[str, float],
        slider: QSlider,
        scale: float = 1.0,
        wrap: bool = False,
    ) -> None:
        """Sync a numerical value from input text to the slider."""
        try:
            f_val = float(val)
            if wrap:
                f_val = (f_val + 180) % 360 - 180

            slider.blockSignals(True)
            slider.setValue(int(round(f_val * scale)))
            slider.blockSignals(False)
        except (ValueError, TypeError):
            # Safe defensive fallback catching ValueError, TypeError
            logging.debug("Suppressed non-critical error", exc_info=True)

    def on_slider_pressed(self) -> None:
        """Prepare for a slider drag operation by saving a geometry snapshot."""
        if not self._is_selection_complete():
            return

        self._slider_dragging = True

        # Capture geometry snapshot to ensure stable axes during rotation/dragging
        self._snapshot_positions = self.mol.GetConformer().GetPositions().copy()

    def on_slider_released(self) -> None:
        """Finalize a slider drag operation."""
        self._slider_dragging = False
        # Snapshot is usually kept until selection changes to preserve turn direction,
        # but subclasses can override this behavior if needed.
        self.main_window.view_3d_manager.draw_molecule_3d(self.mol)
        if hasattr(self.main_window.view_3d_manager, "update_chiral_labels"):
            self.main_window.view_3d_manager.update_chiral_labels()

    def on_slider_value_changed_click(
        self, value: int, input_box: QLineEdit, scale: float = 1.0
    ) -> None:
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

    def on_slider_moved_realtime(
        self, value: int, input_box: QLineEdit, scale: float = 1.0
    ) -> None:
        """Update geometry in real-time as the slider is dragged."""
        if not self._is_selection_complete():
            return

        input_box.blockSignals(True)
        input_box.setText(f"{value / scale:.3f}")
        input_box.blockSignals(False)

        self.apply_geometry_update(float(value / scale))

    def _is_selection_complete(self) -> bool:
        """Must be implemented by subclass to check if enough atoms are picked."""
        raise NotImplementedError

    def apply_geometry_update(self, new_value: float) -> None:  # pylint: disable=arguments-renamed
        """Must be implemented by subclass to perform the actual RdKit/Numpy update."""
        raise NotImplementedError
