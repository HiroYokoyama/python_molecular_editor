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
from typing import Any, Optional

from PyQt6.QtCore import QEvent
from PyQt6.QtGui import QMouseEvent
from pyvistaqt import QtInteractor


class CustomQtInteractor(QtInteractor):
    """PyVista QtInteractor subclass that exposes the main window for event handling."""

    def __init__(
        self,
        parent: Optional[Any] = None,
        main_window: Optional[Any] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(parent, **kwargs)
        self.main_window = main_window

    def wheelEvent(self, event: Any) -> None:
        """
        Override the mouse wheel event.
        """
        # First call the parent class event to perform normal zoom processing
        super().wheelEvent(event)

        # Force focus back to the 2D view after zoom processing
        if self.main_window and hasattr(self.main_window.init_manager, "view_2d"):
            self.main_window.init_manager.view_2d.setFocus()

    def mouseReleaseEvent(self, event: Any) -> None:
        """
        Override the Qt mouse release event to return focus to the 2D view after
        all 3D view operations.
        """
        super().mouseReleaseEvent(event)  # Process parent class event first
        if self.main_window and hasattr(self.main_window.init_manager, "view_2d"):
            self.main_window.init_manager.view_2d.setFocus()

    def mouseDoubleClickEvent(self, event: Any) -> None:
        """Re-dispatch double-clicks as plain presses so fast clicking works.

        Qt turns every second fast click into a double-click event. Forwarding
        it unchanged would reach VTK with a repeat count and be dropped by the
        interactor style, so the click would be lost (e.g. when rapidly
        selecting atoms in measurement mode). Synthesizing a normal press keeps
        press/release pairing intact and makes each fast click act as a click.
        """
        try:
            synthetic_press = QMouseEvent(
                QEvent.Type.MouseButtonPress,
                event.position(),
                event.globalPosition(),
                event.button(),
                event.buttons(),
                event.modifiers(),
            )
            super().mousePressEvent(synthetic_press)
            event.accept()
        except (AttributeError, RuntimeError, TypeError):
            # Safe defensive fallback catching AttributeError, RuntimeError, TypeError
            logging.debug("Suppressed non-critical error", exc_info=True)
