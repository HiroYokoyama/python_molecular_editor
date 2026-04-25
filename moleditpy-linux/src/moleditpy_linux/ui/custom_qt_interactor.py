#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import time
from typing import Any, Optional

from pyvistaqt import QtInteractor


class CustomQtInteractor(QtInteractor):
    def __init__(
        self,
        parent: Optional[Any] = None,
        main_window: Optional[Any] = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(parent, **kwargs)
        self.main_window = main_window
        self._last_click_time = 0.0
        self._click_count = 0
        self._ignore_next_release = False

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
        all 3D view operations. Also filters out "Ghost Release" (release without
        a corresponding press).
        """
        if self._ignore_next_release:
            self._ignore_next_release = False
            event.accept()
            return

        super().mouseReleaseEvent(event)  # Process parent class event first
        if self.main_window and hasattr(self.main_window.init_manager, "view_2d"):
            self.main_window.init_manager.view_2d.setFocus()

    def mousePressEvent(self, event: Any) -> None:
        """
        Custom mouse press handling to track accumulated clicks and filter out
        triple-clicks.
        """
        current_time = time.time()
        # Reset count if too much time has passed (0.5s is standard double-click time)
        if current_time - self._last_click_time > 0.5:
            self._click_count = 0

        self._click_count += 1
        self._last_click_time = current_time

        # If this is the 3rd click (or more), swallow it to prevent
        # the internal state desync that happens with rapid clicking sequences.
        if self._click_count >= 3:
            self._ignore_next_release = True
            event.accept()
            return

        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event: Any) -> None:
        """Ignore mouse double-clicks on the 3D widget to avoid accidental actions.

        Swallow the double-click event so it doesn't trigger selection, editing,
        or camera jumps. We intentionally do not call the superclass handler.
        Crucially, we also flag the NEXT release event to be swallowed, preventing
        a "Ghost Release" (Release without Press) from reaching VTK.
        """
        current_time = time.time()
        self._last_click_time = current_time
        # Set to 2 to ensure the next click counts as 3rd
        if current_time - self._last_click_time < 0.5:
            self._click_count = 2
        else:
            self._click_count = 2  # Force sync

        self._ignore_next_release = True

        try:
            # Accept the event to mark it handled and prevent further processing.
            event.accept()
        except (AttributeError, RuntimeError, TypeError):
            # If event doesn't support accept for some reason, just return.
            return
