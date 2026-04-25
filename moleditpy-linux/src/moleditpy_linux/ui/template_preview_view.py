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

from typing import Any, Optional

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QResizeEvent, QShowEvent
from PyQt6.QtWidgets import QGraphicsScene, QGraphicsView


class TemplatePreviewView(QGraphicsView):
    """Custom view class for template previews."""

    def __init__(self, scene: QGraphicsScene) -> None:
        super().__init__(scene)
        self.original_scene_rect = None
        self.template_data = None  # Store template data for dynamic redrawing
        self.parent_dialog = None  # Reference to parent dialog for redraw access

    def set_template_data(self, template_data: Any, parent_dialog: Any) -> None:
        """Set the template data and parent dialog reference for dynamic redrawing."""
        self.template_data = template_data
        self.parent_dialog = parent_dialog

    def resizeEvent(self, event: Optional[QResizeEvent]) -> None:
        """Handle resize events to ensure the preview is properly fitted."""
        super().resizeEvent(event)
        if self.original_scene_rect and not self.original_scene_rect.isEmpty():
            # Delay the fitInView call to ensure proper widget sizing
            QTimer.singleShot(10, self.refit_view)

    def refit_view(self) -> None:
        """Refit the view to the stored original scene rectangle."""
        try:
            if self.original_scene_rect and not self.original_scene_rect.isEmpty():
                self.fitInView(
                    self.original_scene_rect, Qt.AspectRatioMode.KeepAspectRatio
                )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            print(f"Warning: Failed to refit template preview: {e}")

    def showEvent(self, event: Optional[QShowEvent]) -> None:
        """Handle the show event to ensure the preview fits correctly upon display."""
        super().showEvent(event)
        # Ensure proper fitting when widget becomes visible
        if self.original_scene_rect:
            QTimer.singleShot(50, self.refit_view)

    def redraw_with_current_size(self) -> None:
        """Redraw the template structure scaled to the current view dimensions."""
        if not (self.template_data and self.parent_dialog):
            return
        try:
            self.scene().clear()
            view_size = (self.width(), self.height())
            self.parent_dialog.draw_template_preview(
                self.scene(), self.template_data, view_size
            )
            bounding_rect = self.scene().itemsBoundingRect()
            if (
                not bounding_rect.isEmpty()
                and bounding_rect.width() > 0
                and bounding_rect.height() > 0
            ):
                content_size = max(bounding_rect.width(), bounding_rect.height())
                padding = max(20, content_size * 0.2)
                padded_rect = bounding_rect.adjusted(
                    -padding, -padding, padding, padding
                )
                self.scene().setSceneRect(padded_rect)
                self.original_scene_rect = padded_rect
                QTimer.singleShot(
                    10,
                    lambda: self.fitInView(
                        padded_rect, Qt.AspectRatioMode.KeepAspectRatio
                    ),
                )
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            print(f"Warning: Failed to redraw template preview: {e}")
