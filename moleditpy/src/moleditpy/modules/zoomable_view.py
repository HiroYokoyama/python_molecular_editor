#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import QEvent, QPointF, Qt
from PyQt6.QtWidgets import QGraphicsView


class ZoomableView(QGraphicsView):
    """QGraphicsView with zoom functionality via mouse wheel and panning via middle button or Shift+left drag"""

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self.main_window = parent
        self.setAcceptDrops(False)

        self._is_panning = False
        self._pan_start_pos = QPointF()
        self._pan_start_scroll_h = 0
        self._pan_start_scroll_v = 0

    def wheelEvent(self, event):
        """Event handler for mouse wheel rotation"""
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            zoom_in_factor = 1.1
            zoom_out_factor = 1 / zoom_in_factor

            transform = self.transform()
            current_scale = transform.m11()
            min_scale, max_scale = 0.05, 20.0

            if event.angleDelta().y() > 0:
                if max_scale > current_scale:
                    self.scale(zoom_in_factor, zoom_in_factor)
            else:
                if min_scale < current_scale:
                    self.scale(zoom_out_factor, zoom_out_factor)

            event.accept()
        else:
            super().wheelEvent(event)

    def mousePressEvent(self, event):
        """Start panning (view movement) mode if middle button or Shift+left button is pressed"""
        is_middle_button = event.button() == Qt.MouseButton.MiddleButton
        is_shift_left_button = (
            event.button() == Qt.MouseButton.LeftButton
            and event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        )

        if is_middle_button or is_shift_left_button:
            self._is_panning = True
            self._pan_start_pos = event.pos()  # Record start point in viewport coordinates
            # Record current scrollbar values
            self._pan_start_scroll_h = self.horizontalScrollBar().value()
            self._pan_start_scroll_v = self.verticalScrollBar().value()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        """Handle view movement by updating scrollbars during panning"""
        if self._is_panning:
            delta = event.pos() - self._pan_start_pos  # Calculate mouse movement delta
            # Update scroll position by subtracting movement delta from start position
            self.horizontalScrollBar().setValue(self._pan_start_scroll_h - delta.x())
            self.verticalScrollBar().setValue(self._pan_start_scroll_v - delta.y())
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        """End panning mode and restore cursor when the relevant button is released"""
        # Check if the middle or left button used for panning was released
        is_middle_button_release = event.button() == Qt.MouseButton.MiddleButton
        is_left_button_release = event.button() == Qt.MouseButton.LeftButton

        if self._is_panning and (is_middle_button_release or is_left_button_release):
            self._is_panning = False
            # Restore cursor based on the current scene mode
            current_mode = self.scene().mode if self.scene() else "select"
            if current_mode == "select":
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif current_mode.startswith(("atom", "bond", "template")):
                self.setCursor(Qt.CursorShape.CrossCursor)
            elif current_mode.startswith(("charge", "radical")):
                self.setCursor(Qt.CursorShape.CrossCursor)
            else:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def viewportEvent(self, event):
        """Handle native gestures (like pinch zoom on trackpads)"""
        if event.type() == QEvent.Type.NativeGesture:
            # Detect pinch zoom gestures
            if event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
                # event.value() returns the scale factor delta
                # (positive for zoom-in, negative for zoom-out)
                factor = 1.0 + event.value()

                current_scale = self.transform().m11()
                min_scale, max_scale = 0.05, 20.0

                # Apply scaling if within limits
                self.scale(factor, factor)
                return True

        return super().viewportEvent(event)
