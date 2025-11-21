#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import (
    QGraphicsView
)


from PyQt6.QtCore import (
    Qt, QPointF, QEvent
)

class ZoomableView(QGraphicsView):
    """ マウスホイールでのズームと、中ボタン or Shift+左ドラッグでのパン機能を追加したQGraphicsView """
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
        """ マウスホイールを回した際のイベント """
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
        """ 中ボタン or Shift+左ボタンが押されたらパン（視点移動）モードを開始 """
        is_middle_button = event.button() == Qt.MouseButton.MiddleButton
        is_shift_left_button = (event.button() == Qt.MouseButton.LeftButton and
                                event.modifiers() & Qt.KeyboardModifier.ShiftModifier)

        if is_middle_button or is_shift_left_button:
            self._is_panning = True
            self._pan_start_pos = event.pos() # ビューポート座標で開始点を記録
            # 現在のスクロールバーの位置を記録
            self._pan_start_scroll_h = self.horizontalScrollBar().value()
            self._pan_start_scroll_v = self.verticalScrollBar().value()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        """ パンモード中にマウスを動かしたら、スクロールバーを操作して視点を移動させる """
        if self._is_panning:
            delta = event.pos() - self._pan_start_pos # マウスの移動量を計算
            # 開始時のスクロール位置から移動量を引いた値を新しいスクロール位置に設定
            self.horizontalScrollBar().setValue(self._pan_start_scroll_h - delta.x())
            self.verticalScrollBar().setValue(self._pan_start_scroll_v - delta.y())
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        """ パンに使用したボタンが離されたらパンモードを終了 """
        # パンを開始したボタン（中 or 左）のどちらかが離されたかをチェック
        is_middle_button_release = event.button() == Qt.MouseButton.MiddleButton
        is_left_button_release = event.button() == Qt.MouseButton.LeftButton

        if self._is_panning and (is_middle_button_release or is_left_button_release):
            self._is_panning = False
            # 現在の描画モードに応じたカーソルに戻す
            current_mode = self.scene().mode if self.scene() else 'select'
            if current_mode == 'select':
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif current_mode.startswith(('atom', 'bond', 'template')):
                self.setCursor(Qt.CursorShape.CrossCursor)
            elif current_mode.startswith(('charge', 'radical')):
                self.setCursor(Qt.CursorShape.CrossCursor)
            else:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def viewportEvent(self, event):
        """ Macのトラックパッドなどのネイティブジェスチャー（ピンチ）に対応 """
        if event.type() == QEvent.Type.NativeGesture:
            # ピンチズーム（ZoomNativeGesture）の場合
            if event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
                # event.value() は拡大縮小の倍率差分を返します
                # 例: 拡大時は正の値、縮小時は負の値
                factor = 1.0 + event.value()
                
                # 既存の最大・最小スケール制限を考慮する場合（オプション）
                current_scale = self.transform().m11()
                min_scale, max_scale = 0.05, 20.0
                
                # 制限内であればスケールを適用
                # (厳密に制限を守るならif文でガードしますが、最小実装としては以下で動作します)
                self.scale(factor, factor)
                return True
                
        return super().viewportEvent(event)
