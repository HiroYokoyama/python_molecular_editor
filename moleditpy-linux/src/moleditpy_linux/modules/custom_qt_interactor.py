#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from pyvistaqt import QtInteractor
import time


class CustomQtInteractor(QtInteractor):
    def __init__(self, parent=None, main_window=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.main_window = main_window
        # Track recent clicks so we can detect and swallow triple-clicks
        # Triple-clicks are not a distinct Qt event on all platforms, so we
        # implement a small timing-based counter here and accept the event
        # when 3 rapid clicks are detected to prevent them from reaching
        # the VTK interactor and causing unexpected behaviour in the 3D view.
        self._last_click_time = 0.0
        self._click_count = 0
        self._ignore_next_release = False

    def wheelEvent(self, event):
        """
        マウスホイールイベントをオーバーライドする。
        """
        # 最初に親クラスのイベントを呼び、通常のズーム処理を実行させる
        super().wheelEvent(event)
        
        # ズーム処理の完了後、2Dビューにフォーカスを強制的に戻す
        if self.main_window and hasattr(self.main_window, 'view_2d'):
            self.main_window.view_2d.setFocus()

    def mouseReleaseEvent(self, event):
        """
        Qtのマウスリリースイベントをオーバーライドし、
        3Dビューでの全ての操作完了後に2Dビューへフォーカスを戻す。
        また、Ghost Release（対応するPressがないRelease）をフィルタリングする。
        """
        if self._ignore_next_release:
            self._ignore_next_release = False
            event.accept()
            return

        super().mouseReleaseEvent(event) # 親クラスのイベントを先に処理
        if self.main_window and hasattr(self.main_window, 'view_2d'):
            self.main_window.view_2d.setFocus()

    def mousePressEvent(self, event):
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

    def mouseDoubleClickEvent(self, event):
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
             self._click_count = 2 # Force sync
        
        # We must ignore the release event that follows this double-click event
        self._ignore_next_release = True
        
        try:
            # Accept the event to mark it handled and prevent further processing.
            event.accept()
        except Exception:
            # If event doesn't support accept for some reason, just return.
            return
