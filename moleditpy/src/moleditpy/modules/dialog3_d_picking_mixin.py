#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import numpy as np
from PyQt6.QtCore import QEvent, Qt

try:
    from .constants import pt
except Exception:
    from modules.constants import pt


class Dialog3DPickingMixin:
    """3D原子選択のための共通機能を提供するMixin"""

    def __init__(self):
        """Mixinの初期化"""
        self.picking_enabled = False

    def eventFilter(self, obj, event):
        """3Dビューでのマウスクリックをキャプチャする（元の3D editロジックを正確に再現）"""
        if (
            obj == self.main_window.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonPress
            and event.button() == Qt.MouseButton.LeftButton
        ):
            # Start tracking for smart selection (click vs drag)
            self._mouse_press_pos = event.pos()
            self._mouse_moved = False

            try:
                # VTKイベント座標を取得（元のロジックと同じ）
                interactor = self.main_window.plotter.interactor
                click_pos = interactor.GetEventPosition()
                picker = self.main_window.plotter.picker
                picker.Pick(
                    click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer
                )

                if picker.GetActor() is self.main_window.atom_actor:
                    picked_position = np.array(picker.GetPickPosition())
                    distances = np.linalg.norm(
                        self.main_window.atom_positions_3d - picked_position, axis=1
                    )
                    closest_atom_idx = np.argmin(distances)

                    # 範囲チェックを追加
                    if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                        # クリック閾値チェック（元のロジックと同じ）
                        atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                        if atom:
                            try:
                                atomic_num = atom.GetAtomicNum()
                                vdw_radius = pt.GetRvdw(atomic_num)
                                if vdw_radius < 0.1:
                                    vdw_radius = 1.5
                            except Exception:
                                vdw_radius = 1.5
                            click_threshold = vdw_radius * 1.5

                            if distances[closest_atom_idx] < click_threshold:
                                try:
                                    self.main_window._picking_consumed = True
                                except Exception:
                                    import traceback

                                    traceback.print_exc()
                                self.on_atom_picked(int(closest_atom_idx))

                                # We picked an atom, so stop tracking for background click
                                self._mouse_press_pos = None
                                return True

                # 原子以外をクリックした場合
                # 即時には解除せず、回転操作（ドラッグ）を許可する。
                # 実際の解除は MouseButtonRelease イベントで行う。
                return False

            except Exception as e:
                print(f"Error in eventFilter: {e}")
                # On exception, don't swallow the event either — let the normal
                # event pipeline continue so the UI remains responsive.
                return False

        # Add movement tracking for smart selection
        elif (
            obj == self.main_window.plotter.interactor
            and event.type() == QEvent.Type.MouseMove
        ):
            if hasattr(self, "_mouse_press_pos") and self._mouse_press_pos is not None:
                # Check if moved significantly
                diff = event.pos() - self._mouse_press_pos
                if diff.manhattanLength() > 3:
                    self._mouse_moved = True

        # Add release handling for smart selection
        elif (
            obj == self.main_window.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonRelease
            and event.button() == Qt.MouseButton.LeftButton
        ):
            if hasattr(self, "_mouse_press_pos") and self._mouse_press_pos is not None:
                if not getattr(self, "_mouse_moved", False):
                    # Pure click (no drag) on background -> Clear selection
                    if hasattr(self, "clear_selection"):
                        self.clear_selection()

                # Reset state
                self._mouse_press_pos = None
                self._mouse_moved = False

        return super().eventFilter(obj, event)

    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
        # Ensure the main window flag exists
        try:
            self.main_window._picking_consumed = False
        except Exception:
            import traceback

            traceback.print_exc()

    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, "picking_enabled") and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        try:
            # Clear any leftover flag when picking is disabled
            if hasattr(self.main_window, "_picking_consumed"):
                self.main_window._picking_consumed = False
        except Exception:
            import traceback

            traceback.print_exc()

    def try_alternative_picking(self, x, y):
        """代替のピッキング方法（使用しない）"""

    # ------------------------------------------------------------------
    # Label management (shared across dialogs)
    # ------------------------------------------------------------------

    def clear_atom_labels(self):
        """Remove all label actors from the plotter."""
        if hasattr(self, "selection_labels"):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except Exception:
                    import traceback

                    traceback.print_exc()
            self.selection_labels = []

    # Alias — some dialogs use this name instead.
    clear_selection_labels = clear_atom_labels

    def add_selection_label(self, atom_idx, label_text, color="yellow"):
        """Add a point label at the position of *atom_idx*.

        Parameters
        ----------
        atom_idx : int
            Index into ``self.main_window.atom_positions_3d``.
        label_text : str
            Text shown next to the atom.
        color : str, optional
            Label colour (default ``'yellow'``).
        """
        if not hasattr(self, "selection_labels"):
            self.selection_labels = []

        pos = self.main_window.atom_positions_3d[atom_idx]

        label_actor = self.main_window.plotter.add_point_labels(
            [pos],
            [label_text],
            point_size=20,
            font_size=12,
            text_color=color,
            always_visible=True,
        )
        self.selection_labels.append(label_actor)

    def show_atom_labels_for(self, atoms_and_labels, color="yellow"):
        """Clear existing labels and add new ones for each *(idx, text)* pair.

        Parameters
        ----------
        atoms_and_labels : list[tuple[int, str]]
            Each element is ``(atom_idx, label_text)``.
        color : str, optional
            Label colour (default ``'yellow'``).
        """
        self.clear_atom_labels()

        if not hasattr(self, "selection_labels"):
            self.selection_labels = []

        for atom_idx, label_text in atoms_and_labels:
            pos = self.main_window.atom_positions_3d[atom_idx]
            label_actor = self.main_window.plotter.add_point_labels(
                [pos],
                [label_text],
                point_size=20,
                font_size=12,
                text_color=color,
                always_visible=True,
            )
            self.selection_labels.append(label_actor)
