#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import Qt, QEvent
import numpy as np
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
        if (obj == self.main_window.plotter.interactor and 
            event.type() == QEvent.Type.MouseButtonPress and 
            event.button() == Qt.MouseButton.LeftButton):
            
            try:
                # VTKイベント座標を取得（元のロジックと同じ）
                interactor = self.main_window.plotter.interactor
                click_pos = interactor.GetEventPosition()
                picker = self.main_window.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer)

                if picker.GetActor() is self.main_window.atom_actor:
                    picked_position = np.array(picker.GetPickPosition())
                    distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                    closest_atom_idx = np.argmin(distances)

                    # 範囲チェックを追加
                    if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                        # クリック閾値チェック（元のロジックと同じ）
                        atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                        if atom:
                            atomic_num = atom.GetAtomicNum()
                            vdw_radius = pt.GetRvdw(atomic_num)
                            click_threshold = vdw_radius * 1.5

                            if distances[closest_atom_idx] < click_threshold:
                                # We handled the pick (atom clicked) -> consume the event so
                                # other UI elements (including the VTK interactor observers)
                                # don't also process it. Set a flag on the main window so
                                # the VTK-based handlers can ignore the same logical click
                                # when it arrives via the VTK event pipeline.
                                try:
                                    self.main_window._picking_consumed = True
                                except Exception:
                                    pass
                                self.on_atom_picked(int(closest_atom_idx))
                                return True
                
                # 原子以外をクリックした場合は選択をクリア（Measurementモードと同じロジック）
                if hasattr(self, 'clear_selection'):
                    self.clear_selection()
                # We did not actually pick an atom; do NOT consume the event here so
                # the interactor and CustomInteractorStyle can handle camera rotation
                # and other behaviors. Returning False (or calling the base
                # implementation) allows normal processing to continue.
                return False
                    
            except Exception as e:
                print(f"Error in eventFilter: {e}")
                # On exception, don't swallow the event either — let the normal
                # event pipeline continue so the UI remains responsive.
                return False

        return super().eventFilter(obj, event)
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
        # Ensure the main window flag exists
        try:
            self.main_window._picking_consumed = False
        except Exception:
            pass
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, 'picking_enabled') and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        try:
            # Clear any leftover flag when picking is disabled
            if hasattr(self.main_window, '_picking_consumed'):
                self.main_window._picking_consumed = False
        except Exception:
            pass
    
    def try_alternative_picking(self, x, y):
        """代替のピッキング方法（使用しない）"""
