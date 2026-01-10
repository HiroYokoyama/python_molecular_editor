#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
main_window_edit_3d.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowEdit3d
"""


import numpy as np


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    from . import sip_isdeleted_safe
except Exception:
    from modules import sip_isdeleted_safe

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QGraphicsTextItem
)

from PyQt6.QtGui import (
    QColor, QFont
)


from PyQt6.QtCore import (
    QPointF
)

import pyvista as pv

# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None
    
# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import VDW_RADII
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import VDW_RADII


# --- クラス定義 ---
class MainWindowEdit3d(object):
    """ main_window.py から分離された機能クラス """


    def toggle_measurement_mode(self, checked):
        """測定モードのオン/オフを切り替える"""
        if checked:
            # 測定モードをオンにする時は、3D Dragモードを無効化
            if self.is_3d_edit_mode:
                self.edit_3d_action.setChecked(False)
                self.toggle_3d_edit_mode(False)
            
            # アクティブな3D編集ダイアログを閉じる
            self.close_all_3d_edit_dialogs()
        
        self.measurement_mode = checked
        
        if not checked:
            self.clear_measurement_selection()
        
        # ボタンのテキストとステータスメッセージを更新
        if checked:
            self.statusBar().showMessage("Measurement mode enabled. Click atoms to measure distances/angles/dihedrals.")
        else:
            self.statusBar().showMessage("Measurement mode disabled.")
    


    def close_all_3d_edit_dialogs(self):
        """すべてのアクティブな3D編集ダイアログを閉じる"""
        dialogs_to_close = self.active_3d_dialogs.copy()
        for dialog in dialogs_to_close:
            try:
                dialog.close()
            except Exception:
                pass
        self.active_3d_dialogs.clear()



    def handle_measurement_atom_selection(self, atom_idx):
        """測定用の原子選択を処理する"""
        # 既に選択されている原子の場合は除外
        if atom_idx in self.selected_atoms_for_measurement:
            return
        
        self.selected_atoms_for_measurement.append(atom_idx)
        
        '''
        # 4つ以上選択された場合はクリア
        if len(self.selected_atoms_for_measurement) > 4:
            self.clear_measurement_selection()
            self.selected_atoms_for_measurement.append(atom_idx)
        '''
        
        # 原子にラベルを追加
        self.add_measurement_label(atom_idx, len(self.selected_atoms_for_measurement))
        
        # 測定値を計算して表示
        self.calculate_and_display_measurements()



    def add_measurement_label(self, atom_idx, label_number):
        """原子に数字ラベルを追加する"""
        if not self.current_mol or atom_idx >= self.current_mol.GetNumAtoms():
            return
        
        # 測定ラベルリストを更新
        self.measurement_labels.append((atom_idx, str(label_number)))
        
        # 3Dビューの測定ラベルを再描画
        self.update_measurement_labels_display()
        
        # 2Dビューの測定ラベルも更新
        self.update_2d_measurement_labels()



    def update_measurement_labels_display(self):
        """測定ラベルを3D表示に描画する（原子中心配置）"""
        try:
            # 既存の測定ラベルを削除
            self.plotter.remove_actor('measurement_labels')
        except Exception:
            pass
        
        if not self.measurement_labels or not self.current_mol:
            return
        
        # ラベル位置とテキストを準備
        pts, labels = [], []
        for atom_idx, label_text in self.measurement_labels:
            if atom_idx < len(self.atom_positions_3d):
                coord = self.atom_positions_3d[atom_idx].copy()
                # オフセットを削除して原子中心に配置
                pts.append(coord)
                labels.append(label_text)
        
        if pts and labels:
            # PyVistaのpoint_labelsを使用（赤色固定）
            self.plotter.add_point_labels(
                np.array(pts), 
                labels, 
                font_size=16,
                point_size=0,
                text_color='red',  # 測定時は常に赤色
                name='measurement_labels',
                always_visible=True,
                tolerance=0.01,
                show_points=False
            )



    def clear_measurement_selection(self):
        """測定選択をクリアする"""
        self.selected_atoms_for_measurement.clear()
        
        # 3Dビューのラベルを削除
        self.measurement_labels.clear()
        try:
            self.plotter.remove_actor('measurement_labels')
        except Exception:
            pass
        
        # 2Dビューの測定ラベルも削除
        self.clear_2d_measurement_labels()
        
        # 測定結果のテキストを削除
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
                self.measurement_text_actor = None
            except Exception:
                pass
        
        self.plotter.render()



    def update_2d_measurement_labels(self):
        """2Dビューで測定ラベルを更新表示する"""
        # 既存の2D測定ラベルを削除
        self.clear_2d_measurement_labels()
        
        # 現在の分子から原子-AtomItemマッピングを作成
        if not self.current_mol or not hasattr(self, 'data') or not self.data.atoms:
            return
            
        # RDKit原子インデックスから2D AtomItemへのマッピングを作成
        atom_idx_to_item = {}
        
        # シーンからAtomItemを取得してマッピング
        if hasattr(self, 'scene'):
            for item in self.scene.items():
                if hasattr(item, 'atom_id') and hasattr(item, 'symbol'):  # AtomItemかチェック
                    # 原子IDから対応するRDKit原子インデックスを見つける
                    rdkit_idx = self.find_rdkit_atom_index(item)
                    if rdkit_idx is not None:
                        atom_idx_to_item[rdkit_idx] = item
        
        # 測定ラベルを2Dビューに追加
        if not hasattr(self, 'measurement_label_items_2d'):
            self.measurement_label_items_2d = []
            
        for atom_idx, label_text in self.measurement_labels:
            if atom_idx in atom_idx_to_item:
                atom_item = atom_idx_to_item[atom_idx]
                self.add_2d_measurement_label(atom_item, label_text)



    def add_2d_measurement_label(self, atom_item, label_text):
        """特定のAtomItemに測定ラベルを追加する"""
        # ラベルアイテムを作成
        label_item = QGraphicsTextItem(label_text)
        label_item.setDefaultTextColor(QColor(255, 0, 0))  # 赤色
        label_item.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        
        # Z値を設定して最前面に表示（原子ラベルより上）
        label_item.setZValue(2000)  # より高い値で確実に最前面に配置
        
        # 原子の右上により近く配置
        atom_pos = atom_item.pos()
        atom_rect = atom_item.boundingRect()
        label_pos = QPointF(
            atom_pos.x() + atom_rect.width() / 4 + 2,
            atom_pos.y() - atom_rect.height() / 4 - 8
        )
        label_item.setPos(label_pos)
        
        # シーンに追加
        self.scene.addItem(label_item)
        self.measurement_label_items_2d.append(label_item)



    def clear_2d_measurement_labels(self):
        """2Dビューの測定ラベルを全て削除する"""
        if hasattr(self, 'measurement_label_items_2d'):
            for label_item in self.measurement_label_items_2d:
                try:
                    # Avoid touching partially-deleted wrappers
                    if sip_isdeleted_safe(label_item):
                        continue
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except Exception:
                        # Scene access or removal failed; skip
                        continue
                except Exception:
                    # If sip check itself fails, fall back to best-effort removal
                    try:
                        if label_item.scene():
                            self.scene.removeItem(label_item)
                    except Exception:
                        continue
            self.measurement_label_items_2d.clear()



    def find_rdkit_atom_index(self, atom_item):
        """AtomItemから対応するRDKit原子インデックスを見つける"""
        if not self.current_mol or not atom_item:
            return None
        
        # マッピング辞書を使用（最も確実）
        if hasattr(self, 'atom_id_to_rdkit_idx_map') and atom_item.atom_id in self.atom_id_to_rdkit_idx_map:
            return self.atom_id_to_rdkit_idx_map[atom_item.atom_id]
        
        # マッピングが存在しない場合はNone（外部ファイル読み込み時など）
        return None



    def calculate_and_display_measurements(self):
        """選択された原子に基づいて測定値を計算し表示する"""
        num_selected = len(self.selected_atoms_for_measurement)
        if num_selected < 2:
            return
        
        measurement_text = []
        
        if num_selected >= 2:
            # 距離の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            distance = self.calculate_distance(atom1_idx, atom2_idx)
            measurement_text.append(f"Distance 1-2: {distance:.3f} Å")
        
        if num_selected >= 3:
            # 角度の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1] 
            atom3_idx = self.selected_atoms_for_measurement[2]
            angle = self.calculate_angle(atom1_idx, atom2_idx, atom3_idx)
            measurement_text.append(f"Angle 1-2-3: {angle:.2f}°")
        
        if num_selected >= 4:
            # 二面角の計算
            atom1_idx = self.selected_atoms_for_measurement[0]
            atom2_idx = self.selected_atoms_for_measurement[1]
            atom3_idx = self.selected_atoms_for_measurement[2]
            atom4_idx = self.selected_atoms_for_measurement[3]
            dihedral = self.calculate_dihedral(atom1_idx, atom2_idx, atom3_idx, atom4_idx)
            measurement_text.append(f"Dihedral 1-2-3-4: {dihedral:.2f}°")
        
        # 測定結果を3D画面の右上に表示
        self.display_measurement_text(measurement_text)



    def calculate_distance(self, atom1_idx, atom2_idx):
        """2原子間の距離を計算する"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])
        return np.linalg.norm(pos2 - pos1)



    def calculate_angle(self, atom1_idx, atom2_idx, atom3_idx):
        """3原子の角度を計算する（中央が頂点）"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])  # 頂点
        pos3 = np.array(self.atom_positions_3d[atom3_idx])
        
        # ベクトルを計算
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        # 角度を計算（ラジアンから度に変換）
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        # 数値誤差による範囲外の値をクリップ
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)



    def calculate_dihedral(self, atom1_idx, atom2_idx, atom3_idx, atom4_idx):
        """4原子の二面角を計算する（正しい公式を使用）"""
        pos1 = np.array(self.atom_positions_3d[atom1_idx])
        pos2 = np.array(self.atom_positions_3d[atom2_idx])
        pos3 = np.array(self.atom_positions_3d[atom3_idx])
        pos4 = np.array(self.atom_positions_3d[atom4_idx])
        
        # Vectors between consecutive atoms
        v1 = pos2 - pos1  # 1->2
        v2 = pos3 - pos2  # 2->3 (central bond)
        v3 = pos4 - pos3  # 3->4
        
        # Normalize the central bond vector
        v2_norm = v2 / np.linalg.norm(v2)
        
        # Calculate plane normal vectors
        n1 = np.cross(v1, v2)  # Normal to plane 1-2-3
        n2 = np.cross(v2, v3)  # Normal to plane 2-3-4
        
        # Normalize the normal vectors
        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)
        
        if n1_norm == 0 or n2_norm == 0:
            return 0.0  # Atoms are collinear
        
        n1 = n1 / n1_norm
        n2 = n2 / n2_norm
        
        # Calculate the cosine of the dihedral angle
        cos_angle = np.dot(n1, n2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        # Calculate the sine for proper sign determination
        sin_angle = np.dot(np.cross(n1, n2), v2_norm)
        
        # Calculate the dihedral angle with correct sign
        angle_rad = np.arctan2(sin_angle, cos_angle)
        return np.degrees(angle_rad)



    def display_measurement_text(self, measurement_lines):
        """測定結果のテキストを3D画面の左上に表示する（小さな等幅フォント）"""
        # 既存のテキストを削除
        if self.measurement_text_actor:
            try:
                self.plotter.remove_actor(self.measurement_text_actor)
            except Exception:
                pass
        
        if not measurement_lines:
            self.measurement_text_actor = None
            return
        
        # テキストを結合
        text = '\n'.join(measurement_lines)
        
        # 背景色から適切なテキスト色を決定
        try:
            bg_color_hex = self.settings.get('background_color', '#919191')
            bg_qcolor = QColor(bg_color_hex)
            if bg_qcolor.isValid():
                luminance = bg_qcolor.toHsl().lightness()
                text_color = 'black' if luminance > 128 else 'white'
            else:
                text_color = 'white'
        except Exception:
            text_color = 'white'
        
        # 左上に表示（小さな等幅フォント）
        self.measurement_text_actor = self.plotter.add_text(
            text,
            position='upper_left',
            font_size=10,  # より小さく
            color=text_color,  # 背景に合わせた色
            font='courier',  # 等幅フォント
            name='measurement_display'
        )
        
        self.plotter.render()

    # --- 3D Drag functionality ---
    


    def toggle_atom_selection_3d(self, atom_idx):
        """3Dビューで原子の選択状態をトグルする"""
        if atom_idx in self.selected_atoms_3d:
            self.selected_atoms_3d.remove(atom_idx)
        else:
            self.selected_atoms_3d.add(atom_idx)
        
        # 選択状態のビジュアルフィードバックを更新
        self.update_3d_selection_display()
    


    def clear_3d_selection(self):
        """3Dビューでの原子選択をクリア"""
        self.selected_atoms_3d.clear()
        self.update_3d_selection_display()
    


    def update_3d_selection_display(self):
        """3Dビューでの選択原子のハイライト表示を更新"""
        try:
            # 既存の選択ハイライトを削除
            self.plotter.remove_actor('selection_highlight')
        except Exception:
            pass
        
        if not self.selected_atoms_3d or not self.current_mol:
            self.plotter.render()
            return
        
        # 選択された原子のインデックスリストを作成
        selected_indices = list(self.selected_atoms_3d)
        
        # 選択された原子の位置を取得
        selected_positions = self.atom_positions_3d[selected_indices]
        
        # 原子の半径を少し大きくしてハイライト表示
        selected_radii = np.array([VDW_RADII.get(
            self.current_mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3 
            for i in selected_indices])
        
        # ハイライト用のデータセットを作成
        highlight_source = pv.PolyData(selected_positions)
        highlight_source['radii'] = selected_radii
        
        # 黄色の半透明球でハイライト
        highlight_glyphs = highlight_source.glyph(
            scale='radii', 
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16), 
            orient=False
        )
        
        self.plotter.add_mesh(
            highlight_glyphs, 
            color='yellow', 
            opacity=0.3, 
            name='selection_highlight'
        )
        
        self.plotter.render()
    
    def remove_dialog_from_list(self, dialog):
        """ダイアログをアクティブリストから削除"""
        if dialog in self.active_3d_dialogs:
            self.active_3d_dialogs.remove(dialog)
