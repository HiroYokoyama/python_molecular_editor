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
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton, QMessageBox
)
from PyQt6.QtCore import Qt
import numpy as np

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

class AlignPlaneDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, plane, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.plane = plane
        self.selected_atoms = set()
        
        # 事前選択された原子を追加
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms)
        
        self.init_ui()
        
        # 事前選択された原子にラベルを追加
        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()
    
    def init_ui(self):
        plane_names = {'xy': 'XY', 'xz': 'XZ', 'yz': 'YZ'}
        self.setWindowTitle(f"Align to {plane_names[self.plane]} Plane")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(f"Click atoms in the 3D view to select them for align to the {plane_names[self.plane]} plane. At least 3 atoms are required.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for alignment")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply align")
        self.apply_button.clicked.connect(self.apply_PlaneAlign)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def enable_picking(self):
        """3Dビューでの原子選択を有効にする"""
        self.main_window.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        if hasattr(self, 'picking_enabled') and self.picking_enabled:
            self.main_window.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_PlaneAlign()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
    
    def select_all_atoms(self):
        """Select all atoms in the current molecule and update labels/UI."""
        try:
            # Prefer RDKit molecule if available
            if hasattr(self, 'mol') and self.mol is not None:
                try:
                    n = self.mol.GetNumAtoms()
                    # create a set of indices [0..n-1]
                    self.selected_atoms = set(range(n))
                except Exception:
                    # fallback to main_window data map
                    self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()
            else:
                # fallback to main_window data map
                self.selected_atoms = set(self.main_window.data.atoms.keys()) if hasattr(self.main_window, 'data') else set()

            # Update labels and display
            self.show_atom_labels()
            self.update_display()

        except Exception as e:
            QMessageBox.warning(self, "Warning", f"Failed to select all atoms: {e}")
    
    def update_display(self):
        """表示を更新"""
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select for align (minimum 3 required)")
            self.apply_button.setEnabled(False)
        else:
            atom_list = sorted(self.selected_atoms)
            atom_display = []
            for i, atom_idx in enumerate(atom_list):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_display.append(f"#{i+1}: {symbol}({atom_idx})")
            
            self.selection_label.setText(f"Selected {count} atoms: {', '.join(atom_display)}")
            self.apply_button.setEnabled(count >= 3)
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
            
        if self.selected_atoms:
            sorted_atoms = sorted(self.selected_atoms)
            
            for i, atom_idx in enumerate(sorted_atoms):
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"#{i+1}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color='blue',
                    always_visible=True
                )
                self.selection_labels.append(label_actor)
    
    def clear_atom_labels(self):
        """原子ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except Exception:
                    pass
            self.selection_labels = []
    
    def apply_PlaneAlign(self):
        """alignを適用（回転ベース）"""
        if len(self.selected_atoms) < 3:
            QMessageBox.warning(self, "Warning", "Please select at least 3 atoms for align.")
            return
        try:

            # 選択された原子の位置を取得
            selected_indices = list(self.selected_atoms)
            selected_positions = self.main_window.atom_positions_3d[selected_indices].copy()

            # 重心を計算
            centroid = np.mean(selected_positions, axis=0)

            # 重心を原点に移動
            centered_positions = selected_positions - centroid

            # 主成分分析で最適な平面を見つける
            # 選択された原子の座標の共分散行列を計算
            cov_matrix = np.cov(centered_positions.T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

            # 固有値が最も小さい固有ベクトルが平面の法線方向
            normal_vector = eigenvectors[:, 0]  # 最小固有値に対応する固有ベクトル

            # 目標の平面の法線ベクトルを定義
            if self.plane == 'xy':
                target_normal = np.array([0, 0, 1])  # Z軸方向
            elif self.plane == 'xz':
                target_normal = np.array([0, 1, 0])  # Y軸方向
            elif self.plane == 'yz':
                target_normal = np.array([1, 0, 0])  # X軸方向

            # 法線ベクトルの向きを調整（内積が正になるように）
            if np.dot(normal_vector, target_normal) < 0:
                normal_vector = -normal_vector

            # 回転軸と回転角度を計算
            rotation_axis = np.cross(normal_vector, target_normal)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-10:  # 回転が必要な場合
                rotation_axis = rotation_axis / rotation_axis_norm
                cos_angle = np.dot(normal_vector, target_normal)
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                rotation_angle = np.arccos(cos_angle)

                # Rodrigues回転公式を使用して全分子を回転
                def rodrigues_rotation(v, axis, angle):
                    cos_a = np.cos(angle)
                    sin_a = np.sin(angle)
                    return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)

                # 分子全体を回転させる
                conf = self.mol.GetConformer()
                for i in range(self.mol.GetNumAtoms()):
                    current_pos = np.array(conf.GetAtomPosition(i))
                    # 重心基準で回転
                    centered_pos = current_pos - centroid
                    rotated_pos = rodrigues_rotation(centered_pos, rotation_axis, rotation_angle)
                    new_pos = rotated_pos + centroid
                    conf.SetAtomPosition(i, new_pos.tolist())
                    self.main_window.atom_positions_3d[i] = new_pos

            # 3D表示を更新
            self.main_window.draw_molecule_3d(self.mol)

            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Undo状態を保存
            self.main_window.push_undo_state()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply align: {str(e)}")
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        super().accept()
