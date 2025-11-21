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
    QDialog, QVBoxLayout, QLabel, QGridLayout, QHBoxLayout, QPushButton, QLineEdit, QMessageBox
)
from PyQt6.QtCore import Qt, QEvent
import numpy as np
import pyvista as pv
from rdkit import Chem
try:
    from .constants import VDW_RADII
except Exception:
    from modules.constants import VDW_RADII

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

class MoveGroupDialog(Dialog3DPickingMixin, QDialog):
    """結合している分子グループを選択して並行移動・回転するダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()
        self.group_atoms = set()  # 選択原子に結合している全原子
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Move Group")
        self.setModal(False)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        self.resize(300,400)  # ウィンドウサイズを設定
        layout = QVBoxLayout(self)
        
        # ドラッグ状態管理
        self.is_dragging_group = False
        self.drag_start_pos = None
        self.mouse_moved_during_drag = False  # ドラッグ中にマウスが動いたかを追跡
        
        # Instructions
        instruction_label = QLabel("Click an atom in the 3D view to select its connected molecule group.\n"
                                   "Left-drag: Move the group\n"
                                   "Right-drag: Rotate the group around its center")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected group display
        self.selection_label = QLabel("No group selected")
        layout.addWidget(self.selection_label)
        
        # Translation controls
        trans_group = QLabel("Translation (Å):")
        trans_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(trans_group)
        
        trans_layout = QGridLayout()
        self.x_trans_input = QLineEdit("0.0")
        self.y_trans_input = QLineEdit("0.0")
        self.z_trans_input = QLineEdit("0.0")
        
        # Enterキーでapply_translationを実行
        self.x_trans_input.returnPressed.connect(self.apply_translation)
        self.y_trans_input.returnPressed.connect(self.apply_translation)
        self.z_trans_input.returnPressed.connect(self.apply_translation)
        
        trans_layout.addWidget(QLabel("X:"), 0, 0)
        trans_layout.addWidget(self.x_trans_input, 0, 1)
        trans_layout.addWidget(QLabel("Y:"), 1, 0)
        trans_layout.addWidget(self.y_trans_input, 1, 1)
        trans_layout.addWidget(QLabel("Z:"), 2, 0)
        trans_layout.addWidget(self.z_trans_input, 2, 1)
        
        trans_button_layout = QHBoxLayout()
        reset_trans_button = QPushButton("Reset")
        reset_trans_button.clicked.connect(self.reset_translation_inputs)
        trans_button_layout.addWidget(reset_trans_button)
        
        apply_trans_button = QPushButton("Apply Translation")
        apply_trans_button.clicked.connect(self.apply_translation)
        trans_button_layout.addWidget(apply_trans_button)
        
        trans_layout.addLayout(trans_button_layout, 3, 0, 1, 2)
        
        layout.addLayout(trans_layout)
        
        layout.addSpacing(10)
        
        # Rotation controls
        rot_group = QLabel("Rotation (degrees):")
        rot_group.setStyleSheet("font-weight: bold;")
        layout.addWidget(rot_group)
        
        rot_layout = QGridLayout()
        self.x_rot_input = QLineEdit("0.0")
        self.y_rot_input = QLineEdit("0.0")
        self.z_rot_input = QLineEdit("0.0")
        
        # Enterキーでapply_rotationを実行
        self.x_rot_input.returnPressed.connect(self.apply_rotation)
        self.y_rot_input.returnPressed.connect(self.apply_rotation)
        self.z_rot_input.returnPressed.connect(self.apply_rotation)
        
        rot_layout.addWidget(QLabel("Around X:"), 0, 0)
        rot_layout.addWidget(self.x_rot_input, 0, 1)
        rot_layout.addWidget(QLabel("Around Y:"), 1, 0)
        rot_layout.addWidget(self.y_rot_input, 1, 1)
        rot_layout.addWidget(QLabel("Around Z:"), 2, 0)
        rot_layout.addWidget(self.z_rot_input, 2, 1)
        
        rot_button_layout = QHBoxLayout()
        reset_rot_button = QPushButton("Reset")
        reset_rot_button.clicked.connect(self.reset_rotation_inputs)
        rot_button_layout.addWidget(reset_rot_button)
        
        apply_rot_button = QPushButton("Apply Rotation")
        apply_rot_button.clicked.connect(self.apply_rotation)
        rot_button_layout.addWidget(apply_rot_button)
        
        rot_layout.addLayout(rot_button_layout, 3, 0, 1, 2)
        
        layout.addLayout(rot_layout)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Enable picking to handle atom selection
        self.enable_picking()
    
    def eventFilter(self, obj, event):
        """3Dビューでのマウスイベント処理 - グループが選択されている場合はCustomInteractorStyleに任せる"""
        if obj == self.main_window.plotter.interactor:
            # ダブルクリック/トリプルクリックで状態が混乱するのを防ぐ
            if event.type() == QEvent.Type.MouseButtonDblClick:
                # ダブルクリックは無視し、状態をリセット
                self.is_dragging_group = False
                self.drag_start_pos = None
                self.mouse_moved_during_drag = False
                self.potential_drag = False
                if hasattr(self, 'clicked_atom_for_toggle'):
                    delattr(self, 'clicked_atom_for_toggle')
                return False
            
            if event.type() == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
                # 前回の状態をクリーンアップ（トリプルクリック対策）
                self.is_dragging_group = False
                self.potential_drag = False
                if hasattr(self, 'clicked_atom_for_toggle'):
                    delattr(self, 'clicked_atom_for_toggle')
                # グループが既に選択されている場合は、CustomInteractorStyleに処理を任せる
                if self.group_atoms:
                    return False
                
                # マウスプレス時の処理
                # マウスプレス時の処理
                try:
                    interactor = self.main_window.plotter.interactor
                    click_pos = interactor.GetEventPosition()
                    
                    # まずピッキングしてどの原子がクリックされたか確認
                    picker = self.main_window.plotter.picker
                    picker.Pick(click_pos[0], click_pos[1], 0, self.main_window.plotter.renderer)
                    
                    clicked_atom_idx = None
                    if picker.GetActor() is self.main_window.atom_actor:
                        picked_position = np.array(picker.GetPickPosition())
                        distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                        closest_atom_idx = np.argmin(distances)
                        
                        # 閾値チェック
                        if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                            atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                            if atom:
                                atomic_num = atom.GetAtomicNum()
                                pt = Chem.GetPeriodicTable()
                                vdw_radius = pt.GetRvdw(atomic_num)
                                click_threshold = vdw_radius * 1.5
                                
                                if distances[closest_atom_idx] < click_threshold:
                                    clicked_atom_idx = int(closest_atom_idx)
                    
                    
                    # クリックされた原子の処理
                    if clicked_atom_idx is not None:
                        if self.group_atoms and clicked_atom_idx in self.group_atoms:
                            # 既存のグループ内の原子 - ドラッグ準備（まだドラッグとは確定しない）
                            self.is_dragging_group = False  # まだドラッグ中ではない
                            self.drag_start_pos = click_pos
                            self.drag_atom_idx = clicked_atom_idx
                            self.mouse_moved_during_drag = False
                            self.potential_drag = True  # ドラッグの可能性がある
                            self.clicked_atom_for_toggle = clicked_atom_idx  # トグル用に保存
                            # イベントを消費せず、カメラ操作を許可（閾値超えたらドラッグ開始）
                            return False
                        else:
                            # グループ外の原子 - 新しいグループを選択
                            # 親クラス（Mixin）のon_atom_pickedを手動で呼ぶ
                            self.on_atom_picked(clicked_atom_idx)
                            return True
                    else:
                        # 原子以外をクリック
                        # グループがあっても通常のカメラ操作を許可
                        return False
                    
                except Exception as e:
                    print(f"Error in mouse press: {e}")
                    return False
            
            elif event.type() == QEvent.Type.MouseMove:
                # マウス移動時の処理
                if getattr(self, 'potential_drag', False) and self.drag_start_pos and not self.is_dragging_group:
                    # potential_drag状態：閾値チェック
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]
                        
                        # 閾値を超えたらドラッグ開始
                        drag_threshold = 5  # ピクセル
                        if abs(dx) > drag_threshold or abs(dy) > drag_threshold:
                            # ドラッグ開始を確定
                            self.is_dragging_group = True
                            self.potential_drag = False
                            try:
                                self.main_window.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                            except Exception:
                                pass
                    except Exception:
                        pass
                    
                    # 閾値以下の場合はカメラ操作を許可
                    if not self.is_dragging_group:
                        return False
                
                if self.is_dragging_group and self.drag_start_pos:
                    # ドラッグモード中 - 移動距離を記録するのみ（リアルタイム更新なし）
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        
                        dx = current_pos[0] - self.drag_start_pos[0]
                        dy = current_pos[1] - self.drag_start_pos[1]
                        
                        if abs(dx) > 2 or abs(dy) > 2:
                            self.mouse_moved_during_drag = True
                    except Exception:
                        pass
                    
                    # ドラッグ中はイベントを消費してカメラ回転を防ぐ
                    return True
                
                # ホバー処理（ドラッグ中でない場合）
                if self.group_atoms:
                    try:
                        interactor = self.main_window.plotter.interactor
                        current_pos = interactor.GetEventPosition()
                        picker = self.main_window.plotter.picker
                        picker.Pick(current_pos[0], current_pos[1], 0, self.main_window.plotter.renderer)
                        
                        if picker.GetActor() is self.main_window.atom_actor:
                            picked_position = np.array(picker.GetPickPosition())
                            distances = np.linalg.norm(self.main_window.atom_positions_3d - picked_position, axis=1)
                            closest_atom_idx = np.argmin(distances)
                            
                            if closest_atom_idx in self.group_atoms:
                                self.main_window.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                            else:
                                self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                        else:
                            self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                    except Exception:
                        pass
                
                # ドラッグ中でない場合はカメラ回転を許可
                return False
            
            elif event.type() == QEvent.Type.MouseButtonRelease and event.button() == Qt.MouseButton.LeftButton:
                # マウスリリース時の処理
                if getattr(self, 'potential_drag', False) or (self.is_dragging_group and self.drag_start_pos):
                    try:
                        if self.is_dragging_group and self.mouse_moved_during_drag:
                            # ドラッグが実行された - CustomInteractorStyleに任せる（何もしない）
                            pass
                        else:
                            # マウスが閾値以下の移動 = 単なるクリック
                            # グループ内の原子をクリックした場合は選択/解除をトグル
                            if hasattr(self, 'clicked_atom_for_toggle'):
                                clicked_atom = self.clicked_atom_for_toggle
                                delattr(self, 'clicked_atom_for_toggle')
                                # ドラッグ状態をリセットしてからトグル処理
                                self.is_dragging_group = False
                                self.drag_start_pos = None
                                self.mouse_moved_during_drag = False
                                self.potential_drag = False
                                if hasattr(self, 'last_drag_positions'):
                                    delattr(self, 'last_drag_positions')
                                # トグル処理を実行
                                self.on_atom_picked(clicked_atom)
                                try:
                                    self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                                except Exception:
                                    pass
                                return True
                        
                    except Exception:
                        pass
                    finally:
                        # ドラッグ状態をリセット
                        self.is_dragging_group = False
                        self.drag_start_pos = None
                        self.mouse_moved_during_drag = False
                        self.potential_drag = False
                        # 保存していた位置情報をクリア
                        if hasattr(self, 'last_drag_positions'):
                            delattr(self, 'last_drag_positions')
                        try:
                            self.main_window.plotter.setCursor(Qt.CursorShape.ArrowCursor)
                        except Exception:
                            pass
                    
                    return True  # イベントを消費
                
                # ドラッグ中でない場合は通常のリリース処理
                return False
        
        # その他のイベントは親クラスに渡す
        return False
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときに、その原子が属する連結成分全体を選択（複数グループ対応）"""
        # ドラッグ中は選択を変更しない（ただしリリース時のトグルは許可）
        if getattr(self, 'is_dragging_group', False):
            return
        
        # BFS/DFSで連結成分を探索
        visited = set()
        queue = [atom_idx]
        visited.add(atom_idx)
        
        while queue:
            current_idx = queue.pop(0)
            for bond_idx in range(self.mol.GetNumBonds()):
                bond = self.mol.GetBondWithIdx(bond_idx)
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                
                if begin_idx == current_idx and end_idx not in visited:
                    visited.add(end_idx)
                    queue.append(end_idx)
                elif end_idx == current_idx and begin_idx not in visited:
                    visited.add(begin_idx)
                    queue.append(begin_idx)
        
        # 新しいグループとして追加または解除
        if visited.issubset(self.group_atoms):
            # すでに選択されている - 解除
            self.group_atoms -= visited
        else:
            # 新しいグループを追加
            self.group_atoms |= visited
        
        self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()
    
    def update_display(self):
        if not self.group_atoms:
            self.selection_label.setText("No group selected")
        else:
            atom_info = []
            for atom_idx in sorted(self.group_atoms):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_info.append(f"{symbol}({atom_idx})")
            
            self.selection_label.setText(f"Selected group: {len(self.group_atoms)} atoms - {', '.join(atom_info[:5])}{' ...' if len(atom_info) > 5 else ''}")
    
    def show_atom_labels(self):
        """選択されたグループの原子をハイライト表示（Ctrlクリックと同じスタイル）"""
        self.clear_atom_labels()
        
        if not self.group_atoms:
            return
        
        # 選択された原子のインデックスリストを作成
        selected_indices = list(self.group_atoms)
        
        # 選択された原子の位置を取得
        selected_positions = self.main_window.atom_positions_3d[selected_indices]
        
        # 原子の半径を少し大きくしてハイライト表示
        selected_radii = np.array([VDW_RADII.get(
            self.mol.GetAtomWithIdx(i).GetSymbol(), 0.4) * 1.3 
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
        
        # ハイライトアクターを追加して保存（ピッキング不可に設定）
        self.highlight_actor = self.main_window.plotter.add_mesh(
            highlight_glyphs, 
            color='yellow', 
            opacity=0.3, 
            name='move_group_highlight',
            pickable=False  # ピッキングを無効化
        )
        
        self.main_window.plotter.render()
    
    def clear_atom_labels(self):
        """原子ハイライトをクリア"""
        try:
            self.main_window.plotter.remove_actor('move_group_highlight')
        except Exception:
            pass
        
        if hasattr(self, 'highlight_actor'):
            try:
                self.main_window.plotter.remove_actor(self.highlight_actor)
            except Exception:
                pass
            self.highlight_actor = None
        
        try:
            self.main_window.plotter.render()
        except Exception:
            pass
    
    def reset_translation_inputs(self):
        """Translation入力フィールドをリセット"""
        self.x_trans_input.setText("0.0")
        self.y_trans_input.setText("0.0")
        self.z_trans_input.setText("0.0")
    
    def apply_translation(self):
        """選択したグループを並行移動"""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return
        
        try:
            dx = float(self.x_trans_input.text())
            dy = float(self.y_trans_input.text())
            dz = float(self.z_trans_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid translation values.")
            return
        
        translation_vector = np.array([dx, dy, dz])
        
        conf = self.mol.GetConformer()
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            new_pos = atom_pos + translation_vector
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # ラベルを再描画
        self.main_window.push_undo_state()
    
    def reset_rotation_inputs(self):
        """Rotation入力フィールドをリセット"""
        self.x_rot_input.setText("0.0")
        self.y_rot_input.setText("0.0")
        self.z_rot_input.setText("0.0")
    
    def apply_rotation(self):
        """選択したグループを回転"""
        if not self.group_atoms:
            QMessageBox.warning(self, "Warning", "Please select a group first.")
            return
        
        try:
            rx = float(self.x_rot_input.text())
            ry = float(self.y_rot_input.text())
            rz = float(self.z_rot_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid rotation values.")
            return
        
        # 度をラジアンに変換
        rx_rad = np.radians(rx)
        ry_rad = np.radians(ry)
        rz_rad = np.radians(rz)
        
        # グループの重心を計算
        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.group_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        centroid = np.mean(positions, axis=0)
        
        # 回転行列を作成
        # X軸周り
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx_rad), -np.sin(rx_rad)],
            [0, np.sin(rx_rad), np.cos(rx_rad)]
        ])
        # Y軸周り
        Ry = np.array([
            [np.cos(ry_rad), 0, np.sin(ry_rad)],
            [0, 1, 0],
            [-np.sin(ry_rad), 0, np.cos(ry_rad)]
        ])
        # Z軸周り
        Rz = np.array([
            [np.cos(rz_rad), -np.sin(rz_rad), 0],
            [np.sin(rz_rad), np.cos(rz_rad), 0],
            [0, 0, 1]
        ])
        
        # 合成回転行列 (Z * Y * X)
        R = Rz @ Ry @ Rx
        
        # 各原子を回転
        for atom_idx in self.group_atoms:
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            # 重心を原点に移動
            centered_pos = atom_pos - centroid
            # 回転
            rotated_pos = R @ centered_pos
            # 重心を元に戻す
            new_pos = rotated_pos + centroid
            conf.SetAtomPosition(atom_idx, new_pos.tolist())
            self.main_window.atom_positions_3d[atom_idx] = new_pos
        
        self.main_window.draw_molecule_3d(self.mol)
        self.main_window.update_chiral_labels()
        self.show_atom_labels()  # ラベルを再描画
        self.main_window.push_undo_state()
    
    def clear_selection(self):
        """選択をクリア"""
        self.selected_atoms.clear()
        self.group_atoms.clear()
        self.clear_atom_labels()
        self.update_display()
        # ドラッグ関連のフラグもリセット
        self.is_dragging_group = False
        self.drag_start_pos = None
        if hasattr(self, 'last_drag_positions'):
            delattr(self, 'last_drag_positions')
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().reject()
