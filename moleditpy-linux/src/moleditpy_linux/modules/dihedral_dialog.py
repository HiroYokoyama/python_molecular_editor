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
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton, QLineEdit, QWidget, QRadioButton
)

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

from PyQt6.QtCore import Qt
import numpy as np
from PyQt6.QtWidgets import QMessageBox


class DihedralDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 4:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # central bond start
            self.atom3_idx = preselected_atoms[2]  # central bond end
            self.atom4_idx = preselected_atoms[3]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Dihedral Angle")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click four atoms in order to define a dihedral angle. The rotation will be around the bond between the 2nd and 3rd atoms.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current dihedral angle display
        self.dihedral_label = QLabel("")
        layout.addWidget(self.dihedral_label)
        
        # New dihedral angle input
        dihedral_layout = QHBoxLayout()
        dihedral_layout.addWidget(QLabel("New dihedral angle (degrees):"))
        self.dihedral_input = QLineEdit()
        self.dihedral_input.setPlaceholderText("180.0")
        dihedral_layout.addWidget(self.dihedral_input)
        layout.addLayout(dihedral_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Move:"))
        
        self.move_group_radio = QRadioButton("Atom 1,2,3: Fixed, Atom 4 group: Rotate")
        self.move_group_radio.setChecked(True)
        group_layout.addWidget(self.move_group_radio)
        
        self.move_atom_radio = QRadioButton("Atom 1,2,3: Fixed, Atom 4: Rotate atom only")
        group_layout.addWidget(self.move_atom_radio)
        
        self.both_groups_radio = QRadioButton("Central bond fixed: Both groups rotate equally")
        group_layout.addWidget(self.both_groups_radio)
        
        layout.addWidget(group_box)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker for DihedralDialog
        self.picker_connection = None
        self.enable_picking()
        
        # 事前選択された原子がある場合は初期表示を更新
        if self.atom1_idx is not None:
            self.show_atom_labels()
            self.update_display()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        elif self.atom3_idx is None:
            self.atom3_idx = atom_idx
        elif self.atom4_idx is None:
            self.atom4_idx = atom_idx
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
            self.atom4_idx = None
        
        # 原子ラベルを表示
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_changes()
            event.accept()
        else:
            super().keyPressEvent(event)
    
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
    
    def clear_selection(self):
        """選択をクリア"""
        self.atom1_idx = None
        self.atom2_idx = None  # central bond start
        self.atom3_idx = None  # central bond end
        self.atom4_idx = None
        self.clear_atom_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
        labels = ["1st", "2nd (bond start)", "3rd (bond end)", "4th"]
        colors = ["yellow", "yellow", "yellow", "yellow"]  # 全て黄色に統一
        
        for i, atom_idx in enumerate(selected_atoms):
            if atom_idx is not None:
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"{labels[i]}"
                
                # ラベルを追加
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text], 
                    point_size=20, 
                    font_size=12,
                    text_color=colors[i],
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
    
    def update_display(self):
        """表示を更新"""
        selected_count = sum(x is not None for x in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx])
        
        if selected_count == 0:
            self.selection_label.setText("No atoms selected")
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input when no selection
            try:
                self.dihedral_input.clear()
            except Exception:
                pass
        elif selected_count < 4:
            selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
            
            display_parts = []
            for atom_idx in selected_atoms:
                if atom_idx is not None:
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    display_parts.append(f"{symbol}({atom_idx})")
                else:
                    display_parts.append("?")
            
            self.selection_label.setText(" - ".join(display_parts))
            self.dihedral_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear dihedral input while selection is incomplete
            try:
                self.dihedral_input.clear()
            except Exception:
                pass
        else:
            selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]
            
            display_parts = []
            for atom_idx in selected_atoms:
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                display_parts.append(f"{symbol}({atom_idx})")
            
            self.selection_label.setText(" - ".join(display_parts))
            
            # Calculate current dihedral angle
            current_dihedral = self.calculate_dihedral()
            self.dihedral_label.setText(f"Current dihedral: {current_dihedral:.2f}°")
            self.apply_button.setEnabled(True)
            # Update dihedral input box with current dihedral
            try:
                self.dihedral_input.setText(f"{current_dihedral:.2f}")
            except Exception:
                pass
    
    def calculate_dihedral(self):
        """現在の二面角を計算（正しい公式を使用）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        pos4 = np.array(conf.GetAtomPosition(self.atom4_idx))
        
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
    
    def apply_changes(self):
        """変更を適用"""
        if any(idx is None for idx in [self.atom1_idx, self.atom2_idx, self.atom3_idx, self.atom4_idx]):
            return
        
        try:
            new_dihedral = float(self.dihedral_input.text())
            if new_dihedral < -180 or new_dihedral > 180:
                QMessageBox.warning(self, "Invalid Input", "Dihedral angle must be between -180 and 180 degrees.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Apply the dihedral angle change
        self.adjust_dihedral(new_dihedral)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()

        # Undo状態を保存
        self.main_window.push_undo_state()
    
    def adjust_dihedral(self, new_dihedral_deg):
        """二面角を調整（改善されたアルゴリズム）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        pos4 = np.array(conf.GetAtomPosition(self.atom4_idx))
        
        # Current dihedral angle
        current_dihedral = self.calculate_dihedral()
        
        # Calculate rotation angle needed
        rotation_angle_deg = new_dihedral_deg - current_dihedral
        
        # Handle angle wrapping for shortest rotation
        if rotation_angle_deg > 180:
            rotation_angle_deg -= 360
        elif rotation_angle_deg < -180:
            rotation_angle_deg += 360
        
        rotation_angle_rad = np.radians(rotation_angle_deg)
        
        # Skip if no rotation needed
        if abs(rotation_angle_rad) < 1e-6:
            return
        
        # Rotation axis is the bond between atom2 and atom3
        rotation_axis = pos3 - pos2
        axis_length = np.linalg.norm(rotation_axis)
        
        if axis_length == 0:
            return  # Atoms are at the same position
        
        rotation_axis = rotation_axis / axis_length
        
        # Rodrigues' rotation formula implementation
        def rotate_point_around_axis(point, axis_point, axis_direction, angle):
            """Rotate a point around an axis using Rodrigues' formula"""
            # Translate point so axis passes through origin
            translated_point = point - axis_point
            
            # Apply Rodrigues' rotation formula
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            
            rotated = (translated_point * cos_a + 
                      np.cross(axis_direction, translated_point) * sin_a + 
                      axis_direction * np.dot(axis_direction, translated_point) * (1 - cos_a))
            
            # Translate back to original coordinate system
            return rotated + axis_point
        
        if self.both_groups_radio.isChecked():
            # Both groups rotate equally around the central bond (half angle each in opposite directions)
            half_rotation = rotation_angle_rad / 2
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom2_idx, exclude=self.atom3_idx)
            group4_atoms = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate group1 (atom1 side) by -half_rotation
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, -half_rotation)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Rotate group4 (atom4 side) by +half_rotation
            for atom_idx in group4_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, half_rotation)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.move_group_radio.isChecked():
            # Move the connected group containing atom4
            # Find all atoms connected to atom3 (excluding atom2 side)
            atoms_to_rotate = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate all atoms in the group
            for atom_idx in atoms_to_rotate:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = rotate_point_around_axis(current_pos, pos2, rotation_axis, rotation_angle_rad)
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        else:
            # Move only atom4
            new_pos4 = rotate_point_around_axis(pos4, pos2, rotation_axis, rotation_angle_rad)
            conf.SetAtomPosition(self.atom4_idx, new_pos4.tolist())
            self.main_window.atom_positions_3d[self.atom4_idx] = new_pos4
        
        # Update the 3D view
        self.main_window.draw_molecule_3d(self.mol)
    
    def get_connected_group(self, start_atom, exclude=None):
        """指定された原子から連結されているグループを取得"""
        visited = set()
        to_visit = [start_atom]
        
        while to_visit:
            current = to_visit.pop()
            if current in visited or current == exclude:
                continue
            
            visited.add(current)
            
            # Get neighboring atoms
            atom = self.mol.GetAtomWithIdx(current)
            for bond in atom.GetBonds():
                other_idx = bond.GetOtherAtomIdx(current)
                if other_idx not in visited and other_idx != exclude:
                    to_visit.append(other_idx)
        
        return visited
