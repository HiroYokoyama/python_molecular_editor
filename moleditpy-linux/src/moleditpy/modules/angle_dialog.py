from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton, QLineEdit, QWidget, QRadioButton
)

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMessageBox
import numpy as np


class AngleDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 3:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]  # vertex
            self.atom3_idx = preselected_atoms[2]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Angle")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click three atoms in order: first-vertex-third. The angle around the vertex atom will be adjusted.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current angle display
        self.angle_label = QLabel("")
        layout.addWidget(self.angle_label)
        
        # New angle input
        angle_layout = QHBoxLayout()
        angle_layout.addWidget(QLabel("New angle (degrees):"))
        self.angle_input = QLineEdit()
        self.angle_input.setPlaceholderText("109.5")
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Rotation Options:"))
        
        self.rotate_group_radio = QRadioButton("Atom 1,2: Fixed, Atom 3: Rotate connected group")
        self.rotate_group_radio.setChecked(True)
        group_layout.addWidget(self.rotate_group_radio)

        self.rotate_atom_radio = QRadioButton("Atom 1,2: Fixed, Atom 3: Rotate atom only")
        group_layout.addWidget(self.rotate_atom_radio)
        
        self.both_groups_radio = QRadioButton("Vertex fixed: Both arms rotate equally")
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
        
        # Connect to main window's picker for AngleDialog
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
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
            self.atom3_idx = None
        
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
        self.atom2_idx = None  # vertex atom
        self.atom3_idx = None
        self.clear_selection_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx, self.atom3_idx]
        labels = ["1st", "2nd (vertex)", "3rd"]
        colors = ["yellow", "yellow", "yellow"]  # 全て黄色に統一
        
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
                except:
                    pass
            self.selection_labels = []
    
    def clear_selection_labels(self):
        """選択ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def add_selection_label(self, atom_idx, label_text):
        """選択された原子にラベルを追加"""
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        # 原子の位置を取得
        pos = self.main_window.atom_positions_3d[atom_idx]
        
        # ラベルを追加
        label_actor = self.main_window.plotter.add_point_labels(
            [pos], [label_text], 
            point_size=20, 
            font_size=12,
            text_color='yellow',
            always_visible=True
        )
        self.selection_labels.append(label_actor)
    
    def update_display(self):
        """表示を更新"""
        # 既存のラベルをクリア
        self.clear_selection_labels()
        
        if self.atom1_idx is None:
            self.selection_label.setText("No atoms selected")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear angle input when no selection
            try:
                self.angle_input.clear()
            except Exception:
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.clear()
            except Exception:
                pass
        elif self.atom3_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(f"Selected: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - ?")
            self.angle_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            # Clear angle input while selection is incomplete
            try:
                self.angle_input.clear()
            except Exception:
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            symbol3 = self.mol.GetAtomWithIdx(self.atom3_idx).GetSymbol()
            self.selection_label.setText(f"Angle: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx}) - {symbol3}({self.atom3_idx})")
            
            # Calculate current angle
            current_angle = self.calculate_angle()
            self.angle_label.setText(f"Current angle: {current_angle:.2f}°")
            self.apply_button.setEnabled(True)
            # Update angle input box with current angle
            try:
                self.angle_input.setText(f"{current_angle:.2f}")
            except Exception:
                pass
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2(vertex)")
            self.add_selection_label(self.atom3_idx, "3")
    
    def calculate_angle(self):
        """現在の角度を計算"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))  # vertex
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)
    
    def apply_changes(self):
        """変更を適用"""
        if self.atom1_idx is None or self.atom2_idx is None or self.atom3_idx is None:
            return
        
        try:
            new_angle = float(self.angle_input.text())
            if new_angle < 0 or new_angle >= 360:
                QMessageBox.warning(self, "Invalid Input", "Angle must be between 0 and 360 degrees.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Undo状態を保存
        self.main_window.push_undo_state()
        
        # Apply the angle change
        self.adjust_angle(new_angle)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()
    
    def adjust_angle(self, new_angle_deg):
        """角度を調整（均等回転オプション付き）"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))  # vertex
        pos3 = np.array(conf.GetAtomPosition(self.atom3_idx))
        
        vec1 = pos1 - pos2
        vec2 = pos3 - pos2
        
        # Current angle
        current_angle_rad = np.arccos(np.clip(
            np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)), -1.0, 1.0))
        
        # Target angle
        target_angle_rad = np.radians(new_angle_deg)
        
        # Rotation axis (perpendicular to the plane containing vec1 and vec2)
        rotation_axis = np.cross(vec1, vec2)
        rotation_axis_norm = np.linalg.norm(rotation_axis)
        
        if rotation_axis_norm == 0:
            # Vectors are parallel, cannot rotate
            return
        
        rotation_axis = rotation_axis / rotation_axis_norm
        
        # Total rotation angle needed
        total_rotation_angle = target_angle_rad - current_angle_rad
        
        # Rodrigues' rotation formula
        def rotate_vector(v, axis, angle):
            cos_a = np.cos(angle)
            sin_a = np.sin(angle)
            return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)
        
        if self.both_groups_radio.isChecked():
            # Both arms rotate equally (half angle each in opposite directions)
            half_rotation = total_rotation_angle / 2
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom1_idx, exclude=self.atom2_idx)
            group3_atoms = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            # Rotate group 1 by -half_rotation
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                relative_pos = current_pos - pos2
                rotated_pos = rotate_vector(relative_pos, rotation_axis, -half_rotation)
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Rotate group 3 by +half_rotation
            for atom_idx in group3_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                relative_pos = current_pos - pos2
                rotated_pos = rotate_vector(relative_pos, rotation_axis, half_rotation)
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.rotate_atom_radio.isChecked():
            # Move only the third atom
            new_vec2 = rotate_vector(vec2, rotation_axis, total_rotation_angle)
            new_pos3 = pos2 + new_vec2
            conf.SetAtomPosition(self.atom3_idx, new_pos3.tolist())
            self.main_window.atom_positions_3d[self.atom3_idx] = new_pos3
        else:
            # Rotate the connected group around atom2 (vertex) - default behavior
            atoms_to_move = self.get_connected_group(self.atom3_idx, exclude=self.atom2_idx)
            
            for atom_idx in atoms_to_move:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                # Transform to coordinate system centered at atom2
                relative_pos = current_pos - pos2
                # Rotate around the rotation axis
                rotated_pos = rotate_vector(relative_pos, rotation_axis, total_rotation_angle)
                # Transform back to world coordinates
                new_pos = pos2 + rotated_pos
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
        
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
