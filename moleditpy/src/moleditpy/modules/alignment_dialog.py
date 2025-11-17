from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton
)
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMessageBox
import numpy as np

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

class AlignmentDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, axis, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.axis = axis
        self.selected_atoms = set()
        
        # 事前選択された原子を追加（最大2個まで）
        if preselected_atoms:
            self.selected_atoms.update(preselected_atoms[:2])
        
        self.init_ui()
        
        # 事前選択された原子にラベルを追加
        if self.selected_atoms:
            for i, atom_idx in enumerate(sorted(self.selected_atoms), 1):
                self.add_selection_label(atom_idx, f"Atom {i}")
            self.update_display()
    
    def init_ui(self):
        axis_names = {'x': 'X-axis', 'y': 'Y-axis', 'z': 'Z-axis'}
        self.setWindowTitle(f"Align to {axis_names[self.axis]}")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel(f"Click atoms in the 3D view to select them for alignment to the {axis_names[self.axis]}. Exactly 2 atoms are required. The first atom will be moved to the origin, and the second atom will be positioned on the {axis_names[self.axis]}.")
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
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Alignment")
        self.apply_button.clicked.connect(self.apply_alignment)
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
        # Dialog3DPickingMixinの機能を使用
        super().enable_picking()
    
    def disable_picking(self):
        """3Dビューでの原子選択を無効にする"""
        # Dialog3DPickingMixinの機能を使用
        super().disable_picking()
    
    def on_atom_picked(self, atom_idx):
        """原子がクリックされた時の処理"""
        if self.main_window.current_mol is None:
            return
            
        if atom_idx in self.selected_atoms:
            # 既に選択されている場合は選択解除
            self.selected_atoms.remove(atom_idx)
            self.remove_atom_label(atom_idx)
        else:
            # 2つまでしか選択できない
            if len(self.selected_atoms) < 2:
                self.selected_atoms.add(atom_idx)
                # ラベルの順番を示す
                label_text = f"Atom {len(self.selected_atoms)}"
                self.add_selection_label(atom_idx, label_text)
        
        self.update_display()
    
    def update_display(self):
        """選択状態の表示を更新"""
        if len(self.selected_atoms) == 0:
            self.selection_label.setText("Click atoms to select for alignment (exactly 2 required)")
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 1:
            selected_list = list(self.selected_atoms)
            atom = self.mol.GetAtomWithIdx(selected_list[0])
            self.selection_label.setText(f"Selected 1 atom: {atom.GetSymbol()}{selected_list[0]+1}")
            self.apply_button.setEnabled(False)
        elif len(self.selected_atoms) == 2:
            selected_list = sorted(list(self.selected_atoms))
            atom1 = self.mol.GetAtomWithIdx(selected_list[0])
            atom2 = self.mol.GetAtomWithIdx(selected_list[1])
            self.selection_label.setText(f"Selected 2 atoms: {atom1.GetSymbol()}{selected_list[0]+1}, {atom2.GetSymbol()}{selected_list[1]+1}")
            self.apply_button.setEnabled(True)
    
    def clear_selection(self):
        """選択をクリア"""
        self.clear_selection_labels()
        self.selected_atoms.clear()
        self.update_display()
    
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
    
    def remove_atom_label(self, atom_idx):
        """特定の原子のラベルを削除"""
        # 簡単化のため、全ラベルをクリアして再描画
        self.clear_selection_labels()
        for i, idx in enumerate(sorted(self.selected_atoms), 1):
            if idx != atom_idx:
                self.add_selection_label(idx, f"Atom {i}")
    
    def clear_selection_labels(self):
        """選択ラベルをクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []
    
    def apply_alignment(self):
        """アライメントを適用"""
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Warning", "Please select exactly 2 atoms for alignment.")
            return
        try:

            selected_list = sorted(list(self.selected_atoms))
            atom1_idx, atom2_idx = selected_list[0], selected_list[1]

            conf = self.mol.GetConformer()

            # 原子の現在位置を取得
            pos1 = np.array(conf.GetAtomPosition(atom1_idx))
            pos2 = np.array(conf.GetAtomPosition(atom2_idx))

            # 最初に全分子を移動して、atom1を原点に配置
            translation = -pos1
            for i in range(self.mol.GetNumAtoms()):
                current_pos = np.array(conf.GetAtomPosition(i))
                new_pos = current_pos + translation
                conf.SetAtomPosition(i, new_pos.tolist())

            # atom2の新しい位置を取得（移動後）
            pos2_translated = pos2 + translation

            # atom2を選択した軸上に配置するための回転を計算
            axis_vectors = {
                'x': np.array([1.0, 0.0, 0.0]),
                'y': np.array([0.0, 1.0, 0.0]),
                'z': np.array([0.0, 0.0, 1.0])
            }
            target_axis = axis_vectors[self.axis]
            
            # atom2から原点への方向ベクトル
            current_vector = pos2_translated
            current_length = np.linalg.norm(current_vector)
            
            if current_length > 1e-10:  # ゼロベクトルでない場合
                current_vector_normalized = current_vector / current_length
                
                # 回転軸と角度を計算
                rotation_axis = np.cross(current_vector_normalized, target_axis)
                rotation_axis_length = np.linalg.norm(rotation_axis)
                
                if rotation_axis_length > 1e-10:  # 回転が必要
                    rotation_axis = rotation_axis / rotation_axis_length
                    cos_angle = np.dot(current_vector_normalized, target_axis)
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    rotation_angle = np.arccos(cos_angle)
                    
                    # ロドリゲスの回転公式を使用
                    def rodrigues_rotation(v, k, theta):
                        cos_theta = np.cos(theta)
                        sin_theta = np.sin(theta)
                        return (v * cos_theta + 
                               np.cross(k, v) * sin_theta + 
                               k * np.dot(k, v) * (1 - cos_theta))
                    
                    # 全ての原子に回転を適用
                    for i in range(self.mol.GetNumAtoms()):
                        current_pos = np.array(conf.GetAtomPosition(i))
                        rotated_pos = rodrigues_rotation(current_pos, rotation_axis, rotation_angle)
                        conf.SetAtomPosition(i, rotated_pos.tolist())
            
            # 3D座標を更新
            self.main_window.atom_positions_3d = np.array([
                list(conf.GetAtomPosition(i)) for i in range(self.mol.GetNumAtoms())
            ])
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Undo状態を保存
            self.main_window.push_undo_state()
            
            QMessageBox.information(self, "Success", f"Alignment to {self.axis.upper()}-axis completed.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply alignment: {str(e)}")
    
    def closeEvent(self, event):
        """ダイアログが閉じられる時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().closeEvent(event)
    
    def reject(self):
        """キャンセル時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().reject()
    
    def accept(self):
        """OK時の処理"""
        self.clear_selection_labels()
        self.disable_picking()
        super().accept()
