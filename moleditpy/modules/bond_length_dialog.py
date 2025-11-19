from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton, QLineEdit, QWidget, QRadioButton, QMessageBox
)

from .dialog3_d_picking_mixin import Dialog3DPickingMixin

from PyQt6.QtCore import Qt
import numpy as np

class BondLengthDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.atom1_idx = None
        self.atom2_idx = None
        
        # 事前選択された原子を設定
        if preselected_atoms and len(preselected_atoms) >= 2:
            self.atom1_idx = preselected_atoms[0]
            self.atom2_idx = preselected_atoms[1]
        
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Adjust Bond Length")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click two atoms in the 3D view to select a bond, then specify the new length.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Current distance display
        self.distance_label = QLabel("")
        layout.addWidget(self.distance_label)
        
        # New distance input
        distance_layout = QHBoxLayout()
        distance_layout.addWidget(QLabel("New distance (Å):"))
        self.distance_input = QLineEdit()
        self.distance_input.setPlaceholderText("1.54")
        distance_layout.addWidget(self.distance_input)
        layout.addLayout(distance_layout)
        
        # Movement options
        group_box = QWidget()
        group_layout = QVBoxLayout(group_box)
        group_layout.addWidget(QLabel("Movement Options:"))
        
        self.atom1_fix_group_radio = QRadioButton("Atom 1: Fixed, Atom 2: Move connected group")
        self.atom1_fix_group_radio.setChecked(True)
        group_layout.addWidget(self.atom1_fix_group_radio)

        self.atom1_fix_radio = QRadioButton("Atom 1: Fixed, Atom 2: Move atom only")
        group_layout.addWidget(self.atom1_fix_radio)
        
        self.both_groups_radio = QRadioButton("Both groups: Move towards center equally")
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
        
        # Connect to main window's picker
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
        else:
            # Reset and start over
            self.atom1_idx = atom_idx
            self.atom2_idx = None
        
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
        self.atom2_idx = None
        self.clear_selection_labels()
        self.update_display()
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()
        
        # 新しいラベルを表示
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        selected_atoms = [self.atom1_idx, self.atom2_idx]
        labels = ["1st", "2nd"]
        colors = ["yellow", "yellow"]
        
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
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # Clear distance input when no selection
            try:
                self.distance_input.clear()
            except Exception:
                pass
        elif self.atom2_idx is None:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            self.selection_label.setText(f"First atom: {symbol1} (index {self.atom1_idx})")
            self.distance_label.setText("")
            self.apply_button.setEnabled(False)
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            # Clear distance input while selection is incomplete
            try:
                self.distance_input.clear()
            except Exception:
                pass
        else:
            symbol1 = self.mol.GetAtomWithIdx(self.atom1_idx).GetSymbol()
            symbol2 = self.mol.GetAtomWithIdx(self.atom2_idx).GetSymbol()
            self.selection_label.setText(f"Bond: {symbol1}({self.atom1_idx}) - {symbol2}({self.atom2_idx})")
            
            # Calculate current distance
            conf = self.mol.GetConformer()
            pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
            pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
            current_distance = np.linalg.norm(pos2 - pos1)
            self.distance_label.setText(f"Current distance: {current_distance:.3f} Å")
            self.apply_button.setEnabled(True)
            # Update the distance input box to show current distance
            try:
                self.distance_input.setText(f"{current_distance:.3f}")
            except Exception:
                pass
            # ラベル追加
            self.add_selection_label(self.atom1_idx, "1")
            self.add_selection_label(self.atom2_idx, "2")
    
    def apply_changes(self):
        """変更を適用"""
        if self.atom1_idx is None or self.atom2_idx is None:
            return
        
        try:
            new_distance = float(self.distance_input.text())
            if new_distance <= 0:
                QMessageBox.warning(self, "Invalid Input", "Distance must be positive.")
                return
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number.")
            return
        
        # Undo状態を保存
        self.main_window.push_undo_state()
        
        # Apply the bond length change
        self.adjust_bond_length(new_distance)
        
        # キラルラベルを更新
        self.main_window.update_chiral_labels()
    
    def adjust_bond_length(self, new_distance):
        """結合長を調整"""
        conf = self.mol.GetConformer()
        pos1 = np.array(conf.GetAtomPosition(self.atom1_idx))
        pos2 = np.array(conf.GetAtomPosition(self.atom2_idx))
        
        # Direction vector from atom1 to atom2
        direction = pos2 - pos1
        current_distance = np.linalg.norm(direction)
        
        if current_distance == 0:
            return
        
        direction = direction / current_distance
        
        if self.both_groups_radio.isChecked():
            # Both groups move towards center equally
            bond_center = (pos1 + pos2) / 2
            half_distance = new_distance / 2
            
            # New positions for both atoms
            new_pos1 = bond_center - direction * half_distance
            new_pos2 = bond_center + direction * half_distance
            
            # Get both connected groups
            group1_atoms = self.get_connected_group(self.atom1_idx, exclude=self.atom2_idx)
            group2_atoms = self.get_connected_group(self.atom2_idx, exclude=self.atom1_idx)
            
            # Calculate displacements
            displacement1 = new_pos1 - pos1
            displacement2 = new_pos2 - pos2
            
            # Move group 1
            for atom_idx in group1_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement1
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
            
            # Move group 2
            for atom_idx in group2_atoms:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement2
                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                self.main_window.atom_positions_3d[atom_idx] = new_pos
                
        elif self.atom1_fix_radio.isChecked():
            # Move only the second atom
            new_pos2 = pos1 + direction * new_distance
            conf.SetAtomPosition(self.atom2_idx, new_pos2.tolist())
            self.main_window.atom_positions_3d[self.atom2_idx] = new_pos2
        else:
            # Move the connected group (default behavior)
            new_pos2 = pos1 + direction * new_distance
            atoms_to_move = self.get_connected_group(self.atom2_idx, exclude=self.atom1_idx)
            displacement = new_pos2 - pos2
            
            for atom_idx in atoms_to_move:
                current_pos = np.array(conf.GetAtomPosition(atom_idx))
                new_pos = current_pos + displacement
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
