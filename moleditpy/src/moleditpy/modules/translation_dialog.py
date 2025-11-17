from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QGridLayout, QLineEdit, QCheckBox, QPushButton, QHBoxLayout
)
from PyQt6.QtWidgets import QMessageBox
from PyQt6.QtCore import Qt
import numpy as np

from .dialog3_d_picking_mixin import Dialog3DPickingMixin

class TranslationDialog(Dialog3DPickingMixin, QDialog):
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()  # 複数原子選択用
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Translation")
        self.setModal(False)  # モードレスにしてクリックを阻害しない
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)  # 常に前面表示
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Click atoms in the 3D view to select them. The centroid of selected atoms will be moved to the target coordinates, translating the entire molecule.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)
        
        # Coordinate inputs
        coord_layout = QGridLayout()
        coord_layout.addWidget(QLabel("Target X:"), 0, 0)
        self.x_input = QLineEdit("0.0")
        coord_layout.addWidget(self.x_input, 0, 1)
        
        coord_layout.addWidget(QLabel("Target Y:"), 1, 0)
        self.y_input = QLineEdit("0.0")
        coord_layout.addWidget(self.y_input, 1, 1)
        
        coord_layout.addWidget(QLabel("Target Z:"), 2, 0)
        self.z_input = QLineEdit("0.0")
        coord_layout.addWidget(self.z_input, 2, 1)
        
        layout.addLayout(coord_layout)

        # Translation target toggle: Entire molecule (default) or Selected atoms only
        self.translate_selected_only_checkbox = QCheckBox("Translate selected atoms only")
        self.translate_selected_only_checkbox.setToolTip(
            "When checked, only the atoms you selected will be moved so their centroid matches the target.\n"
            "When unchecked (default), the entire molecule will be translated so the selected atoms' centroid moves to the target."
        )
        self.translate_selected_only_checkbox.setChecked(False)  # default: entire molecule
        layout.addWidget(self.translate_selected_only_checkbox)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
    
        # Select all atoms button
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for translation")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)
        
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Translation")
        self.apply_button.clicked.connect(self.apply_translation)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
        
        # Connect to main window's picker
        self.picker_connection = None
        self.enable_picking()
    
    def on_atom_picked(self, atom_idx):
        """原子がピックされたときの処理"""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理"""
        if event.key() == Qt.Key.Key_Return or event.key() == Qt.Key.Key_Enter:
            if self.apply_button.isEnabled():
                self.apply_translation()
            event.accept()
        else:
            super().keyPressEvent(event)
    
    def update_display(self):
        """表示を更新"""
        if not self.selected_atoms:
            self.selection_label.setText("No atoms selected")
            self.apply_button.setEnabled(False)
        else:
            # 分子の有効性チェック
            if not self.mol or self.mol.GetNumConformers() == 0:
                self.selection_label.setText("Error: No valid molecule or conformer")
                self.apply_button.setEnabled(False)
                return
            
            try:
                conf = self.mol.GetConformer()
                # 選択原子の重心を計算
                centroid = self.calculate_centroid()
                
                # 選択原子の情報を表示
                atom_info = []
                for atom_idx in sorted(self.selected_atoms):
                    symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    atom_info.append(f"{symbol}({atom_idx})")
                
                self.selection_label.setText(
                    f"Selected atoms: {', '.join(atom_info)}\n"
                    f"Centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})"
                )
                self.apply_button.setEnabled(True)
            except Exception as e:
                self.selection_label.setText(f"Error accessing atom data: {str(e)}")
                self.apply_button.setEnabled(False)

        # Update the coordinate input fields when selection changes
        # If there are selected atoms, fill the inputs with the computed centroid (or single atom pos).
        # If no atoms selected, clear or reset to 0.0.
        try:
            if self.selected_atoms:
                # Use the centroid we just computed if available; otherwise compute now.
                try:
                    coords = centroid
                except NameError:
                    coords = self.calculate_centroid()

                # Format with reasonable precision
                self.x_input.setText(f"{coords[0]:.4f}")
                self.y_input.setText(f"{coords[1]:.4f}")
                self.z_input.setText(f"{coords[2]:.4f}")
            else:
                # No selection: reset fields to default
                self.x_input.setText("0.0")
                self.y_input.setText("0.0")
                self.z_input.setText("0.0")
        except Exception:
            # Be tolerant: do not crash the UI if inputs cannot be updated
            pass
    
    def calculate_centroid(self):
        """選択原子の重心を計算"""
        if not self.selected_atoms:
            return np.array([0.0, 0.0, 0.0])
        
        conf = self.mol.GetConformer()
        positions = []
        for atom_idx in self.selected_atoms:
            pos = conf.GetAtomPosition(atom_idx)
            positions.append([pos.x, pos.y, pos.z])
        
        return np.mean(positions, axis=0)
    
    def apply_translation(self):
        """平行移動を適用"""
        if not self.selected_atoms:
            QMessageBox.warning(self, "Warning", "Please select at least one atom.")
            return

        # 分子の有効性チェック
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Warning", "No valid molecule or conformer available.")
            return

        try:
            target_x = float(self.x_input.text())
            target_y = float(self.y_input.text())
            target_z = float(self.z_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Please enter valid coordinates.")
            return

        try:
            # 選択原子の重心を計算
            current_centroid = self.calculate_centroid()
            target_pos = np.array([target_x, target_y, target_z])

            # 移動ベクトルを計算
            translation_vector = target_pos - current_centroid

            conf = self.mol.GetConformer()

            if self.translate_selected_only_checkbox.isChecked():
                # Move only the selected atoms: shift selected atoms by translation_vector
                for i in range(self.mol.GetNumAtoms()):
                    if i in self.selected_atoms:
                        atom_pos = np.array(conf.GetAtomPosition(i))
                        new_pos = atom_pos + translation_vector
                        conf.SetAtomPosition(i, new_pos.tolist())
                        # Update 3d positions for this atom only
                        try:
                            self.main_window.atom_positions_3d[i] = new_pos
                        except Exception:
                            pass
                    else:
                        # leave other atoms unchanged
                        continue
            else:
                # Default: translate entire molecule so centroid moves to target
                for i in range(self.mol.GetNumAtoms()):
                    atom_pos = np.array(conf.GetAtomPosition(i))
                    new_pos = atom_pos + translation_vector
                    conf.SetAtomPosition(i, new_pos.tolist())
                    self.main_window.atom_positions_3d[i] = new_pos


            # 3D表示を更新
            self.main_window.draw_molecule_3d(self.mol)

            # キラルラベルを更新
            self.main_window.update_chiral_labels()

            # Apply後に選択解除
            self.clear_selection()

            # Undo状態を保存
            self.main_window.push_undo_state()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply translation: {str(e)}")
    
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
    
    def show_atom_labels(self):
        """選択された原子にラベルを表示"""
        # 既存のラベルをクリア
        self.clear_atom_labels()

        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []

        if self.selected_atoms:
            positions = []
            labels = []

            for i, atom_idx in enumerate(sorted(self.selected_atoms)):
                pos = self.main_window.atom_positions_3d[atom_idx]
                positions.append(pos)
                labels.append(f"S{i+1}")

            # 重心位置も表示
            if len(self.selected_atoms) > 1:
                centroid = self.calculate_centroid()
                positions.append(centroid)
                labels.append("CEN")

            # ラベルを追加
            if positions:
                label_actor = self.main_window.plotter.add_point_labels(
                    positions, labels,
                    point_size=20,
                    font_size=12,
                    text_color='cyan',
                    always_visible=True
                )
                # add_point_labelsがリストを返す場合も考慮
                if isinstance(label_actor, list):
                    self.selection_labels.extend(label_actor)
                else:
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
        # ラベル消去後に再描画を強制
        try:
            self.main_window.plotter.render()
        except Exception:
            pass
    
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
    
    def accept(self):
        """OK時の処理"""
        self.clear_atom_labels()
        self.disable_picking()
        try:
            self.main_window.draw_molecule_3d(self.mol)
        except Exception:
            pass
        super().accept()
