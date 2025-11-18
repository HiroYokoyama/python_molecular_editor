from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QHBoxLayout, QPushButton, QMessageBox
)
from PyQt6.QtCore import Qt
import numpy as np

try:
    from .dialog3_d_picking_mixin import Dialog3DPickingMixin
except Exception:
    from modules.dialog3_d_picking_mixin import Dialog3DPickingMixin

class PlanarizeDialog(Dialog3DPickingMixin, QDialog):

    """選択原子群を最適フィット平面へ投影して planarize するダイアログ
    AlignPlane を参考にした選択UIを持ち、Apply ボタンで選択原子を平面へ直交射影する。
    """
    def __init__(self, mol, main_window, preselected_atoms=None, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = set()

        if preselected_atoms:
            # 事前選択された原子を追加
            self.selected_atoms.update(preselected_atoms)

        self.init_ui()

        if self.selected_atoms:
            self.show_atom_labels()
            self.update_display()

    def init_ui(self):
        self.setWindowTitle("Planarize")
        self.setModal(False)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        layout = QVBoxLayout(self)

        instruction_label = QLabel("Click atoms in the 3D view to select them for planarization (minimum 3 required).")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)


        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        button_layout.addWidget(self.clear_button)
    
        # Select All Atoms ボタンを追加
        self.select_all_button = QPushButton("Select All Atoms")
        self.select_all_button.setToolTip("Select all atoms in the molecule for planarization")
        self.select_all_button.clicked.connect(self.select_all_atoms)
        button_layout.addWidget(self.select_all_button)

        self.apply_button = QPushButton("Apply planarize")
        self.apply_button.clicked.connect(self.apply_planarize)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
    
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        button_layout.addWidget(close_button)

        button_layout.addStretch()
    
        layout.addLayout(button_layout)

        # enable picking
        self.picker_connection = None
        self.enable_picking()

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.add(atom_idx)
        self.show_atom_labels()
        self.update_display()

    def clear_selection(self):
        self.selected_atoms.clear()
        self.clear_atom_labels()
        self.update_display()

    def update_display(self):
        count = len(self.selected_atoms)
        if count == 0:
            self.selection_label.setText("Click atoms to select for planarize (minimum 3 required)")
            self.apply_button.setEnabled(False)
        else:
            atom_list = sorted(self.selected_atoms)
            atom_display = []
            for i, atom_idx in enumerate(atom_list):
                symbol = self.mol.GetAtomWithIdx(atom_idx).GetSymbol()
                atom_display.append(f"#{i+1}: {symbol}({atom_idx})")
            self.selection_label.setText(f"Selected {count} atoms: {', '.join(atom_display)}")
            self.apply_button.setEnabled(count >= 3)

    def select_all_atoms(self):
        """Select all atoms in the current molecule (or fallback) and update labels/UI."""
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
        self.clear_atom_labels()
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        if self.selected_atoms:
            for i, atom_idx in enumerate(sorted(self.selected_atoms)):
                pos = self.main_window.atom_positions_3d[atom_idx]
                label_text = f"#{i+1}"
                label_actor = self.main_window.plotter.add_point_labels(
                    [pos], [label_text],
                    point_size=20,
                    font_size=12,
                    text_color='cyan',
                    always_visible=True
                )
                self.selection_labels.append(label_actor)

    def clear_atom_labels(self):
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except:
                    pass
            self.selection_labels = []

    def apply_planarize(self):
        if not self.selected_atoms or len(self.selected_atoms) < 3:
            QMessageBox.warning(self, "Warning", "Please select at least 3 atoms for planarize.")
            return

        try:
            selected_indices = list(sorted(self.selected_atoms))
            selected_positions = self.main_window.atom_positions_3d[selected_indices].copy()

            centroid = np.mean(selected_positions, axis=0)
            centered_positions = selected_positions - centroid

            # SVDによる最小二乗平面の法線取得
            u, s, vh = np.linalg.svd(centered_positions, full_matrices=False)
            normal = vh[-1]
            norm = np.linalg.norm(normal)
            if norm == 0:
                QMessageBox.warning(self, "Warning", "Cannot determine fit plane (degenerate positions).")
                return
            normal = normal / norm

            # 各点を重心を通る平面へ直交射影
            projections = centered_positions - np.outer(np.dot(centered_positions, normal), normal)
            new_positions = projections + centroid

            # 分子座標を更新
            conf = self.mol.GetConformer()
            for i, new_pos in zip(selected_indices, new_positions):
                conf.SetAtomPosition(int(i), new_pos.tolist())
                self.main_window.atom_positions_3d[int(i)] = new_pos

            # 3Dビュー更新
            self.main_window.draw_molecule_3d(self.mol)
            self.main_window.update_chiral_labels()
            self.main_window.push_undo_state()

            QMessageBox.information(self, "Success", f"Planarized {len(selected_indices)} atoms to best-fit plane.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to planarize: {e}")

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
