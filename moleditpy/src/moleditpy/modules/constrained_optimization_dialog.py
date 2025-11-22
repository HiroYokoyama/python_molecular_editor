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
    QDialog, QVBoxLayout, QLabel, QFormLayout, QComboBox, QLineEdit, QTableWidget, QTableWidgetItem, QHBoxLayout,
    QPushButton, QMessageBox, QAbstractItemView
)

from .dialog3_d_picking_mixin import Dialog3DPickingMixin

from PyQt6.QtCore import Qt
from rdkit.Chem import AllChem, rdMolTransforms


class ConstrainedOptimizationDialog(Dialog3DPickingMixin, QDialog):
    """制約付き最適化ダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        QDialog.__init__(self, parent)
        Dialog3DPickingMixin.__init__(self)
        self.mol = mol
        self.main_window = main_window
        self.selected_atoms = []  # 順序が重要なのでリストを使用
        self.constraints = []  # (type, atoms_indices, value)
        self.constraint_labels = [] # 3Dラベルアクター
        self.init_ui()
        self.enable_picking()

        # MainWindowから既存の制約を読み込む
        if self.main_window.constraints_3d:
            self.constraint_table.blockSignals(True) # 読み込み中のシグナルをブロック
            try:
                # self.constraints には (Type, (Idx...), Value, Force) のタプル形式で読み込む
                for const_data in self.main_window.constraints_3d:
                    # 後方互換性のため、3要素または4要素の制約に対応
                    if len(const_data) == 4:
                        const_type, atom_indices, value, force_const = const_data
                    else:
                        const_type, atom_indices, value = const_data
                        force_const = 1.0e5  # デフォルト値
                    
                    # タプル化して内部リストに追加
                    self.constraints.append((const_type, tuple(atom_indices), value, force_const))
                    
                    row_count = self.constraint_table.rowCount()
                    self.constraint_table.insertRow(row_count)
                    
                    value_str = ""
                    if const_type == "Distance":
                        value_str = f"{value:.3f}"
                    else:
                        value_str = f"{value:.2f}"

                    # カラム 0 (Type)
                    item_type = QTableWidgetItem(const_type)
                    item_type.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 0, item_type)

                    # カラム 1 (Atom Indices)
                    item_indices = QTableWidgetItem(str(atom_indices))
                    item_indices.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 1, item_indices)

                    # カラム 2 (Value)
                    item_value = QTableWidgetItem(value_str)
                    item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 2, item_value)
                    
                    # カラム 3 (Force)
                    item_force = QTableWidgetItem(f"{force_const:.2e}")
                    item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.constraint_table.setItem(row_count, 3, item_force)
            finally:
                self.constraint_table.blockSignals(False)

            # <<< MainWindowの現在の最適化設定を読み込み、デフォルトにする >>>
        try:
            # (修正) None の場合に備えてフォールバックを追加
            current_method_str = self.main_window.optimization_method or "MMFF_RDKIT"
            current_method = current_method_str.upper()
            
            # (修正) 比較順序を厳密化
            
            # 1. UFF_RDKIT
            if current_method == "UFF_RDKIT":
                self.ff_combo.setCurrentText("UFF")
            
            # 2. MMFF94_RDKIT (MMFF94)
            elif current_method == "MMFF94_RDKIT":
                self.ff_combo.setCurrentText("MMFF94")

            # 3. MMFF_RDKIT (MMFF94s) - これがデフォルトでもある
            elif current_method == "MMFF_RDKIT":
                self.ff_combo.setCurrentText("MMFF94s")

            # 4. (古い設定ファイルなどからのフォールバック)
            elif "UFF" in current_method:
                self.ff_combo.setCurrentText("UFF")
            elif "MMFF94S" in current_method:
                self.ff_combo.setCurrentText("MMFF94s")
            elif "MMFF94" in current_method: # MMFF94_RDKITも含むが、先で処理済み
                 self.ff_combo.setCurrentText("MMFF94")

            # 5. デフォルト
            else:
                self.ff_combo.setCurrentText("MMFF94s")
                
        except Exception as e:
            print(f"Could not set default force field: {e}")

    def init_ui(self):
        self.setWindowTitle("Constrained Optimization")
        self.setModal(False)
        self.resize(450, 500)
        self.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.WindowStaysOnTopHint)
        layout = QVBoxLayout(self)

        # 1. 説明
        instruction_label = QLabel("Select 2-4 atoms to add a constraint. Select constraints in the table to remove them.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # 2. 最適化方法とForce Constant
        form_layout = QFormLayout()
        self.ff_combo = QComboBox()
        self.ff_combo.addItems(["MMFF94s", "MMFF94", "UFF"])
        form_layout.addRow("Force Field:", self.ff_combo)
        
        # Force Constant設定
        self.force_const_input = QLineEdit("1.0e5")
        self.force_const_input.setToolTip("Force constant for constraints (default: 1.0e5)")
        form_layout.addRow("Force Constant:", self.force_const_input)
        
        layout.addLayout(form_layout)
        
        # 3. 選択中の原子
        self.selection_label = QLabel("Selected atoms: None")
        layout.addWidget(self.selection_label)

        # 4. 制約の表
        self.constraint_table = QTableWidget()
        self.constraint_table.setColumnCount(4)
        self.constraint_table.setHorizontalHeaderLabels(["Type", "Atom Indices", "Value (Å or °)", "Force"])
        self.constraint_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        # 編集トリガーをダブルクリックなどに変更
        self.constraint_table.setEditTriggers(QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.SelectedClicked | QTableWidget.EditTrigger.EditKeyPressed)
        self.constraint_table.itemSelectionChanged.connect(self.show_constraint_labels)
        self.constraint_table.cellChanged.connect(self.on_cell_changed)

        self.constraint_table.setStyleSheet("""
            QTableWidget QLineEdit {
                background-color: white;
                color: black;
                border: none;
            }
        """)

        layout.addWidget(self.constraint_table)

        # 5. ボタン (Add / Remove)
        button_layout = QHBoxLayout()
        self.add_button = QPushButton("Add Constraint")
        self.add_button.clicked.connect(self.add_constraint)
        self.add_button.setEnabled(False)
        button_layout.addWidget(self.add_button)
        
        self.remove_button = QPushButton("Remove Selected")
        self.remove_button.clicked.connect(self.remove_constraint)
        self.remove_button.setEnabled(False)
        button_layout.addWidget(self.remove_button)
        layout.addLayout(button_layout)

        self.remove_all_button = QPushButton("Remove All")
        self.remove_all_button.clicked.connect(self.remove_all_constraints)
        button_layout.addWidget(self.remove_all_button)
        
        # 6. メインボタン (Optimize / Close)
        main_buttons = QHBoxLayout()
        main_buttons.addStretch()
        self.optimize_button = QPushButton("Optimize")
        self.optimize_button.clicked.connect(self.apply_optimization)
        main_buttons.addWidget(self.optimize_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        main_buttons.addWidget(close_button)
        layout.addLayout(main_buttons)

    def on_atom_picked(self, atom_idx):
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            if len(self.selected_atoms) >= 4:
                self.selected_atoms.pop(0)  # 4つまで
            self.selected_atoms.append(atom_idx)
        
        self.show_selection_labels()
        self.update_selection_display()

    def update_selection_display(self):
        self.show_selection_labels()
        n = len(self.selected_atoms)

        atom_str = ", ".join(map(str, self.selected_atoms))
        prefix = ""
        can_add = False

        if n == 0:
            prefix = "Selected atoms: None"
            atom_str = ""  # atom_str を空にする
        elif n == 1:
            prefix = "Selected atoms: "
        elif n == 2:
            prefix = "Selected atoms: <b>Distance</b> "
            can_add = True
        elif n == 3:
            prefix = "Selected atoms: <b>Angle</b> "
            can_add = True
        elif n == 4:
            prefix = "Selected atoms: <b>Torsion</b> "
            can_add = True
        else: # n > 4
            prefix = "Selected atoms (max 4): "

        # ラベルテキストを設定
        if n == 0:
            self.selection_label.setText(prefix)
        else:
            self.selection_label.setText(f"{prefix}[{atom_str}]")

        # ボタンのテキストは常に固定
        self.add_button.setText("Add Constraint")
        # ボタンの有効状態を設定
        self.add_button.setEnabled(can_add)

    def add_constraint(self):
        n = len(self.selected_atoms)
        conf = self.mol.GetConformer()
        
        # Force Constantを取得
        try:
            force_const = float(self.force_const_input.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid Force Constant. Using default 1.0e5.")
            force_const = 1.0e5
        
        if n == 2:
            constraint_type = "Distance"
            value = conf.GetAtomPosition(self.selected_atoms[0]).Distance(conf.GetAtomPosition(self.selected_atoms[1]))
            value_str = f"{value:.3f}"
        elif n == 3:
            constraint_type = "Angle"
            value = rdMolTransforms.GetAngleDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        elif n == 4:
            constraint_type = "Torsion"
            value = rdMolTransforms.GetDihedralDeg(conf, *self.selected_atoms)
            value_str = f"{value:.2f}"
        else:
            return

        atom_indices = tuple(self.selected_atoms)
        
        # 既存の制約と重複チェック (原子インデックスが同じもの)
        for const in self.constraints:
            if const[0] == constraint_type and const[1] == atom_indices:
                QMessageBox.warning(self, "Warning", "This exact constraint already exists.")
                return

        self.constraints.append((constraint_type, atom_indices, value, force_const))
        
        # 表を更新
        # 表を更新
        row_count = self.constraint_table.rowCount()
        self.constraint_table.insertRow(row_count)

        # --- カラム 0 (Type) ---
        item_type = QTableWidgetItem(constraint_type)
        # 編集不可フラグを設定 (ItemIsEnabled | ItemIsSelectable)
        item_type.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_type.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 0, item_type)

        # --- カラム 1 (Atom Indices) ---
        item_indices = QTableWidgetItem(str(atom_indices))
        # 編集不可フラグを設定
        item_indices.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
        item_indices.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        self.constraint_table.setItem(row_count, 1, item_indices)

        # --- カラム 2 (Value) ---
        item_value = QTableWidgetItem(value_str)
        item_value.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # 編集可能フラグはデフォルトで有効 (ItemIsEnabled | ItemIsSelectable | ItemIsEditable)
        self.constraint_table.setItem(row_count, 2, item_value)
        
        # --- カラム 3 (Force) ---
        item_force = QTableWidgetItem(f"{force_const:.2e}")
        item_force.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        # 編集可能
        self.constraint_table.setItem(row_count, 3, item_force)

        # 選択をクリア
        self.selected_atoms.clear()
        self.update_selection_display()

    def remove_constraint(self):
        selected_rows = sorted(list(set(index.row() for index in self.constraint_table.selectedIndexes())), reverse=True)
        if not selected_rows:
            return

        self.constraint_table.blockSignals(True) 
        
        for row in selected_rows:
            self.constraints.pop(row)
            self.constraint_table.removeRow(row)
            
        self.constraint_table.blockSignals(False) 
            
        self.clear_constraint_labels()

    def remove_all_constraints(self):
        """全ての制約をクリアする"""
        if not self.constraints:
            return
            
        # 内部リストをクリア
        self.constraints.clear()
        
        # テーブルの行を全て削除
        self.constraint_table.blockSignals(True) 
        self.constraint_table.setRowCount(0)
        self.constraint_table.blockSignals(False) 
            
        # 3Dラベルをクリア
        self.clear_constraint_labels()
        
        # 選択ボタンを無効化
        self.remove_button.setEnabled(False)

    def show_constraint_labels(self):
        self.clear_constraint_labels()
        selected_items = self.constraint_table.selectedItems()
        if not selected_items:
            self.remove_button.setEnabled(False)
            return
            
        self.remove_button.setEnabled(True)
        
        # 選択された行の制約を取得 (最初の選択行のみ)
        try:
            row = selected_items[0].row()
            constraint_type, atom_indices, value, force_const = self.constraints[row]
        except (IndexError, TypeError, ValueError):
            # 古い形式の制約の場合は3要素でunpack
            try:
                constraint_type, atom_indices, value = self.constraints[row]
            except (IndexError, TypeError):
                return
        
        labels = []
        if constraint_type == "Distance":
            labels = ["A1", "A2"]
        elif constraint_type == "Angle":
            labels = ["A1", "A2 (V)", "A3"]
        elif constraint_type == "Torsion":
            labels = ["A1", "A2", "A3", "A4"]
        
        positions = []
        texts = []
        for i, atom_idx in enumerate(atom_indices):
            positions.append(self.main_window.atom_positions_3d[atom_idx])
            texts.append(labels[i])
        
        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions, texts,
                point_size=20, font_size=12, text_color='cyan', always_visible=True
            )
            self.constraint_labels.append(label_actor)

    def clear_constraint_labels(self):
        for label_actor in self.constraint_labels:
            try:
                self.main_window.plotter.remove_actor(label_actor)
            except Exception:
                pass
        self.constraint_labels = []

    def apply_optimization(self):
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No valid 3D molecule found.")
            return

        ff_name = self.ff_combo.currentText()
        conf = self.mol.GetConformer()
        
        try:
            if ff_name.startswith("MMFF"):
                props = AllChem.MMFFGetMoleculeProperties(self.mol, mmffVariant=ff_name)
                ff = AllChem.MMFFGetMoleculeForceField(self.mol, props, confId=0)
                add_dist_constraint = ff.MMFFAddDistanceConstraint
                add_angle_constraint = ff.MMFFAddAngleConstraint
                add_torsion_constraint = ff.MMFFAddTorsionConstraint
            else: # UFF
                ff = AllChem.UFFGetMoleculeForceField(self.mol, confId=0)
                add_dist_constraint = ff.UFFAddDistanceConstraint
                add_angle_constraint = ff.UFFAddAngleConstraint
                add_torsion_constraint = ff.UFFAddTorsionConstraint

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to initialize force field {ff_name}: {e}")
            return

        # 制約を追加
        try:
            for constraint in self.constraints:
                # 後方互換性のため、4要素または3要素の制約に対応
                if len(constraint) == 4:
                    const_type, atoms, value, force_const = constraint
                else:
                    const_type, atoms, value = constraint
                    force_const = 1.0e5  # デフォルト値
                
                if const_type == "Distance":
                    # C++ signature: (self, idx1, idx2, bool relative, minLen, maxLen, forceConst)
                    add_dist_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        False, 
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
                elif const_type == "Angle":
                    # C++ signature: (self, idx1, idx2, idx3, bool relative, minDeg, maxDeg, forceConst)
                    add_angle_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        int(atoms[2]),
                        False,  
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
                elif const_type == "Torsion":
                    # C++ signature: (self, idx1, idx2, idx3, idx4, bool relative, minDeg, maxDeg, forceConst)
                    add_torsion_constraint(
                        int(atoms[0]), 
                        int(atoms[1]), 
                        int(atoms[2]), 
                        int(atoms[3]),
                        False, 
                        float(value), 
                        float(value), 
                        float(force_const)
                    )
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to add constraints: {e}")
            print(e)
            return

        # 最適化の実行
        try:
            self.main_window.statusBar().showMessage(f"Running constrained {ff_name} optimization...")
            ff.Minimize(maxIts=20000)
            
            # 最適化後の座標をメインウィンドウの numpy 配列に反映
            for i in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                self.main_window.atom_positions_3d[i] = [pos.x, pos.y, pos.z]
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            self.main_window.update_chiral_labels()
            self.main_window.push_undo_state()
            self.main_window.statusBar().showMessage("Constrained optimization finished.")

            try:
                constrained_method_name = f"Constrained_{ff_name}"
                self.main_window.last_successful_optimization_method = constrained_method_name
            except Exception as e:
                print(f"Failed to set last_successful_optimization_method: {e}")

            # (修正) 最適化成功時にも制約リストをMainWindowに保存 (reject と同じロジック)
            try:
                # JSON互換のため、タプルをリストに変換して保存
                json_safe_constraints = []
                for const in self.constraints:
                    # 4要素の制約（Type, Indices, Value, Force）
                    if len(const) == 4:
                        json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                    else:
                        # 古い形式の場合は3要素にデフォルトのForceを追加
                        json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
                
                # 変更があった場合のみ MainWindow を更新
                if self.main_window.constraints_3d != json_safe_constraints:
                    self.main_window.constraints_3d = json_safe_constraints
                    self.main_window.has_unsaved_changes = True # 制約の変更も「未保存」扱い
                    self.main_window.update_window_title()
                    
            except Exception as e:
                print(f"Failed to save constraints post-optimization: {e}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Optimization failed: {e}")

    def closeEvent(self, event):
        self.reject()
        event.accept()
    
    def reject(self):
        self.clear_constraint_labels()
        self.clear_selection_labels()
        self.disable_picking()

        # ダイアログを閉じる際に現在の制約リストをMainWindowに保存
        try:
            # JSON互換のため、タプルをリストに変換して保存
            json_safe_constraints = []
            for const in self.constraints:
                # (Type, (Idx...), Value, Force) -> [Type, [Idx...], Value, Force]
                if len(const) == 4:
                    json_safe_constraints.append([const[0], list(const[1]), const[2], const[3]])
                else:
                    # 古い形式の場合は3要素にデフォルトのForceを追加
                    json_safe_constraints.append([const[0], list(const[1]), const[2], 1.0e5])
            
            # 変更があった場合のみ MainWindow を更新
            if self.main_window.constraints_3d != json_safe_constraints:
                self.main_window.constraints_3d = json_safe_constraints
                self.main_window.has_unsaved_changes = True # 制約の変更も「未保存」扱い
                self.main_window.update_window_title()
                
        except Exception as e:
            print(f"Failed to save constraints to main window: {e}")

        super().reject()

    def clear_selection(self):
        """選択をクリア (原子以外をクリックした時にMixinから呼ばれる)"""
        self.selected_atoms.clear()
        self.clear_selection_labels()
        self.update_selection_display()

    def show_selection_labels(self):
        """選択された原子にラベルを表示"""
        self.clear_selection_labels()
        
        if not hasattr(self, 'selection_labels'):
            self.selection_labels = []
        
        if not hasattr(self.main_window, 'atom_positions_3d') or self.main_window.atom_positions_3d is None:
            return  # 3D座標データがない場合は何もしない
            
        max_idx = len(self.main_window.atom_positions_3d) - 1
        positions = []
        texts = []
        
        for i, atom_idx in enumerate(self.selected_atoms):
            if atom_idx is not None and 0 <= atom_idx <= max_idx:
                positions.append(self.main_window.atom_positions_3d[atom_idx])
                texts.append(f"A{i+1}")
            elif atom_idx is not None:
                # インデックスが無効な場合はログ（デバッグ用）
                print(f"Warning: Invalid atom index {atom_idx} in show_selection_labels")
        
        if positions:
            label_actor = self.main_window.plotter.add_point_labels(
                positions, texts,
                point_size=20, font_size=12, text_color='yellow', always_visible=True
            )
            # add_point_labelsがリストを返す場合も考慮
            if isinstance(label_actor, list):
                self.selection_labels.extend(label_actor)
            else:
                self.selection_labels.append(label_actor)

    def clear_selection_labels(self):
        """選択ラベル(A1, A2...)をクリア"""
        if hasattr(self, 'selection_labels'):
            for label_actor in self.selection_labels:
                try:
                    self.main_window.plotter.remove_actor(label_actor)
                except Exception:
                    pass
            self.selection_labels = []

    def on_cell_changed(self, row, column):
        """テーブルのセルが編集されたときに内部データを更新する"""
        
        # "Value" 列 (カラムインデックス 2) と "Force" 列 (カラムインデックス 3) のみ対応
        if column not in [2, 3]:
            return

        try:
            # 変更されたアイテムからテキストを取得
            item = self.constraint_table.item(row, column)
            if not item:
                return
            
            new_value_str = item.text()
            new_value = float(new_value_str)
            
            # 内部の constraints リストを更新
            old_constraint = self.constraints[row]
            
            # 後方互換性のため、3要素または4要素の制約に対応
            if len(old_constraint) == 4:
                if column == 2:  # Value列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], new_value, old_constraint[3])
                elif column == 3:  # Force列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], old_constraint[2], new_value)
            else:
                # 古い3要素形式の場合
                if column == 2:  # Value列
                    self.constraints[row] = (old_constraint[0], old_constraint[1], new_value, 1.0e5)
                elif column == 3:  # Force列（新規追加）
                    self.constraints[row] = (old_constraint[0], old_constraint[1], old_constraint[2], new_value)

        except (ValueError, TypeError):
            # 不正な値（数値以外）が入力された場合
            # 元の値をテーブルに戻す
            self.constraint_table.blockSignals(True)
            
            if column == 2:  # Value列
                old_value = self.constraints[row][2]
                if self.constraints[row][0] == "Distance":
                    item.setText(f"{old_value:.3f}")
                else:
                    item.setText(f"{old_value:.2f}")
            elif column == 3:  # Force列
                old_force = self.constraints[row][3] if len(self.constraints[row]) == 4 else 1.0e5
                item.setText(f"{old_force:.2e}")
            
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.constraint_table.blockSignals(False)
            
            QMessageBox.warning(self, "Invalid Value", "Please enter a valid floating-point number.")
        except IndexError:
            # constraints リストとテーブルが同期していない場合（通常発生しない）
            pass
    
    def keyPressEvent(self, event):
        """キーボードイベントを処理 (Delete/Backspaceで制約を削除, Enterで最適化)"""
        key = event.key()
        
        # DeleteキーまたはBackspaceキーが押されたかチェック
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            # テーブルがフォーカスを持っているか、またはアイテムが選択されているか確認
            if self.constraint_table.hasFocus() or len(self.constraint_table.selectedIndexes()) > 0:
                self.remove_constraint()
                event.accept()
                return

        # Enter/Returnキーが押されたかチェック (最適化を実行)
        if key == Qt.Key.Key_Return or key == Qt.Key.Key_Enter:
            # テーブルが編集中でないことを確認（セルの編集中にEnterを押した場合）
            if self.constraint_table.state() != QAbstractItemView.State.EditingState:
                if self.optimize_button.isEnabled():
                    self.apply_optimization()
                event.accept()
                return
            
        # それ以外のキーはデフォルトの処理
        super().keyPressEvent(event)
