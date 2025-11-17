from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QButtonGroup, QRadioButton,
    QHBoxLayout, QPushButton, QMessageBox
)
from rdkit import Chem

class MirrorDialog(QDialog):
    """分子の鏡像を作成するダイアログ"""
    
    def __init__(self, mol, main_window, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.main_window = main_window
        self.init_ui()
    
    def init_ui(self):
        self.setWindowTitle("Mirror Molecule")
        self.setMinimumSize(300, 200)
        
        layout = QVBoxLayout(self)
        
        # 説明テキスト
        info_label = QLabel("Select the mirror plane to create molecular mirror image:")
        layout.addWidget(info_label)
        
        # ミラー平面選択のラジオボタン
        self.plane_group = QButtonGroup(self)
        
        self.xy_radio = QRadioButton("XY plane (Z = 0)")
        self.xz_radio = QRadioButton("XZ plane (Y = 0)")
        self.yz_radio = QRadioButton("YZ plane (X = 0)")
        
        self.xy_radio.setChecked(True)  # デフォルト選択
        
        self.plane_group.addButton(self.xy_radio, 0)
        self.plane_group.addButton(self.xz_radio, 1)
        self.plane_group.addButton(self.yz_radio, 2)
        
        layout.addWidget(self.xy_radio)
        layout.addWidget(self.xz_radio)
        layout.addWidget(self.yz_radio)
        
        layout.addSpacing(20)
        
        # ボタン
        button_layout = QHBoxLayout()
        
        apply_button = QPushButton("Apply Mirror")
        apply_button.clicked.connect(self.apply_mirror)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)

        button_layout.addWidget(apply_button)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)
    
    def apply_mirror(self):
        """選択された平面に対してミラー変換を適用"""
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No 3D coordinates available.")
            return
        
        # 選択された平面を取得
        plane_id = self.plane_group.checkedId()
        
        try:
            conf = self.mol.GetConformer()
            
            # 各原子の座標を変換
            for atom_idx in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(atom_idx)
                
                if plane_id == 0:  # XY平面（Z軸に対してミラー）
                    new_pos = [pos.x, pos.y, -pos.z]
                elif plane_id == 1:  # XZ平面（Y軸に対してミラー）
                    new_pos = [pos.x, -pos.y, pos.z]
                elif plane_id == 2:  # YZ平面（X軸に対してミラー）
                    new_pos = [-pos.x, pos.y, pos.z]
                
                # 新しい座標を設定
                from rdkit.Geometry import Point3D
                conf.SetAtomPosition(atom_idx, Point3D(new_pos[0], new_pos[1], new_pos[2]))
            
            # 3Dビューを更新
            self.main_window.draw_molecule_3d(self.mol)
            
            # ミラー変換後にキラルタグを強制的に再計算
            try:
                if self.mol.GetNumConformers() > 0:
                    # 既存のキラルタグをクリア
                    for atom in self.mol.GetAtoms():
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                    # 3D座標から新しいキラルタグを計算
                    Chem.AssignAtomChiralTagsFromStructure(self.mol, confId=0)
            except Exception as e:
                print(f"Error updating chiral tags: {e}")
            
            # キラルラベルを更新（鏡像変換でキラリティが変わる可能性があるため）
            self.main_window.update_chiral_labels()
            
            self.main_window.push_undo_state()
            
            plane_names = ["XY", "XZ", "YZ"]
            self.main_window.statusBar().showMessage(f"Molecule mirrored across {plane_names[plane_id]} plane.")
        

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply mirror transformation: {str(e)}")
