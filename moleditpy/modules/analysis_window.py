from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QGridLayout, QLabel, QLineEdit, QPushButton, QApplication
)
from PyQt6.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rd_inchi

class AnalysisWindow(QDialog):
    def __init__(self, mol, parent=None, is_xyz_derived=False):
        super().__init__(parent)
        self.mol = mol
        self.is_xyz_derived = is_xyz_derived  # XYZ由来かどうかのフラグ
        self.setWindowTitle("Molecule Analysis")
        self.setMinimumWidth(400)
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)
        grid_layout = QGridLayout()
        
        # --- 分子特性を計算 ---
        try:
            # RDKitのモジュールをインポート


            if self.is_xyz_derived:
                # XYZ由来の場合：元のXYZファイルの原子情報から直接計算
                # （結合推定の影響を受けない）
                
                # XYZファイルから読み込んだ元の原子情報を取得
                if hasattr(self.mol, '_xyz_atom_data'):
                    xyz_atoms = self.mol._xyz_atom_data
                else:
                    # フォールバック: RDKitオブジェクトから取得
                    xyz_atoms = [(atom.GetSymbol(), 0, 0, 0) for atom in self.mol.GetAtoms()]
                
                # 原子数と元素種を集計
                atom_counts = {}
                total_atoms = len(xyz_atoms)
                num_heavy_atoms = 0
                
                for symbol, x, y, z in xyz_atoms:
                    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                    if symbol != 'H':  # 水素以外
                        num_heavy_atoms += 1
                
                # 化学式を手動で構築（元素順序を考慮）
                element_order = ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
                formula_parts = []
                
                # 定義された順序で元素を追加
                remaining_counts = atom_counts.copy()
                for element in element_order:
                    if element in remaining_counts:
                        count = remaining_counts[element]
                        if count == 1:
                            formula_parts.append(element)
                        else:
                            formula_parts.append(f"{element}{count}")
                        del remaining_counts[element]
                
                # 残りの元素をアルファベット順で追加
                for element in sorted(remaining_counts.keys()):
                    count = remaining_counts[element]
                    if count == 1:
                        formula_parts.append(element)
                    else:
                        formula_parts.append(f"{element}{count}")
                
                mol_formula = ''.join(formula_parts)
                
                # 分子量と精密質量をRDKitから取得
                
                mol_wt = 0.0
                exact_mw = 0.0
                pt = Chem.GetPeriodicTable()
                
                for symbol, count in atom_counts.items():
                    try:
                        # RDKitの周期表から原子量と精密質量を取得
                        atomic_num = pt.GetAtomicNumber(symbol)
                        atomic_weight = pt.GetAtomicWeight(atomic_num)
                        exact_mass = pt.GetMostCommonIsotopeMass(atomic_num)
                        
                        mol_wt += atomic_weight * count
                        exact_mw += exact_mass * count
                    except (ValueError, RuntimeError):
                        # 認識されない元素の場合はスキップ
                        print(f"Warning: Unknown element {symbol}, skipping in mass calculation")
                        continue
                
                # 表示するプロパティを辞書にまとめる（XYZ元データから計算）
                properties = {
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Total Atoms:": str(total_atoms),
                }
                
                # 注意メッセージを追加
                note_label = QLabel("<i>Note: SMILES and structure-dependent properties are not available for XYZ-derived structures due to potential bond estimation inaccuracies.</i>")
                note_label.setWordWrap(True)
                main_layout.addWidget(note_label)
                
            else:
                # 通常の分子（MOLファイルや2Dエディタ由来）の場合：全てのプロパティを計算
                
                # SMILES生成用に、一時的に水素原子を取り除いた分子オブジェクトを作成
                mol_for_smiles = Chem.RemoveHs(self.mol)
                # 水素を取り除いた分子からSMILESを生成（常に簡潔な表記になる）
                smiles = Chem.MolToSmiles(mol_for_smiles, isomericSmiles=True)

                # 各種プロパティを計算
                mol_formula = rdMolDescriptors.CalcMolFormula(self.mol)
                mol_wt = Descriptors.MolWt(self.mol)
                exact_mw = Descriptors.ExactMolWt(self.mol)
                num_heavy_atoms = self.mol.GetNumHeavyAtoms()
                num_rings = rdMolDescriptors.CalcNumRings(self.mol)
                log_p = Descriptors.MolLogP(self.mol)
                tpsa = Descriptors.TPSA(self.mol)
                num_h_donors = rdMolDescriptors.CalcNumHBD(self.mol)
                num_h_acceptors = rdMolDescriptors.CalcNumHBA(self.mol)
                
                # InChIを生成
                try:
                    inchi = Chem.MolToInchi(self.mol)
                except Exception:
                    inchi = "N/A"

                # InChIKeyを生成（RDKitのinchi APIが無い場合に備えてフォールバック）
                try:
                    # Prefer Chem.MolToInchiKey when available
                    inchi_key = None
                    try:
                        inchi_key = Chem.MolToInchiKey(self.mol)
                    except Exception:
                        # Fallback to rdkit.Chem.inchi if present
                        try:
                            inchi_key = rd_inchi.MolToInchiKey(self.mol)
                        except Exception:
                            inchi_key = None

                    if not inchi_key:
                        inchi_key = "N/A"
                except Exception:
                    inchi_key = "N/A"

                # 表示するプロパティを辞書にまとめる
                properties = {
                    "SMILES:": smiles,
                    "InChI:": inchi,
                    "InChIKey:": inchi_key,
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Ring Count:": str(num_rings),
                    "LogP (o/w):": f"{log_p:.3f}",
                    "TPSA (Å²):": f"{tpsa:.2f}",
                    "H-Bond Donors:": str(num_h_donors),
                    "H-Bond Acceptors:": str(num_h_acceptors),
                }
        except Exception as e:
            main_layout.addWidget(QLabel(f"Error calculating properties: {e}"))
            return

        # --- 計算結果をUIに表示 ---
        row = 0
        for label_text, value_text in properties.items():
            label = QLabel(f"<b>{label_text}</b>")
            value = QLineEdit(value_text)
            value.setReadOnly(True)
            
            copy_btn = QPushButton("Copy")
            copy_btn.clicked.connect(lambda _, v=value: self.copy_to_clipboard(v.text()))

            grid_layout.addWidget(label, row, 0)
            grid_layout.addWidget(value, row, 1)
            grid_layout.addWidget(copy_btn, row, 2)
            row += 1
            
        main_layout.addLayout(grid_layout)
        
        # --- OKボタン ---
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        main_layout.addWidget(ok_button, 0, Qt.AlignmentFlag.AlignCenter)
        
        self.setLayout(main_layout)

    def copy_to_clipboard(self, text):
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        if self.parent() and hasattr(self.parent(), 'statusBar'):
            self.parent().statusBar().showMessage(f"Copied '{text}' to clipboard.", 2000)
