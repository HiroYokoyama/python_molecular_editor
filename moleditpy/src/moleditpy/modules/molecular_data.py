from rdkit import Chem
import traceback

try:
    from .constants import ANGSTROM_PER_PIXEL
except Exception:
    from modules.constants import ANGSTROM_PER_PIXEL

class MolecularData:
    def __init__(self):
        self.atoms = {}
        self.bonds = {}
        self._next_atom_id = 0
        self.adjacency_list = {} 

    def add_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self._next_atom_id
        self.atoms[atom_id] = {'symbol': symbol, 'pos': pos, 'item': None, 'charge': charge, 'radical': radical}
        self.adjacency_list[atom_id] = [] 
        self._next_atom_id += 1
        return atom_id

    def add_bond(self, id1, id2, order=1, stereo=0):
        # 立体結合の場合、IDの順序は方向性を意味するため、ソートしない。
        # 非立体結合の場合は、キーを正規化するためにソートする。
        if stereo == 0:
            if id1 > id2: id1, id2 = id2, id1

        bond_data = {'order': order, 'stereo': stereo, 'item': None}
        
        # 逆方向のキーも考慮して、新規結合かどうかをチェック
        is_new_bond = (id1, id2) not in self.bonds and (id2, id1) not in self.bonds
        if is_new_bond:
            if id1 in self.adjacency_list and id2 in self.adjacency_list:
                self.adjacency_list[id1].append(id2)
                self.adjacency_list[id2].append(id1)

        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)].update(bond_data)
            return (id1, id2), 'updated'
        else:
            self.bonds[(id1, id2)] = bond_data
            return (id1, id2), 'created'

    def remove_atom(self, atom_id):
        if atom_id in self.atoms:
            try:
                # Safely get neighbors before deleting the atom's own entry
                neighbors = self.adjacency_list.get(atom_id, [])
                for neighbor_id in neighbors:
                    if neighbor_id in self.adjacency_list and atom_id in self.adjacency_list[neighbor_id]:
                        self.adjacency_list[neighbor_id].remove(atom_id)

                # Now, safely delete the atom's own entry from the adjacency list
                if atom_id in self.adjacency_list:
                    del self.adjacency_list[atom_id]

                del self.atoms[atom_id]
                
                # Remove bonds involving this atom
                bonds_to_remove = [key for key in self.bonds if atom_id in key]
                for key in bonds_to_remove:
                    del self.bonds[key]
                    
            except Exception as e:
                print(f"Error removing atom {atom_id}: {e}")
                
                traceback.print_exc()

    def remove_bond(self, id1, id2):
        try:
            # 方向性のある立体結合(順方向/逆方向)と、正規化された非立体結合のキーを探す
            key_to_remove = None
            if (id1, id2) in self.bonds:
                key_to_remove = (id1, id2)
            elif (id2, id1) in self.bonds:
                key_to_remove = (id2, id1)

            if key_to_remove:
                if id1 in self.adjacency_list and id2 in self.adjacency_list[id1]:
                    self.adjacency_list[id1].remove(id2)
                if id2 in self.adjacency_list and id1 in self.adjacency_list[id2]:
                    self.adjacency_list[id2].remove(id1)
                del self.bonds[key_to_remove]
                
        except Exception as e:
            print(f"Error removing bond {id1}-{id2}: {e}")
            
            traceback.print_exc()


    def to_rdkit_mol(self, use_2d_stereo=True):
        """
        use_2d_stereo: Trueなら2D座標からE/Zを推定（従来通り）。FalseならE/Zラベル優先、ラベルがない場合のみ2D座標推定。
        3D変換時はuse_2d_stereo=Falseで呼び出すこと。
        """
        if not self.atoms:
            return None
        mol = Chem.RWMol()

        # --- Step 1: atoms ---
        atom_id_to_idx_map = {}
        for atom_id, data in self.atoms.items():
            atom = Chem.Atom(data['symbol'])
            atom.SetFormalCharge(data.get('charge', 0))
            atom.SetNumRadicalElectrons(data.get('radical', 0))
            atom.SetIntProp("_original_atom_id", atom_id)
            idx = mol.AddAtom(atom)
            atom_id_to_idx_map[atom_id] = idx

        # --- Step 2: bonds & stereo info保存（ラベル情報はここで保持） ---
        bond_stereo_info = {}  # bond_idx -> {'type': int, 'atom_ids': (id1,id2), 'bond_data': bond_data}
        for (id1, id2), bond_data in self.bonds.items():
            if id1 not in atom_id_to_idx_map or id2 not in atom_id_to_idx_map:
                continue
            idx1, idx2 = atom_id_to_idx_map[id1], atom_id_to_idx_map[id2]

            order_val = float(bond_data['order'])
            order = {1.0: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
                     2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE}.get(order_val, Chem.BondType.SINGLE)

            bond_idx = mol.AddBond(idx1, idx2, order) - 1

            # stereoラベルがあれば、bond_idxに対して詳細を保持（あとで使う）
            if 'stereo' in bond_data and bond_data['stereo'] in [1, 2, 3, 4]:
                bond_stereo_info[bond_idx] = {
                    'type': int(bond_data['stereo']),
                    'atom_ids': (id1, id2),
                    'bond_data': bond_data
                }

        # --- Step 3: sanitize ---
        final_mol = mol.GetMol()
        try:
            Chem.SanitizeMol(final_mol)
        except Exception:
            return None

        # --- Step 4: add 2D conformer ---
        # Convert from scene pixels to angstroms when creating RDKit conformer.
        conf = Chem.Conformer(final_mol.GetNumAtoms())
        conf.Set3D(False)
        for atom_id, data in self.atoms.items():
            if atom_id in atom_id_to_idx_map:
                idx = atom_id_to_idx_map[atom_id]
                pos = data.get('pos')
                if pos:
                    ax = pos.x() * ANGSTROM_PER_PIXEL
                    ay = -pos.y() * ANGSTROM_PER_PIXEL  # Y座標を反転（画面座標系→化学座標系）
                    conf.SetAtomPosition(idx, (ax, ay, 0.0))
        final_mol.AddConformer(conf)

        # --- Step 5: E/Zラベル優先の立体設定 ---
        # まず、E/Zラベルがあるbondを記録
        ez_labeled_bonds = set()
        for bond_idx, info in bond_stereo_info.items():
            if info['type'] in [3, 4]:
                ez_labeled_bonds.add(bond_idx)

        # 2D座標からE/Zを推定するのは、use_2d_stereo=True かつE/Zラベルがないbondのみ
        if use_2d_stereo:
            Chem.SetDoubleBondNeighborDirections(final_mol, final_mol.GetConformer(0))
        else:
            # 3D変換時: E/Zラベルがある場合は座標ベースの推定を完全に無効化
            if ez_labeled_bonds:
                # E/Zラベルがある場合は、すべての結合のBondDirをクリアして座標ベースの推定を無効化
                for b in final_mol.GetBonds():
                    b.SetBondDir(Chem.BondDir.NONE)
            else:
                # E/Zラベルがない場合のみ座標ベースの推定を実行
                Chem.SetDoubleBondNeighborDirections(final_mol, final_mol.GetConformer(0))

        # ヘルパー: 重原子優先で近傍を選ぶ
        def pick_preferred_neighbor(atom, exclude_idx):
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                if nbr.GetAtomicNum() > 1:
                    return nbr.GetIdx()
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() != exclude_idx:
                    return nbr.GetIdx()
            return None

        # --- Step 6: ラベルベースで上書き（E/Z を最優先） ---
        for bond_idx, info in bond_stereo_info.items():
            stereo_type = info['type']
            bond = final_mol.GetBondWithIdx(bond_idx)

            # 単結合の wedge/dash ラベル（1/2）がある場合
            if stereo_type in [1, 2]:
                if stereo_type == 1:
                    bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif stereo_type == 2:
                    bond.SetBondDir(Chem.BondDir.BEGINDASH)
                continue

            # 二重結合の E/Z ラベル（3/4）
            if stereo_type in [3, 4]:
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue

                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()

                bond_data = info.get('bond_data', {}) or {}
                stereo_atoms_specified = bond_data.get('stereo_atoms')

                if stereo_atoms_specified:
                    try:
                        a1_id, a2_id = stereo_atoms_specified
                        neigh1_idx = atom_id_to_idx_map.get(a1_id)
                        neigh2_idx = atom_id_to_idx_map.get(a2_id)
                    except Exception:
                        neigh1_idx = None
                        neigh2_idx = None
                else:
                    neigh1_idx = pick_preferred_neighbor(final_mol.GetAtomWithIdx(begin_atom_idx), end_atom_idx)
                    neigh2_idx = pick_preferred_neighbor(final_mol.GetAtomWithIdx(end_atom_idx), begin_atom_idx)

                if neigh1_idx is None or neigh2_idx is None:
                    continue

                bond.SetStereoAtoms(neigh1_idx, neigh2_idx)
                if stereo_type == 3:
                    bond.SetStereo(Chem.BondStereo.STEREOZ)
                elif stereo_type == 4:
                    bond.SetStereo(Chem.BondStereo.STEREOE)

                # 座標ベースでつけられた隣接単結合の BondDir（wedge/dash）がラベルと矛盾する可能性があるので消す
                b1 = final_mol.GetBondBetweenAtoms(begin_atom_idx, neigh1_idx)
                b2 = final_mol.GetBondBetweenAtoms(end_atom_idx, neigh2_idx)
                if b1 is not None:
                    b1.SetBondDir(Chem.BondDir.NONE)
                if b2 is not None:
                    b2.SetBondDir(Chem.BondDir.NONE)

        # Step 7: 最終化（キャッシュ更新 + 立体割当の再実行）
        final_mol.UpdatePropertyCache(strict=False)
        
        # 3D変換時（use_2d_stereo=False）でE/Zラベルがある場合は、force=Trueで強制適用
        if not use_2d_stereo and ez_labeled_bonds:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=True)
        else:
            Chem.AssignStereochemistry(final_mol, cleanIt=False, force=False)
        return final_mol

    def to_mol_block(self):
        try:
            mol = self.to_rdkit_mol()
            if mol:
                return Chem.MolToMolBlock(mol, includeStereo=True)
        except Exception:
            pass
        if not self.atoms: return None
        atom_map = {old_id: new_id for new_id, old_id in enumerate(self.atoms.keys())}
        num_atoms, num_bonds = len(self.atoms), len(self.bonds)
        mol_block = "\n  MoleditPy\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
        for old_id, atom in self.atoms.items():
            # Convert scene pixel coordinates to angstroms when emitting MOL block
            x_px = atom['item'].pos().x()
            y_px = -atom['item'].pos().y()
            x, y = x_px * ANGSTROM_PER_PIXEL, y_px * ANGSTROM_PER_PIXEL
            z, symbol = 0.0, atom['symbol']
            charge = atom.get('charge', 0)

            chg_code = 0
            if charge == 3: chg_code = 1
            elif charge == 2: chg_code = 2
            elif charge == 1: chg_code = 3
            elif charge == -1: chg_code = 5
            elif charge == -2: chg_code = 6
            elif charge == -3: chg_code = 7

            mol_block += f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0{chg_code:3d}  0  0  0  0  0  0  0\n"

        for (id1, id2), bond in self.bonds.items():
            idx1, idx2, order = atom_map[id1] + 1, atom_map[id2] + 1, bond['order']
            stereo_code = 0
            bond_stereo = bond.get('stereo', 0)
            if bond_stereo == 1:
                stereo_code = 1
            elif bond_stereo == 2:
                stereo_code = 6

            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}{stereo_code:3d}  0  0  0\n"
            
        mol_block += "M  END\n"
        return mol_block
