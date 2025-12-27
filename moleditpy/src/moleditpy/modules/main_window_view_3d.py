#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
main_window_view_3d.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowView3d
"""


import numpy as np
import vtk
import logging


# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QGraphicsView
)

from PyQt6.QtGui import (
    QColor, QTransform
)


from PyQt6.QtCore import (
    Qt, QRectF
)

import pyvista as pv

# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        logging.warning("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None

# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import CPK_COLORS_PV, DEFAULT_CPK_COLORS, VDW_RADII, pt
    from .template_preview_item import TemplatePreviewItem
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import CPK_COLORS_PV, DEFAULT_CPK_COLORS, VDW_RADII, pt
    from modules.template_preview_item import TemplatePreviewItem

# --- クラス定義 ---
class MainWindowView3d(object):
    """ main_window.py から分離された機能クラス """


    def set_3d_style(self, style_name):
        """3D表示スタイルを設定し、ビューを更新する"""
        if self.current_3d_style == style_name:
            return

        # 描画モード変更時に測定モードと3D編集モードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化
        
        # 3D原子選択をクリア
        self.clear_3d_selection()

        self.current_3d_style = style_name
        self.statusBar().showMessage(f"3D style set to: {style_name}")
        
        # 現在表示中の分子があれば、新しいスタイルで再描画する
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)



    def draw_molecule_3d(self, mol):
        """Dispatch to custom style or standard drawing."""
        mw = self
        
        if hasattr(mw, 'plugin_manager') and hasattr(mw.plugin_manager, 'custom_3d_styles'):
             if hasattr(self, 'current_3d_style') and self.current_3d_style in mw.plugin_manager.custom_3d_styles:
                 handler = mw.plugin_manager.custom_3d_styles[self.current_3d_style]['callback']
                 try:
                     handler(mw, mol)
                     return
                 except Exception as e:
                     logging.error(f"Error in custom 3d style '{self.current_3d_style}': {e}")
        
        self.draw_standard_3d_style(mol)

    def draw_standard_3d_style(self, mol, style_override=None):
        """3D 分子を描画し、軸アクターの参照をクリアする（軸の再制御は apply_3d_settings に任せる）"""
        
        current_style = style_override if style_override else self.current_3d_style

        # 測定選択をクリア（分子が変更されたため）
        if hasattr(self, 'measurement_mode'):
            self.clear_measurement_selection()
        
        # 色情報追跡のための辞書を初期化
        if not hasattr(self, '_3d_color_map'):
            self._3d_color_map = {}
        self._3d_color_map.clear()
        
        # 1. カメラ状態とクリア
        camera_state = self.plotter.camera.copy()

        # **残留防止のための強制削除**
        if self.axes_actor is not None:
            try:
                self.plotter.remove_actor(self.axes_actor)
            except Exception:
                pass 
            self.axes_actor = None

        self.plotter.clear()
            
        # 2. 背景色の設定
        self.plotter.set_background(self.settings.get('background_color', '#4f4f4f'))

        # 3. mol が None または原子数ゼロの場合は、背景と軸のみで終了
        if mol is None or mol.GetNumAtoms() == 0:
            self.atom_actor = None
            self.current_mol = None
            self.plotter.render()
            return
            
        # 4. ライティングの設定
        is_lighting_enabled = self.settings.get('lighting_enabled', True)

        if is_lighting_enabled:
            light = pv.Light(
                position=(1, 1, 2),
                light_type='cameralight',
                intensity=self.settings.get('light_intensity', 1.2)
            )
            self.plotter.add_light(light)
            
        # 5. 分子描画ロジック
        # Optionally kekulize aromatic systems for 3D visualization.
        mol_to_draw = mol
        if self.settings.get('display_kekule_3d', False):
            try:
                # Operate on a copy to avoid mutating the original molecule
                mol_to_draw = Chem.Mol(mol)
                Chem.Kekulize(mol_to_draw, clearAromaticFlags=True)
            except Exception as e:
                # Kekulize failed; keep original and warn user
                try:
                    self.statusBar().showMessage(f"Kekulize failed: {e}")
                except Exception:
                    pass
                mol_to_draw = mol

        # Use the original molecule's conformer (positions) to ensure coordinates
        # are preserved even when we create a kekulized copy for bond types.
        conf = mol.GetConformer()

        # Use the kekulized molecule's atom ordering for color/size decisions
        self.atom_positions_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_to_draw.GetNumAtoms())])

        # Use the possibly-kekulized molecule for symbol/bond types
        sym = [a.GetSymbol() for a in mol_to_draw.GetAtoms()]
        col = np.array([CPK_COLORS_PV.get(s, [0.5, 0.5, 0.5]) for s in sym])

        # Apply plugin color overrides
        if hasattr(self, '_plugin_color_overrides') and self._plugin_color_overrides:
            for atom_idx, hex_color in self._plugin_color_overrides.items():
                if 0 <= atom_idx < len(col):
                    try:
                        c = QColor(hex_color)
                        col[atom_idx] = [c.redF(), c.greenF(), c.blueF()]
                    except Exception:
                        pass

        # スタイルに応じて原子の半径を設定（設定から読み込み）
        if current_style == 'cpk':
            atom_scale = self.settings.get('cpk_atom_scale', 1.0)
            resolution = self.settings.get('cpk_resolution', 32)
            rad = np.array([pt.GetRvdw(pt.GetAtomicNumber(s)) * atom_scale for s in sym])
        elif current_style == 'wireframe':
            # Wireframeでは原子を描画しないので、この設定は実際には使用されない
            resolution = self.settings.get('wireframe_resolution', 6)
            rad = np.array([0.01 for s in sym])  # 極小値（使用されない）
        elif current_style == 'stick':
            atom_radius = self.settings.get('stick_bond_radius', 0.15)  # Use bond radius for atoms
            resolution = self.settings.get('stick_resolution', 16)
            rad = np.array([atom_radius for s in sym])
        else:  # ball_and_stick
            atom_scale = self.settings.get('ball_stick_atom_scale', 1.0)
            resolution = self.settings.get('ball_stick_resolution', 16)
            rad = np.array([VDW_RADII.get(s, 0.4) * atom_scale for s in sym])

        self.glyph_source = pv.PolyData(self.atom_positions_3d)
        self.glyph_source['colors'] = col
        self.glyph_source['radii'] = rad

        # メッシュプロパティを共通で定義
        mesh_props = dict(
            smooth_shading=True,
            specular=self.settings.get('specular', 0.2),
            specular_power=self.settings.get('specular_power', 20),
            lighting=is_lighting_enabled,
        )

        # Wireframeスタイルの場合は原子を描画しない
        if current_style != 'wireframe':
            # Stickモードで末端二重結合・三重結合の原子を分裂させるための処理
            if current_style == 'stick':
                # 末端原子（次数1）で多重結合を持つものを検出
                split_atoms = []  # (atom_idx, bond_order, offset_vecs)
                skip_atoms = set()  # スキップする原子のインデックス
                
                for i in range(mol_to_draw.GetNumAtoms()):
                    atom = mol_to_draw.GetAtomWithIdx(i)
                    if atom.GetDegree() == 1:  # 末端原子
                        bonds = atom.GetBonds()
                        if len(bonds) == 1:
                            bond = bonds[0]
                            bond_type = bond.GetBondType()
                            
                            if bond_type in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
                                # 多重結合を持つ末端原子を発見
                                # 結合のもう一方の原子を取得
                                other_idx = bond.GetBeginAtomIdx() if bond.GetEndAtomIdx() == i else bond.GetEndAtomIdx()
                                
                                # 結合ベクトルを計算
                                pos_i = np.array(conf.GetAtomPosition(i))
                                pos_other = np.array(conf.GetAtomPosition(other_idx))
                                bond_vec = pos_i - pos_other
                                bond_length = np.linalg.norm(bond_vec)
                                
                                if bond_length > 0:
                                    bond_unit = bond_vec / bond_length
                                    
                                    # 二重結合の場合は実際の描画と同じオフセット方向を使用
                                    if bond_type == Chem.BondType.DOUBLE:
                                        offset_dir1 = self._calculate_double_bond_offset(mol_to_draw, bond, conf)
                                    else:
                                        # 三重結合の場合は結合描画と同じロジック
                                        v_arb = np.array([0, 0, 1])
                                        if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0):
                                            v_arb = np.array([0, 1, 0])
                                        offset_dir1 = np.cross(bond_unit, v_arb)
                                        offset_dir1 /= np.linalg.norm(offset_dir1)
                                    
                                    # 二重/三重結合描画のオフセット値と半径を取得（結合描画と完全に一致させる）
                                    try:
                                        cyl_radius = self.settings.get('stick_bond_radius', 0.15)
                                        if bond_type == Chem.BondType.DOUBLE:
                                            radius_factor = self.settings.get('stick_double_bond_radius_factor', 0.60)
                                            offset_factor = self.settings.get('stick_double_bond_offset_factor', 1.5)
                                            # 二重結合：s_double / 2 を使用
                                            offset_distance = cyl_radius * offset_factor / 2
                                        else:  # TRIPLE
                                            radius_factor = self.settings.get('stick_triple_bond_radius_factor', 0.40)
                                            offset_factor = self.settings.get('stick_triple_bond_offset_factor', 1.0)
                                            # 三重結合：s_triple をそのまま使用（/ 2 なし）
                                            offset_distance = cyl_radius * offset_factor
                                        
                                        # 結合描画と同じ計算
                                        sphere_radius = cyl_radius * radius_factor
                                    except Exception:
                                        sphere_radius = 0.09  # デフォルト値
                                        offset_distance = 0.15  # デフォルト値
                                    
                                    if bond_type == Chem.BondType.DOUBLE:
                                        # 二重結合：2個に分裂
                                        offset_vecs = [
                                            offset_dir1 * offset_distance,
                                            -offset_dir1 * offset_distance
                                        ]
                                        split_atoms.append((i, 2, offset_vecs))
                                    else:  # TRIPLE
                                        # 三重結合：3個に分裂（中心 + 両側2つ）
                                        # 結合描画と同じ配置
                                        offset_vecs = [
                                            np.array([0, 0, 0]),  # 中心
                                            offset_dir1 * offset_distance,  # +side
                                            -offset_dir1 * offset_distance  # -side
                                        ]
                                        split_atoms.append((i, 3, offset_vecs))
                                    
                                    skip_atoms.add(i)
                
                # 分裂させる原子がある場合、新しい位置リストを作成
                if split_atoms:
                    new_positions = []
                    new_colors = []
                    new_radii = []
                    
                    # 通常の原子を追加（スキップリスト以外）
                    for i in range(len(self.atom_positions_3d)):
                        if i not in skip_atoms:
                            new_positions.append(self.atom_positions_3d[i])
                            new_colors.append(col[i])
                            new_radii.append(rad[i])
                    
                    # 分裂した原子を追加
                    # 上記で計算されたsphere_radiusを使用（結合描画のradius_factorを適用済み）
                    for atom_idx, bond_order, offset_vecs in split_atoms:
                        pos = self.atom_positions_3d[atom_idx]
                        # この原子の結合から半径を取得（上記ループで計算済み）
                        # 簡便のため、最後に計算されたsphere_radiusを使用
                        for offset_vec in offset_vecs:
                            new_positions.append(pos + offset_vec)
                            new_colors.append(col[atom_idx])
                            new_radii.append(sphere_radius)
                    
                    # PolyDataを新しい位置で作成
                    glyph_source = pv.PolyData(np.array(new_positions))
                    glyph_source['colors'] = np.array(new_colors)
                    glyph_source['radii'] = np.array(new_radii)
                else:
                    glyph_source = self.glyph_source
            else:
                glyph_source = self.glyph_source
            
            glyphs = glyph_source.glyph(scale='radii', geom=pv.Sphere(radius=1.0, theta_resolution=resolution, phi_resolution=resolution), orient=False)

            if is_lighting_enabled:
                self.atom_actor = self.plotter.add_mesh(glyphs, scalars='colors', rgb=True, **mesh_props)
            else:
                self.atom_actor = self.plotter.add_mesh(
                    glyphs, scalars='colors', rgb=True, 
                    style='surface', show_edges=True, edge_color='grey',
                    **mesh_props
                )
                self.atom_actor.GetProperty().SetEdgeOpacity(0.3)
            
            # 原子の色情報を記録
            for i, atom_color in enumerate(col):
                atom_rgb = [int(c * 255) for c in atom_color]
                self._3d_color_map[f'atom_{i}'] = atom_rgb


        # ボンドの描画（ball_and_stick、wireframe、stickで描画）
        if current_style in ['ball_and_stick', 'wireframe', 'stick']:
            # スタイルに応じてボンドの太さと解像度を設定（設定から読み込み）
            if current_style == 'wireframe':
                cyl_radius = self.settings.get('wireframe_bond_radius', 0.01)
                bond_resolution = self.settings.get('wireframe_resolution', 6)
            elif current_style == 'stick':
                cyl_radius = self.settings.get('stick_bond_radius', 0.15)
                bond_resolution = self.settings.get('stick_resolution', 16)
            else:  # ball_and_stick
                cyl_radius = self.settings.get('ball_stick_bond_radius', 0.1)
                bond_resolution = self.settings.get('ball_stick_resolution', 16)
            
            # Ball and Stick用の共通色
            bs_bond_rgb = [127, 127, 127]
            if current_style == 'ball_and_stick':
                try:
                    bs_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                    q = QColor(bs_hex)
                    bs_bond_rgb = [q.red(), q.green(), q.blue()]
                except Exception:
                    pass

            # バッチ処理用のリスト
            all_points = []
            all_lines = []
            all_radii = []
            all_colors = [] # Cell data (one per line segment)
            
            current_point_idx = 0
            bond_counter = 0

            for bond in mol_to_draw.GetBonds():
                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()
                sp = np.array(conf.GetAtomPosition(begin_atom_idx))
                ep = np.array(conf.GetAtomPosition(end_atom_idx))
                bt = bond.GetBondType()
                d = ep - sp
                h = np.linalg.norm(d)
                if h == 0: continue

                # ボンドの色
                begin_color = col[begin_atom_idx]
                end_color = col[end_atom_idx]
                begin_color_rgb = [int(c * 255) for c in begin_color]
                end_color_rgb = [int(c * 255) for c in end_color]

                # Check for plugin override
                bond_idx = bond.GetIdx()
                # Override handling: if set, force both ends and uniform color to this value
                if hasattr(self, '_plugin_bond_color_overrides') and bond_idx in self._plugin_bond_color_overrides:
                     try:
                         # Expecting hex string
                         hex_c = self._plugin_bond_color_overrides[bond_idx]
                         c_obj = QColor(hex_c)
                         ov_rgb = [c_obj.red(), c_obj.green(), c_obj.blue()]
                         begin_color_rgb = ov_rgb
                         end_color_rgb = ov_rgb
                         # Also override uniform color in case style uses it
                         # We need to use a local variable for this iteration instead of the global bs_bond_rgb
                         # But wait, bs_bond_rgb is defined outside loop.
                         # We can define local_bs_bond_rgb
                     except Exception:
                         pass
                
                # Determine effective uniform color for this bond
                local_bs_bond_rgb = begin_color_rgb if (hasattr(self, '_plugin_bond_color_overrides') and bond_idx in self._plugin_bond_color_overrides) else bs_bond_rgb


                # セグメント追加用ヘルパー関数
                def add_segment(p1, p2, radius, color_rgb):
                    nonlocal current_point_idx
                    all_points.append(p1)
                    all_points.append(p2)
                    all_lines.append([2, current_point_idx, current_point_idx + 1])
                    all_radii.append(radius)
                    all_radii.append(radius)
                    all_colors.append(color_rgb)
                    current_point_idx += 2

                QApplication.processEvents()

                # Get CPK bond color setting once for all bond types
                use_cpk_bond = self.settings.get('ball_stick_use_cpk_bond_color', False)
                # If overwritten, treat as if we want to show that color (effectively behave like CPK_Split but with same color, or Uniform).
                # To be robust, if overwritten, we can force "use_cpk_bond" logic but with our same colors?
                # Actually, if overridden, we probably want the whole bond to be that color. 
                
                is_overridden = hasattr(self, '_plugin_bond_color_overrides') and bond_idx in self._plugin_bond_color_overrides

                if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                    if current_style == 'ball_and_stick' and not use_cpk_bond and not is_overridden:
                        # 単一セグメント (Uniform color) - Default behavior
                        add_segment(sp, ep, cyl_radius, local_bs_bond_rgb)
                        self._3d_color_map[f'bond_{bond_counter}'] = local_bs_bond_rgb
                    else:
                        # 分割セグメント (CPK split colors OR Overridden uniform)
                        # If overridden, begin/end are same, so this produces a uniform looking bond split in middle
                        mid_point = (sp + ep) / 2
                        add_segment(sp, mid_point, cyl_radius, begin_color_rgb)
                        add_segment(mid_point, ep, cyl_radius, end_color_rgb)
                        self._3d_color_map[f'bond_{bond_counter}_start'] = begin_color_rgb
                        self._3d_color_map[f'bond_{bond_counter}_end'] = end_color_rgb

                else:
                    # 多重結合のパラメータ計算
                    v1 = d / h
                    # モデルごとの半径ファクターを適用
                    if current_style == 'ball_and_stick':
                        double_radius_factor = self.settings.get('ball_stick_double_bond_radius_factor', 0.8)
                        triple_radius_factor = self.settings.get('ball_stick_triple_bond_radius_factor', 0.75)
                    elif current_style == 'wireframe':
                        double_radius_factor = self.settings.get('wireframe_double_bond_radius_factor', 0.8)
                        triple_radius_factor = self.settings.get('wireframe_triple_bond_radius_factor', 0.75)
                    elif current_style == 'stick':
                        double_radius_factor = self.settings.get('stick_double_bond_radius_factor', 0.60)
                        triple_radius_factor = self.settings.get('stick_triple_bond_radius_factor', 0.40)
                    else:
                        double_radius_factor = 1.0
                        triple_radius_factor = 0.75
                    
                    # 設定からオフセットファクターを取得（モデルごと）
                    if current_style == 'ball_and_stick':
                        double_offset_factor = self.settings.get('ball_stick_double_bond_offset_factor', 2.0)
                        triple_offset_factor = self.settings.get('ball_stick_triple_bond_offset_factor', 2.0)
                    elif current_style == 'wireframe':
                        double_offset_factor = self.settings.get('wireframe_double_bond_offset_factor', 3.0)
                        triple_offset_factor = self.settings.get('wireframe_triple_bond_offset_factor', 3.0)
                    elif current_style == 'stick':
                        double_offset_factor = self.settings.get('stick_double_bond_offset_factor', 1.5)
                        triple_offset_factor = self.settings.get('stick_triple_bond_offset_factor', 1.0)
                    else:
                        double_offset_factor = 2.0
                        triple_offset_factor = 2.0

                    if bt == Chem.rdchem.BondType.DOUBLE:
                        r = cyl_radius * double_radius_factor
                        off_dir = self._calculate_double_bond_offset(mol_to_draw, bond, conf)
                        s_double = cyl_radius * double_offset_factor
                        
                        p1_start = sp + off_dir * (s_double / 2)
                        p1_end = ep + off_dir * (s_double / 2)
                        p2_start = sp - off_dir * (s_double / 2)
                        p2_end = ep - off_dir * (s_double / 2)

                        if current_style == 'ball_and_stick' and not use_cpk_bond and not is_overridden:
                            add_segment(p1_start, p1_end, r, local_bs_bond_rgb)
                            add_segment(p2_start, p2_end, r, local_bs_bond_rgb)
                            self._3d_color_map[f'bond_{bond_counter}_1'] = local_bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2'] = local_bs_bond_rgb
                        else:
                            mid1 = (p1_start + p1_end) / 2
                            mid2 = (p2_start + p2_end) / 2
                            add_segment(p1_start, mid1, r, begin_color_rgb)
                            add_segment(mid1, p1_end, r, end_color_rgb)
                            add_segment(p2_start, mid2, r, begin_color_rgb)
                            add_segment(mid2, p2_end, r, end_color_rgb)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_end'] = end_color_rgb

                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        r = cyl_radius * triple_radius_factor
                        v_arb = np.array([0, 0, 1])
                        if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                        off_dir = np.cross(v1, v_arb)
                        off_dir /= np.linalg.norm(off_dir)
                        s_triple = cyl_radius * triple_offset_factor

                        # Center
                        if current_style == 'ball_and_stick' and not use_cpk_bond and not is_overridden:
                            add_segment(sp, ep, r, local_bs_bond_rgb)
                            self._3d_color_map[f'bond_{bond_counter}_1'] = local_bs_bond_rgb
                        else:
                            mid = (sp + ep) / 2
                            add_segment(sp, mid, r, begin_color_rgb)
                            add_segment(mid, ep, r, end_color_rgb)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                        
                        # Sides
                        for sign in [1, -1]:
                            offset = off_dir * s_triple * sign
                            p_start = sp + offset
                            p_end = ep + offset
                            
                            if current_style == 'ball_and_stick' and not use_cpk_bond and not is_overridden:
                                add_segment(p_start, p_end, r, local_bs_bond_rgb)
                                suffix = '_2' if sign == 1 else '_3'
                                self._3d_color_map[f'bond_{bond_counter}{suffix}'] = local_bs_bond_rgb
                            else:
                                mid = (p_start + p_end) / 2
                                add_segment(p_start, mid, r, begin_color_rgb)
                                add_segment(mid, p_end, r, end_color_rgb)
                                suffix = '_2' if sign == 1 else '_3'
                                self._3d_color_map[f'bond_{bond_counter}{suffix}_start'] = begin_color_rgb
                                self._3d_color_map[f'bond_{bond_counter}{suffix}_end'] = end_color_rgb

                bond_counter += 1

            # ジオメトリの生成と描画
            if all_points:
                # Create PolyData
                bond_pd = pv.PolyData(np.array(all_points), lines=np.hstack(all_lines))
                # lines needs to be a flat array with padding indicating number of points per cell
                # all_lines is [[2, i, j], [2, k, l], ...], flatten it
                
                # Add data
                bond_pd.point_data['radii'] = np.array(all_radii)
                
                # Convert colors to 0-1 range for PyVista if needed, but add_mesh with rgb=True expects uint8 if using direct array?
                # Actually pyvista scalars usually prefer float 0-1 or uint8 0-255. 
                # Let's use uint8 0-255 and rgb=True.
                bond_pd.cell_data['colors'] = np.array(all_colors, dtype=np.uint8)
                
                # Tube filter
                # n_sides (resolution) corresponds to theta_resolution in Cylinder
                tube = bond_pd.tube(scalars='radii', absolute=True, radius_factor=1.0, n_sides=bond_resolution, capping=True)
                
                # Add to plotter
                self.plotter.add_mesh(tube, scalars='colors', rgb=True, **mesh_props)

        # Aromatic ring circles display
        if self.settings.get('display_aromatic_circles_3d', False):
            try:
                ring_info = mol_to_draw.GetRingInfo()
                aromatic_rings = []
                
                # Find aromatic rings
                for ring in ring_info.AtomRings():
                    # Check if all atoms in ring are aromatic
                    is_aromatic = all(mol_to_draw.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                    if is_aromatic:
                        aromatic_rings.append(ring)
                
                # Draw circles for aromatic rings
                for ring in aromatic_rings:
                    # Get atom positions
                    ring_positions = [self.atom_positions_3d[idx] for idx in ring]
                    ring_positions_np = np.array(ring_positions)
                    
                    # Calculate ring center
                    center = np.mean(ring_positions_np, axis=0)
                    
                    # Calculate ring normal using PCA or cross product
                    # Use first 3 atoms to get two vectors
                    if len(ring) >= 3:
                        v1 = ring_positions_np[1] - ring_positions_np[0]
                        v2 = ring_positions_np[2] - ring_positions_np[0]
                        normal = np.cross(v1, v2)
                        normal_length = np.linalg.norm(normal)
                        if normal_length > 0:
                            normal = normal / normal_length
                        else:
                            normal = np.array([0, 0, 1])
                    else:
                        normal = np.array([0, 0, 1])
                    
                    # Calculate ring radius (average distance from center)
                    distances = [np.linalg.norm(pos - center) for pos in ring_positions_np]
                    ring_radius = np.mean(distances) * 0.55  # Slightly smaller
                    
                    # Get bond radius from current style settings for torus thickness
                    if current_style == 'stick':
                        bond_radius = self.settings.get('stick_bond_radius', 0.15)
                    elif current_style == 'ball_and_stick':
                        bond_radius = self.settings.get('ball_stick_bond_radius', 0.1)
                    elif current_style == 'wireframe':
                        bond_radius = self.settings.get('wireframe_bond_radius', 0.01)
                    else:
                        bond_radius = 0.1  # Default
                    # Apply user-defined thickness factor (default 0.6)
                    thickness_factor = self.settings.get('aromatic_torus_thickness_factor', 0.6)
                    tube_radius = bond_radius * thickness_factor
                    theta = np.linspace(0, 2.2 * np.pi, 64)
                    circle_x = ring_radius * np.cos(theta)
                    circle_y = ring_radius * np.sin(theta)
                    circle_z = np.zeros_like(theta)
                    circle_points = np.c_[circle_x, circle_y, circle_z]
                    
                    # Create line from points
                    circle_line = pv.Spline(circle_points, n_points=64).tube(radius=tube_radius, n_sides=16)
                    
                    # Rotate torus to align with ring plane
                    # Default torus is in XY plane (normal = [0, 0, 1])
                    default_normal = np.array([0, 0, 1])
                    
                    # Calculate rotation axis and angle
                    if not np.allclose(normal, default_normal) and not np.allclose(normal, -default_normal):
                        axis = np.cross(default_normal, normal)
                        axis_length = np.linalg.norm(axis)
                        if axis_length > 0:
                            axis = axis / axis_length
                            angle = np.arccos(np.clip(np.dot(default_normal, normal), -1.0, 1.0))
                            angle_deg = np.degrees(angle)
                            
                            # Rotate torus
                            circle_line = circle_line.rotate_vector(axis, angle_deg, point=[0, 0, 0])
                    
                    # Translate to ring center
                    circle_line = circle_line.translate(center)
                    
                    # Get torus color from bond color settings
                    # Calculate most common atom type in ring for CPK color
                    from collections import Counter
                    atom_symbols = [mol_to_draw.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
                    most_common_symbol = Counter(atom_symbols).most_common(1)[0][0] if atom_symbols else None
                    
                    if current_style == 'ball_and_stick':
                        # Check if using CPK bond colors
                        use_cpk = self.settings.get('ball_stick_use_cpk_bond_color', False)
                        if use_cpk:
                            # Use CPK color of most common atom type in ring
                            if most_common_symbol:
                                cpk_color = CPK_COLORS_PV.get(most_common_symbol, [0.5, 0.5, 0.5])
                                torus_color = cpk_color
                            else:
                                torus_color = [0.5, 0.5, 0.5]
                        else:
                            # Use Ball & Stick bond color setting
                            bond_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                            q = QColor(bond_hex)
                            torus_color = [q.red() / 255.0, q.green() / 255.0, q.blue() / 255.0]
                    else:
                        # For Wireframe and Stick, use CPK color of most common atom
                        if most_common_symbol:
                            cpk_color = CPK_COLORS_PV.get(most_common_symbol, [0.5, 0.5, 0.5])
                            torus_color = cpk_color
                        else:
                            torus_color = [0.5, 0.5, 0.5]
                    
                    self.plotter.add_mesh(circle_line, color=torus_color, **mesh_props)
                    
            except Exception as e:
                logging.error(f"Error rendering aromatic circles: {e}")

        if getattr(self, 'show_chiral_labels', False):
            try:
                # 3D座標からキラル中心を計算
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if chiral_centers:
                    pts, labels = [], []
                    z_off = 0
                    for idx, lbl in chiral_centers:
                        coord = self.atom_positions_3d[idx].copy(); coord[2] += z_off
                        pts.append(coord); labels.append(lbl if lbl is not None else '?')
                    try: self.plotter.remove_actor('chiral_labels')
                    except Exception: pass
                    self.plotter.add_point_labels(np.array(pts), labels, font_size=20, point_size=0, text_color='blue', name='chiral_labels', always_visible=True, tolerance=0.01, show_points=False)
            except Exception as e: self.statusBar().showMessage(f"3D chiral label drawing error: {e}")

        # E/Zラベルも表示
        if getattr(self, 'show_chiral_labels', False):
            try:
                # If we drew a kekulized molecule use it for E/Z detection so
                # E/Z labels reflect Kekulé rendering; pass mol_to_draw as the
                # molecule to scan for bond stereochemistry.
                self.show_ez_labels_3d(mol)
            except Exception as e: 
                self.statusBar().showMessage(f"3D E/Z label drawing error: {e}")

        self.plotter.camera = camera_state

        # Ensure the underlying VTK camera's parallel/projection flag matches
        # the saved application setting. draw_molecule_3d restores a PyVista
        # camera object which may not propagate the ParallelProjection flag
        # to the VTK renderer camera; enforce it here to guarantee the
        # projection mode selected in settings actually takes effect.
        try:
            proj_mode = self.settings.get('projection_mode', 'Perspective')
            if hasattr(self.plotter, 'renderer') and hasattr(self.plotter.renderer, 'GetActiveCamera'):
                vcam = self.plotter.renderer.GetActiveCamera()
                if vcam:
                    if proj_mode == 'Orthographic':
                        vcam.SetParallelProjection(True)
                    else:
                        vcam.SetParallelProjection(False)
                    try:
                        # Force a render so the change is visible immediately
                        self.plotter.render()
                    except Exception:
                        pass
        except Exception:
            pass
        
        # AtomIDまたは他の原子情報が表示されている場合は再表示
        if hasattr(self, 'atom_info_display_mode') and self.atom_info_display_mode is not None:
            self.show_all_atom_info()
        
        # メニューテキストと状態を現在の分子の種類に応じて更新
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()



    def _calculate_double_bond_offset(self, mol, bond, conf):
        """
        二重結合のオフセット方向を計算する。
        結合している原子の他の結合を考慮して、平面的になるようにする。
        """
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        
        begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
        end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
        
        bond_vec = end_pos - begin_pos
        bond_length = np.linalg.norm(bond_vec)
        if bond_length == 0:
            # フォールバック: Z軸基準
            return np.array([0, 0, 1])
        
        bond_unit = bond_vec / bond_length
        
        # 両端の原子の隣接原子を調べる
        begin_neighbors = []
        end_neighbors = []
        
        for neighbor in begin_atom.GetNeighbors():
            if neighbor.GetIdx() != bond.GetEndAtomIdx():
                neighbor_pos = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                begin_neighbors.append(neighbor_pos)
        
        for neighbor in end_atom.GetNeighbors():
            if neighbor.GetIdx() != bond.GetBeginAtomIdx():
                neighbor_pos = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                end_neighbors.append(neighbor_pos)
        
        # 平面の法線ベクトルを計算
        normal_candidates = []
        
        # 開始原子の隣接原子から平面を推定
        if len(begin_neighbors) >= 1:
            for neighbor_pos in begin_neighbors:
                vec_to_neighbor = neighbor_pos - begin_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # bond_vec と neighbor_vec の外積が平面の法線
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)
        
        # 終了原子の隣接原子から平面を推定
        if len(end_neighbors) >= 1:
            for neighbor_pos in end_neighbors:
                vec_to_neighbor = neighbor_pos - end_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # bond_vec と neighbor_vec の外積が平面の法線
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)
        
        # 複数の法線ベクトルがある場合は平均を取る
        if normal_candidates:
            # 方向を統一するため、最初のベクトルとの内積が正になるように調整
            reference_normal = normal_candidates[0]
            aligned_normals = []
            
            for normal in normal_candidates:
                if np.dot(normal, reference_normal) < 0:
                    normal = -normal
                aligned_normals.append(normal)
            
            avg_normal = np.mean(aligned_normals, axis=0)
            norm_length = np.linalg.norm(avg_normal)
            if norm_length > 1e-6:
                avg_normal /= norm_length
                
                # 法線ベクトルと結合ベクトルに垂直な方向を二重結合のオフセット方向とする
                offset_dir = np.cross(bond_unit, avg_normal)
                offset_length = np.linalg.norm(offset_dir)
                if offset_length > 1e-6:
                    return offset_dir / offset_length
        
        # フォールバック: 結合ベクトルに垂直な任意の方向
        v_arb = np.array([0, 0, 1])
        if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0):
            v_arb = np.array([0, 1, 0])
        
        off_dir = np.cross(bond_unit, v_arb)
        off_dir /= np.linalg.norm(off_dir)
        return off_dir



    def show_ez_labels_3d(self, mol):
        """3DビューでE/Zラベルを表示する（RDKitのステレオ化学判定を使用）"""
        if not mol:
            return
        
        try:
            # 既存のE/Zラベルを削除
            self.plotter.remove_actor('ez_labels')
        except Exception:
            pass
        
        pts, labels = [], []
        
        # 3D座標が存在するかチェック
        if mol.GetNumConformers() == 0:
            return
            
        conf = mol.GetConformer()
        
        # 二重結合でRDKitが判定したE/Z立体化学を表示
        
        try:
            # 3D座標からステレオ化学を再計算 (molに対して行う)
            # これにより、2Dでの描画状態に関わらず、現在の3D座標に基づいたE/Z判定が行われる
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        except Exception:
            pass

        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                new_stereo = bond.GetStereo()
                
                if new_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                    # 結合の中心座標を計算
                    begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
                    end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
                    center_pos = (begin_pos + end_pos) / 2
                    
                    # ラベルの決定
                    label = 'E' if new_stereo == Chem.BondStereo.STEREOE else 'Z'
                    
                    # 2Dとの不一致チェック
                    # main_window_compute.py で保存された2D由来の立体化学プロパティを取得
                    try:
                        old_stereo = bond.GetIntProp("_original_2d_stereo")
                    except KeyError:
                        old_stereo = Chem.BondStereo.STEREONONE

                    # 2D側でもE/Zが指定されていて、かつ3Dと異なる場合は「?」にする
                    if old_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                        if old_stereo != new_stereo:
                            label = '?'
                            
                    pts.append(center_pos)
                    labels.append(label)
        
        if pts and labels:
            self.plotter.add_point_labels(
                np.array(pts), 
                labels, 
                font_size=18,
                point_size=0,
                text_color='darkgreen',  # 暗い緑色
                name='ez_labels',
                always_visible=True,
                tolerance=0.01,
                show_points=False
            )




    def toggle_chiral_labels_display(self, checked):
        """Viewメニューのアクションに応じてキラルラベル表示を切り替える"""
        self.show_chiral_labels = checked
        
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol) 
        
        if checked:
            self.statusBar().showMessage("Chiral labels: will be (re)computed after Convert→3D.")
        else:
            self.statusBar().showMessage("Chiral labels disabled.")




    def update_chiral_labels(self):
        """分子のキラル中心を計算し、2Dビューの原子アイテムにR/Sラベルを設定/解除する
        ※ 可能なら 3D（self.current_mol）を優先して計算し、なければ 2D から作った RDKit 分子を使う。
        """
        # まず全てのアイテムからラベルをクリア
        for atom_data in self.data.atoms.values():
            if atom_data.get('item'):
                atom_data['item'].chiral_label = None

        if not self.show_chiral_labels:
            self.scene.update()
            return

        # 3D の RDKit Mol（コンフォマーを持つもの）を使う
        mol_for_chirality = None
        if getattr(self, 'current_mol', None) is not None:
            mol_for_chirality = self.current_mol
        else:
            return

        if mol_for_chirality is None or mol_for_chirality.GetNumAtoms() == 0:
            self.scene.update()
            return

        try:
            # --- 重要：3D コンフォマーがあるなら、それを使って原子のキラルタグを割り当てる ---
            if mol_for_chirality.GetNumConformers() > 0:
                # confId=0（最初のコンフォマー）を指定して、原子のキラリティータグを3D座標由来で設定
                try:
                    Chem.AssignAtomChiralTagsFromStructure(mol_for_chirality, confId=0)
                except Exception:
                    # 古い RDKit では関数が無い場合があるので（念のため保護）
                    pass

            # RDKit の通常の stereochemistry 割当（念のため）
            #Chem.AssignStereochemistry(mol_for_chirality, cleanIt=True, force=True, flagPossibleStereoCenters=True)

            # キラル中心の取得（(idx, 'R'/'S'/'?') のリスト）
            chiral_centers = Chem.FindMolChiralCenters(mol_for_chirality, includeUnassigned=True)

            # RDKit atom index -> エディタ側 atom_id へのマッピング
            rdkit_idx_to_my_id = {}
            for atom in mol_for_chirality.GetAtoms():
                if atom.HasProp("_original_atom_id"):
                    rdkit_idx_to_my_id[atom.GetIdx()] = atom.GetIntProp("_original_atom_id")

            # 見つかったキラル中心を対応する AtomItem に設定
            for idx, label in chiral_centers:
                if idx in rdkit_idx_to_my_id:
                    atom_id = rdkit_idx_to_my_id[idx]
                    if atom_id in self.data.atoms and self.data.atoms[atom_id].get('item'):
                        # 'R' / 'S' / '?'
                        self.data.atoms[atom_id]['item'].chiral_label = label

        except Exception as e:
            self.statusBar().showMessage(f"Update chiral labels error: {e}")

        # 最後に 2D シーンを再描画
        self.scene.update()



    def toggle_atom_info_display(self, mode):
        """原子情報表示モードを切り替える"""
        # 現在の表示をクリア
        self.clear_all_atom_info_labels()
        
        # 同じモードが選択された場合はOFFにする
        if self.atom_info_display_mode == mode:
            self.atom_info_display_mode = None
            # 全てのアクションのチェックを外す
            self.show_atom_id_action.setChecked(False)
            self.show_rdkit_id_action.setChecked(False)
            self.show_atom_coords_action.setChecked(False)
            self.show_atom_symbol_action.setChecked(False)
            self.statusBar().showMessage("Atom info display disabled.")
        else:
            # 新しいモードを設定
            self.atom_info_display_mode = mode
            # 該当するアクションのみチェック
            self.show_atom_id_action.setChecked(mode == 'id')
            self.show_rdkit_id_action.setChecked(mode == 'rdkit_id')
            self.show_atom_coords_action.setChecked(mode == 'coords')
            self.show_atom_symbol_action.setChecked(mode == 'symbol')
            
            mode_names = {'id': 'Atom ID', 'rdkit_id': 'RDKit Index', 'coords': 'Coordinates', 'symbol': 'Element Symbol'}
            self.statusBar().showMessage(f"Displaying: {mode_names[mode]}")
            
            # すべての原子に情報を表示
            self.show_all_atom_info()



    def is_xyz_derived_molecule(self):
        """現在の分子がXYZファイル由来かどうかを判定"""
        if not self.current_mol:
            return False
        try:
            # 最初の原子がxyz_unique_idプロパティを持っているかチェック
            if self.current_mol.GetNumAtoms() > 0:
                return self.current_mol.GetAtomWithIdx(0).HasProp("xyz_unique_id")
        except Exception:
            pass
        return False



    def has_original_atom_ids(self):
        """現在の分子がOriginal Atom IDsを持っているかどうかを判定"""
        if not self.current_mol:
            return False
        try:
            # いずれかの原子が_original_atom_idプロパティを持っているかチェック
            for atom_idx in range(self.current_mol.GetNumAtoms()):
                atom = self.current_mol.GetAtomWithIdx(atom_idx)
                if atom.HasProp("_original_atom_id"):
                    return True
        except Exception:
            pass
        return False



    def update_atom_id_menu_text(self):
        """原子IDメニューのテキストを現在の分子の種類に応じて更新"""
        if hasattr(self, 'show_atom_id_action'):
            if self.is_xyz_derived_molecule():
                self.show_atom_id_action.setText("Show XYZ Unique ID")
            else:
                self.show_atom_id_action.setText("Show Original ID / Index")



    def update_atom_id_menu_state(self):
        """原子IDメニューの有効/無効状態を更新"""
        if hasattr(self, 'show_atom_id_action'):
            has_original_ids = self.has_original_atom_ids()
            has_xyz_ids = self.is_xyz_derived_molecule()
            
            # Original IDまたはXYZ IDがある場合のみ有効化
            self.show_atom_id_action.setEnabled(has_original_ids or has_xyz_ids)
            
            # 現在選択されているモードが無効化される場合は解除
            if not (has_original_ids or has_xyz_ids) and self.atom_info_display_mode == 'id':
                self.atom_info_display_mode = None
                self.show_atom_id_action.setChecked(False)
                self.clear_all_atom_info_labels()




    def show_all_atom_info(self):
        """すべての原子に情報を表示"""
        if self.atom_info_display_mode is None or not hasattr(self, 'atom_positions_3d') or self.atom_positions_3d is None:
            return
        
        # 既存のラベルをクリア
        self.clear_all_atom_info_labels()

        # ラベルを表示するためにタイプ別に分けてリストを作る
        rdkit_positions = []
        rdkit_texts = []
        id_positions = []
        id_texts = []
        xyz_positions = []
        xyz_texts = []
        other_positions = []
        other_texts = []

        for atom_idx, pos in enumerate(self.atom_positions_3d):
            # default: skip if no display mode
            if self.atom_info_display_mode is None:
                continue

            if self.atom_info_display_mode == 'id':
                # Original IDがある場合は優先表示、なければXYZのユニークID、最後にRDKitインデックス
                try:
                    if self.current_mol:
                        atom = self.current_mol.GetAtomWithIdx(atom_idx)
                        if atom.HasProp("_original_atom_id"):
                            original_id = atom.GetIntProp("_original_atom_id")
                            # プレフィックスを削除して数値だけ表示
                            id_positions.append(pos)
                            id_texts.append(str(original_id))
                        elif atom.HasProp("xyz_unique_id"):
                            unique_id = atom.GetIntProp("xyz_unique_id")
                            xyz_positions.append(pos)
                            xyz_texts.append(str(unique_id))
                        else:
                            rdkit_positions.append(pos)
                            rdkit_texts.append(str(atom_idx))
                    else:
                        rdkit_positions.append(pos)
                        rdkit_texts.append(str(atom_idx))
                except Exception:
                    rdkit_positions.append(pos)
                    rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == 'rdkit_id':
                rdkit_positions.append(pos)
                rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == 'coords':
                other_positions.append(pos)
                other_texts.append(f"({pos[0]:.2f},{pos[1]:.2f},{pos[2]:.2f})")

            elif self.atom_info_display_mode == 'symbol':
                if self.current_mol:
                    symbol = self.current_mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    other_positions.append(pos)
                    other_texts.append(symbol)
                else:
                    other_positions.append(pos)
                    other_texts.append("?")

            else:
                continue

        # 色の定義（暗めの青/緑/赤）
        rdkit_color = '#003366'   # 暗めの青
        id_color = '#006400'      # 暗めの緑
        xyz_color = '#8B0000'     # 暗めの赤
        other_color = 'black'

        # それぞれのグループごとにラベルを追加し、参照をリストで保持する
        self.current_atom_info_labels = []
        try:
            if rdkit_positions:
                a = self.plotter.add_point_labels(
                    np.array(rdkit_positions), rdkit_texts,
                    point_size=12, font_size=18, text_color=rdkit_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_rdkit'
                )
                self.current_atom_info_labels.append(a)

            if id_positions:
                a = self.plotter.add_point_labels(
                    np.array(id_positions), id_texts,
                    point_size=12, font_size=18, text_color=id_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_id'
                )
                self.current_atom_info_labels.append(a)

            if xyz_positions:
                a = self.plotter.add_point_labels(
                    np.array(xyz_positions), xyz_texts,
                    point_size=12, font_size=18, text_color=xyz_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_xyz'
                )
                self.current_atom_info_labels.append(a)

            if other_positions:
                a = self.plotter.add_point_labels(
                    np.array(other_positions), other_texts,
                    point_size=12, font_size=18, text_color=other_color,
                    always_visible=True, tolerance=0.01, show_points=False,
                    name='atom_labels_other'
                )
                self.current_atom_info_labels.append(a)
        except Exception as e:
            print(f"Error adding atom info labels: {e}")

        # 右上に凡例を表示（既存の凡例は消す）
        try:
            # 古い凡例削除
            if hasattr(self, 'atom_label_legend_names') and self.atom_label_legend_names:
                for nm in self.atom_label_legend_names:
                    try:
                        self.plotter.remove_actor(nm)
                    except Exception:
                        pass
            self.atom_label_legend_names = []

            # 凡例テキストを右上に縦並びで追加（背景なし、太字のみ）
            legend_entries = []
            if rdkit_positions:
                legend_entries.append(('RDKit', rdkit_color, 'legend_rdkit'))
            if id_positions:
                legend_entries.append(('ID', id_color, 'legend_id'))
            if xyz_positions:
                legend_entries.append(('XYZ', xyz_color, 'legend_xyz'))
            # Do not show 'Other' in the legend per UI requirement
            # (other_positions are still labeled in-scene but not listed in the legend)

            # 左下に凡例ラベルを追加（背景なし、太字のみ）
            # Increase spacing to avoid overlapping when short labels like 'RDKit' and 'ID' appear
            spacing = 30
            for i, (label_text, label_color, label_name) in enumerate(legend_entries):
                # 左下基準でy座標を上げる
                # Add a small horizontal offset for very short adjacent labels so they don't visually collide
                y = 0.0 + i * spacing
                x_offset = 0.0
                # If both RDKit and ID are present, nudge the second entry slightly to the right to avoid overlap
                try:
                    if label_text == 'ID' and any(e[0] == 'RDKit' for e in legend_entries):
                        x_offset = 0.06
                except Exception:
                    x_offset = 0.0
                try:
                    actor = self.plotter.add_text(
                        label_text,
                        position=(0.0 + x_offset, y),
                        font_size=12,
                        color=label_color,
                        name=label_name,
                        font='arial'
                    )
                    self.atom_label_legend_names.append(label_name)
                    # 太字のみ設定（背景は設定しない）
                    try:
                        if hasattr(actor, 'GetTextProperty'):
                            tp = actor.GetTextProperty()
                            try:
                                tp.SetBold(True)
                            except Exception:
                                pass
                    except Exception:
                        pass
                except Exception:
                    continue

        except Exception:
            pass



    def clear_all_atom_info_labels(self):
        """すべての原子情報ラベルをクリア"""
        # Remove label actors (may be a single actor, a list, or None)
        try:
            if hasattr(self, 'current_atom_info_labels') and self.current_atom_info_labels:
                if isinstance(self.current_atom_info_labels, (list, tuple)):
                    for a in list(self.current_atom_info_labels):
                        try:
                            self.plotter.remove_actor(a)
                        except Exception:
                            pass
                else:
                    try:
                        self.plotter.remove_actor(self.current_atom_info_labels)
                    except Exception:
                        pass
        except Exception:
            pass
        finally:
            self.current_atom_info_labels = None

        # Remove legend text actors if present
        try:
            if hasattr(self, 'atom_label_legend_names') and self.atom_label_legend_names:
                for nm in list(self.atom_label_legend_names):
                    try:
                        self.plotter.remove_actor(nm)
                    except Exception:
                        pass
        except Exception:
            pass
        finally:
            self.atom_label_legend_names = []



    def setup_3d_hover(self):
        """3Dビューでの表示を設定（常時表示に変更）"""
        if self.atom_info_display_mode is not None:
            self.show_all_atom_info()



    def zoom_in(self):
        """ ビューを 20% 拡大する """
        self.view_2d.scale(1.2, 1.2)



    def zoom_out(self):
        """ ビューを 20% 縮小する """
        self.view_2d.scale(1/1.2, 1/1.2)
        


    def reset_zoom(self):
        """ ビューの拡大率をデフォルト (75%) にリセットする """
        transform = QTransform()
        transform.scale(0.75, 0.75)
        self.view_2d.setTransform(transform)



    def fit_to_view(self):
        """ シーン上のすべてのアイテムがビューに収まるように調整する """
        if not self.scene.items():
            self.reset_zoom()
            return
            
        # 合計の表示矩形（目に見えるアイテムのみ）を計算
        visible_items_rect = QRectF()
        for item in self.scene.items():
            if item.isVisible() and not isinstance(item, TemplatePreviewItem):
                if visible_items_rect.isEmpty():
                    visible_items_rect = item.sceneBoundingRect()
                else:
                    visible_items_rect = visible_items_rect.united(item.sceneBoundingRect())

        if visible_items_rect.isEmpty():
            self.reset_zoom()
            return

        # 少し余白を持たせる（パディング）
        padding_factor = 1.10  # 10% の余裕
        cx = visible_items_rect.center().x()
        cy = visible_items_rect.center().y()
        w = visible_items_rect.width() * padding_factor
        h = visible_items_rect.height() * padding_factor
        padded = QRectF(cx - w / 2.0, cy - h / 2.0, w, h)

        # フィット時にマウス位置に依存するアンカーが原因でジャンプすることがあるため
        # 一時的にトランスフォームアンカーをビュー中心にしてから fitInView を呼ぶ
        try:
            old_ta = self.view_2d.transformationAnchor()
            old_ra = self.view_2d.resizeAnchor()
        except Exception:
            old_ta = old_ra = None

        try:
            self.view_2d.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
            self.view_2d.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
            self.view_2d.fitInView(padded, Qt.AspectRatioMode.KeepAspectRatio)
        finally:
            # 元のアンカーを復元
            try:
                if old_ta is not None:
                    self.view_2d.setTransformationAnchor(old_ta)
                if old_ra is not None:
                    self.view_2d.setResizeAnchor(old_ra)
            except Exception:
                pass



    def update_cpk_colors_from_settings(self):
        """Update global CPK_COLORS and CPK_COLORS_PV from saved settings overrides.

        This modifies the in-memory CPK_COLORS mapping (not persisted until settings are saved).
        Only keys present in self.settings['cpk_colors'] are changed; other elements keep the defaults.
        """
        try:
            # Overridden CPK settings are stored in self.settings['cpk_colors'].
            # To ensure that 2D modules (e.g., atom_item.py) which imported the
            # `CPK_COLORS` mapping from `modules.constants` at import time see
            # updates, mutate the mapping in-place on the constants module
            # instead of rebinding a new local variable here.
            overrides = self.settings.get('cpk_colors', {}) or {}

            # Import the constants module so we can update mappings directly
            try:
                from . import constants as constants_mod
            except Exception:
                import modules.constants as constants_mod

            # Reset constants.CPK_COLORS to defaults but keep the same dict
            constants_mod.CPK_COLORS.clear()
            for k, v in DEFAULT_CPK_COLORS.items():
                constants_mod.CPK_COLORS[k] = QColor(v) if not isinstance(v, QColor) else v

            # Apply overrides from settings
            for k, hexv in overrides.items():
                if isinstance(hexv, str) and hexv:
                    constants_mod.CPK_COLORS[k] = QColor(hexv)

            # Rebuild the PV representation in-place too
            constants_mod.CPK_COLORS_PV.clear()
            for k, c in constants_mod.CPK_COLORS.items():
                constants_mod.CPK_COLORS_PV[k] = [c.redF(), c.greenF(), c.blueF()]
        except Exception as e:
            print(f"Failed to update CPK colors from settings: {e}")




    def apply_3d_settings(self, redraw=True):
        # Projection mode
        proj_mode = self.settings.get('projection_mode', 'Perspective')
        if hasattr(self.plotter, 'renderer') and hasattr(self.plotter.renderer, 'GetActiveCamera'):
            cam = self.plotter.renderer.GetActiveCamera()
            if cam:
                if proj_mode == 'Orthographic':
                    cam.SetParallelProjection(True)
                else:
                    cam.SetParallelProjection(False)
        """3Dビューの視覚設定を適用する"""
        if not hasattr(self, 'plotter'):
            return  
        
        # レンダラーのレイヤー設定を有効化（テキストオーバーレイ用）
        renderer = self.plotter.renderer
        if renderer and hasattr(renderer, 'SetNumberOfLayers'):
            try:
                renderer.SetNumberOfLayers(2)  # レイヤー0:3Dオブジェクト、レイヤー1:2Dオーバーレイ
            except Exception:
                pass  # PyVistaのバージョンによってはサポートされていない場合がある  

        # --- 3D軸ウィジェットの設定 ---
        show_axes = self.settings.get('show_3d_axes', True) 

        # ウィジェットがまだ作成されていない場合は作成する
        if self.axes_widget is None and hasattr(self.plotter, 'interactor'):
            axes = vtk.vtkAxesActor()
            self.axes_widget = vtk.vtkOrientationMarkerWidget()
            self.axes_widget.SetOrientationMarker(axes)
            self.axes_widget.SetInteractor(self.plotter.interactor)
            # 左下隅に設定 (幅・高さ20%)
            self.axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)

        # 設定に応じてウィジェットを有効化/無効化
        if self.axes_widget:
            if show_axes:
                self.axes_widget.On()
                self.axes_widget.SetInteractive(False)  
            else:
                self.axes_widget.Off()  

        if redraw:
            self.draw_molecule_3d(self.current_mol)

        # 設定変更時にカメラ位置をリセットしない（初回のみリセット）
        if not getattr(self, '_camera_initialized', False):
            try:
                self.plotter.reset_camera()
            except Exception:
                pass
            self._camera_initialized = True
        
        # 強制的にプロッターを更新
        try:
            self.plotter.render()
            if hasattr(self.plotter, 'update'):
                self.plotter.update()
        except Exception:
            pass




    def update_bond_color_override(self, bond_idx, hex_color):
        """Plugin API helper to override bond color."""
        if not hasattr(self, '_plugin_bond_color_overrides'):
            self._plugin_bond_color_overrides = {}
            
        if hex_color is None:
            if bond_idx in self._plugin_bond_color_overrides:
                del self._plugin_bond_color_overrides[bond_idx]
        else:
            self._plugin_bond_color_overrides[bond_idx] = hex_color

        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)

    def update_atom_color_override(self, atom_index, color_hex):
        """Plugin helper to update specific atom color override."""
        if not hasattr(self, '_plugin_color_overrides'):
            self._plugin_color_overrides = {}
        
        if color_hex is None:
            if atom_index in self._plugin_color_overrides:
                del self._plugin_color_overrides[atom_index]
        else:
            self._plugin_color_overrides[atom_index] = color_hex
            
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)




