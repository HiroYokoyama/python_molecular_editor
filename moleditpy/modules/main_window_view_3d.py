#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_window_view_3d.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowView3d
"""


import numpy as np
import vtk


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
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
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

    def __init__(self, main_window):
        """ クラスの初期化 """
        self.mw = main_window


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
        """3D 分子を描画し、軸アクターの参照をクリアする（軸の再制御は apply_3d_settings に任せる）"""
        
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

        # スタイルに応じて原子の半径を設定（設定から読み込み）
        if self.current_3d_style == 'cpk':
            atom_scale = self.settings.get('cpk_atom_scale', 1.0)
            resolution = self.settings.get('cpk_resolution', 32)
            rad = np.array([pt.GetRvdw(pt.GetAtomicNumber(s)) * atom_scale for s in sym])
        elif self.current_3d_style == 'wireframe':
            # Wireframeでは原子を描画しないので、この設定は実際には使用されない
            resolution = self.settings.get('wireframe_resolution', 6)
            rad = np.array([0.01 for s in sym])  # 極小値（使用されない）
        elif self.current_3d_style == 'stick':
            atom_radius = self.settings.get('stick_atom_radius', 0.15)
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
        if self.current_3d_style != 'wireframe':
            glyphs = self.glyph_source.glyph(scale='radii', geom=pv.Sphere(radius=1.0, theta_resolution=resolution, phi_resolution=resolution), orient=False)

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
        if self.current_3d_style in ['ball_and_stick', 'wireframe', 'stick']:
            # スタイルに応じてボンドの太さと解像度を設定（設定から読み込み）
            if self.current_3d_style == 'wireframe':
                cyl_radius = self.settings.get('wireframe_bond_radius', 0.01)
                bond_resolution = self.settings.get('wireframe_resolution', 6)
            elif self.current_3d_style == 'stick':
                cyl_radius = self.settings.get('stick_bond_radius', 0.15)
                bond_resolution = self.settings.get('stick_resolution', 16)
            else:  # ball_and_stick
                cyl_radius = self.settings.get('ball_stick_bond_radius', 0.1)
                bond_resolution = self.settings.get('ball_stick_resolution', 16)
            
            bond_counter = 0  # 結合の個別識別用
            
            # Ball and Stick用のシリンダーリストを準備（高速化のため）
            if self.current_3d_style == 'ball_and_stick':
                bond_cylinders = []
                # Compute the configured grey/uniform bond color for Ball & Stick
                try:
                    bs_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                    q = QColor(bs_hex)
                    bs_bond_rgb = [q.red(), q.green(), q.blue()]
                except Exception:
                    bs_bond_rgb = [127, 127, 127]
            
            for bond in mol_to_draw.GetBonds():
                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()
                sp = np.array(conf.GetAtomPosition(begin_atom_idx))
                ep = np.array(conf.GetAtomPosition(end_atom_idx))
                bt = bond.GetBondType()
                c = (sp + ep) / 2
                d = ep - sp
                h = np.linalg.norm(d)
                if h == 0: continue

                # ボンドの色を原子の色から決定（各半分で異なる色）
                begin_color = col[begin_atom_idx]
                end_color = col[end_atom_idx]
                
                # 結合の色情報を記録
                begin_color_rgb = [int(c * 255) for c in begin_color]
                end_color_rgb = [int(c * 255) for c in end_color]

                # UI応答性維持のためイベント処理
                QApplication.processEvents()
                if bt == Chem.rdchem.BondType.SINGLE or bt == Chem.rdchem.BondType.AROMATIC:
                    if self.current_3d_style == 'ball_and_stick':
                        # Ball and stickは全結合をまとめて処理（高速化）
                        cyl = pv.Cylinder(center=c, direction=d, radius=cyl_radius, height=h, resolution=bond_resolution)
                        bond_cylinders.append(cyl)
                        self._3d_color_map[f'bond_{bond_counter}'] = bs_bond_rgb  # グレー (configurable)
                    else:
                        # その他（stick, wireframe）は中央で色が変わる2つの円柱
                        mid_point = (sp + ep) / 2
                        
                        # 前半（開始原子の色）
                        cyl1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                        actor1 = self.plotter.add_mesh(cyl1, color=begin_color, **mesh_props)
                        self._3d_color_map[f'bond_{bond_counter}_start'] = begin_color_rgb
                        
                        # 後半（終了原子の色）
                        cyl2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=cyl_radius, height=h/2, resolution=bond_resolution)
                        actor2 = self.plotter.add_mesh(cyl2, color=end_color, **mesh_props)
                        self._3d_color_map[f'bond_{bond_counter}_end'] = end_color_rgb
                else:
                    v1 = d / h
                    # モデルごとの半径ファクターを適用
                    if self.current_3d_style == 'ball_and_stick':
                        double_radius_factor = self.settings.get('ball_stick_double_bond_radius_factor', 0.8)
                        triple_radius_factor = self.settings.get('ball_stick_triple_bond_radius_factor', 0.75)
                    elif self.current_3d_style == 'wireframe':
                        double_radius_factor = self.settings.get('wireframe_double_bond_radius_factor', 1.0)
                        triple_radius_factor = self.settings.get('wireframe_triple_bond_radius_factor', 0.75)
                    elif self.current_3d_style == 'stick':
                        double_radius_factor = self.settings.get('stick_double_bond_radius_factor', 0.60)
                        triple_radius_factor = self.settings.get('stick_triple_bond_radius_factor', 0.40)
                    else:
                        double_radius_factor = 1.0
                        triple_radius_factor = 0.75
                    r = cyl_radius * 0.8  # fallback, will be overridden below
                    # 設定からオフセットファクターを取得（モデルごと）
                    if self.current_3d_style == 'ball_and_stick':
                        double_offset_factor = self.settings.get('ball_stick_double_bond_offset_factor', 2.0)
                        triple_offset_factor = self.settings.get('ball_stick_triple_bond_offset_factor', 2.0)
                    elif self.current_3d_style == 'wireframe':
                        double_offset_factor = self.settings.get('wireframe_double_bond_offset_factor', 3.0)
                        triple_offset_factor = self.settings.get('wireframe_triple_bond_offset_factor', 3.0)
                    elif self.current_3d_style == 'stick':
                        double_offset_factor = self.settings.get('stick_double_bond_offset_factor', 1.5)
                        triple_offset_factor = self.settings.get('stick_triple_bond_offset_factor', 1.0)
                    else:
                        double_offset_factor = 2.0
                        triple_offset_factor = 2.0
                    s = cyl_radius * 2.0  # デフォルト値

                    if bt == Chem.rdchem.BondType.DOUBLE:
                        r = cyl_radius * double_radius_factor
                        # 二重結合の場合、結合している原子の他の結合を考慮してオフセット方向を決定
                        off_dir = self._calculate_double_bond_offset(mol_to_draw, bond, conf)
                        # 設定から二重結合のオフセットファクターを適用
                        s_double = cyl_radius * double_offset_factor
                        c1, c2 = c + off_dir * (s_double / 2), c - off_dir * (s_double / 2)
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickは全結合をまとめて処理（高速化）
                            cyl1 = pv.Cylinder(center=c1, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c2, direction=d, radius=r, height=h, resolution=bond_resolution)
                            bond_cylinders.extend([cyl1, cyl2])
                            self._3d_color_map[f'bond_{bond_counter}_1'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2'] = bs_bond_rgb
                        else:
                            # その他（stick, wireframe）は中央で色が変わる
                            mid_point = (sp + ep) / 2
                            
                            # 第一の結合線（前半・後半）
                            cyl1_1 = pv.Cylinder(center=(sp + mid_point) / 2 + off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl1_2 = pv.Cylinder(center=(mid_point + ep) / 2 + off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl1_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                            
                            # 第二の結合線（前半・後半）
                            cyl2_1 = pv.Cylinder(center=(sp + mid_point) / 2 - off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl2_2 = pv.Cylinder(center=(mid_point + ep) / 2 - off_dir * (s_double / 2), direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl2_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl2_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_2_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_end'] = end_color_rgb
                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        r = cyl_radius * triple_radius_factor
                        # 三重結合
                        v_arb = np.array([0, 0, 1])
                        if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0): v_arb = np.array([0, 1, 0])
                        off_dir = np.cross(v1, v_arb)
                        off_dir /= np.linalg.norm(off_dir)
                        
                        # 設定から三重結合のオフセットファクターを適用
                        s_triple = cyl_radius * triple_offset_factor
                        
                        if self.current_3d_style == 'ball_and_stick':
                            # Ball and stickは全結合をまとめて処理（高速化）
                            cyl1 = pv.Cylinder(center=c, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl2 = pv.Cylinder(center=c + off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            cyl3 = pv.Cylinder(center=c - off_dir * s_triple, direction=d, radius=r, height=h, resolution=bond_resolution)
                            bond_cylinders.extend([cyl1, cyl2, cyl3])
                            self._3d_color_map[f'bond_{bond_counter}_1'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2'] = bs_bond_rgb
                            self._3d_color_map[f'bond_{bond_counter}_3'] = bs_bond_rgb
                        else:
                            # その他（stick, wireframe）は中央で色が変わる
                            mid_point = (sp + ep) / 2
                            
                            # 中央の結合線（前半・後半）
                            cyl1_1 = pv.Cylinder(center=(sp + mid_point) / 2, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl1_2 = pv.Cylinder(center=(mid_point + ep) / 2, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl1_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl1_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_1_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_1_end'] = end_color_rgb
                            
                            # 上側の結合線（前半・後半）
                            cyl2_1 = pv.Cylinder(center=(sp + mid_point) / 2 + off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl2_2 = pv.Cylinder(center=(mid_point + ep) / 2 + off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl2_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl2_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_2_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_2_end'] = end_color_rgb
                            
                            # 下側の結合線（前半・後半）
                            cyl3_1 = pv.Cylinder(center=(sp + mid_point) / 2 - off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            cyl3_2 = pv.Cylinder(center=(mid_point + ep) / 2 - off_dir * s_triple, direction=d, radius=r, height=h/2, resolution=bond_resolution)
                            self.plotter.add_mesh(cyl3_1, color=begin_color, **mesh_props)
                            self.plotter.add_mesh(cyl3_2, color=end_color, **mesh_props)
                            self._3d_color_map[f'bond_{bond_counter}_3_start'] = begin_color_rgb
                            self._3d_color_map[f'bond_{bond_counter}_3_end'] = end_color_rgb

                bond_counter += 1
            
            # Ball and Stick用：全結合をまとめて一括描画（高速化）
            if self.current_3d_style == 'ball_and_stick' and bond_cylinders:
                # 全シリンダーを結合してMultiBlockを作成
                combined_bonds = pv.MultiBlock(bond_cylinders)
                combined_mesh = combined_bonds.combine()
                
                # 一括でグレーで描画
                # Use the configured Ball & Stick bond color (hex) for the combined bonds
                try:
                    bs_hex = self.settings.get('ball_stick_bond_color', '#7F7F7F')
                    q = QColor(bs_hex)
                    # Use normalized RGB for pyvista (r,g,b) floats in [0,1]
                    bond_color = (q.redF(), q.greenF(), q.blueF())
                    bond_actor = self.plotter.add_mesh(combined_mesh, color=bond_color, **mesh_props)
                except Exception:
                    bond_actor = self.plotter.add_mesh(combined_mesh, color='grey', **mesh_props)
                
                # まとめて色情報を記録
                self._3d_color_map['bonds_combined'] = bs_bond_rgb

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
                self.show_ez_labels_3d(mol, scan_mol=mol_to_draw)
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



    def show_ez_labels_3d(self, mol, scan_mol=None):
        """3DビューでE/Zラベルを表示する（RDKitのステレオ化学判定を使用）"""
        if not mol:
            return
        
        try:
            # 既存のE/Zラベルを削除
            self.plotter.remove_actor('ez_labels')
        except:
            pass
        
        pts, labels = [], []
        
        # 3D座標が存在するかチェック
        if mol.GetNumConformers() == 0:
            return
            
        conf = mol.GetConformer()
        
        # RDKitに3D座標からステレオ化学を計算させる
        try:
            # 3D座標からステレオ化学を再計算
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        except:
            pass
        
        # 二重結合でRDKitが判定したE/Z立体化学を表示
        # `scan_mol` is used for stereochemistry detection (bond types); default
        # to the provided molecule if not supplied.
        if scan_mol is None:
            scan_mol = mol

        for bond in scan_mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                    # 結合の中心座標を計算
                    # Use positions from the original molecule's conformer; `bond` may
                    # come from `scan_mol` which can be kekulized but position indices
                    # correspond to the original `mol`.
                    begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
                    end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
                    center_pos = (begin_pos + end_pos) / 2
                    
                    # RDKitの判定結果を使用
                    label = 'E' if stereo == Chem.BondStereo.STEREOE else 'Z'
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
                    except:
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
                        except:
                            pass
                else:
                    try:
                        self.plotter.remove_actor(self.current_atom_info_labels)
                    except:
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
                    except:
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
            except:
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



