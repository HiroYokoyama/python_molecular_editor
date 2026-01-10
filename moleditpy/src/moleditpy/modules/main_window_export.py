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
main_window_export.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowExport
"""


import numpy as np
import math
import os


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QFileDialog, QMessageBox
)

from PyQt6.QtGui import (
    QBrush, QColor, QPainter, QImage
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
    from .atom_item import AtomItem
    from .bond_item import BondItem
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem


# --- クラス定義 ---
class MainWindowExport(object):
    """ main_window.py から分離された機能クラス """


    def export_stl(self):
        """STLファイルとしてエクスポート（色なし）"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから直接データを取得（色情報なし）
            combined_mesh = self.export_from_3d_view_no_color()
            
            if combined_mesh is None or combined_mesh.n_points == 0:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            if not file_path.lower().endswith('.stl'):
                file_path += '.stl'
            
            combined_mesh.save(file_path, binary=True)
            self.statusBar().showMessage(f"STL exported to {file_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting STL: {e}")



    def export_obj_mtl(self):
        """OBJ/MTLファイルとしてエクスポート（表示中のモデルベース、色付き）"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as OBJ/MTL (with colors)", default_dir, "OBJ Files (*.obj);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから表示中のメッシュデータを色情報とともに取得
            meshes_with_colors = self.export_from_3d_view_with_colors()
            
            if not meshes_with_colors:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            # ファイル拡張子を確認・追加
            if not file_path.lower().endswith('.obj'):
                file_path += '.obj'
            
            # OBJ+MTL形式で保存（オブジェクトごとに色分け）
            mtl_path = file_path.replace('.obj', '.mtl')
            
            self.create_multi_material_obj(meshes_with_colors, file_path, mtl_path)
            
            self.statusBar().showMessage(f"OBJ+MTL files with individual colors exported to {file_path} and {mtl_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting OBJ/MTL: {e}")



    def create_multi_material_obj(self, meshes_with_colors, obj_path, mtl_path):
        """複数のマテリアルを持つOBJファイルとMTLファイルを作成（改良版）"""
        try:
            
            # MTLファイルを作成
            with open(mtl_path, 'w') as mtl_file:
                mtl_file.write(f"# Material file for {os.path.basename(obj_path)}\n")
                mtl_file.write("# Generated with individual object colors\n\n")
                
                for i, mesh_data in enumerate(meshes_with_colors):
                    color = mesh_data['color']
                    material_name = f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    
                    mtl_file.write(f"newmtl {material_name}\n")
                    mtl_file.write(f"Ka 0.2 0.2 0.2\n")  # Ambient
                    mtl_file.write(f"Kd {color[0]/255.0:.3f} {color[1]/255.0:.3f} {color[2]/255.0:.3f}\n")  # Diffuse
                    mtl_file.write(f"Ks 0.5 0.5 0.5\n")  # Specular
                    mtl_file.write(f"Ns 32.0\n")          # Specular exponent
                    mtl_file.write(f"illum 2\n")          # Illumination model
                    mtl_file.write(f"\n")
            
            # OBJファイルを作成
            with open(obj_path, 'w') as obj_file:
                obj_file.write(f"# OBJ file with multiple materials\n")
                obj_file.write(f"# Generated with individual object colors\n")
                obj_file.write(f"mtllib {os.path.basename(mtl_path)}\n\n")
                
                vertex_offset = 1  # OBJファイルの頂点インデックスは1から始まる
                
                for i, mesh_data in enumerate(meshes_with_colors):
                    mesh = mesh_data['mesh']
                    material_name = f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    
                    obj_file.write(f"# Object {i}: {mesh_data['name']}\n")
                    obj_file.write(f"# Color: RGB({mesh_data['color'][0]}, {mesh_data['color'][1]}, {mesh_data['color'][2]})\n")
                    obj_file.write(f"o object_{i}\n")
                    obj_file.write(f"usemtl {material_name}\n")
                    
                    # 頂点を書き込み
                    points = mesh.points
                    for point in points:
                        obj_file.write(f"v {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
                    
                    # 面を書き込み
                    faces_written = 0
                    for j in range(mesh.n_cells):
                        cell = mesh.get_cell(j)
                        if cell.type == 5:  # VTK_TRIANGLE
                            points_in_cell = cell.point_ids
                            v1 = points_in_cell[0] + vertex_offset
                            v2 = points_in_cell[1] + vertex_offset
                            v3 = points_in_cell[2] + vertex_offset
                            obj_file.write(f"f {v1} {v2} {v3}\n")
                            faces_written += 1
                        elif cell.type == 6:  # VTK_TRIANGLE_STRIP
                            # Triangle strips share vertices between adjacent triangles
                            # For n points, we get (n-2) triangles
                            points_in_cell = cell.point_ids
                            n_points = len(points_in_cell)
                            for k in range(n_points - 2):
                                if k % 2 == 0:
                                    # Even triangles: use points k, k+1, k+2
                                    v1 = points_in_cell[k] + vertex_offset
                                    v2 = points_in_cell[k+1] + vertex_offset
                                    v3 = points_in_cell[k+2] + vertex_offset
                                else:
                                    # Odd triangles: reverse winding to maintain consistent orientation
                                    v1 = points_in_cell[k+1] + vertex_offset
                                    v2 = points_in_cell[k] + vertex_offset
                                    v3 = points_in_cell[k+2] + vertex_offset
                                obj_file.write(f"f {v1} {v2} {v3}\n")
                                faces_written += 1
                        elif cell.type == 9:  # VTK_QUAD
                            points_in_cell = cell.point_ids
                            v1 = points_in_cell[0] + vertex_offset
                            v2 = points_in_cell[1] + vertex_offset
                            v3 = points_in_cell[2] + vertex_offset
                            v4 = points_in_cell[3] + vertex_offset
                            obj_file.write(f"f {v1} {v2} {v3} {v4}\n")
                            faces_written += 1
                    
                    
                    vertex_offset += mesh.n_points
                    obj_file.write(f"\n")
                
        except Exception as e:
            raise Exception(f"Failed to create multi-material OBJ: {e}")



    def export_color_stl(self):
        """カラーSTLファイルとしてエクスポート"""
        if not self.current_mol:
            self.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return
            
        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except Exception:
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export as Color STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )
        
        if not file_path:
            return
            
        try:
            
            # 3Dビューから直接データを取得
            combined_mesh = self.export_from_3d_view()
            
            if combined_mesh is None or combined_mesh.n_points == 0:
                self.statusBar().showMessage("No 3D geometry to export.")
                return
            
            # STL形式で保存
            if not file_path.lower().endswith('.stl'):
                file_path += '.stl'
            combined_mesh.save(file_path, binary=True)
            self.statusBar().showMessage(f"STL exported to {file_path}")
                
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting STL: {e}")
    


    def export_from_3d_view(self):
        """現在の3Dビューから直接メッシュデータを取得"""
        try:
            
            # PyVistaプロッターから全てのアクターを取得
            combined_mesh = pv.PolyData()
            
            # プロッターのレンダラーからアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得する複数の方法を試行
                    mesh = None
                    
                    # 方法1: mapperのinputから取得 (Improved)
                    mapper = None
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, 'GetMapper'):
                        mapper = actor.GetMapper()
                    
                    if mapper is not None:
                        if hasattr(mapper, 'input') and mapper.input is not None:
                            mesh = mapper.input
                        elif hasattr(mapper, 'GetInput') and mapper.GetInput() is not None:
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, 'GetInputAsDataSet'):
                            mesh = mapper.GetInputAsDataSet()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # 元のメッシュを変更しないようにコピーを作成
                        mesh_copy = mesh.copy()
                        
                        # コピーしたメッシュにカラー情報を追加
                        if hasattr(actor, 'prop') and hasattr(actor.prop, 'color'):
                            color = actor.prop.color
                            # RGB値を0-255の範囲に変換
                            rgb = np.array([int(c * 255) for c in color], dtype=np.uint8)
                            
                            # Blender対応のPLY形式用カラー属性を設定
                            mesh_copy.point_data['diffuse_red'] = np.full(mesh_copy.n_points, rgb[0], dtype=np.uint8)
                            mesh_copy.point_data['diffuse_green'] = np.full(mesh_copy.n_points, rgb[1], dtype=np.uint8) 
                            mesh_copy.point_data['diffuse_blue'] = np.full(mesh_copy.n_points, rgb[2], dtype=np.uint8)
                            
                            # 標準的なPLY形式もサポート
                            mesh_copy.point_data['red'] = np.full(mesh_copy.n_points, rgb[0], dtype=np.uint8)
                            mesh_copy.point_data['green'] = np.full(mesh_copy.n_points, rgb[1], dtype=np.uint8) 
                            mesh_copy.point_data['blue'] = np.full(mesh_copy.n_points, rgb[2], dtype=np.uint8)
                            
                            # 従来の colors 配列も保持（STL用）
                            mesh_colors = np.tile(rgb, (mesh_copy.n_points, 1))
                            mesh_copy.point_data['colors'] = mesh_colors
                        
                        # メッシュを結合
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)
                            
                except Exception:
                    continue
            
            return combined_mesh
            
        except Exception:
            return None



    def export_from_3d_view_no_color(self):
        """現在の3Dビューから直接メッシュデータを取得（色情報なし）"""
        try:
            
            # PyVistaプロッターから全てのアクターを取得
            combined_mesh = pv.PolyData()
            
            # プロッターのレンダラーからアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得する複数の方法を試行
                    mesh = None
                    
                    # 方法1: mapperのinputから取得 (Improved)
                    mapper = None
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, 'GetMapper'):
                        mapper = actor.GetMapper()
                    
                    if mapper is not None:
                        if hasattr(mapper, 'input') and mapper.input is not None:
                            mesh = mapper.input
                        elif hasattr(mapper, 'GetInput') and mapper.GetInput() is not None:
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, 'GetInputAsDataSet'):
                            mesh = mapper.GetInputAsDataSet()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    # 方法3: Removed unsafe fallback
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # 元のメッシュを変更しないようにコピーを作成（色情報は追加しない）
                        mesh_copy = mesh.copy()
                        
                        # メッシュを結合
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)
                            
                except Exception:
                    continue
            
            return combined_mesh
            
        except Exception:
            return None



    def export_from_3d_view_with_colors(self):
        """現在の3Dビューから直接メッシュデータを色情報とともに取得"""
        try:
            
            meshes_with_colors = []
            
            # PyVistaプロッターから全てのアクターを取得
            renderer = self.plotter.renderer
            actors = renderer.actors
            
            actor_count = 0
            
            for actor_name, actor in actors.items():
                try:
                    # VTKアクターからポリデータを取得
                    mesh = None
                    
                    # 方法1: mapperのinputから取得 (Improved)
                    mapper = None
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, 'GetMapper'):
                        mapper = actor.GetMapper()
                    
                    if mapper is not None:
                        if hasattr(mapper, 'input') and mapper.input is not None:
                            mesh = mapper.input
                        elif hasattr(mapper, 'GetInput') and mapper.GetInput() is not None:
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, 'GetInputAsDataSet'):
                            mesh = mapper.GetInputAsDataSet()
                    
                    # 方法2: PyVistaプロッターの内部データから取得
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]
                    
                    if mesh is not None and hasattr(mesh, 'n_points') and mesh.n_points > 0:
                        # PyVistaメッシュに変換（必要な場合）
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, 'extract_surface'):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)
                        
                        # アクターから色情報を取得
                        color = [128, 128, 128]  # デフォルト色（グレー）
                        
                        try:
                            # VTKアクターのプロパティから色を取得
                            if hasattr(actor, 'prop') and actor.prop is not None:
                                vtk_color = actor.prop.GetColor()
                                color = [int(c * 255) for c in vtk_color]
                            elif hasattr(actor, 'GetProperty'):
                                prop = actor.GetProperty()
                                if prop is not None:
                                    vtk_color = prop.GetColor()
                                    color = [int(c * 255) for c in vtk_color]
                        except Exception:
                            # 色取得に失敗した場合はデフォルト色をそのまま使用
                            pass
                        
                        # メッシュのコピーを作成
                        mesh_copy = mesh.copy()

                        # もしメッシュに頂点ごとの色情報が含まれている場合、
                        # それぞれの色ごとにサブメッシュに分割して個別マテリアルを作る。
                        # これにより、glyphs（すべての原子が一つのメッシュにまとめられる場合）でも
                        # 各原子の色を保持してOBJ/MTLへ出力できる。
                        try:
                            colors = None
                            pd = mesh_copy.point_data
                            # 優先的にred/green/blue配列を使用
                            if 'red' in pd and 'green' in pd and 'blue' in pd:
                                r = np.asarray(pd['red']).reshape(-1)
                                g = np.asarray(pd['green']).reshape(-1)
                                b = np.asarray(pd['blue']).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # diffuse_* のキーもサポート
                            elif 'diffuse_red' in pd and 'diffuse_green' in pd and 'diffuse_blue' in pd:
                                r = np.asarray(pd['diffuse_red']).reshape(-1)
                                g = np.asarray(pd['diffuse_green']).reshape(-1)
                                b = np.asarray(pd['diffuse_blue']).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # 単一の colors 配列があればそれを使う
                            elif 'colors' in pd:
                                colors = np.asarray(pd['colors'])
                            
                            # cell_dataのcolorsも確認（Tubeフィルタなどはcell_dataに色を持つ場合がある）
                            if colors is None and 'colors' in mesh_copy.cell_data:
                                try:
                                    # cell_dataをpoint_dataに変換
                                    temp_mesh = mesh_copy.cell_data_to_point_data()
                                    if 'colors' in temp_mesh.point_data:
                                        colors = np.asarray(temp_mesh.point_data['colors'])
                                except Exception:
                                    pass

                            if colors is not None and colors.size > 0:
                                # 整数に変換。colors が 0-1 の float の場合は 255 倍して正規化する。
                                colors_arr = np.asarray(colors)
                                # 期待形状に整形
                                if colors_arr.ndim == 1:
                                    # 1次元の場合は単一チャンネルとして扱う
                                    colors_arr = colors_arr.reshape(-1, 1)

                                # float かどうか判定して正規化
                                if np.issubdtype(colors_arr.dtype, np.floating):
                                    # 値の最大が1付近なら0-1レンジとみなして255倍
                                    if colors_arr.max() <= 1.01:
                                        colors_int = np.clip((colors_arr * 255.0).round(), 0, 255).astype(np.int32)
                                    else:
                                        # 既に0-255レンジのfloatならそのまま丸める
                                        colors_int = np.clip(colors_arr.round(), 0, 255).astype(np.int32)
                                else:
                                    colors_int = np.clip(colors_arr, 0, 255).astype(np.int32)
                                # Ensure shape is (n_points, 3)
                                if colors_int.ndim == 1:
                                    # 単一値が入っている場合は同一RGBとして扱う
                                    colors_int = np.vstack([colors_int, colors_int, colors_int]).T

                                # 一意な色ごとにサブメッシュを抽出して追加
                                unique_colors, inverse = np.unique(colors_int, axis=0, return_inverse=True)
                                
                                split_success = False
                                if unique_colors.shape[0] > 1:
                                    for uc_idx, uc in enumerate(unique_colors):
                                        point_inds = np.where(inverse == uc_idx)[0]
                                        if point_inds.size == 0:
                                            continue
                                        try:
                                            # Use temp_mesh if available (has point data), else mesh_copy
                                            target_mesh = temp_mesh if 'temp_mesh' in locals() else mesh_copy
                                            
                                            # extract_points with adjacent_cells=False to avoid pulling in neighbors
                                            submesh = target_mesh.extract_points(point_inds, adjacent_cells=False) 
                                            
                                        except Exception:
                                            # extract_points が利用できない場合はスキップ
                                            continue
                                        if submesh is None or getattr(submesh, 'n_points', 0) == 0:
                                            continue
                                        
                                        color_rgb = [int(uc[0]), int(uc[1]), int(uc[2])]
                                        meshes_with_colors.append({
                                            'mesh': submesh,
                                            'color': color_rgb,
                                            'name': f'{actor_name}_color_{uc_idx}',
                                            'type': 'display_actor',
                                            'actor_name': actor_name
                                        })
                                        split_success = True
                                    
                                    if split_success:
                                        actor_count += 1
                                        # 分割に成功したので以下の通常追加は行わない
                                        continue
                                    # If splitting failed (no submeshes added), fall through to default
                                else:
                                    # 色が1色のみの場合は、その色を使用してメッシュ全体を出力
                                    uc = unique_colors[0]
                                    color = [int(uc[0]), int(uc[1]), int(uc[2])]
                                    # ここでは continue せず、下のデフォルト追加処理に任せる（colorを更新したため）
                        except Exception:
                            # 分割処理に失敗した場合はフォールバックで単体メッシュを追加
                            pass
                        
                        meshes_with_colors.append({
                            'mesh': mesh_copy,
                            'color': color,
                            'name': f'actor_{actor_count}_{actor_name}',
                            'type': 'display_actor',
                            'actor_name': actor_name
                        })
                        
                        actor_count += 1
                            
                except Exception as e:
                    continue
            
            
            return meshes_with_colors
            
        except Exception as e:
            print(f"Error in export_from_3d_view_with_colors: {e}")
            return []



    def export_2d_png(self):
        if not self.data.atoms:
            self.statusBar().showMessage("Nothing to export.")
            return

        # default filename: based on current file, append -2d for 2D exports
        default_name = "untitled-2d"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}-2d"
        except Exception:
            default_name = "untitled-2d"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
        except Exception:
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(self, "Export 2D as PNG", default_path, "PNG Files (*.png)")
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for a white background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.Yes)

        if reply == QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = (reply == QMessageBox.StandardButton.Yes)

        QApplication.processEvents()

        items_to_restore = {}
        original_background = self.scene.backgroundBrush()

        try:
            all_items = list(self.scene.items())
            for item in all_items:
                is_mol_part = isinstance(item, (AtomItem, BondItem))
                if not (is_mol_part and item.isVisible()):
                    items_to_restore[item] = item.isVisible()
                    item.hide()

            molecule_bounds = QRectF()
            for item in self.scene.items():
                if isinstance(item, (AtomItem, BondItem)) and item.isVisible():
                    molecule_bounds = molecule_bounds.united(item.sceneBoundingRect())

            if molecule_bounds.isEmpty() or not molecule_bounds.isValid():
                self.statusBar().showMessage("Error: Could not determine molecule bounds for export.")
                return

            if is_transparent:
                self.scene.setBackgroundBrush(QBrush(Qt.BrushStyle.NoBrush))
            else:
                self.scene.setBackgroundBrush(QBrush(QColor("#FFFFFF")))

            rect_to_render = molecule_bounds.adjusted(-20, -20, 20, 20)

            w = max(1, int(math.ceil(rect_to_render.width())))
            h = max(1, int(math.ceil(rect_to_render.height())))

            if w <= 0 or h <= 0:
                self.statusBar().showMessage("Error: Invalid image size calculated.")
                return

            image = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
            if is_transparent:
                image.fill(Qt.GlobalColor.transparent)
            else:
                image.fill(Qt.GlobalColor.white)

            painter = QPainter()
            ok = painter.begin(image)
            if not ok or not painter.isActive():
                self.statusBar().showMessage("Failed to start QPainter for image rendering.")
                return

            try:
                painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                target_rect = QRectF(0, 0, w, h)
                source_rect = rect_to_render
                self.scene.render(painter, target_rect, source_rect)
            finally:
                painter.end()

            saved = image.save(filePath, "PNG")
            if saved:
                self.statusBar().showMessage(f"2D view exported to {filePath}")
            else:
                self.statusBar().showMessage("Failed to save image. Check file path or permissions.")

        except Exception as e:
            self.statusBar().showMessage(f"An unexpected error occurred during 2D export: {e}")

        finally:
            for item, was_visible in items_to_restore.items():
                item.setVisible(was_visible)
            self.scene.setBackgroundBrush(original_background)
            if self.view_2d:
                self.view_2d.viewport().update()



    def export_3d_png(self):
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to export.", 2000)
            return

        # default filename: match XYZ/MOL naming (use base name without suffix)
        default_name = "untitled"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}"
        except Exception:
            default_name = "untitled"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
        except Exception:
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(self, "Export 3D as PNG", default_path, "PNG Files (*.png)")
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(self, 'Choose Background',
                                     'Do you want a transparent background?\n(Choose "No" for current background)',
                                     QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                                     QMessageBox.StandardButton.Yes)

        if reply == QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = (reply == QMessageBox.StandardButton.Yes)

        try:
            self.plotter.screenshot(filePath, transparent_background=is_transparent)
            self.statusBar().showMessage(f"3D view exported to {filePath}", 3000)
        except Exception as e:
            self.statusBar().showMessage(f"Error exporting 3D PNG: {e}")


