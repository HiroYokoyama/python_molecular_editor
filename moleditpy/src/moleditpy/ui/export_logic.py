#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations
import math
import os
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

# PyQt6 Modules
import pyvista as pv
from PyQt6.QtCore import QRectF, QSize, Qt
from PyQt6.QtGui import QBrush, QImage, QPainter
from PyQt6.QtSvg import QSvgGenerator
from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .atom_item import AtomItem
    from .bond_item import BondItem
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.atom_item import AtomItem
    from moleditpy.ui.bond_item import BondItem


# --- Class Definition ---
class ExportManager:
    """Independent manager for export logic, ported from MainWindowExport mixin."""

    _cls: Optional[type[ExportManager]] = None

    def __init__(self, host: Any) -> None:
        self.host = host

    def __getattr__(self, name: str) -> Any:
        """Delegate back to host for attributes not found on this manager."""
        return getattr(self.host, name)

    def export_stl(self) -> None:
        if not self.current_mol:
            self.host.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self.host, "Export as STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )

        if not file_path:
            return

        try:
            # Get 3D data from view (no color)
            combined_mesh = self.export_from_3d_view_no_color()

            if combined_mesh is None or combined_mesh.n_points == 0:
                self.host.statusBar().showMessage("No 3D geometry to export.")
                return

            if not file_path.lower().endswith(".stl"):
                file_path += ".stl"

            combined_mesh.save(file_path, binary=True)
            self.host.statusBar().showMessage(f"STL exported to {file_path}")

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(f"Error exporting STL: {e}")

    def export_obj_mtl(self) -> None:
        """Export as OBJ/MTL (with colors)."""
        if not self.current_mol:
            self.host.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self.host,
            "Export as OBJ/MTL (with colors)",
            default_dir,
            "OBJ Files (*.obj);;All Files (*)",
        )

        if not file_path:
            return

        try:
            # Get mesh data with colors from 3D view
            meshes_with_colors = self.export_from_3d_view_with_colors()

            if not meshes_with_colors:
                self.host.statusBar().showMessage("No 3D geometry to export.")
                return

            # Ensure file extension
            if not file_path.lower().endswith(".obj"):
                file_path += ".obj"

            # Save as OBJ+MTL with material per object
            mtl_path = file_path.replace(".obj", ".mtl")

            self.create_multi_material_obj(meshes_with_colors, file_path, mtl_path)

            self.host.statusBar().showMessage(
                f"OBJ+MTL files with individual colors exported to {file_path} and {mtl_path}"
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(f"Error exporting OBJ/MTL: {e}")

    def create_multi_material_obj(self, meshes_with_colors: List[Dict[str, Any]], obj_path: str, mtl_path: str) -> None:
        """Create multi-material OBJ/MTL files.
        meshes_with_colors: list of dicts with 'mesh', 'color', 'name'
        """
        try:
            # Create MTL file
            with open(mtl_path, "w") as mtl_file:
                mtl_file.write(f"# Material file for {os.path.basename(obj_path)}\n")
                mtl_file.write("# Generated with individual object colors\n\n")

                for i, mesh_data in enumerate(meshes_with_colors):
                    color = mesh_data["color"]
                    material_name = (
                        f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    )

                    mtl_file.write(f"newmtl {material_name}\n")
                    mtl_file.write("Ka 0.2 0.2 0.2\n")  # Ambient
                    mtl_file.write(
                        f"Kd {color[0] / 255.0:.3f} {color[1] / 255.0:.3f} {color[2] / 255.0:.3f}\n"
                    )  # Diffuse
                    mtl_file.write("Ks 0.5 0.5 0.5\n")  # Specular
                    mtl_file.write("Ns 32.0\n")  # Specular exponent
                    mtl_file.write("illum 2\n")  # Illumination model
                    mtl_file.write("\n")

            # Create OBJ file
            with open(obj_path, "w") as obj_file:
                obj_file.write("# OBJ file with multiple materials\n")
                obj_file.write("# Generated with individual object colors\n")
                obj_file.write(f"mtllib {os.path.basename(mtl_path)}\n\n")

                vertex_offset = 1  # OBJ indices start at 1

                for i, mesh_data in enumerate(meshes_with_colors):
                    mesh = mesh_data["mesh"]
                    material_name = (
                        f"material_{i}_{mesh_data['name'].replace(' ', '_')}"
                    )

                    obj_file.write(f"# Object {i}: {mesh_data['name']}\n")
                    obj_file.write(
                        f"# Color: RGB({mesh_data['color'][0]}, {mesh_data['color'][1]}, {mesh_data['color'][2]})\n"
                    )
                    obj_file.write(f"o object_{i}\n")
                    obj_file.write(f"usemtl {material_name}\n")

                    # Write vertices
                    points = mesh.points
                    for point in points:
                        obj_file.write(
                            f"v {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n"
                        )

                    # Write faces
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
                                    v2 = points_in_cell[k + 1] + vertex_offset
                                    v3 = points_in_cell[k + 2] + vertex_offset
                                else:
                                    # Odd triangles: reverse winding to maintain consistent orientation
                                    v1 = points_in_cell[k + 1] + vertex_offset
                                    v2 = points_in_cell[k] + vertex_offset
                                    v3 = points_in_cell[k + 2] + vertex_offset
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
                    obj_file.write("\n")

        except (AttributeError, RuntimeError, ValueError) as e:
            raise Exception(f"Failed to create multi-material OBJ: {e}")

    def export_color_stl(self) -> None:
        """Export as Color STL."""
        if not self.current_mol:
            self.host.statusBar().showMessage("Error: Please generate a 3D structure first.")
            return

        # prefer same directory as current file when available
        default_dir = ""
        try:
            if self.current_file_path:
                default_dir = os.path.dirname(self.current_file_path)
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_dir = ""

        file_path, _ = QFileDialog.getSaveFileName(
            self.host, "Export as Color STL", default_dir, "STL Files (*.stl);;All Files (*)"
        )

        if not file_path:
            return

        try:
            # Get 3D data from view
            combined_mesh = self.export_from_3d_view()

            if combined_mesh is None or combined_mesh.n_points == 0:
                self.host.statusBar().showMessage("No 3D geometry to export.")
                return

            # Save as STL
            if not file_path.lower().endswith(".stl"):
                file_path += ".stl"
            combined_mesh.save(file_path, binary=True)
            self.host.statusBar().showMessage(f"STL exported to {file_path}")

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(f"Error exporting STL: {e}")

    def export_from_3d_view(self) -> Optional[pv.PolyData]:
        """Get mesh data from 3D view."""
        try:
            # Get all actors from PyVista plotter
            combined_mesh = pv.PolyData()

            # Get actors from renderer
            renderer = self.plotter.renderer
            actors = renderer.actors

            for actor_name, actor in actors.items():
                try:
                    # Attempt to get polydata from VTK actor
                    mesh = None

                    # Method 1: Get from mapper input
                    mapper = None
                    if hasattr(actor, "mapper") and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, "GetMapper"):
                        mapper = actor.GetMapper()

                    if mapper is not None:
                        if hasattr(mapper, "input") and mapper.input is not None:
                            mesh = mapper.input
                        elif (
                            hasattr(mapper, "GetInput")
                            and mapper.GetInput() is not None
                        ):
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, "GetInputAsDataSet"):
                            mesh = mapper.GetInputAsDataSet()

                    # Method 2: Get from PyVista plotter internal data
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]

                    if (
                        mesh is not None
                        and hasattr(mesh, "n_points")
                        and mesh.n_points > 0
                    ):
                        # Convert to PyVista mesh (if necessary)
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, "extract_surface"):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)

                        # Create copy of mesh
                        mesh_copy = mesh.copy()

                        # Add color info
                        if hasattr(actor, "prop") and hasattr(actor.prop, "color"):
                            color = actor.prop.color
                            # Convert RGB to 0-255
                            rgb = np.array(
                                [int(c * 255) for c in color], dtype=np.uint8
                            )

                            # Set PLY color attributes for Blender
                            mesh_copy.point_data["diffuse_red"] = np.full(
                                mesh_copy.n_points, rgb[0], dtype=np.uint8
                            )
                            mesh_copy.point_data["diffuse_green"] = np.full(
                                mesh_copy.n_points, rgb[1], dtype=np.uint8
                            )
                            mesh_copy.point_data["diffuse_blue"] = np.full(
                                mesh_copy.n_points, rgb[2], dtype=np.uint8
                            )

                            # Support standard PLY
                            mesh_copy.point_data["red"] = np.full(
                                mesh_copy.n_points, rgb[0], dtype=np.uint8
                            )
                            mesh_copy.point_data["green"] = np.full(
                                mesh_copy.n_points, rgb[1], dtype=np.uint8
                            )
                            mesh_copy.point_data["blue"] = np.full(
                                mesh_copy.n_points, rgb[2], dtype=np.uint8
                            )

                            # Keep 'colors' array for STL
                            mesh_colors = np.tile(rgb, (mesh_copy.n_points, 1))
                            mesh_copy.point_data["colors"] = mesh_colors

                        # Combine meshes
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)

                except (AttributeError, RuntimeError, ValueError, TypeError):
                    continue

            return combined_mesh

        except (AttributeError, RuntimeError, ValueError, TypeError):
            return None

    def export_from_3d_view_no_color(self) -> Optional[pv.PolyData]:
        """Get mesh data from 3D view (no color)."""
        try:
            # Get actors from PyVista plotter
            combined_mesh = pv.PolyData()

            # Get actors from renderer
            renderer = self.plotter.renderer
            actors = renderer.actors

            for actor_name, actor in actors.items():
                try:
                    # Attempt to get polydata from VTK actor
                    mesh = None

                    # Method 1: Get from mapper input
                    mapper = None
                    if hasattr(actor, "mapper") and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, "GetMapper"):
                        mapper = actor.GetMapper()

                    if mapper is not None:
                        if hasattr(mapper, "input") and mapper.input is not None:
                            mesh = mapper.input
                        elif (
                            hasattr(mapper, "GetInput")
                            and mapper.GetInput() is not None
                        ):
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, "GetInputAsDataSet"):
                            mesh = mapper.GetInputAsDataSet()

                    # Method 2: Get from PyVista plotter internal data
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]

                    # Method 3: Removed unsafe fallback

                    if (
                        mesh is not None
                        and hasattr(mesh, "n_points")
                        and mesh.n_points > 0
                    ):
                        # Convert to PyVista mesh if necessary
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, "extract_surface"):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)

                        # Create copy of mesh
                        mesh_copy = mesh.copy()

                        # Combine meshes
                        if combined_mesh.n_points == 0:
                            combined_mesh = mesh_copy.copy()
                        else:
                            combined_mesh = combined_mesh.merge(mesh_copy)

                except (AttributeError, RuntimeError, ValueError, TypeError):
                    continue

            return combined_mesh

        except (AttributeError, RuntimeError, ValueError, TypeError):
            return None

    def export_from_3d_view_with_colors(self) -> List[Dict[str, Any]]:
        """Get mesh data with colors from 3D view."""
        try:
            meshes_with_colors = []

            # Get actors from PyVista plotter
            renderer = self.plotter.renderer
            actors = renderer.actors

            actor_count = 0

            for actor_name, actor in actors.items():
                try:
                    # Get polydata from VTK actor
                    mesh = None

                    # Method 1: Get from mapper input (Improved)
                    mapper = None
                    if hasattr(actor, "mapper") and actor.mapper is not None:
                        mapper = actor.mapper
                    elif hasattr(actor, "GetMapper"):
                        mapper = actor.GetMapper()

                    if mapper is not None:
                        if hasattr(mapper, "input") and mapper.input is not None:
                            mesh = mapper.input
                        elif (
                            hasattr(mapper, "GetInput")
                            and mapper.GetInput() is not None
                        ):
                            mesh = mapper.GetInput()
                        elif hasattr(mapper, "GetInputAsDataSet"):
                            mesh = mapper.GetInputAsDataSet()

                    # Method 2: Get from PyVista plotter internal data
                    if mesh is None and actor_name in self.plotter.mesh:
                        mesh = self.plotter.mesh[actor_name]

                    if (
                        mesh is not None
                        and hasattr(mesh, "n_points")
                        and mesh.n_points > 0
                    ):
                        # Convert to PyVista mesh if necessary
                        if not isinstance(mesh, pv.PolyData):
                            if hasattr(mesh, "extract_surface"):
                                mesh = mesh.extract_surface()
                            else:
                                mesh = pv.wrap(mesh)

                        # Get color from actor
                        color = [128, 128, 128]  # Default color (gray)

                        try:
                            # Get color from properties
                            if hasattr(actor, "prop") and actor.prop is not None:
                                vtk_color = actor.prop.GetColor()
                                color = [int(c * 255) for c in vtk_color]
                            elif hasattr(actor, "GetProperty"):
                                prop = actor.GetProperty()
                                if prop is not None:
                                    vtk_color = prop.GetColor()
                                    color = [int(c * 255) for c in vtk_color]
                        except (AttributeError, RuntimeError, TypeError):
                            # Use default color on failure to avoid console noise during complex mesh export
                            pass

                        # Create mesh copy
                        mesh_copy = mesh.copy()

                        # Split into submeshes by color if vertex colors exist
                        # This allows glyphs (where all atoms are in one mesh)
                        # to retain individual atom colors when exporting to OBJ/MTL.
                        try:
                            colors = None
                            pd = mesh_copy.point_data
                            # Use red/green/blue arrays if available
                            if "red" in pd and "green" in pd and "blue" in pd:
                                r = np.asarray(pd["red"]).reshape(-1)
                                g = np.asarray(pd["green"]).reshape(-1)
                                b = np.asarray(pd["blue"]).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # Support diffuse_* keys
                            elif (
                                "diffuse_red" in pd
                                and "diffuse_green" in pd
                                and "diffuse_blue" in pd
                            ):
                                r = np.asarray(pd["diffuse_red"]).reshape(-1)
                                g = np.asarray(pd["diffuse_green"]).reshape(-1)
                                b = np.asarray(pd["diffuse_blue"]).reshape(-1)
                                colors = np.vstack([r, g, b]).T
                            # Use 'colors' array if available
                            elif "colors" in pd:
                                colors = np.asarray(pd["colors"])

                            # Check cell_data colors
                            if colors is None and "colors" in mesh_copy.cell_data:
                                try:
                                    # Convert cell_data to point_data
                                    temp_mesh = mesh_copy.cell_data_to_point_data()
                                    if "colors" in temp_mesh.point_data:
                                        colors = np.asarray(
                                            temp_mesh.point_data["colors"]
                                        )
                                except (
                                    AttributeError,
                                    RuntimeError,
                                    ValueError,
                                    TypeError,
                                ):
                                    # Fail silently and fall through to default mesh addition
                                    pass
                            if colors is not None and colors.size > 0:
                                # Normalize float colors to 0-255
                                colors_arr = np.asarray(colors)
                                # Reshape to expected shape
                                if colors_arr.ndim == 1:
                                    # Treat as a single channel if 1D
                                    colors_arr = colors_arr.reshape(-1, 1)

                                # Normalize if float
                                if np.issubdtype(colors_arr.dtype, np.floating):
                                    # Assume 0-1 range if max <= 1.01
                                    if colors_arr.max() <= 1.01:
                                        colors_int = np.clip(
                                            (colors_arr * 255.0).round(), 0, 255
                                        ).astype(np.int32)
                                    else:
                                        # Round if already 0-255 float
                                        colors_int = np.clip(
                                            colors_arr.round(), 0, 255
                                        ).astype(np.int32)
                                else:
                                    colors_int = np.clip(colors_arr, 0, 255).astype(
                                        np.int32
                                    )
                                # Ensure shape is (n_points, 3)
                                if colors_int.ndim == 1:
                                    # Treat single color value as grayscale/same RGB
                                    colors_int = np.vstack(
                                        [colors_int, colors_int, colors_int]
                                    ).T

                                # Extract submeshes per unique color
                                unique_colors, inverse = np.unique(
                                    colors_int, axis=0, return_inverse=True
                                )

                                split_success = False
                                if unique_colors.shape[0] > 1:
                                    for uc_idx, uc in enumerate(unique_colors):
                                        point_inds = np.where(inverse == uc_idx)[0]
                                        if point_inds.size == 0:
                                            continue
                                        try:
                                            # Use temp_mesh if point data available
                                            target_mesh = (
                                                temp_mesh
                                                if "temp_mesh" in locals()
                                                else mesh_copy
                                            )

                                            # extract_points with adjacent_cells=False to avoid pulling in neighbors
                                            submesh = target_mesh.extract_points(
                                                point_inds, adjacent_cells=False
                                            )

                                        except (
                                            AttributeError,
                                            RuntimeError,
                                            TypeError,
                                        ):
                                            # Skip if extraction unavailable
                                            continue
                                        if (
                                            submesh is None
                                            or getattr(submesh, "n_points", 0) == 0
                                        ):
                                            continue

                                        color_rgb = [int(uc[0]), int(uc[1]), int(uc[2])]
                                        meshes_with_colors.append(
                                            {
                                                "mesh": submesh,
                                                "color": color_rgb,
                                                "name": f"{actor_name}_color_{uc_idx}",
                                                "type": "display_actor",
                                                "actor_name": actor_name,
                                            }
                                        )
                                        split_success = True

                                    if split_success:
                                        actor_count += 1
                                        # Skip default addition on success
                                        continue
                                    # If splitting failed (no submeshes added), fall through to default
                                else:
                                    # Case: single color
                                    uc = unique_colors[0]
                                    color = [int(uc[0]), int(uc[1]), int(uc[2])]
                                    # Do not continue here; let the default addition handle it (color has been updated)
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            # Fallback: add single mesh on failure
                            pass

                        meshes_with_colors.append(
                            {
                                "mesh": mesh_copy,
                                "color": color,
                                "name": f"actor_{actor_count}_{actor_name}",
                                "type": "display_actor",
                                "actor_name": actor_name,
                            }
                        )

                        actor_count += 1

                except (AttributeError, RuntimeError, ValueError):
                    continue

            return meshes_with_colors

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error in export_from_3d_view_with_colors: {e}")
            return []

    def export_2d_png(self) -> None:
        if not self.data.atoms:
            self.host.statusBar().showMessage("Nothing to export.")
            return

        # default filename: based on current file, append -2d for 2D exports
        default_name = "untitled-2d"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}-2d"
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_name = "untitled-2d"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(
                    os.path.dirname(self.current_file_path), default_name
                )
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(
            self.host, "Export 2D as PNG", default_path, "PNG Files (*.png)"
        )
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(
            self.host,
            "Choose Background",
            'Do you want a transparent background?\n(Choose "No" to use the current background color)',
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )

        if reply == QMessageBox.StandardButton.Cancel:
            self.host.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = reply == QMessageBox.StandardButton.Yes

        QApplication.processEvents()

        items_to_restore = {}
        original_background = None
        try:
            original_background = self.scene.backgroundBrush()
        except (AttributeError, RuntimeError, ValueError, TypeError):
            # Minimal risk; keep default brush
            pass

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
                self.host.statusBar().showMessage(
                    "Error: Could not determine molecule bounds for export."
                )
                return

            # Handle transparency
            if is_transparent:
                self.scene.setBackgroundBrush(QBrush(Qt.BrushStyle.NoBrush))

            rect_to_render = molecule_bounds.adjusted(-20, -20, 20, 20)

            w = max(1, int(math.ceil(rect_to_render.width())))
            h = max(1, int(math.ceil(rect_to_render.height())))

            if w <= 0 or h <= 0:
                self.host.statusBar().showMessage("Error: Invalid image size calculated.")
                return

            image = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
            # Initialize transparent image
            image.fill(Qt.GlobalColor.transparent)

            painter = QPainter()
            ok = painter.begin(image)
            if not ok or not painter.isActive():
                self.host.statusBar().showMessage(
                    "Failed to start QPainter for image rendering."
                )
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
                self.host.statusBar().showMessage(f"2D view exported to {filePath}")
            else:
                self.host.statusBar().showMessage(
                    "Failed to save image. Check file path or permissions."
                )

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(
                f"An unexpected error occurred during 2D export: {e}"
            )

        finally:
            for item, was_visible in items_to_restore.items():
                item.setVisible(was_visible)
            self.scene.setBackgroundBrush(original_background)
            if self.view_2d:
                self.view_2d.viewport().update()

    def export_2d_svg(self) -> None:
        """Export 2D drawing as SVG."""
        if not self.data.atoms:
            self.host.statusBar().showMessage("Nothing to export.")
            return

        # default filename
        default_name = "untitled-2d"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}-2d"
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_name = "untitled-2d"

        # prefer same directory
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(
                    os.path.dirname(self.current_file_path), default_name
                )
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(
            self.host, "Export 2D as SVG", default_path, "SVG Files (*.svg)"
        )
        if not filePath:
            return

        if not (filePath.lower().endswith(".svg")):
            filePath += ".svg"

        # Ask about transparency
        reply = QMessageBox.question(
            self.host,
            "Choose Background",
            'Do you want a transparent background?\n(Choose "No" to use the current background color)',
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )

        if reply == QMessageBox.StandardButton.Cancel:
            self.host.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = reply == QMessageBox.StandardButton.Yes

        original_background = None
        try:
            # 1. Hide non-molecular items
            items_to_restore = {}
            original_background = self.scene.backgroundBrush()

            all_items = list(self.scene.items())
            for item in all_items:
                is_mol_part = isinstance(item, (AtomItem, BondItem))
                if not (is_mol_part and item.isVisible()):
                    # Hide non-atom/bond items for consistency with PNG export
                    items_to_restore[item] = item.isVisible()
                    item.hide()

            # 2. Calculate bounds
            molecule_bounds = QRectF()
            for item in self.scene.items():
                if isinstance(item, (AtomItem, BondItem)) and item.isVisible():
                    molecule_bounds = molecule_bounds.united(item.sceneBoundingRect())

            if molecule_bounds.isEmpty() or not molecule_bounds.isValid():
                self.host.statusBar().showMessage(
                    "Error: Could not determine molecule bounds for export."
                )
                # Restore
                for item, was_visible in items_to_restore.items():
                    item.setVisible(was_visible)
                return

            if is_transparent:
                self.scene.setBackgroundBrush(QBrush(Qt.BrushStyle.NoBrush))

            # Margin
            rect_to_render = molecule_bounds.adjusted(-20, -20, 20, 20)

            width = int(rect_to_render.width())
            height = int(rect_to_render.height())

            # 3. Setup QSvgGenerator
            generator = QSvgGenerator()
            generator.setFileName(filePath)
            generator.setSize(QSize(width, height))
            generator.setViewBox(rect_to_render)
            generator.setTitle("MoleditPy Molecule")

            # 4. Render
            painter = QPainter()
            painter.begin(generator)
            try:
                self.scene.render(painter, rect_to_render, rect_to_render)
            finally:
                painter.end()

            self.host.statusBar().showMessage(f"2D view exported to {filePath}")

        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(
                f"An unexpected error occurred during SVG export: {e}"
            )

        finally:
            # Restore
            for item, was_visible in items_to_restore.items():
                item.setVisible(was_visible)
            if original_background is not None:
                self.scene.setBackgroundBrush(original_background)
            if self.view_2d:
                self.view_2d.viewport().update()

    def export_3d_png(self) -> None:
        """Export 3D view as PNG."""
        if not self.current_mol:
            self.host.statusBar().showMessage("No 3D molecule to export.", 2000)
            return

        # Default filename: {name}.png
        default_name = "untitled"
        try:
            if self.current_file_path:
                base = os.path.basename(self.current_file_path)
                name = os.path.splitext(base)[0]
                default_name = f"{name}"
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_name = "untitled"

        # prefer same directory as current file when available
        default_path = default_name
        try:
            if self.current_file_path:
                default_path = os.path.join(
                    os.path.dirname(self.current_file_path), default_name
                )
        except (AttributeError, RuntimeError, ValueError, TypeError):
            default_path = default_name

        filePath, _ = QFileDialog.getSaveFileName(
            self.host, "Export 3D as PNG", default_path, "PNG Files (*.png)"
        )
        if not filePath:
            return

        if not (filePath.lower().endswith(".png")):
            filePath += ".png"

        reply = QMessageBox.question(
            self.host,
            "Choose Background",
            'Do you want a transparent background?\n(Choose "No" for current background)',
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )

        if reply == QMessageBox.StandardButton.Cancel:
            self.host.statusBar().showMessage("Export cancelled.", 2000)
            return

        is_transparent = reply == QMessageBox.StandardButton.Yes

        try:
            self.plotter.screenshot(filePath, transparent_background=is_transparent)
            self.host.statusBar().showMessage(f"3D view exported to {filePath}", 3000)
        except (AttributeError, RuntimeError, ValueError) as e:
            self.host.statusBar().showMessage(f"Error exporting 3D PNG: {e}")


ExportManager._cls = ExportManager
