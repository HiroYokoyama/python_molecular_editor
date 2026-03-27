#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""


import logging
import contextlib

import numpy as np
import vtk

# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem

# PyQt6 Modules
import pyvista as pv
from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QColor, QTransform
from PyQt6.QtWidgets import QGraphicsView

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except (ImportError, AttributeError, TypeError):
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import CPK_COLORS_PV, DEFAULT_CPK_COLORS, VDW_RADII, pt
    from .template_preview_item import TemplatePreviewItem
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.utils.constants import CPK_COLORS_PV, VDW_RADII, pt
    from moleditpy.ui.template_preview_item import TemplatePreviewItem


# --- Class Definition ---
class MainWindowView3d:
    _cls = None  # Class-level reference for plugin patching accessibility

    def set_3d_style(self, style_name):
        """Set 3D display style and update view"""
        current_stored_style = getattr(self, "current_3d_style", None)
        if current_stored_style == style_name:
            return

        # Reset measurement and 3D edit modes on style change
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # Disable measurement mode

        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # Disable 3D edit mode

        # Clear 3D selection
        self.clear_3d_selection()

        self.current_3d_style = style_name
        self.statusBar().showMessage(f"3D style set to: {style_name}")

        # Redraw if molecule is displayed
        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)

    def draw_molecule_3d(self, mol):
        """Dispatch to custom style or standard drawing."""
        mw = self

        if hasattr(mw, "plugin_manager") and hasattr(
            mw.plugin_manager, "custom_3d_styles"
        ):
            if (
                hasattr(self, "current_3d_style")
                and self.current_3d_style in mw.plugin_manager.custom_3d_styles
            ):
                handler = mw.plugin_manager.custom_3d_styles[self.current_3d_style][
                    "callback"
                ]
                try:
                    handler(mw, mol)
                    return
                except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                    logging.error(
                        f"Error in custom 3d style '{self.current_3d_style}': {e}"
                    )

        self.draw_standard_3d_style(mol)

    def draw_standard_3d_style(self, mol, style_override=None):
        """Draw 3D molecule and clear axis actor reference (re-controlled by apply_3d_settings)"""

        # Re-entrancy guard: prevent overlapping draws that cause ghost
        # bonds when slider events are processed mid-render.
        if getattr(self, "_drawing_3d", False):
            return
        self._drawing_3d = True

        try:
            MainWindowView3d._draw_standard_3d_style_body(self, mol, style_override)
        finally:
            self._drawing_3d = False

    def _prepare_3d_kekule_mol(self, mol):
        """Optionally kekulize aromatic systems for 3D visualization."""
        if self.settings.get("display_kekule_3d", False):
            try:
                mol_to_draw = Chem.Mol(mol)
                Chem.Kekulize(mol_to_draw, clearAromaticFlags=True)
                return mol_to_draw
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                # Kekulize failed; keep original and warn user
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    self.statusBar().showMessage(f"Kekulize failed: {e}")
        return mol

    def _draw_standard_3d_style_body(self, mol, style_override=None):
        current_style = (
            style_override
            if style_override
            else getattr(self, "current_3d_style", "Ball and Stick")
        )

        # Clear measurement selection (molecule changed)
        if hasattr(self, "measurement_mode"):
            self.clear_measurement_selection()

        # Initialize 3D color map
        if not hasattr(self, "_3d_color_map"):
            self._3d_color_map = {}
        self._3d_color_map.clear()

        # 1. Camera state and clear
        camera_state = self.plotter.camera.copy()

        # **Force removal to prevent residue**
        # Pylint: access-member-before-definition fix
        old_axes_actor = getattr(self, "axes_actor", None)
        if old_axes_actor is not None:
            with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                self.plotter.remove_actor(old_axes_actor)
            self.axes_actor = None

        self.plotter.clear()

        # 2. Background color setting
        self.plotter.set_background(self.settings.get("background_color", "#4f4f4f"))

        # 3. End with background and axes if mol is None or empty
        if mol is None or mol.GetNumAtoms() == 0:
            self.atom_actor = None
            self.current_mol = None
            self.plotter.render()
            return

        # 4. Lighting setting
        is_lighting_enabled = self.settings.get("lighting_enabled", True)

        if is_lighting_enabled:
            light = pv.Light(
                position=(1, 1, 2),
                light_type="cameralight",
                intensity=self.settings.get("light_intensity", 1.2),
            )
            self.plotter.add_light(light)

        # 5. Molecule drawing logic
        # Optionally kekulize aromatic systems for 3D visualization.
        mol_to_draw = mol
        if self.settings.get("display_kekule_3d", False):
            try:
                mol_to_draw = Chem.Mol(mol)
                Chem.Kekulize(mol_to_draw, clearAromaticFlags=True)
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                # Kekulize failed; keep original and warn user
                with contextlib.suppress(AttributeError, RuntimeError, TypeError):
                    self.statusBar().showMessage(f"Kekulize failed: {e}")
                mol_to_draw = mol

        # Use the original molecule's conformer (positions) to ensure coordinates
        # are preserved even when we create a kekulized copy for bond types.
        conf = mol.GetConformer()

        # Use the kekulized molecule's atom ordering for color/size decisions
        self.atom_positions_3d = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol_to_draw.GetNumAtoms())]
        )

        # Use the possibly-kekulized molecule for symbol/bond types
        sym = [a.GetSymbol() for a in mol_to_draw.GetAtoms()]
        col = np.array([CPK_COLORS_PV.get(s, [0.5, 0.5, 0.5]) for s in sym])

        # Apply plugin color overrides
        if hasattr(self, "_plugin_color_overrides") and self._plugin_color_overrides:
            for atom_idx, hex_color in self._plugin_color_overrides.items():
                if 0 <= atom_idx < len(col):
                    try:
                        c = QColor(hex_color)
                        col[atom_idx] = [c.redF(), c.greenF(), c.blueF()]
                    except (
                        AttributeError,
                        RuntimeError,
                        TypeError,
                        ValueError,
                        KeyError,
                    ):
                        # Suppress traceback
                        pass

        # Define common mesh properties
        mesh_props = dict(
            smooth_shading=True,
            specular=self.settings.get("specular", 0.2),
            specular_power=self.settings.get("specular_power", 20),
            lighting=is_lighting_enabled,
        )

        # --- Mod: Extract variables for delegates ---
        self._add_3d_atom_glyphs(
            mol_to_draw, conf, sym, col, current_style, is_lighting_enabled, mesh_props
        )
        self._add_3d_bond_cylinders(mol_to_draw, conf, col, current_style, mesh_props)
        self._add_3d_aromatic_rings(mol_to_draw, current_style, mesh_props)
        self._add_3d_labels(mol, mol_to_draw)
        self.plotter.camera = camera_state

        # Update projection mode and force render
        settings = getattr(self, "settings", {})
        proj_mode = settings.get("projection_mode", "Perspective")
        if hasattr(self.plotter, "renderer") and hasattr(
            self.plotter.renderer, "GetActiveCamera"
        ):
            vcam = self.plotter.renderer.GetActiveCamera()
            if vcam:
                vcam.SetParallelProjection(proj_mode == "Orthographic")
                try:
                    # Force a render so the change is visible immediately
                    self.plotter.render()
                except (AttributeError, RuntimeError, TypeError):
                    # Suppress non-critical 3D rendering errors
                    pass

        # Re-display if AtomID or other atom info is shown
        if (
            hasattr(self, "atom_info_display_mode")
            and self.atom_info_display_mode is not None
        ):
            self.show_all_atom_info()

        # Update menu text and state depending on molecule type
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()

    def _add_3d_atom_glyphs(
        self,
        mol_to_draw,
        conf,
        sym,
        col,
        current_style,
        is_lighting_enabled,
        mesh_props,
    ):
        # Set atom radii based on style
        if current_style == "cpk":
            atom_scale = self.settings.get("cpk_atom_scale", 1.0)
            resolution = self.settings.get("cpk_resolution", 32)

            # Safe VDW lookup to handle custom elements like 'Bq'
            def get_safe_rvdw(s):
                try:
                    r = pt.GetRvdw(pt.GetAtomicNumber(s))
                    return r if r > 0.1 else 1.5
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    return 1.5

            rad = np.array([get_safe_rvdw(s) * atom_scale for s in sym])
        elif current_style == "wireframe":
            # Atoms not drawn in wireframe mode
            resolution = self.settings.get("wireframe_resolution", 6)
            rad = np.array([0.01 for s in sym])  # Minimal value (not used)
        elif current_style == "stick":
            atom_radius = self.settings.get(
                "stick_bond_radius", 0.15
            )  # Use bond radius for atoms
            resolution = self.settings.get("stick_resolution", 16)
            rad = np.array([atom_radius for s in sym])
        else:  # ball_and_stick
            atom_scale = self.settings.get("ball_stick_atom_scale", 1.0)
            resolution = self.settings.get("ball_stick_resolution", 16)
            rad = np.array([VDW_RADII.get(s, 0.4) * atom_scale for s in sym])

        self.glyph_source = pv.PolyData(self.atom_positions_3d)
        self.glyph_source["colors"] = col
        self.glyph_source["radii"] = rad

        # Do not draw atoms in wireframe mode
        if current_style != "wireframe":
            # Split terminal multiple bond atoms in Stick mode
            if current_style == "stick":
                # Detect terminal atoms with multiple bonds
                split_atoms = []  # (atom_idx, bond_order, offset_vecs)
                skip_atoms = set()  # Indices of atoms to skip

                for i in range(mol_to_draw.GetNumAtoms()):
                    atom = mol_to_draw.GetAtomWithIdx(i)
                    if atom.GetDegree() == 1:  # Terminal atom
                        bonds = atom.GetBonds()
                        if len(bonds) == 1:
                            bond = bonds[0]
                            bond_type = bond.GetBondType()

                            if bond_type in [
                                Chem.BondType.DOUBLE,
                                Chem.BondType.TRIPLE,
                            ]:
                                # Found terminal atom with multiple bond
                                # Get adjacent atom
                                other_idx = (
                                    bond.GetBeginAtomIdx()
                                    if bond.GetEndAtomIdx() == i
                                    else bond.GetEndAtomIdx()
                                )

                                # Calculate bond vector
                                pos_i = np.array(conf.GetAtomPosition(i))
                                pos_other = np.array(conf.GetAtomPosition(other_idx))
                                bond_vec = pos_i - pos_other
                                bond_length = np.linalg.norm(bond_vec)

                                if bond_length > 0:
                                    bond_unit = bond_vec / bond_length

                                    # Use same offset direction as bond drawing for double bonds
                                    if bond_type == Chem.BondType.DOUBLE:
                                        offset_dir1 = (
                                            self._calculate_double_bond_offset(
                                                mol_to_draw, bond, conf
                                            )
                                        )
                                    else:
                                        # Use same logic as bond drawing for triple bonds
                                        v_arb = np.array([0, 0, 1])
                                        if np.allclose(
                                            np.abs(np.dot(bond_unit, v_arb)), 1.0
                                        ):
                                            v_arb = np.array([0, 1, 0])
                                        offset_dir1 = np.cross(bond_unit, v_arb)
                                        offset_dir1 /= np.linalg.norm(offset_dir1)

                                    # Get offset and radius for multiple bonds (matching bond drawing)
                                    try:
                                        cyl_radius = self.settings.get(
                                            "stick_bond_radius", 0.15
                                        )
                                        if bond_type == Chem.BondType.DOUBLE:
                                            radius_factor = self.settings.get(
                                                "stick_double_bond_radius_factor", 0.60
                                            )
                                            offset_factor = self.settings.get(
                                                "stick_double_bond_offset_factor", 1.5
                                            )
                                            # Double bond: use s_double / 2
                                            offset_distance = (
                                                cyl_radius * offset_factor / 2
                                            )
                                        else:  # TRIPLE
                                            radius_factor = self.settings.get(
                                                "stick_triple_bond_radius_factor", 0.40
                                            )
                                            offset_factor = self.settings.get(
                                                "stick_triple_bond_offset_factor", 1.0
                                            )
                                            # Triple bond: use s_triple without division
                                            offset_distance = cyl_radius * offset_factor

                                        # Same calculation as bond drawing
                                        sphere_radius = cyl_radius * radius_factor
                                    except (
                                        AttributeError,
                                        RuntimeError,
                                        TypeError,
                                        ValueError,
                                    ):
                                        sphere_radius = 0.09  # Default values
                                        offset_distance = 0.15  # Default values

                                    if bond_type == Chem.BondType.DOUBLE:
                                        # Double bond: split into 2
                                        offset_vecs = [
                                            offset_dir1 * offset_distance,
                                            -offset_dir1 * offset_distance,
                                        ]
                                        split_atoms.append(
                                            (i, 2, offset_vecs, sphere_radius)
                                        )
                                    else:  # TRIPLE
                                        # Triple bond: split into 3 (center + two sides)
                                        # Same arrangement as bond drawing
                                        offset_vecs = [
                                            np.array([0, 0, 0]),  # Center
                                            offset_dir1 * offset_distance,  # +side
                                            -offset_dir1 * offset_distance,  # -side
                                        ]
                                        split_atoms.append(
                                            (i, 3, offset_vecs, sphere_radius)
                                        )

                                    skip_atoms.add(i)

                # Create new positions if atoms are split
                if split_atoms:
                    new_positions = []
                    new_colors = []
                    new_radii = []

                    # Add normal atoms (excluding skip list)
                    for i in range(len(self.atom_positions_3d)):
                        if i not in skip_atoms:
                            new_positions.append(self.atom_positions_3d[i])
                            new_colors.append(col[i])
                            new_radii.append(rad[i])

                    # Add split atoms
                    # Use calculated sphere_radius (with radius_factor applied)
                    for atom_idx, bond_order, offset_vecs, s_radius in split_atoms:
                        pos = self.atom_positions_3d[atom_idx]
                        # Get radius from bond (calculated above)
                        for offset_vec in offset_vecs:
                            new_positions.append(pos + offset_vec)
                            new_colors.append(col[atom_idx])
                            new_radii.append(s_radius)

                    # Create PolyData at new positions
                    glyph_source = pv.PolyData(np.array(new_positions))
                    glyph_source["colors"] = np.array(new_colors)
                    glyph_source["radii"] = np.array(new_radii)
                else:
                    glyph_source = self.glyph_source
            else:
                glyph_source = self.glyph_source

            glyphs = glyph_source.glyph(
                scale="radii",
                geom=pv.Sphere(
                    radius=1.0, theta_resolution=resolution, phi_resolution=resolution
                ),
                orient=False,
            )

            if is_lighting_enabled:
                self.atom_actor = self.plotter.add_mesh(
                    glyphs, scalars="colors", rgb=True, **mesh_props
                )
            else:
                self.atom_actor = self.plotter.add_mesh(
                    glyphs,
                    scalars="colors",
                    rgb=True,
                    style="surface",
                    show_edges=True,
                    edge_color="grey",
                    **mesh_props,
                )
                self.atom_actor.GetProperty().SetEdgeOpacity(0.3)

            # Record atom color info
            for i, atom_color in enumerate(col):
                atom_rgb = [int(c * 255) for c in atom_color]
                self._3d_color_map[f"atom_{i}"] = atom_rgb

    def _add_3d_bond_cylinders(self, mol_to_draw, conf, col, current_style, mesh_props):
        # Draw bonds (ball_and_stick, wireframe, stick)
        if current_style in ["ball_and_stick", "wireframe", "stick"]:
            # Set bond radius and resolution based on style
            if current_style == "wireframe":
                cyl_radius = self.settings.get("wireframe_bond_radius", 0.01)
                bond_resolution = self.settings.get("wireframe_resolution", 6)
            elif current_style == "stick":
                cyl_radius = self.settings.get("stick_bond_radius", 0.15)
                bond_resolution = self.settings.get("stick_resolution", 16)
            else:  # ball_and_stick
                cyl_radius = self.settings.get("ball_stick_bond_radius", 0.1)
                bond_resolution = self.settings.get("ball_stick_resolution", 16)

            # Common color for Ball and Stick
            bs_bond_rgb = [127, 127, 127]
            if current_style == "ball_and_stick":
                try:
                    bs_hex = self.settings.get("ball_stick_bond_color", "#7F7F7F")
                    q = QColor(bs_hex)
                    bs_bond_rgb = [q.red(), q.green(), q.blue()]
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    # Suppress non-critical bond drawing color fallback noise
                    pass

            # Lists for batch processing
            all_points = []
            all_lines = []
            all_radii = []
            all_colors = []  # Cell data (one per line segment)

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
                if h == 0:
                    continue

                # Bond colors
                begin_color = col[begin_atom_idx]
                end_color = col[end_atom_idx]
                begin_color_rgb = [int(c * 255) for c in begin_color]
                end_color_rgb = [int(c * 255) for c in end_color]

                # Check for plugin override
                bond_idx = bond.GetIdx()
                # Override handling: if set, force both ends and uniform color to this value
                if (
                    hasattr(self, "_plugin_bond_color_overrides")
                    and bond_idx in self._plugin_bond_color_overrides
                ):
                    try:
                        # Expecting hex string
                        hex_c = self._plugin_bond_color_overrides[bond_idx]
                        c_obj = QColor(hex_c)
                        ov_rgb = [c_obj.red(), c_obj.green(), c_obj.blue()]
                        begin_color_rgb = ov_rgb
                        end_color_rgb = ov_rgb
                    except (
                        AttributeError,
                        RuntimeError,
                        TypeError,
                        ValueError,
                        KeyError,
                    ):
                        # Suppress traceback
                        pass

                # Determine effective uniform color for this bond
                local_bs_bond_rgb = (
                    begin_color_rgb
                    if (
                        hasattr(self, "_plugin_bond_color_overrides")
                        and bond_idx in self._plugin_bond_color_overrides
                    )
                    else bs_bond_rgb
                )

                # Helper to add segments
                def add_segment(p1, p2, radius, color_rgb):
                    nonlocal current_point_idx
                    all_points.append(p1)
                    all_points.append(p2)
                    all_lines.append([2, current_point_idx, current_point_idx + 1])
                    all_radii.append(radius)
                    all_radii.append(radius)
                    all_colors.append(color_rgb)
                    current_point_idx += 2

                # Get CPK bond color setting once for all bond types
                use_cpk_bond = self.settings.get("ball_stick_use_cpk_bond_color", False)
                # If overwritten, treat as if we want to show that color (effectively behave like CPK_Split but with same color, or Uniform).
                # To be robust, if overwritten, we can force "use_cpk_bond" logic but with our same colors?
                # Actually, if overridden, we probably want the whole bond to be that color.

                is_overridden = (
                    hasattr(self, "_plugin_bond_color_overrides")
                    and bond_idx in self._plugin_bond_color_overrides
                )

                if (
                    bt == Chem.rdchem.BondType.SINGLE
                    or bt == Chem.rdchem.BondType.AROMATIC
                ):
                    if (
                        current_style == "ball_and_stick"
                        and not use_cpk_bond
                        and not is_overridden
                    ):
                        # Single segment (Uniform color)
                        add_segment(sp, ep, cyl_radius, local_bs_bond_rgb)
                        self._3d_color_map[f"bond_{bond_counter}"] = local_bs_bond_rgb
                    else:
                        # Split segments (CPK split colors OR Overridden uniform)
                        # If overridden, begin/end are same, so this produces a uniform looking bond split in middle
                        mid_point = (sp + ep) / 2
                        add_segment(sp, mid_point, cyl_radius, begin_color_rgb)
                        add_segment(mid_point, ep, cyl_radius, end_color_rgb)
                        self._3d_color_map[f"bond_{bond_counter}_start"] = (
                            begin_color_rgb
                        )
                        self._3d_color_map[f"bond_{bond_counter}_end"] = end_color_rgb

                else:
                    # Calculate multiple bond parameters
                    v1 = d / h
                    # Apply radius factor per model
                    if current_style == "ball_and_stick":
                        double_radius_factor = self.settings.get(
                            "ball_stick_double_bond_radius_factor", 0.8
                        )
                        triple_radius_factor = self.settings.get(
                            "ball_stick_triple_bond_radius_factor", 0.75
                        )
                    elif current_style == "wireframe":
                        double_radius_factor = self.settings.get(
                            "wireframe_double_bond_radius_factor", 0.8
                        )
                        triple_radius_factor = self.settings.get(
                            "wireframe_triple_bond_radius_factor", 0.75
                        )
                    elif current_style == "stick":
                        double_radius_factor = self.settings.get(
                            "stick_double_bond_radius_factor", 0.60
                        )
                        triple_radius_factor = self.settings.get(
                            "stick_triple_bond_radius_factor", 0.40
                        )
                    else:
                        double_radius_factor = 1.0
                        triple_radius_factor = 0.75

                    # Get offset factor from settings (per model)
                    if current_style == "ball_and_stick":
                        double_offset_factor = self.settings.get(
                            "ball_stick_double_bond_offset_factor", 2.0
                        )
                        triple_offset_factor = self.settings.get(
                            "ball_stick_triple_bond_offset_factor", 2.0
                        )
                    elif current_style == "wireframe":
                        double_offset_factor = self.settings.get(
                            "wireframe_double_bond_offset_factor", 3.0
                        )
                        triple_offset_factor = self.settings.get(
                            "wireframe_triple_bond_offset_factor", 3.0
                        )
                    elif current_style == "stick":
                        double_offset_factor = self.settings.get(
                            "stick_double_bond_offset_factor", 1.5
                        )
                        triple_offset_factor = self.settings.get(
                            "stick_triple_bond_offset_factor", 1.0
                        )
                    else:
                        double_offset_factor = 2.0
                        triple_offset_factor = 2.0

                    if bt == Chem.rdchem.BondType.DOUBLE:
                        r = cyl_radius * double_radius_factor
                        off_dir = self._calculate_double_bond_offset(
                            mol_to_draw, bond, conf
                        )
                        s_double = cyl_radius * double_offset_factor

                        p1_start = sp + off_dir * (s_double / 2)
                        p1_end = ep + off_dir * (s_double / 2)
                        p2_start = sp - off_dir * (s_double / 2)
                        p2_end = ep - off_dir * (s_double / 2)

                        if (
                            current_style == "ball_and_stick"
                            and not use_cpk_bond
                            and not is_overridden
                        ):
                            add_segment(p1_start, p1_end, r, local_bs_bond_rgb)
                            add_segment(p2_start, p2_end, r, local_bs_bond_rgb)
                            self._3d_color_map[f"bond_{bond_counter}_1"] = (
                                local_bs_bond_rgb
                            )
                            self._3d_color_map[f"bond_{bond_counter}_2"] = (
                                local_bs_bond_rgb
                            )
                        else:
                            mid1 = (p1_start + p1_end) / 2
                            mid2 = (p2_start + p2_end) / 2
                            add_segment(p1_start, mid1, r, begin_color_rgb)
                            add_segment(mid1, p1_end, r, end_color_rgb)
                            add_segment(p2_start, mid2, r, begin_color_rgb)
                            add_segment(mid2, p2_end, r, end_color_rgb)
                            self._3d_color_map[f"bond_{bond_counter}_1_start"] = (
                                begin_color_rgb
                            )
                            self._3d_color_map[f"bond_{bond_counter}_1_end"] = (
                                end_color_rgb
                            )
                            self._3d_color_map[f"bond_{bond_counter}_2_start"] = (
                                begin_color_rgb
                            )
                            self._3d_color_map[f"bond_{bond_counter}_2_end"] = (
                                end_color_rgb
                            )

                    elif bt == Chem.rdchem.BondType.TRIPLE:
                        r = cyl_radius * triple_radius_factor
                        v_arb = np.array([0, 0, 1])
                        if np.allclose(np.abs(np.dot(v1, v_arb)), 1.0):
                            v_arb = np.array([0, 1, 0])
                        off_dir = np.cross(v1, v_arb)
                        off_dir /= np.linalg.norm(off_dir)
                        s_triple = cyl_radius * triple_offset_factor

                        # Center
                        if (
                            current_style == "ball_and_stick"
                            and not use_cpk_bond
                            and not is_overridden
                        ):
                            add_segment(sp, ep, r, local_bs_bond_rgb)
                            self._3d_color_map[f"bond_{bond_counter}_1"] = (
                                local_bs_bond_rgb
                            )
                        else:
                            mid = (sp + ep) / 2
                            add_segment(sp, mid, r, begin_color_rgb)
                            add_segment(mid, ep, r, end_color_rgb)
                            self._3d_color_map[f"bond_{bond_counter}_1_start"] = (
                                begin_color_rgb
                            )
                            self._3d_color_map[f"bond_{bond_counter}_1_end"] = (
                                end_color_rgb
                            )

                        # Sides
                        for sign in [1, -1]:
                            offset = off_dir * s_triple * sign
                            p_start = sp + offset
                            p_end = ep + offset

                            if (
                                current_style == "ball_and_stick"
                                and not use_cpk_bond
                                and not is_overridden
                            ):
                                add_segment(p_start, p_end, r, local_bs_bond_rgb)
                                suffix = "_2" if sign == 1 else "_3"
                                self._3d_color_map[f"bond_{bond_counter}{suffix}"] = (
                                    local_bs_bond_rgb
                                )
                            else:
                                mid = (p_start + p_end) / 2
                                add_segment(p_start, mid, r, begin_color_rgb)
                                add_segment(mid, p_end, r, end_color_rgb)
                                suffix = "_2" if sign == 1 else "_3"
                                self._3d_color_map[
                                    f"bond_{bond_counter}{suffix}_start"
                                ] = begin_color_rgb
                                self._3d_color_map[
                                    f"bond_{bond_counter}{suffix}_end"
                                ] = end_color_rgb

                bond_counter += 1

            # Generate and draw geometry
            if all_points:
                # Create PolyData
                bond_pd = pv.PolyData(np.array(all_points), lines=np.hstack(all_lines))
                # lines needs to be a flat array with padding indicating number of points per cell
                # all_lines is [[2, i, j], [2, k, l], ...], flatten it

                # Add data
                bond_pd.point_data["radii"] = np.array(all_radii)

                # Convert colors to 0-1 range for PyVista if needed, but add_mesh with rgb=True expects uint8 if using direct array?
                # Actually pyvista scalars usually prefer float 0-1 or uint8 0-255.
                # Let's use uint8 0-255 and rgb=True.
                bond_pd.cell_data["colors"] = np.array(all_colors, dtype=np.uint8)

                # Tube filter
                # n_sides corresponds to theta_resolution in Cylinder
                tube = bond_pd.tube(
                    scalars="radii",
                    absolute=True,
                    radius_factor=1.0,
                    n_sides=bond_resolution,
                    capping=True,
                )

                # Add to plotter
                self.plotter.add_mesh(tube, scalars="colors", rgb=True, **mesh_props)

    def _add_3d_aromatic_rings(self, mol_to_draw, current_style, mesh_props):
        # Aromatic ring circles display
        if self.settings.get("display_aromatic_circles_3d", False):
            try:
                ring_info = mol_to_draw.GetRingInfo()
                aromatic_rings = []

                # Find aromatic rings
                for ring in ring_info.AtomRings():
                    # Check if all atoms in ring are aromatic
                    is_aromatic = all(
                        mol_to_draw.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring
                    )
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
                    distances = [
                        np.linalg.norm(pos - center) for pos in ring_positions_np
                    ]
                    ring_radius = np.mean(distances) * 0.55  # Slightly smaller

                    # Get bond radius from current style settings for torus thickness
                    if current_style == "stick":
                        bond_radius = self.settings.get("stick_bond_radius", 0.15)
                    elif current_style == "ball_and_stick":
                        bond_radius = self.settings.get("ball_stick_bond_radius", 0.1)
                    elif current_style == "wireframe":
                        bond_radius = self.settings.get("wireframe_bond_radius", 0.01)
                    else:
                        bond_radius = 0.1  # Default
                    # Apply user-defined thickness factor (default 0.6)
                    thickness_factor = self.settings.get(
                        "aromatic_torus_thickness_factor", 0.6
                    )
                    tube_radius = bond_radius * thickness_factor
                    theta = np.linspace(0, 2.2 * np.pi, 64)
                    circle_x = ring_radius * np.cos(theta)
                    circle_y = ring_radius * np.sin(theta)
                    circle_z = np.zeros_like(theta)
                    circle_points = np.c_[circle_x, circle_y, circle_z]

                    # Create line from points
                    circle_line = pv.Spline(circle_points, n_points=64).tube(
                        radius=tube_radius, n_sides=16
                    )

                    # Rotate torus to align with ring plane
                    # Default torus is in XY plane (normal = [0, 0, 1])
                    default_normal = np.array([0, 0, 1])

                    # Calculate rotation axis and angle
                    if not np.allclose(normal, default_normal) and not np.allclose(
                        normal, -default_normal
                    ):
                        axis = np.cross(default_normal, normal)
                        axis_length = np.linalg.norm(axis)
                        if axis_length > 0:
                            axis = axis / axis_length
                            angle = np.arccos(
                                np.clip(np.dot(default_normal, normal), -1.0, 1.0)
                            )
                            angle_deg = np.degrees(angle)

                            # Rotate torus
                            circle_line = circle_line.rotate_vector(
                                axis, angle_deg, point=[0, 0, 0]
                            )

                    # Translate to ring center
                    circle_line = circle_line.translate(center)

                    # Get torus color from bond color settings
                    # Calculate most common atom type in ring for CPK color
                    from collections import Counter

                    atom_symbols = [
                        mol_to_draw.GetAtomWithIdx(idx).GetSymbol() for idx in ring
                    ]
                    most_common_symbol = (
                        Counter(atom_symbols).most_common(1)[0][0]
                        if atom_symbols
                        else None
                    )

                    if current_style == "ball_and_stick":
                        # Check if using CPK bond colors
                        use_cpk = self.settings.get(
                            "ball_stick_use_cpk_bond_color", False
                        )
                        if use_cpk:
                            # Use CPK color of most common atom type in ring
                            if most_common_symbol:
                                cpk_color = CPK_COLORS_PV.get(
                                    most_common_symbol, [0.5, 0.5, 0.5]
                                )
                                torus_color = cpk_color
                            else:
                                torus_color = [0.5, 0.5, 0.5]
                        else:
                            # Use Ball & Stick bond color setting
                            bond_hex = self.settings.get(
                                "ball_stick_bond_color", "#7F7F7F"
                            )
                            q = QColor(bond_hex)
                            torus_color = [
                                q.red() / 255.0,
                                q.green() / 255.0,
                                q.blue() / 255.0,
                            ]
                    else:
                        # For Wireframe and Stick, use CPK color of most common atom
                        if most_common_symbol:
                            cpk_color = CPK_COLORS_PV.get(
                                most_common_symbol, [0.5, 0.5, 0.5]
                            )
                            torus_color = cpk_color
                        else:
                            torus_color = [0.5, 0.5, 0.5]

                    self.plotter.add_mesh(circle_line, color=torus_color, **mesh_props)

            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                logging.error(f"Error rendering aromatic circles: {e}")

    def _add_3d_labels(self, mol, mol_to_draw):
        if getattr(self, "show_chiral_labels", False):
            try:
                # Calculate chiral centers from 3D coordinates
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if chiral_centers:
                    pts, labels = [], []
                    z_off = 0
                    for idx, lbl in chiral_centers:
                        coord = self.atom_positions_3d[idx].copy()
                        coord[2] += z_off
                        pts.append(coord)
                        labels.append(lbl if lbl is not None else "?")
                    try:
                        self.plotter.remove_actor("chiral_labels")
                    except (AttributeError, RuntimeError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress non-critical 3D label update errors
                    self.plotter.add_point_labels(
                        np.array(pts),
                        labels,
                        font_size=20,
                        point_size=0,
                        text_color="blue",
                        name="chiral_labels",
                        always_visible=True,
                        tolerance=0.01,
                        show_points=False,
                    )
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                self.statusBar().showMessage(f"3D chiral label drawing error: {e}")

        # Also display E/Z labels
        if getattr(self, "show_chiral_labels", False):
            try:
                # If we drew a kekulized molecule use it for E/Z detection so
                # E/Z labels reflect Kekulﾃｩ rendering; pass mol_to_draw as the
                # molecule to scan for bond stereochemistry.
                self.show_ez_labels_3d(mol)
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                self.statusBar().showMessage(f"3D E/Z label drawing error: {e}")

    def _calculate_double_bond_offset(self, mol, bond, conf):
        """
        Calculate double bond offset direction.
        Consider other bonds of connected atoms to keep it planar.
        """
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())

        begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
        end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))

        bond_vec = end_pos - begin_pos
        bond_length = np.linalg.norm(bond_vec)
        if bond_length == 0:
            # Fallback: Z-axis reference
            return np.array([0, 0, 1])

        bond_unit = bond_vec / bond_length

        # Check neighbors of both terminal atoms
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

        # Calculate plane normal vector
        normal_candidates = []

        # Estimate plane from start atom's neighbors
        if len(begin_neighbors) >= 1:
            for neighbor_pos in begin_neighbors:
                vec_to_neighbor = neighbor_pos - begin_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # Normal to the plane is the cross product of bond_vec and neighbor_vec
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)

        # Estimate plane from neighbors of the end atom
        if len(end_neighbors) >= 1:
            for neighbor_pos in end_neighbors:
                vec_to_neighbor = neighbor_pos - end_pos
                if np.linalg.norm(vec_to_neighbor) > 1e-6:
                    # Cross product of bond_vec and neighbor_vec is the plane normal
                    normal = np.cross(bond_vec, vec_to_neighbor)
                    norm_length = np.linalg.norm(normal)
                    if norm_length > 1e-6:
                        normal_candidates.append(normal / norm_length)

        # If multiple normal vectors exist, take the average
        if normal_candidates:
            # Adjust so the dot product with the first vector is positive for consistent direction
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

                # Set offset direction perpendicular to normal and bond vectors
                offset_dir = np.cross(bond_unit, avg_normal)
                offset_length = np.linalg.norm(offset_dir)
                if offset_length > 1e-6:
                    return offset_dir / offset_length

        # Fallback: Arbitrary direction perpendicular to the bond vector
        v_arb = np.array([0, 0, 1])
        if np.allclose(np.abs(np.dot(bond_unit, v_arb)), 1.0):
            v_arb = np.array([0, 1, 0])

        off_dir = np.cross(bond_unit, v_arb)
        off_dir /= np.linalg.norm(off_dir)
        return off_dir

    def show_ez_labels_3d(self, mol):
        """Display E/Z labels in 3D view (using RDKit stereochemistry determination)"""
        if not mol:
            return

        # Remove existing E/Z labels
        if (
            hasattr(self.plotter, "renderer")
            and "ez_labels" in self.plotter.renderer.actors
        ):
            try:
                self.plotter.remove_actor("ez_labels")
            except (AttributeError, RuntimeError, TypeError):
                # Ignore label removal failure on stale plotter
                pass

        pts, labels = [], []

        # Check if 3D coordinates exist
        if mol.GetNumConformers() == 0:
            return

        conf = mol.GetConformer()

        # Display E/Z stereochemistry determined by RDKit for double bonds

        try:
            # Recalculate stereochemistry from 3D coordinates (on mol).
            # This ensures determination is based on actual 3D positions regardless of 2D state.
            Chem.AssignStereochemistry(
                mol, cleanIt=True, force=True, flagPossibleStereoCenters=True
            )
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Suppress non-critical stereochemistry assignment noise during 3D label update
            pass

        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                new_stereo = bond.GetStereo()

                if new_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                    # Calculate bond center coordinates
                    begin_pos = np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx()))
                    end_pos = np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
                    center_pos = (begin_pos + end_pos) / 2

                    # Determine label
                    label = "E" if new_stereo == Chem.BondStereo.STEREOE else "Z"

                    # Check for discrepancy with 2D
                    # Get 2D-derived stereochemistry property saved in main_window_compute.py
                    try:
                        old_stereo = bond.GetIntProp("_original_2d_stereo")
                    except KeyError:
                        old_stereo = Chem.BondStereo.STEREONONE

                    # Set to "?" if 2D also has E/Z specified but it differs from 3D
                    if old_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                        if old_stereo != new_stereo:
                            label = "?"

                    pts.append(center_pos)
                    labels.append(label)

        if pts and labels:
            self.plotter.add_point_labels(
                np.array(pts),
                labels,
                font_size=18,
                point_size=0,
                text_color="darkgreen",  # Dark green color
                name="ez_labels",
                always_visible=True,
                tolerance=0.01,
                show_points=False,
            )

    def toggle_chiral_labels_display(self, checked):
        """Toggle chiral label display based on View menu action"""
        self.show_chiral_labels = checked

        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)

        if checked:
            self.statusBar().showMessage(
                "Chiral labels: will be (re)computed after Convert→3D."
            )
        else:
            self.statusBar().showMessage("Chiral labels disabled.")

    def update_chiral_labels(self):
        """Calculate chiral centers and set/clear R/S labels on 2D AtomItems.
        Prefer 3D (self.current_mol) if available; otherwise use RDKit mol from 2D.
        """
        # First clear labels from all items
        for atom_data in self.data.atoms.values():
            if atom_data.get("item"):
                atom_data["item"].chiral_label = None

        if not self.show_chiral_labels:
            self.scene.update()
            return

        # Use 3D RDKit Mol with conformer
        mol_for_chirality = None
        if getattr(self, "current_mol", None) is not None:
            mol_for_chirality = self.current_mol
        else:
            return

        if mol_for_chirality is None or mol_for_chirality.GetNumAtoms() == 0:
            self.scene.update()
            return

        try:
            # --- IMPORTANT: if 3D conformer exists, use it to assign chiral tags ---
            if mol_for_chirality.GetNumConformers() > 0:
                # Set chiral tags from 3D coordinates using confId=0
                try:
                    Chem.AssignAtomChiralTagsFromStructure(mol_for_chirality, confId=0)
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    # Guard for older RDKit versions or invalid coordinates
                    pass

            # Get chiral centers (list of (idx, 'R'/'S'/'?'))
            chiral_centers = Chem.FindMolChiralCenters(
                mol_for_chirality, includeUnassigned=True
            )

            # Mapping from RDKit atom index to editor atom_id
            rdkit_idx_to_my_id = {}
            for atom in mol_for_chirality.GetAtoms():
                if atom.HasProp("_original_atom_id"):
                    rdkit_idx_to_my_id[atom.GetIdx()] = atom.GetIntProp(
                        "_original_atom_id"
                    )

            # Set found chiral centers to corresponding AtomItems
            for idx, label in chiral_centers:
                if idx in rdkit_idx_to_my_id:
                    atom_id = rdkit_idx_to_my_id[idx]
                    if atom_id in self.data.atoms and self.data.atoms[atom_id].get(
                        "item"
                    ):
                        # 'R' / 'S' / '?'
                        self.data.atoms[atom_id]["item"].chiral_label = label

        except (AttributeError, RuntimeError, TypeError, ValueError) as e:
            self.statusBar().showMessage(f"Update chiral labels error: {e}")

        # Finally redraw 2D scene
        self.scene.update()

    def toggle_atom_info_display(self, mode):
        """Toggle atom info display mode"""
        # Clear current display
        self.clear_all_atom_info_labels()

        # Turn OFF if the same mode is selected
        if self.atom_info_display_mode == mode:
            self.atom_info_display_mode = None
            # Uncheck all actions
            self.show_atom_id_action.setChecked(False)
            self.show_rdkit_id_action.setChecked(False)
            self.show_atom_coords_action.setChecked(False)
            self.show_atom_symbol_action.setChecked(False)
            self.statusBar().showMessage("Atom info display disabled.")
        else:
            # Set new mode
            self.atom_info_display_mode = mode
            # Check only the relevant action
            self.show_atom_id_action.setChecked(mode == "id")
            self.show_rdkit_id_action.setChecked(mode == "rdkit_id")
            self.show_atom_coords_action.setChecked(mode == "coords")
            self.show_atom_symbol_action.setChecked(mode == "symbol")

            mode_names = {
                "id": "Atom ID",
                "rdkit_id": "RDKit Index",
                "coords": "Coordinates",
                "symbol": "Element Symbol",
            }
            self.statusBar().showMessage(f"Displaying: {mode_names[mode]}")

            # Display info for all atoms
            self.show_all_atom_info()

    def is_xyz_derived_molecule(self):
        """Determine if the current molecule is derived from an XYZ file"""
        if not self.current_mol:
            return False
        if not self.current_mol or self.current_mol.GetNumAtoms() == 0:
            return False

        try:
            # Check if the first atom has xyz_unique_id property
            return self.current_mol.GetAtomWithIdx(0).HasProp("xyz_unique_id")
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Suppress non-critical property access noise
            return False

        return False

    def has_original_atom_ids(self):
        """Determine if the current molecule has Original Atom IDs"""
        if not self.current_mol:
            return False
        if not self.current_mol:
            return False
        try:
            # Check if any atom has _original_atom_id property
            for atom_idx in range(self.current_mol.GetNumAtoms()):
                atom = self.current_mol.GetAtomWithIdx(atom_idx)
                if atom.HasProp("_original_atom_id"):
                    return True
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Suppress non-critical property access noise
            pass
        return False

        return False

    def update_atom_id_menu_text(self):
        """Update Atom ID menu text based on molecule type"""
        if hasattr(self, "show_atom_id_action"):
            if self.is_xyz_derived_molecule():
                self.show_atom_id_action.setText("Show XYZ Unique ID")
            else:
                self.show_atom_id_action.setText("Show Original ID / Index")

    def update_atom_id_menu_state(self):
        """Update Atom ID menu enabled/disabled state"""
        if hasattr(self, "show_atom_id_action"):
            has_original_ids = self.has_original_atom_ids()
            has_xyz_ids = self.is_xyz_derived_molecule()

            # Enable only if Original ID or XYZ ID exists
            self.show_atom_id_action.setEnabled(has_original_ids or has_xyz_ids)

            # Disable selection if the currently selected mode is no longer valid
            if (
                not (has_original_ids or has_xyz_ids)
                and self.atom_info_display_mode == "id"
            ):
                self.atom_info_display_mode = None
                self.show_atom_id_action.setChecked(False)
                self.clear_all_atom_info_labels()

    def show_all_atom_info(self):
        """Display info for all atoms"""
        if (
            self.atom_info_display_mode is None
            or not hasattr(self, "atom_positions_3d")
            or self.atom_positions_3d is None
        ):
            return

        # Clear existing labels
        self.clear_all_atom_info_labels()

        # Group by type to create lists for label display
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

            if self.atom_info_display_mode == "id":
                # Display Original ID if available, otherwise XYZ unique ID, finally RDKit index
                try:
                    if self.current_mol:
                        atom = self.current_mol.GetAtomWithIdx(atom_idx)
                        if atom.HasProp("_original_atom_id"):
                            original_id = atom.GetIntProp("_original_atom_id")
                            # Remove prefix and display only the number
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
                except (AttributeError, RuntimeError, TypeError, ValueError):
                    rdkit_positions.append(pos)
                    rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == "rdkit_id":
                rdkit_positions.append(pos)
                rdkit_texts.append(str(atom_idx))

            elif self.atom_info_display_mode == "coords":
                other_positions.append(pos)
                other_texts.append(f"({pos[0]:.2f},{pos[1]:.2f},{pos[2]:.2f})")

            elif self.atom_info_display_mode == "symbol":
                if self.current_mol:
                    symbol = self.current_mol.GetAtomWithIdx(atom_idx).GetSymbol()
                    other_positions.append(pos)
                    other_texts.append(symbol)
                else:
                    other_positions.append(pos)
                    other_texts.append("?")

            else:
                continue

        # Color definitions (dark blue/green/red)
        rdkit_color = "#003366"  # Dark blue
        id_color = "#006400"  # Dark green
        xyz_color = "#8B0000"  # Dark red
        other_color = "black"

        # Add labels for each group and keep references in a list
        self.current_atom_info_labels = []
        try:
            if rdkit_positions:
                a = self.plotter.add_point_labels(
                    np.array(rdkit_positions),
                    rdkit_texts,
                    point_size=12,
                    font_size=18,
                    text_color=rdkit_color,
                    always_visible=True,
                    tolerance=0.01,
                    show_points=False,
                    name="atom_labels_rdkit",
                )
                self.current_atom_info_labels.append(a)

            if id_positions:
                a = self.plotter.add_point_labels(
                    np.array(id_positions),
                    id_texts,
                    point_size=12,
                    font_size=18,
                    text_color=id_color,
                    always_visible=True,
                    tolerance=0.01,
                    show_points=False,
                    name="atom_labels_id",
                )
                self.current_atom_info_labels.append(a)

            if xyz_positions:
                a = self.plotter.add_point_labels(
                    np.array(xyz_positions),
                    xyz_texts,
                    point_size=12,
                    font_size=18,
                    text_color=xyz_color,
                    always_visible=True,
                    tolerance=0.01,
                    show_points=False,
                    name="atom_labels_xyz",
                )
                self.current_atom_info_labels.append(a)

            if other_positions:
                a = self.plotter.add_point_labels(
                    np.array(other_positions),
                    other_texts,
                    point_size=12,
                    font_size=18,
                    text_color=other_color,
                    always_visible=True,
                    tolerance=0.01,
                    show_points=False,
                    name="atom_labels_other",
                )
                self.current_atom_info_labels.append(a)
        except (AttributeError, RuntimeError, TypeError, ValueError) as e:
            print(f"Error adding atom info labels: {e}")

        # Display legend in the top-right (remove existing legends)
        try:
            if (
                hasattr(self, "atom_label_legend_names")
                and self.atom_label_legend_names
            ):
                for nm in self.atom_label_legend_names:
                    try:
                        self.plotter.remove_actor(nm)
                    except (AttributeError, RuntimeError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress non-critical 3D rendering/actor update errors
            self.atom_label_legend_names = []

            legend_entries = []
            if rdkit_positions:
                legend_entries.append(("RDKit", rdkit_color, "legend_rdkit"))
            if id_positions:
                legend_entries.append(("ID", id_color, "legend_id"))
            if xyz_positions:
                legend_entries.append(("XYZ", xyz_color, "legend_xyz"))

            # Add legend labels at bottom-left (no background, bold only)
            # Increase spacing to avoid overlapping when short labels like 'RDKit' and 'ID' appear
            spacing = 30
            for i, (label_text, label_color, label_name) in enumerate(legend_entries):
                # Incremental y-coordinate from bottom-left
                # Add a small horizontal offset for very short adjacent labels so they don't visually collide
                y = 0.0 + i * spacing
                x_offset = 0.0
                # If both RDKit and ID are present, nudge the second entry slightly to the right to avoid overlap
                try:
                    if label_text == "ID" and any(
                        e[0] == "RDKit" for e in legend_entries
                    ):
                        x_offset = 0.06
                except (AttributeError, RuntimeError, TypeError):
                    x_offset = 0.0
                try:
                    actor = self.plotter.add_text(
                        label_text,
                        position=(0.0 + x_offset, y),
                        font_size=12,
                        color=label_color,
                        name=label_name,
                        font="arial",
                    )
                    self.atom_label_legend_names.append(label_name)
                    # Set bold only (no background)
                    try:
                        if hasattr(actor, "GetTextProperty"):
                            tp = actor.GetTextProperty()
                            try:
                                tp.SetBold(True)
                            except (AttributeError, RuntimeError, TypeError):
                                # Suppress non-critical actor property update noise
                                pass
                    except (AttributeError, RuntimeError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress non-critical 3D rendering/actor update errors
                except (AttributeError, RuntimeError, TypeError):
                    continue

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress legend addition errors

    def clear_all_atom_info_labels(self):
        """Clear all atom info labels"""
        # Remove label actors (may be a single actor, a list, or None)
        try:
            if (
                hasattr(self, "current_atom_info_labels")
                and self.current_atom_info_labels
            ):
                if isinstance(self.current_atom_info_labels, (list, tuple)):
                    for a in list(self.current_atom_info_labels):
                        try:
                            self.plotter.remove_actor(a)
                        except (AttributeError, RuntimeError, TypeError):
                            # Suppress non-critical actor removal noise
                            pass
                else:
                    try:
                        self.plotter.remove_actor(self.current_atom_info_labels)
                    except (AttributeError, RuntimeError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress non-critical 3D rendering/actor update errors
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical 3D state update errors
        finally:
            self.current_atom_info_labels = None

        # Remove legend text actors if present
        try:
            if (
                hasattr(self, "atom_label_legend_names")
                and self.atom_label_legend_names
            ):
                for nm in list(self.atom_label_legend_names):
                    try:
                        self.plotter.remove_actor(nm)
                    except (AttributeError, RuntimeError, TypeError) as e:
                        logging.debug(
                            f"Suppressed exception: {e}"
                        )  # Suppress non-critical 3D rendering/actor update errors
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical 3D state update errors
        finally:
            self.atom_label_legend_names = []

    def setup_3d_hover(self):
        """Configure 3D view display (changed to persistent)"""
        if self.atom_info_display_mode is not None:
            self.show_all_atom_info()

    def zoom_in(self):
        """Zoom in by 20%"""
        self.view_2d.scale(1.2, 1.2)

    def zoom_out(self):
        """Zoom out by 20%"""
        self.view_2d.scale(1 / 1.2, 1 / 1.2)

    def reset_zoom(self):
        """Reset zoom to default (75%)"""
        transform = QTransform()
        transform.scale(0.75, 0.75)
        self.view_2d.setTransform(transform)

    def fit_to_view(self):
        """Fit all items in the scene into the view"""
        if not self.scene.items():
            self.reset_zoom()
            return

        # Calculate total bounding rect of visible items
        visible_items_rect = QRectF()
        for item in self.scene.items():
            if item.isVisible() and not isinstance(item, TemplatePreviewItem):
                if visible_items_rect.isEmpty():
                    visible_items_rect = item.sceneBoundingRect()
                else:
                    visible_items_rect = visible_items_rect.united(
                        item.sceneBoundingRect()
                    )

        if visible_items_rect.isEmpty():
            self.reset_zoom()
            return

        # Add some padding
        padding_factor = 1.10  # 10% margin
        cx = visible_items_rect.center().x()
        cy = visible_items_rect.center().y()
        w = visible_items_rect.width() * padding_factor
        h = visible_items_rect.height() * padding_factor
        padded = QRectF(cx - w / 2.0, cy - h / 2.0, w, h)

        # Temporarily set anchor to center before calling fitInView
        try:
            old_ta = self.view_2d.transformationAnchor()
            old_ra = self.view_2d.resizeAnchor()
        except (AttributeError, RuntimeError, TypeError):
            old_ta = old_ra = None

        try:
            self.view_2d.setTransformationAnchor(
                QGraphicsView.ViewportAnchor.AnchorViewCenter
            )
            self.view_2d.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
            self.view_2d.fitInView(padded, Qt.AspectRatioMode.KeepAspectRatio)
        finally:
            # Restore original anchor
            try:
                if old_ta is not None:
                    self.view_2d.setTransformationAnchor(old_ta)
                if old_ra is not None:
                    self.view_2d.setResizeAnchor(old_ra)
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress non-critical 3D view/actor cleanup errors

    def apply_3d_settings(self, redraw=True):
        """Apply 3D view visual settings"""
        # Projection mode
        proj_mode = self.settings.get("projection_mode", "Perspective")
        if hasattr(self.plotter, "renderer") and hasattr(
            self.plotter.renderer, "GetActiveCamera"
        ):
            cam = self.plotter.renderer.GetActiveCamera()
            if cam:
                if proj_mode == "Orthographic":
                    cam.SetParallelProjection(True)
                else:
                    cam.SetParallelProjection(False)

        if not hasattr(self, "plotter"):
            return

        # Enable renderer layers (for text overlay)
        renderer = self.plotter.renderer
        if renderer and hasattr(renderer, "SetNumberOfLayers"):
            try:
                renderer.SetNumberOfLayers(2)  # Layer 0: 3D, Layer 1: 2D Overlay
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # May not be supported depending on PyVista version

        # --- 3D Axis Widget Settings ---
        show_axes = self.settings.get("show_3d_axes", True)

        # Create widget if not already created
        if getattr(self, "axes_widget", None) is None and hasattr(
            self.plotter, "interactor"
        ):
            axes = vtk.vtkAxesActor()
            self.axes_widget = vtk.vtkOrientationMarkerWidget()
            self.axes_widget.SetOrientationMarker(axes)
            self.axes_widget.SetInteractor(self.plotter.interactor)
            # Position at bottom-left corner (20% width/height)
            self.axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)

        # Enable/disable widget based on settings
        if self.axes_widget:
            if show_axes:
                self.axes_widget.On()
                self.axes_widget.SetInteractive(False)
            else:
                self.axes_widget.Off()

        if redraw:
            self.draw_molecule_3d(self.current_mol)

        # Do not reset camera on settings change (reset only once)
        if not getattr(self, "_camera_initialized", False):
            try:
                self.plotter.reset_camera()
            except (AttributeError, RuntimeError, TypeError) as e:
                logging.debug(
                    f"Suppressed exception: {e}"
                )  # Suppress non-critical 3D view/actor cleanup errors

            self._camera_initialized = True

        # Force plotter update
        try:
            self.plotter.render()
            if hasattr(self.plotter, "update"):
                self.plotter.update()
        except (AttributeError, RuntimeError, TypeError) as e:
            logging.debug(
                f"Suppressed exception: {e}"
            )  # Suppress non-critical 3D state update errors

    def update_bond_color_override(self, bond_idx, hex_color):
        """Plugin API helper to override bond color."""
        if not hasattr(self, "_plugin_bond_color_overrides"):
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
        if not hasattr(self, "_plugin_color_overrides"):
            self._plugin_color_overrides = {}

        if color_hex is None:
            if atom_index in self._plugin_color_overrides:
                del self._plugin_color_overrides[atom_index]
        else:
            self._plugin_color_overrides[atom_index] = color_hex

        if self.current_mol:
            self.draw_molecule_3d(self.current_mol)


# Set class-level marker for plugin compatibility
MainWindowView3d._cls = MainWindowView3d
