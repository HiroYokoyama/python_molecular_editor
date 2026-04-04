#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

DEFAULT_SETTINGS = {
    # --- 3D Scene Visuals ---
    "background_color": "#919191",
    "projection_mode": "Perspective",
    "lighting_enabled": True,
    "specular": 0.20,
    "specular_power": 20,
    "light_intensity": 1.0,
    "show_3d_axes": True,
    "show_chiral_labels": False,
    # --- 3D Model Parameters (Ball and Stick) ---
    "ball_stick_atom_scale": 1.0,
    "ball_stick_bond_radius": 0.1,
    "ball_stick_resolution": 16,
    "ball_stick_bond_color": "#7F7F7F",
    "ball_stick_use_cpk_bond_color": False,
    "ball_stick_double_bond_offset_factor": 2.0,
    "ball_stick_triple_bond_offset_factor": 2.0,
    "ball_stick_double_bond_radius_factor": 0.8,
    "ball_stick_triple_bond_radius_factor": 0.75,
    # --- 3D Model Parameters (CPK) ---
    "cpk_atom_scale": 1.0,
    "cpk_resolution": 32,
    "cpk_colors": {},
    # --- 3D Model Parameters (Wireframe) ---
    "wireframe_bond_radius": 0.02,
    "wireframe_resolution": 6,
    "wireframe_double_bond_offset_factor": 3.0,
    "wireframe_triple_bond_offset_factor": 3.0,
    "wireframe_double_bond_radius_factor": 0.8,
    "wireframe_triple_bond_radius_factor": 0.75,
    # --- 3D Model Parameters (Stick) ---
    "stick_bond_radius": 0.15,
    "stick_resolution": 16,
    "stick_double_bond_offset_factor": 1.5,
    "stick_triple_bond_offset_factor": 1.0,
    "stick_double_bond_radius_factor": 0.6,
    "stick_triple_bond_radius_factor": 0.4,
    # --- Other 3D Rendering Options ---
    "aromatic_torus_thickness_factor": 0.6,
    "display_aromatic_circles_3d": False,
    "display_kekule_3d": False,
    # --- 3D Conversion and Optimization ---
    "3d_conversion_mode": "rdkit",
    "optimization_method": "MMFF_RDKIT",
    "optimize_intermolecular_interaction_rdkit": True,
    "skip_chemistry_checks": False,
    "always_ask_charge": False,
    # --- Interaction and Selection ---
    "use_high_fidelity_selection": True,
    "selection_color": "#FFD700",
    "high_quality_meshing": True,
    "atom_label_color": "#000000",
    "bond_color": "#808080",
    # --- 2D Settings ---
    "bond_width_2d": 2.0,
    "bond_spacing_double_2d": 3.5,
    "bond_spacing_triple_2d": 3.5,
    "atom_font_size_2d": 20,
    "background_color_2d": "#FFFFFF",
    "bond_color_2d": "#222222",
    "atom_use_bond_color_2d": False,
    "bond_cap_style_2d": "Round",
    "bond_wedge_width_2d": 6.0,
    "bond_dash_count_2d": 8,
    "atom_font_family_2d": "Arial",
    # --- Application Session / UI State ---
    "theme": "light",
    "window_size": [1200, 800],
    "window_position": [100, 100],
    "splitter_sizes": [600, 600],
    "last_dir": "",
}
