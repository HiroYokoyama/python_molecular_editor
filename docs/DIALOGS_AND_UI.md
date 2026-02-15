# Dialogs and UI Components Documentation

This document provides a detailed overview of the various dialogs and custom UI components used in **MoleditPy**.

## General Dialogs

### SettingsDialog (`settings_dialog.py`)
The central configuration hub for the application.
- **Tabs**:
    - **2D Settings**: Atom/Bond colors, fonts, bond width/spacing.
    - **3D Scene**: Background color, projection mode (Perspective/Orthographic), lighting, fog.
    - **Ball & Stick**: Sphere scale, bond radius, resolution.
    - **CPK / Space Filling**: Sphere scale, resolution.
    - **Wireframe / Stick**: Bond widths.
    - **Other**: Optimization defaults, startup preferences.
- **Key Features**: Updates global `settings` dictionary and triggers immediate redraws in both 2D and 3D views.

### AboutDialog (`about_dialog.py`)
Displays application version, author, license, and logo.

### PeriodicTableDialog (`periodic_table_dialog.py`)
A grid-based dialog for selecting chemical elements.
- **Usage**: Used for changing atom types or selecting elements for specific operations.

## Geometry & Manipulation Dialogs

### TranslationDialog (`translation_dialog.py`)
Allows precise translation of the molecule or selected atoms in 3D space.
- **Features**:
    - Select atoms by clicking in the 3D view.
    - Translate by specific dX, dY, dZ values.
    - "Center" functionality to move the molecule to the origin.

### MoveGroupDialog (`move_group_dialog.py`)
Advanced tool for moving and rotating substructures or connected groups.
- **Features**:
    - **Group Selection**: Click an atom to select its connected group (BFS/DFS). Ctrl+Click to add/remove groups.
    - **Interaction**: Drag atoms to move the group; Right-drag to rotate the group around its geometric center.
    - **Visuals**: Highlights selected atoms and displays labels.

### AlignPlaneDialog (`align_plane_dialog.py`)
Aligns a selected set of atoms (minimum 3) to a specific Cartesian plane (XY, XZ, or YZ).
- **Algorithm**: Uses PCA (Principal Component Analysis) to find the best-fit plane of the selected atoms and rotates the molecule so that this plane aligns with the target global plane.

### AlignmentDialog (`alignment_dialog.py`)
Aligns a specific bond axis to a global axis (X, Y, or Z).
- **Usage**: Requires exactly 2 selected atoms.
- **Effect**: Translates the first atom to the origin and rotates the molecule so the second atom lies on the specified axis.

### MirrorDialog (`mirror_dialog.py`)
Reflects the molecule across a selected plane (XY, XZ, YZ).
- **Usage**: Useful for generating enantiomers.

### PlanarizeDialog (`planarize_dialog.py`)
Flattens selected atoms onto their best-fit plane.
- **Usage**: Useful for correcting distorted aromatic rings or planar groups.

## Geometric Adjustment Dialogs

### BondLengthDialog (`bond_length_dialog.py`)
Adjusts the distance between two bonded atoms.
- **Modes**:
    - Move Atom 2 only (Atom 1 fixed).
    - Move Atom 2's connected group (Atom 1 fixed).
    - Move both groups towards/away from the bond center.

### AngleDialog (`angle_dialog.py`)
Adjusts the bond angle defined by three atoms.
- **Modes**:
    - Rotate Atom 3 only.
    - Rotate Atom 3's connected group.
    - Rotate both arms (Atom 1 group and Atom 3 group) equally.

### DihedralDialog (`dihedral_dialog.py`)
Adjusts the torsion (dihedral) angle defined by four atoms.
- **Modes**:
    - Rotate the group attached to Atom 4 (Atom 1-2-3 fixed).
    - Rotate the group attached to Atom 1 (Atom 2-3-4 fixed).

## Analysis & Calculation

### AnalysisWindow (`analysis_window.py`)
Displays detailed physicochemical properties of the current molecule.
- **Properties**: Molecular Formula, Weight, Exact Mass, LogP, TPSA, Ring Count, Rotatable Bonds, H-Bond Donors/Acceptors, SMILES, InChI.
- **Logic**: Handles XYZ-derived structures (coordinates only) differently from MOL-derived structures (with bond topology).

### ConstrainedOptimizationDialog (`constrained_optimization_dialog.py`)
Setup for performing geometry optimization with distance constraints.
- **Features**:
    - Select pairs of atoms to constrain.
    - Specify target distances.
    - Run optimization (MMFF/UFF) while respecting these constraints (via adding force field constraints).

### CalculationWorker (`calculation_worker.py`)
*Technically a background worker, but critical for UI responsiveness.*
- **Role**: Runs valid 3D embedding and optimization tasks in a separate thread.
- **Key Methods**: `run_conversion`, `run_optimization`.

## User Interface Components

### PluginManagerWindow (`plugin_manager_window.py`)
Manages external Python plugins.
- **Features**:
    - Drag & Drop installation of `.py` or `.zip` plugins.
    - Enable/Disable/Remove plugins.
    - Reload plugins dynamically.

### ColorSettingsDialog (`color_settings_dialog.py`)
Specialized dialog for customizing CPK element colors.
- **Features**:
    - Grid of elements to pick colors.
    - Persistent overrides saved to settings.
    - "Reset All" functionality.

### UserTemplateDialog (`user_template_dialog.py`)
Manages user-defined structural templates.
- **Features**:
    - Save current structure as a template.
    - Grid view of saved templates.
    - Apply template to current scene.

## Custom Graphics & Interaction

### MoleculeScene (`molecule_scene.py`)
The core 2D editing canvas (`QGraphicsScene` subclass).
- **Role**: Manages all 2D items (atoms, bonds) and handles mouse events for drawing, selecting, and modifying structure.

### ZoomableView (`zoomable_view.py`)
The widget (`QGraphicsView` subclass) that displays the `MoleculeScene`.
- **Features**: Handles zooming (wheel events) and panning interactions.

### AtomItem / BondItem (`atom_item.py`, `bond_item.py`)
Custom `QGraphicsItem` implementations for high-performance 2D rendering.
- **AtomItem**: Renders symbol, implicit hydrogens (with subscript support), charge, and radical state. Handles hit detection.
- **BondItem**: Renders lines for bond orders. Handles complex visualizations like wedges/dashes for stereochemistry and "inside" double bonds for rings.

### CustomQtInteractor (`custom_qt_interactor.py`)
A QWidget wrapper for the VTK render window.
- **Role**: Embeds the VTK 3D view into the PyQt application.

### CustomInteractorStyle (`custom_interactor_style.py`)
Custom VTK interaction style for the 3D viewer.
- **Role**: Intercepts mouse events to enable "Atom Picking" and specific editing interactions (like dragging atoms) while maintaining standard trackball camera controls.

## Helpers & Mixins

### Dialog3DPickingMixin (`dialog3_d_picking_mixin.py`)
A mixin class used by various 3D dialogs (like `TranslationDialog`, `AlignPlaneDialog`).
- **Role**: Provides common functionality to handle "picking" atoms in the 3D view and updating the dialog's state (e.g., "Selected: Atom 5") when the user clicks in the 3D window.
