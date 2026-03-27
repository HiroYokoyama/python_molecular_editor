# Dialogs and UI Components Documentation

This document provides a detailed overview of the various dialogs and custom UI components used in **MoleditPy**. All UI-related files are located in `moleditpy/src/moleditpy/ui/`.

## General Dialogs

### SettingsDialog (`settings_dialog.py`)
The central configuration hub.
- **Architecture**: UI is separated into a base class and specialized tab modules in the `settings_tabs/` subdirectory.
- **Tabs**:
    - **2D Settings**: Colors, fonts, and bond visual styles.
    - **3D Scene**: Backgrounds, projection modes, and lighting effects.
    - **3D Styles**: Detailed configuration for Ball & Stick, CPK, Wireframe, and Stick models.

### AboutDialog (`about_dialog.py`)
Displays application version, credentials, and legal information.

### PeriodicTableDialog (`periodic_table_dialog.py`)
Interactive grid for element selection. Used for one-shot atom changes or global element settings.

---

## Geometry & Manipulation Dialogs

These dialogs provide precision tools for 3D molecular editing. Many utilize the `Dialog3DPickingMixin` for interactive atom selection.

| Dialog | Functional File | Logic/Helper | Description |
| :--- | :--- | :--- | :--- |
| **Translation** | `translation_dialog.py` | | Move atoms or molecules to specific coordinates. |
| **Move Group** | `move_group_dialog.py` | | BFS-based connected group movement and rotation. |
| **Align Plane** | `align_plane_dialog.py` | | PCA-based alignment of atoms to Cartesian planes. |
| **Alignment** | `alignment_dialog.py` | | Align bond axes to Cartesian X, Y, or Z axes. |
| **Mirror** | `mirror_dialog.py` | | Generate enantiomers via plane reflection. |
| **Planarize** | `planarize_dialog.py` | | Project atoms onto their best-fit plane. |

### Bond & Angle Adjustment
- **BondLengthDialog (`bond_length_dialog.py`)**: Precise adjustment of inter-atomic distances.
- **AngleDialog (`angle_dialog.py`)**: Adjustment of bond angles with various "arm" rotation modes.
- **DihedralDialog (`dihedral_dialog.py`)**: Torsion angle adjustment.

---

## Analysis & Calculation UI

- **AnalysisWindow (`analysis_window.py`)**: Physicochemical property reporting for graph-based and coordinate-based structures.
- **ConstrainedOptimizationDialog (`constrained_optimization_dialog.py`)**: interface for setting up 3D minimizations with fixed distance/angle constraints.

---

## Core Visual Components

### MoleculeScene (`ui/molecule_scene.py`)
The 2D editing canvas.
- **Architecture**: Uses mixins from `molecular_scene_handler.py` (Template, Keyboard, and Query logic).
- **Items**: Manages `AtomItem` and `BondItem` instances, handling all click/drag interaction logic.

### 3D Viewer (`ui/view_3d.py`)
High-performance rendering engine.
- **Logic**: Delegated to `view_3d_logic.py` (View3DManager).
- **Styles**: Supports Ball & Stick, Stick, Wireframe, and CPK visualizations with VTK/PyVista.

### Interaction Wrappers
- **ZoomableView (`zoomable_view.py`)**: `QGraphicsView` wrapper providing smooth panning and zooming for the 2D scene.
- **CustomQtInteractor (`custom_qt_interactor.py`)**: Embeds the VTK render window into the Qt layout.
- **CustomInteractorStyle (`custom_interactor_style.py`)**: Implements the "Direct 3D Edit" mode where atoms can be dragged in 3D to modify conformations.

---

## Managers & Mixins

### DialogManager (`ui/dialog_manager.py` / `dialog_logic.py`)
Centralized registry that tracks open dialog instances, ensuring only one instance of each tool (e.g., Settings) is active at a time.

### UI Manager (`ui/ui_manager.py`)
Manages high-level UI states, mode switching (Select/Atom/Bond/Template), and the visibility of 2D/3D panels.
