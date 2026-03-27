# Core Data Structures & Scene Graph

This document details the internal data models and the visual scene graph architecture used in **MoleditPy**.

## Overview

The application strictly separates the chemical data model (`MolecularData`) from the visual representation in the 2D editor (`MoleculeScene`, `AtomItem`, `BondItem`).

| Component | Responsibility | Key Class | Location |
| :--- | :--- | :--- | :--- |
| **Data Model** | Pure chemical state (atoms, bonds, connectivity). | `MolecularData` | `core/molecular_data.py` |
| **App State** | Undo/Redo and project serialization. | `MainWindowAppState` | `core/app_state.py` |
| **Scene Graph** | 2D visualization and interaction logic. | `MoleculeScene` | `ui/molecule_scene.py` |
| **Visual Items** | Rendering individual atoms and bonds. | `AtomItem`, `BondItem` | `ui/` package |

---

## 1. Data Model (`core/molecular_data.py`)

The `MolecularData` class is the source of truth for the chemistry being edited.

### Structure
- **`atoms`**: A dictionary mapping unique `atom_id` (int) to a dictionary of properties:
    - `"symbol"`: Element symbol (str, e.g., "C").
    - `"pos"`: QPointF coordinates (scene space).
    - `"charge"`: Formal charge (int).
    - `"radical"`: Number of radical electrons (int).
    - `"item"`: Reference to the corresponding `AtomItem` (for sync).
- **`bonds`**: A dictionary mapping a tuple `(id1, id2)` to bond properties:
    - Order of IDs in the key matters for storage but uniqueness is managed.
    - `"order"`: Bond order (1, 2, 3, 1.5).
    - `"stereo"`: Stereochemistry flag (0=None, 1=Wedge, 2=Dash, 3=Z, 4=E).
- **`adjacency_list`**: A dictionary mapping `atom_id` to a list of connected neighbor IDs for O(1) traversal.

### Key Features

#### `to_rdkit_mol(use_2d_stereo=True)`
This is the **critical bridge** between the 2D drawing and RDKit's chemical logic.
1.  **Reconstruction**: Builds an `RDKit.RWMol` from the `atoms` and `bonds` dictionaries.
2.  **Stereochemistry**:
    - **Wedge/Dash**: Converts `stereo=1` (Wedge) and `stereo=2` (Dash) properties into `Chem.BondDir.BEGINWEDGE` / `BEGINDASH`.
    - **E/Z Isomerism**:
        - If `stereo` is 3 (Z) or 4 (E), it explicitly sets `Chem.BondStereo.STEREOZ` or `STEREOE`.
        - It typically prioritizes explicit user labels over coordinate-based inference if labels exist.

---

## 2. The Scene (`ui/molecule_scene.py`)

`MoleculeScene` inherits from `QGraphicsScene` and several specialized mixin classes defined in `molecular_scene_handler.py`. This decomposition keeps the core scene logic manageable and separates concerns like template handling and keyboard events.

### Mixins
- **`TemplateMixin`**: Manages the template preview system and fragment insertion logic.
- **`KeyboardMixin`**: Centralizes keyboard event handling and shortcuts.
- **`SceneQueryMixin`**: Utility methods for spatial queries and basic item CRUD operations.

### Interaction Modes
The scene behaves differently based on `self.mode`:
- **`select`**: Standard selection, moving atoms, right-click to delete.
- **`atom_[Symbol]`**: e.g., `atom_C`, `atom_O`. Clicking creates/changes atoms. Dragging from an atom creates a bond + new atom.
- **`bond_[Order]_[Stereo]`**: e.g., `bond_1_0` (Single), `bond_1_1` (Wedge), `bond_2_0` (Double). Clicking existing bonds changes their type. Dragging creates specific bonds.
- **`charge_plus` / `charge_minus` / `radical`**: Clicking atoms modifies their electronic state.
- **`template_[Type]`**: Shows a ghost preview and places a chemical fragment (e.g., Benzene) on click.

### Event Handling
- **`mousePressEvent`**:
    - Detects clicks on atoms/nodes.
    - Initiates "drag line" for bonding if in atom/bond mode.
    - Handles Right-Click for deletion (context-sensitive).
- **`mouseMoveEvent`**:
    - Updates the "temp line" (dotted red line) during bond creation dragging.
    - Updates template previews.
- **`mouseReleaseEvent`**:
    - Finalizes actions (creates atoms/bonds, applying templates).
    - Updates the data model and visual items.
    - Triggers Undo State capture.

---

## 3. Visual Items

These classes inherit from `QGraphicsItem` and handle the actual rendering code (painting).

### `AtomItem` (`atom_item.py`)
Renders a single atom.
- **`paint()`**:
    - Draws the atom symbol centered or aligned based on connectivity.
    - **Implicit Hydrogens**: Calculates and renders "H", "CH₂", etc., automatically based on valency if configured.
    - **Charge/Radical**: Draws standard chemistry notations (+, -, dots) adjacent to the symbol.
    - **Highlighting**: Draws selection or hover rectangles.
- **`shape()`**: Defines a precise collision circle for mouse detection, separate from the visual bounding box.

### `BondItem` (`bond_item.py`)
Renders a bond between two `AtomItem`s.
- **Dynamic Updates**: Tracks the positions of its connected atoms (`atom1`, `atom2`) and updates its line coordinates automatically.
- **Rendering Logic**:
    - **Single**: Simple line.
    - **Double**: Two parallel lines.
        - **Ring Detection**: If part of a ring, it draws one center line and one shorter "inner" line (classic chemical style).
    - **Triple**: Three parallel lines.
    - **Wedge/Dash**: Renders a standard filled wedge triangle or a series of hash lines.
- **E/Z Support**: Can render "E" or "Z" labels directly on the bond for absolute clarity in stereochemical editing.
