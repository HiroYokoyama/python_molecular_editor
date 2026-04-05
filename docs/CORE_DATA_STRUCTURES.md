# Core Data Structures & Scene Graph

This document details the internal data models and the visual scene graph architecture used in **MoleditPy**.

## Overview

The application maintains a strict separation between the persistent chemical data model (`MolecularData`) and the interactive visual representation (`MoleculeScene`).

| Component | Responsibility | Key Class | Location |
| :--- | :--- | :--- | :--- |
| **Data Model** | Pure chemical state (atoms, bonds, connectivity). | `MolecularData` | `core/molecular_data.py` |
| **State Manager**| Undo/Redo tracking and state serialization. | `StateManager` | `ui/app_state.py` |
| **Persistent IO** | Project saving and loading logic. | `IOManager` | `ui/io_logic.py` |
| **Scene Graph** | 2D visualization and interaction logic. | `MoleculeScene` | `ui/molecule_scene.py` |
| **Visual Items** | Rendering individual atoms and bonds. | `AtomItem`, `BondItem` | `ui/` package |

---

## 1. Data Model (`core/molecular_data.py`)

`MolecularData` is the source of truth for the chemistry being edited. It is designed to be independent of the UI framework.

### Internal Storage
- **`atoms`**: A dictionary mapping a unique `atom_id` (int) to a property dictionary:
    - `"symbol"`: Element symbol (e.g., "C", "N", "O").
    - `"pos"`: Coordinates (stored as `PointTuple` for UI independence).
    - `"charge"`: Formal charge (int).
    - `"radical"`: Number of radical electrons (int).
- **`bonds`**: A dictionary mapping an atom ID pair `(id1, id2)` to bond properties:
    - `"order"`: Bond order (1, 2, 3, 1.5).
    - `"stereo"`: Stereochemistry flag (0=None, 1=Wedge, 2=Dash, 3=Z, 4=E).
- **`adjacency_list`**: A dictionary mapping `atom_id` to a set of neighbor IDs for high-performance topological traversal.

### Key Operations
- **`to_rdkit_mol()`**: The bridge to RDKit. It handles the reconstruction of a 3D-aware RDKit molecule object, preserving atom properties and explicit stereochemical directions (Wedges/Dashes).
- **`from_rdkit_mol()`**: populates the `atoms` and `bonds` dictionaries from an existing RDKit molecule.

---

## 2. Interactive Scene (`ui/molecule_scene.py`)

`MoleculeScene` manages the 2D canvas using the Qt Graphics View Framework. It is composed of multiple functional handlers (located in `ui/molecular_scene_handler.py`):

- **Template Handling**: Manages ghost-previews of chemical fragments (rings) and handles their merger with existing structures.
- **Keyboard Interaction**: Maps shortcuts (e.g., `1`, `2`, `3` for bond orders, `C`, `O` for elements) to local editing actions.
- **Query Logic**: Provides spatial methods (e.g., `get_atom_at(pos)`) used by mouse event handlers.

---

## 3. Visual Items (`ui/atom_item.py`, `ui/bond_item.py`)

The visual representation layer handles custom painting and user feedback.

### AtomItem
Renders the element symbol and related notations.
- **Dynamic Labeling**: Automatically calculates and renders implicit hydrogens (e.g., "H", "CH₂") based on the atom's valency and current bonding state.
- **Electronic State**: Renders formal charges and radical dots using standard chemical notation.
- **Syncing**: Holds a reference back to its `atom_id` in the `MolecularData` model.

### BondItem
Renders a line or wedge between two `AtomItem` objects.
- **Order Rendering**: Renders single, double, triple, and aromatic (dashed) bonds. 
- **Stereo Rendering**: Implements filled wedges and hashed dashes for stereochemical visualization.
- **Dynamic Layout**: Automatically updates its start and end points as the connected `AtomItem`s are moved on the canvas.
