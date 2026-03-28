# Dialogs and UI Components

MoleditPy contains over a dozen specialized dialogs for chemical editing and application configuration. These are managed centrally by the `DialogManager` (`ui/dialog_logic.py`).

## Overview

All dialogs follow a standard pattern:
- **UI Logic Separation**: The visual form (typically PyQt6) is kept as simple as possible, delegating complex calculations to the `MainWindow` or specialized manager modules.
- **Session Persistence**: Dialog states (last used values, window positions) are often serialized.

| Component | Responsibility |
| :--- | :--- |
| **`DialogManager`**| Central hub for instantiating and showing application-wide dialogs. |
| **`MainWindowUiManager`**| Manages the high-level layout, toolbar visibility, and coordinate overlays. |

---

## 1. Scientific Editing Dialogs

These dialogs are used to modify the 3D geometry of the current molecule.

| Dialog | Purpose |
| :--- | :--- |
| **Bond Length / Angle / Dihedral** | Set precise numerical values for internal coordinates. |
| **Translation** | Translate the molecule or selected atoms to coordinate targets. |
| **Align Plane / Align Axis** | Reorient the entire molecule relative to the global unit axes. |
| **Mirror Image** | Generate a mirror image across a specific plane. |
| **Planarize** | Project a set of atoms onto their best-fit plane. |
| **Constrained Optimization** | UI for adding/removing force-field constraints (MMFF94/UFF). |

---

## 2. Interactive Components

| Component | Description |
| :--- | :--- |
| **Periodic Table** | A visual element selector for atom placement and element switching. |
| **User Templates** | Interface for saving and retrieving custom chemical fragments. |
| **Analysis Window** | Real-time display of molecular properties (MW, LogP). |

---

## 3. Global Configuration

| Dialog | Description |
| :--- | :--- |
| **Settings Dialog** | Tabbed interface for Appearance (colors), Rendering (lighting), and IO settings. |
| **About Dialog** | Displays version, DOI, and license metadata. |
