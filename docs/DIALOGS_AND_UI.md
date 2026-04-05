# Dialogs and UI Components

MoleditPy contains over a dozen specialized dialogs for chemical editing and application configuration. These are managed centrally by the `DialogManager` (`ui/dialog_logic.py`).

## Overview

Post-v3.0, all dialogs have been refactored to prioritize **State Isolation** and **Diagnostic Visibility**:

- **Scientific Base Classes**: Most editing dialogs inherit from `BasePickingDialog` or `Dialog3DPickingMixin` to ensure consistent 3D interaction logic and synchronized lifecycle management.
- **UI Logic Delegation**: Presentation logic is strictly separated from scientific state; dialogs delegate compute tasks back to the host's specialized Managers (e.g., `ComputeManager`).
- **Never Hide Interaction Errors**: To prevent silent failures during 3D picking or coordinate transforms, all interaction methods use explicit `try-except` blocks or `hasattr` checks that **never fail silently**. Every deviation is logged via `logging.error` to ensure developers maintain absolute visibility into the UI/Logic interface.

| Component | Responsibility | Location |
| :--- | :--- | :--- |
| **`DialogManager`**| Unified registry for instantiating and lifecycle-tracking all application-wide dialogs. | `ui/dialog_logic.py` |
| **`UIManager`**| Manages the high-level layout, toolbar state, and coordinate system synchronization. | `ui/ui_manager.py` |

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
