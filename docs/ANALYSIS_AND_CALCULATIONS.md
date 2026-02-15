# Analysis and Calculations

This document details the architecture for computational tasks, including 2D-to-3D conversion, property analysis, and thread management.

## Overview

Computational tasks in MoleditPy are designed to keep the main UI thread responsive. Heavy calculations are offloaded to background threads using a worker-signal pattern.

| Component | Responsibility | Key Class/File |
| :--- | :--- | :--- |
| **Worker** | Executes heavy calculations (embedding, optimization) in a background thread. | `CalculationWorker` (`calculation_worker.py`) |
| **Compute Mixin** | Manages the worker thread, UI state (progress bars, buttons), and results within the Main Window. | `MainWindowCompute` (`main_window_compute.py`) |
| **Analysis** | Calculates physiochemical properties (MW, LogP, etc.) and displays them. | `AnalysisWindow` (`analysis_window.py`) |

---

## 1. Calculation Worker (`calculation_worker.py`)

The `CalculationWorker` class inherits from `QObject` and is designed to run in a structured `QThread`.

### Workflow
1.  **Instantiation**: Created by `MainWindowCompute` when a calculation is requested.
2.  **Input**: Receives a `mol_block` (text) and an `options` dictionary via the `start_work` signal.
3.  **Processing** (`run_calculation`):
    - **Parsing**: Converts the MOL block to an RDKit molecule.
    - **Stereo Handling**: Respects explicit E/Z stereochemistry labels from the 2D editor, applying them *after* hydrogen addition to prevent RDKit from overriding them.
    - **Embedding**: Generates 3D coordinates.
        - *Modes*: `rdkit` (ETKDGv2), `direct` (2D coords + Z-offset), `fallback` (tries RDKit then others).
    - **Optimization**: Minimizes energy using force fields (MMFF94, MMFF94s, UFF).
4.  **Output**: Emits `finished` with the resulting RDKit molecule, or `error` if failed.

### Key Logic
- **Explicit Stereo Preservation**: The worker aggressively enforces stereochemistry defined in the 2D editor (e.g., specific E/Z configurations) by setting `BondStereo` tags and using distance constraints during embedding if necessary.
- **Explicit Stereo Preservation**: The worker aggressively enforces stereochemistry defined in the 2D editor (e.g., specific E/Z configurations) by setting `BondStereo` tags and using distance constraints during embedding if necessary.
- **Worker ID & Halting**:
    - Each calculation is assigned a unique `worker_id` by the main window.
    - The worker periodically checks a shared `halt_ids` set (passed by reference).
    - If its `worker_id` is found in `halt_ids`, the worker raises a `RuntimeError("Halted")`, allowing cleaner termination than forcefully killing the thread.

---

## 2. Main Window Integration (`main_window_compute.py`)

This mixin integrates the computation logic into the main application window.

### Thread Management
- **Lifecycle**: Creates a new `QThread` and `CalculationWorker` pair for each run.
- **Cleanup**: Connects `finished` signals to automatically quit the thread and delete objects (`deleteLater`).
- **Safety**: Maintains a list of active threads to prevent garbage collection during execution.

### User Interaction
- **Menus**: dynamic right-click menus on the "Convert" and "Optimize" buttons allow users to choose specific algorithms (e.g., "MMFF94s only", "Direct conversion").
- **Visual Feedback**: Updates the status bar and displays a "Calculating..." text actor in the 3D view.
- **Halting**: Transforms the "Convert" button into a "Halt" button during execution, creating a responsive UX.

---

## 3. Property Analysis (`analysis_window.py`)

The `AnalysisWindow` provides a summary of molecular properties.

### Modes
1.  **Standard Mode (MOL-derived)**:
    - Used when a chemical graph is available (connectivity is known).
    - Calculates: SMILES, InChI, InChIKey, Formula, MW, Exact Mass, LogP, TPSA, H-Bond Donors/Acceptors, Ring Count.
    - Uses `rdkit.Chem.Descriptors` and `rdMolDescriptors`.
2.  **XYZ Mode (Coordinate-derived)**:
    - Used when loading raw XYZ files where bonding is uncertain.
    - Calculates: Formula, MW, Exact Mass, Atom Counts.
    - **Limitation**: Does not calculate connectivity-dependent properties (LogP, SMILES) to avoid errors from incorrect bond guessing.
