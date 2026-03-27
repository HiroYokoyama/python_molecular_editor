# Analysis and Calculations

This document details the architecture for computational tasks, including 2D-to-3D conversion, property analysis, and thread management.

## Overview

Computational tasks in MoleditPy are designed to keep the main UI thread responsive. Heavy calculations are offloaded to background threads using a worker-signal pattern managed by specialized logic managers.

| Component | Responsibility | Key Class/File | Location |
| :--- | :--- | :--- | :--- |
| **Worker** | Executes heavy calculations (embedding, optimization) in a background thread. | `CalculationWorker` | `core/calculation_worker.py` |
| **Compute Engine**| Orchestrates 3D embedding algorithms and library coordination. | `ComputeEngine` | `core/compute_engine.py` |
| **Compute Manager**| Manages the worker thread, UI state, and results within the UI layer. | `ComputeManager` | `ui/compute_logic.py` |
| **Analysis** | Calculates physiochemical properties (MW, LogP, etc.) and displays them. | `AnalysisWindow` | `ui/analysis_window.py` |

---

## 1. Calculation Worker (`core/calculation_worker.py`)

The `CalculationWorker` class inherits from `QObject` and is designed to run in a structured `QThread`.

### Workflow
1.  **Instantiation**: Created by `ComputeManager` (as part of the `MainWindow` state) when a calculation is requested.
2.  **Input**: Receives a `mol_block` (text) and an `options` dictionary via the `start_work` signal.
3.  **Processing** (`run_calculation`):
    - **Parsing**: Converts the MOL block to an RDKit molecule.
    - **Stereo Handling**: Respects explicit E/Z stereochemistry labels from the 2D editor, applying them *after* hydrogen addition to prevent RDKit from overriding them.
    - **Embedding**: Generates 3D coordinates.
        - *Modes*: `rdkit` (ETKDGv2), `direct` (2D coords + Z-offset), `fallback` (tries RDKit then others).
    - **Optimization**: Minimizes energy using force fields (MMFF94, MMFF94s, UFF).
        - *Robustness*: RDKit optimization failures are treated as non-fatal warnings, ensuring the pipeline continues even if specific algorithms fail.
4.  **Output**: Emits `finished` with the resulting RDKit molecule (and original worker ID), or `error` if failed.

---

## 2. Stability and Isolation

To ensure application stability, heavy or potentially unstable third-party library calls are isolated.

### Open Babel Subprocess Isolation
Calculations involving Open Babel's `make3D()` are executed within a dedicated Python subprocess.
- **Reasoning**: Open Babel's C++ core can occasionally trigger hard crashes (segfaults) that would otherwise terminate the entire PyQt application.
- **Implementation**: The `CalculationWorker` spawns a subprocess using `sys.executable` and passes the molecular data via stdin.
- **Timeout**: A 20-second timeout is enforced to prevent zombie processes in case of hung calculations.

### Non-Fatal Optimization
When performing Force Field (FF) optimizations (MMFF94, UFF):
- The worker treats failures to converge or missing parameters as non-fatal.
- If optimization fails, the worker provides the unoptimized (but valid) coordinate set to the user rather than crashing the pipeline.

### Worker ID & Halting
- Each calculation is assigned a unique `worker_id` by the `ComputeManager`.
- The worker periodically checks a shared `halt_ids` set.
- If its `worker_id` is found in `halt_ids`, the worker terminates gracefully, allowing the user to cancel long-running tasks.

---

## 3. UI Integration (`ui/compute_logic.py`)

The `ComputeManager` handles the high-level orchestration of calculations from the UI perspective.

### Manager Responsibilities
- **Thread Management**: Creates and cleans up `QThread` and `CalculationWorker` pairs.
- **User Feedback**: Updates the status bar and manages the "Calculating..." overlay in the 3D viewer.
- **Result Processing**: Delegates the final result (RDKit molecule) back to the `MainWindow` for rendering and state updates.
- **Method Selection**: Provides right-click menus for choosing specific optimization methods (GAFF, MMFF94, etc.).

---

## 4. Property Analysis (`ui/analysis_window.py`)

The `AnalysisWindow` provides a summary of molecular properties.

### Modes
1.  **Standard Mode (MOL-derived)**: Calculates SMILES, InChI, InChIKey, Formula, MW, Exact Mass, LogP, TPSA, H-Bond counts, and Ring Count.
2.  **XYZ Mode (Coordinate-derived)**: Used for raw XYZ files. Calculates Formula and MW but suppresses connectivity-dependent properties (LogP, SMILES) to avoid inaccuracies from guessed bonds.
