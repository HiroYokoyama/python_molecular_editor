# Analysis and Calculations

This document details the architecture for computational tasks, including 2D-to-3D conversion, property analysis, and thread management.

## Overview

Computational tasks in MoleditPy are designed to keep the main UI thread responsive by offloading heavy calculations to specialized background workers managed by a centralized manager.

| Component | Responsibility | Key Class | Location |
| :--- | :--- | :--- | :--- |
| **Worker** | Executes heavy calculations (conformer generation, energy minimization). | `CalculationWorker` | `ui/calculation_worker.py` |
| **Compute Manager**| Orchestrates 3D embedding, thread lifecycles, UI feedback, and results synchronization. | `ComputeManager` | `ui/compute_logic.py` |
| **Analysis** | Calculates physiochemical properties (MW, LogP, etc.). | `AnalysisWindow` | `ui/analysis_window.py` |

---

## 1. Calculation Worker (`ui/calculation_worker.py`)

The `CalculationWorker` inherits from `QObject` and is designed to run within a dedicated `QThread`, ensuring that 3D operations do not freeze the interface.

### Operations
1.  **3D Embedding**: Uses RDKit's ETKDGv2 algorithm for primary conformer generation.
2.  **Stereo Enforcement**: Explicitly applies E/Z and Wedge/Dash constraints to the 3D conformer.
3.  **Fallback Logic**: Automatically switches to Open Babel's `make3D()` (running in an isolated subprocess) if RDKit embedding fails.
4.  **Optimization**: Minimizes the molecule's potential energy using force fields (MMFF94, UFF).

---

## 2. Stability and Isolation

Computational chemistry libraries (C++ extensions) can sometimes cause hard crashes. MoleditPy uses multiple isolation strategies:

### Subprocess Isolation
Calculations involving Open Babel's `make3D()` are executed within a dedicated Python subprocess. This prevents a potential Open Babel segmentation fault from terminating the core application.

### Exception Safety
All background computations and their UI feedback loops are protected by a global exception hook in `main.py`. This ensures that even unexpected calculation failures are captured in the terminal logs (`ERROR`) without a silent application exit.

---

## 3. UI Synchronization (`ui/compute_logic.py`)

The `ComputeManager` bridges the gap between the scientific logic and the 3D viewer.

### Key Features
- **Progress Tracking**: Updates the status bar and viewer overlays during long-running tasks.
- **Cancellation**: Provides a unique `worker_id` system to halt specifically requested threads.
- **Final Projection**: Transfers the result (RDKit Conformer) into the `View3DManager` for visual rendering.
