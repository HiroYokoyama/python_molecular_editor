# Installation

It is **strongly recommended** to install the required packages using the provided `requirements.txt` file. This ensures that you are using a set of libraries (including PyVista and PyQt) that are tested and confirmed to be compatible.

Recent updates to some core dependencies (specifically PyQt 6.10.0 and newer) have introduced regressions that are known to cause a **`Segmentation fault`** (crash) when used with PyVista.

## Recommended Setup (Using `requirements.txt`)

Please follow these steps to ensure a stable environment.

1.  **(Optional but Recommended)** Create and activate a new virtual environment:

    ```bash
    # For standard Python venv
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate

    # Or for Conda
    # conda create -n moleditpy python=3.10
    # conda activate moleditpy
    ```

2.  Install all dependencies from the `requirements.txt` file:

    ```bash
    pip install -r requirements.txt
    ```

    This file locks `pyqt6` to a stable version (like `6.9.1`) that avoids this critical bug.

    For Linux users, please use `requirements-linux.txt` file:

    ```bash
    pip install -r requirements-linux.txt
    ```

### :warning: Avoid Manual Installation

**Do not** install packages manually (e.g., `pip install pyvista pyqt6`). This will likely fetch the latest, incompatible versions and cause the application to crash.
