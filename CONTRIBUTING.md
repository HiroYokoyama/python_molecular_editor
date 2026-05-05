# Contributing to MoleditPy

Thank you for your interest in contributing to **MoleditPy**! We welcome contributions from the community to make this molecular editor more robust and feature-rich for scientific research.

To maintain the high quality and stability of the application (Pylint > 9.0, Core Logic Coverage > 80%), please review the following guidelines.

## 1. How to Report Bugs

Since MoleditPy is designed to be **fail-safe** (it avoids crashing via extensive exception handling), bugs may appear as "silent failures" or console errors rather than application crashes.

When opening an Issue, please include:
1.  **Steps to Reproduce:** A detailed checklist of actions taken.
2.  **Expected Behavior vs. Actual Behavior.**
3.  **Logs / Traceback:** This is critical. Please attach the contents of the console output or the `app_errors.log` file. Even if the app did not crash, the error details are likely recorded there.
4.  **Environment:** OS, Python version, and MoleditPy version, installed plugins.

## 2. Development Setup

1.  Clone the repository.
2.  Install the package in editable mode (recommended for development) or standard mode:
    ```bash
    # Standard install
    pip install .

    # Or for active development (changes reflect immediately)
    # pip install -e .
    ```
3.  Install development tools (pylint, pytest, coverage):
    ```bash
    pip install pylint pytest coverage
    ```

## 3. Coding Standards

We maintain a strict code quality standard. Before submitting a Pull Request (PR), ensure your code meets the following:

* **Pylint Score:** The project maintains a score of **> 9.0/10**.
    * Run `pylint src/moleditpy` and fix warnings.
    * Do not suppress warnings locally without a valid reason.
* **Style:** Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/).
* **Type Hinting:** Use Python type hints for all function arguments and return values in Core Logic.

## 4. Testing & Stability Policy

MoleditPy uses a hybrid approach to quality assurance: **Automated Unit Tests** for logic and **Manual Checklists** for UI.

### A. Core Molecular Logic (Automated)
* Files in `core/` related to calculation, parsing, or data structure (e.g., `molecular_data.py`, `calculation_worker.py`) **MUST** have unit tests.
* **Target Coverage:** Maintain or exceed **80%** coverage for logic modules.
* Any new public method added to `plugin_interface.py` (`PluginContext` or `Plugin3DController`) **MUST** have a corresponding test in `tests/unit/test_plugin_interface.py`.

### B. Error Handling — No Hiding, No Crashing

MoleditPy's stability policy is: **never hide errors, never crash.**

| Situation | Rule |
|-----------|------|
| UI slot / event handler | Wrap in `try/except`, log the traceback, show a status message or dialog |
| Plugin `initialize()` / `run()` | Same — always log; optionally show a warning to the user |
| Internal helper (non-UI) | Let exceptions propagate; do NOT swallow them |
| C-extension boundary (RDKit, PyVista, PyQt6) | Broad `except Exception` is acceptable **only here**; must still log |

**Critical rule:** an `except` block that contains only `pass` is **never** acceptable. You must at minimum write to the log:

```python
# WRONG
try:
    do_something()
except Exception:
    pass

# CORRECT
try:
    do_something()
except Exception:
    import traceback
    print(traceback.format_exc())   # or use the app logger
```

### C. GUI & Interaction (Manual + Defensive)
* **Avoid Broad Exceptions:** Minimize the use of `try-except Exception` blocks. Always prefer granular, specific exception types (e.g., `AttributeError`, `ValueError`, `RuntimeError`, `OSError`) to prevent masking unexpected bugs. Use broad `Exception` only at the highest level of a UI slot or when interfacing with unpredictable C-extensions.
* **Manual Testing:** For UI changes, you must verify the behavior using the **Manual Test Checklist**. Please confirm in your PR description that you have manually verified the fix.

## 5. Pull Request Process

1.  **Branching:** Create a new branch for your feature or fix (e.g., `feature/new-bond-algo` or `fix/dialog-layout`).
2.  **Verification:**
    * Run tests: `pytest`
    * Check linting: `pylint src/moleditpy`
3.  **Description:** Describe your changes clearly. If it's a UI change, screenshots are highly appreciated.
4.  **Review:** Wait for a maintainer to review your code. We may ask for changes to meet the Pylint score or logging requirements.

## 6. License

By contributing, you agree that your contributions will be licensed under the project's [LICENSE](./LICENSE).

---
*Happy Coding!*
