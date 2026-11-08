# Contributing to MoleditPy

Thank you for your interest in contributing to **MoleditPy**! We welcome contributions from the community to make this molecular editor more robust and feature-rich for scientific research.

To maintain the high quality and stability of the application (Pylint > 9.0, Full Application Coverage ~ 80%), please review the following guidelines.

## 1. How to Report Bugs

Since MoleditPy is designed to be **fail-safe** (it avoids crashing via extensive exception handling), bugs may appear as "silent failures" or console errors rather than application crashes.

When opening an Issue, please include:
1.  **Steps to Reproduce:** A detailed checklist of actions taken.
2.  **Expected Behavior vs. Actual Behavior.**
3.  **Logs / Traceback:** This is critical. Before sharing, remove usernames, local paths, project data, plugin details, and secrets. If the logs cannot be safely sanitized, ask a maintainer for a private submission option. Please attach the console/terminal output, copy the traceback from the error dialog's "Show Details" button, or attach `~/.moleditpy/moleditpy.log` (if "Save log to file" is enabled in Settings ▸ Other). Even if the app did not crash, the error details are recorded there.
4.  **Environment:** OS, Python version, and MoleditPy version, installed plugins.

## 2. Development Setup

1.  Clone the repository.
2.  Install the package in editable mode (recommended for development) or standard mode:
    ```bash
    # Standard install
    pip install ./moleditpy

    # Or for active development (changes reflect immediately)
    pip install -e moleditpy/
    ```
3.  Install development tools and test dependencies:
    ```bash
    pip install pytest pytest-qt pytest-cov pytest-timeout pylint
    ```

## 3. Coding Standards

We maintain a strict code quality standard. Before submitting a Pull Request (PR), ensure your code meets the following:

* **Pylint Score:** The project maintains a score of **> 9.0/10**.
    * Run `pylint moleditpy/src/moleditpy/` and fix warnings.
    * Do not suppress warnings locally without a valid reason.
* **Style:** Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/).
* **Type Hinting:** Use Python type hints for all function and method arguments and return values.

## 4. Testing & Stability Policy

MoleditPy uses a hybrid approach to quality assurance: **Automated Unit Tests** for logic and **Manual Checklists** for UI.

### A. Automated Tests & Coverage
* Core files in `moleditpy/src/moleditpy/core/` related to molecular geometry or data structures (e.g., `molecular_data.py`), or files related to calculation/parsing (e.g., `calculation_worker.py` and `io_logic.py` in `moleditpy/src/moleditpy/ui/`) **MUST** have unit tests.
* **Target Coverage:** Maintain or exceed **80%** full application coverage.
* Any new public method added to `plugin_interface.py` (`PluginContext` or `Plugin3DController`) **MUST** have a corresponding test in `tests/unit/plugins/test_plugin_interface.py`.

### B. Error Handling — No Hiding, No Crashing

MoleditPy's stability policy is: **never hide errors, never crash.**

| Situation | Rule |
|-----------|------|
| UI slot / event handler | Wrap in `try/except`, log the traceback, show a status message or dialog |
| Plugin `initialize()` / `run()` | Same — always log; optionally show a warning to the user |
| Internal helper (non-UI) | Let exceptions propagate; do NOT swallow them |
| C-extension boundary (RDKit, PyVista, PySide6) | Broad `except Exception` is acceptable **only here**; must still log |

**Critical rule:** an `except` block that contains only `pass` is **never** acceptable. You must at minimum write to the log:

```python
# WRONG
try:
    do_something()
except Exception:
    pass

# CORRECT
import logging

try:
    do_something()
except Exception as e:
    logging.exception("Failed to do something: %s", e)
```

### C. GUI & Interaction (Manual + Defensive)
* **Avoid Broad Exceptions:** Minimize the use of `try-except Exception` blocks. Always prefer granular, specific exception types (e.g., `AttributeError`, `ValueError`, `RuntimeError`, `OSError`) to prevent masking unexpected bugs. Use broad `Exception` only at the highest level of a UI slot or when interfacing with unpredictable C-extensions.
* **Manual Testing:** For UI changes, you must verify the behavior using the **Manual Test Checklist**. Please confirm in your PR description that you have manually verified the fix.

## 5. Core Feature Acceptance Policy

MoleditPy follows the design philosophy:

**"Simplicity at the Core. Limitless Potential in the Extensions."**

We welcome contributions that improve MoleditPy. However, not every useful feature needs to become part of the core application.

Our goal is to keep the core simple, stable, safe, and maintainable while allowing users to customize and extend their workflows through plugins.

### What Belongs in the Core?

We generally welcome contributions that:

- Fix bugs or improve stability, performance, and reliability.
- Improve essential molecular editing, visualization, and file handling functionality.
- Improve the usability, accessibility, or safety of existing features.
- Refactor code to improve maintainability without introducing unnecessary complexity.
- Improve testing, documentation, or development infrastructure.
- Enhance the plugin system or provide broadly useful APIs for extensions.

Core features should provide a clear benefit to the application's fundamental functionality and remain safe and predictable for general users.

### What Belongs in a Plugin?

Features that are useful primarily for specific users or workflows are generally better suited as plugins.

Examples include:

- Specialized scientific calculations, analyses, or visualizations.
- Integrations with external scientific software.
- Workflow-specific automation and productivity tools.
- Optional UI customizations or specialized shortcuts.

Such features can still be valuable contributions to the MoleditPy ecosystem, even if they are not included in the core.

### How Are New Features Evaluated?

When reviewing feature proposals and pull requests, maintainers consider:

1. **General usefulness:** Does the feature improve MoleditPy's fundamental functionality?
2. **Safety and usability:** Is the feature appropriate and predictable for general users?
3. **Maintainability:** Can the feature be maintained without introducing unnecessary long-term complexity?
4. **Extensibility:** Could the feature be implemented as a plugin instead?

These are guidelines rather than strict rules. Features do not need to be used by every user to justify inclusion in the core.

Even well-implemented and thoroughly tested features may be better suited as plugins.

### Plugin Flexibility and Compatibility

Plugin developers are free to access internal core APIs when necessary. However, internal APIs are not guaranteed to remain stable across MoleditPy versions.

Plugins that depend on internal implementation details may require updates when the core changes.

Contributors are also welcome to propose improvements to the public plugin API when they would provide broadly reusable functionality.

### Before Submitting a Feature PR

For significant new features, we encourage contributors to open an Issue to discuss whether the feature belongs in the core or would be better implemented as a plugin.

A decision not to merge a feature into the core does not necessarily mean the feature is unwanted or poorly implemented. It may simply be better suited to an optional extension.

**We value all contributions and aim to keep MoleditPy both easy to use and flexible enough to support diverse scientific workflows.**

## 6. Pull Request Process

1.  **Branching:** Create a new branch for your feature or fix (e.g., `feature/new-bond-algo` or `fix/dialog-layout`).
2.  **Verification:**
    * Run tests: Run the unified test runner using `python tests/run_all_tests.py` (or specify `--unit` / `--integration` as needed, and add `--headless` or prefix with `MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen` if running in a headless environment like CI or WSL).
    * Check linting: `pylint moleditpy/src/moleditpy/` (or run `python tests/run_all_tests.py --pylint`).
3.  **Description:** Describe your changes clearly. If it's a UI change, screenshots are highly appreciated.
4.  **Review:** Wait for a maintainer to review your code. We may ask for changes to meet the Pylint score or logging requirements.

## 7. License

By contributing, you agree that your contributions will be licensed under the project's [LICENSE](./LICENSE).

---
*Happy Coding!*
