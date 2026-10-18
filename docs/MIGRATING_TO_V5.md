# Migrating to MoleditPy 5

MoleditPy 5.0.0a1 uses one `moleditpy` package on Windows, macOS and Linux, Python 3.9–3.14, PySide6-Essentials, RDKit and PyVista. The core does not install, import or bundle Open Babel. The separate Linux distribution and its build/publish pipeline have been removed.

```sh
python -m pip uninstall moleditpy-linux openbabel-wheel
python -m pip install --upgrade --pre 'moleditpy>=5.0.0a1,<6'
```

Use a fresh virtual environment when migrating from PyQt6 to avoid loading two Qt bindings. Open Babel can also be uninstalled if no independently installed plugin or other application needs it.

Fallback conversion tries RDKit and then Direct. Direct preserves the sketch and adds missing hydrogens. Saved Open Babel conversion settings migrate to Fallback; retired built-in Open Babel force fields migrate to RDKit MMFF94s. Extra force fields can be registered by an optional optimizer plugin using the existing callback API.

Plugin authors must change PyQt6 imports to PySide6, `pyqtSignal`/`pyqtSlot`/`pyqtProperty` to `Signal`/`Slot`/`Property`, and use Shiboken6 wrapper validity checks. Never mix Qt bindings in one process. Qt worker callbacks that update widgets need a GUI-thread QObject receiver and typed slots.

The core is LGPL-3.0-or-later. Independently developed plugins may choose any license; existing plugin licenses remain unchanged. Read [LICENSING.md](../LICENSING.md) and [DISTRIBUTION.md](../DISTRIBUTION.md) for the exact permission and frozen-library replacement/source instructions.

MoleditPy 5.0.0a1 has been published from the migration branch as a GitHub prerelease and on PyPI. Only the main package was published; existing Linux-package downloads remain available. The shortcut helper and migrated plugin PRs remain unmerged, and their new package versions have not been published.

The app sets `QT_API=pyside6` inside its launch process; it does not change persistent operating-system environment variables.
