# Main application audit — October 2026

Assisted-by: GPT-6.1 Sol

The source review covered all 79 Python files under `moleditpy/src/moleditpy`
at `48c4e88d` (35,149 lines), including the entry points, core model and
stereochemistry, every manager and dialog, scene and interactor code, plugins,
settings tabs, and utilities. The added project validator was reviewed with
its callers and regression tests. This records the scope of the read-through;
it does not claim that no further defects exist.

## Confirmed fixes

- Invalid 3D label preferences and built-in setting types could raise during
  startup, rendering, or opening settings. Invalid values now use defaults;
  unknown plugin preferences remain available.
- Axis-alignment dialogs lacked the molecule invalidation hook. Plugin molecule
  setters also bypassed invalidation. Replacement now closes the old tools.
- XYZ coordinates could contain NaN/infinity, and the charge dialog accepted
  fractional values or raised on infinity. Inputs are validated before use.
- A plugin installed from its own destination folder deleted its source.
  Installations also reject copying a directory into itself.
- User-template serialization, previews, and placement dropped R/S annotations.
  They now retain labels and refresh the cached preview when labels change.
- OBJ color splitting could reuse the previous actor's converted mesh.
  Each actor now supplies its own geometry.
- XYZ exports lost original pseudo-atom labels. File and plugin exports retain
  those labels and their coordinates.
- Structurally corrupt projects could clear the document before failing.
  JSON and raw geometry are checked before clearing document or plugin state.
- Closing during conversion could outlive the 200 ms thread wait. Shutdown now
  cancels work and resumes closing after native thread cleanup, without blocking
  the GUI. A controlled thread test uses events rather than timing sleeps.
- Extreme finite geometry input could overflow slider conversion and leave
  signals blocked. Values are clamped before integer conversion and the prior
  signal-blocking state is restored.

Every confirmed defect has a regression. Existing test assertions are retained;
two raw-file test fixtures now contain valid atom records instead of placeholder
strings that the real loader could never restore.

## Dependency decision

[PR #168](https://github.com/HiroYokoyama/python_molecular_editor/pull/168)
groups RDKit 2026.9.1 with PySide6-Essentials 6.12.0. These changes require
separate decisions.

- **RDKit:** permit the 2026.9 series (`<2026.10`). An isolated Python 3.13
  environment using RDKit 2026.9.1 and the supported PySide6 6.11.2 passed all
  2,942 then-existing headless tests (3 skipped) and all 35 native GUI tests.
  The final suite also includes the new audit regressions.
- **PySide6:** retain `<6.12` on Python >=3.10 and the 6.9 range on Python 3.9.
  On 2026-10-24, official tracker issues
  [PYSIDE-3474](https://bugreports.qt.io/browse/PYSIDE-3474) and
  [PYSIDE-3424](https://bugreports.qt.io/browse/PYSIDE-3424) were still Reported,
  without a resolution or fix version. PyPI still supplied 6.12.0. Keep the
  native reference-count regression enabled and validate corrected wheels
  before removing this exclusion.
- **Python 3.15:** retain `<3.15` for this patch. Python 3.15 is released and
  RDKit supplies matching wheels, but the supported Qt series excludes that
  interpreter and VTK 9.7.1 has no CPython 3.15 wheel on PyPI. Do not advertise
  support before installation, full tests, and native rendering work on the
  supported operating systems.

Primary references checked on 2026-10-24:
[RDKit release](https://github.com/rdkit/rdkit/releases/tag/Release_2026_09_1),
[Qt 6.12 release notes](https://doc.qt.io/qtforpython-6/release_notes/pyside6_release_notes.html),
[Qt development notes](https://wiki.qt.io/Qt_for_Python_Development_Notes),
[Python releases](https://www.python.org/downloads/), and PyPI metadata for
[RDKit](https://pypi.org/project/rdkit/),
[PySide6-Essentials](https://pypi.org/project/PySide6-Essentials/), and
[VTK](https://pypi.org/project/vtk/).
