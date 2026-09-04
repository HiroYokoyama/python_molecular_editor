# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
* **Non-UTF-8 XYZ/MOL Import**: XYZ and MOL file imports now fall back through common non-UTF-8 encodings (UTF-8 BOM, Shift-JIS/CP932, EUC-JP) instead of failing outright with `UnicodeDecodeError` on files containing Japanese (or other non-UTF-8) comment/title text.

### Changed
* **Lint Toolchain & Type Annotations**: Added a per-module mypy override for `moleditpy.ui.custom_qt_interactor` to allow untyped `pyvistaqt` interactor subclassing without local errors ([#130](https://github.com/HiroYokoyama/python_molecular_editor/pull/130)).

### Documentation
* **AI Usage Documentation**: Added Gemini 3.7 Flash to Google AI models list in `docs/AI_USAGE.md` ([#131](https://github.com/HiroYokoyama/python_molecular_editor/pull/131)).
* **Contributing & Logging Guidelines**: Clarified log and traceback submission channels (console output, Error Dialog details, and optional `~/.moleditpy/moleditpy.log`) in `CONTRIBUTING.md` ([#131](https://github.com/HiroYokoyama/python_molecular_editor/pull/131)).
* **Test Inventory**: Documented previously unlisted test files across `tests/unit/README.md`, `tests/integration/README.md`, and `tests/e2e/README.md` ([#131](https://github.com/HiroYokoyama/python_molecular_editor/pull/131)).

---

## [4.9.0] - 2026-08-30

### Added
* **2D Cleanup Settings Tab**: Introduced a dedicated "2D Cleanup" settings tab exposing RDKit depiction layout options for `Clean Up 2D` (`Ctrl+J`), including:
  * Prefer CoordGen algorithm
  * Canonical orientation
  * Use ring templates
  * Straighten bonds
  * Avoid clashes ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).
* **CI Codecov Integration**: Added coverage uploading to Codecov in the GitHub Actions workflow ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).

### Changed
* **2D Coordinate Normalization**: `optimize_2d_coords` now normalizes bond lengths to a consistent target across both CoordGen (~1.0 Å) and the default generator (~1.5 Å) ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).
* **Default 2D Style Tuning**: Updated default 2D rendering parameters to heavier, clearer defaults (bond width 2.0 pt, atom font size 10 pt) ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).

### Fixed
* **Atom Font Initialization**: Fixed default atom font application upon atom creation in the 2D canvas ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).
* **RDKit Type Annotations**: Silenced platform-specific RDKit overload type warnings ([#129](https://github.com/HiroYokoyama/python_molecular_editor/pull/129)).

---

## [4.8.3] - 2026-08-29

### Changed
* **Lint Toolchain Alignment**: Pinned exact toolchain versions in `requirements-lint.txt` (`mypy`, `ruff`, `pylint`) and targeted mypy at Python 3.13 ([#128](https://github.com/HiroYokoyama/python_molecular_editor/pull/128)).
* **CI Linting Workflows**: Standardized CI checks to run `ruff format --check`, `ruff check`, and multi-platform `mypy` (win32, linux, darwin) alongside pylint ([#128](https://github.com/HiroYokoyama/python_molecular_editor/pull/128)).

### Fixed
* **PyInstaller Frozen Packaging**: Bundled `rdkit.ForceField` in frozen application distributions to ensure energy minimization and conformation generation work in standalone binaries ([#128](https://github.com/HiroYokoyama/python_molecular_editor/pull/128)).
* **Type Safety Cleanups**: Resolved type errors across core modules, eliminated unneeded `# type: ignore` directives, and refreshed coverage and pylint reports ([#128](https://github.com/HiroYokoyama/python_molecular_editor/pull/128)).

---

## [4.8.2] - 2026-08-28

### Fixed
* **Direct 2D-to-3D Coordinate Conversion**: Centred, rescaled, and de-stretched direct 2D-to-3D conversion conformers. Previously, canvas pixel coordinates mapped directly to 3D space, causing structures drawn off-centre to appear far from the 3D origin ([#126](https://github.com/HiroYokoyama/python_molecular_editor/pull/126)).
* **3D Bounding Box Accuracy**: Bounding box calculations in direct conversion now accurately reflect actual conformer geometry ([#126](https://github.com/HiroYokoyama/python_molecular_editor/pull/126)).

### Documentation
* Regenerated test assertion catalog and coverage reports ([#126](https://github.com/HiroYokoyama/python_molecular_editor/pull/126)).

---

## [4.8.1] - 2026-08-22

### Added
* **Plugin Menu Header Slot**: Enabled plugins to request placement in the top-level Plugin menu header slot ([#125](https://github.com/HiroYokoyama/python_molecular_editor/pull/125)).

### Changed
* **Valence-Aware Template Fusion**: Planned fused template bond orders in a single valence-aware pass, preventing fused ring placements from violating atomic valence limits ([#125](https://github.com/HiroYokoyama/python_molecular_editor/pull/125)).
* **Plugin Menu Lifecycle**: Improved plugin menu tear-down to properly destroy menu objects and prevent menu bar leaking on uninstalls or reloads ([#125](https://github.com/HiroYokoyama/python_molecular_editor/pull/125)).

### Fixed
* **Template Ghost Parity**: Synchronized template ghost preview graphics with actual placement geometry, bond trimming, and Kekulé rotations ([#125](https://github.com/HiroYokoyama/python_molecular_editor/pull/125)).
* **Menu Dividers and Clear All Zoom**: Fixed plugin menu separators, ordering inconsistencies, and viewport zoom resets on "Clear All" ([#125](https://github.com/HiroYokoyama/python_molecular_editor/pull/125)).

---

## [4.8.0] - 2026-08-21

### Added
* **About Image 3D Build**: Added shortcut to construct the 2D molecule structure in 3D directly from the About dialog image ([#124](https://github.com/HiroYokoyama/python_molecular_editor/pull/124)).

### Changed
* **Template Atom Replacement**: Clicking an existing atom when placing a user template now replaces it with the template's first atom (preserving element, charge, and radical) instead of retaining the original element ([#124](https://github.com/HiroYokoyama/python_molecular_editor/pull/124)).
* **Toolbar Mode Clarity**: Highlighted active template buttons across multiple toolbars and allowed exiting template mode via `Escape` ([#124](https://github.com/HiroYokoyama/python_molecular_editor/pull/124)).

### Fixed
* **Ghost Rendering Parity**: Rendered 2D template ghost previews using native `AtomItem` and `BondItem` widgets, matching real placement font settings and implicit hydrogen labels ([#124](https://github.com/HiroYokoyama/python_molecular_editor/pull/124)).
* **Test Isolation**: Isolated test suite executions from user configurations in `~/.moleditpy` ([#124](https://github.com/HiroYokoyama/python_molecular_editor/pull/124)).

---

## [4.7.2] - 2026-08-21

### Added
* **macOS Application Packaging**: Added automated build and publishing of macOS `.app` bundles in the release workflow ([#123](https://github.com/HiroYokoyama/python_molecular_editor/pull/123)).

### Fixed
* **OBJ/MTL Export UTF-8 Encoding**: Explicitly specified UTF-8 encoding when writing OBJ and MTL text files, preventing broken material references on non-UTF8 system locales (e.g. cp932) ([#123](https://github.com/HiroYokoyama/python_molecular_editor/pull/123)).
* **Export Error Handling**: Surfaced OBJ/MTL file write errors (read-only directories, permissions) via status bar alerts instead of unhandled exceptions ([#123](https://github.com/HiroYokoyama/python_molecular_editor/pull/123)).
* **Release Tag Validation**: Guarded release CI workflow from building missing tags ([#123](https://github.com/HiroYokoyama/python_molecular_editor/pull/123)).

---

## [4.7.1] - 2026-08-15

### Added
* **Plugin File Context API**: Added `PluginContext.set_current_file()` allowing plugins to update the active document filename after open/save operations ([#122](https://github.com/HiroYokoyama/python_molecular_editor/pull/122)).
* **Plugin Metadata**: Documented `PLUGIN_OPTIONAL_DEPENDENCIES` convention for plugin authors ([#122](https://github.com/HiroYokoyama/python_molecular_editor/pull/122)).

### Changed
* **Release Guard**: Prevented release packaging workflows when version strings disagree with `main` ([#122](https://github.com/HiroYokoyama/python_molecular_editor/pull/122)).
* **Zenodo Metadata**: Normalized legacy license identifiers prior to Zenodo metadata submission ([#122](https://github.com/HiroYokoyama/python_molecular_editor/pull/122)).

### Fixed
* **Linux E2E Package Name**: Resolved e2e test failures where `moleditpy` package name was hardcoded on Linux ([#122](https://github.com/HiroYokoyama/python_molecular_editor/pull/122)).

---

## [4.7.0] - 2026-08-11

### Changed
* **Zoom Gesture Performance**: Coalesced zoom re-indexing to execute once per gesture, improving responsiveness during pinch and wheel operations ([#121](https://github.com/HiroYokoyama/python_molecular_editor/pull/121)).
* **Code Cleanups**: Cleaned all remaining MyPy type issues and pruned redundant `hasattr` checks ([#121](https://github.com/HiroYokoyama/python_molecular_editor/pull/121)).

### Fixed
* **Headless .pmeraw Warning**: Removed blocking confirmation dialog on `.pmeraw` file loads under `MOLEDITPY_HEADLESS` ([#121](https://github.com/HiroYokoyama/python_molecular_editor/pull/121)).
* **Zoomed Hit Detection**: Tightened atom and bond hit-detection boxes at high zoom levels ([#121](https://github.com/HiroYokoyama/python_molecular_editor/pull/121)).

---

## [4.6.1] - 2026-08-08

### Changed
* **Chain Tool Refinements**: Refactored chain drawing logic for tighter endpoint snapping using `find_atom_near` ([#120](https://github.com/HiroYokoyama/python_molecular_editor/pull/120)).

### Fixed
* **Test Teardown Race**: Ensured full-GUI VTK plotter guard cleanup outranks `pytest-qt` teardown hooks ([#120](https://github.com/HiroYokoyama/python_molecular_editor/pull/120)).

---

## [4.6.0] - 2026-08-07

### Added
* **Interactive Chain Drawing Tool**: Added interactive continuous aliphatic chain tool with real-time cursor preview and automatic joining to hovered atoms ([#118](https://github.com/HiroYokoyama/python_molecular_editor/pull/118)).
* **Live Atom Count Badge**: Added real-time floating badge displaying the number of atoms currently being drawn by the chain tool ([#118](https://github.com/HiroYokoyama/python_molecular_editor/pull/118)).

---

## [4.5.2] - 2026-08-05

### Added
* **Full-GUI Test Tier**: Introduced `tests/gui/` full-GUI test suite executing real VTK/PyVista rendering and end-to-end benzene synthesis pipelines in headless CI ([#116](https://github.com/HiroYokoyama/python_molecular_editor/pull/116)).

### Changed
* **OS Theme Query Caching**: Cached system theme detection queries to avoid repetitive IPC latency ([#117](https://github.com/HiroYokoyama/python_molecular_editor/pull/117)).

### Fixed
* **First Heteroatom UI Stutter**: Resolved initial UI freeze when placing the first N/O heteroatom by warming up `QFontMetricsF` during startup rather than on the first canvas click ([#117](https://github.com/HiroYokoyama/python_molecular_editor/pull/117)).
* **CI Dev Branch Triggers**: Configured GitHub Actions CI workflows to run across all `dev-*` branch pushes ([#117](https://github.com/HiroYokoyama/python_molecular_editor/pull/117)).

---

## [4.5.1] - 2026-07-30

### Changed
* **Move Selected Atoms Dialog**: Improved atom transformation handling and expanded unit test coverage ([#113](https://github.com/HiroYokoyama/python_molecular_editor/pull/113)).

---

## [4.5.0] - 2026-07-28

### Added
* **Real-Time 2D Atom Dragging**: Atoms and connected bonds now update dynamically in real time during drag operations ([#112](https://github.com/HiroYokoyama/python_molecular_editor/pull/112)).
* **Large Molecule Drag Optimization**: Implemented automatic frame-skipping during real-time dragging for molecules exceeding 300 atoms ([#112](https://github.com/HiroYokoyama/python_molecular_editor/pull/112)).

### Fixed
* **Move Group Dialog Crash**: Resolved crash when opening Move Group dialog with preselected atoms ([#112](https://github.com/HiroYokoyama/python_molecular_editor/pull/112)).
* **Orphaned Drag Gestures**: Fixed dangling mouse drag states on abrupt focus losses ([#112](https://github.com/HiroYokoyama/python_molecular_editor/pull/112)).

---

## [4.4.1] - 2026-07-24

### Added
* **Mouse Sensitivity Setting**: Added configurable 3D mouse sensitivity slider in user settings ([#111](https://github.com/HiroYokoyama/python_molecular_editor/pull/111)).

### Documentation
* Removed deprecated measure mode references and updated user manual to version 4.4 ([#111](https://github.com/HiroYokoyama/python_molecular_editor/pull/111)).

---

## [4.4.0] - 2026-07-21

### Added
* **Interactor Style Watchdog**: Added 2-second background watchdog (`UIManager._check_interactor_style`) ensuring `CustomInteractorStyle` remains active throughout complex view operations ([#110](https://github.com/HiroYokoyama/python_molecular_editor/pull/110)).
* **Mouse State Self-Healing**: Automatically force-reset stuck rotation/pan flags when `QApplication.mouseButtons()` detects no buttons held ([#110](https://github.com/HiroYokoyama/python_molecular_editor/pull/110)).

### Changed
* **Rapid Click Handling**: Replaced click-swallowing filter with synthetic single-press redispatching, making rapid atom selection responsive and precise ([#110](https://github.com/HiroYokoyama/python_molecular_editor/pull/110)).

---

## [4.3.2] - 2026-07-17

### Added
* **Multi-Resolution Icons**: Added smaller application and document file icons for crisp taskbar and file manager display ([#108](https://github.com/HiroYokoyama/python_molecular_editor/pull/108), [#109](https://github.com/HiroYokoyama/python_molecular_editor/pull/109)).

---

## [4.3.1] - 2026-07-15

### Added
* **Corrupted-File Robustness Battery**: Added a 34-case automated test suite for empty, malformed, or truncated project and geometry files ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).

### Changed
* **Constrained Optimization Dialog**: Closing the dialog now cleanly cancels in-flight calculations; synchronized 3D constraints before undo stack pushes ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).

### Fixed
* **Extended Element Radii**: Fixed covalent-radius fallback dropping bonds for elements beyond vanadium (I, Br, Fe, etc.) during XYZ connectivity estimation ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).
* **OBJ Material Matching**: Fixed case-sensitive `.OBJ` extension causing MTL material references to be lost ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).
* **Stereo Undo Checkpoint**: Added missing undo state checkpoint on stereo wedge/hash direction flips (`W`/`D`) ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).
* **Startup File Exception Deferral**: Deferred startup file loading out of the `MainWindow` constructor so corrupted initial files do not crash the app window ([#107](https://github.com/HiroYokoyama/python_molecular_editor/pull/107)).

---

## [4.3.0] - 2026-07-12

### Added
* **Centralized Error & Exception Dialogs**: Introduced thread-safe error handler presenting unhandled exceptions and critical log events in a GUI details dialog ([#106](https://github.com/HiroYokoyama/python_molecular_editor/pull/106)).

### Changed
* **Camera Zoom Preservation**: Undo and redo actions now retain the user's active camera zoom instead of resetting the 3D viewport ([#106](https://github.com/HiroYokoyama/python_molecular_editor/pull/106)).

### Fixed
* **Plugin Action Group Leaks**: Tagged plugin actions to detach properly from `QActionGroup` parent containers on menu rebuild ([#106](https://github.com/HiroYokoyama/python_molecular_editor/pull/106)).

---

## [4.2.4] - 2026-07-11

### Fixed
* **Plugin Menu Cleanup**: Fixed `rebuild_plugin_menus` failing to clear stale plugin actions from standalone toolbars ([#105](https://github.com/HiroYokoyama/python_molecular_editor/pull/105)).
* **Menu Iteration Refactoring**: Extracted reusable menu iteration logic into `_get_plugin_target_menus()` ([#105](https://github.com/HiroYokoyama/python_molecular_editor/pull/105)).

---

## [4.2.3] - 2026-07-04

### Changed
* **PluginContext API Naming**: Renamed `PluginContext.fit_3d_view()` to `fit_2d_view()` to accurately represent 2D scene fitting, updating manual and unit tests ([#104](https://github.com/HiroYokoyama/python_molecular_editor/pull/104)).

---

## [4.2.2] - 2026-07-04

### Added
* **Plugin Optimizer Integration**: Allowed plugin-registered 3D optimizers to participate in 2D→3D conversion via background MMFF pre-optimization followed by plugin refinement ([#103](https://github.com/HiroYokoyama/python_molecular_editor/pull/103)).

### Fixed
* **Persisted Optimization Defaults**: Restored persisted plugin optimization methods at startup ([#103](https://github.com/HiroYokoyama/python_molecular_editor/pull/103)).
* **Menu Rebuild Idempotency**: Fixed plugin optimization actions disappearing from menus after plugin reload or uninstall ([#103](https://github.com/HiroYokoyama/python_molecular_editor/pull/103)).

---

## [4.2.1] - 2026-07-03

### Added
* **Dynamic Plugin Optimization UI**: Integrated plugin 3D optimization methods dynamically into both the top-level Optimization menu and right-click context menus ([#102](https://github.com/HiroYokoyama/python_molecular_editor/pull/102)).

---

## [4.2.0] - 2026-07-02

### Changed
* **Python 3.9 Compatibility**: Pinned `pyqt6 < 6.10` on Python 3.9 platforms ([#101](https://github.com/HiroYokoyama/python_molecular_editor/pull/101)).

### Fixed
* **macOS Dark Mode Support**: Improved OS dark-mode detection and theme synchronization on macOS ([#101](https://github.com/HiroYokoyama/python_molecular_editor/pull/101)).
* **Charge Validation & IO Robustness**: Resolved charge calculation anomalies and edge cases in `IOManager` ([#101](https://github.com/HiroYokoyama/python_molecular_editor/pull/101)).

---

## [4.1.4] - 2026-06-30

### Added
* **Headless E2E Test Suite**: Added comprehensive `tests/e2e/` test suite validating the complete 2D→3D pipeline with real RDKit (no mock chemistry) ([#96](https://github.com/HiroYokoyama/python_molecular_editor/pull/96)).
* **Windows Installer Upgrades**: Added automatic uninstall of previous installations and prompt to close running instances during setup ([#96](https://github.com/HiroYokoyama/python_molecular_editor/pull/96)).

### Fixed
* **Plugin Hot Reloading**: Plugin hot-reload now updates toolbars and menus immediately and purges stale `sys.modules` cache entries ([#96](https://github.com/HiroYokoyama/python_molecular_editor/pull/96)).

---

## [4.1.3] - 2026-06-29

### Added
* **Box Selection in Move Group**: Added rectangle box selection mode to `MoveSelectedAtomsDialog` via PyVista picking ([#95](https://github.com/HiroYokoyama/python_molecular_editor/pull/95)).

### Changed
* **Single-Click Reset**: Single clicks in Box Selection mode now clear active atom selections ([#95](https://github.com/HiroYokoyama/python_molecular_editor/pull/95)).

---

## [4.1.2] - 2026-06-26

### Fixed
* Addressed UI event handling issues and stabilized 3D interactor state transitions ([#94](https://github.com/HiroYokoyama/python_molecular_editor/pull/94)).

---

## [4.1.1] - 2026-06-24

### Fixed
* Fixed plugin manager dialog initialization and 3D viewport update synchronization ([#93](https://github.com/HiroYokoyama/python_molecular_editor/pull/93)).

---

## [4.1.0] - 2026-06-24

### Added
* **PluginContext Extension**: Added canvas, viewport, and structure helper methods to `PluginContext` ([#92](https://github.com/HiroYokoyama/python_molecular_editor/pull/92)).

---

## [4.0.2] - 2026-06-21

### Fixed
* Resolved runtime exception handling in template preview rendering loop.

---

## [4.0.1] - 2026-06-20

### Fixed
* Fixed Python 3.9 compatibility issues in plugin menu manager ([#90](https://github.com/HiroYokoyama/python_molecular_editor/pull/90)).

---

## [4.0.0] - 2026-06-20

### Added
* **Plugin Architecture V4**: Major release upgrading the plugin interface to V4 (`PluginContext` V4 API):
  * `load_from_smiles()`
  * `to_xyz_block()`
  * Dynamic canvas and view manipulation helpers
  * Full type annotations across plugin contracts ([#89](https://github.com/HiroYokoyama/python_molecular_editor/pull/89)).
* **Plugin Documentation V4**: Published `docs/PLUGIN_DEVELOPMENT_MANUAL_V4.md` and archived V3 documentation ([#89](https://github.com/HiroYokoyama/python_molecular_editor/pull/89)).

### Changed
* **Code Quality & Linting**: Eliminated all unused imports and flake8 warnings across `moleditpy/src/` ([#89](https://github.com/HiroYokoyama/python_molecular_editor/pull/89)).

---

## [3.6.6] - 2026-06-11

### Changed
* **XYZ Coordinate Formatting**: Standardized exported XYZ coordinates into aligned columns with 8 decimal places (`symbol:<5` and `coords:>15.8f`) ([#84](https://github.com/HiroYokoyama/python_molecular_editor/pull/84), [#85](https://github.com/HiroYokoyama/python_molecular_editor/pull/85)).

---

## [3.6.4] - 2026-06-02

### Added
* **Plugin Metadata Specifications**: Documented `PLUGIN_TAGS`, `PLUGIN_SUPPORTED_MOLEDITPY_VERSION`, and optional dependencies in plugin manuals ([#83](https://github.com/HiroYokoyama/python_molecular_editor/pull/83), [#82](https://github.com/HiroYokoyama/python_molecular_editor/pull/82)).

---

## [3.6.0] - 2026-05-30

### Added
* **Modeless Dialog Architecture**: Enhanced non-blocking dialog management across UI components ([#75](https://github.com/HiroYokoyama/python_molecular_editor/pull/75)).
* **Manual CI Dispatch**: Added `workflow_dispatch` trigger to GitHub Actions tests workflow ([#75](https://github.com/HiroYokoyama/python_molecular_editor/pull/75)).

### Fixed
* **Windows Taskbar Icon**: Forced window icon refresh after initialization on Windows ([#75](https://github.com/HiroYokoyama/python_molecular_editor/pull/75)).
* **PyInstaller Asset Bundling**: Corrected asset bundling path for Linux packaging ([#75](https://github.com/HiroYokoyama/python_molecular_editor/pull/75)).
* **Thread Flakiness in Tests**: Mocked `QThread.start` in dummy compute tests to eliminate background segfaults ([#75](https://github.com/HiroYokoyama/python_molecular_editor/pull/75)).

---

## [3.5.0] - 2026-05-29

### Fixed
* **3D Export Keys and Colors**: Resolved `KeyError` during 3D scene exports and corrected atom color mapping ([#72](https://github.com/HiroYokoyama/python_molecular_editor/pull/72)).

---

## [3.4.0] - 2026-05-10

### Changed
* **Dependency Upgrades**: Updated RDKit dependency constraints to `<2026.4` across Linux and Windows configurations ([#68](https://github.com/HiroYokoyama/python_molecular_editor/pull/68), [#69](https://github.com/HiroYokoyama/python_molecular_editor/pull/69), [#70](https://github.com/HiroYokoyama/python_molecular_editor/pull/70)).

---

## [3.3.0] - 2026-05-05

### Fixed
* **3D Atom Picking Deadlock**: Avoided VTK cell picking during atom clicks to prevent deadlocks and interactor freezing ([#65](https://github.com/HiroYokoyama/python_molecular_editor/pull/65)).
* **3D Export Coloring**: Fixed missing color attributes when exporting 3D geometries ([#65](https://github.com/HiroYokoyama/python_molecular_editor/pull/65)).

---

## [3.2.0] - 2026-04-25

### Added
* **Plugin Manual V3**: Added comprehensive Plugin Development Manual V3 and expanded test suites ([#58](https://github.com/HiroYokoyama/python_molecular_editor/pull/58), [#60](https://github.com/HiroYokoyama/python_molecular_editor/pull/60), [#61](https://github.com/HiroYokoyama/python_molecular_editor/pull/61), [#63](https://github.com/HiroYokoyama/python_molecular_editor/pull/63)).

---

## [3.1.0] - 2026-04-14

### Added
* **Non-Destructive 2D Imports**: SMILES, InChI, and MOL/SDF 2D imports now append structures alongside existing molecules (placed with an 80 px offset) rather than overwriting the canvas ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).
* **Import Undo Support**: Pushed 2D structure imports to the undo stack and marked unsaved state ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).

### Changed
* **String Importer Refactoring**: Extracted shared placement helpers in `StringImporterManager` to unify SMILES and InChI code paths ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).

### Fixed
* **Close Unsaved Changes Crash**: Fixed `AttributeError` crash when exiting with unsaved changes ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).
* **3D Projection Mode Setting**: Fixed perspective/orthographic setting defaulting incorrectly to perspective on reload ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).
* **Template Unsaved Flag**: Fixed saving user templates (`.pmetmplt`) inadvertently clearing project dirty flags ([#57](https://github.com/HiroYokoyama/python_molecular_editor/pull/57)).

---

## [3.0.6] - 2026-04-13

### Added
* **Headless Plugin Installation**: Added command-line option for unattended plugin installation ([#56](https://github.com/HiroYokoyama/python_molecular_editor/pull/56)).

### Fixed
* Fixed 3D viewport zoom resetting unexpectedly during 3D edit tool operations ([#56](https://github.com/HiroYokoyama/python_molecular_editor/pull/56)).

---

## [3.0.5] - 2026-04-11

### Added
* **CLI Version Argument**: Added `--version` CLI flag to output application version ([#55](https://github.com/HiroYokoyama/python_molecular_editor/pull/55)).

### Fixed
* Fixed file launching crashes when operating in safe mode ([#55](https://github.com/HiroYokoyama/python_molecular_editor/pull/55)).

---

## [3.0.4] - 2026-04-09

### Fixed
* Fixed camera orientation resetting after loading `.pemprj` projects or changing initial display styles ([#54](https://github.com/HiroYokoyama/python_molecular_editor/pull/54)).

---

## [3.0.3] - 2026-04-09

### Fixed
* Fixed translate and move dialog coordinate transformation logic ([#53](https://github.com/HiroYokoyama/python_molecular_editor/pull/53)).

---

## [3.0.0] - 2026-04-02

### Added
* **Plugin Architecture V3**: Major architectural overhaul introducing the V3 plugin system with decoupled UI menus and life-cycle handlers ([#49](https://github.com/HiroYokoyama/python_molecular_editor/pull/49)).
* **Automated CI Workflows**: Established unified unit and integration testing pipelines on GitHub Actions ([#49](https://github.com/HiroYokoyama/python_molecular_editor/pull/49)).

---

## [2.8.0] - 2026-03-03

### Added
* **Interactive Angle & Bond Length Sliders**: Added real-time slider controls with automatic angle range wrapping (`[-180°, 180°]`) and stable rotation axes during slider drags in geometry transformation dialogs.
* **Folder Comparison Tooling**: Added PowerShell comparison tooling for Linux/Windows package syncing.

### Changed
* **Windows Installer Build**: Updated InnoSetup script and distribution packaging.

---

## [2.7.0] - 2026-02-26

### Added
* **Developer Documentation Suite**: Published comprehensive architecture documentation across `docs/`:
  * `docs/ARCHITECTURE.md`
  * `docs/CORE_DATA_STRUCTURES.md`
  * `docs/DIALOGS_AND_UI.md`
  * `docs/ANALYSIS_AND_CALCULATIONS.md`
* **Static Analysis Tooling**: Integrated `.pylintrc` and code quality enforcement.

### Fixed
* Resolved multi-fragment calculation crashes and calculation worker thread synchronization ([#37](https://github.com/HiroYokoyama/python_molecular_editor/pull/37)).

---

## [2.6.3] - 2026-02-22

### Changed
* **Python 3.14 Compatibility**: Updated dependency versions and package build requirements to support Python 3.14 ([#33](https://github.com/HiroYokoyama/python_molecular_editor/pull/33)).

---

## [2.6.2] - 2026-02-22

### Changed
* **Calculation Engine Refactoring**: Refactored calculation worker and compute logic for background thread safety ([#32](https://github.com/HiroYokoyama/python_molecular_editor/pull/32)).
* **Chiral Tests**: Enabled chiral stereochemistry validation tests on Python 3.14 ([#32](https://github.com/HiroYokoyama/python_molecular_editor/pull/32)).

---

## [2.6.0] - 2026-02-15

### Added
* **Dedicated Geometry & Alignment Dialogs**:
  * `AlignmentDialog` / `AlignPlaneDialog`: Align structures to specific coordinate axes or planes.
  * `AngleDialog` & `BondLengthDialog`: Modify bond lengths, angles, and dihedrals with numerical precision.
  * `ColorSettingsDialog`: Interactive CPK color customizer.
* **Stereochemistry and Geometry Tests**: Added comprehensive unit tests for E/Z stereochemistry, MOL block parsing priorities, and atom item implicit hydrogen rendering ([#28](https://github.com/HiroYokoyama/python_molecular_editor/pull/28), [#29](https://github.com/HiroYokoyama/python_molecular_editor/pull/29), [#30](https://github.com/HiroYokoyama/python_molecular_editor/pull/30)).
* **Test Utility Suite**: Added `tests/utils/generate_assertion_catalog.py` and coverage report generator `tests/utils/print_cov.py`.

### Changed
* **Package Structure**: Moved application assets to the package root for uniform wheel and frozen distribution packaging.

### Fixed
* **Windows Taskbar Icon**: Set explicit application-level `.ico` icon for Windows taskbar consistency.

---

## [2.5.0] - 2026-02-01

### Added
* **Dynamic 2D Atom Styling**: Added support for dynamic font family and font size adjustments on 2D `AtomItem` and `BondItem` widgets.
* **CPK / Bond Color Toggles**: Added settings to render 2D atom text in uniform bond colors or CPK element colors.

---

## [2.4.0] - 2026-01-16

### Added
* **Plugin Document Reset Handlers**: Allowed plugins to register callbacks for document reset and file opening events.
* **Plugin File Openers**: Added priority-sorted plugin file opening registrations with custom keyboard shortcuts.

---

## [2.2.0] - 2025-12-25

### Added
* **Design Principles**: Published core architectural principles and guidelines in `DESIGN_PRINCIPLES.md` ([#23](https://github.com/HiroYokoyama/python_molecular_editor/pull/23)).
* **Modularization**: Expanded modular component separation between 2D canvas, 3D viewport, and calculation backends.

---

## [2.1.0] - 2025-12-14

### Added
* **Installer Documentation**: Added comprehensive Windows installer documentation in English and Japanese.
* **Icon Assets**: Added high-resolution `.ico` and `.icns` icons and icon generator scripts.

---

## [2.0.0] - 2025-12-12

### Added
* **Major Version 2.0.0 Release**:
  * **Aromatic Ring Display Styles**: Added support for displaying aromatic rings as single bonds, toroids (rings), or alternating Kekulé structures.
  * **3D View Settings Dialog**: Configurable atom sphere radii, bond cylinder thickness, rendering quality, and multiple-bond spacing.
  * **Customizable CPK Colors**: Interactive CPK color customizer.
  * **Plugin System V1**: Initial plugin system supporting external tools (such as MS spectrum simulation and Hello World demo).
  * **Interactive 3D Dragging**: Comprehensive geometric manipulation directly in the 3D PyVista viewport.
  * **Constrained Force Field Optimization**: MMFF94 and UFF optimization under fixed distance, angle, and dihedral constraints.

---

## [1.18.0] - 2025-11-22

### Added
* **GPL-v3 Licensing**: Officially transitioned to the GNU General Public License v3.0 (GPL-v3).
* **Installer CLI**: Introduced `moleditpy-installer` CLI command to create desktop and start menu shortcuts across Windows, macOS, and Linux.

---

## [1.16.0] - 2025-11-18

### Changed
* **Package Modularization**: Major refactoring transitioning the monolithic source code into the modular `moleditpy/` subpackage structure ([#22](https://github.com/HiroYokoyama/python_molecular_editor/pull/22)).

---

## [1.14.0] - 2025-11-15

### Added
* **3D Geometric Tools Suite**:
  * **Translation**: Move molecules or selected atom subsets to exact coordinates.
  * **Planarization**: Flatten selections onto XY, XZ, YZ, or best-fit planes.
  * **Mirror Reflection**: Invert molecular geometry across principal planes.
  * **Numerical Geometry Control**: Set exact bond lengths, bond angles, and torsional dihedral angles.
* **Native Project File Format (`.pmeprj`)**: Save and load complete sessions including 2D sketch, 3D conformers, and active geometric constraints.
* **Linux Specific Package (`moleditpy-linux`)**: Created standalone Linux package configuration to resolve Open Babel library collisions.

---

## [1.9.0] - 2025-10-14

### Added
* **Ring Template System**: Instant placement of 3- to 9-membered rings and benzene with live preview and automatic valence-aware Kekulé adjustment upon ring fusion.
* **Formal Charges & Radicals**: Added support for setting atomic formal charges (`+`/`-`) and radical states (`.`) via toolbar and keyboard shortcuts.
* **Periodic Table Dialog**: Full periodic table element selector.
* **Clipboard Support**: Full `Cut` (`Ctrl+X`), `Copy` (`Ctrl+C`), and `Paste` (`Ctrl+V`) for molecular fragments.

---

## [1.2.0] - 2025-10-08

### Added
* **Cross-Platform Installation Guides**: Added Linux installation notes with `conda-forge` environment instructions and dependency pinning.

---

## [1.0.0] - 2025-10-07

### Added
* **Initial Stable Release (v1.0.0)**:
  * 2D molecular drawing on interactive `QGraphicsScene` canvas with `AtomItem` and `BondItem`.
  * 3D molecular visualization using PyVista and VTK.
  * RDKit conformer generation and 2D-to-3D structure conversion.
  * File import and export for standard chemical formats (MOL, SDF, XYZ).
  * State serialization with full Undo/Redo history.

---

## [0.3.0] - 2025-10-04

### Added
* Early prototype 2D chemical drawing canvas, preliminary 3D viewer integration, and initial project README.

---

## [0.1.0] - 2025-10-01

### Added
* Initial repository creation and proof-of-concept molecular editor.
