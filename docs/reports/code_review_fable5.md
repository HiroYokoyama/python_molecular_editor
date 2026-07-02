# MoleditPy Package Code Review (fable5)

- **Date:** 2026-07-02
- **Reviewer:** Claude Code (Fable 5)
- **Scope:** `python_molecular_editor/moleditpy/src/moleditpy/` (~29.6k lines). Full read of `core/`, `plugins/`, `utils/`, `main.py`, `app_state.py`, `io_logic.py`, plus targeted reads of `calculation_worker.py`, `main_window.py` and pattern sweeps across `ui/`. Docs cross-checked against code.

Severity: 🔴 bug (wrong behavior possible), 🟡 latent bug / robustness, 🔵 inconsistency / stale doc / cosmetic.

---

## 1. Core data model (`core/`)

### 🔴 1.1 `MolecularData.add_bond` can create duplicate bonds under reversed keys
`core/molecular_data.py:82-106`. Stereo bonds keep caller key order `(id1, id2)`; non-stereo keys are sorted. The `is_new_bond` check correctly considers both key directions, but the update/insert step does **not**:

```python
if (id1, id2) in self.bonds:
    self.bonds[(id1, id2)].update(bond_data)
else:
    self.bonds[(id1, id2)] = bond_data   # reversed key (id2, id1) may still exist!
```

If a bond exists as `(2, 5)` and `add_bond(5, 2, stereo=1)` is called (or a stereo bond `(5, 2)` exists and `add_bond(2, 5, stereo=0)` sorts to `(2, 5)`), **two entries for the same atom pair** end up in `self.bonds`, while adjacency is not duplicated — the model is now internally inconsistent, and `to_mol_block()` / `to_rdkit_mol()` would emit/attempt the bond twice (RDKit `AddBond` on a duplicate raises → sanitize failure → silent fallback).

The 2D scene code works around this by always calling `remove_bond` before `add_bond` (`ui/molecule_scene.py:654-671`), but the data-layer API itself (exposed to plugins and tests) has no guard. Fix: when the exact key is absent, also check the reversed key and update/replace it (popping the old key if the direction changed).

### 🔴 1.2 `to_mol_block()` manual fallback breaks on float bond order and missing positions
`core/molecular_data.py:339-390` (the fallback used when RDKit sanitization fails):
- `f"{idx1:3d}{idx2:3d}{order:3d}..."` — `order` is typed `Union[int, float]` and `to_rdkit_mol` explicitly supports `1.5` (aromatic). A float raises `ValueError: Unknown format code 'd'` — and this fallback path is *not* wrapped in try/except, so the caller crashes exactly in the case the fallback was built for.
- Atoms with no/invalid `pos` are `continue`-skipped from the atom block, but `num_atoms` in the counts line and `atom_map` indices still include them → the emitted MOL block has fewer atom lines than declared and bond indices pointing at shifted/nonexistent atoms.

Fix: cast `int(round(order))` (or map 1.5 → 4 "aromatic" per V2000) and build `atom_map` only from atoms actually written.

### 🟡 1.3 `resolve_2d_overlaps` recursive union-find
`core/mol_geometry.py:620-624`. `find_set` uses recursive path compression; for very large flat overlap chains this can hit Python's recursion limit. Iterative compression is a trivial fix. Low likelihood, but this is the "pure logic" layer meant to be robust.

---

## 2. Plugin system (`plugins/`)

### 🔵 2.1 Docstring vs. implementation: `get_selected_atom_indices`
`plugins/plugin_interface.py:134-139` says *"currently selected in the 2D **or 3D** view"*. The implementation (`plugin_manager.py:576-618`) only reads the 2D scene selection. Either implement 3D selection or fix the docstring — plugin authors will rely on this.

### 🔵 2.2 Dead/misleading `except ImportError` in `get_selected_atom_indices`
`plugins/plugin_manager.py:611-614`:
```python
except (ImportError):  # [OPTIONAL DEP] importlib.metadata unavailable (<3.8); silently skip.
    pass
```
Nothing in the `try` block imports anything; the comment about `importlib.metadata` belongs to some other code. Dead handler + stale comment — remove.

### 🟡 2.3 `install_plugin` ZIP "nested" heuristic misfires on single-file ZIPs
`plugins/plugin_manager.py:150-184`. A ZIP whose only top-level entry is one file (e.g. `myplugin.py`) yields `roots == {"myplugin.py"}` → `is_nested=True` → extracted directly into the plugin root instead of a wrapper folder (Case B). Only `__init__.py` is special-cased. Works by accident for plain `.py` plugins, but a ZIP containing a single non-`__init__` data or py file bypasses the intended flat-ZIP handling. Guard should be "the single root is a directory", e.g. check `name.endswith('/')` or use `zipfile.Path`.

### 🟡 2.4 `discover_plugins` doesn't skip dot-directories
`plugins/plugin_manager.py:251-253`. Comment says "Exclude hidden directories" but only `__*` is excluded; `.git`, `.vscode`, etc. inside the plugin dir are scanned, and any `.py` inside them is loaded/executed as a plugin.

### 🔵 2.5 `mark_project_modified` violates the project's own error-handling standard
`plugins/plugin_interface.py:394-400` uses `except Exception: pass`. CLAUDE.md says "never bare pass" and ARCHITECTURE.md advertises a "Never Hide Errors" policy; there are several other `contextlib.suppress(Exception)` / commented-`pass` blocks in `ui/` too. Either relax the documented policy or add `logging.debug(..., exc_info=True)` here.

---

## 3. I/O layer (`ui/io_logic.py`, `ui/app_state.py`)

### 🔴 3.1 `save_3d_as_mol` never appends the `.mol` extension
`ui/io_logic.py:877-901`. Every other save function (`save_as_mol`, `save_as_xyz`, `save_project_as`, `save_raw_data`, `save_as_json`) appends the extension when missing; `save_3d_as_mol` does not — typing a bare name saves an extensionless file that the app's own openers then won't match by extension.

### 🟡 3.2 `.pmeraw` uses `pickle.load` on user-supplied files
`ui/io_logic.py:989-990`. Opening an untrusted `.pmeraw` project executes arbitrary code by design of pickle. At minimum this deserves a warning in docs/README ("only open .pmeraw files you created"); better, migrate the format to the JSON `.pmeprj` path entirely (it already exists) and deprecate raw pickle load.

### 🟡 3.3 `prompt_for_charge` silently swallows invalid input
`ui/io_logic.py:319-323`. Typing garbage into the charge dialog returns `(0, True, False)` with no feedback — the user believes their charge was applied. Should re-prompt or show a validation message.

### 🟡 3.4 `_process` inner try/except is a no-op
`ui/io_logic.py:170-181`: `except (RuntimeError, ValueError, TypeError): raise` — catch-and-reraise with nothing else; remove.

### 🔵 3.5 `load_json_data` / `load_raw_data` skip the unsaved-changes check
`open_project_file` gates on `check_unsaved_changes()`, but calling `load_json_data(path)` / `load_raw_data(path)` directly (menu actions, plugin file openers, drag-drop) goes straight to `clear_all(skip_check=True)` — unsaved work is discarded without prompting on those paths.

### 🔵 3.6 `save_project` (Ctrl+S) doesn't refresh the undo snapshot
`ui/io_logic.py:361-397` — `save_project_as`/`save_raw_data`/`save_as_json` all call `host.save_state_snapshot()` after a successful save; the plain `save_project` path does not. Whatever "compare against saved state" logic uses `saved_state` will diverge after a Ctrl+S save.

---

## 4. Calculation worker (`ui/calculation_worker.py`)

### 🟡 4.1 Null-check happens after the call that would crash
Line 657-659:
```python
rd_mol = Chem.AddHs(Chem.MolFromMolBlock(out_mol_block, removeHs=False))
if not rd_mol:
    raise ValueError("Open Babel produced invalid MOL.")
```
If `MolFromMolBlock` returns `None`, `Chem.AddHs(None)` raises a Boost `ArgumentError` first, so the intended "invalid MOL" error message can never be produced. Split the two calls and check between them.

### 🔵 4.2 No-op line in the embedded Open Babel subprocess script
Line 631: `if "BABEL_DATADIR" in os.environ: os.environ["BABEL_DATADIR"] = os.environ["BABEL_DATADIR"]` assigns a variable to itself. Child processes inherit the environment anyway — delete the line (and the "Inherit babel variables" comment).

---

## 5. Misc application code

### 🟡 5.1 macOS dark-mode detection is a stub
`utils/system_utils.py:48-50` — on Darwin the function unconditionally returns `"light"`, so macOS users never get automatic dark theme. `defaults read -g AppleInterfaceStyle` (returns `Dark`/error) would implement it in the same subprocess style already used for GNOME.

### 🔵 5.2 `VDW_RADII` are not van der Waals radii
`utils/constants.py:208` — `pt.GetRvdw(i) * 0.3` bakes a display scale factor into a constant named `VDW_RADII`. Anyone reusing this for chemistry (e.g. a plugin) gets radii 3.3× too small. Rename (e.g. `VDW_DISPLAY_RADII`) or store the scale separately.

### 🔵 5.3 Hardcoded AppUserModelID version
`main.py:112` — `myappid = "hyoko.moleditpy.1.0"` while the app is at 4.x. Harmless for taskbar grouping but stale; consider deriving from `VERSION` major.

---

## 6. Documentation / config staleness

| # | Where | Issue |
|---|-------|-------|
| 🔵 6.1 | `python_molecular_editor/CLAUDE.md` | Says "MoleditPy … (v3.6.0)"; `pyproject.toml` is at **4.1.4**. |
| 🔵 6.2 | `DEV_MAIN/CLAUDE.md` (workspace root) | Refers to `context.register_reset_handler`; the real API is `register_document_reset_handler` (interface + V4 manual agree). |
| 🔵 6.3 | `docs/ARCHITECTURE.md` "Never Hide Errors" section | Claims every defensive `try/except` includes `logging.error()` or re-raise; the codebase contains many `contextlib.suppress(...)` and commented `pass` blocks (see 2.5). The policy text overstates current practice. |
| 🔵 6.4 | `plugins/plugin_manager.py:611` comment | Stale importlib.metadata comment (see 2.2). |

Checked and **found consistent** (no action needed): all `context.*` methods in `PLUGIN_DEVELOPMENT_MANUAL_V4.md` exist in `PluginContext`; all cross-manager attribute references used by the plugin API (`view_3d_manager`, `init_manager.scene`, `ui_manager.enter_3d_viewer_mode`, `compute_manager.check_chemistry_problems_fallback`, `MainWindow.scene/plotter` proxies, `string_importer_manager.load_from_smiles`, color-override methods in `view_3d_logic.py`) resolve to real definitions.

---

## Suggested priority

1. **1.1 / 1.2** — data-model correctness (duplicate bonds, broken MOL fallback); these are in the "pure, tested core" and cheap to fix + unit-test.
2. **3.1, 4.1** — small user-facing/robustness fixes, one-liners.
3. **2.3, 2.4, 3.2** — plugin-install and pickle hardening.
4. **Docs**: bump CLAUDE.md version string, fix `register_reset_handler` name in the workspace CLAUDE.md, tone down or enforce the "Never Hide Errors" claim.

*Not covered in depth:* `view_3d_logic.py`, `main_window_init.py`, `molecular_scene_handler.py`, `edit_actions_logic.py` internals (≈7k lines of Qt/VTK UI code) beyond pattern sweeps and cross-reference checks — a follow-up pass on those files could be worthwhile.

---

# Follow-up pass (2026-07-02, second review of remaining UI modules)

Scope: `compute_logic.py`, `string_importers.py`, `edit_actions_logic.py` (history core), `export_logic.py` (extension/save paths), plus a late-binding-lambda scan across all of `ui/` (clean — all loop lambdas correctly bind via default args).

## Fixed on dev-4.2.0

- 🔴 **`compute_logic._remove_calculating_text`**: accessed `plotter.renderer` in the `if` condition *before* the `contextlib.suppress` block → `AttributeError` crash when the plotter is None (headless / 3D disabled). Also had a duplicated `and plotter.renderer and plotter.renderer` condition, and cleared the actor attribute from `host.__dict__` although it is stored on the ComputeManager itself (stale actor reference survived). All three fixed.
- 🟡 **`compute_logic.trigger_conversion` / `on_calculation_finished`**: unguarded `plotter.add_text(...)` / `plotter.reset_camera()` — same None-plotter crash family. Guarded.
- 🟡 **`io_logic` deferred QTimer lambdas** (4 sites): `lambda: ...plotter.view_isometric()` / `.render()` fire 50–100 ms later with no None check — if the plotter is gone by then, the exception lands in the Qt event loop. Replaced with guarded `_plotter_view_isometric` / `_plotter_render` helpers.
- 🔵 **`compute_logic._handle_chemistry_problems`**: `test_trigger_conversion_chemistry_problems` fails when run in isolation (pre-existing, order-dependent — a mocked host makes `QMessageBox.critical` raise `TypeError`). This is a test-fixture issue, not a production bug; production code intentionally left unchanged (no test-only suppression added). **Fixed in the test** (commit `762ad0b`): it now patches `QMessageBox.critical` like its sibling tests.

## Noted, not changed

- 🔵 `edit_actions_logic.push_undo_state`: the duplicate-state comparison covers atoms/bonds/3D mol but not `constraints_3d` / viewer-mode — an edit that only changes constraints does not create an undo entry.
- 🔵 The undo stack is uncapped; each entry deep-copies the full state including the 3D binary. Long sessions with large molecules grow memory without bound. Consider a max depth (e.g. 100).
- 🔵 `string_importers.load_from_smiles` / `load_from_inchi` are near-identical ~35-line bodies — could share one `_load_from_mol(mol, label)` helper.
- ✅ Export extension handling (`.stl/.obj/.png/.svg`) is consistent; undo/redo stack logic itself is sound; no late-binding lambda captures anywhere in `ui/`.
