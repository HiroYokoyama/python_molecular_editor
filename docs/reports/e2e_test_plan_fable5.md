# E2E Test Expansion Plan (fable5) — executable by a smaller model

**Goal:** extend `python_molecular_editor/tests/e2e/` from 14 ethane-conversion tests to full user-workflow coverage (project round-trip, import/export, undo/redo, editing, plugins).

**Ground rules (read first, follow always):**
1. Work ONLY in `python_molecular_editor/tests/e2e/`. NEVER modify `moleditpy/src/**` or `moleditpy-linux/**`. If a test seems to need a production change, STOP and report instead.
2. After EVERY new test file: run it, then run the whole e2e suite, then commit. Never batch multiple files into one commit.
3. Reuse the existing fixtures — do not write a new MainWindow fixture. `tests/e2e/conftest.py` provides:
   - `app` (session QApplication), `window` (real MainWindow, headless, VTK/PyVista mocked, dummy PluginManager), `qtbot`.
   - Package indirection: import the app package via `from conftest import _PKG` style used in `test_ethane_gui.py` (copy its import pattern exactly — it selects `moleditpy` vs `moleditpy_linux`).
4. Mark every test that uses `window` with `@pytest.mark.gui` (see `pytest.ini`, `--strict-markers` is on).
5. All modal dialogs MUST be monkeypatched before triggering them, or the test hangs forever:
   - `QFileDialog.getSaveFileName` / `getOpenFileName` → return `(str(tmp_path/"name.ext"), None)`
   - `QMessageBox.question` → return `QMessageBox.StandardButton.No` (or Yes when the flow needs it)
   - `QMessageBox.critical` / `warning` / `information` → `lambda *a, **k: None`
   - Patch at the module that *uses* it, e.g. `monkeypatch.setattr(f"{PKG}.ui.io_logic.QFileDialog.getSaveFileName", ...)`.
6. Run command (from `python_molecular_editor/`, Windows Git Bash):
   `MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python -m pytest tests/e2e -q -p no:cacheprovider`
7. Commit message format: one line `e2e: add <area> workflow tests`, then the trailer lines used by earlier commits on this branch (Co-Authored-By + Claude-Session).
8. Useful internal APIs (all real, no mocks needed):
   - Draw atoms/bonds: `window.scene.create_atom("C", QPointF(x, y))` → returns atom_id; `window.scene.create_bond(window.scene.atom_items[id1], window.scene.atom_items[id2], bond_order=1, bond_stereo=0)`
   - Data model: `window.data.atoms`, `window.data.bonds`, `window.data.next_atom_id`
   - State: `window.state_manager.create_json_data()`, `.load_from_json_data(d)`, `.get_current_state()`, `.set_state_from_data(d)`
   - Undo: `window.edit_actions_manager.push_undo_state()`, `.undo()`, `.redo()`, `.undo_stack`, `.redo_stack`
   - IO: `window.io_manager.save_project()`, `.open_project_file(path)`, `.load_mol_file(path)`, `.load_xyz_file(path)`, `.save_as_mol()`, `.save_3d_as_mol()`, `.save_as_xyz()`
   - Import: `window.string_importer_manager.load_from_smiles("CCO")`
   - How ethane is drawn + converted: copy `_draw_ethane(window)` and the conversion-wait pattern (`qtbot.waitUntil(lambda: window.current_mol is not None, timeout=...)`) from `test_ethane_gui.py` verbatim.

---

## Step 0 — Baseline (do before anything)

```bash
cd python_molecular_editor
MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python -m pytest tests/e2e -q -p no:cacheprovider
```
Expected: all existing tests pass. If not, STOP and report.

## Step 1 — `tests/e2e/test_project_roundtrip.py`

Purpose: the single most valuable missing E2E — save a project, reload it, verify nothing is lost.

Tests (each `@pytest.mark.gui`):
1. `test_pmeprj_roundtrip(window, qtbot, tmp_path, monkeypatch)`
   - Draw ethane (copy `_draw_ethane`). Record `atoms_before = {id: (d["symbol"], d["charge"]) for ...}`, bond count, `next_atom_id`.
   - `json_data = window.state_manager.create_json_data()`; write with `json.dump` to `tmp_path/"t.pmeprj"`.
   - Clear: `window.edit_actions_manager.clear_all(skip_check=True)`; assert `window.data.atoms == {}`.
   - Load via `window.io_manager.open_project_file(str(path))` (monkeypatch `QMessageBox.question` → No first, in case of unsaved-changes prompt).
   - Assert atoms/bonds/next_atom_id match the recorded values (compare by symbol/charge/bond order, not object identity; positions with `pytest.approx`).
2. `test_pmeprj_save_via_dialog_appends_extension(window, qtbot, tmp_path, monkeypatch)`
   - Draw ethane; monkeypatch `io_logic.QFileDialog.getSaveFileName` → `(str(tmp_path/"noext"), None)`; call `window.io_manager.save_project_as()`; assert `(tmp_path/"noext.pmeprj").exists()` and `window.state_manager.has_unsaved_changes is False`.
3. `test_pmeraw_roundtrip(window, qtbot, tmp_path, monkeypatch)` — same as (1) but through `save_raw_data`/`load_raw_data` with the save dialog patched.
4. `test_load_prompts_on_unsaved_changes(window, qtbot, tmp_path, monkeypatch)`
   - Draw ethane so `has_unsaved_changes` True; patch `app_state.QMessageBox.question` → Cancel; call `open_project_file(str(existing_pmeprj))`; assert document unchanged (atoms still present).

Pitfalls: `create_json_data` needs no dialog; `open_project_file` may fire deferred `QTimer.singleShot` plotter calls — they are already None-guarded, but call `qtbot.wait(150)` after load before asserting.

## Step 2 — `tests/e2e/test_export_files.py`

1. `test_save_as_mol_2d(window, tmp_path, monkeypatch)` — draw ethane, patch save dialog → `tmp_path/"out"` (NO extension), call `save_as_mol()`; assert `out.mol` exists; parse with `Chem.MolFromMolFile(str(p))`; assert 2 heavy atoms; assert `"MoleditPy"` in file text (header line).
2. `test_save_3d_as_mol_appends_extension(window, qtbot, tmp_path, monkeypatch)` — draw + convert to 3D (copy conversion pattern), patch dialog → bare name; `save_3d_as_mol()`; assert `.mol` created and RDKit-parseable with a 3D conformer (`mol.GetConformer().Is3D()` may be False for RDKit-written blocks — assert conformer exists and z-coords not all zero instead; if all-zero on some RDKit versions, relax to conformer-exists).
3. `test_save_as_xyz(window, qtbot, tmp_path, monkeypatch)` — after conversion, patch dialog → bare name; assert `.xyz` created; first line == atom count (int), second line contains `chrg =` and `mult =`; each remaining line has 4 fields.

## Step 3 — `tests/e2e/test_import_workflows.py`

1. `test_smiles_import_creates_scene_atoms(window)` — `load_from_smiles("CCO")`; assert `len(window.data.atoms) == 3`; bond count 2; `has_unsaved_changes` True.
2. `test_invalid_smiles_no_crash(window)` — `load_from_smiles("not_a_smiles")`; assert atom count unchanged (0) and app alive (`window.isVisible() in (True, False)` — i.e., no exception raised).
3. `test_mol_file_import(window, tmp_path)` — write a benzene MOL via RDKit (`Chem.MolToMolBlock(Chem.MolFromSmiles("c1ccccc1"))` after `Compute2DCoords`), save to tmp file, `window.io_manager.load_mol_file(str(p))`; assert 6 atoms, 6 bonds; kekulized orders are ints in `{1, 2}`.
4. `test_xyz_load_enters_viewer_mode(window, qtbot, tmp_path, monkeypatch)` — write a valid water XYZ (3 atoms, header count + comment line); patch `QMessageBox.question` → No and set `window.init_manager.settings["skip_chemistry_checks"] = True` (avoids the charge dialog entirely); call `load_xyz_for_3d_viewing(str(p))`; `qtbot.wait(150)`; assert `window.current_mol` is not None with 3 atoms and `window.ui_manager.is_2d_editable is False`.

## Step 4 — `tests/e2e/test_undo_redo_workflow.py`

1. `test_draw_undo_redo_cycle(window)` — reset history; draw atom A (push), atom B (push); `undo()` → 1 atom; `undo()` → 0 atoms (stack floor: undo() only acts when len>1, so first push is the empty baseline — draw AFTER `reset_undo_stack()` and push after each atom); `redo()` twice → 2 atoms again. Assert exact atom counts at each step.
2. `test_undo_restores_bond_data(window)` — draw ethane (bond order 1), push; change bond via `window.data.add_bond(a, b, order=2)` + push; `undo()`; assert bond order back to 1.
3. `test_undo_stack_capped_e2e(window)` — from `moleditpy.utils.constants import UNDO_STACK_MAX_DEPTH` (use `_PKG` indirection); loop `UNDO_STACK_MAX_DEPTH + 5` times: move one atom (`window.data.set_atom_pos(aid, (i, 0))`) + push; assert `len(undo_stack) == UNDO_STACK_MAX_DEPTH`.

## Step 5 — `tests/e2e/test_plugin_workflow.py`

Uses the REAL `PluginManager` (not the dummy in the window fixture) against `tmp_path` — no `window` needed for 1–3.
1. `test_install_and_discover_py_plugin(tmp_path)` — write a minimal plugin (`PLUGIN_NAME`, `initialize(context)` registering a menu action), `pm = PluginManager(main_window=None)`; `pm.plugin_dir = str(tmp_path/"plugins")`; `pm.install_plugin(src)`; `pm.discover_plugins()`; assert plugin listed with status `"Loaded"` and 1 menu action.
2. `test_plugin_save_load_state_roundtrip(window, monkeypatch)` — attach a real PluginManager to `window.plugin_manager` with a save handler returning `{"k": 1}` and a load handler capturing input; `create_json_data()` → assert `d["plugins"]["Name"] == {"k": 1}`; `load_from_json_data(d)` → assert load handler received `{"k": 1}`.
3. `test_zip_single_file_install(tmp_path)` — zip containing only `solo.py`; install; assert extracted into wrapper folder `<zipname>/solo.py` (regression for the review fix).

## Step 6 — Wire-up + final gate

1. No runner changes needed — `run_all_tests.py` picks up `tests/e2e/` wholesale. Confirm:
   `MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python tests/run_all_tests.py --e2e --no-cov --no-report`
2. Full gate before final commit:
   - `python -m pytest tests/e2e -q` (all new + old pass)
   - `python -m pytest tests/unit -q` (unchanged, still green)
   - `python -m mypy moleditpy/src/moleditpy/ --config-file mypy.ini` → 0 (tests aren't type-checked, but confirms no accidental src edits)
   - `python -m ruff check tests/e2e/`
3. Update `tests/e2e/README.md`: add one table row per new file (one line each).

## Failure protocol (for the executing model)

- A test hangs > 60 s → a real dialog opened; kill, find the unpatched `QFileDialog`/`QMessageBox`/`QInputDialog` in the code path, patch it, retry.
- `AttributeError` on `window.X` → check the attribute name in `moleditpy/src/moleditpy/ui/main_window.py` proxies before assuming a fixture problem.
- Conversion never finishes → the `window` fixture disables the worker thread (`init_worker_thread` is patched out); use the same conversion trigger as `test_ethane_gui.py::test_ethane_2d_to_3d_via_button`, not a bespoke one.
- Any need to touch `moleditpy/src/**` → STOP, report the blocker, do not work around it in production code.

**Estimated size:** 5 files, ~17 tests, ~600 lines of test code. Priority order = step order (1 and 2 carry the most user-facing risk coverage).
