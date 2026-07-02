# MoleditPy — Final Code Quality Report (fable5)

- **Date:** 2026-07-02 (final revision — supersedes earlier versions of this file)
- **Author:** Claude Code (Fable 5)
- **Scope:** `python_molecular_editor/moleditpy/src/moleditpy/` (~30k lines, 66 modules), branch **`dev-4.2.0`** (39 commits ahead of `main`, pushed to origin)
- **Companions:** `code_review_fable5.md` (original findings), `e2e_test_plan_fable5.md` (executed in full)

---

## 1. Final metrics (verified on the pushed tree)

| Check | Result | Session start |
|---|---|---|
| **mypy** (py3.9 target) | **0 errors** in 66 files | could not run → 1111 errors |
| **ruff check + format** | **fully clean** | format never applied; 2 F401 |
| **pylint** (full) | **9.33 / 10** (target > 9.0) | ~9.39 (−0.06 = line-length from ignore/logging comments) |
| **pylint** (errors-only) | **10.00 / 10** | — |
| **Tests** | **1354 passed** (1255 unit + 64 integration + 32 E2E + new additions; 1 platform skip) | 1244, no E2E workflows |
| **Test side effects** | **Zero** — full suite leaves `~/.moleditpy` (430 files) and repo tree byte-identical (audited) | E2E tests wrote to real user settings |
| **App startup** | Verified importable after every commit | — |

## 2. Everything fixed on `dev-4.2.0`

### Correctness bugs (each with regression tests where testable)
- **Core:** `add_bond` duplicate entries on reversed stereo/non-stereo keys; manual MOL-block fallback crash on aromatic order + atom-count desync; recursive union-find.
- **I/O & state:** `save_3d_as_mol` missing `.mol` extension; Ctrl+S not refreshing the saved-state snapshot; direct project loads bypassing the unsaved-changes check; stale 2D→3D atom map (`set_atom_id_to_rdkit_idx_map` wrote to an attribute nothing read); undo history uncapped (now `UNDO_STACK_MAX_DEPTH` = 100).
- **Dialogs:** charge dialog silently accepting invalid input → inline validation in the same dialog (always-visible error label; the initial `setVisible` toggle didn't render reliably and was replaced with `setText` + `adjustSize`).
- **3D/compute:** multiple None-plotter crash paths guarded (overlay removal, deferred QTimer lambdas, `reset_camera`); Open Babel null-check ordering; dead code removed (unreachable constants-import fallback, nonexistent `ColorSettingsDialog.refresh_ui()` call, always-true guard, no-op env line).
- **Plugins:** ZIP single-file misdetection; dot-directories scanned for plugins; **`register_optimization_method` was registered/validated but never invoked** — now wired: appears in the Optimize 3D right-click menu as `<name> (Plugin)`, runs synchronously, success redraws + pushes undo, failures logged and status-reported. Manual updated accordingly.
- **Platform:** macOS dark-mode detection implemented; `VDW_RADII` → `VDW_DISPLAY_RADII` (deprecated alias kept); versioned AppUserModelID.

### Error-handling hardening (owner-directed)
- New `utils/suppress_log.py`: `contextlib.suppress` semantics + DEBUG-level traceback logging. **All 61 suppress sites converted; all 90 silent `except: pass` handlers now log.** Nothing in the app swallows an exception invisibly.
- Plugin touchpoints stay **broadly defensive** (plugins have full app access): menu/toolbar/export callbacks, drop handlers, custom 3D styles, save/load/reset handlers, file openers (widened to `except Exception`), and the new optimization dispatch — all catch everything, log full tracebacks, never crash the app.

### Typing (mypy 1111 → 0)
~440 genuine annotation fixes + 180 per-line code-specific `# type: ignore[...]` for Qt/dynamic idioms. No blanket ignores, no config-level disabling. Known caveat: UI widget attrs are pragmatically `Any`.

### Test infrastructure
- **18 new E2E workflow tests** (project round-trip, export, import, undo/redo, plugins) on the real MainWindow, plus fixture repairs.
- **E2E home-dir sandbox**: fixtures redirect `HOME`/`USERPROFILE` to tmp — this fixed a real incident where a test persisted `skip_chemistry_checks=true` into the user's actual settings.json (breaking XYZ charge prompting in the real app). Post-fix audit: zero file modifications anywhere.
- One order-dependent unit test fixed (mock dialog patch); deferred-QTimer teardown pattern documented in tests.

### Documentation
- CLAUDE.md version claim, workspace API-name staleness, ARCHITECTURE.md error-policy claim, plugin manual (`register_optimization_method` real behaviour), E2E README table, and **README security note (EN/JA) for `.pmeraw` pickle loading**.

## 3. Verdict

**A− for its category.** Correctness and observability are now strong: two full review passes with all findings fixed, four independent quality gates green, no silent error paths, workflow-level E2E coverage that has already caught real regressions twice, and side-effect-free tests (audited, not assumed).

Compared to general OSS: comfortably top-tier vs. the median Python project; exceptional vs. scientific GUI tools; the remaining gaps vs. flagship infrastructure are refactor-scale, not defect-scale.

## 4. Remaining items (known, deliberate, non-blocking)

| Item | Class |
|---|---|
| `view_3d_logic.py` ~2000 lines; ~300 pylint complexity warnings; 92 duplicate-code blocks; 114 f-string logging calls | Refactor-scale style debt — the largest remaining quality item |
| UI attrs typed `Any` instead of precise widget types | Typing depth; tighten opportunistically |
| `.pmeraw` pickle loading | Kept by owner decision; **now documented in README (EN/JA)** |
| Undo comparison ignores `constraints_3d`; SMILES/InChI importer duplication | Accepted as-is by owner |
| E2E tests use mocked VTK (real rendering untested headlessly); GUI suite runs in CI only | Inherent to headless testing |
| Bus factor of one | Structural |

## 5. Process record

39 commits, every one gated on package-import smoke + mypy + tests. Four self-inflicted regressions during the campaign — missing `Any` imports (startup crash), gesture-handling `isinstance` change, obabel subprocess script corruption, e2e settings pollution — **all four caught by the gates or the new E2E suite within the same session and fixed**, which is the strongest evidence the safety net now works as designed.
