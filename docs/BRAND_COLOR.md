# Brand Color

## MoleditPy Blue — `#3577F7`

This is the project's brand blue, taken from the application icon
(`moleditpy/src/moleditpy/assets/icon.png` / `icon.ico`). It is used wherever
the app needs a single accent color that should read as "MoleditPy", rather
than an arbitrary UI blue.

| Value | Where |
|---|---|
| `#3577F7` | RGB `(53, 119, 247)` |

### Where it's used

- **Nitrogen's default CPK color** — `moleditpy/src/moleditpy/utils/constants.py`,
  `CPK_COLORS["N"]`. Nitrogen's standard CPK blue happens to sit close to the
  brand color, so the default was tightened to the exact brand hex instead of
  an independently-chosen blue.
- **The file icon** — `icon_images/file-icon/file_icon.png` (and every size,
  `.ico` and `.icns` derived from it by `scripts/generate_file_icon.py`). The
  "MoleditPy File" bar and the molecule's nitrogen atoms used an
  independently-picked blue (`#1981dc`) and now use this exact hex. The
  *application* icon (`icon.png` / `icon.ico`, no bar) is unchanged — only
  the file-association icon carries the bar.
- **moleditpy_job_manager** — the blue rack unit in its favicon
  (`FAVICON_SVG` in `job_manager/web_monitor.py`) uses this value, so the
  plugin's icon reads as part of the MoleditPy family rather than a
  generic blue. The plugin's UI accent color (`CY_ACCENT` in `theme.py`) is
  a separate, deliberately-chosen blue and is not tied to this value.

### Adding it elsewhere

Do not eyeball a "close enough" blue for a new MoleditPy-branded accent —
use `#3577F7` exactly, so every place that claims to be "the brand color"
actually matches. Derived shades (e.g. a darker variant for hover/focus
states) should be computed from this value, not from a different blue.
