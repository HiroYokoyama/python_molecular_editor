#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate every "file icon" asset from scratch, in one deterministic pass.

The file icon (what Explorer/Finder show for .pmeprj / .pmeraw, distinct from
the application icon) used to be built in two disconnected steps: a
matplotlib script drew the paper/bar background, and the application icon's
molecule artwork was pasted on top by hand in GIMP (``file_icon.xcf``), then
exported once and copied around by hand into every size, .ico and .icns. Two
manual steps drift -- the bar's blue and the app's brand blue had quietly
diverged. This script does both steps itself, so re-running it after either
the background design or the application icon changes reproduces every
derived file identically, on any platform.

Run from the repo root::

    python scripts/generate_file_icon.py

Requires matplotlib and Pillow (dev-only tools; not runtime dependencies of
the application itself).
"""

from __future__ import annotations

import io
import os
import struct

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICON_DIR = os.path.join(BASE_DIR, "icon_images", "file-icon")

#: The application icon (molecule + "MoleditPy" script text), pasted onto
#: the generated background. This is the single source both this file icon
#: and the app icon itself are derived from.
APP_ICON_PNG = os.path.join(BASE_DIR, "icon_images", "icon.png")

MASTER_PNG = os.path.join(ICON_DIR, "file_icon.png")
SMALL_DIR = os.path.join(ICON_DIR, "small")
ICO_PATH = os.path.join(ICON_DIR, "file_icon.ico")
ICNS_PATH = os.path.join(ICON_DIR, "file_icon.icns")

#: Where the built .ico is actually consumed by the app. moleditpy-linux's
#: copy is synced separately by scripts/sync_linux_version.py, which mirrors
#: any binary file under assets/ -- no separate step needed here for that.
APP_ASSET_ICO = os.path.join(
    BASE_DIR, "moleditpy", "src", "moleditpy", "assets", "file_icon.ico"
)

#: Also bundled by the installer, which has no automated pipeline of its own
#: pointed at this repo -- keep this list in sync by hand if that changes.
INSTALLER_DATA_DIR = os.path.join(
    os.path.dirname(BASE_DIR),
    "python_molecular_editor_installer",
    "moleditpy-installer",
    "moleditpy_installer",
    "data",
)

#: MoleditPy Blue -- see docs/BRAND_COLOR.md. The single place this script's
#: output should ever get its blue from.
BRAND_BLUE = "#3577F7"

#: Sizes exported as standalone PNGs (used by the Linux .desktop icon theme).
SMALL_SIZES = (16, 22, 24, 32, 48, 64, 128)

#: Sizes Pillow's ICNS writer always bakes in (32, 64, 128, 256, 512, 1024 --
#: 256 and 512 each cover two of its internal type codes). This is NOT a
#: `sizes=` argument to Image.save: Pillow's ICNS plugin ignores that kwarg
#: entirely and hardcodes this exact size set, falling back to a plain
#: (non-LANCZOS) im.resize() for any size not supplied via append_images.
#: Passing every size here explicitly is what makes each entry a proper
#: high-quality downscale of the master instead of that fallback.
ICNS_SIZES = (32, 64, 128, 256, 512, 1024)

#: The background is drawn in a 100x100-ish data-unit space (matplotlib
#: "axes" units); everything below is in that space, not pixels.
DOC_W = 70
DOC_H = 90
DOC_X = (100 - DOC_W) / 2
DOC_Y = (100 - DOC_H) / 2
FOLD_SIZE = 20
BAR_H = 16
BAR_Y = DOC_Y + 12

#: Margin (in the same data units) kept clear around the molecule artwork:
#: below the folded corner at the top, above the bar at the bottom, and in
#: from both sides.
CONTENT_MARGIN = 6

#: Extra shrink applied after fitting the molecule to its available box, so
#: it doesn't run edge-to-edge against the margin on its long axis.
MOLECULE_SCALE = 0.8

#: How the leftover vertical space (after fitting+shrinking) is split between
#: the top and bottom of the molecule's box: 0.5 centers it; higher pushes it
#: down (more slack kept above it, less below), which reads better than dead
#: center since the fold notch already crowds the top-right corner.
VERTICAL_BIAS = 0.95

#: "MoleditPy File" text size, in points (matplotlib fontsize).
BAR_FONT_SIZE = 52


#: The data-coordinate window the background is drawn in. Kept as a fixed,
#: known range (rather than relying on matplotlib's autoscale/tight-bbox
#: cropping) so a data coordinate maps to a pixel coordinate by one constant
#: scale factor -- which _paste_molecule depends on to place the artwork.
AXES_RANGE = (-5, 105)


def _draw_background(size_px: int) -> plt.Figure:
    """The paper/bar background, ported from the former create-icon-bkg.py."""
    dpi = 100
    fig_size_inch = size_px / dpi
    fig = plt.figure(figsize=(fig_size_inch, fig_size_inch), dpi=dpi)
    # Axes fill the entire figure with zero margin, and bbox_inches="tight" is
    # never used when saving: both would crop the canvas to whatever content
    # happens to be drawn, decoupling pixel position from data coordinates.
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(*AXES_RANGE)
    ax.set_ylim(*AXES_RANGE)
    ax.set_aspect("equal")
    ax.axis("off")

    p1 = (DOC_X, DOC_Y)
    p2 = (DOC_X + DOC_W, DOC_Y)
    p3 = (DOC_X + DOC_W, DOC_Y + DOC_H - FOLD_SIZE)
    p4 = (DOC_X + DOC_W - FOLD_SIZE, DOC_Y + DOC_H)
    p5 = (DOC_X, DOC_Y + DOC_H)
    paper_coords = [p1, p2, p3, p4, p5]

    shadow_coords = [(p[0] + 3, p[1] - 3) for p in paper_coords]
    shadow_poly = patches.Polygon(
        shadow_coords, closed=True, facecolor="black", alpha=0.2, zorder=0,
        joinstyle="round",
    )
    shadow_poly.set_path_effects(
        [path_effects.withStroke(linewidth=15, foreground="black", alpha=0.1)]
    )
    ax.add_patch(shadow_poly)

    paper_poly = patches.Polygon(
        paper_coords, closed=True, facecolor="#F9F9F9", edgecolor="#DDDDDD",
        linewidth=1, zorder=1, joinstyle="round",
    )
    ax.add_patch(paper_poly)

    hex_size = 8
    h_step = hex_size * np.sqrt(3)
    v_step = hex_size * 1.5
    # Padded by a full step (not just hex_size) past each edge: the pattern
    # is clipped to paper_poly regardless, so overshooting is free, but
    # undershooting left a bare strip along the top edge where the last row's
    # center was far enough below y=DOC_Y+DOC_H that the hexagon didn't reach it.
    for row in range(
        int((DOC_Y - v_step) // v_step), int((DOC_Y + DOC_H + v_step) // v_step) + 1
    ):
        for col in range(
            int((DOC_X - h_step) // h_step), int((DOC_X + DOC_W + h_step) // h_step) + 1
        ):
            x_pos = col * h_step
            y_pos = row * v_step
            if row % 2 == 1:
                x_pos += h_step / 2
            angles = np.linspace(0, 2 * np.pi, 7)
            poly = patches.Polygon(
                np.column_stack((x_pos + hex_size * np.cos(angles), y_pos + hex_size * np.sin(angles))),
                closed=True, edgecolor="#39CCCC", facecolor="none", linewidth=1,
                alpha=0.15, zorder=2,
            )
            poly.set_clip_path(paper_poly)
            ax.add_patch(poly)

    fold_coords = [p3, (DOC_X + DOC_W - FOLD_SIZE, DOC_Y + DOC_H - FOLD_SIZE), p4]
    fold_shadow = patches.Polygon(
        fold_coords, closed=True, fc="black", alpha=0.1, zorder=2.5, joinstyle="round",
    )
    fold_shadow.set_path_effects(
        [path_effects.withStroke(linewidth=5, foreground="black", alpha=0.1)]
    )
    ax.add_patch(fold_shadow)
    fold_poly = patches.Polygon(
        fold_coords, closed=True, facecolor="#EEEEEE", edgecolor="#CCCCCC",
        linewidth=1, zorder=3, joinstyle="round",
    )
    ax.add_patch(fold_poly)

    bar_rect = patches.Rectangle(
        (DOC_X, BAR_Y), DOC_W, BAR_H, facecolor=BRAND_BLUE, alpha=1.0, zorder=2
    )
    bar_rect.set_clip_path(paper_poly)
    ax.add_patch(bar_rect)

    ax.text(
        DOC_X + DOC_W / 2, BAR_Y + BAR_H / 2, "MoleditPy File",
        ha="center", va="center", fontsize=BAR_FONT_SIZE, color="white",
        fontweight="bold", fontname="DejaVu Sans", zorder=3,
    )
    return fig


def _figure_to_image(fig: plt.Figure) -> Image.Image:
    buf = io.BytesIO()
    # No bbox_inches="tight": that crops to whichever pixels ended up
    # non-transparent, which is exactly what must NOT happen here -- the
    # output size has to stay the figure's declared size so a data
    # coordinate maps to a pixel coordinate by one constant scale factor.
    fig.savefig(buf, transparent=True, dpi=fig.dpi)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGBA")


def _paste_molecule(background: Image.Image) -> Image.Image:
    """Scale the application icon to fit between the top of the paper and
    the bar, centered, and paste it onto the background."""
    molecule = Image.open(APP_ICON_PNG).convert("RGBA")
    mol_arr = np.array(molecule)
    ys, xs = np.where(mol_arr[:, :, 3] > 10)
    mol_box = molecule.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))

    W, H = background.size
    px_per_unit = W / (AXES_RANGE[1] - AXES_RANGE[0])

    def x_to_px(x_data: float) -> float:
        return (x_data - AXES_RANGE[0]) * px_per_unit

    def y_to_px(y_data: float) -> float:
        # matplotlib's y axis increases upward; image rows increase downward.
        return (AXES_RANGE[1] - y_data) * px_per_unit

    avail_x0 = x_to_px(DOC_X + CONTENT_MARGIN)
    avail_x1 = x_to_px(DOC_X + DOC_W - CONTENT_MARGIN)
    avail_y0 = y_to_px(DOC_Y + DOC_H - CONTENT_MARGIN)  # near the top edge
    avail_y1 = y_to_px(BAR_Y + BAR_H + CONTENT_MARGIN)  # just above the bar
    avail_w = avail_x1 - avail_x0
    avail_h = avail_y1 - avail_y0

    scale = min(avail_w / mol_box.width, avail_h / mol_box.height) * MOLECULE_SCALE
    new_size = (max(1, round(mol_box.width * scale)), max(1, round(mol_box.height * scale)))
    mol_resized = mol_box.resize(new_size, Image.LANCZOS)

    paste_x = round(avail_x0 + (avail_w - new_size[0]) / 2)
    paste_y = round(avail_y0 + (avail_h - new_size[1]) * VERTICAL_BIAS)

    out = background.copy()
    out.alpha_composite(mol_resized, (paste_x, paste_y))
    return out


def _write_single_frame_ico(master: Image.Image, path: str) -> None:
    """A one-frame ICO whose directory claims 256x256 but whose PNG payload
    is the full-resolution master (Windows scales the declared size)."""
    buf = io.BytesIO()
    master.save(buf, format="PNG")
    png_bytes = buf.getvalue()
    header = struct.pack("<HHH", 0, 1, 1)
    entry = struct.pack("<BBBBHHII", 0, 0, 0, 0, 1, 32, len(png_bytes), 6 + 16)
    with open(path, "wb") as f:
        f.write(header)
        f.write(entry)
        f.write(png_bytes)


def main() -> int:
    background = _figure_to_image(_draw_background(size_px=1024))
    master = _paste_molecule(background)
    master.save(MASTER_PNG)

    os.makedirs(SMALL_DIR, exist_ok=True)
    small_paths = []
    for size in SMALL_SIZES:
        path = os.path.join(SMALL_DIR, f"file_icon_{size}.png")
        master.resize((size, size), Image.LANCZOS).save(path)
        small_paths.append(path)

    _write_single_frame_ico(master, ICO_PATH)
    icns_sizes = [s for s in ICNS_SIZES if s <= master.size[0]]
    master.save(
        ICNS_PATH,
        append_images=[master.resize((s, s), Image.LANCZOS) for s in icns_sizes],
    )

    with open(ICO_PATH, "rb") as src, open(APP_ASSET_ICO, "wb") as dst:
        dst.write(src.read())

    installer_updated = []
    if os.path.isdir(INSTALLER_DATA_DIR):
        with open(ICO_PATH, "rb") as src:
            ico_bytes = src.read()
        with open(os.path.join(INSTALLER_DATA_DIR, "file_icon.ico"), "wb") as dst:
            dst.write(ico_bytes)
        installer_updated.append("file_icon.ico")

        master.save(os.path.join(INSTALLER_DATA_DIR, "file_icon.png"))
        installer_updated.append("file_icon.png")

        with open(ICNS_PATH, "rb") as src:
            icns_bytes = src.read()
        with open(os.path.join(INSTALLER_DATA_DIR, "file_icon.icns"), "wb") as dst:
            dst.write(icns_bytes)
        installer_updated.append("file_icon.icns")

        for size in SMALL_SIZES:
            name = f"file_icon_{size}.png"
            master.resize((size, size), Image.LANCZOS).save(
                os.path.join(INSTALLER_DATA_DIR, name)
            )
            installer_updated.append(name)

    print(f"Wrote {MASTER_PNG}")
    print(f"Wrote {len(small_paths)} PNG sizes under {SMALL_DIR}")
    print(f"Wrote {ICO_PATH} and {ICNS_PATH}")
    print(f"Copied ICO to {APP_ASSET_ICO}")
    if installer_updated:
        print(f"Updated {len(installer_updated)} files in {INSTALLER_DATA_DIR}")
    else:
        print(f"Installer repo not found beside this one ({INSTALLER_DATA_DIR}); skipped")
    print(
        "Note: moleditpy-linux/.../assets/file_icon.ico is synced separately "
        "by scripts/sync_linux_version.py."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
