#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Regenerate every derived "file icon" asset from the master PNG.

The file icon (the one Explorer/Finder show for .pmeprj / .pmeraw, distinct
from the application icon) exists as a dozen exported copies: a handful of
raster sizes, a Windows .ico, and a macOS .icns. Exporting each by hand from
the GIMP source (``file_icon.xcf``) is how they drift, so this script treats
``icon_images/file-icon/file_icon.png`` as the single source of truth and
regenerates every other file from it.

Run after editing the master PNG (e.g. after exporting a new version from
``file_icon.xcf`` in GIMP)::

    python scripts/generate_file_icon.py

The .xcf itself is not touched here -- it is the editable design source and
has no automated pipeline of its own.
"""

from __future__ import annotations

import os
import struct

from PIL import Image

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICON_DIR = os.path.join(BASE_DIR, "icon_images", "file-icon")
MASTER_PNG = os.path.join(ICON_DIR, "file_icon.png")
SMALL_DIR = os.path.join(ICON_DIR, "small")
ICO_PATH = os.path.join(ICON_DIR, "file_icon.ico")
ICNS_PATH = os.path.join(ICON_DIR, "file_icon.icns")

#: Where the built .ico is actually consumed by the app (and copied into
#: moleditpy-linux by scripts/sync_linux_version.py, which mirrors any
#: binary file in assets/ -- no separate step needed here for that copy).
APP_ASSET_ICO = os.path.join(
    BASE_DIR, "moleditpy", "src", "moleditpy", "assets", "file_icon.ico"
)

#: Sizes exported as standalone PNGs (used by the Linux .desktop icon theme).
SMALL_SIZES = (16, 22, 24, 32, 48, 64, 128)

#: Sizes baked into the macOS .icns.
ICNS_SIZES = ((16, 16), (32, 32), (128, 128), (256, 256), (512, 512))


def _write_single_frame_ico(master: Image.Image, path: str) -> None:
    """Write a one-frame ICO whose directory claims 256x256 but whose PNG
    payload is the full-resolution master, matching how this asset has always
    shipped (Windows scales the declared size; the extra resolution is used
    directly by anything that reads the PNG frame instead of trusting the
    directory, e.g. some installers)."""
    import io

    buf = io.BytesIO()
    master.save(buf, format="PNG")
    png_bytes = buf.getvalue()

    header = struct.pack("<HHH", 0, 1, 1)
    # width/height 0 means "256" in the ICO directory format.
    entry = struct.pack(
        "<BBBBHHII", 0, 0, 0, 0, 1, 32, len(png_bytes), 6 + 16
    )
    with open(path, "wb") as f:
        f.write(header)
        f.write(entry)
        f.write(png_bytes)


def main() -> int:
    master = Image.open(MASTER_PNG).convert("RGBA")

    os.makedirs(SMALL_DIR, exist_ok=True)
    for size in SMALL_SIZES:
        resized = master.resize((size, size), Image.LANCZOS)
        resized.save(os.path.join(SMALL_DIR, f"file_icon_{size}.png"))

    _write_single_frame_ico(master, ICO_PATH)

    master.save(
        ICNS_PATH,
        sizes=[s for s in ICNS_SIZES if s[0] <= master.size[0]],
    )

    with open(ICO_PATH, "rb") as src, open(APP_ASSET_ICO, "wb") as dst:
        dst.write(src.read())

    print(f"Regenerated {len(SMALL_SIZES)} PNG sizes, {ICO_PATH}, {ICNS_PATH}")
    print(f"Copied ICO to {APP_ASSET_ICO}")
    print(
        "Note: moleditpy-linux/.../assets/file_icon.ico is synced separately "
        "by scripts/sync_linux_version.py, and python_molecular_editor_installer's "
        "bundled copy must be updated by hand in that repo."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
