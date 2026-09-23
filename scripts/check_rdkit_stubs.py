#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Report whether the installed rdkit ships type stubs mypy can use.

mypy.ini skips rdkit's bundled ``rdkit-stubs`` (``follow_imports = skip`` in
the ``[mypy-rdkit.*]`` section) because some of the generated ``.pyi`` files
are not valid Python, and a single unparseable stub stops mypy before it checks
any of our code. Run this after upgrading rdkit::

    python scripts/check_rdkit_stubs.py

Exit status 0 means every stub parses: rdkit type checking can be turned back
on by deleting that one ``follow_imports = skip`` line, then running mypy on
all three platforms and fixing what it reports. Exit status 1 lists the broken
stubs; leave the skip in place.
"""

from __future__ import annotations

import ast
import sys
from importlib.util import find_spec
from pathlib import Path


def main() -> int:
    """Parse every rdkit stub and report the ones that fail."""
    spec = find_spec("rdkit")
    if spec is None or spec.origin is None:
        print("rdkit is not installed.")
        return 1
    import rdkit

    stubs = Path(spec.origin).parent.parent / "rdkit-stubs"
    print(f"rdkit {rdkit.__version__}")
    if not stubs.is_dir():
        print("No rdkit-stubs shipped: rdkit is untyped, the skip is a no-op.")
        return 1

    files = sorted(stubs.rglob("*.pyi"))
    broken = []
    for path in files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            broken.append(f"  {path.relative_to(stubs.parent)}:{exc.lineno}: {exc.msg}")

    if broken:
        print(f"{len(broken)} of {len(files)} stubs do not parse:")
        print("\n".join(broken))
        print("Keep 'follow_imports = skip' in mypy.ini.")
        return 1
    print(f"All {len(files)} stubs parse. rdkit type checking can be re-enabled:")
    print("delete 'follow_imports = skip' under [mypy-rdkit.*] in mypy.ini.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
