#!/usr/bin/env python3
"""
fix_test_imports.py -- Auto-fix stale imports in test files using AST discovery.

Renames import targets and all usages in test files to match current dev source:
  * Old class/function name -> new name (rename throughout the file)
  * Old module -> new module (when the function/class moved)

Discovered renames (from scan_test_disconnections.py):
  MainWindowMainInit -> MainInitManager  (in main_window_init.py)
  _set_mol_prop_safe: import from io_logic instead of molecular_parsers

Usage:
    python scripts/fix_test_imports.py           # dry-run (show changes only)
    python scripts/fix_test_imports.py --apply   # apply changes to test files
"""

import re
import argparse
from pathlib import Path

ROOT  = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests" / "unit"

# ---------------------------------------------------------------------------
# Rename rules: (old_module, old_name, new_module, new_name, reason)
# ---------------------------------------------------------------------------
RENAME_RULES = [
    (
        "moleditpy.ui.main_window_init",
        "MainWindowMainInit",
        "moleditpy.ui.main_window_init",
        "MainInitManager",
        "Renamed: MainWindowMainInit -> MainInitManager",
    ),
    (
        "moleditpy.ui.molecular_parsers",
        "_set_mol_prop_safe",
        "moleditpy.ui.io_logic",
        "_set_mol_prop_safe",
        "Moved: _set_mol_prop_safe now lives in io_logic",
    ),
]


def fix_file(fpath: Path, rules: list, apply: bool) -> list:
    """Return list of (description, old_text, new_text) changes."""
    text = fpath.read_text(encoding="utf-8", errors="replace")
    new_text = text
    changes = []

    for old_mod, old_name, new_mod, new_name, reason in rules:
        # 1) Fix the import line
        #    "from old_mod import old_name" -> "from new_mod import new_name"
        import_pattern = re.compile(
            rf"(from\s+{re.escape(old_mod)}\s+import\s+)([\w, ]+)"
        )
        def fix_import(m, old_name=old_name, old_mod=old_mod,
                       new_mod=new_mod, new_name=new_name):
            prefix = m.group(1)
            names_str = m.group(2)
            names = [n.strip() for n in names_str.split(",")]
            if old_name not in names:
                return m.group(0)
            # Replace just this name in the list
            new_names = [new_name if n == old_name else n for n in names]
            new_prefix = f"from {new_mod} import " if new_mod != old_mod else prefix
            return new_prefix + ", ".join(new_names)

        new_text2 = import_pattern.sub(fix_import, new_text)
        if new_text2 != new_text:
            changes.append((f"import: {old_mod}.{old_name} -> {new_mod}.{new_name}",
                             reason))
            new_text = new_text2

        # 2) Replace all usages of old_name with new_name in the rest of the file
        #    (only if names differ)
        if old_name != new_name:
            usage_pattern = re.compile(rf"\b{re.escape(old_name)}\b")
            new_text2 = usage_pattern.sub(new_name, new_text)
            if new_text2 != new_text:
                changes.append((f"usage: {old_name} -> {new_name}", reason))
                new_text = new_text2

    if changes and apply:
        fpath.write_text(new_text, encoding="utf-8")

    return changes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true",
                    help="Apply fixes (default: dry-run)")
    args = ap.parse_args()

    test_files = sorted(TESTS.glob("test_*.py"))
    total = 0

    for fpath in test_files:
        changes = fix_file(fpath, RENAME_RULES, args.apply)
        if not changes:
            continue
        rel = fpath.relative_to(ROOT)
        status = "[FIXED]" if args.apply else "[DRY-RUN]"
        print(f"\n{rel}:")
        for desc, reason in changes:
            total += 1
            print(f"  {status} {desc}")
            print(f"           Reason: {reason}")

    print(f"\nTotal changes: {total}", end="")
    if not args.apply and total:
        print(" (dry-run -- run with --apply to apply)")
    else:
        print()


if __name__ == "__main__":
    main()
