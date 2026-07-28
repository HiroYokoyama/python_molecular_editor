import os
import shutil
import re
import sys
import filecmp

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_MAIN = os.path.join(BASE_DIR, "moleditpy", "src", "moleditpy")
SRC_LINUX = os.path.join(BASE_DIR, "moleditpy-linux", "src", "moleditpy_linux")
MAIN_TOML = os.path.join(BASE_DIR, "moleditpy", "pyproject.toml")
LINUX_TOML = os.path.join(BASE_DIR, "moleditpy-linux", "pyproject.toml")

# `moleditpy` only where it names the package: `moleditpy.x`, `import moleditpy`,
# `from moleditpy`. Leaves `~/.moleditpy` and other bare strings alone.
PACKAGE_REF = re.compile(
    r"\bmoleditpy(?=\.)|(?<=\bimport )moleditpy\b|(?<=\bfrom )moleditpy\b"
)
# Renamed above but not module paths: the log file and the Windows app id.
LITERAL_RESTORES = (
    ("moleditpy_linux.log", "moleditpy.log"),
    ("hyoko.moleditpy_linux.", "hyoko.moleditpy."),
)
OBABEL_LINE = re.compile(
    r"^(\s*)OBABEL_AVAILABLE = (?:True|False|importlib\..*)$", re.MULTILINE
)
DEPENDENCIES_BLOCK = re.compile(
    r"^dependencies\s*=\s*\[(.*?)\]", re.MULTILINE | re.DOTALL
)
DEP_NAME = re.compile(r"^[A-Za-z0-9._-]+")


def transform(content: str, is_root_init: bool) -> str:
    """Rewrite one main-package source file into its Linux equivalent."""
    content = PACKAGE_REF.sub("moleditpy_linux", content)
    for renamed, original in LITERAL_RESTORES:
        content = content.replace(renamed, original)
    content = content.replace('version("MoleditPy")', 'version("MoleditPy-linux")')
    if is_root_init:
        content = OBABEL_LINE.sub(r"\1OBABEL_AVAILABLE = False", content)
    return content


def read_dependencies(path: str) -> dict:
    """Map package name -> full requirement string from a pyproject."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            block = DEPENDENCIES_BLOCK.search(f.read())
    except OSError:
        return {}
    if not block:
        return {}
    deps = {}
    for raw in re.findall(r"[\"']([^\"']+)[\"']", block.group(1)):
        name = DEP_NAME.match(raw.strip())
        if name:
            deps[name.group(0).lower()] = raw.strip()
    return deps


def check_shared_pins() -> list:
    """Requirements present in both pyprojects must carry identical specifiers."""
    main_deps = read_dependencies(MAIN_TOML)
    linux_deps = read_dependencies(LINUX_TOML)
    return [
        f"  {name}: main has {main_deps[name]!r}, linux has {linux_deps[name]!r}"
        for name in sorted(set(main_deps) & set(linux_deps))
        if main_deps[name] != linux_deps[name]
    ]


def sync_version(dry_run: bool, prefix: str) -> bool:
    for path in (MAIN_TOML, LINUX_TOML):
        if not os.path.exists(path):
            print(f"Error: pyproject.toml not found at {path}")
            return False

    with open(MAIN_TOML, "r", encoding="utf-8") as f:
        main_toml = f.read()

    version_match = re.search(
        r'^version\s*=\s*["\']([^"\']+)["\']', main_toml, re.MULTILINE
    )
    if not version_match:
        print(f"Error: no version found in {MAIN_TOML}")
        return False
    new_version = version_match.group(1)

    with open(LINUX_TOML, "r", encoding="utf-8") as f:
        linux_toml = f.read()

    new_linux_toml, count = re.subn(
        r'^version\s*=\s*["\'][^"\']+["\']',
        f'version = "{new_version}"',
        linux_toml,
        flags=re.MULTILINE,
    )
    if count == 0:
        print(f"Error: no version found in {LINUX_TOML}")
        return False

    if new_linux_toml != linux_toml:
        if not dry_run:
            with open(LINUX_TOML, "w", encoding="utf-8") as f:
                f.write(new_linux_toml)
        print(f"{prefix}Synchronized pyproject.toml version to {new_version}")
    return True


def sync_linux(dry_run=False, verbose=False):
    """Return the number of files that differ, or None if the sync could not run."""
    prefix = "[DRY RUN] " if dry_run else ""
    print(f"{prefix}Synchronizing Linux version from {SRC_MAIN} to {SRC_LINUX}...")

    if not os.path.exists(SRC_MAIN):
        print(f"Error: Main source directory not found at {SRC_MAIN}")
        return None

    if not dry_run:
        os.makedirs(SRC_LINUX, exist_ok=True)

    synced_files_count = 0
    valid_dest_files = set()

    # 2. Copy all files from main to linux
    for root, dirs, files in os.walk(SRC_MAIN):
        # Skip pycache
        if "__pycache__" in dirs:
            dirs.remove("__pycache__")

        # Calculate relative path
        rel_path = os.path.relpath(root, SRC_MAIN)
        dest_dir = os.path.normpath(os.path.join(SRC_LINUX, rel_path))

        if not dry_run:
            os.makedirs(dest_dir, exist_ok=True)

        for file in files:
            if file.endswith(".pyc") or file.endswith(".pyo"):
                continue

            src_file = os.path.join(root, file)
            dest_file = os.path.join(dest_dir, file)
            valid_dest_files.add(os.path.abspath(dest_file))

            # For Python files, perform transformation
            if file.endswith(".py"):
                try:
                    with open(src_file, "r", encoding="utf-8") as f:
                        content = f.read()

                    content = transform(
                        content,
                        is_root_init=(file == "__init__.py" and rel_path == "."),
                    )

                    # Check if file needs to be updated
                    needs_update = True
                    if os.path.exists(dest_file):
                        with open(dest_file, "r", encoding="utf-8") as f:
                            existing_content = f.read()
                        if existing_content == content:
                            needs_update = False

                    if needs_update:
                        if not dry_run:
                            with open(dest_file, "w", encoding="utf-8") as f:
                                f.write(content)
                        else:
                            print(f"{prefix}Would modify {dest_file}")
                        synced_files_count += 1

                except UnicodeDecodeError:
                    print(
                        f"Warning: Could not decode {src_file} as UTF-8. Copying as binary."
                    )
                    if not os.path.exists(dest_file) or not filecmp.cmp(
                        src_file, dest_file, shallow=False
                    ):
                        if not dry_run:
                            shutil.copy2(src_file, dest_file)
                        else:
                            print(f"{prefix}Would copy {src_file} -> {dest_file}")
                        synced_files_count += 1
            else:
                # Binary files (icons, assets, etc.)
                if not os.path.exists(dest_file) or not filecmp.cmp(
                    src_file, dest_file, shallow=False
                ):
                    if not dry_run:
                        shutil.copy2(src_file, dest_file)
                    else:
                        print(f"{prefix}Would copy {src_file} -> {dest_file}")
                    synced_files_count += 1

    # Cleanup orphan files in destination
    for root, dirs, files in os.walk(SRC_LINUX, topdown=False):
        for file in files:
            dest_file = os.path.abspath(os.path.join(root, file))
            if dest_file not in valid_dest_files:
                is_pyc = (
                    dest_file.endswith(".pyc")
                    or dest_file.endswith(".pyo")
                    or "__pycache__" in dest_file
                )
                if not dry_run:
                    os.remove(dest_file)
                else:
                    if verbose or not is_pyc:
                        print(f"{prefix}Would remove orphan file {dest_file}")
                if not is_pyc:
                    synced_files_count += 1
        # remove empty directories
        if os.path.isdir(root) and not os.listdir(root):
            if not dry_run:
                os.rmdir(root)
            else:
                is_pycache = "__pycache__" in root
                if verbose or not is_pycache:
                    print(f"{prefix}Would remove empty directory {root}")

    if not sync_version(dry_run, prefix):
        return None

    mismatches = check_shared_pins()
    if mismatches:
        print("Warning: shared dependencies differ between pyproject.toml files:")
        print("\n".join(mismatches))

    print(
        f"{prefix}Synchronization complete! {synced_files_count} files synced/modified."
    )
    return synced_files_count


if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        print("Usage: python sync_linux_version.py [OPTIONS]")
        print(
            "\nSynchronizes the main MoleditPy source code to the Linux-specific package."
        )
        print("It transforms package names, synchronizes the pyproject.toml version,")
        print("and intelligently only modifies or copies files that have changed.")
        print("\nOptions:")
        print("  -h, --help       Show this help message and exit")
        print(
            "  --dry, --dry-run Simulate the synchronization process without modifying any files"
        )
        print(
            "  --check          Like --dry-run, but exit 1 if the Linux tree is out of date"
        )
        print(
            "  -v, --verbose    Show verbose output (e.g., include .pyc files in orphan cleanup logs)"
        )
        sys.exit(0)

    is_check = "--check" in sys.argv
    is_dry_run = is_check or "--dry" in sys.argv or "--dry-run" in sys.argv
    is_verbose = "-v" in sys.argv or "--verbose" in sys.argv

    pending = sync_linux(dry_run=is_dry_run, verbose=is_verbose)
    if pending is None:
        sys.exit(1)
    if is_check:
        if pending or check_shared_pins():
            print(
                "Linux package is out of date. Run: python scripts/sync_linux_version.py"
            )
            sys.exit(1)
        print("Linux package is up to date.")
    sys.exit(0)
