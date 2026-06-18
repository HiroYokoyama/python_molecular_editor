import os
import shutil
import re
import sys
import filecmp

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_MAIN = os.path.join(BASE_DIR, "moleditpy", "src", "moleditpy")
SRC_LINUX = os.path.join(BASE_DIR, "moleditpy-linux", "src", "moleditpy_linux")


def sync_linux(dry_run=False, verbose=False):
    prefix = "[DRY RUN] " if dry_run else ""
    print(f"{prefix}Synchronizing Linux version from {SRC_MAIN} to {SRC_LINUX}...")

    if not os.path.exists(SRC_MAIN):
        print(f"Error: Main source directory not found at {SRC_MAIN}")
        return

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
        dest_dir = os.path.join(SRC_LINUX, rel_path)

        if not dry_run:
            os.makedirs(dest_dir, exist_ok=True)

        for file in files:
            if file.endswith(".pyc") or file.endswith(".pyo") or "__pycache__" in root:
                continue

            src_file = os.path.join(root, file)
            dest_file = os.path.join(dest_dir, file)
            valid_dest_files.add(os.path.abspath(dest_file))

            # For Python files, perform transformation
            if file.endswith(".py"):
                try:
                    with open(src_file, "r", encoding="utf-8") as f:
                        content = f.read()

                    # Transformation 1: Rename package
                    content = content.replace("moleditpy.", "moleditpy_linux.")
                    content = content.replace(
                        "from moleditpy import", "from moleditpy_linux import"
                    )
                    content = content.replace(
                        "import moleditpy", "import moleditpy_linux"
                    )

                    # Transformation 2: Fix package name for version checking
                    content = content.replace(
                        'version("MoleditPy")', 'version("MoleditPy-linux")'
                    )

                    # Transformation 3: Explicit Open Babel disablement
                    if file == "__init__.py" and rel_path == ".":
                        content = re.sub(
                            r"OBABEL_AVAILABLE = (True|False|importlib.*None)",
                            "OBABEL_AVAILABLE = False",
                            content,
                            flags=re.DOTALL,
                        )

                    # Transformation 4: Fix double mentions
                    content = content.replace(
                        "moleditpy_linux_linux", "moleditpy_linux"
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
                if not dry_run:
                    os.remove(dest_file)
                else:
                    is_pyc = (
                        dest_file.endswith(".pyc")
                        or dest_file.endswith(".pyo")
                        or "__pycache__" in dest_file
                    )
                    if verbose or not is_pyc:
                        print(f"{prefix}Would remove orphan file {dest_file}")
        # remove empty directories
        if not os.listdir(root):
            if not dry_run:
                os.rmdir(root)
            else:
                is_pycache = "__pycache__" in root
                if verbose or not is_pycache:
                    print(f"{prefix}Would remove empty directory {root}")

    # 3. Synchronize version in pyproject.toml
    main_toml_path = os.path.join(BASE_DIR, "moleditpy", "pyproject.toml")
    linux_toml_path = os.path.join(BASE_DIR, "moleditpy-linux", "pyproject.toml")

    if os.path.exists(main_toml_path) and os.path.exists(linux_toml_path):
        with open(main_toml_path, "r", encoding="utf-8") as f:
            main_toml = f.read()

        # Extract version
        version_match = re.search(
            r'^version\s*=\s*["\']([^"\']+)["\']', main_toml, re.MULTILINE
        )
        if version_match:
            new_version = version_match.group(1)

            with open(linux_toml_path, "r", encoding="utf-8") as f:
                linux_toml = f.read()

            new_linux_toml = re.sub(
                r'^version\s*=\s*["\'][^"\']+["\']',
                f'version = "{new_version}"',
                linux_toml,
                flags=re.MULTILINE,
            )

            if new_linux_toml != linux_toml:
                if not dry_run:
                    with open(linux_toml_path, "w", encoding="utf-8") as f:
                        f.write(new_linux_toml)
                print(f"{prefix}Synchronized pyproject.toml version to {new_version}")

    print(
        f"{prefix}Synchronization complete! {synced_files_count} files synced/modified."
    )


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
            "  -v, --verbose    Show verbose output (e.g., include .pyc files in orphan cleanup logs)"
        )
        sys.exit(0)

    is_dry_run = "--dry" in sys.argv or "--dry-run" in sys.argv
    is_verbose = "-v" in sys.argv or "--verbose" in sys.argv
    sync_linux(dry_run=is_dry_run, verbose=is_verbose)
