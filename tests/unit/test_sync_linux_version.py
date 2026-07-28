"""Tests for scripts/sync_linux_version.py — the main -> Linux package transform."""

import os
import sys
import importlib
import textwrap
import pytest

SCRIPTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "scripts",
)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

sync = importlib.import_module("sync_linux_version")


# =============================================================================
# transform() — package renaming
# =============================================================================


def test_dotted_module_references_are_renamed():
    assert sync.transform("from moleditpy.ui import x", False) == (
        "from moleditpy_linux.ui import x"
    )
    assert sync.transform("import moleditpy.core.mol", False) == (
        "import moleditpy_linux.core.mol"
    )


def test_bare_import_forms_are_renamed():
    assert sync.transform("import moleditpy", False) == "import moleditpy_linux"
    assert sync.transform("from moleditpy import x", False) == (
        "from moleditpy_linux import x"
    )


def test_module_path_inside_a_string_is_renamed():
    """main_window_init compares against a module name at runtime; it must track."""
    src = 'if "moleditpy.utils.constants" in name:'
    assert sync.transform(src, False) == 'if "moleditpy_linux.utils.constants" in name:'


def test_config_directory_is_not_renamed():
    """~/.moleditpy is shared between both builds."""
    src = 'os.path.join(os.path.expanduser("~"), ".moleditpy", "settings.json")'
    assert sync.transform(src, False) == src


def test_log_filename_is_preserved():
    src = 'log_path = os.path.join(log_dir, "moleditpy.log")'
    assert sync.transform(src, False) == src


def test_log_path_inside_ui_text_is_preserved():
    src = '"Save application log to ~/.moleditpy/moleditpy.log (rotated)."'
    assert sync.transform(src, False) == src


def test_windows_app_id_is_preserved():
    src = 'myappid = f"hyoko.moleditpy.{major}"'
    assert sync.transform(src, False) == src


def test_no_double_linux_suffix_is_produced():
    """The old three-rule chain produced moleditpy_linux_linux and patched it after."""
    out = sync.transform("import moleditpy.ui\nfrom moleditpy import y\n", False)
    assert "moleditpy_linux_linux" not in out
    assert out == "import moleditpy_linux.ui\nfrom moleditpy_linux import y\n"


def test_already_renamed_text_is_left_alone():
    src = "import moleditpy_linux.ui"
    assert sync.transform(src, False) == src


def test_version_lookup_package_is_renamed():
    assert sync.transform('version("MoleditPy")', False) == 'version("MoleditPy-linux")'


# =============================================================================
# transform() — Open Babel disablement
# =============================================================================


ROOT_INIT = textwrap.dedent("""\
    import importlib.util

    try:
        OBABEL_AVAILABLE = importlib.util.find_spec("openbabel") is not None
    except ImportError:
        OBABEL_AVAILABLE = False
    """)


def test_obabel_probe_is_disabled():
    out = sync.transform(ROOT_INIT, True)
    assert "find_spec" not in out
    assert out.count("OBABEL_AVAILABLE = False") == 2


def test_obabel_disablement_preserves_indentation():
    out = sync.transform(ROOT_INIT, True)
    assert "    OBABEL_AVAILABLE = False\n" in out


def test_obabel_rewrite_does_not_swallow_later_lines():
    """A greedy DOTALL `.*None` used to delete everything up to the last None."""
    src = ROOT_INIT + "\nDEFAULT_THING = None\nTRAILER = 1\n"
    out = sync.transform(src, True)
    assert "except ImportError:" in out
    assert "DEFAULT_THING = None" in out
    assert "TRAILER = 1" in out
    assert len(out.splitlines()) == len(src.splitlines())


def test_obabel_rewrite_only_applies_to_the_root_init():
    src = "OBABEL_AVAILABLE = True\n"
    assert sync.transform(src, False) == src
    assert sync.transform(src, True) == "OBABEL_AVAILABLE = False\n"


# =============================================================================
# Dependency pin comparison
# =============================================================================


def _write_toml(path, name, version, deps):
    body = ",\n".join(f'  "{d}"' for d in deps)
    path.write_text(
        f'[project]\nname = "{name}"\nversion = "{version}"\n'
        f"dependencies = [\n{body}\n]\n",
        encoding="utf-8",
    )


def test_read_dependencies_maps_name_to_requirement(tmp_path):
    toml = tmp_path / "pyproject.toml"
    _write_toml(toml, "X", "1.0", ["numpy", "pyvista < 0.49"])
    deps = sync.read_dependencies(str(toml))
    assert deps == {"numpy": "numpy", "pyvista": "pyvista < 0.49"}


def test_read_dependencies_on_missing_file_is_empty():
    assert sync.read_dependencies(r"G:\nope\pyproject.toml") == {}


def test_shared_pins_matching_reports_nothing(tmp_path, monkeypatch):
    main, linux = tmp_path / "a.toml", tmp_path / "b.toml"
    _write_toml(
        main, "MoleditPy", "4.5.0", ["numpy", "rdkit < 2026.4", "openbabel-wheel"]
    )
    _write_toml(linux, "MoleditPy-linux", "4.5.0", ["numpy", "rdkit < 2026.4"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux))
    assert sync.check_shared_pins() == []


def test_shared_pin_drift_is_reported(tmp_path, monkeypatch):
    main, linux = tmp_path / "a.toml", tmp_path / "b.toml"
    _write_toml(main, "MoleditPy", "4.5.0", ["rdkit < 2027.1"])
    _write_toml(linux, "MoleditPy-linux", "4.5.0", ["rdkit < 2026.4"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux))
    (message,) = sync.check_shared_pins()
    assert "rdkit" in message and "2027.1" in message and "2026.4" in message


def test_real_pyprojects_share_identical_pins():
    """The shipped pair must not drift; only openbabel is main-only."""
    assert sync.check_shared_pins() == []


# =============================================================================
# sync_version() failure paths
# =============================================================================


def test_sync_version_fails_when_a_pyproject_is_missing(tmp_path, monkeypatch):
    main = tmp_path / "a.toml"
    _write_toml(main, "MoleditPy", "4.5.0", ["numpy"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(tmp_path / "missing.toml"))
    assert sync.sync_version(dry_run=True, prefix="") is False


def test_sync_version_fails_when_version_is_absent(tmp_path, monkeypatch):
    main, linux = tmp_path / "a.toml", tmp_path / "b.toml"
    main.write_text('[project]\nname = "MoleditPy"\n', encoding="utf-8")
    _write_toml(linux, "MoleditPy-linux", "4.4.1", ["numpy"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux))
    assert sync.sync_version(dry_run=True, prefix="") is False


def test_sync_version_writes_the_main_version(tmp_path, monkeypatch):
    main, linux = tmp_path / "a.toml", tmp_path / "b.toml"
    _write_toml(main, "MoleditPy", "4.5.0", ["numpy"])
    _write_toml(linux, "MoleditPy-linux", "4.4.1", ["numpy"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux))
    assert sync.sync_version(dry_run=False, prefix="") is True
    assert 'version = "4.5.0"' in linux.read_text(encoding="utf-8")
    assert 'name = "MoleditPy-linux"' in linux.read_text(encoding="utf-8")


def test_sync_version_dry_run_does_not_write(tmp_path, monkeypatch):
    main, linux = tmp_path / "a.toml", tmp_path / "b.toml"
    _write_toml(main, "MoleditPy", "4.5.0", ["numpy"])
    _write_toml(linux, "MoleditPy-linux", "4.4.1", ["numpy"])
    monkeypatch.setattr(sync, "MAIN_TOML", str(main))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux))
    sync.sync_version(dry_run=True, prefix="[DRY RUN] ")
    assert 'version = "4.4.1"' in linux.read_text(encoding="utf-8")


# =============================================================================
# sync_linux() end to end against a temporary workspace
# =============================================================================


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    """A miniature main/linux pair with the module globals pointed at it."""
    src_main = tmp_path / "main" / "moleditpy"
    src_linux = tmp_path / "linux" / "moleditpy_linux"
    src_main.mkdir(parents=True)
    (src_main / "__init__.py").write_text(ROOT_INIT, encoding="utf-8")
    (src_main / "main.py").write_text("import moleditpy.ui\n", encoding="utf-8")
    (src_main / "assets").mkdir()
    (src_main / "assets" / "icon.png").write_bytes(b"\x89PNG\r\n")

    main_toml = tmp_path / "main.toml"
    linux_toml = tmp_path / "linux.toml"
    _write_toml(main_toml, "MoleditPy", "4.5.0", ["numpy"])
    _write_toml(linux_toml, "MoleditPy-linux", "4.5.0", ["numpy"])

    monkeypatch.setattr(sync, "SRC_MAIN", str(src_main))
    monkeypatch.setattr(sync, "SRC_LINUX", str(src_linux))
    monkeypatch.setattr(sync, "MAIN_TOML", str(main_toml))
    monkeypatch.setattr(sync, "LINUX_TOML", str(linux_toml))
    return src_main, src_linux


def test_missing_main_source_returns_none(workspace, monkeypatch, tmp_path):
    monkeypatch.setattr(sync, "SRC_MAIN", str(tmp_path / "gone"))
    assert sync.sync_linux() is None


def test_sync_creates_the_transformed_tree(workspace):
    src_main, src_linux = workspace
    assert sync.sync_linux() == 3
    assert (src_linux / "main.py").read_text(encoding="utf-8") == (
        "import moleditpy_linux.ui\n"
    )
    assert "find_spec" not in (src_linux / "__init__.py").read_text(encoding="utf-8")
    assert (src_linux / "assets" / "icon.png").read_bytes() == b"\x89PNG\r\n"


def test_second_sync_reports_no_changes(workspace):
    sync.sync_linux()
    assert sync.sync_linux() == 0


def test_dry_run_leaves_the_tree_untouched(workspace):
    _src_main, src_linux = workspace
    assert sync.sync_linux(dry_run=True) == 3
    assert not src_linux.exists()


def test_orphan_files_are_removed_and_counted(workspace):
    _src_main, src_linux = workspace
    sync.sync_linux()
    stale = src_linux / "removed_module.py"
    stale.write_text("# gone from main\n", encoding="utf-8")
    assert sync.sync_linux() == 1
    assert not stale.exists()


def test_pycache_is_not_counted_as_a_change(workspace):
    _src_main, src_linux = workspace
    sync.sync_linux()
    cache = src_linux / "__pycache__"
    cache.mkdir()
    (cache / "main.cpython-313.pyc").write_bytes(b"\x00")
    assert sync.sync_linux() == 0
    assert not cache.exists()


def test_source_pycache_is_not_copied(workspace):
    src_main, src_linux = workspace
    cache = src_main / "__pycache__"
    cache.mkdir()
    (cache / "main.cpython-313.pyc").write_bytes(b"\x00")
    sync.sync_linux()
    assert not (src_linux / "__pycache__").exists()


def test_sync_fails_when_the_version_cannot_be_read(workspace, monkeypatch, tmp_path):
    monkeypatch.setattr(sync, "MAIN_TOML", str(tmp_path / "absent.toml"))
    assert sync.sync_linux() is None
