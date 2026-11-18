"""Tests for running a plugin package's ``__main__`` from the command line."""

import io
import json
import sys
import uuid
from pathlib import Path
from typing import Iterator, Optional
from unittest.mock import patch

import pytest

from moleditpy.main import main
from moleditpy.plugins import plugin_runner
from moleditpy.plugins.plugin_manager import PluginManager

MAIN_SOURCE = """
import json, os, sys
with open(os.environ["RUNNER_RECORD"], "w", encoding="utf-8") as fh:
    json.dump({"argv": sys.argv, "name": __name__}, fh)
code = os.environ.get("RUNNER_EXIT")
if code is not None:
    sys.exit(int(code) if code.lstrip("-").isdigit() else code)
"""


@pytest.fixture
def manager(tmp_path: Path) -> PluginManager:
    pm = PluginManager()
    pm.plugin_dir = str(tmp_path / "plugins")
    pm.disabled_plugins = set()
    Path(pm.plugin_dir).mkdir()
    return pm


@pytest.fixture
def record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    path = tmp_path / "record.json"
    monkeypatch.setenv("RUNNER_RECORD", str(path))
    monkeypatch.delenv("RUNNER_EXIT", raising=False)
    return path


@pytest.fixture
def name() -> Iterator[str]:
    """A fresh package name, so no test sees another's cached import."""
    value = f"runner_plugin_{uuid.uuid4().hex[:8]}"
    yield value
    for module in [m for m in sys.modules if m == value or m.startswith(value + ".")]:
        del sys.modules[module]


def make_package(
    base: Path, name: str, main_source: Optional[str] = MAIN_SOURCE
) -> Path:
    folder = base / name
    folder.mkdir(parents=True)
    (folder / "__init__.py").write_text("", encoding="utf-8")
    if main_source is not None:
        (folder / "__main__.py").write_text(main_source, encoding="utf-8")
    return folder


def test_a_package_with_main_is_runnable(manager, name):
    folder = make_package(Path(manager.plugin_dir), name)
    assert plugin_runner.runnable_plugins(manager) == {name: [str(folder)]}


def test_a_package_inside_a_category_folder_is_found(manager, name):
    folder = make_package(Path(manager.plugin_dir) / "Utility", name)
    assert plugin_runner.runnable_plugins(manager) == {name: [str(folder)]}


def test_a_package_without_main_and_a_single_file_are_not_runnable(manager, name):
    make_package(Path(manager.plugin_dir), name, main_source=None)
    (Path(manager.plugin_dir) / "single.py").write_text("", encoding="utf-8")
    assert plugin_runner.runnable_plugins(manager) == {}


def test_a_disabled_plugin_is_not_runnable(manager, name):
    folder = make_package(Path(manager.plugin_dir), name)
    manager.disabled_plugins = {manager.plugin_path_key(str(folder / "__init__.py"))}
    assert plugin_runner.runnable_plugins(manager) == {}


def test_hidden_and_non_identifier_folders_are_skipped(manager):
    make_package(Path(manager.plugin_dir), ".hidden")
    make_package(Path(manager.plugin_dir), "not-a-name")
    assert plugin_runner.runnable_plugins(manager) == {}


def test_a_missing_plugin_folder_lists_nothing(tmp_path):
    pm = PluginManager()
    pm.plugin_dir = str(tmp_path / "absent")
    out = io.StringIO()
    assert plugin_runner.list_plugins(pm, out) == 0
    assert "No plugins are installed." in out.getvalue()


def plugin_lines(manager):
    out = io.StringIO()
    assert plugin_runner.list_plugins(manager, out) == 0
    return {
        line.split()[0]: line
        for line in out.getvalue().splitlines()
        if line.startswith("  ")
    }


def test_list_plugins_shows_every_plugin_and_marks_the_runnable(manager, name):
    base = Path(manager.plugin_dir)
    make_package(base, name)
    plain = make_package(base / "Utility", "plain_package", main_source=None)
    (plain / "__init__.py").write_text(
        'PLUGIN_NAME = "Plain Tool"\nPLUGIN_VERSION = "2.1"\n', encoding="utf-8"
    )
    (base / "single.py").write_text("", encoding="utf-8")
    off = make_package(base, "off_package")
    manager.disabled_plugins = {manager.plugin_path_key(str(off / "__init__.py"))}

    lines = plugin_lines(manager)
    assert set(lines) == {name, "plain_package", "single", "off_package"}
    assert "[run]" in lines[name]
    assert "[run" not in lines["plain_package"]
    assert "Plain Tool 2.1" in lines["plain_package"]
    assert "Utility" in lines["plain_package"]
    assert "[run" not in lines["single"]
    assert "[disabled]" in lines["off_package"]


def test_run_passes_arguments_and_restores_state(manager, name, record):
    make_package(Path(manager.plugin_dir), name)
    saved_argv, saved_path = sys.argv[:], sys.path[:]
    assert plugin_runner.run_plugin(name, ["gui", "--port", "1"], manager) == 0
    seen = json.loads(record.read_text(encoding="utf-8"))
    # As with ``python -m``, argv[0] is the path of the module being run.
    assert seen["argv"][0].endswith("__main__.py")
    assert seen["argv"][1:] == ["gui", "--port", "1"]
    assert seen["name"] == "__main__"
    assert sys.argv == saved_argv
    assert sys.path == saved_path


@pytest.mark.parametrize(
    ("requested", "status"), [("3", 3), ("0", 0), ("failed badly", 1)]
)
def test_the_plugins_exit_request_becomes_the_status(
    manager, name, record, monkeypatch, capsys, requested, status
):
    make_package(Path(manager.plugin_dir), name)
    monkeypatch.setenv("RUNNER_EXIT", requested)
    assert plugin_runner.run_plugin(name, [], manager) == status
    if status == 1:
        assert "failed badly" in capsys.readouterr().err


def test_an_unknown_plugin_is_an_error(manager, capsys):
    assert plugin_runner.run_plugin("nothing_here", [], manager) == 2
    assert "no runnable plugin named 'nothing_here'" in capsys.readouterr().err


def test_two_plugins_with_one_name_are_an_error(manager, name, capsys):
    make_package(Path(manager.plugin_dir) / "A", name)
    make_package(Path(manager.plugin_dir) / "B", name)
    assert plugin_runner.run_plugin(name, [], manager) == 2
    assert "more than one plugin" in capsys.readouterr().err


# Command line --------------------------------------------------------------


def run_main(argv, manager):
    with (
        patch.object(sys, "argv", ["moleditpy", *argv]),
        patch.object(plugin_runner, "PluginManager", return_value=manager),
        pytest.raises(SystemExit) as exit_info,
    ):
        main()
    return exit_info.value.code


def test_run_plugin_flag_hands_later_options_to_the_plugin(
    manager, name, record, capsys
):
    make_package(Path(manager.plugin_dir), name)
    argv = ["-r", name, "serve", "--config", "x.json", "--safe"]
    assert run_main(argv, manager) == 0
    assert json.loads(record.read_text(encoding="utf-8"))["argv"][1:] == [
        "serve",
        "--config",
        "x.json",
        "--safe",
    ]
    assert "MoleditPy - A Python-based" in capsys.readouterr().out


def test_list_plugins_flag_prints_the_banner_and_works_in_safe_mode(
    manager, name, capsys
):
    make_package(Path(manager.plugin_dir), name)
    assert run_main(["--safe", "--list-plugins"], manager) == 0
    out = capsys.readouterr().out
    assert out.startswith("-----")
    assert "MoleditPy - A Python-based" in out
    assert name in out


def test_run_plugin_without_a_name_is_an_error(manager, capsys):
    assert run_main(["--run-plugin"], manager) == 2
    assert "needs a plugin name" in capsys.readouterr().err


def test_safe_mode_cannot_run_a_plugin(manager, name, capsys):
    make_package(Path(manager.plugin_dir), name)
    assert run_main(["--safe", "-r", name], manager) == 2
    assert "--safe" in capsys.readouterr().err


def test_a_relative_argument_resolves_from_the_callers_folder(
    manager, name, record, tmp_path, monkeypatch
):
    make_package(Path(manager.plugin_dir), name)
    caller = tmp_path / "caller"
    caller.mkdir()
    monkeypatch.chdir(caller)
    (Path(manager.plugin_dir) / name / "__main__.py").write_text(
        MAIN_SOURCE.replace('"argv": sys.argv', '"argv": [os.getcwd()]'),
        encoding="utf-8",
    )
    assert plugin_runner.run_plugin(name, ["input.inp"], manager) == 0
    assert json.loads(record.read_text(encoding="utf-8"))["argv"] == [str(caller)]


@pytest.mark.parametrize(
    "argv",
    [["-reverse"], ["mol.mol", "-reverse"], ["-reverse", "--safe", "mol.mol"]],
)
def test_qt_reverse_is_not_read_as_run_plugin(argv):
    from moleditpy.main import _split_qt_options

    ours, qt = _split_qt_options(argv)
    assert qt == ["-reverse"]
    assert "-reverse" not in ours


def test_reverse_after_run_plugin_belongs_to_the_plugin():
    from moleditpy.main import _split_qt_options

    argv = ["--safe", "-r", "tool", "-reverse"]
    assert _split_qt_options(argv) == (argv, [])


def test_reverse_reaches_qt_and_opens_the_editor(monkeypatch):
    seen = {}

    class StopAtQt(Exception):
        pass

    def fake_qapplication(argv):
        seen["argv"] = argv
        raise StopAtQt

    monkeypatch.setattr(sys, "argv", ["moleditpy", "-reverse"])
    monkeypatch.setattr("moleditpy.main.QApplication", fake_qapplication)
    with pytest.raises(StopAtQt):
        main()
    assert seen["argv"][1:] == ["-reverse"]
