"""Runner failures must stay visible instead of being retried into a pass."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import Mock

import pytest


@pytest.fixture
def runner(monkeypatch: pytest.MonkeyPatch):
    for name in ("PYTEST_ADDOPTS", "NO_COLOR"):
        monkeypatch.setenv(name, os.environ.get(name, ""))
    path = Path(__file__).resolve().parents[2] / "run_all_tests.py"
    spec = importlib.util.spec_from_file_location("test_suite_runner", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("returncode", [1, 2, 3, 4, 5])
def test_pytest_failure_is_not_retried(runner, monkeypatch, returncode: int) -> None:
    run_once = Mock(return_value=(returncode, "1 failed", 0.1))
    monkeypatch.setattr(runner, "_run_once", run_once)
    assert runner.run_suite("UNIT", "unused", enable_cov=False) == returncode
    run_once.assert_called_once()


def test_known_native_crash_can_still_retry(runner, monkeypatch) -> None:
    run_once = Mock(side_effect=[(-11, "", 0.1), (0, "1 passed", 0.1)])
    monkeypatch.setattr(runner, "_run_once", run_once)
    assert runner.run_suite("UNIT", "unused", enable_cov=False) == 0
    assert run_once.call_count == 2


def test_fail_fast_disables_native_crash_retry(runner, monkeypatch) -> None:
    run_once = Mock(return_value=(-11, "", 0.1))
    monkeypatch.setattr(runner, "_run_once", run_once)
    assert runner.run_suite("UNIT", "unused", enable_cov=False, exitfirst=True) == -11
    run_once.assert_called_once()


@pytest.mark.parametrize(
    "summary", ["=== 23 passed, 2 skipped in 1.0s ===", "23 passed, 2 skipped in 1.0s"]
)
def test_counts_ignore_status_words_in_test_ids(runner, summary: str) -> None:
    output = "tests/test_example.py::test_status[1 passed-1 failed] PASSED\n" + summary
    assert runner._parse_test_counts(output) == {
        "passed": 23,
        "failed": 0,
        "skipped": 2,
    }


@pytest.mark.parametrize("returncode, summary", [(0, "1 passed"), (1, "1 failed")])
def test_full_gui_cli_exit_status(
    tmp_path: Path, returncode: int, summary: str
) -> None:
    source = Path(__file__).resolve().parents[2] / "run_all_tests.py"
    copied_runner = tmp_path / "tests/run_all_tests.py"
    copied_runner.parent.mkdir()
    copied_runner.write_bytes(source.read_bytes())
    (tmp_path / "pytest.py").write_text(
        f"print({summary!r})\nraise SystemExit({returncode})\n", encoding="utf-8"
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = str(tmp_path)
    result = subprocess.run(
        [sys.executable, str(copied_runner), "--full-gui", "--no-cov", "--no-report"],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == returncode, result.stdout + result.stderr
    status = "PASSED" if returncode == 0 else "FAILED"
    assert f"FULL_GUI    : {status}" in result.stdout
