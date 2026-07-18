"""Run unit + integration + e2e + GUI tests with combined coverage and generate coverage_report.md."""

import subprocess
import json
import re
import sys
import os
import platform
import argparse

# Adjusted for tests/utils/ location
UTILS = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(UTILS)  # tests/
ROOT = os.path.dirname(BASE)  # project root
SRC = os.path.join(ROOT, "moleditpy", "src")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--skip-run",
    action="store_true",
    help="Skip running tests and use existing .coverage data",
)
args_cov = parser.parse_known_args()[0]

env = os.environ.copy()
pp = env.get("PYTHONPATH", "")
env["PYTHONPATH"] = f"{SRC}{os.pathsep}{pp}".strip(os.pathsep)

devnull = subprocess.DEVNULL

# Per-suite test counts parsed from each pytest summary line.
suite_counts = {}
_counts_cache = os.path.join(BASE, ".test_counts.json")


def _parse_counts(text):
    """Extract passed/failed/skipped/error counts from a pytest summary line."""

    def grab(word):
        m = re.search(rf"(\d+) {word}", text)
        return int(m.group(1)) if m else 0

    return {
        "passed": grab("passed"),
        "failed": grab("failed"),
        "skipped": grab("skipped"),
        "errors": grab("error"),
    }


def _run_counting(label, cmd, run_env):
    """Run a pytest suite, record its test counts, and return the exit code."""
    result = subprocess.run(cmd, env=run_env, cwd=ROOT, capture_output=True, text=True)
    counts = _parse_counts(result.stdout + result.stderr)
    suite_counts[label] = counts
    extra = ""
    if counts["skipped"]:
        extra += f", {counts['skipped']} skipped"
    if counts["failed"]:
        extra += f", {counts['failed']} FAILED"
    print(f"  -> {label}: {counts['passed']} passed{extra}")
    return result.returncode


if not args_cov.skip_run:
    # Ensure we are in ROOT when running coverage
    print("Running Unit Tests with coverage (headless)...")
    gui_env = env.copy()
    gui_env["MOLEDITPY_HEADLESS"] = "1"
    gui_env["QT_QPA_PLATFORM"] = "offscreen"
    _run_counting(
        "Unit",
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/unit",
            "-q",
            "--cov=moleditpy",
            "--cov-report=",
            "--tb=short",
        ],
        gui_env,
    )

    print("Running Integration Tests with coverage (headless)...")
    _run_counting(
        "Integration",
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/integration",
            "-q",
            "--cov=moleditpy",
            "--cov-append",
            "--cov-report=",
            "--tb=short",
        ],
        gui_env,
    )

    print("Running E2E Tests with coverage (headless)...")
    # On Linux, e2e imports moleditpy_linux; everywhere else it uses moleditpy
    e2e_cov_source = "moleditpy_linux" if platform.system() == "Linux" else "moleditpy"
    # Auto-sync Linux package before running e2e coverage
    if platform.system() == "Linux":
        sync_script = os.path.join(ROOT, "scripts", "sync_linux_version.py")
        if os.path.exists(sync_script):
            subprocess.run(
                [sys.executable, sync_script],
                cwd=ROOT,
                stdout=devnull,
                stderr=devnull,
            )
    _run_counting(
        "E2E",
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/e2e",
            "-q",
            "--timeout=60",
            f"--cov={e2e_cov_source}",
            "--cov-append",
            "--cov-report=",
            "--tb=short",
        ],
        gui_env,
    )

    print("Running GUI Tests with coverage (headless)...")
    _run_counting(
        "GUI",
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/gui",
            "-q",
            "--timeout=60",
            "--cov=moleditpy",
            "--cov-append",
            "--cov-report=",
            "--tb=short",
        ],
        gui_env,
    )

    # Cache counts so a later --skip-run invocation can still report them.
    try:
        with open(_counts_cache, "w", encoding="utf-8") as fh:
            json.dump(suite_counts, fh)
    except OSError:
        pass
else:
    print("Skipping test execution, using existing coverage data...")
    try:
        with open(_counts_cache, "r", encoding="utf-8") as fh:
            suite_counts = json.load(fh)
    except (OSError, ValueError):
        suite_counts = {}

# Step 4: Generate reports
print("Generating coverage reports...")

cvrc = os.path.join(ROOT, ".coveragerc")
subprocess.run(
    [
        sys.executable,
        "-m",
        "coverage",
        "json",
        "--rcfile",
        cvrc,
        "-o",
        os.path.join(BASE, "cov_full.json"),
    ],
    env=env,
    cwd=ROOT,
    stdout=devnull,
    stderr=devnull,
    check=True,
)

# Generate HTML report for full application coverage
subprocess.run(
    [
        sys.executable,
        "-m",
        "coverage",
        "html",
        "--rcfile",
        cvrc,
        "-d",
        os.path.join(BASE, "coverage_html"),
    ],
    env=env,
    cwd=ROOT,
    stdout=devnull,
    stderr=devnull,
)


# Step 5: Build merged markdown report
def get_stats(path):
    with open(path, "r") as f:
        data = json.load(f)
    return data["totals"], data["files"]


def build_table(files, totals, title):
    lines = []
    lines.append(f"### {title}")
    lines.append("")
    lines.append("| File | Stmts | Miss | Cover |")
    lines.append("| :--- | :--- | :--- | :--- |")

    files_by_dir = {}
    for fname, info in sorted(files.items()):
        s = info["summary"]
        if s["num_statements"] == 0:
            continue
        d = os.path.dirname(fname) or "."
        files_by_dir.setdefault(d, []).append((fname, s))

    for d in sorted(files_by_dir.keys()):
        for fname, s in files_by_dir[d]:
            p = s["percent_covered"]
            ml = s["missing_lines"]
            lines.append(
                f"| {fname:<55} | {s['num_statements']:>6} | {ml:>6} | {p:>6.1f}% |"
            )

    lines.append(
        f"| **TOTAL** | **{totals['num_statements']}** | **{totals['missing_lines']}** | **{totals['percent_covered']:.2f}%** |"
    )
    lines.append("")
    return lines


full_totals, full_files = get_stats(os.path.join(BASE, "cov_full.json"))

markdown_lines = []
markdown_lines.append("# MoleditPy Coverage Report")
markdown_lines.append("")
markdown_lines.append(
    f"- **Overall Project Coverage**: **{full_totals['percent_covered']:.2f}%**"
)
markdown_lines.append("")

markdown_lines.extend(build_table(full_files, full_totals, "Coverage Breakdown"))


def _status_line(label, key):
    counts = suite_counts.get(key)
    if not counts:
        return f"- **{label}**: PASSED"
    passed = counts.get("passed", 0)
    failed = counts.get("failed", 0)
    skipped = counts.get("skipped", 0)
    status = "FAILED" if failed else "PASSED"
    detail = f"{passed} passed"
    if skipped:
        detail += f", {skipped} skipped"
    if failed:
        detail += f", {failed} failed"
    return f"- **{label}**: {status} ({detail})"


total_passed = sum(c.get("passed", 0) for c in suite_counts.values())
total_skipped = sum(c.get("skipped", 0) for c in suite_counts.values())

markdown_lines.append("## Test Suite Status")
if total_passed:
    total_line = f"- **Total tests passed**: {total_passed}"
    if total_skipped:
        total_line += f" ({total_skipped} skipped)"
    markdown_lines.append(total_line)
markdown_lines.append(_status_line("Unit tests", "Unit"))
markdown_lines.append(_status_line("Integration tests", "Integration"))
markdown_lines.append(_status_line("E2E tests", "E2E"))
markdown_lines.append(_status_line("GUI tests", "GUI"))
markdown_lines.append("")
markdown_lines.append("[View Detailed HTML Report](coverage_html/index.html)")

output_file = os.path.join(BASE, "coverage_report.md")
with open(output_file, "w", encoding="utf-8") as f:
    f.write("\n".join(markdown_lines))

print(f"Coverage report written to: {output_file}")

# Step 6: Update Assertion Catalog automatically
print("Updating Assertion Catalog...")
subprocess.run(
    [sys.executable, os.path.join(UTILS, "generate_assertion_catalog.py")], cwd=ROOT
)
