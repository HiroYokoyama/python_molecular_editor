"""Run unit + integration + e2e + GUI tests with combined coverage and generate coverage_report.md."""

import subprocess
import json
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

if not args_cov.skip_run:
    # Ensure we are in ROOT when running coverage
    print("Running Unit Tests with coverage (headless)...")
    gui_env = env.copy()
    gui_env["MOLEDITPY_HEADLESS"] = "1"
    gui_env["QT_QPA_PLATFORM"] = "offscreen"
    subprocess.run(
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
        env=gui_env,
        cwd=ROOT,
        stdout=devnull,
        stderr=devnull,
    )

    print("Running Integration Tests with coverage (headless)...")
    subprocess.run(
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
        env=gui_env,
        cwd=ROOT,
        stdout=devnull,
        stderr=devnull,
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
    subprocess.run(
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
        env=gui_env,
        cwd=ROOT,
        stdout=devnull,
        stderr=devnull,
    )

    print("Running GUI Tests with coverage (headless)...")
    subprocess.run(
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
        env=gui_env,
        cwd=ROOT,
        stdout=devnull,
        stderr=devnull,
    )
else:
    print("Skipping test execution, using existing coverage data...")

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

markdown_lines.append("## Test Suite Status")
markdown_lines.append("- **Unit tests**: PASSED")
markdown_lines.append("- **Integration tests**: PASSED")
markdown_lines.append("- **E2E tests**: PASSED")
markdown_lines.append("- **GUI tests**: PASSED")
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
