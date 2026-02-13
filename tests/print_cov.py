"""Run unit + integration + GUI tests with combined coverage and generate coverage_report.md."""
import subprocess, json, sys, os

BASE = os.path.dirname(os.path.abspath(__file__))  # tests/
ROOT = os.path.dirname(BASE)  # project root
SRC = os.path.join(ROOT, "moleditpy", "src")

env = os.environ.copy()
pp = env.get("PYTHONPATH", "")
env["PYTHONPATH"] = f"{SRC}{os.pathsep}{pp}".strip(os.pathsep)

# Suppress pytest output by redirecting to devnull
devnull = subprocess.DEVNULL

# Step 1: Run unit tests with coverage
print("Running Unit Tests with coverage...")
r1 = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/unit", "-q",
     "--cov=moleditpy", "--cov-report=", "--tb=short"],
    env=env, cwd=ROOT, stdout=devnull, stderr=devnull
)
unit_ok = r1.returncode == 0

# Step 2: Run integration tests, appending to existing coverage data
print("Running Integration Tests with coverage...")
r_int = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/integration", "-q",
     "--cov=moleditpy", "--cov-append", "--cov-report=", "--tb=short"],
    env=env, cwd=ROOT, stdout=devnull, stderr=devnull
)
int_ok = r_int.returncode == 0

# Step 3: Run GUI tests, appending to existing coverage data
print("Running GUI Tests with coverage (headless)...")
gui_env = env.copy()
gui_env["MOLEDITPY_HEADLESS"] = "1"
gui_env["QT_QPA_PLATFORM"] = "offscreen"
r2 = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/gui", "-q", "--timeout=60",
     "--cov=moleditpy", "--cov-append", "--cov-report=", "--tb=short"],
    env=gui_env, cwd=ROOT, stdout=devnull, stderr=devnull
)
gui_ok = r2.returncode == 0

# Step 4: Generate JSON report from combined .coverage
print("Generating combined coverage report...")
subprocess.run(
    [sys.executable, "-m", "coverage", "json", "-o", os.path.join(BASE, "cov_combined.json")],
    env=env, cwd=ROOT, stdout=devnull, stderr=devnull
)
subprocess.run(
    [sys.executable, "-m", "coverage", "html", "-d", os.path.join(BASE, "coverage_html")],
    env=env, cwd=ROOT, stdout=devnull, stderr=devnull
)

# Step 5: Build the markdown report
with open(os.path.join(BASE, "cov_combined.json"), "r") as f:
    data = json.load(f)

totals = data["totals"]
pct = totals["percent_covered"]
covered = totals["covered_lines"]
stmts = totals["num_statements"]
missing = stmts - covered

markdown_lines = []
markdown_lines.append("# Combined Coverage Report")
markdown_lines.append("")
markdown_lines.append("**Suites**: Unit + Integration + GUI (Headless)")
markdown_lines.append("")
markdown_lines.append(f"- **Total Statements**: {stmts}")
markdown_lines.append(f"- **Total Covered**: {covered}")
markdown_lines.append(f"- **Overall Coverage**: **{pct:.2f}%**")
markdown_lines.append("")
markdown_lines.append("## File Breakdown")
markdown_lines.append("")
markdown_lines.append("| File | Stmts | Miss | Cover |")
markdown_lines.append("| :--- | :--- | :--- | :--- |")

# Group by directory for readability
files_by_dir = {}
for fname, info in sorted(data["files"].items()):
    s = info["summary"]
    if s["num_statements"] == 0:
        continue
    parts = fname.replace("\\", "/").split("/")
    if len(parts) > 1:
        d = "/".join(parts[:-1])
    else:
        d = "."
    files_by_dir.setdefault(d, []).append((fname, s))

for d in sorted(files_by_dir.keys()):
    entries = files_by_dir[d]
    for fname, s in entries:
        p = s["percent_covered"]
        ml = s["missing_lines"]
        markdown_lines.append(f"| {fname:<55} | {s['num_statements']:>6} | {ml:>6} | {p:>6.1f}% |")

markdown_lines.append(f"| **TOTAL** | **{stmts}** | **{missing}** | **{pct:.2f}%** |")
markdown_lines.append("")
markdown_lines.append("## Test Suite Status")
markdown_lines.append(f"- **Unit tests**: {'PASSED' if unit_ok else 'FAILED'}")
markdown_lines.append(f"- **Integration tests**: {'PASSED' if int_ok else 'FAILED'}")
markdown_lines.append(f"- **GUI tests**: {'PASSED' if gui_ok else 'FAILED'}")
markdown_lines.append(f"- **Overall Coverage**: **{pct:.2f}%**")
markdown_lines.append("")
markdown_lines.append("[View HTML Report](coverage_html/index.html)")
markdown_lines.append("")

# Write to coverage_report.md
output_file = os.path.join(BASE, "coverage_report.md")
with open(output_file, "w", encoding="utf-8") as f:
    f.write("\n".join(markdown_lines))

print(f"Coverage report written to: {output_file}")
print(f"Overall Coverage: {pct:.2f}%")
print(f"Unit: {'PASSED' if unit_ok else 'FAILED'} | Integration: {'PASSED' if int_ok else 'FAILED'} | GUI: {'PASSED' if gui_ok else 'FAILED'}")
