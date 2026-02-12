"""Run unit + integration + GUI tests with combined coverage and print a comprehensive report."""
import subprocess, json, sys, os

BASE = os.path.dirname(os.path.abspath(__file__))  # tests/
ROOT = os.path.dirname(BASE)  # project root
SRC = os.path.join(ROOT, "moleditpy", "src")

env = os.environ.copy()
pp = env.get("PYTHONPATH", "")
env["PYTHONPATH"] = f"{SRC}{os.pathsep}{pp}".strip(os.pathsep)

# Step 1: Run unit tests with coverage
print("=" * 60)
print("STEP 1: Running Unit Tests with coverage...")
print("=" * 60)
r1 = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/unit", "-q",
     "--cov=moleditpy", "--cov-report=", "--tb=short"],
    env=env, cwd=ROOT
)
unit_ok = r1.returncode == 0

# Step 2: Run integration tests, appending to existing coverage data
print()
print("=" * 60)
print("STEP 2: Running Integration Tests with coverage...")
print("=" * 60)
r_int = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/integration", "-q",
     "--cov=moleditpy", "--cov-append", "--cov-report=", "--tb=short"],
    env=env, cwd=ROOT
)
int_ok = r_int.returncode == 0

# Step 3: Run GUI tests, appending to existing coverage data
print()
print("=" * 60)
print("STEP 3: Running GUI Tests with coverage (headless)...")
print("=" * 60)
gui_env = env.copy()
gui_env["MOLEDITPY_HEADLESS"] = "1"
gui_env["QT_QPA_PLATFORM"] = "offscreen"
r2 = subprocess.run(
    [sys.executable, "-m", "pytest", "tests/gui", "-q", "--timeout=60",
     "--cov=moleditpy", "--cov-append", "--cov-report=", "--tb=short"],
    env=gui_env, cwd=ROOT
)
gui_ok = r2.returncode == 0

# Step 4: Generate JSON report from combined .coverage
print()
print("=" * 60)
print("STEP 4: Generating combined coverage report...")
print("=" * 60)
subprocess.run(
    [sys.executable, "-m", "coverage", "json", "-o", os.path.join(BASE, "cov_combined.json")],
    env=env, cwd=ROOT
)
subprocess.run(
    [sys.executable, "-m", "coverage", "html", "-d", os.path.join(BASE, "coverage_html")],
    env=env, cwd=ROOT
)

# Step 5: Print the report
with open(os.path.join(BASE, "cov_combined.json"), "r") as f:
    data = json.load(f)

totals = data["totals"]
pct = totals["percent_covered"]
covered = totals["covered_lines"]
stmts = totals["num_statements"]
missing = stmts - covered

print("# Combined Coverage Report")
print("\n**Suites**: Unit + Integration + GUI (Headless)")
print(f"\n- **Total Statements**: {stmts}")
print(f"- **Total Covered**: {covered}")
print(f"- **Overall Coverage**: **{pct:.2f}%**")

print("\n## File Breakdown")
print("\n| File | Stmts | Miss | Cover |")
print("| :--- | :--- | :--- | :--- |")

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
        print(f"| {fname:<55} | {s['num_statements']:>6} | {ml:>6} | {p:>6.1f}% |")

print(f"| **TOTAL** | **{stmts}** | **{missing}** | **{pct:.1f}%** |")

print("\n## Test Suite Status")
print(f"- **Unit tests**: {'✅ PASSED' if unit_ok else '❌ FAILED'}")
print(f"- **Integration tests**: {'✅ PASSED' if int_ok else '❌ FAILED'}")
print(f"- **GUI tests**: {'✅ PASSED' if gui_ok else '❌ FAILED'}")
print(f"\n[View HTML Report](coverage_html/index.html)")
