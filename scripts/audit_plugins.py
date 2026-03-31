"""
Plugin Audit Tool for MoleditPy
Checks all plugin entry-point files for broken run() stubs.

Rules:
  - If run() is present, the framework auto-registers a menu entry.
    run() just must not be a stub (pass).
  - initialize(context) is optional — used for advanced registration
    (sub-menus, export actions, 3D styles, etc.)
  - A plugin is unreachable only if BOTH run() and initialize() are absent.

Usage: python tools/audit_plugins.py [plugin_dir]
"""
import ast
import os
import sys

PLUGIN_DIR = sys.argv[1] if len(sys.argv) > 1 else "tmp/plugins"

# Sub-package helper dirs (not plugin entry points)
SKIP_DIRS = {
    'gaussian_fchk_mo_analyzer', 'orca_result_analyzer',
    'orca_input_generator_pro', 'reaction_sketcher',
    'nmr_predicator_nmrshiftdb2', '__pycache__',
}

# Keywords that indicate run() does real work
RUN_ACTION_KEYWORDS = [
    'run_plugin', '_launch_fn', 'dialog', 'context',
    'launch', 'export', 'open_settings', 'identify',
    'resolve', 'plugin_manager', 'QFileDialog', 'win',
]


def get_func_node(tree, name):
    for n in tree.body:
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return n
    return None


def check_file(path, plugin_dir):
    rel = os.path.relpath(path, plugin_dir)
    with open(path, encoding='utf-8', errors='replace') as f:
        src = f.read()

    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        return rel, [f"SYNTAX ERROR: {e}"]

    issues = []
    run_node  = get_func_node(tree, 'run')
    init_node = get_func_node(tree, 'initialize')

    if run_node:
        seg = ast.get_source_segment(src, run_node) or ''
        body = run_node.body

        is_pure_pass = len(body) == 1 and isinstance(body[0], ast.Pass)
        # host-redirect then pass: if hasattr(mw,'host'): ... \n pass
        has_host_then_pass = (
            len(body) == 2
            and isinstance(body[1], ast.Pass)
        )
        if is_pure_pass or has_host_then_pass:
            issues.append("run() is a stub (pass) -- does nothing")
        elif not any(kw in seg for kw in RUN_ACTION_KEYWORDS):
            issues.append("run() has no recognizable action -- review manually")

    elif not init_node:
        issues.append("No run() and no initialize() -- plugin is unreachable")

    return rel, issues


results = {}
for root, dirs, files in os.walk(PLUGIN_DIR):
    dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
    for fname in sorted(files):
        if not fname.endswith('.py') or fname.startswith('_'):
            continue
        fpath = os.path.join(root, fname)
        rel, issues = check_file(fpath, PLUGIN_DIR)
        results[rel] = issues

broken = {k: v for k, v in results.items() if v}
ok     = [k for k, v in results.items() if not v]

print("=" * 60)
print(f"Plugin Audit -- {PLUGIN_DIR}")
print("=" * 60)
if broken:
    print(f"\n[ISSUES] {len(broken)} file(s):\n")
    for path, issues in sorted(broken.items()):
        print(f"  {path}")
        for iss in issues:
            print(f"    ! {iss}")
else:
    print("\nNo issues found.")

print(f"\n[OK] {len(ok)} file(s) passed.\n")
for p in ok:
    print(f"  {p}")
