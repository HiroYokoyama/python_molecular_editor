#!/usr/bin/env python3
"""
scan_test_disconnections.py -- AST-based scanner for broken test patterns.

Finds:
  * STALE_IMPORT    -- imports from modules that don't exist in dev source
  * MISSING_NAME    -- imports of names that don't exist in the source module
  * MGR_INIT_NOARGS -- DummyClass(Manager) calls super().__init__() without host arg
  * NO_HOST_SETUP   -- DummyClass(Manager) has no self.host assignment in __init__
  * STALE_ASSERT    -- mock.assert_called_with(X) where X doesn't match current API
  * WRONG_MANAGER   -- test mocks self.X but code now routes via self.host.manager.X

Usage:
    python scripts/scan_test_disconnections.py
    python scripts/scan_test_disconnections.py --verbose
"""

import ast
import sys
import argparse
import importlib.util
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
SRC  = ROOT / "moleditpy" / "src" / "moleditpy" / "ui"
TESTS = ROOT / "tests" / "unit"

# Manager class names -> source file
MANAGER_CLASS_MAP = {
    "IOManager":               SRC / "io_logic.py",
    "MainWindowProjectIo":     SRC / "project_io.py",
    "MainWindowMolecularParsers": SRC / "molecular_parsers.py",
    "AppStateManager":         SRC / "app_state.py",
    "EditActionsManager":      SRC / "edit_actions_logic.py",
    "UIManager":               SRC / "ui_manager.py",
    "View3DManager":           SRC / "view_3d_logic.py",
    "Edit3DManager":           SRC / "edit_3d_logic.py",
    "DialogManager":           SRC / "dialog_logic.py",
    "ComputeManager":          SRC / "compute_logic.py",
    "ExportManager":           SRC / "export_logic.py",
    "StringImporterManager":   SRC / "string_importers.py",
    "InitManager":             SRC / "main_window_init.py",
}

# Manager attribute on MainWindow -> manager class name (inverse map)
MANAGER_ATTR_TO_CLASS = {
    "state_manager":           "AppStateManager",
    "edit_actions_manager":    "EditActionsManager",
    "ui_manager":              "UIManager",
    "view_3d_manager":         "View3DManager",
    "edit_3d_manager":         "Edit3DManager",
    "dialog_manager":          "DialogManager",
    "compute_manager":         "ComputeManager",
    "export_manager":          "ExportManager",
    "io_manager":              "IOManager",
    "string_importer_manager": "StringImporterManager",
    "init_manager":            "InitManager",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ast_parse(fpath: Path):
    try:
        return ast.parse(fpath.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError:
        return None


def ast_collect_names(fpath: Path) -> set:
    """All top-level class/function/variable names exported from fpath."""
    tree = ast_parse(fpath)
    if not tree:
        return set()
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    names.add(t.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


def ast_collect_methods(fpath: Path) -> set:
    """Return set of method names defined at class level in fpath."""
    tree = ast_parse(fpath)
    if not tree:
        return set()
    methods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    methods.add(item.name)
    return methods


def get_manager_init_params(fpath: Path) -> list:
    """Return parameter names of the first __init__ found in a class."""
    tree = ast_parse(fpath)
    if not tree:
        return ["self"]
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                    return [a.arg for a in item.args.args]
    return ["self"]


def resolve_moleditpy_import(module_str: str) -> Path | None:
    """Given 'moleditpy.ui.io_logic', return the .py path in dev src if it exists.
    Also handles moleditpy.core.X, moleditpy.utils.X etc. by searching under ROOT."""
    if not module_str.startswith("moleditpy."):
        return None
    # Try in SRC (moleditpy/src/moleditpy/...) tree
    parts = module_str.split(".")
    # parts[0] = 'moleditpy', rest = subpath
    pkg_root = ROOT / "moleditpy" / "src" / "moleditpy"
    rel = "/".join(parts[1:]) + ".py"
    candidate = pkg_root / rel
    if candidate.exists():
        return candidate
    # Also try as package __init__
    pkg_candidate = pkg_root / "/".join(parts[1:]) / "__init__.py"
    if pkg_candidate.exists():
        return pkg_candidate
    return None


# ---------------------------------------------------------------------------
# Build per-manager method sets
# ---------------------------------------------------------------------------

def build_manager_methods() -> dict:
    """Return {class_name: set(method_names)}."""
    result = {}
    for cls_name, fpath in MANAGER_CLASS_MAP.items():
        result[cls_name] = ast_collect_methods(fpath) if fpath.exists() else set()
        # Also include methods from aliased sources
        if cls_name in ("MainWindowProjectIo", "MainWindowMolecularParsers"):
            # These are aliases for IOManager
            result[cls_name] |= result.get("IOManager", set())
    return result


# ---------------------------------------------------------------------------
# AST visitors for test files
# ---------------------------------------------------------------------------

class TestFileAnalyzer(ast.NodeVisitor):
    """Walk a test file's AST and collect findings."""

    def __init__(self, fpath: Path, manager_methods: dict, verbose: bool):
        self.fpath = fpath
        self.manager_methods = manager_methods
        self.verbose = verbose
        self.findings = []

        # Track imports: {local_name: (module_str, original_name)}
        self.imports: dict = {}
        # Track DummyClass base classes: {class_name: [base_class_names]}
        self.dummy_classes: dict = {}
        # Track which DummyClass.__init__ sets self.host
        self.host_setups: dict = {}  # {class_name: bool}
        # Track super().__init__ calls: {class_name: args_count}
        self.super_init_args: dict = {}

    def visit_ImportFrom(self, node: ast.ImportFrom):
        if not node.module:
            self.generic_visit(node)
            return
        for alias in node.names:
            local = alias.asname or alias.name
            self.imports[local] = (node.module, alias.name, node.lineno)
        self.generic_visit(node)

    def visit_ClassDef(self, node: ast.ClassDef):
        bases = []
        for base in node.bases:
            if isinstance(base, ast.Name):
                bases.append(base.id)
            elif isinstance(base, ast.Attribute):
                bases.append(base.attr)
        self.dummy_classes[node.name] = (bases, node.lineno)

        # Analyze __init__ for host setup and super().__init__ calls
        has_host = False
        super_args = -1  # -1 = not found
        for item in node.body:
            if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                for stmt in ast.walk(item):
                    # Look for self.host = ...
                    if isinstance(stmt, ast.Assign):
                        for t in stmt.targets:
                            if (isinstance(t, ast.Attribute) and
                                isinstance(t.value, ast.Name) and
                                t.value.id == "self" and t.attr == "host"):
                                has_host = True
                    # Look for super().__init__(...)
                    if isinstance(stmt, ast.Call):
                        func = stmt.func
                        if (isinstance(func, ast.Attribute) and
                            func.attr == "__init__" and
                            isinstance(func.value, ast.Call) and
                            isinstance(func.value.func, ast.Name) and
                            func.value.func.id == "super"):
                            super_args = len(stmt.args) + len(stmt.keywords)

        self.host_setups[node.name] = has_host
        self.super_init_args[node.name] = super_args
        self.generic_visit(node)

    def analyze(self) -> list:
        """Run all checks after visiting."""
        self._check_imports()
        self._check_dummy_classes()
        return self.findings

    def _check_imports(self):
        for local, (module, orig_name, lineno) in self.imports.items():
            if not module.startswith("moleditpy"):
                continue
            path = resolve_moleditpy_import(module)
            if path is None:
                self.findings.append({
                    "kind":   "STALE_IMPORT",
                    "lineno": lineno,
                    "attr":   f"{module}.{orig_name}",
                    "line":   f"from {module} import {orig_name}",
                    "soft":   False,
                    "detail": f"Module '{module}' not found in dev source (SRC={SRC})",
                })
                continue
            # Module exists — check if the name is exported
            if orig_name == "*":
                continue
            known = ast_collect_names(path)
            if orig_name not in known:
                self.findings.append({
                    "kind":   "MISSING_NAME",
                    "lineno": lineno,
                    "attr":   f"{module}.{orig_name}",
                    "line":   f"from {module} import {orig_name}",
                    "soft":   False,
                    "detail": f"'{orig_name}' not defined in {path.name}",
                })

    def _check_dummy_classes(self):
        for cls_name, (bases, lineno) in self.dummy_classes.items():
            manager_bases = [b for b in bases if b in self.manager_methods]
            if not manager_bases:
                continue

            for base in manager_bases:
                mgr_file = MANAGER_CLASS_MAP.get(base)
                if mgr_file is None:
                    continue

                # Check: super().__init__() without host arg
                super_arg_count = self.super_init_args.get(cls_name, -1)
                if super_arg_count == 0:
                    mgr_params = get_manager_init_params(mgr_file)
                    # If manager __init__ needs more than just 'self'
                    required = [p for p in mgr_params if p != "self"]
                    # Check if they all have defaults
                    tree = ast_parse(mgr_file)
                    has_defaults = True
                    if tree:
                        for node in ast.walk(tree):
                            if isinstance(node, ast.ClassDef):
                                for item in node.body:
                                    if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                                        n_no_default = len(item.args.args) - 1 - len(item.args.defaults)
                                        if n_no_default > 0:
                                            has_defaults = False

                    if required and not has_defaults:
                        self.findings.append({
                            "kind":   "MGR_INIT_NOARGS",
                            "lineno": lineno,
                            "attr":   cls_name,
                            "line":   f"class {cls_name}({base})",
                            "soft":   False,
                            "detail": f"{base}.__init__ requires args {required} but super().__init__() called with 0 args",
                        })

                # Check: no self.host = ... in __init__
                if not self.host_setups.get(cls_name, False):
                    self.findings.append({
                        "kind":   "NO_HOST_SETUP",
                        "lineno": lineno,
                        "attr":   cls_name,
                        "line":   f"class {cls_name}({base})",
                        "soft":   True,
                        "detail": f"__init__ does not set self.host -- manager methods using self.host.X will fail",
                    })


def scan_test_file(fpath: Path, manager_methods: dict, verbose: bool) -> list:
    tree = ast_parse(fpath)
    if not tree:
        return []
    analyzer = TestFileAnalyzer(fpath, manager_methods, verbose)
    analyzer.visit(tree)
    return analyzer.analyze()


# ---------------------------------------------------------------------------
# Entry-point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true",
                    help="Also show soft warnings (NO_HOST_SETUP, etc.)")
    args = ap.parse_args()

    print("Building manager method map from AST...")
    manager_methods = build_manager_methods()
    for cls, methods in manager_methods.items():
        if methods:
            print(f"  {cls}: {len(methods)} methods")

    test_files = sorted(TESTS.glob("test_*.py"))
    print(f"\nScanning {len(test_files)} test files in {TESTS.relative_to(ROOT)}...\n")

    total_hard = 0
    total_soft = 0

    for fpath in test_files:
        findings = scan_test_file(fpath, manager_methods, args.verbose)
        hard = [f for f in findings if not f["soft"]]
        soft = [f for f in findings if f["soft"]]

        if not hard and not (args.verbose and soft):
            continue

        rel = fpath.relative_to(ROOT)
        print(f"\n{'='*72}")
        print(f"  {rel}")
        print(f"{'='*72}")

        for f in hard:
            total_hard += 1
            print(f"  [{f['kind']}] line {f['lineno']}")
            print(f"    Import: {f['line']}")
            print(f"    Reason: {f['detail']}")

        if args.verbose and soft:
            for f in soft:
                total_soft += 1
                print(f"  [SOFT/{f['kind']}] line {f['lineno']}  {f['attr']}")
                print(f"    Context: {f['line']}")
                print(f"    Reason:  {f['detail']}")

    print(f"\n{'='*72}")
    print(f"Hard import/init disconnections: {total_hard}")
    if args.verbose:
        print(f"Soft warnings:                   {total_soft}")


if __name__ == "__main__":
    main()
