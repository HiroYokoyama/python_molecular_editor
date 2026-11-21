#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: LGPL-3.0-or-later
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

# Run a plugin's command line from the command line, without the editor.
#
# For a package plugin, ``moleditpy --run-plugin NAME ARGS...`` behaves like
# ``python -m NAME ARGS...`` with the plugin's parent folder put on sys.path;
# the package opts in by shipping a ``__main__.py``. For a single-file plugin it
# behaves like ``python NAME.py ARGS...``; the file opts in with a top-level
# ``if __name__ == "__main__":`` block. The working directory is left alone, so
# relative paths in ARGS mean what they meant to the caller. No QApplication or
# main window is created; the plugin builds whatever it needs, and its exit
# status becomes the process's.

import ast
import logging
import os
import runpy
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, TextIO

from .plugin_manager import PluginManager


def _is_main_guard(test: ast.expr) -> bool:
    if not (
        isinstance(test, ast.Compare)
        and len(test.ops) == 1
        and isinstance(test.ops[0], ast.Eq)
    ):
        return False
    sides = [test.left, test.comparators[0]]
    names = [s for s in sides if isinstance(s, ast.Name) and s.id == "__name__"]
    values = [s for s in sides if isinstance(s, ast.Constant) and s.value == "__main__"]
    return len(names) == 1 and len(values) == 1


def has_main_guard(path: str) -> bool:
    """Whether the file has a top-level ``if __name__ == "__main__":`` block.

    The file is parsed, not imported, so no plugin code runs.
    """
    try:
        with open(path, "rb") as fh:
            tree = ast.parse(fh.read(), filename=path)
    except (OSError, SyntaxError, ValueError) as error:
        logging.error("Could not parse plugin file %s: %s", path, error)
        return False
    return any(
        isinstance(node, ast.If) and _is_main_guard(node.test) for node in tree.body
    )


@dataclass(frozen=True)
class InstalledPlugin:
    """One plugin found in the plugin folder, read without importing it."""

    name: str
    path: str
    display_name: str
    version: str
    package: bool
    enabled: bool
    has_main: bool

    @property
    def runnable(self) -> bool:
        """Whether ``--run-plugin`` can run this plugin."""
        return self.enabled and self.has_main and self.name.isidentifier()


def installed_plugins(manager: Optional[PluginManager] = None) -> List[InstalledPlugin]:
    """Return every installed plugin, the way the editor discovers them.

    Package plugins (folders with ``__init__.py``) and single-file plugins are
    both included, in and below category folders, enabled or not. Metadata is
    read by parsing, so no plugin code runs.
    """
    manager = manager or PluginManager()
    found: List[InstalledPlugin] = []
    if not os.path.isdir(manager.plugin_dir):
        return found

    def add(entry: str, name: str, path: str, package: bool, has_main: bool) -> None:
        info = manager.get_plugin_info_safe(entry, fallback_name=name)
        found.append(
            InstalledPlugin(
                name=name,
                path=path,
                display_name=str(info.get("name", name)),
                version=str(info.get("version", "Unknown")),
                package=package,
                enabled=not manager.is_plugin_disabled(entry),
                has_main=has_main,
            )
        )

    for root, dirs, files in os.walk(manager.plugin_dir):
        dirs[:] = sorted(
            d for d in dirs if not d.startswith("__") and not d.startswith(".")
        )
        if "__init__.py" in files:
            dirs[:] = []
            add(
                os.path.join(root, "__init__.py"),
                os.path.basename(root),
                root,
                True,
                "__main__.py" in files,
            )
            continue
        for filename in sorted(files):
            if filename.endswith(".py") and not filename.startswith("__"):
                path = os.path.join(root, filename)
                add(
                    path,
                    os.path.splitext(filename)[0],
                    path,
                    False,
                    has_main_guard(path),
                )
    return found


def runnable_plugins(manager: Optional[PluginManager] = None) -> Dict[str, List[str]]:
    """Map each runnable plugin name to the folders or files that provide it.

    More than one folder under one name is reported, not resolved.
    """
    runnable: Dict[str, List[str]] = {}
    for plugin in installed_plugins(manager):
        if plugin.runnable:
            runnable.setdefault(plugin.name, []).append(plugin.path)
    return runnable


def list_plugins(
    manager: Optional[PluginManager] = None, out: Optional[TextIO] = None
) -> int:
    """Print every installed plugin and whether it can be run; return 0."""
    manager = manager or PluginManager()
    out = out or sys.stdout
    plugins = installed_plugins(manager)
    print(f"Plugin folder: {manager.plugin_dir}", file=out)
    if not plugins:
        print("No plugins are installed.", file=out)
        return 0
    print("(run: can be started with --run-plugin NAME)\n", file=out)
    for plugin in sorted(plugins, key=lambda p: (p.name.lower(), p.path)):
        flags = []
        if plugin.runnable:
            flags.append("run")
        if not plugin.enabled:
            flags.append("disabled")
        label = f"[{', '.join(flags)}]" if flags else ""
        where = os.path.relpath(plugin.path, manager.plugin_dir)
        print(
            f"  {plugin.name:<28} {label:<16} {plugin.display_name} "
            f"{plugin.version}  ({where})",
            file=out,
        )
    return 0


def _exit_status(code: object) -> int:
    if code is None:
        return 0
    if isinstance(code, int):
        return code
    print(code, file=sys.stderr)
    return 1


def run_plugin(
    name: str, argv: Sequence[str], manager: Optional[PluginManager] = None
) -> int:
    """Run plugin *name*'s ``__main__`` with *argv* and return its exit status."""
    paths = runnable_plugins(manager).get(name, [])
    if not paths:
        print(
            f"Error: no runnable plugin named '{name}'. A plugin can be run when "
            "it is enabled and is a package with a __main__.py or a single file "
            'with an if __name__ == "__main__": block; see --list-plugins.',
            file=sys.stderr,
        )
        return 2
    if len(paths) > 1:
        print(
            f"Error: more than one plugin is named '{name}':\n  " + "\n  ".join(paths),
            file=sys.stderr,
        )
        return 2

    path = paths[0]
    single_file = os.path.isfile(path)
    saved_argv, saved_path = sys.argv[:], sys.path[:]
    # As with ``python FILE``, argv[0] is the file and its folder is importable.
    sys.argv = [path if single_file else name, *argv]
    sys.path.insert(0, os.path.dirname(path))
    try:
        if single_file:
            runpy.run_path(path, run_name="__main__")
        else:
            runpy.run_module(name, run_name="__main__", alter_sys=True)
    except SystemExit as exit_request:
        return _exit_status(exit_request.code)
    finally:
        sys.argv, sys.path[:] = saved_argv, saved_path
    return 0
