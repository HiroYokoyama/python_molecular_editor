"""Run a plugin package's ``__main__`` from the command line, without the editor.

``moleditpy --run-plugin NAME ARGS...`` behaves like ``python -m NAME ARGS...``
run from the plugin's own folder: a package plugin opts in by shipping a
``__main__.py``. No QApplication or main window is created; the plugin builds
whatever it needs, and its exit status becomes the process's.
"""

import os
import runpy
import sys
from typing import Dict, List, Optional, Sequence, TextIO

from .plugin_manager import PluginManager


def runnable_plugins(manager: Optional[PluginManager] = None) -> Dict[str, List[str]]:
    """Map each runnable plugin name to the package folders that provide it.

    A runnable plugin is an enabled package plugin (a folder with
    ``__init__.py``) that also has ``__main__.py``. Folders are found the way
    the editor discovers plugins, category folders included. More than one
    folder under one name is reported, not resolved.
    """
    manager = manager or PluginManager()
    found: Dict[str, List[str]] = {}
    if not os.path.isdir(manager.plugin_dir):
        return found
    for root, dirs, files in os.walk(manager.plugin_dir):
        dirs[:] = sorted(
            d for d in dirs if not d.startswith("__") and not d.startswith(".")
        )
        if "__init__.py" not in files:
            continue
        dirs[:] = []
        if "__main__.py" not in files:
            continue
        if manager.is_plugin_disabled(os.path.join(root, "__init__.py")):
            continue
        name = os.path.basename(root)
        if name.isidentifier():
            found.setdefault(name, []).append(root)
    return found


def list_plugins(
    manager: Optional[PluginManager] = None, out: Optional[TextIO] = None
) -> int:
    """Print the runnable plugins, one per line, and return an exit status."""
    out = out or sys.stdout
    plugins = runnable_plugins(manager)
    if not plugins:
        print("No plugin can be run from the command line.", file=out)
        return 0
    for name in sorted(plugins):
        print(f"{name}\t{plugins[name][0]}", file=out)
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
    folders = runnable_plugins(manager).get(name, [])
    if not folders:
        print(
            f"Error: no runnable plugin named '{name}'. A plugin can be run when "
            "it is an enabled package with a __main__.py; see --list-plugins.",
            file=sys.stderr,
        )
        return 2
    if len(folders) > 1:
        print(
            f"Error: more than one plugin is named '{name}':\n  "
            + "\n  ".join(folders),
            file=sys.stderr,
        )
        return 2

    parent = os.path.dirname(folders[0])
    saved_argv, saved_path = sys.argv[:], sys.path[:]
    sys.argv = [name, *argv]
    sys.path.insert(0, parent)
    try:
        runpy.run_module(name, run_name="__main__", alter_sys=True)
    except SystemExit as exit_request:
        return _exit_status(exit_request.code)
    finally:
        sys.argv, sys.path[:] = saved_argv, saved_path
    return 0
