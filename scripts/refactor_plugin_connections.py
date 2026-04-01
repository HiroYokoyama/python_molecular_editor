#!/usr/bin/env python3
"""Codemod: fix plugin callback connection wiring for tmp/plugins.

This rewrites menu registrations inside initialize(context) from:
    context.add_menu_action(..., run)
into:
    context.add_menu_action(..., lambda: run(context.get_main_window()))

It also handles add_plugin_menu/register_menu_action when they pass bare `run`.
Only the callback expression token is replaced, so formatting stays mostly intact.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

TARGET_METHODS = {"add_menu_action", "add_plugin_menu", "register_menu_action"}


@dataclass
class Edit:
    start: int
    end: int
    replacement: str


def line_starts(text: str) -> List[int]:
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)
    return starts


def to_offset(starts: List[int], line: int, col: int) -> int:
    return starts[line - 1] + col


def callback_nodes_for_call(call: ast.Call) -> Iterable[ast.AST]:
    func = call.func
    if not isinstance(func, ast.Attribute):
        return []

    method = func.attr
    if method not in TARGET_METHODS:
        return []

    if method in {"add_menu_action", "add_plugin_menu"}:
        if len(call.args) >= 2:
            return [call.args[1]]
        return []

    # register_menu_action supports both:
    # (path, callback, ...)
    # (path, text, callback, ...)
    nodes: List[ast.AST] = []
    if len(call.args) >= 2:
        nodes.append(call.args[1])
    if len(call.args) >= 3:
        nodes.append(call.args[2])
    return nodes


def collect_edits(text: str) -> List[Edit]:
    tree = ast.parse(text)
    starts = line_starts(text)
    edits: List[Edit] = []

    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        if node.name != "initialize":
            continue
        if not node.args.args:
            continue

        context_name = node.args.args[0].arg

        for sub in ast.walk(node):
            if not isinstance(sub, ast.Call):
                continue
            if not isinstance(sub.func, ast.Attribute):
                continue
            owner = sub.func.value
            if not isinstance(owner, ast.Name):
                continue
            if owner.id != context_name:
                continue

            for cb in callback_nodes_for_call(sub):
                if isinstance(cb, ast.Name) and cb.id == "run":
                    s = to_offset(starts, cb.lineno, cb.col_offset)
                    e = to_offset(starts, cb.end_lineno, cb.end_col_offset)
                    repl = f"lambda: run({context_name}.get_main_window())"
                    edits.append(Edit(start=s, end=e, replacement=repl))

    return edits


def apply_edits(text: str, edits: List[Edit]) -> str:
    if not edits:
        return text
    out = text
    for ed in sorted(edits, key=lambda x: x.start, reverse=True):
        out = out[: ed.start] + ed.replacement + out[ed.end :]
    return out


def process_file(path: Path) -> Tuple[bool, int]:
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="cp932")

    edits = collect_edits(text)
    if not edits:
        return False, 0

    new_text = apply_edits(text, edits)
    if new_text != text:
        path.write_text(new_text, encoding="utf-8")
        return True, len(edits)

    return False, 0


def main() -> int:
    root = Path("tmp/plugins")
    changed_files = 0
    changed_refs = 0

    for py in root.rglob("*.py"):
        changed, edits = process_file(py)
        if changed:
            changed_files += 1
            changed_refs += edits
            print(f"updated: {py} ({edits} replacement(s))")

    print(f"done: files={changed_files}, replacements={changed_refs}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
