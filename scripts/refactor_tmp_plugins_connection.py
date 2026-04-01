#!/usr/bin/env python3
"""Refactor tmp/plugins for robust context/main-window connection handling.

What this codemod does:
1) Normalize top-level `run(...)` entrypoints so they accept either
   - PluginContext
   - MainWindow
   - host wrapper (legacy)
2) Fix menu registration callbacks in initialize(context):
   bare `run` callback -> `lambda: run(context)`

Scope: tmp/plugins/**/*.py
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

ROOT = Path("tmp/plugins")
MENU_METHODS = {"add_menu_action", "add_plugin_menu", "register_menu_action"}


@dataclass
class Edit:
    start: int
    end: int
    replacement: str


def build_line_starts(text: str) -> List[int]:
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)
    return starts


def off(starts: List[int], line: int, col: int) -> int:
    return starts[line - 1] + col


def callback_nodes_for_call(call: ast.Call) -> Iterable[ast.AST]:
    if not isinstance(call.func, ast.Attribute):
        return []
    m = call.func.attr
    if m not in MENU_METHODS:
        return []

    if m in {"add_menu_action", "add_plugin_menu"}:
        return [call.args[1]] if len(call.args) >= 2 else []

    nodes: List[ast.AST] = []
    if len(call.args) >= 2:
        nodes.append(call.args[1])
    if len(call.args) >= 3:
        nodes.append(call.args[2])
    return nodes


def has_run_resolver(func_node: ast.FunctionDef, source_segment: str) -> bool:
    needles = ["get_main_window", ".host", "PluginContext"]
    return any(n in source_segment for n in needles)


def collect_edits(text: str) -> List[Edit]:
    tree = ast.parse(text)
    starts = build_line_starts(text)
    edits: List[Edit] = []

    # 1) Menu callback rewrites inside initialize(context)
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or node.name != "initialize":
            continue
        if not node.args.args:
            continue

        ctx_name = node.args.args[0].arg

        for sub in ast.walk(node):
            if not isinstance(sub, ast.Call):
                continue
            if not isinstance(sub.func, ast.Attribute):
                continue
            owner = sub.func.value
            if not isinstance(owner, ast.Name) or owner.id != ctx_name:
                continue

            for cb in callback_nodes_for_call(sub):
                if isinstance(cb, ast.Name) and cb.id == "run":
                    s = off(starts, cb.lineno, cb.col_offset)
                    e = off(starts, cb.end_lineno, cb.end_col_offset)
                    edits.append(Edit(s, e, f"lambda: run({ctx_name})"))

    # 2) Add run resolver to top-level run(...) when needed
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or node.name != "run":
            continue
        if len(node.args.args) != 1:
            continue

        arg_name = node.args.args[0].arg

        seg_start = off(starts, node.lineno, node.col_offset)
        seg_end = off(starts, node.end_lineno, node.end_col_offset)
        fn_seg = text[seg_start:seg_end]
        if has_run_resolver(node, fn_seg):
            continue

        insert_before_line = node.body[0].lineno
        if (
            isinstance(node.body[0], ast.Expr)
            and isinstance(node.body[0].value, ast.Constant)
            and isinstance(node.body[0].value.value, str)
        ):
            if len(node.body) > 1:
                insert_before_line = node.body[1].lineno
            else:
                # def + docstring only: insert before function end
                insert_before_line = node.end_lineno

        insert_pos = off(starts, insert_before_line, 0)
        indent = " " * (node.col_offset + 4)
        snippet = (
            f"{indent}# Accept either PluginContext or MainWindow entrypoint\n"
            f"{indent}if hasattr({arg_name}, 'get_main_window'):\n"
            f"{indent}    {arg_name} = {arg_name}.get_main_window()\n"
            f"{indent}if hasattr({arg_name}, 'host'):\n"
            f"{indent}    {arg_name} = {arg_name}.host\n\n"
        )
        edits.append(Edit(insert_pos, insert_pos, snippet))

    return edits


def apply_edits(text: str, edits: List[Edit]) -> str:
    out = text
    for e in sorted(edits, key=lambda x: x.start, reverse=True):
        out = out[: e.start] + e.replacement + out[e.end :]
    return out


def process(path: Path) -> Tuple[bool, int]:
    try:
        src = path.read_text(encoding="utf-8")
        encoding = "utf-8"
    except UnicodeDecodeError:
        src = path.read_text(encoding="cp932")
        encoding = "cp932"

    edits = collect_edits(src)
    if not edits:
        return False, 0

    dst = apply_edits(src, edits)
    if dst == src:
        return False, 0

    # Write back as utf-8 to keep code path uniform in repo
    path.write_text(dst, encoding="utf-8")
    return True, len(edits)


def main() -> int:
    changed = 0
    edit_count = 0
    for py in ROOT.rglob("*.py"):
        ok, n = process(py)
        if ok:
            changed += 1
            edit_count += n
            print(f"updated: {py} ({n} edit(s))")

    print(f"done: files={changed}, edits={edit_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
