#!/usr/bin/env python3
"""
harden_hasattr.py  --  Systematically adds logging to hasattr() blocks.

Usage:
    python scripts/harden_hasattr.py          # Add logs
    python scripts/harden_hasattr.py --strip  # Remove added logs (letter removal)
"""

import ast
import argparse
import sys
from pathlib import Path

# Descriptive marker for easy finding and removal
MARKER = "# [REPORT ERROR MISSING ATTRIBUTE]"
OLD_MARKERS = ["# [M]"]

def contains_hasattr(node):
    """Recursively check if an expression is a pure hasattr call or a negated one."""
    if isinstance(node, ast.Call):
        func = node.func
        if isinstance(func, ast.Name) and func.id == "hasattr":
            # extract obj and attr for logging
            obj_name = "object"
            attr_repr = "'unknown'"
            if len(node.args) >= 2:
                arg0 = node.args[0]
                if isinstance(arg0, ast.Name):
                    obj_name = arg0.id
                elif isinstance(arg0, ast.Attribute):
                    if isinstance(arg0.value, ast.Name) and arg0.value.id == "self":
                        obj_name = f"self.{arg0.attr}"
                
                arg1 = node.args[1]
                if isinstance(arg1, (ast.Constant, ast.Str)):
                    val = arg1.value if isinstance(arg1, ast.Constant) else arg1.s
                    attr_repr = f"'{val}'"
                elif isinstance(arg1, ast.Name):
                    attr_repr = f"{{{arg1.id}}}"
            return True, obj_name, attr_repr, False # False = NOT negated
    
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        found, obj, attr, _ = contains_hasattr(node.operand)
        if found:
            return True, obj, attr, True # True = negated
            
    return False, None, None, False

def harden_file(fpath: Path):
    if not fpath.exists():
        return
    
    text = fpath.read_text(encoding="utf-8")
    lines = text.splitlines()
    
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        print(f"Error parsing {fpath}: {e}")
        return

    # To simplify insertion without breaking line numbers/comments,
    # we'll work backwards through the lines.
    
    patches = []

    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            found, obj_name, attr_repr, negated = contains_hasattr(node.test)
            
            # Check for [SAFE] comment to skip
            if_line = lines[node.lineno - 1]
            if "[SAFE]" in if_line:
                continue

            if found and not negated and not node.orelse:
                # Find the last line of the body to determine indentation
                last_node = node.body[-1]
                body_line = lines[last_node.lineno - 1]
                
                if_line = lines[node.lineno - 1]
                base_indent = if_line[:len(if_line) - len(if_line.lstrip())]
                
                # Check for one-liner: if hasattr(...): body
                if last_node.lineno == node.lineno:
                    # One-liner: use base_indent + 4 spaces
                    indent = base_indent + "    "
                else:
                    indent = body_line[:len(body_line) - len(body_line.lstrip())]
                
                end_lineno = getattr(last_node, "end_lineno", last_node.lineno)
                patches.append((end_lineno, base_indent, indent, obj_name, attr_repr))

    if not patches:
        return

    # Check if logging is imported
    has_logging = any(isinstance(node, (ast.Import, ast.ImportFrom)) and 
                      any(alias.name == "logging" for alias in node.names) 
                      for node in ast.walk(tree))
    
    # Sort patches in reverse to avoid shifting lines
    patches.sort(key=lambda x: x[0], reverse=True)
    
    modified = False
    new_lines = list(lines)
    for lineno, b_indent, indent, obj, attr_repr in patches:
        if lineno < len(new_lines) and new_lines[lineno].strip().startswith("else"):
            continue
            
        else_line = f"{b_indent}else:  {MARKER}"
        log_line = f"{indent}logging.error(f\"REPORT ERROR: Missing attribute {attr_repr} on {obj}\")"
        
        new_lines.insert(lineno, else_line)
        new_lines.insert(lineno + 1, log_line)
        modified = True

    if modified:
        if not has_logging:
            # Insert import logging at the top (after software title / docstring)
            insert_pos = 0
            in_docstring = False
            for i, line in enumerate(new_lines):
                if line.startswith("#!") or line.startswith("# -*-"):
                    insert_pos = i + 1
                    continue
                if '"""' in line or "'''" in line:
                    if not in_docstring:
                        in_docstring = True
                        continue
                    else:
                        in_docstring = False
                        insert_pos = i + 1
                        continue # check for __future__ on same line
                if not in_docstring:
                    if line.startswith("from __future__"):
                        insert_pos = i + 1
                        continue
                    if line.strip() != "":
                        break
                    
            new_lines.insert(insert_pos, f"import logging  {MARKER}")
        
        fpath.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
        print(f"Hardened {fpath}")

def strip_markers(fpath: Path):
    if not fpath.exists():
        return
        
    text = fpath.read_text(encoding="utf-8")
    lines = text.splitlines()
    
    new_lines = []
    skip_next = False
    modified = False
    
    all_markers = [MARKER] + OLD_MARKERS
    
    for line in lines:
        if skip_next:
            skip_next = False
            continue
            
        found_marker = False
        for m in all_markers:
            if m in line:
                if "else:" in line:
                    found_marker = True
                    skip_next = True
                    break
                elif "import logging" in line:
                    found_marker = True
                    break
        
        if found_marker:
            modified = True
            continue
            
        new_lines.append(line)
        
    if modified:
        fpath.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
        print(f"Stripped {fpath}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--strip", action="store_true", help="Remove added logging markers (letter removal)")
    args = parser.parse_args()
    
    src_dir = Path("moleditpy/src/moleditpy/ui")
    files = list(src_dir.glob("*.py"))
    
    for f in files:
        if args.strip:
            strip_markers(f)
        else:
            harden_file(f)

if __name__ == "__main__":
    main()
