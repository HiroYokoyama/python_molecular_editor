import ast
import os
import hashlib
from collections import defaultdict

def normalize_ast(node, fuzzy=False):
    """Normalize AST for comparison."""
    # Create a copy of the node to avoid modifying the original tree (though we are using ast.parse again later)
    # Actually ast.walk visits all children.
    # To normalize, we strip location info.
    for n in ast.walk(node):
        for attr in ('lineno', 'col_offset', 'end_lineno', 'end_col_offset'):
            if hasattr(n, attr):
                setattr(n, attr, 0)
        
        # For fuzzy matching, neutralize names but keep structure
        if fuzzy:
            if isinstance(n, ast.Name):
                n.id = 'var'
            elif isinstance(n, ast.arg):
                n.arg = 'arg'
            elif isinstance(n, ast.Attribute):
                n.attr = 'attr'
            elif isinstance(n, ast.keyword):
                n.arg = 'kwarg'
            # Constants can also be neutralized if very fuzzy
            # elif isinstance(n, ast.Constant):
            #     n.value = None
    
    return ast.dump(node)

def find_duplicates(directory, fuzzy=False, min_lines=2):
    """Find duplicated functions/methods by AST structure."""
    duplicates = defaultdict(list)
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        source = f.read()
                        tree = ast.parse(source)
                    
                    for node in ast.walk(tree):
                        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                            # Only consider functions with a minimum "complexity" or line count
                            # Simple heuristic for line count: 
                            if hasattr(node, 'end_lineno') and hasattr(node, 'lineno'):
                                lines = node.end_lineno - node.lineno
                            else:
                                lines = len(node.body)
                                
                            if lines < min_lines:
                                continue
                                
                            norm = normalize_ast(node, fuzzy=fuzzy)
                            h = hashlib.md5(norm.encode()).hexdigest()
                            duplicates[h].append({
                                'file': path,
                                'name': node.name,
                                'line': getattr(node, 'lineno', 0),
                                'lines': lines
                            })
                except Exception as e:
                    print(f"Error parsing {path}: {e}")
                    
    return duplicates

def report_duplicates(duplicates):
    """Print a summary of duplicates found."""
    found_any = False
    sorted_dups = sorted(duplicates.items(), key=lambda x: len(x[1]), reverse=True)
    
    for h, instances in sorted_dups:
        if len(instances) > 1:
            found_any = True
            # Check if all instances are the same function (e.g. property setters/getters might look similar)
            # Actually, the name check is not enough because we want to find refactored logic.
            
            print(f"[{len(instances)} Matches] AST Hash: {h[:8]}...")
            for inst in instances:
                print(f"  - {inst['file']}:{inst['line']} -> {inst['name']} ({inst['lines']} lines)")
            print("-" * 40)
            
    if not found_any:
        print("No duplicates found with the current settings.")

if __name__ == "__main__":
    import sys
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    fuzzy_mode = "--fuzzy" in sys.argv
    
    print(f"Scanning directory: {target_dir}")
    print(f"Fuzzy mode: {'ENABLED' if fuzzy_mode else 'DISABLED'}")
    print("=" * 60)
    
    dups = find_duplicates(target_dir, fuzzy=fuzzy_mode)
    report_duplicates(dups)
