import ast
import os
import sys

def extract_assertions(filepath):
    """Parse a test file and extract function names, docstrings, and assertions."""
    with open(filepath, "r", encoding="utf-8") as f:
        try:
            tree = ast.parse(f.read())
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
            return []

    results = []
    
    # helper to process function nodes
    def process_function(func_node, parent_name=None):
        if not func_node.name.startswith("test_"):
            return

        display_name = f"{parent_name}.{func_node.name}" if parent_name else func_node.name
        docstring = ast.get_docstring(func_node) or "No description provided."
        doc_lines = docstring.strip().split("\n")
        summary = doc_lines[0].strip() if doc_lines else "No description."
        
        assertions = []
        for subnode in ast.walk(func_node):
            if isinstance(subnode, ast.Assert):
                try:
                    if hasattr(ast, "unparse"):
                        assertions.append(f"assert {ast.unparse(subnode.test)}")
                    else:
                        assertions.append("assert [expression]")
                except:
                    assertions.append("assert [complex expression]")
        
        results.append({
            "name": display_name,
            "description": summary,
            "assertions": assertions
        })

    for node in tree.body:
        # Top-level functions
        if isinstance(node, ast.FunctionDef):
            process_function(node)
        # Class methods
        elif isinstance(node, ast.ClassDef):
            for subnode in node.body:
                if isinstance(subnode, ast.FunctionDef):
                    process_function(subnode, parent_name=node.name)

    return results

def generate_catalog():
    # Adjusted for tests/utils/ location
    utils_dir = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.dirname(utils_dir)
    root_dir = os.path.dirname(tests_dir)
    output_file = os.path.join(tests_dir, "assertion_catalog.md")
    
    subdirs = ["unit", "integration", "gui"]
    
    markdown_lines = ["# Test Assertions Catalog", ""]
    
    for subdir in subdirs:
        target_path = os.path.join(tests_dir, subdir)
        if not os.path.exists(target_path):
            continue
            
        test_files = [f for f in os.listdir(target_path) if f.startswith("test_") and f.endswith(".py")]
        
        for test_file in sorted(test_files):
            rel_path = f"tests/{subdir}/{test_file}"
            markdown_lines.append(f"## {rel_path.replace('\\', '/')}")
            markdown_lines.append("")
            
            file_results = extract_assertions(os.path.join(target_path, test_file))
            for res in file_results:
                markdown_lines.append(f"### {res['name']}")
                markdown_lines.append(f"_{res['description']}_")
                markdown_lines.append("")
                for ass in res['assertions']:
                    markdown_lines.append(f"- {ass}")
                markdown_lines.append("")
                
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(markdown_lines))
    
    print(f"Assertion catalog updated: {output_file}")

if __name__ == "__main__":
    generate_catalog()
