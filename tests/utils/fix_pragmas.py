
import os
import re

TARGET_DIR = r"e:\Research\Calculation\moleditpy\DEV_MAIN\python_molecular_editor\moleditpy\src\moleditpy"

def fix_pragmas_robust(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    new_lines = []
    modified = False
    i = 0
    
    pattern_except = re.compile(r'^(\s*)except(\s+Exception.*)?\s*:')
    pattern_try = re.compile(r'^(\s*)try:\s*(#\s*pragma:\s*no\s*cover)?')

    # Check for global traceback import
    global_has_traceback = any(re.match(r'^import\s+traceback', l) for l in lines)

    while i < len(lines):
        line = lines[i]
        match_try = pattern_try.match(line)
        match_except = pattern_except.match(line)
        
        target_match = match_try or match_except
        
        if target_match:
            indent = target_match.group(1)
            j = i + 1
            body_lines = []
            while j < len(lines):
                if not lines[j].strip():
                    body_lines.append(lines[j])
                    j += 1
                    continue
                next_indent = len(lines[j]) - len(lines[j].lstrip())
                if next_indent > len(indent):
                    body_lines.append(lines[j])
                    j += 1
                else: break
            
            content = [l.strip() for l in body_lines if l.strip() and not l.strip().startswith('#')]
            
            if match_try:
                if content and "# pragma: no cover" in line:
                    line = line.split('#')[0].rstrip() + "\n"
                    modified = True
            
            if match_except:
                exempt_statements = {"pass", "import traceback", "traceback.print_exc()"}
                is_exempt = len(content) > 0 and all(stmt in exempt_statements for stmt in content)
                is_traceback_block = "traceback.print_exc()" in content
                
                # Check for missing import traceback in the block OR globally
                if is_traceback_block and "import traceback" not in content and not global_has_traceback:
                    # Insert import traceback at the top of the block
                    body_lines_temp = [l for l in body_lines if l.strip()]
                    # Find first line of body
                    first_line_idx = -1
                    for idx, bl in enumerate(body_lines):
                        if bl.strip():
                            first_line_idx = idx
                            break
                    if first_line_idx != -1:
                        # Add import traceback before the content
                        body_lines.insert(first_line_idx, indent + "    import traceback\n")
                        modified = True
                        # Re-calculate content
                        content = [l.strip() for l in body_lines if l.strip() and not l.strip().startswith('#')]
                
                has_pragma = "# pragma: no cover" in line
                if is_exempt and not has_pragma:
                    line = line.split('#')[0].rstrip() + "  # pragma: no cover\n"
                    modified = True
                elif not is_exempt and has_pragma:
                    line = line.split('#')[0].rstrip() + "\n"
                    modified = True

                if is_traceback_block:
                    # Clean up internal blank lines
                    temp_body = [l for l in body_lines if l.strip()]
                    if temp_body != body_lines:
                        body_lines = temp_body
                        modified = True
                    
                    if j < len(lines):
                        next_line = lines[j].strip()
                        if next_line and not next_line.startswith(('except', 'else', 'finally')):
                            if not (body_lines and not body_lines[-1].strip()):
                                body_lines.append("\n")
                                modified = True
            
            new_lines.append(line)
            new_lines.extend(body_lines)
            i = j
            continue

        new_lines.append(line)
        i += 1

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        return True
    return False

def run_main():
    print(f"Scanning {TARGET_DIR}...")
    count = 0
    for root, dirs, files in os.walk(TARGET_DIR):
        for file in files:
            if file.endswith(".py"):
                filepath = os.path.join(root, file)
                if fix_pragmas_robust(filepath):
                    print(f"Modified: {file}")
                    count += 1
    print(f"Finished. Modified {count} files.")

if __name__ == "__main__":
    run_main()
