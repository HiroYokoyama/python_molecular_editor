
import os
import re

TARGET_DIR = r"e:\Research\Calculation\moleditpy\DEV_MAIN\python_molecular_editor\moleditpy\src\moleditpy"

def fix_pragmas_in_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    modified = False
    temp_lines = list(lines)
    
    # Pattern for try or except lines
    # We now target any 'except Exception' line to ensure it has or lacks pragma correctly.
    # We also target 'try' lines to REMOVE pragmas if they have them, as try blocks usually contain logic.
    pattern_except = re.compile(r'^(\s*)except\s+Exception.*:')
    pattern_try = re.compile(r'^(\s*)try:\s*#\s*pragma:\s*no\s*cover')

    for i in range(len(temp_lines)):
        line = temp_lines[i]
        
        # 1. Handle 'try' lines (Remove pragma if body has logic, which is almost always)
        match_try = pattern_try.match(line)
        if match_try:
            indent = match_try.group(1)
            body_lines = []
            j = i + 1
            while j < len(temp_lines):
                next_line = temp_lines[j]
                if not next_line.strip():
                    j += 1
                    continue
                next_indent = len(next_line) - len(next_line.lstrip())
                if next_indent > len(indent):
                    body_lines.append(next_line)
                    j += 1
                else:
                    break
            
            content_lines = [l.strip() for l in body_lines if l.strip() and not l.strip().startswith('#')]
            if content_lines: # If try block has ANY code, remove pragma
                new_line = line.split('#')[0].rstrip() + "\n"
                if new_line != line:
                    temp_lines[i] = new_line
                    modified = True

        # 2. Handle 'except Exception' lines
        match_except = pattern_except.match(line)
        if match_except:
            indent = match_except.group(1)
            body_lines = []
            j = i + 1
            while j < len(temp_lines):
                next_line = temp_lines[j]
                if not next_line.strip():
                    j += 1
                    continue
                next_indent = len(next_line) - len(next_line.lstrip())
                if next_indent > len(indent):
                    body_lines.append(next_line)
                    j += 1
                else:
                    break
            
            content_lines = [l.strip() for l in body_lines if l.strip() and not l.strip().startswith('#')]
            
            # Logic: Should have pragma if block ONLY contains pass/import traceback/traceback.print_exc()
            exempt_statements = {"pass", "import traceback", "traceback.print_exc()"}
            is_exempt = len(content_lines) > 0 and all(stmt in exempt_statements for stmt in content_lines)
            
            has_pragma = "# pragma: no cover" in line
            
            if is_exempt and not has_pragma:
                # ADD pragma
                new_line = line.split('#')[0].rstrip() + "  # pragma: no cover\n"
                if new_line != line:
                    temp_lines[i] = new_line
                    modified = True
            elif not is_exempt and has_pragma:
                # REMOVE pragma
                new_line = line.split('#')[0].rstrip() + "\n"
                if new_line != line:
                    temp_lines[i] = new_line
                    modified = True

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(temp_lines)
        return True
    return False

def main():
    print(f"Scanning {TARGET_DIR}...")
    count = 0
    for root, dirs, files in os.walk(TARGET_DIR):
        for file in files:
            if file.endswith(".py"):
                filepath = os.path.join(root, file)
                if fix_pragmas_in_file(filepath):
                    print(f"Fixed: {file}")
                    count += 1
    print(f"Finished. Fixed {count} files.")

if __name__ == "__main__":
    main()
