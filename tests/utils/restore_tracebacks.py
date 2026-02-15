
import os
import re

TARGET_DIR = r"e:\Research\Calculation\moleditpy\DEV_MAIN\python_molecular_editor\moleditpy\src\moleditpy"

def restore_traceback_pattern(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    original = content
    
    # Pattern to find:
    # except Exception...: # pragma: no cover
    #     import traceback
    #     pass
    
    # We'll use a regex that matches the indentation and handles optional "as e"
    # and optional multiple lines/newline between import and pass.
    
    # Regex explanation:
    # (^\s+except\s+Exception.*:\s*#\s*pragma:\s*no\s*cover\s*\n)  <- Group 1: except line
    # (\s+)                                                       <- Group 2: indentation
    # (import\s+traceback\s*\n)                                   <- Group 3: import line
    # (\s*\n)*                                                   <- Optional blank lines
    # (\s+)                                                       <- Group 4: same indentation as pass
    # pass(?=\n|$)                                                <- The pass keyword
    
    pattern = re.compile(
        r'(^\s+except\s+Exception.*:\s*#\s*pragma:\s*no\s*cover\s*\n)'
        r'(\s+)(import\s+traceback\s*\n)'
        r'(?:\s*\n)*'
        r'\s+pass(?=\n|$)',
        re.MULTILINE
    )
    
    def replace_func(match):
        except_line = match.group(1)
        indent = match.group(2)
        import_line = match.group(3)
        return f"{except_line}{indent}{import_line}{indent}traceback.print_exc()"

    content = pattern.sub(replace_func, content)
    
    # Also handle the cases where "import traceback" is NOT there but we want it
    # (If my previous scripts replaced the whole block with just "pass")
    
    # Pattern: except Exception...: # pragma: no cover
    #     pass
    pattern_no_import = re.compile(
        r'(^\s+except\s+Exception.*:\s*#\s*pragma:\s*no\s*cover\s*\n)'
        r'(\s+)pass(?=\n|$)',
        re.MULTILINE
    )
    
    def replace_no_import(match):
        except_line = match.group(1)
        indent = match.group(2)
        return f"{except_line}{indent}import traceback\n{indent}traceback.print_exc()"

    content = pattern_no_import.sub(replace_no_import, content)

    if content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

def main():
    print(f"Scanning {TARGET_DIR}...")
    count = 0
    for root, dirs, files in os.walk(TARGET_DIR):
        for file in files:
            if file.endswith(".py"):
                filepath = os.path.join(root, file)
                if restore_traceback_pattern(filepath):
                    print(f"Restored: {file}")
                    count += 1
    print(f"Finished. Restored {count} files.")

if __name__ == "__main__":
    main()
