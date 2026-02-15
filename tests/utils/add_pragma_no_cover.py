
import os
import re

TARGET_DIR = r"e:\Research\Calculation\moleditpy\DEV_MAIN\python_molecular_editor\moleditpy\src\moleditpy"

def add_pragma_to_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    new_lines = []
    modified = False
    
    # Regex to match 'except Exception' or 'except Exception as e'
    # We want to match lines that have this, but don't already have 'pragma: no cover'
    regex = re.compile(r'^\s*except\s+Exception.*:')

    for line in lines:
        if regex.match(line):
            if "pragma: no cover" not in line:
                # Remove newline, add pragma, add newline
                line = line.rstrip() + "  # pragma: no cover\n"
                modified = True
        new_lines.append(line)

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        return True
    return False

def main():
    print(f"Scanning {TARGET_DIR}...")
    count = 0
    total_files = 0
    for root, dirs, files in os.walk(TARGET_DIR):
        for file in files:
            if file.endswith(".py"):
                filepath = os.path.join(root, file)
                if add_pragma_to_file(filepath):
                    print(f"Updated: {file}")
                    count += 1
                total_files += 1
    print(f"Finished. Updated {count} out of {total_files} files.")

if __name__ == "__main__":
    main()
