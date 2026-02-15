
import os
import re

# Target the entire package source to be thorough
TARGET_DIR = r"e:\Research\Calculation\moleditpy\DEV_MAIN\python_molecular_editor\moleditpy\src\moleditpy"

def silence_tracebacks_in_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Pattern 1: import traceback; traceback.print_exc() -> pass
    content = re.sub(r'import\s+traceback\s*;\s*traceback\.print_exc\(\)', 'pass', content)

    # Standard replacement for calls
    content = content.replace('traceback.print_exc()', 'pass')

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

def main():
    print(f"Scanning {TARGET_DIR}...")
    count = 0
    for root, dirs, files in os.walk(TARGET_DIR):
        for file in files:
            if file.endswith(".py"):
                filepath = os.path.join(root, file)
                with open(filepath, 'r', encoding='utf-8') as f:
                    if 'traceback.print_exc' in f.read():
                        silence_tracebacks_in_file(filepath)
                        print(f"Processed: {file}")
                        count += 1
    print(f"Finished. Processed {count} files.")

if __name__ == "__main__":
    main()
