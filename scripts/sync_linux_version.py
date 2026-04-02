import os
import shutil
import re

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_MAIN = os.path.join(BASE_DIR, "moleditpy", "src", "moleditpy")
SRC_LINUX = os.path.join(BASE_DIR, "moleditpy-linux", "src", "moleditpy_linux")

def sync_linux():
    print(f"Synchronizing Linux version from {SRC_MAIN} to {SRC_LINUX}...")

    if not os.path.exists(SRC_MAIN):
        print(f"Error: Main source directory not found at {SRC_MAIN}")
        return

    # 1. Clean destination (preserve __pycache__ if desired, but here we clean all)
    if os.path.exists(SRC_LINUX):
        # We preserve the directory itself but clear its contents
        for item in os.listdir(SRC_LINUX):
            item_path = os.path.join(SRC_LINUX, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
            else:
                os.remove(item_path)
    else:
        os.makedirs(SRC_LINUX, exist_ok=True)

    # 2. Copy all files from main to linux
    for root, dirs, files in os.walk(SRC_MAIN):
        # Skip pycache
        if "__pycache__" in dirs:
            dirs.remove("__pycache__")
        
        # Calculate relative path
        rel_path = os.path.relpath(root, SRC_MAIN)
        dest_dir = os.path.join(SRC_LINUX, rel_path)
        os.makedirs(dest_dir, exist_ok=True)

        for file in files:
            if file.endswith(".pyc") or file.endswith(".pyo") or "__pycache__" in root:
                continue
            
            src_file = os.path.join(root, file)
            dest_file = os.path.join(dest_dir, file)
            
            # For Python files, perform transformation
            if file.endswith(".py"):
                try:
                    with open(src_file, "r", encoding="utf-8") as f:
                        content = f.read()
                    
                    # Transformation 1: Rename package
                    content = content.replace("moleditpy.", "moleditpy_linux.")
                    content = content.replace("from moleditpy import", "from moleditpy_linux import")
                    content = content.replace("import moleditpy", "import moleditpy_linux")
                    
                    # Transformation 2: Explicit Open Babel disablement
                    if file == "__init__.py" and rel_path == ".":
                        content = re.sub(
                            r"OBABEL_AVAILABLE = (True|False|importlib.*None)",
                            "OBABEL_AVAILABLE = False",
                            content,
                            flags=re.DOTALL
                        )

                    # Transformation 3: Fix double mentions
                    content = content.replace("moleditpy_linux_linux", "moleditpy_linux")

                    with open(dest_file, "w", encoding="utf-8") as f:
                        f.write(content)
                except UnicodeDecodeError:
                    print(f"Warning: Could not decode {src_file} as UTF-8. Copying as binary.")
                    shutil.copy2(src_file, dest_file)
            else:
                # Binary files (icons, assets, etc.)
                shutil.copy2(src_file, dest_file)

    print("Synchronization complete!")

if __name__ == "__main__":
    sync_linux()
