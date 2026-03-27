import sys
import os

# Set base path
base_path = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(base_path, "moleditpy", "src")
sys.path.insert(0, src_path)

print(f"DEBUG: sys.path[0] = {sys.path[0]}")

try:
    print("DEBUG: Importing moleditpy.ui.main_window...")
    from moleditpy.ui.main_window import MainWindow
    print("Import successful!")
except Exception as e:
    import traceback
    traceback.print_exc()
