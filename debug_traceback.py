import sys
import os
from PyQt6.QtWidgets import QApplication

# Set base path
base_path = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(base_path, "moleditpy", "src")
sys.path.insert(0, src_path)

# Mock QApplication
app = QApplication(sys.argv)

try:
    from moleditpy.ui.main_window import MainWindow
    print("DEBUG: Instantiating MainWindow...")
    window = MainWindow()
    
    print("DEBUG: Creating an atom...")
    from PyQt6.QtCore import QPointF
    window.scene.create_atom("C", QPointF(0, 0))
    window.push_undo_state()
    
    print("DEBUG: Calling undo...")
    window.undo()
    
    print("DEBUG: Calling redo...")
    window.redo()
    
    print("Success!")
except Exception as e:
    import traceback
    traceback.print_exc()
