import pytest
import os
import sys
from PyQt6.QtWidgets import QApplication

# Path setup
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'moleditpy', 'src'))
if os.path.isdir(src_path) and src_path not in sys.path:
    sys.path.insert(0, src_path)

@pytest.fixture(scope="session")
def app():
    """QApplication session-wide instance for integration tests."""
    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)
    return q_app
