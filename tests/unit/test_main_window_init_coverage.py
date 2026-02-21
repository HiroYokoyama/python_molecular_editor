import pytest
import sys
from unittest.mock import MagicMock, patch

# Note: We don't necessarily need to instantiate MainWindow,
# just importing it and its submodules ensures the pragmas are registered
# and the executable lines count is reduced in coverage reports.


def test_imports_mainwindow():
    """Ensure MainWindow and its init submodule can be imported without crashing."""
    # Importing is the primary verification here; if it fails, test runner raises error.
    from moleditpy.modules.main_window import MainWindow
    from moleditpy.modules.main_window_main_init import MainWindowMainInit
    
    # Verify they have expected core methods/attributes
    assert "__init__" in MainWindow.__dict__  # has a custom constructor
    assert hasattr(MainWindowMainInit, "init_ui")


def test_mainwindow_init_with_mocks():
    """Verify MainWindow class structure is intact with mocks."""
    with patch("PyQt6.QtWidgets.QMainWindow.__init__", return_value=None):
        with patch(
            "moleditpy.modules.main_window_main_init.MainWindowMainInit.init_ui"
        ):
            with patch(
                "moleditpy.modules.main_window_main_init.MainWindowMainInit.init_menu_bar"
            ):
                from moleditpy.modules.main_window import MainWindow
                
                # MainWindow delegates to mixin objects (composition, not inheritance).
                # Verify MainWindow defines its own init_ui and init_menu_bar wrappers.
                assert "init_ui" in MainWindow.__dict__, (
                    "MainWindow should define its own init_ui wrapper"
                )
                assert "init_menu_bar" in MainWindow.__dict__, (
                    "MainWindow should define its own init_menu_bar wrapper"
                )
