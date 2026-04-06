from unittest.mock import patch

def test_imports_mainwindow():
    """Ensure MainWindow and its init submodule can be imported without crashing."""
    from moleditpy.ui.main_window_init import MainInitManager
    assert hasattr(MainInitManager, "init_ui")

def test_mainwindow_init_with_mocks():
    """Verify MainWindow instantiates MainInitManager during initialization."""
    # Patch the MainInitManager class in the main_window module
    with patch("moleditpy.ui.main_window.MainInitManager") as MockInitManager:
        from moleditpy.ui.main_window import MainWindow

        # Instantiate MainWindow
        mw = MainWindow()

        # Verify MainInitManager was instantiated
        MockInitManager.assert_called_once()

        # Verify the instance is assigned to mw.init_manager
        assert mw.init_manager == MockInitManager.return_value

        # Verify it was called with mw as the first argument (host)
        args, kwargs = MockInitManager.call_args
        assert args[0] == mw
