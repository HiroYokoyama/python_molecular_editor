from unittest.mock import patch


def test_imports_mainwindow():
    """Ensure MainWindow and its init submodule can be imported without crashing."""
    from moleditpy.ui.main_window_init import MainInitManager

    assert hasattr(MainInitManager, "init_ui")


def test_mainwindow_init_with_mocks(app):
    """Verify MainWindow instantiates MainInitManager during initialization."""
    # Patch the MainInitManager class in the main_window module
    with patch("moleditpy.ui.main_window.MainInitManager") as MockInitManager:
        from moleditpy.ui.main_window import MainWindow

        # Instantiate MainWindow. Close + delete it before teardown: an un-closed
        # real QMainWindow GC'd after QApplication teardown crashes on Windows
        # when this file runs in isolation.
        mw = MainWindow()
        try:
            # Verify MainInitManager was instantiated
            MockInitManager.assert_called_once()

            # Verify the instance is assigned to mw.init_manager
            assert mw.init_manager == MockInitManager.return_value

            # Verify it was called with mw as the first argument (host)
            args, kwargs = MockInitManager.call_args
            assert args[0] == mw
        finally:
            mw.close()
            mw.deleteLater()


def test_save_settings_logs_unwritable_settings_file(tmp_path, caplog):
    """A settings file that cannot be written is logged, not raised on close."""
    from types import SimpleNamespace

    from moleditpy.ui.main_window_init import MainInitManager

    fake = SimpleNamespace(
        settings={"bond_color": "#000000"},
        settings_dirty=True,
        settings_dir=str(tmp_path),
        settings_file=str(tmp_path / "settings.json"),
        host=SimpleNamespace(initial_settings={}),
    )
    with patch("builtins.open", side_effect=PermissionError("read-only")):
        MainInitManager.save_settings(fake)

    assert fake.settings_dirty is True
    assert "Error saving settings" in caplog.text
