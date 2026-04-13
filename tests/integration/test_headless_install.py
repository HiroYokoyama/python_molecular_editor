import os
import sys
import textwrap
from unittest.mock import patch, MagicMock
import pytest
from moleditpy.main import main
from moleditpy.plugins.plugin_manager import PluginManager


def test_headless_install_success(tmp_path, capsys):
    """Test successful headless plugin installation."""
    # 1. Setup a dummy plugin file
    plugin_src = tmp_path / "test_cli_plugin.py"
    plugin_src.write_text(
        textwrap.dedent("""
        PLUGIN_NAME = "CLI Test Plugin"
        PLUGIN_VERSION = "0.0.1"
        PLUGIN_AUTHOR = "Test Bot"
        PLUGIN_DESCRIPTION = "Integration test for CLI"
    """)
    )

    # 2. Setup a dummy plugin directory
    test_plugin_dir = tmp_path / "moledit_plugins"
    test_plugin_dir.mkdir()

    # 3. Patch PluginManager class in the source module
    with (
        patch("moleditpy.plugins.plugin_manager.PluginManager") as MockPM,
        patch.object(sys, "argv", ["moleditpy", "--install-plugin", str(plugin_src)]),
        patch("builtins.input", return_value="y"),
    ):
        # Setup the mock instance
        mock_pm = MockPM.return_value
        mock_pm.plugin_dir = str(test_plugin_dir)
        mock_pm.get_plugin_info_safe.return_value = {
            "name": "CLI Test Plugin",
            "author": "Test Bot",
            "version": "0.0.1",
            "description": "Integration test for CLI",
        }
        mock_pm._compute_sha256.return_value = "fake-sha256"
        mock_pm.install_plugin.return_value = (True, "Success")

        # Run main() which handles the install-plugin arg
        with pytest.raises(SystemExit) as cm:
            main()
        assert cm.value.code == 0

    # 5. Verify stdout
    captured = capsys.readouterr()
    assert "PLUGIN INSTALLATION (HEADLESS)" in captured.out
    assert "CLI Test Plugin" in captured.out
    assert "Success:" in captured.out


def test_headless_install_abort(tmp_path, capsys):
    """Test aborted headless plugin installation (answering 'n')."""
    plugin_src = tmp_path / "abort_plugin.py"
    plugin_src.write_text("PLUGIN_NAME = 'Abort Plugin'")

    # Patch PluginManager class in the source module
    with (
        patch("moleditpy.plugins.plugin_manager.PluginManager") as MockPM,
        patch.object(sys, "argv", ["moleditpy", "--install-plugin", str(plugin_src)]),
        patch("builtins.input", return_value="n"),
    ):
        mock_pm = MockPM.return_value
        mock_pm.plugin_dir = str(tmp_path)
        mock_pm.get_plugin_info_safe.return_value = {"name": "Abort"}
        mock_pm._compute_sha256.return_value = "fake-sha256"

        with pytest.raises(SystemExit) as cm:
            main()
        assert cm.value.code == 0

    # Verify install_plugin was NOT called
    mock_pm.install_plugin.assert_not_called()

    captured = capsys.readouterr()
    assert "Installation aborted." in captured.out


def test_headless_install_invalid_path(tmp_path, capsys):
    """Test error message when plugin path is invalid."""
    with patch.object(
        sys, "argv", ["moleditpy", "--install-plugin", "non_existent_file.py"]
    ):
        with pytest.raises(SystemExit) as cm:
            main()
        assert cm.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Plugin path not found" in captured.out
