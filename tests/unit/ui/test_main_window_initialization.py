from unittest.mock import patch

import json
from unittest.mock import MagicMock

import pytest


@pytest.mark.parametrize(
    "key,value",
    [
        ("optimization_method", []),
        ("3d_conversion_mode", {}),
        ("atom_font_family_2d", []),
        ("atom_font_size_2d", "bad"),
        ("atom_font_size_2d", 1.5),
        ("ball_stick_resolution", 0),
        ("bond_dash_count_2d", 0),
        ("specular", float("nan")),
        ("light_intensity", float("inf")),
        ("bond_width_2d", 10**400),
        ("mouse_rotation_sensitivity", -1),
        ("lighting_enabled", []),
        ("cpk_colors", []),
    ],
)
def test_load_settings_replaces_invalid_builtin_values(tmp_path, key, value):
    from moleditpy.ui.main_window_init import MainInitManager
    from moleditpy.utils.default_settings import DEFAULT_SETTINGS

    manager = object.__new__(MainInitManager)
    manager.host = MagicMock()
    manager.settings_file = str(tmp_path / "settings.json")
    payload = {key: value, "plugin_setting": ["retained"], "bond_color_2d": "#123456"}
    (tmp_path / "settings.json").write_text(json.dumps(payload), encoding="utf-8")
    manager.load_settings()
    assert manager.settings[key] == DEFAULT_SETTINGS[key]
    assert manager.settings["plugin_setting"] == ["retained"]
    assert manager.settings["bond_color_2d"] == "#123456"


def test_default_settings_and_color_overrides_have_independent_values():
    from moleditpy.utils.default_settings import DEFAULT_SETTINGS, normalize_settings

    settings = normalize_settings({"cpk_colors": {"C": "#123456", "N": []}})
    assert settings["cpk_colors"] == {"C": "#123456"}
    settings["cpk_colors"]["O"] = "#654321"
    assert DEFAULT_SETTINGS["cpk_colors"] == {}
    assert normalize_settings({"specular_power": 32.0})["specular_power"] == 32


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
