"""The icon painter must ask the OS for its theme once, not once per icon."""

from unittest.mock import MagicMock, patch

from PyQt6.QtGui import QColor

from moleditpy.ui.main_window_init import _UNQUERIED, MainInitManager


def _manager(settings=None):
    manager = MagicMock(spec=MainInitManager)
    manager.settings = settings if settings is not None else {}
    manager._os_dark_pref = _UNQUERIED
    return manager


class TestOsThemeIsQueriedOnce:
    def test_repeated_icons_fork_only_one_helper(self, app):
        """detect_system_dark_mode forks `defaults`/`gsettings`; 15+ per launch hung CI."""
        manager = _manager()
        with patch(
            "moleditpy.ui.main_window_init.detect_system_dark_mode", return_value=True
        ) as detect:
            for _ in range(15):
                MainInitManager._get_icon_foreground_color(manager)
        detect.assert_called_once()

    def test_unknown_preference_is_not_re_queried(self, app):
        """None is a real answer ("OS could not say"), so it must also be cached."""
        manager = _manager()
        with patch(
            "moleditpy.ui.main_window_init.detect_system_dark_mode", return_value=None
        ) as detect:
            for _ in range(5):
                MainInitManager._get_icon_foreground_color(manager)
        detect.assert_called_once()


class TestColourIsStillCorrect:
    def test_dark_os_gives_white_icons(self, app):
        manager = _manager()
        with patch(
            "moleditpy.ui.main_window_init.detect_system_dark_mode", return_value=True
        ):
            assert MainInitManager._get_icon_foreground_color(manager) == QColor(
                "#FFFFFF"
            )

    def test_light_os_gives_black_icons(self, app):
        manager = _manager()
        with patch(
            "moleditpy.ui.main_window_init.detect_system_dark_mode", return_value=False
        ):
            assert MainInitManager._get_icon_foreground_color(manager) == QColor(
                "#000000"
            )

    def test_explicit_setting_wins_without_asking_the_os(self, app):
        manager = _manager({"icon_foreground": "#FF0000"})
        with patch("moleditpy.ui.main_window_init.detect_system_dark_mode") as detect:
            assert MainInitManager._get_icon_foreground_color(manager) == QColor(
                "#FF0000"
            )
        detect.assert_not_called()

    def test_falls_back_to_background_luminance(self, app):
        """With no OS answer, a dark canvas still has to yield light icons."""
        manager = _manager({"background_color": "#101010"})
        with patch(
            "moleditpy.ui.main_window_init.detect_system_dark_mode", return_value=None
        ):
            assert MainInitManager._get_icon_foreground_color(manager) == QColor(
                "#FFFFFF"
            )
