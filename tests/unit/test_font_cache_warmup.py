"""The startup font-cache warm-up that keeps the first heteroatom responsive."""

from unittest.mock import MagicMock, patch

from PyQt6.QtGui import QFont

from moleditpy.ui.atom_item import SUBSCRIPT_MAP, _WARM_UP_TEXT, warm_font_cache


class TestSubscriptMap:
    def test_maps_every_digit(self):
        assert "0123456789".translate(SUBSCRIPT_MAP) == "₀₁₂₃₄₅₆₇₈₉"

    def test_labels_render_as_subscripts(self):
        assert "H" + "2".translate(SUBSCRIPT_MAP) == "H₂"


class TestWarmUpText:
    def test_covers_every_subscript_digit(self):
        """One warm-up must serve any implicit-H count, not just H₂."""
        for digit in "₀₁₂₃₄₅₆₇₈₉":
            assert digit in _WARM_UP_TEXT, f"{digit} would still cost a cache miss"

    def test_covers_latin_shaping_too(self):
        assert "C" in _WARM_UP_TEXT and "H" in _WARM_UP_TEXT


class TestWarmFontCache:
    def test_measures_the_warm_up_text(self, app):
        with patch("moleditpy.ui.atom_item.QFontMetricsF") as fm:
            warm_font_cache()
        fm.return_value.boundingRect.assert_called_once_with(_WARM_UP_TEXT)

    def test_uses_the_font_it_is_given(self, app):
        font = QFont("Times New Roman", 14)
        with patch("moleditpy.ui.atom_item.QFontMetricsF") as fm:
            warm_font_cache(font)
        assert fm.call_args[0][0] is font

    def test_never_raises(self, app):
        """A cosmetic optimisation must not be able to break startup."""
        with patch("moleditpy.ui.atom_item.QFontMetricsF", side_effect=RuntimeError):
            warm_font_cache()


class TestInitManagerHook:
    def test_warms_with_the_configured_font(self, app):
        from moleditpy.ui.main_window_init import MainInitManager

        manager = MagicMock(spec=MainInitManager)
        manager.scene = MagicMock()
        manager.scene.get_setting.side_effect = lambda key, default: {
            "atom_font_family_2d": "Times New Roman",
            "atom_font_size_2d": 14,
            "atom_font_bold_2d": False,
        }.get(key, default)

        with patch("moleditpy.ui.main_window_init.warm_font_cache") as warm:
            MainInitManager._warm_font_cache(manager)

        font = warm.call_args[0][0]
        assert font.family() == "Times New Roman"
        assert font.pointSize() == 14
        assert font.weight() == QFont.Weight.Normal

    def test_falls_back_before_the_scene_exists(self, app):
        """Warm-up can run before the scene is built, so a missing scene must
        fall back to the default font rather than raise.

        This previously asserted a scene that exists but lacks get_setting();
        that path was a hasattr() guard for test fakes only — the app's only
        scene is MoleculeScene, which always has it."""
        from moleditpy.ui.main_window_init import MainInitManager

        manager = MagicMock(spec=MainInitManager)
        manager.scene = None

        with patch("moleditpy.ui.main_window_init.warm_font_cache") as warm:
            MainInitManager._warm_font_cache(manager)

        warm.assert_called_once_with(None)
