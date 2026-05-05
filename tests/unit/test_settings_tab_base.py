from unittest.mock import MagicMock
from PyQt6.QtWidgets import QFrame, QSlider, QLabel, QWidget
from moleditpy.ui.settings_tabs.settings_tab_base import SettingsTabBase
from moleditpy.utils.default_settings import DEFAULT_SETTINGS


class ConcreteTab(SettingsTabBase):
    """Minimal concrete subclass for testing the base class."""

    def update_ui(self, settings_dict):
        pass

    def get_settings(self):
        return {}


def test_init_stores_default_settings(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    assert tab.default_settings is DEFAULT_SETTINGS


def test_create_separator_returns_hline(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    sep = tab._create_separator()
    assert isinstance(sep, QFrame)
    assert sep.frameShape() == QFrame.Shape.HLine
    assert sep.frameShadow() == QFrame.Shadow.Sunken


def test_create_slider_range(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    slider, label = tab._create_slider(10, 200, 10.0)
    assert isinstance(slider, QSlider)
    assert slider.minimum() == 10
    assert slider.maximum() == 200


def test_create_slider_float_label_updates(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    slider, label = tab._create_slider(0, 100, 10.0)
    slider.setValue(50)
    assert label.text() == "5.00"


def test_create_slider_int_label_updates(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    slider, label = tab._create_slider(3, 20, 1.0, is_int=True)
    slider.setValue(10)
    assert label.text() == "10"


def test_wrap_layout_returns_widget_with_children(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    slider, label = tab._create_slider(0, 100, 1.0)
    container = tab._wrap_layout(slider, label)
    assert isinstance(container, QWidget)


def test_reset_to_defaults_calls_update_ui(app):
    tab = ConcreteTab(DEFAULT_SETTINGS)
    tab.update_ui = MagicMock()
    tab.reset_to_defaults()
    tab.update_ui.assert_called_once_with(DEFAULT_SETTINGS)


def test_update_ui_not_implemented(app):
    import pytest

    # The base class raises NotImplementedError; ConcreteTab overrides it
    class RawBase(SettingsTabBase):
        pass

    base = RawBase(DEFAULT_SETTINGS)
    with pytest.raises(NotImplementedError):
        base.update_ui({})


def test_get_settings_not_implemented(app):
    import pytest

    class RawBase(SettingsTabBase):
        pass

    base = RawBase(DEFAULT_SETTINGS)
    with pytest.raises(NotImplementedError):
        base.get_settings()
