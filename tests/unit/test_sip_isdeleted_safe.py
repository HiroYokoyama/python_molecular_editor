"""PySide wrapper lifetime checks used by core UI and legacy plugins."""

from unittest.mock import patch
from PySide6.QtCore import QObject
from shiboken6 import delete

from moleditpy.utils.sip_isdeleted_safe import sip_isdeleted_safe


def test_live_and_deleted_wrappers():
    obj = QObject()
    assert not sip_isdeleted_safe(obj)
    delete(obj)
    assert sip_isdeleted_safe(obj)
    assert sip_isdeleted_safe(None)


def test_plain_python_objects_are_safe():
    assert not sip_isdeleted_safe(object())


def test_invalid_checker_is_conservative():
    with patch("moleditpy.utils.sip_isdeleted_safe.isValid", side_effect=RuntimeError):
        assert not sip_isdeleted_safe(object())
