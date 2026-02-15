import pytest
import sys
from unittest.mock import patch, MagicMock
from moleditpy.modules import sip_isdeleted_safe

def test_sip_isdeleted_safe_none():
    """Test safe check with None."""
    assert sip_isdeleted_safe(None) is False

def test_sip_isdeleted_safe_valid_obj():
    """Test safe check with a valid object (mocked)."""
    obj = MagicMock()
    # If sip.isdeleted returns False (not deleted)
    with patch("moleditpy.modules._sip_isdeleted", return_value=False):
        assert sip_isdeleted_safe(obj) is False

def test_sip_isdeleted_safe_deleted_obj():
    """Test safe check with a deleted object."""
    obj = MagicMock()
    # If sip.isdeleted returns True (deleted)
    with patch("moleditpy.modules._sip_isdeleted", return_value=True):
        assert sip_isdeleted_safe(obj) is True

def test_sip_isdeleted_safe_exception():
    """Test safe check when an exception occurs."""
    obj = MagicMock()
    with patch("moleditpy.modules._sip_isdeleted", side_effect=Exception("SIP error")):
        assert sip_isdeleted_safe(obj) is False

def test_sip_isdeleted_safe_no_sip():
    """Test safe check when _sip_isdeleted is None (sip import failed)."""
    with patch("moleditpy.modules._sip_isdeleted", None):
        assert sip_isdeleted_safe(MagicMock()) is False
