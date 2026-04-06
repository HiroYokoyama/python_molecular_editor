from unittest.mock import patch, MagicMock
from moleditpy.utils.sip_isdeleted_safe import sip_isdeleted_safe


def test_sip_isdeleted_safe_valid_obj():
    """Test safe check with a valid object (mocked)."""
    obj = MagicMock()
    # If sip.isdeleted returns False (not deleted)
    with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", return_value=False):
        assert sip_isdeleted_safe(obj) is False


def test_sip_isdeleted_safe_deleted_obj():
    """Test safe check with a deleted object."""
    obj = MagicMock()
    # If sip.isdeleted returns True (deleted)
    with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", return_value=True):
        assert sip_isdeleted_safe(obj) is True


def test_sip_isdeleted_safe_exception():
    """Test safe check when an exception occurs."""
    obj = MagicMock()
    with patch(
        "moleditpy.utils.sip_isdeleted_safe._sip_isdeleted",
        side_effect=RuntimeError("SIP error"),
    ):
        assert sip_isdeleted_safe(obj) is False


def test_sip_isdeleted_safe_no_sip():
    """Test safe check when _sip_isdeleted is None (sip import failed)."""
    with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", None):
        # MagicMock should return False (not deleted) when sip unavailable
        result = sip_isdeleted_safe(MagicMock())
        assert result is False
