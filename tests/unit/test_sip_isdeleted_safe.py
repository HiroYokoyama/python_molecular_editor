"""
Unit tests for utils/sip_isdeleted_safe.py.

Covers:
  - None input -> True (treated as deleted)
  - _sip_isdeleted unavailable -> False (conservative)
  - sip reports deleted -> True
  - sip reports not deleted -> False
  - sip raises exception -> False (conservative)
"""

import os
import sys
from unittest.mock import MagicMock, patch

_src = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "moleditpy", "src")
)
if os.path.isdir(_src) and _src not in sys.path:
    sys.path.insert(0, _src)

from moleditpy.utils.sip_isdeleted_safe import sip_isdeleted_safe


class TestSipIsDeletedSafe:
    def test_none_is_treated_as_deleted(self):
        """None input is treated as a deleted object, returning True."""
        assert sip_isdeleted_safe(None) is True

    def test_returns_false_when_sip_unavailable(self):
        """Returns False (not deleted) when the sip checker is not available."""
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", None):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_true_when_sip_reports_deleted(self):
        """Returns True when sip confirms the object is deleted."""
        mock_checker = MagicMock(return_value=True)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is True

    def test_returns_false_when_sip_reports_not_deleted(self):
        """Returns False when sip confirms the object is still alive."""
        mock_checker = MagicMock(return_value=False)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_runtime_error(self):
        """Returns False conservatively when sip raises RuntimeError."""
        mock_checker = MagicMock(side_effect=RuntimeError("deleted"))
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_attribute_error(self):
        """Returns False conservatively when sip raises AttributeError."""
        mock_checker = MagicMock(side_effect=AttributeError)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_type_error(self):
        """Returns False conservatively when sip raises TypeError."""
        mock_checker = MagicMock(side_effect=TypeError)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_fallback_to_standard_sip(self):
        """Falls back to standard sip.isdeleted when PyQt6.sip is unavailable."""
        import sys
        import importlib
        from unittest.mock import patch, MagicMock

        # Mock standard sip module
        mock_sip = MagicMock()
        mock_sip.isdeleted = MagicMock(return_value=True)

        with patch.dict(
            sys.modules, {"PyQt6": None, "PyQt6.sip": None, "sip": mock_sip}
        ):
            import moleditpy.utils.sip_isdeleted_safe as sds

            importlib.reload(sds)
            assert sds.sip_isdeleted_safe(object()) is True

        # Restore standard state
        importlib.reload(sds)

    def test_fallback_when_both_sip_packages_missing(self):
        """Returns False and sets _sip_isdeleted=None when both sip packages are absent."""
        import sys
        import importlib
        from unittest.mock import patch

        with patch.dict(sys.modules, {"PyQt6": None, "PyQt6.sip": None, "sip": None}):
            import moleditpy.utils.sip_isdeleted_safe as sds

            importlib.reload(sds)
            assert sds.sip_isdeleted_safe(object()) is False
            assert sds._sip_isdeleted is None

        # Restore standard state
        importlib.reload(sds)
