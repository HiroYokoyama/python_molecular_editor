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
        assert sip_isdeleted_safe(None) is True

    def test_returns_false_when_sip_unavailable(self):
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", None):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_true_when_sip_reports_deleted(self):
        mock_checker = MagicMock(return_value=True)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is True

    def test_returns_false_when_sip_reports_not_deleted(self):
        mock_checker = MagicMock(return_value=False)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_runtime_error(self):
        mock_checker = MagicMock(side_effect=RuntimeError("deleted"))
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_attribute_error(self):
        mock_checker = MagicMock(side_effect=AttributeError)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False

    def test_returns_false_on_type_error(self):
        mock_checker = MagicMock(side_effect=TypeError)
        with patch("moleditpy.utils.sip_isdeleted_safe._sip_isdeleted", mock_checker):
            assert sip_isdeleted_safe(object()) is False
