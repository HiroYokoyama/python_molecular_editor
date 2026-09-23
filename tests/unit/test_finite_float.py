"""Unit tests for utils.finite_float."""

import pytest

from moleditpy.utils.finite_float import finite_float


@pytest.mark.parametrize("text,value", [("1.5", 1.5), (" -2 ", -2.0), ("0", 0.0)])
def test_parses_finite_numbers(text, value):
    assert finite_float(text) == value


@pytest.mark.parametrize("text", ["nan", "NaN", "inf", "-Infinity", "abc", ""])
def test_rejects_non_finite_and_non_numbers(text):
    with pytest.raises(ValueError):
        finite_float(text)
