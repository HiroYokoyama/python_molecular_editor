"""
Unit tests for utils/constants.py.
"""

import importlib
import builtins
from unittest.mock import patch, MagicMock, mock_open

# We want to test different import and file states for constants._get_version()
# using importlib.reload.


def test_constants_version_from_metadata():
    """Test when importlib.metadata successfully finds the package version."""
    mock_version = MagicMock(return_value="1.2.3")

    with patch("importlib.metadata.version", mock_version):
        import moleditpy.utils.constants as constants

        importlib.reload(constants)
        assert constants.VERSION == "1.2.3"
        assert mock_version.call_count >= 1


def test_constants_version_metadata_package_not_found():
    """Test fallback to pyproject.toml when PackageNotFoundError is raised."""
    from importlib.metadata import PackageNotFoundError

    def raise_not_found(*args, **kwargs):
        raise PackageNotFoundError()

    # Mock os.path.exists and builtins.open to return a mock pyproject.toml
    mock_exists = MagicMock(side_effect=lambda path: "pyproject.toml" in path)
    toml_data = 'version = "2.3.4"\n'
    m_open = mock_open(read_data=toml_data)

    with (
        patch("importlib.metadata.version", side_effect=raise_not_found),
        patch("os.path.exists", mock_exists),
        patch("builtins.open", m_open),
    ):
        import moleditpy.utils.constants as constants

        importlib.reload(constants)
        assert constants.VERSION == "2.3.4"


def test_constants_version_import_error_fallback():
    """Test fallback to pyproject.toml when importlib.metadata cannot be imported (ImportError)."""
    # Force importlib.metadata to raise ImportError when imported inside _get_version
    # We can do this by patching sys.modules or importlib.import_module.
    # But since _get_version does `from importlib.metadata import version`,
    # let's patch builtins.__import__ or mock sys.modules.

    original_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "importlib.metadata" or (
            args and args[0] == "metadata" and name == "importlib"
        ):
            raise ImportError("Mocked import error")
        return original_import(name, *args, **kwargs)

    mock_exists = MagicMock(side_effect=lambda path: "pyproject.toml" in path)
    toml_data = 'version = "3.4.5"\n'
    m_open = mock_open(read_data=toml_data)

    with (
        patch("builtins.__import__", side_effect=mock_import),
        patch("os.path.exists", mock_exists),
        patch("builtins.open", m_open),
    ):
        import moleditpy.utils.constants as constants

        importlib.reload(constants)
        assert constants.VERSION == "3.4.5"


def test_constants_version_all_fail_returns_unknown():
    """Test that it falls back to 'Unknown' when all lookup methods fail."""
    from importlib.metadata import PackageNotFoundError

    def raise_not_found(*args, **kwargs):
        raise PackageNotFoundError()

    mock_exists = MagicMock(return_value=False)

    with (
        patch("importlib.metadata.version", side_effect=raise_not_found),
        patch("os.path.exists", mock_exists),
    ):
        import moleditpy.utils.constants as constants

        importlib.reload(constants)
        assert constants.VERSION == "Unknown"


def test_constants_version_file_exception_returns_unknown():
    """Test that any reading file exception is handled and returns 'Unknown'."""
    from importlib.metadata import PackageNotFoundError

    def raise_not_found(*args, **kwargs):
        raise PackageNotFoundError()

    mock_exists = MagicMock(side_effect=lambda path: "pyproject.toml" in path)
    m_open = MagicMock(side_effect=PermissionError("Permission Denied"))

    with (
        patch("importlib.metadata.version", side_effect=raise_not_found),
        patch("os.path.exists", mock_exists),
        patch("builtins.open", m_open),
    ):
        import moleditpy.utils.constants as constants

        importlib.reload(constants)
        assert constants.VERSION == "Unknown"


def test_restore_state():
    """Helper to restore constants state after reloading it under mocks."""
    import moleditpy.utils.constants as constants

    importlib.reload(constants)
