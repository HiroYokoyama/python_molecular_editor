"""Tests for PluginManager — metadata extraction, registration, and discovery."""

import os
import textwrap
from moleditpy.plugins.plugin_manager import PluginManager
from unittest.mock import MagicMock


# =============================================================================
# get_plugin_info_safe (AST-based metadata extraction)
# =============================================================================


def test_plugin_info_extracts_all_fields(tmp_path):
    """Verify AST parser extracts PLUGIN_NAME, VERSION, AUTHOR, DESCRIPTION."""
    plugin_file = tmp_path / "test_plugin.py"
    plugin_file.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "My Cool Plugin"
        PLUGIN_VERSION = "1.2.3"
        PLUGIN_AUTHOR = "Test Author"
        PLUGIN_DESCRIPTION = "Does amazing things"
        PLUGIN_CATEGORY = "Analysis"

        def run(ctx):
            pass
    """)
    )

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))

    assert info["name"] == "My Cool Plugin"
    assert info["version"] == "1.2.3"
    assert info["author"] == "Test Author"
    assert info["description"] == "Does amazing things"
    assert info["category"] == "Analysis"


def test_plugin_info_version_tuple(tmp_path):
    """Version tuples like (1, 0, 0) should be joined as '1.0.0'."""
    plugin_file = tmp_path / "tuple_ver.py"
    plugin_file.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "Tuple Version"
        PLUGIN_VERSION = (2, 3, 1)
    """)
    )

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))
    assert info["version"] == "2.3.1"


def test_plugin_info_fallback_dunder(tmp_path):
    """__version__ and __author__ should be used when PLUGIN_ variants are missing."""
    plugin_file = tmp_path / "dunder_plugin.py"
    plugin_file.write_text(
        textwrap.dedent("""\
        __version__ = "0.5.0"
        __author__ = "Dunder Author"
        PLUGIN_NAME = "Dunder Plugin"
    """)
    )

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))
    assert info["version"] == "0.5.0"
    assert info["author"] == "Dunder Author"


def test_plugin_info_docstring_fallback(tmp_path):
    """Module docstring should be used as description if PLUGIN_DESCRIPTION is missing."""
    plugin_file = tmp_path / "docstring_plugin.py"
    plugin_file.write_text(
        textwrap.dedent('''\
        """This is the plugin description from docstring."""
        PLUGIN_NAME = "Docstring Plugin"
    ''')
    )

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))
    assert info["description"] == "This is the plugin description from docstring."


def test_plugin_info_missing_file(tmp_path):
    """Non-existent file should return defaults, not crash."""
    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(tmp_path / "nonexistent.py"))
    assert info["name"] == "nonexistent.py"
    assert info["version"] == "Unknown"


def test_plugin_info_syntax_error(tmp_path):
    """File with syntax errors should return defaults, not crash."""
    plugin_file = tmp_path / "broken.py"
    plugin_file.write_text("def foo(:\n    pass\n")

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))
    assert info["version"] == "Unknown"


def test_plugin_info_empty_file(tmp_path):
    """Empty file should return defaults."""
    plugin_file = tmp_path / "empty.py"
    plugin_file.write_text("")

    pm = PluginManager()
    info = pm.get_plugin_info_safe(str(plugin_file))
    assert info["name"] == "empty.py"


# =============================================================================
# Registration methods
# =============================================================================


def test_register_menu_action():
    """register_menu_action should store action metadata."""
    pm = PluginManager()
    pm.register_menu_action(
        "TestPlugin", "Tools/Test", lambda: None, "Test", None, None
    )
    assert len(pm.menu_actions) == 1
    assert pm.menu_actions[0]["plugin"] == "TestPlugin"


def test_register_toolbar_action():
    pm = PluginManager()
    pm.register_toolbar_action("TestPlugin", lambda: None, "Test", None, "Tooltip")
    assert len(pm.toolbar_actions) == 1


def test_register_export_action():
    pm = PluginManager()
    pm.register_export_action("TestPlugin", "Export as PDF", lambda: None)
    assert len(pm.export_actions) == 1
    assert pm.export_actions[0]["label"] == "Export as PDF"


def test_register_optimization_method():
    pm = PluginManager()
    pm.register_optimization_method("TestPlugin", "MMFF94", lambda: None)
    assert "MMFF94" in pm.optimization_methods


def test_register_file_opener():
    pm = PluginManager()
    pm.register_file_opener("TestPlugin", ".cif", lambda: None, priority=10)
    assert ".cif" in pm.file_openers


def test_register_file_opener_priority():
    """Higher priority opener should replace lower priority one in sorted order."""
    pm = PluginManager()
    cb_low = lambda: "low"
    cb_high = lambda: "high"
    pm.register_file_opener("Plugin1", ".xyz", cb_low, priority=0)
    pm.register_file_opener("Plugin2", ".xyz", cb_high, priority=10)
    # It's a list, sorted by priority descending
    assert pm.file_openers[".xyz"][0]["callback"] == cb_high


def test_register_analysis_tool():
    pm = PluginManager()
    pm.register_analysis_tool("TestPlugin", "NMR Analysis", lambda: None)
    assert len(pm.analysis_tools) == 1


def test_register_save_handler():
    pm = PluginManager()
    pm.register_save_handler("TestPlugin", lambda: None)
    assert "TestPlugin" in pm.save_handlers


def test_register_load_handler():
    pm = PluginManager()
    pm.register_load_handler("TestPlugin", lambda: None)
    assert "TestPlugin" in pm.load_handlers


def test_register_3d_style():
    pm = PluginManager()
    pm.register_3d_style("TestPlugin", "Wireframe", lambda: None)
    assert "Wireframe" in pm.custom_3d_styles


def test_register_document_reset_handler():
    pm = PluginManager()
    pm.register_document_reset_handler("TestPlugin", lambda: None)
    assert len(pm.document_reset_handlers) == 1


def test_invoke_document_reset_handlers():
    """All registered reset handlers should be called."""
    pm = PluginManager()
    called = []
    pm.register_document_reset_handler("P1", lambda: called.append("P1"))
    pm.register_document_reset_handler("P2", lambda: called.append("P2"))

    pm.invoke_document_reset_handlers()
    assert called == ["P1", "P2"]


# =============================================================================
# Discovery
# =============================================================================


def test_discover_plugins_empty_dir(tmp_path):
    """Empty plugin directory should return no plugins."""
    pm = PluginManager()
    pm.plugin_dir = str(tmp_path)
    plugins = pm.discover_plugins()
    assert plugins == []


def test_discover_plugins_single_file(tmp_path):
    """Single .py file in plugin dir should be discovered."""
    plugin_file = tmp_path / "hello_plugin.py"
    plugin_file.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "Hello"
        PLUGIN_VERSION = "1.0"

        def register(ctx):
            pass
    """)
    )

    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(tmp_path)
    plugins = pm.discover_plugins()

    assert len(plugins) >= 1
    names = [p["name"] for p in plugins]
    assert "Hello" in names


def test_discover_plugins_package(tmp_path):
    """Folder with __init__.py should be discovered as package plugin."""
    pkg_dir = tmp_path / "my_package"
    pkg_dir.mkdir()
    init_file = pkg_dir / "__init__.py"
    init_file.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "Package Plugin"
        PLUGIN_VERSION = "2.0"

        def register(ctx):
            pass
    """)
    )

    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(tmp_path)
    plugins = pm.discover_plugins()

    assert len(plugins) >= 1
    names = [p["name"] for p in plugins]
    assert "Package Plugin" in names


def test_discover_plugins_ignores_dunder_files(tmp_path):
    """Files starting with __ should be ignored, and root __init__.py shouldn't count as a plugin."""
    cat_dir = tmp_path / "category"
    cat_dir.mkdir()
    (cat_dir / "__init__.py").write_text("# Not a plugin")

    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(tmp_path)
    plugins = pm.discover_plugins()

    (cat_dir / "__trash.py").write_text("# ignore me")
    assert len(plugins) == 1  # The 'category' package itself
    names = [p["name"] for p in plugins]
    assert "__trash" not in names


def test_ensure_plugin_dir_creates_directory(tmp_path):
    """ensure_plugin_dir should create the directory if it doesn't exist."""
    pm = PluginManager()
    new_dir = str(tmp_path / "new_plugins")
    pm.plugin_dir = new_dir
    pm.ensure_plugin_dir()
    assert os.path.isdir(new_dir)


# =============================================================================
# install_plugin + full discover flow (merged from test_plugin_manager_redundant)
# =============================================================================


def test_install_and_discover_single_file(tmp_path):
    """install_plugin copies a .py file; discover_plugins returns its metadata."""
    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()
    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(plugin_dir)

    source = tmp_path / "test_plugin.py"
    source.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "Test Plugin"
        PLUGIN_VERSION = "1.0"
        PLUGIN_DESCRIPTION = "A simple test plugin"

        def initialize(context):
            pass
    """)
    )

    success, msg = pm.install_plugin(str(source))
    assert success, f"Installation failed: {msg}"
    assert (plugin_dir / "test_plugin.py").exists()

    plugins = pm.discover_plugins()
    assert len(plugins) == 1
    p = plugins[0]
    assert p["name"] == "Test Plugin"
    assert p["version"] == "1.0"
    assert p["description"] == "A simple test plugin"
    assert p["status"] == "Loaded"


def test_install_plugin_registers_menu_action(tmp_path):
    """A plugin that calls context.add_menu_action in initialize() should register it."""
    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()
    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(plugin_dir)

    source = tmp_path / "action_plugin.py"
    source.write_text(
        textwrap.dedent("""\
        PLUGIN_NAME = "Action Plugin"

        def initialize(context):
            context.add_menu_action(
                path="Test > Action",
                callback=lambda: None,
                text="Test Action",
                icon=None,
                shortcut=None
            )
    """)
    )

    pm.install_plugin(str(source))
    pm.discover_plugins()

    assert len(pm.menu_actions) == 1
    assert pm.menu_actions[0]["plugin"] == "Action Plugin"
    assert pm.menu_actions[0]["text"] == "Test Action"


def test_install_zip_extracts_and_discovers(tmp_path):
    """install_plugin accepts a .zip file; extracted package is discovered."""
    import zipfile

    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()
    pm = PluginManager(main_window=MagicMock())
    pm.plugin_dir = str(plugin_dir)

    zip_source = tmp_path / "zip_source"
    zip_source.mkdir()
    (zip_source / "__init__.py").write_text(
        "PLUGIN_NAME='Zipped Plugin'", encoding="utf-8"
    )

    zip_file = tmp_path / "plugin.zip"
    with zipfile.ZipFile(zip_file, "w") as zf:
        zf.write(zip_source / "__init__.py", "MyPlugin/__init__.py")

    success, msg = pm.install_plugin(str(zip_file))
    assert success, f"ZIP installation failed: {msg}"
    assert (plugin_dir / "MyPlugin" / "__init__.py").exists()

    plugins = pm.discover_plugins()
    assert any(p["name"] == "Zipped Plugin" for p in plugins)
