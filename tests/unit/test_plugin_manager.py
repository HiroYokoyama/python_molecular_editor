"""Tests for PluginManager — metadata extraction, registration, and discovery."""

import os
import sys
import ast
import zipfile
import textwrap
from moleditpy.plugins.plugin_manager import PluginManager
from unittest.mock import MagicMock, patch


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


# =============================================================================
# SHA-256 Calculation
# =============================================================================


def test_compute_sha256(tmp_path):
    """Test SHA-256 calculation for files and directories."""
    pm = PluginManager()

    # 1. Non-existent path
    assert pm.compute_sha256("/non/existent/path") == "N/A"

    # 2. Single file
    f = tmp_path / "test.txt"
    f.write_text("hello-sha256", encoding="utf-8")
    sha1 = pm.compute_sha256(str(f))
    assert sha1 != "N/A"
    assert len(sha1) == 64

    # 3. Directory
    d = tmp_path / "sha_dir"
    d.mkdir()
    (d / "a.txt").write_text("data-a", encoding="utf-8")
    (d / "b.txt").write_text("data-b", encoding="utf-8")
    sha_dir = pm.compute_sha256(str(d))
    assert sha_dir != "N/A"
    assert len(sha_dir) == 64

    # 4. Error handling (OSError)
    with patch("builtins.open", side_effect=OSError):
        assert pm._sha256_for_file(str(f)) == "N/A"
        assert pm._sha256_for_directory(str(d)) == "N/A"


# =============================================================================
# Extended tests (merged from test_plugin_manager_extended.py)
# =============================================================================


class TestPluginManagerExtended:
    def test_imports_fallback(self):
        code = """
import sys
sys.modules['moleditpy.plugins.plugin_interface'] = None
try:
    from moleditpy.plugins.plugin_manager import PluginManager
except BaseException:
    pass
"""
        exec(code)

    def test_get_set_main_window(self):
        pm = PluginManager()
        pm.set_main_window("mw")
        assert pm.get_main_window() == "mw"

    @patch("os.makedirs", side_effect=OSError("test mkdir err"))
    @patch("os.path.exists", return_value=False)
    @patch("logging.error")
    def test_ensure_plugin_dir_error(self, mock_log, mock_exists, mock_makedirs):
        pm = PluginManager()
        pm.ensure_plugin_dir()
        mock_log.assert_called_with("Error creating plugin directory: test mkdir err")

    @patch("moleditpy.plugins.plugin_manager.QDesktopServices.openUrl")
    def test_open_plugin_folder(self, mock_open_url):
        pm = PluginManager()
        pm.ensure_plugin_dir = MagicMock()
        pm.open_plugin_folder()
        mock_open_url.assert_called_once()

    @patch("shutil.rmtree")
    @patch("shutil.copytree")
    @patch("os.remove")
    @patch("os.path.exists")
    @patch("os.path.isdir")
    def test_install_plugin_folder(
        self, mock_isdir, mock_exists, mock_remove, mock_copytree, mock_rmtree
    ):
        pm = PluginManager()
        pm.ensure_plugin_dir = MagicMock()

        mock_exists.return_value = True
        mock_isdir.side_effect = [True, True]
        success, msg = pm.install_plugin("/some/folder")
        assert success
        assert "package" in msg
        mock_rmtree.assert_called()

        mock_isdir.side_effect = [True, False]
        success, msg = pm.install_plugin("/some/folder")
        mock_remove.assert_called()

    @patch("os.path.isdir", side_effect=[False, True])
    @patch("os.path.exists", return_value=True)
    @patch("shutil.copy2")
    @patch("shutil.rmtree")
    def test_install_plugin_file_existing_dir(
        self, mock_rmtree, mock_copy2, mock_exists, mock_isdir
    ):
        pm = PluginManager()
        pm.ensure_plugin_dir = MagicMock()
        pm.install_plugin("plugin.py")
        mock_rmtree.assert_called()

    @patch("os.path.basename", side_effect=RuntimeError("Install err"))
    def test_install_plugin_exception(self, mock_basename):
        pm = PluginManager()
        pm.ensure_plugin_dir = MagicMock()
        success, msg = pm.install_plugin("any")
        assert not success
        assert "Install err" in msg

    def test_zip_extraction(self, tmpdir):
        pm = PluginManager()
        pm.plugin_dir = str(tmpdir.mkdir("plugins"))
        pm.main_window = MagicMock()
        pm.discover_plugins = MagicMock()

        # Case A: nested folder
        zip1 = str(tmpdir.join("single.zip"))
        with zipfile.ZipFile(zip1, "w") as zf:
            zf.writestr("TopFolder/test.py", "print(1)")
        os.makedirs(os.path.join(pm.plugin_dir, "TopFolder"))
        success, _ = pm.install_plugin(zip1)
        assert success

        zip1a = str(tmpdir.join("single2.zip"))
        with zipfile.ZipFile(zip1a, "w") as zf:
            zf.writestr("Top2/test.py", "print(1)")
        with open(os.path.join(pm.plugin_dir, "Top2"), "w") as f:
            f.write("hello")
        success, _ = pm.install_plugin(zip1a)
        assert success

        zip_init = str(tmpdir.join("init_only.zip"))
        with zipfile.ZipFile(zip_init, "w") as zf:
            zf.writestr("__init__.py/file.txt", "print")
        success, _ = pm.install_plugin(zip_init)
        assert success

        # Flat ZIP: more than 1 root -> is_nested = False
        zip2 = str(tmpdir.join("flat.zip"))
        with zipfile.ZipFile(zip2, "w") as zf:
            zf.writestr("f1.py", "")
            zf.writestr("f2.py", "")
        os.makedirs(os.path.join(pm.plugin_dir, "flat"))
        success, _ = pm.install_plugin(zip2)
        assert success

        zip3 = str(tmpdir.join("flat2.zip"))
        with zipfile.ZipFile(zip3, "w") as zf:
            zf.writestr("f1.py", "")
            zf.writestr("f2.py", "")
        with open(os.path.join(pm.plugin_dir, "flat2"), "w") as f:
            f.write("hi")
        success, _ = pm.install_plugin(zip3)
        assert success

    def test_discover_plugins_not_exists(self):
        pm = PluginManager()
        pm.plugin_dir = "/non/existent/dir"
        with patch("os.path.exists", side_effect=[True, False]):
            assert pm.discover_plugins() == []

    @patch("importlib.util.spec_from_file_location")
    def test_load_single_plugin_exceptions_and_stub(self, mock_spec):
        pm = PluginManager()
        mock_spec.side_effect = RuntimeError("Load err")
        with patch("logging.error") as mock_log:
            pm._load_single_plugin("test.py", "test", "cat")
            mock_log.assert_called()

        mock_spec.side_effect = None
        mock_spec.return_value = None
        pm._load_single_plugin("__init__.py", "plugin", "Cat.SubCat")
        assert "Cat" in sys.modules
        assert "Cat.SubCat" in sys.modules
        sys.modules.pop("Cat", None)
        sys.modules.pop("Cat.SubCat", None)

    def test_plugin_init_exceptions(self, tmpdir):
        pm = PluginManager()
        pm.plugin_dir = str(tmpdir.mkdir("plugins"))

        badinit = tmpdir.join("plugins").join("badinit.py")
        badinit.write("def initialize(ctx):\n    raise RuntimeError('INIT_ERR')")
        pm.discover_plugins()
        plugin = next((p for p in pm.plugins if "badinit" in p["filepath"]), None)
        assert plugin and "Error (Init): INIT_ERR" in plugin["status"]

        badauto = tmpdir.mkdir("plugins2").join("badauto.py")
        badauto.write("def autorun(mw):\n    raise ValueError('AUTO_ERR')")
        pm.plugin_dir = str(tmpdir.join("plugins2"))
        pm.discover_plugins("MW")
        plugin2 = next((p for p in pm.plugins if "badauto" in p["filepath"]), None)
        assert plugin2 and "Error (Autorun): AUTO_ERR" in plugin2["status"]

    @patch("importlib.util.spec_from_file_location")
    @patch("importlib.util.module_from_spec")
    def test_load_plugin_version_tuple(self, mock_mod, mock_spec):
        pm = PluginManager()
        module = MagicMock()
        module.PLUGIN_NAME = "VPlugin"
        module.PLUGIN_VERSION = (3, 1, 4)
        del module.initialize
        del module.autorun
        mock_mod.return_value = module

        spec = MagicMock()
        spec.name = "VPlugin"
        mock_spec.return_value = spec

        pm._load_single_plugin("a.py", "VPlugin", "")
        assert pm.plugins[0]["version"] == "3.1.4"

    @patch("moleditpy.plugins.plugin_manager.QMessageBox.critical")
    def test_run_plugin_exceptions(self, mock_crit):
        pm = PluginManager()
        module = MagicMock()
        module.run.side_effect = RuntimeError("Run Err")
        pm.run_plugin(module, "mw")
        mock_crit.assert_called_once()

    def test_register_drop_handler(self):
        pm = PluginManager()
        pm.register_drop_handler("A", "c1", 1)
        pm.register_drop_handler("B", "c2", 5)
        assert pm.drop_handlers[0]["plugin"] == "B"

    def test_manager_api_helpers(self):
        pm = PluginManager()

        pm.register_file_opener("P1", "txt", lambda x: None)
        mw = MagicMock()
        pm.main_window = mw
        pm.show_status_message("test", 1000)
        mw.statusBar().showMessage.assert_called_with("test", 1000)

        pm.push_undo_checkpoint()
        mw.state_manager.push_undo_state.assert_called()

        pm.refresh_3d_view()
        mw.plotter.render.assert_called()

        pm.reset_3d_camera()
        mw.plotter.reset_camera.assert_called()
        mw.plotter.render.assert_called()

    def test_get_selected_atom_indices_complex(self):
        pm = PluginManager()

        pm.main_window = None
        assert pm.get_selected_atom_indices() == []

        pm.main_window = MagicMock()

        pm.main_window.scene = MagicMock()
        item1 = MagicMock()
        item1.atom_id = 99
        item2 = MagicMock()
        del item2.atom_id
        pm.main_window.scene.selectedItems.return_value = [item1, item2]

        mol = MagicMock()
        mol.GetNumAtoms.return_value = 1
        atom = MagicMock()
        atom.HasProp.return_value = True
        atom.GetIntProp.return_value = 99
        mol.GetAtomWithIdx.return_value = atom

        pm.main_window.current_mol = mol

        indices = pm.get_selected_atom_indices()
        assert indices == [0]

        atom.HasProp.return_value = False
        assert pm.get_selected_atom_indices() == []

        atom.HasProp.return_value = True
        atom.GetIntProp.side_effect = TypeError("err")
        assert pm.get_selected_atom_indices() == []

    def test_register_get_window(self):
        pm = PluginManager()
        pm.register_window("P1", "w1", "WIN")
        assert pm.plugin_windows["P1"]["w1"] == "WIN"
        assert pm.get_window("P1", "w1") == "WIN"

    def test_invoke_document_reset_handlers_logs_error(self):
        pm = PluginManager()

        def bad_cb():
            raise RuntimeError("Reset err")

        pm.register_document_reset_handler("P1", bad_cb)
        with patch("logging.error") as mock_log:
            pm.invoke_document_reset_handlers()
            mock_log.assert_called()

    def test_get_plugin_info_safe_exceptions_and_ast(self, tmpdir):
        pm = PluginManager()

        f = tmpdir.join("a.py")
        f.write("""
__version__: str = "1.0"
a = (1, "x")
PLUGIN_VERSION = (1, ("x",))
""")

        with patch("ast.parse") as mock_parse:

            class FakeNode:
                def __init__(self):
                    self.targets = []
                    self.target = ast.Name(id="PLUGIN_VERSION", ctx=ast.Store())

                    self.value = ast.Tuple(
                        elts=[ast.Name(id="x", ctx=ast.Load())], ctx=ast.Load()
                    )

            class FakeTree:
                body = [FakeNode()]

            mock_parse.return_value = FakeTree()
            info = pm.get_plugin_info_safe(str(f))
            assert info["version"] == "Unknown"

        with patch("builtins.open", side_effect=OSError("Err")):
            pm.get_plugin_info_safe(str(f))
