import os
import sys
import shutil
import zipfile
import ast
import logging
from unittest.mock import MagicMock, patch

import pytest
from moleditpy.plugins.plugin_manager import PluginManager


class TestPluginManagerExtended:
    def test_imports_fallback(self):
        # To truly hit lines 29-31 we would need to reload the module or execute it locally:
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

        # line 126: if top_folder == "__init__.py" -> is_nested = False
        zip_init = str(tmpdir.join("init_only.zip"))
        with zipfile.ZipFile(zip_init, "w") as zf:
            zf.writestr("__init__.py/file.txt", "print")  # __init__.py folder!
        success, _ = pm.install_plugin(zip_init)
        assert success

        # Flat ZIP: more than 1 root -> is_nested = False
        zip2 = str(tmpdir.join("flat.zip"))
        with zipfile.ZipFile(zip2, "w") as zf:
            zf.writestr("f1.py", "")
            zf.writestr("f2.py", "")
        # Create an existing directory to test line 144-155 shutil.rmtree
        os.makedirs(os.path.join(pm.plugin_dir, "flat"))
        success, _ = pm.install_plugin(zip2)
        assert success

        # Flat ZIP existing dest is file
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
        # Instead of mocking spec, actually write a bad plugin!
        pm = PluginManager()
        pm.plugin_dir = str(tmpdir.mkdir("plugins"))

        # Bad plugin for init
        badinit = tmpdir.join("plugins").join("badinit.py")
        badinit.write("def initialize(ctx):\n    raise RuntimeError('INIT_ERR')")
        pm.discover_plugins()
        plugin = next((p for p in pm.plugins if "badinit" in p["filepath"]), None)
        assert plugin and "Error (Init): INIT_ERR" in plugin["status"]

        # Bad plugin for autorun
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
        # Fallback will parse strings normally
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

    # API MISSING METHODS Tests
    def test_manager_api_helpers(self):
        pm = PluginManager()

        # 460
        pm.register_file_opener("P1", "txt", lambda x: None)
        # 501-503 show_status_message success
        mw = MagicMock()
        pm.main_window = mw
        pm.show_status_message("test", 1000)
        mw.statusBar().showMessage.assert_called_with("test", 1000)

        # 508 push_undo_checkpoint
        pm.push_undo_checkpoint()
        mw.state_manager.push_undo_state.assert_called()

        # 517 refresh_3d_view
        pm.refresh_3d_view()
        mw.plotter.render.assert_called()

        # 526-527 reset_3d_camera
        pm.reset_3d_camera()
        mw.plotter.reset_camera.assert_called()
        mw.plotter.render.assert_called()

    def test_get_selected_atom_indices_complex(self):
        pm = PluginManager()

        pm.main_window = MagicMock()

        # 533 if not self.main_window: return []
        pm.main_window = None
        assert pm.get_selected_atom_indices() == []

        pm.main_window = MagicMock()

        # 545-562 success branch
        pm.main_window.scene = MagicMock()
        item1 = MagicMock()
        item1.atom_id = 99
        # include an item with no atom_id (e.g. bond item)
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

        # test HasProp = False
        atom.HasProp.return_value = False
        assert pm.get_selected_atom_indices() == []

        # TypeError in GetIntProp -> 565-566
        atom.HasProp.return_value = True
        atom.GetIntProp.side_effect = TypeError("err")
        assert pm.get_selected_atom_indices() == []

    def test_register_get_window(self):
        pm = PluginManager()
        pm.register_window("P1", "w1", "WIN")
        assert pm.plugin_windows["P1"]["w1"] == "WIN"
        assert pm.get_window("P1", "w1") == "WIN"

    def test_invoke_document_reset_handlers(self):
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
            # manually cause ast tuple unpacking error lines 636-639, 641-648
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
            # It will fail down in the tuple parsing
            info = pm.get_plugin_info_safe(str(f))
            assert info["version"] == "Unknown"

        # 679-682 exception
        with patch("builtins.open", side_effect=OSError("Err")):
            pm.get_plugin_info_safe(str(f))
