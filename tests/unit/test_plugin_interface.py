import pytest
from unittest.mock import MagicMock, patch
from moleditpy.plugins.plugin_interface import PluginContext, Plugin3DController


class TestPluginInterface:
    @pytest.fixture
    def mock_manager(self):
        return MagicMock()

    @pytest.fixture
    def mock_main_window(self):
        mw = MagicMock()
        mw.view_3d_manager.current_mol = "mock_molecule"
        return mw

    def test_plugin_context_init(self, mock_manager):
        """Test PluginContext initialization."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        assert ctx._manager == mock_manager
        assert ctx._plugin_name == "TestPlugin"

    def test_add_menu_action(self, mock_manager):
        """Test add_menu_action delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.add_menu_action("File/Test", callback, "Test Action", "icon.png", "Ctrl+T")
        mock_manager.register_menu_action.assert_called_once_with(
            "TestPlugin", "File/Test", callback, "Test Action", "icon.png", "Ctrl+T"
        )

    def test_add_toolbar_action(self, mock_manager):
        """Test add_toolbar_action delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.add_toolbar_action(callback, "Toolbar Action", "icon.png", "Tooltip")
        mock_manager.register_toolbar_action.assert_called_once_with(
            "TestPlugin", callback, "Toolbar Action", "icon.png", "Tooltip"
        )

    def test_register_drop_handler(self, mock_manager):
        """Test register_drop_handler delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_drop_handler(callback, 5)
        mock_manager.register_drop_handler.assert_called_once_with(
            "TestPlugin", callback, 5
        )

    def test_get_3d_controller(self, mock_manager, mock_main_window):
        """Test get_3d_controller returns a controller linked to main window."""
        mock_manager.get_main_window.return_value = mock_main_window
        ctx = PluginContext(mock_manager, "TestPlugin")
        controller = ctx.get_3d_controller()
        assert isinstance(controller, Plugin3DController)
        assert controller._mw == mock_main_window

    def test_get_main_window(self, mock_manager, mock_main_window):
        """Test get_main_window delegation."""
        mock_manager.get_main_window.return_value = mock_main_window
        ctx = PluginContext(mock_manager, "TestPlugin")
        assert ctx.get_main_window() == mock_main_window

    def test_current_molecule_property(self, mock_manager, mock_main_window):
        """Test current_molecule getter and setter."""
        mock_manager.get_main_window.return_value = mock_main_window
        ctx = PluginContext(mock_manager, "TestPlugin")

        # Test getter
        assert ctx.current_molecule == "mock_molecule"

        # Test setter
        ctx.current_molecule = "new_molecule"
        assert mock_main_window.view_3d_manager.current_mol == "new_molecule"
        mock_main_window.view_3d_manager.draw_molecule_3d.assert_called_once_with("new_molecule")

    def test_current_molecule_no_window(self, mock_manager):
        """Test current_molecule when main window is None."""
        mock_manager.get_main_window.return_value = None
        ctx = PluginContext(mock_manager, "TestPlugin")

        assert ctx.current_molecule is None

        # Setter should safely do nothing
        ctx.current_molecule = "fail"
        # No error raised is the pass condition

    def test_add_export_action(self, mock_manager):
        """Test add_export_action delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.add_export_action("Export Plugin", callback)
        mock_manager.register_export_action.assert_called_once_with(
            "TestPlugin", "Export Plugin", callback
        )

    def test_register_optimization_method(self, mock_manager):
        """Test register_optimization_method delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_optimization_method("My Opt", callback)
        mock_manager.register_optimization_method.assert_called_once_with(
            "TestPlugin", "My Opt", callback
        )

    def test_register_file_opener(self, mock_manager):
        """Test register_file_opener delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_file_opener(".ext", callback, 10)
        mock_manager.register_file_opener.assert_called_once_with(
            "TestPlugin", ".ext", callback, 10
        )

    def test_add_analysis_tool(self, mock_manager):
        """Test add_analysis_tool delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.add_analysis_tool("Analyze This", callback)
        mock_manager.register_analysis_tool.assert_called_once_with(
            "TestPlugin", "Analyze This", callback
        )

    def test_register_save_handler(self, mock_manager):
        """Test register_save_handler delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_save_handler(callback)
        mock_manager.register_save_handler.assert_called_once_with(
            "TestPlugin", callback
        )

    def test_register_load_handler(self, mock_manager):
        """Test register_load_handler delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_load_handler(callback)
        mock_manager.register_load_handler.assert_called_once_with(
            "TestPlugin", callback
        )

    def test_register_3d_context_menu(self, mock_manager, capsys):
        """Test deprecated register_3d_context_menu."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        ctx.register_3d_context_menu(MagicMock(), "Label")
        captured = capsys.readouterr()
        assert "deprecated" in captured.out or "deprecated" in captured.err

    def test_register_3d_style(self, mock_manager):
        """Test register_3d_style delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_3d_style("My Style", callback)
        mock_manager.register_3d_style.assert_called_once_with(
            "TestPlugin", "My Style", callback
        )

    def test_register_document_reset_handler(self, mock_manager):
        """Test register_document_reset_handler delegation."""
        ctx = PluginContext(mock_manager, "TestPlugin")
        callback = MagicMock()
        ctx.register_document_reset_handler(callback)
        mock_manager.register_document_reset_handler.assert_called_once_with(
            "TestPlugin", callback
        )

    def test_3d_controller_set_atom_color(self, mock_main_window):
        """Test Plugin3DController.set_atom_color."""
        controller = Plugin3DController(mock_main_window)
        # Mock the view_3d and plotter
        mock_main_window.main_window_view_3d = MagicMock()
        mock_main_window.plotter = MagicMock()

        controller.set_atom_color(1, "#FF0000")

        mock_main_window.main_window_view_3d.update_atom_color_override.assert_called_once_with(
            1, "#FF0000"
        )
        mock_main_window.plotter.render.assert_called_once()

    def test_3d_controller_set_bond_color(self, mock_main_window):
        """Test Plugin3DController.set_bond_color."""
        controller = Plugin3DController(mock_main_window)
        # Mock the view_3d and plotter
        mock_main_window.main_window_view_3d = MagicMock()
        mock_main_window.plotter = MagicMock()

        controller.set_bond_color(2, "#00FF00")

        mock_main_window.main_window_view_3d.update_bond_color_override.assert_called_once_with(
            2, "#00FF00"
        )
        mock_main_window.plotter.render.assert_called_once()
