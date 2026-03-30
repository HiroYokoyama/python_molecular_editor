import pytest
import os
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QAction, QActionGroup, QColor, QBrush
from PyQt6.QtWidgets import QMenu, QToolBar, QApplication, QFileDialog

def test_plugin_menu_actions_population(window, qtbot):
    """Test that menu actions registered by plugins are correctly added to the menu bar."""
    # 1. Setup mock plugin data
    # Create the top-level Plugins menu if it doesn't exist
    menu_bar = window.menuBar()
    
    # Ensure existing menu to test path merging
    existing_menu = menu_bar.addMenu("Analysis")
    
    window.plugin_manager.menu_actions = [
        {
            "plugin": "TestPlugin",
            "path": "Plugins/SubMenu/TestAction",
            "callback": MagicMock(),
            "text": "Click Me",
            "icon": None,
            "shortcut": "Ctrl+Shift+P"
        },
        {
            "plugin": "TestPlugin",
            "path": "Analysis/ExtraTool",
            "callback": MagicMock(),
            "text": "Plugin Analysis",
            "icon": None,
            "shortcut": None
        }
    ]
    
    # 2. Trigger update
    window.init_manager.update_plugin_menu(window.init_manager.plugin_menu)
    
    # 3. Verify top-level "Plugins" menu was found/created
    plugins_action = next((a for a in menu_bar.actions() if "Plugins" in a.text().replace("&", "")), None)
    assert plugins_action is not None
    plugins_menu = plugins_action.menu()
    assert plugins_menu is not None
    
    # 4. Verify SubMenu
    sub_menu_action = next((a for a in plugins_menu.actions() if "SubMenu" in a.text().replace("&", "")), None)
    assert sub_menu_action is not None
    sub_menu = sub_menu_action.menu()
    assert sub_menu is not None
    
    # 5. Verify Action
    test_action = next((a for a in sub_menu.actions() if "Click Me" in a.text().replace("&", "")), None)
    assert test_action is not None
    assert test_action.shortcut().toString() == "Ctrl+Shift+P"
    
    # 6. Verify merging into existing Analysis menu
    analysis_action = next((a for a in menu_bar.actions() if "Analysis" in a.text().replace("&", "")), None)
    analysis_menu = analysis_action.menu()
    plugin_analysis_action = next((a for a in analysis_menu.actions() if "Plugin Analysis" in a.text().replace("&", "")), None)
    assert plugin_analysis_action is not None

def test_plugin_toolbar_actions_visibility(window, qtbot):
    """Test that the plugin toolbar is shown/hidden and populated correctly."""
    # Ensure toolbar exists
    if not hasattr(window, "plugin_toolbar"):
        window.init_manager.plugin_toolbar = QToolBar("Plugins", window)
        window.addToolBar(window.init_manager.plugin_toolbar)
        
    # 1. Empty plugins -> Toolbar hidden
    window.plugin_manager.toolbar_actions = []
    window.init_manager._add_plugin_toolbar_actions()
    assert window.init_manager.plugin_toolbar.isHidden()
    
    # 2. With plugins -> Toolbar shown
    mock_cb = MagicMock()
    window.plugin_manager.toolbar_actions = [
        {"text": "ToolBtn", "callback": mock_cb, "icon": None, "tooltip": "Hint"}
    ]
    window.init_manager._add_plugin_toolbar_actions()
    assert not window.init_manager.plugin_toolbar.isHidden()
    
    # Verify action
    toolbar_actions = window.init_manager.plugin_toolbar.actions()
    assert len(toolbar_actions) == 1
    assert toolbar_actions[0].text() == "ToolBtn"
    assert toolbar_actions[0].toolTip() == "Hint"

def test_ui_sync_after_reset(window, qtbot):
    """Test that UI elements (background, checked states) sync correctly after settings reset."""
    # 1. Setup stale state
    window.init_manager.settings["background_color_2d"] = "#FF0000"
    window.init_manager.settings["optimization_method"] = "UFF_RDKIT"
    
    # Ensure scene exists and is a real MoleculeScene
    assert window.init_manager.scene is not None
    
    # 2. Simulate new settings arrival
    window.init_manager.settings["background_color_2d"] = "#0000FF" 
    window.init_manager.settings["optimization_method"] = "MMFF_RDKIT"
    window.init_manager.settings["3d_conversion_mode"] = "rdkit"
    
    # 3. Trigger UI refresh
    window.init_manager._refresh_ui_after_reset()
    
    # 4. Verify Scene Background
    actual_bg = window.init_manager.scene.backgroundBrush().color().name().upper()
    assert actual_bg == "#0000FF", f"Scene background should be #0000FF, but got {actual_bg}"
    
    # 5. Verify Menu Checkstates (using existing window.opt3d_actions if available)
    if hasattr(window, "opt3d_actions"):
        if "MMFF_RDKIT" in window.opt3d_actions:
            assert window.opt3d_actions["MMFF_RDKIT"].isChecked(), "MMFF action should be checked"
    
    if hasattr(window, "conv_actions"):
        if "rdkit" in window.conv_actions:
            assert window.conv_actions["rdkit"].isChecked(), "RDKit conversion should be checked"

def test_custom_3d_style_integration(window, qtbot):
    """Test that custom 3D styles from plugins appear in the style menu."""
    # Ensure style_button and menu exist
    if not hasattr(window, "style_button") or not window.init_manager.style_button.menu():
        # Fallback for headless/mocked plotter where style_button might not be fully initialized
        window.init_manager.style_button = QToolButton(window)
        window.init_manager.style_button.setMenu(QMenu(window))
    
    style_menu = window.init_manager.style_button.menu()
    style_group = QActionGroup(window)
    style_group.setExclusive(True)
    
    # Add a dummy existing style to check grouping
    existing_style = QAction("Stick", window, checkable=True)
    existing_style.setActionGroup(style_group)
    style_menu.addAction(existing_style)
    
    window.plugin_manager.custom_3d_styles = {"Vantablack": {"plugin": "CoolPlugin", "callback": MagicMock()}}
    window.init_manager._update_style_menu_with_plugins()
    
    style_action = next((a for a in style_menu.actions() if "Vantablack" in a.text()), None)
    assert style_action is not None, f"Vantablack style not found in {[a.text() for a in style_menu.actions()]}"
    
    # Verify it joins a group if one exists in the menu
    assert style_action.actionGroup() is not None

def test_integrate_plugin_export_actions(window, qtbot):
    """Test that export actions are added to File/Export menu."""
    # Find existing menus
    file_action = next((a for a in window.menuBar().actions() if "File" in a.text().replace("&", "")), None)
    assert file_action is not None
    file_menu = file_action.menu()
    
    export_action = next((a for a in file_menu.actions() if a.menu() and "Export" in a.text().replace("&", "")), None)
    assert export_action is not None
    export_menu = export_action.menu()
    
    # Mock Toolbar Export Button if missing
    if not hasattr(window.init_manager, 'export_button'):
        window.init_manager.export_button = QToolButton(window)
        window.init_manager.export_button.setMenu(QMenu(window))
    btn_menu = window.init_manager.export_button.menu()
    
    window.plugin_manager.export_actions = [
        {"plugin": "Xerox", "label": "Export to Fax", "callback": MagicMock()}
    ]
    
    window.init_manager._integrate_plugin_export_actions()
    
    assert any("Export to Fax" in a.text() for a in btn_menu.actions()), "Missing in button menu"
    assert any("Export to Fax" in a.text() for a in export_menu.actions()), "Missing in Export menu"

def test_integrate_plugin_analysis_tools(window, qtbot):
    """Test integration into the Analysis menu."""
    analysis_action = next((a for a in window.menuBar().actions() if "Analysis" in a.text().replace("&", "")), None)
    assert analysis_action is not None
    analysis_menu = analysis_action.menu()
    
    window.plugin_manager.analysis_tools = [
        {"plugin": "StatBot", "label": "Calculate Karma", "callback": MagicMock()}
    ]
    
    window.init_manager._integrate_plugin_analysis_tools()
    
    found = any("Calculate Karma" in a.text() for a in analysis_menu.actions())
    assert found, f"Karma tool not found in {[a.text() for a in analysis_menu.actions()]}"

def test_integrate_plugin_file_openers_ui(window, qtbot):
    """Test integration of plugin openers into the Import menu."""
    window.import_menu = QMenu("Import", window)
    
    # Register opener for .fake extension
    cb = MagicMock()
    window.plugin_manager.file_openers = {
        ".fake": [{"plugin": "FakePlugin", "callback": cb, "priority": 10}]
    }
    
    window.init_manager._integrate_plugin_file_openers()
    
    # Verify action existence
    import_action = next((a for a in window.import_menu.actions() if "FakePlugin" in a.text()), None)
    assert import_action is not None
    assert "Import .fake" in import_action.text()
    
    # We can't easily test the QFileDialog trigger without mocking QFileDialog.getOpenFileName
    with patch.object(QFileDialog, "getOpenFileName", return_value=("data.fake", "All Files (*.)")):
        # Mock update_window_title to avoid side effects
        window.state_manager.update_window_title = MagicMock()
        import_action.trigger()
        cb.assert_called_once_with("data.fake")
        # Attributes are still proxied on host for now, but we check via io_manager if available
        assert window.current_file_path == "data.fake"
