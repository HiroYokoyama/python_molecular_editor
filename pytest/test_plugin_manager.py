
import pytest
import os
import shutil
import zipfile
from unittest import mock
from PyQt6.QtWidgets import QApplication

# Adjust import path to finding modules
import sys
sys.path.append(os.path.abspath("src"))

from moleditpy.modules.plugin_manager import PluginManager

@pytest.fixture
def plugin_manager(tmp_path):
    """Fixture to provide a PluginManager instance using a temporary directory."""
    # Create the temp plugin directory
    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()
    
    # Initialize PluginManager
    pm = PluginManager()
    # Override the plugin directory to use the temp one
    pm.plugin_dir = str(plugin_dir)
    return pm

def test_init(plugin_manager):
    """Test initialization of PluginManager."""
    assert os.path.exists(plugin_manager.plugin_dir)
    assert plugin_manager.plugins == []

def test_install_and_discover_single_file(plugin_manager, tmp_path):
    """Test installing and discovering a single-file plugin."""
    # Create a dummy plugin file
    dummy_plugin_content = """
PLUGIN_NAME = "Test Plugin"
PLUGIN_VERSION = "1.0"
PLUGIN_DESCRIPTION = "A simple test plugin"

def initialize(context):
    pass
"""
    source_file = tmp_path / "test_plugin.py"
    source_file.write_text(dummy_plugin_content, encoding="utf-8")
    
    # Install the plugin
    success, msg = plugin_manager.install_plugin(str(source_file))
    assert success, f"Installation failed: {msg}"
    assert (tmp_path / "plugins" / "test_plugin.py").exists()
    
    # Discover plugins
    plugins = plugin_manager.discover_plugins()
    assert len(plugins) == 1
    p = plugins[0]
    assert p['name'] == "Test Plugin"
    assert p['version'] == "1.0"
    assert p['description'] == "A simple test plugin"
    assert p['status'] == "Loaded"

def test_plugin_registration(plugin_manager, tmp_path):
    """Test that a plugin can register actions via the context."""
    dummy_plugin_content = """
PLUGIN_NAME = "Action Plugin"

def initialize(context):
    context.add_menu_action(
        path="Test > Action",
        callback=lambda: None,
        text="Test Action",
        icon=None,
        shortcut=None
    )
"""
    source_file = tmp_path / "action_plugin.py"
    source_file.write_text(dummy_plugin_content, encoding="utf-8")
    
    plugin_manager.install_plugin(str(source_file))
    plugin_manager.discover_plugins()
    
    assert len(plugin_manager.menu_actions) == 1
    action = plugin_manager.menu_actions[0]
    assert action['plugin'] == "Action Plugin"
    assert action['text'] == "Test Action"

def test_install_zip(plugin_manager, tmp_path):
    """Test installing a plugin from a ZIP file."""
    # Create a directory structure to zip
    zip_source = tmp_path / "zip_source"
    zip_source.mkdir()
    (zip_source / "__init__.py").write_text("PLUGIN_NAME='Zipped Plugin'", encoding="utf-8")
    
    # Create zip file
    zip_file = tmp_path / "plugin.zip"
    with zipfile.ZipFile(zip_file, 'w') as zf:
        zf.write(zip_source / "__init__.py", "MyPlugin/__init__.py")
        
    # Install
    success, msg = plugin_manager.install_plugin(str(zip_file))
    assert success
    
    # Check extraction
    extracted_dir = tmp_path / "plugins" / "MyPlugin"
    assert extracted_dir.exists()
    assert (extracted_dir / "__init__.py").exists()
    
    # Discover
    plugins = plugin_manager.discover_plugins()
    assert len(plugins) == 1
    assert plugins[0]['name'] == 'Zipped Plugin'

def test_ast_metadata_parsing(plugin_manager, tmp_path):
    """Test the safe metadata extraction using AST."""
    source_file = tmp_path / "ast_test.py"
    source_file.write_text("""
PLUGIN_NAME = "AST Plugin"
PLUGIN_VERSION = (1, 2, 3)
__author__ = "Test Author"
'''
Docstring description.
This should be extracted.
'''
""", encoding="utf-8")
    
    info = plugin_manager.get_plugin_info_safe(str(source_file))
    assert info['name'] == "AST Plugin"
    assert info['version'] == "1.2.3"
    assert info['author'] == "Test Author"
    assert info['description'] == "Docstring description."
