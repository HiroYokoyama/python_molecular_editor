import os
import shutil
import hashlib
from unittest.mock import MagicMock, patch

import pytest
from PyQt6.QtCore import Qt, QUrl, QMimeData
from PyQt6.QtGui import QDragEnterEvent, QDropEvent
from PyQt6.QtWidgets import QMessageBox

from moleditpy.plugins.plugin_manager_window import PluginManagerWindow


@pytest.fixture
def mock_plugin_manager():
    manager = MagicMock()
    manager.plugin_dir = "/fake/plugins"
    manager.plugins = [
        {
            "name": "Test Plugin 1",
            "version": "1.0",
            "author": "Author A",
            "status": "Loaded",
            "filepath": "/fake/plugins/plugin1.py",
            "description": "Desc 1",
        },
        {
            "name": "Test Plugin 2",
            "version": "2.0",
            "author": "Author B",
            "status": "Error",
            "filepath": "/fake/plugins/plugin2.py",
            "description": "Desc 2",
        },
        {
            "name": "Package Plugin",
            "version": "1.1",
            "author": "Author C",
            "status": "No Entry Point",
            "filepath": "/fake/plugins/pkg_plugin/__init__.py",
            "description": "Desc 3",
        },
    ]
    return manager


def test_init_and_refresh(mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Plugin Manager"
    assert window.table.rowCount() == 3

    # Check rows are inserted with color correctly
    assert window.table.item(0, 0).text() == "Loaded"
    assert window.table.item(0, 0).foreground().color() == Qt.GlobalColor.darkGreen
    assert window.table.item(1, 0).text() == "Error"
    assert window.table.item(1, 0).foreground().color() == Qt.GlobalColor.red
    assert window.table.item(2, 0).text() == "No Entry Point"
    assert window.table.item(2, 0).foreground().color() == Qt.GlobalColor.gray

    # Check relative path resolution
    assert window.table.item(0, 4).text() == "plugin1.py"


def test_refresh_relative_path_error(mock_plugin_manager, qtbot):
    # If relpath throws exception it shows basename
    mock_plugin_manager.plugins = [{"filepath": "plugin3.py"}]
    with patch("os.path.relpath", side_effect=ValueError):
        window = PluginManagerWindow(mock_plugin_manager)
        qtbot.addWidget(window)
        assert window.table.item(0, 4).text() == "plugin3.py"


def test_update_button_state(mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    assert not window.btn_remove.isEnabled()
    window.table.selectRow(0)
    window.table.itemSelectionChanged.emit()
    assert window.btn_remove.isEnabled()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_on_reload_main_window_present(mock_info, mock_plugin_manager, qtbot):
    mock_plugin_manager.main_window = MagicMock()
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.on_reload()
    mock_plugin_manager.discover_plugins.assert_called_with(
        mock_plugin_manager.main_window
    )
    mock_info.assert_called_once()

    # test silent reload
    mock_info.reset_mock()
    window.on_reload(silent=True)
    mock_info.assert_not_called()


def test_on_reload_no_main_window(mock_plugin_manager, qtbot):
    mock_plugin_manager.main_window = None
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    with patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information"):
        window.on_reload()
    mock_plugin_manager.discover_plugins.assert_called_with()


@patch("moleditpy.plugins.plugin_manager_window.QDesktopServices.openUrl")
def test_explore_plugins_online(mock_open_url, mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    children = window.findChildren(type(window.btn_remove))
    btn = [b for b in children if "Explore" in b.text()][0]
    btn.click()
    mock_open_url.assert_called_once()
    assert (
        "https://hiroyokoyama.github.io/moleditpy-plugins/explorer/"
        in mock_open_url.call_args[0][0].url()
    )


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.warning")
def test_on_remove_plugin_no_selection(mock_warn, mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.on_remove_plugin()
    mock_warn.assert_called_with(window, "Warning", "Please select a plugin to remove.")


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
@patch("os.path.exists", return_value=True)
@patch("os.remove")
def test_on_remove_plugin_single_file(
    mock_remove, mock_exists, mock_info, mock_question, mock_plugin_manager, qtbot
):
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(0)
    window.on_remove_plugin()

    mock_remove.assert_called_with("/fake/plugins/plugin1.py")
    mock_info.assert_called()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
@patch("os.path.exists", return_value=True)
@patch("shutil.rmtree")
def test_on_remove_plugin_package(
    mock_rmtree, mock_exists, mock_info, mock_question, mock_plugin_manager, qtbot
):
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(2)
    window.on_remove_plugin()

    mock_rmtree.assert_called_with("/fake/plugins/pkg_plugin")
    mock_info.assert_called()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.critical")
@patch("os.path.exists", return_value=True)
@patch("os.remove", side_effect=RuntimeError("Remove error"))
def test_on_remove_plugin_error(
    mock_remove, mock_exists, mock_critical, mock_question, mock_plugin_manager, qtbot
):
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(0)
    window.on_remove_plugin()

    mock_critical.assert_called_with(
        window, "Error", "Failed to delete plugin: Remove error"
    )


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.warning")
@patch("os.path.exists", return_value=False)
def test_on_remove_plugin_not_exists(
    mock_exists, mock_warn, mock_plugin_manager, qtbot
):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(0)
    window.on_remove_plugin()

    # Needs to match the call
    mock_warn.assert_called()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_show_plugin_details(mock_info, mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    item = window.table.item(0, 0)
    window.show_plugin_details(item)
    mock_info.assert_called_once()
    assert "Test Plugin 1" in mock_info.call_args[0][2]


def test_drag_enter_event(mock_plugin_manager, qtbot):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    event = MagicMock(spec=QDragEnterEvent)
    mime_data = MagicMock(spec=QMimeData)
    mime_data.hasUrls.return_value = True
    event.mimeData.return_value = mime_data

    window.dragEnterEvent(event)
    event.accept.assert_called_once()

    mime_data.hasUrls.return_value = False
    event.reset_mock()
    window.dragEnterEvent(event)
    event.ignore.assert_called_once()


@patch("os.path.isfile")
@patch("os.path.isdir")
@patch("os.path.exists")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_drop_event_valid_files(
    mock_info,
    mock_question,
    mock_exists,
    mock_isdir,
    mock_isfile,
    mock_plugin_manager,
    qtbot,
):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager._compute_sha256 = MagicMock(return_value="abc")

    mock_plugin_manager.get_plugin_info_safe.return_value = {
        "name": "Dropped",
        "version": "1.0",
        "author": "me",
        "description": "desc",
    }
    mock_plugin_manager.install_plugin.return_value = (True, "Installed fine")

    mock_question.return_value = QMessageBox.StandardButton.Yes

    # Make dropping a .py file valid
    mock_isfile.return_value = True
    mock_isdir.return_value = False

    event = MagicMock(spec=QDropEvent)
    mime_data = MagicMock(spec=QMimeData)
    url_mock = MagicMock()
    url_mock.toLocalFile.return_value = "/some/file.py"
    mime_data.urls.return_value = [url_mock]
    event.mimeData.return_value = mime_data

    window.dropEvent(event)

    mock_plugin_manager.install_plugin.assert_called_once_with("/some/file.py")
    mock_info.assert_called_once()
    assert "Installed fine" in mock_info.call_args[0][2]


@patch("os.path.isfile")
@patch("os.path.isdir")
@patch("os.path.exists")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_drop_event_init_py_package(
    mock_info,
    mock_question,
    mock_exists,
    mock_isdir,
    mock_isfile,
    mock_plugin_manager,
    qtbot,
):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager._compute_sha256 = MagicMock(return_value="abc")

    mock_plugin_manager.get_plugin_info_safe.return_value = {
        "name": "Pkg",
        "version": "1.0",
        "author": "me",
        "description": "",
    }
    mock_plugin_manager.install_plugin.return_value = (False, "Error inst")

    mock_question.return_value = QMessageBox.StandardButton.Yes

    # __init__.py converts to folder
    def isfile_side(x):
        return True

    def isdir_side(x):
        return False

    def exists_side(x):
        return True

    mock_isfile.side_effect = isfile_side
    mock_isdir.side_effect = isdir_side
    mock_exists.side_effect = exists_side

    event = MagicMock(spec=QDropEvent)
    mime_data = MagicMock(spec=QMimeData)
    url_mock = MagicMock()
    url_mock.toLocalFile.return_value = "/some/folder/__init__.py"
    mime_data.urls.return_value = [url_mock]
    event.mimeData.return_value = mime_data

    window.dropEvent(event)

    mock_plugin_manager.install_plugin.assert_called_once_with("/some/folder")
    assert "Error inst" in mock_info.call_args[0][2]


@patch("os.path.isfile")
@patch("os.path.isdir")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_drop_event_zip_file(
    mock_info, mock_question, mock_isdir, mock_isfile, mock_plugin_manager, qtbot
):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager._compute_sha256 = MagicMock(return_value="abc")
    mock_plugin_manager.install_plugin.return_value = (True, "")
    mock_question.return_value = QMessageBox.StandardButton.Yes

    mock_isfile.return_value = True
    mock_isdir.return_value = False

    event = MagicMock(spec=QDropEvent)
    url_mock = MagicMock()
    url_mock.toLocalFile.return_value = "/some/file.zip"
    event.mimeData.return_value = MagicMock(urls=lambda: [url_mock])

    window.dropEvent(event)
    mock_plugin_manager.install_plugin.assert_called_with("/some/file.zip")


@patch("os.path.isdir")
@patch("os.path.isfile")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_drop_event_pure_folder(
    mock_info, mock_question, mock_isfile, mock_isdir, mock_plugin_manager, qtbot
):
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager._compute_sha256 = MagicMock(return_value="abc")
    mock_plugin_manager.install_plugin.return_value = (True, "Installed folder")
    mock_question.return_value = QMessageBox.StandardButton.Yes

    mock_isfile.return_value = False
    mock_isdir.return_value = True

    event = MagicMock(spec=QDropEvent)
    url_mock = MagicMock()
    url_mock.toLocalFile.return_value = "/some/plugin_folder"
    event.mimeData.return_value = MagicMock(urls=lambda: [url_mock])

    window.dropEvent(event)
    mock_plugin_manager.install_plugin.assert_called_with("/some/plugin_folder")
