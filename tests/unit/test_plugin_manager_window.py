"""Unit tests for PluginManagerWindow drag-and-drop and UI behavior."""

import logging
from unittest.mock import MagicMock, patch

import pytest
from PyQt6.QtCore import Qt, QMimeData
from PyQt6.QtGui import QDragEnterEvent, QDropEvent
from PyQt6.QtWidgets import QCheckBox, QMessageBox, QWidget

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
    """PluginManagerWindow initialises with correct title, row count, and status colours."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    assert window.windowTitle() == "Plugin Manager"
    assert window.table.rowCount() == 3
    assert window.table.columnCount() == 7

    # Check status column (index 1) colors
    assert window.table.item(0, 1).text() == "Loaded"
    assert window.table.item(0, 1).foreground().color() == Qt.GlobalColor.darkGreen
    assert window.table.item(1, 1).text() == "Error"
    assert window.table.item(1, 1).foreground().color() == Qt.GlobalColor.red
    assert window.table.item(2, 1).text() == "No Entry Point"
    assert window.table.item(2, 1).foreground().color() == Qt.GlobalColor.gray

    # Check relative path resolution (column 5)
    assert window.table.item(0, 5).text() == "plugin1.py"
    # Column 0 is the Enabled checkbox
    checkbox = window._status_checkbox(0)
    assert isinstance(checkbox, QCheckBox)
    assert checkbox.isChecked()


def test_status_checkbox_is_horizontal_cell_widget(mock_plugin_manager, qtbot):
    """The status toggle is a clean centered QCheckBox in the Enabled column."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    checkbox = window._status_checkbox(0)
    assert isinstance(checkbox, QCheckBox)
    assert checkbox.text() == ""
    assert checkbox.toolTip() == "Enable or disable this plugin"


def test_refresh_relative_path_error(mock_plugin_manager, qtbot):
    """When os.path.relpath raises ValueError, the filepath column falls back to basename."""
    # If relpath throws exception it shows basename
    mock_plugin_manager.plugins = [{"filepath": "plugin3.py"}]
    with patch("os.path.relpath", side_effect=ValueError):
        window = PluginManagerWindow(mock_plugin_manager)
        qtbot.addWidget(window)
        assert window.table.item(0, 5).text() == "plugin3.py"


def test_update_button_state(mock_plugin_manager, qtbot):
    """Remove button is disabled initially and enabled after a row is selected."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    assert not window.btn_remove.isEnabled()
    window.table.selectRow(0)
    window.table.itemSelectionChanged.emit()
    assert window.btn_remove.isEnabled()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_on_reload_main_window_present(mock_info, mock_plugin_manager, qtbot):
    """on_reload discovers plugins with the main window and shows an info message."""
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
    """on_reload calls discover_plugins without arguments when main_window is None."""
    mock_plugin_manager.main_window = None
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    with patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information"):
        window.on_reload()
    mock_plugin_manager.discover_plugins.assert_called_with()


@patch("moleditpy.plugins.plugin_manager_window.QDesktopServices.openUrl")
def test_explore_plugins_online(mock_open_url, mock_plugin_manager, qtbot):
    """The Explore button opens the plugin explorer URL in the system browser."""
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
    """on_remove_plugin shows a warning when no plugin row is selected."""
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
    """on_remove_plugin removes a single-file plugin by calling os.remove."""
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
    """on_remove_plugin removes a package plugin by calling shutil.rmtree on its directory."""
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(2)
    window.on_remove_plugin()

    mock_rmtree.assert_called_with("/fake/plugins/pkg_plugin")
    mock_info.assert_called()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("os.path.exists", return_value=True)
@patch("os.remove", side_effect=PermissionError("Remove error"))
def test_on_remove_plugin_error(
    mock_remove, mock_exists, mock_question, mock_plugin_manager, qtbot, caplog
):
    """on_remove_plugin logs the error when os.remove raises OSError."""
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(0)
    with caplog.at_level(logging.ERROR):
        window.on_remove_plugin()

    assert "Failed to delete plugin: Remove error" in caplog.text


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.warning")
@patch("os.path.exists", return_value=False)
def test_on_remove_plugin_not_exists(
    mock_exists, mock_warn, mock_plugin_manager, qtbot
):
    """on_remove_plugin shows a warning when the plugin file does not exist."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.table.selectRow(0)
    window.on_remove_plugin()

    # Needs to match the call
    mock_warn.assert_called()


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_show_plugin_details(mock_info, mock_plugin_manager, qtbot):
    """show_plugin_details displays an information dialog containing the plugin name."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    item = window.table.item(0, 0)
    window.show_plugin_details(item)
    mock_info.assert_called_once()
    assert "Test Plugin 1" in mock_info.call_args[0][2]


def test_drag_enter_event(mock_plugin_manager, qtbot):
    """dragEnterEvent accepts URLs and ignores non-URL mime data."""
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
    """Dropping a .py file installs it and shows a success info message."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager.compute_sha256 = MagicMock(return_value="abc")

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
    """Dropping an __init__.py file installs its parent folder as a package."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager.compute_sha256 = MagicMock(return_value="abc")

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
    """Dropping a .zip file installs it via install_plugin."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager.compute_sha256 = MagicMock(return_value="abc")
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
    """Dropping a directory installs it directly via install_plugin."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    mock_plugin_manager.compute_sha256 = MagicMock(return_value="abc")
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


def test_close_persists_disabled_paths_and_reloads_once(mock_plugin_manager, qtbot):
    """Persists disabled plugin paths and reloads only once when closing."""
    mock_plugin_manager.main_window = None
    mock_plugin_manager.plugin_path_key.side_effect = lambda filepath: filepath.rsplit(
        "/", 1
    )[-1]
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window._status_checkbox(0).setChecked(False)
    window.done(0)
    window.done(0)

    mock_plugin_manager.save_disabled_plugins.assert_called_once_with({"plugin1.py"})
    mock_plugin_manager.discover_plugins.assert_called_once_with()
    mock_plugin_manager.rebuild_plugin_menus.assert_not_called()


def test_close_with_main_window_reloads_and_rebuilds_menus(mock_plugin_manager, qtbot):
    """Persists preferences and rebuilds plugin menus when closing."""
    mock_main_window = MagicMock()
    mock_plugin_manager.main_window = mock_main_window
    mock_plugin_manager.plugin_path_key.side_effect = lambda filepath: filepath.rsplit(
        "/", 1
    )[-1]
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window._status_checkbox(0).setChecked(False)
    window.done(0)

    mock_plugin_manager.save_disabled_plugins.assert_called_once_with({"plugin1.py"})
    mock_plugin_manager.discover_plugins.assert_called_once_with(mock_main_window)
    mock_plugin_manager.rebuild_plugin_menus.assert_called_once_with()


def test_close_hides_the_dialog_before_rediscovering(mock_plugin_manager, qtbot):
    """The dialog must be closed before discovery re-runs plugin initialize().

    Discovery executes every plugin's initialize(); one that raises its own
    dialog would sit behind this window while it is still modal and visible.
    """
    mock_plugin_manager.main_window = MagicMock()
    visible_during_discovery = []
    mock_plugin_manager.discover_plugins.side_effect = lambda *a, **k: (
        visible_during_discovery.append(window.isVisible())
    )
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)
    window.show()

    window.done(0)

    assert visible_during_discovery == [False]


def test_on_reload_persists_checkbox_changes(mock_plugin_manager, qtbot):
    """Saves checkbox changes before reloading plugins without a main window."""
    mock_plugin_manager.main_window = None
    mock_plugin_manager.plugin_path_key.side_effect = lambda filepath: filepath.rsplit(
        "/", 1
    )[-1]
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window._status_checkbox(0).setChecked(False)
    window.on_reload(silent=True)

    mock_plugin_manager.save_disabled_plugins.assert_called_with({"plugin1.py"})
    mock_plugin_manager.discover_plugins.assert_called_with()


def test_on_reload_with_main_window_persists_and_rebuilds_menus(
    mock_plugin_manager, qtbot
):
    """Saves checkbox changes and rebuilds menus during a main-window reload."""
    mock_main_window = MagicMock()
    mock_plugin_manager.main_window = mock_main_window
    mock_plugin_manager.plugin_path_key.side_effect = lambda filepath: filepath.rsplit(
        "/", 1
    )[-1]
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window._status_checkbox(0).setChecked(False)
    window.on_reload(silent=True)

    mock_plugin_manager.save_disabled_plugins.assert_called_with({"plugin1.py"})
    mock_plugin_manager.discover_plugins.assert_called_with(mock_main_window)
    mock_plugin_manager.rebuild_plugin_menus.assert_called_with()


def test_status_preserves_status_colors(mock_plugin_manager, qtbot):
    """Applies status foreground colors to the status column item."""
    mock_plugin_manager.plugins[1]["status"] = "Error"
    mock_plugin_manager.plugins[2]["status"] = "Disabled"
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    assert window.table.item(0, 1).foreground().color() == Qt.GlobalColor.darkGreen
    assert window.table.item(1, 1).foreground().color() == Qt.GlobalColor.red
    assert window.table.item(2, 1).foreground().color() == Qt.GlobalColor.gray


def test_search_filter_by_name(mock_plugin_manager, qtbot):
    """Typing in search box filters plugins by name."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Initially all 3 rows visible
    assert not window.table.isRowHidden(0)
    assert not window.table.isRowHidden(1)
    assert not window.table.isRowHidden(2)

    # Search for "Plugin 2"
    window.search_input.setText("Plugin 2")
    assert window.table.isRowHidden(0)
    assert not window.table.isRowHidden(1)
    assert window.table.isRowHidden(2)

    # Clear search
    window.search_input.setText("")
    assert not window.table.isRowHidden(0)
    assert not window.table.isRowHidden(1)
    assert not window.table.isRowHidden(2)


def test_search_filter_by_author_and_description(mock_plugin_manager, qtbot):
    """Search box filters plugins by author or description."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Filter by Author C
    window.search_input.setText("Author C")
    assert window.table.isRowHidden(0)
    assert window.table.isRowHidden(1)
    assert not window.table.isRowHidden(2)

    # Filter by description substring
    window.search_input.setText("Desc 1")
    assert not window.table.isRowHidden(0)
    assert window.table.isRowHidden(1)
    assert window.table.isRowHidden(2)


def test_search_filter_by_location_and_status(mock_plugin_manager, qtbot):
    """Search box filters plugins by location relative path or status."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Filter by status "Error"
    window.search_input.setText("Error")
    assert window.table.isRowHidden(0)
    assert not window.table.isRowHidden(1)
    assert window.table.isRowHidden(2)

    # Filter by relative path part
    window.search_input.setText("pkg_plugin")
    assert window.table.isRowHidden(0)
    assert window.table.isRowHidden(1)
    assert not window.table.isRowHidden(2)


def test_checkbox_toggle_and_save_with_search_filter(mock_plugin_manager, qtbot):
    """Disabling a plugin while search filter is active persists correctly."""
    mock_plugin_manager.main_window = None
    mock_plugin_manager.plugin_path_key.side_effect = lambda filepath: filepath.rsplit(
        "/", 1
    )[-1]
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Filter to only show Package Plugin (row 2)
    window.search_input.setText("Package")
    assert window.table.isRowHidden(0)
    assert window.table.isRowHidden(1)
    assert not window.table.isRowHidden(2)

    # Uncheck row 2's enabled checkbox
    checkbox = window._status_checkbox(2)
    assert checkbox is not None
    checkbox.setChecked(False)

    # Done saves the disabled preference for the right plugin
    window.done(0)
    mock_plugin_manager.save_disabled_plugins.assert_called_once_with({"__init__.py"})


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
def test_show_plugin_details_with_search_filter(mock_info, mock_plugin_manager, qtbot):
    """Double-clicking a filtered row shows correct plugin metadata."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    window.search_input.setText("Plugin 2")
    item = window.table.item(1, 0)
    window.show_plugin_details(item)

    mock_info.assert_called_once()
    assert "Test Plugin 2" in mock_info.call_args[0][2]
    assert "Author B" in mock_info.call_args[0][2]


@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.question")
@patch("moleditpy.plugins.plugin_manager_window.QMessageBox.information")
@patch("os.path.exists", return_value=True)
@patch("os.remove")
def test_on_remove_plugin_with_search_filter(
    mock_remove, mock_exists, mock_info, mock_question, mock_plugin_manager, qtbot
):
    """Removing a plugin when filtered resolves the correct file path."""
    mock_question.return_value = QMessageBox.StandardButton.Yes
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Filter to show only row 1
    window.search_input.setText("Plugin 2")
    window.table.selectRow(1)
    window.on_remove_plugin()

    mock_remove.assert_called_with("/fake/plugins/plugin2.py")
    mock_info.assert_called()


def test_status_checkbox_direct_and_none(mock_plugin_manager, qtbot):
    """_status_checkbox handles direct QCheckBox and missing/empty widgets."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # When a cell contains directly a QCheckBox (not wrapped in a QWidget container)
    direct_cb = QCheckBox()
    window.table.setCellWidget(0, 0, direct_cb)
    assert window._status_checkbox(0) is direct_cb

    # When a cell contains a QWidget with no QCheckBox child
    empty_widget = QWidget()
    window.table.setCellWidget(0, 0, empty_widget)
    assert window._status_checkbox(0) is None

    # When a cell widget is completely None
    window.table.setCellWidget(0, 0, None)
    assert window._status_checkbox(0) is None


def test_plugin_for_row_boundary_conditions(mock_plugin_manager, qtbot):
    """_plugin_for_row returns None for negative or out-of-range rows and index mappings."""
    window = PluginManagerWindow(mock_plugin_manager)
    qtbot.addWidget(window)

    # Negative row or beyond rowCount
    assert window._plugin_for_row(-1) is None
    assert window._plugin_for_row(999) is None

    # Item with data beyond plugin list length
    item = window.table.item(0, 0)
    item.setData(Qt.ItemDataRole.UserRole, 9999)
    assert window._plugin_for_row(0) is None
