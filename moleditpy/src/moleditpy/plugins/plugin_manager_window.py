#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import os
import shutil
from typing import Any, Optional, cast


from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices, QDragEnterEvent, QDropEvent
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)


class PluginManagerWindow(QDialog):
    """Dialog for browsing, installing, and removing plugins."""

    def __init__(self, plugin_manager: Any, parent: Optional[QWidget] = None) -> None:
        """Initialize plugin management dialog window."""
        super().__init__(parent)
        self.btn_remove: Any = None
        self.table: Any = None
        self.search_input: Any = None
        self.plugin_manager = plugin_manager
        self._preferences_applied = False
        self.setWindowTitle("Plugin Manager")
        self.resize(850, 520)
        self.setAcceptDrops(True)  # Enable drag & drop for the whole window

        self.init_ui()
        self.refresh_plugin_list()

    def init_ui(self) -> None:
        """Build plugin manager UI with table, buttons, search box, and drop support."""
        layout = QVBoxLayout(self)

        lbl_info = QLabel(
            "Drag & Drop .py or .zip files to install plugins. Use checkboxes to enable/disable."
        )
        lbl_info.setStyleSheet("color: gray; font-style: italic;")
        layout.addWidget(lbl_info)

        # Search Box
        search_layout = QHBoxLayout()
        search_label = QLabel("Search:")
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText(
            "Search plugins by name, author, location, or description..."
        )
        self.search_input.setClearButtonEnabled(True)
        self.search_input.textChanged.connect(self.filter_plugins)
        search_layout.addWidget(search_label)
        search_layout.addWidget(self.search_input)
        layout.addLayout(search_layout)

        # Plugin Table
        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels(
            [
                "Enabled",
                "Status",
                "Name",
                "Version",
                "Author",
                "Location",
                "Description",
            ]
        )
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive
        )
        self.table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self.table.horizontalHeader().setSectionResizeMode(
            5, QHeaderView.ResizeMode.Interactive
        )
        self.table.horizontalHeader().setSectionResizeMode(
            6, QHeaderView.ResizeMode.Stretch
        )  # Description stretches
        self.table.setColumnWidth(0, 70)  # Enabled checkbox
        self.table.setColumnWidth(1, 100)  # Status
        self.table.setColumnWidth(2, 180)  # Name
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self.update_button_state)
        self.table.itemDoubleClicked.connect(self.show_plugin_details)
        layout.addWidget(self.table)

        # Buttons
        btn_layout = QHBoxLayout()

        btn_reload = QPushButton("Reload Plugins")
        # Not connect(self.on_reload): clicked(bool) would land in *silent*.
        btn_reload.clicked.connect(lambda: self.on_reload())
        btn_layout.addWidget(btn_reload)

        btn_folder = QPushButton("Open Plugin Folder")
        btn_folder.clicked.connect(self.plugin_manager.open_plugin_folder)
        btn_layout.addWidget(btn_folder)

        self.btn_remove = QPushButton("Remove Plugin")
        self.btn_remove.clicked.connect(self.on_remove_plugin)
        self.btn_remove.setEnabled(False)
        btn_layout.addWidget(self.btn_remove)

        btn_explore = QPushButton("Explore Plugins Online")
        btn_explore.clicked.connect(
            lambda: QDesktopServices.openUrl(
                QUrl("https://hiroyokoyama.github.io/moleditpy-plugins/explorer/")
            )
        )
        btn_layout.addWidget(btn_explore)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addStretch()
        btn_layout.addWidget(btn_close)

        layout.addLayout(btn_layout)

    def _create_checkbox_container(
        self, checked: bool, tooltip: str = "Enable or disable this plugin"
    ) -> QWidget:
        """Create a centered QCheckBox container; read it back via _status_checkbox."""
        container = QWidget()
        box_layout = QHBoxLayout(container)
        box_layout.setContentsMargins(0, 0, 0, 0)
        box_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        checkbox = QCheckBox()
        checkbox.setChecked(checked)
        checkbox.setToolTip(tooltip)
        box_layout.addWidget(checkbox)
        return container

    def refresh_plugin_list(self) -> None:
        """Repopulate the plugin table from the current plugin registry."""
        self.table.setRowCount(0)
        plugins = self.plugin_manager.plugins

        self.table.setRowCount(len(plugins))
        for row, p in enumerate(plugins):
            # Column 0: Enabled Checkbox (Clean native check, centered)
            is_enabled = not p.get("disabled", False)
            self.table.setCellWidget(
                row, 0, self._create_checkbox_container(is_enabled)
            )

            # Empty item behind the widget so the column still selects the row.
            self.table.setItem(row, 0, QTableWidgetItem())

            # Column 1: Status (colored text)
            status = str(p.get("status", "Unknown"))
            status_item = QTableWidgetItem(status)
            status_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            if status.startswith("Error"):
                status_item.setForeground(Qt.GlobalColor.red)
            elif status == "Loaded":
                status_item.setForeground(Qt.GlobalColor.darkGreen)
            elif status in ("No Entry Point", "Disabled"):
                status_item.setForeground(Qt.GlobalColor.gray)
            self.table.setItem(row, 1, status_item)

            # Column 2: Name
            self.table.setItem(row, 2, QTableWidgetItem(str(p.get("name", "Unknown"))))
            # Column 3: Version
            self.table.setItem(row, 3, QTableWidgetItem(str(p.get("version", ""))))
            # Column 4: Author
            self.table.setItem(row, 4, QTableWidgetItem(str(p.get("author", ""))))

            # Column 5: Location (Relative Path)
            full_path = p.get("filepath", "")
            rel_path = ""
            if full_path:
                try:
                    rel_path = os.path.relpath(
                        full_path, self.plugin_manager.plugin_dir
                    )
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    rel_path = os.path.basename(full_path)
            self.table.setItem(row, 5, QTableWidgetItem(str(rel_path)))

            # Column 6: Description
            self.table.setItem(row, 6, QTableWidgetItem(str(p.get("description", ""))))

        if self.search_input and self.search_input.text():
            self.filter_plugins(self.search_input.text())

    def filter_plugins(self, text: str) -> None:
        """Filter table rows based on search text."""
        query = text.strip().lower()
        for row in range(self.table.rowCount()):
            if not query:
                self.table.setRowHidden(row, False)
                continue

            # Check status, name, author, location, and description
            status = (
                self.table.item(row, 1).text().lower()
                if self.table.item(row, 1)
                else ""
            )
            name = (
                self.table.item(row, 2).text().lower()
                if self.table.item(row, 2)
                else ""
            )
            author = (
                self.table.item(row, 4).text().lower()
                if self.table.item(row, 4)
                else ""
            )
            location = (
                self.table.item(row, 5).text().lower()
                if self.table.item(row, 5)
                else ""
            )
            description = (
                self.table.item(row, 6).text().lower()
                if self.table.item(row, 6)
                else ""
            )

            match = (
                query in name
                or query in author
                or query in location
                or query in description
                or query in status
            )
            self.table.setRowHidden(row, not match)

        self._drop_hidden_selection()

    def _drop_hidden_selection(self) -> None:
        """Clear the selection when the filter hides the row it points at."""
        row = self.table.currentRow()
        if row >= 0 and self.table.isRowHidden(row):
            self.table.clearSelection()
            self.table.setCurrentCell(-1, -1)
            self.update_button_state()

    def _status_checkbox(self, row: int) -> Optional[QCheckBox]:
        """Return the visible status checkbox for a plugin table row."""
        container = self.table.cellWidget(row, 0)
        if isinstance(container, QCheckBox):
            return container
        if isinstance(container, QWidget):
            checkbox = container.findChild(QCheckBox)
            return checkbox
        return None

    def _save_checkbox_preferences(self) -> set[str]:
        """Collect and persist checkbox choices, hidden rows included."""
        disabled_paths = set()
        for row in range(self.table.rowCount()):
            plugin = self._plugin_for_row(row)
            if plugin is None:
                continue
            status_checkbox = self._status_checkbox(row)
            if status_checkbox is not None and not status_checkbox.isChecked():
                filepath = plugin.get("filepath")
                if filepath:
                    disabled_paths.add(self.plugin_manager.plugin_path_key(filepath))

        self.plugin_manager.save_disabled_plugins(disabled_paths)
        return disabled_paths

    def _apply_plugin_preferences(self) -> None:
        """Persist checkbox choices and reload plugins."""
        if self._preferences_applied:
            return
        self._preferences_applied = True

        self._save_checkbox_preferences()
        if self.plugin_manager.main_window:
            self.plugin_manager.discover_plugins(self.plugin_manager.main_window)
            self.plugin_manager.rebuild_plugin_menus()
        else:
            self.plugin_manager.discover_plugins()

    def done(self, result: int) -> None:
        """Apply plugin preferences for every way the dialog can close.

        Close first: applying them re-runs every plugin's ``initialize()``.
        """
        super().done(result)
        self._apply_plugin_preferences()

    def update_button_state(self) -> None:
        """Enable or disable the Remove button based on table selection."""
        has_selection = self.table.currentRow() >= 0
        if hasattr(self, "btn_remove"):
            self.btn_remove.setEnabled(has_selection)

    def on_reload(self, silent: bool = False) -> None:
        """Reload all plugins from disk, rebuild main-window UI, and refresh the table."""
        self._save_checkbox_preferences()
        if self.plugin_manager.main_window:
            self.plugin_manager.discover_plugins(self.plugin_manager.main_window)
            self.plugin_manager.rebuild_plugin_menus()
            self.refresh_plugin_list()

            if not silent:
                QMessageBox.information(self, "Reloaded", "Plugins have been reloaded.")
        else:
            self.plugin_manager.discover_plugins()
            self.refresh_plugin_list()

    def _plugin_for_row(self, row: int) -> Optional[dict[str, Any]]:
        """Retrieve the plugin dict a table row was built from.

        Row index is the plugin index; filtering hides rows, never reorders.
        """
        if row < 0 or row >= min(
            self.table.rowCount(), len(self.plugin_manager.plugins)
        ):
            return None
        return cast(dict[str, Any], self.plugin_manager.plugins[row])

    def on_remove_plugin(self) -> None:
        """Delete the selected plugin file or folder and reload."""
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.warning(self, "Warning", "Please select a plugin to remove.")
            return

        plugin = self._plugin_for_row(row)
        if plugin is not None:
            filepath = plugin.get("filepath")

            if filepath and os.path.exists(filepath):
                # Check if it is a package plugin (based on __init__.py)
                is_package = os.path.basename(filepath) == "__init__.py"
                target_path = os.path.dirname(filepath) if is_package else filepath

                msg = f"Are you sure you want to remove '{plugin.get('name', 'Unknown')}'?"
                if is_package:
                    msg += f"\n\nThis will delete the entire folder:\n{target_path}"
                else:
                    msg += f"\n\nFile: {filepath}"

                msg += "\nThis cannot be undone."

                reply = QMessageBox.question(
                    self,
                    "Remove Plugin",
                    msg,
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if reply == QMessageBox.StandardButton.Yes:
                    try:
                        if is_package:
                            shutil.rmtree(target_path)
                        else:
                            os.remove(target_path)

                        self.on_reload(silent=True)  # Reload list and plugins
                        QMessageBox.information(
                            self,
                            "Success",
                            f"Removed '{plugin.get('name', 'Unknown')}'.",
                        )
                    except OSError as e:
                        logging.exception("Failed to delete plugin: %s", e)
            else:
                QMessageBox.warning(
                    self, "Error", f"Plugin file not found:\n{filepath}"
                )

    def show_plugin_details(self, item: QTableWidgetItem) -> None:
        """Show a message box with full metadata for the double-clicked plugin."""
        row = item.row()
        plugin = self._plugin_for_row(row)
        if plugin is not None:
            msg = (
                f"Name: {plugin.get('name', 'Unknown')}\n"
                f"Version: {plugin.get('version', 'Unknown')}\n"
                f"Author: {plugin.get('author', 'Unknown')}\n"
                f"Status: {plugin.get('status', 'Unknown')}\n"
                f"Location: {plugin.get('filepath', 'Unknown')}\n\n"
                f"Description:\n{plugin.get('description', 'No description available.')}"
            )
            QMessageBox.information(self, "Plugin Details", msg)

    # --- Drag & Drop Support ---
    def dragEnterEvent(self, event: Optional[QDragEnterEvent]) -> None:
        """Accept drag events carrying file URLs."""
        if event is None:
            return
        if event.mimeData().hasUrls():  # type: ignore[union-attr]
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event: Optional[QDropEvent]) -> None:
        """Install dropped plugin .py files, folders, or zip archives."""
        if event is None:
            return
        files_installed = []
        errors = []
        for url in event.mimeData().urls():  # type: ignore[union-attr]
            file_path = url.toLocalFile()

            is_valid = False
            is_zip = False
            is_folder = False

            if os.path.isfile(file_path):
                # Special handling: If user drops __init__.py,
                # assume they want to install the package (folder)
                if os.path.basename(file_path) == "__init__.py":
                    file_path = os.path.dirname(file_path)
                    is_valid = True
                    is_folder = True
                elif file_path.endswith(".py"):
                    is_valid = True
                elif file_path.endswith(".zip"):
                    is_valid = True
                    is_zip = True

            if os.path.isdir(file_path):
                is_valid = True
                is_folder = True

            if is_valid:
                sha256_value = self.plugin_manager.compute_sha256(file_path)
                # Extract info and confirm
                info = {
                    "name": os.path.basename(file_path),
                    "version": "Unknown",
                    "author": "Unknown",
                    "description": "",
                }

                if is_folder:
                    info["description"] = "Folder Plugin / Category"
                    # Try to parse __init__.py if it exists
                    init_path = os.path.join(file_path, "__init__.py")
                    if os.path.exists(init_path):
                        folder_note = info["description"]
                        info = self.plugin_manager.get_plugin_info_safe(
                            init_path, fallback_name=os.path.basename(file_path)
                        )
                        declared = info["description"]
                        info["description"] = f"{folder_note} (Package: {info['name']})"
                        if declared:
                            info["description"] += f" - {declared}"

                elif is_zip:
                    info["description"] = "ZIP Package Plugin"
                elif file_path.endswith(".py"):
                    info = self.plugin_manager.get_plugin_info_safe(file_path)

                msg = (
                    f"Do you want to install this plugin?\n\n"
                    f"Name: {info['name']}\n"
                    f"Author: {info['author']}\n"
                    f"Version: {info['version']}\n"
                    f"Description: {info['description']}\n\n"
                    f"File: {os.path.basename(file_path)}\n"
                    f"SHA-256: {sha256_value}"
                )

                reply = QMessageBox.question(
                    self,
                    "Install Plugin?",
                    msg,
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )

                if reply == QMessageBox.StandardButton.Yes:
                    success, msg = self.plugin_manager.install_plugin(file_path)
                    if success:
                        files_installed.append(msg)
                    else:
                        errors.append(msg)

        if files_installed or errors:
            self.refresh_plugin_list()
            summary = ""
            if files_installed:
                summary += "Installed:\n" + "\n".join(files_installed) + "\n\n"
            if errors:
                summary += "Errors:\n" + "\n".join(errors)

            QMessageBox.information(self, "Plugin Installation", summary)
