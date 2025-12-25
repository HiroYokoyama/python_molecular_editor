#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import os
import sys
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget, 
    QTableWidgetItem, QHeaderView, QLabel, QFileDialog, QMessageBox, QAbstractItemView
)
from PyQt6.QtCore import Qt, QMimeData, QUrl
from PyQt6.QtGui import QDragEnterEvent, QDropEvent, QDesktopServices

class PluginManagerWindow(QDialog):
    def __init__(self, plugin_manager, parent=None):
        super().__init__(parent)
        self.plugin_manager = plugin_manager
        self.setWindowTitle("Plugin Manager")
        self.resize(800, 500)
        self.setAcceptDrops(True) # Enable drag & drop for the whole window
        
        self.init_ui()
        self.refresh_plugin_list()

    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Header / Instruction
        lbl_info = QLabel("Drag & Drop .py files here to install plugins.")
        lbl_info.setStyleSheet("color: gray; font-style: italic;")
        layout.addWidget(lbl_info)

        # Plugin Table
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(["Status", "Name", "Version", "Author", "Location", "Description"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Interactive) 
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch) # Description stretches
        self.table.setColumnWidth(1, 200) # Make Name column wider
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self.update_button_state)
        self.table.itemDoubleClicked.connect(self.show_plugin_details)
        layout.addWidget(self.table)

        # Buttons
        btn_layout = QHBoxLayout()
        
        btn_reload = QPushButton("Reload Plugins")
        btn_reload.clicked.connect(self.on_reload)
        btn_layout.addWidget(btn_reload)
        
        btn_folder = QPushButton("Open Plugin Folder")
        btn_folder.clicked.connect(self.plugin_manager.open_plugin_folder)
        btn_layout.addWidget(btn_folder)
        
        self.btn_remove = QPushButton("Remove Plugin")
        self.btn_remove.clicked.connect(self.on_remove_plugin)
        self.btn_remove.setEnabled(False)
        btn_layout.addWidget(self.btn_remove)

        btn_explore = QPushButton("Explore Plugins Online")
        btn_explore.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://hiroyokoyama.github.io/moleditpy-plugins/explorer/")))
        btn_layout.addWidget(btn_explore)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_layout.addStretch()
        btn_layout.addWidget(btn_close)
        
        layout.addLayout(btn_layout)

    def refresh_plugin_list(self):
        self.table.setRowCount(0)
        plugins = self.plugin_manager.plugins
        
        self.table.setRowCount(len(plugins))
        for row, p in enumerate(plugins):
            status_item = QTableWidgetItem(str(p.get('status', 'Unknown')))
            status_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.table.setItem(row, 0, status_item)
            self.table.setItem(row, 1, QTableWidgetItem(str(p.get('name', 'Unknown'))))
            self.table.setItem(row, 2, QTableWidgetItem(str(p.get('version', ''))))
            self.table.setItem(row, 3, QTableWidgetItem(str(p.get('author', ''))))
            
            # Location (Relative Path)
            full_path = p.get('filepath', '')
            rel_path = ""
            if full_path:
                try:
                    rel_path = os.path.relpath(full_path, self.plugin_manager.plugin_dir)
                except Exception:
                    rel_path = os.path.basename(full_path)
            self.table.setItem(row, 4, QTableWidgetItem(str(rel_path)))
            
            self.table.setItem(row, 5, QTableWidgetItem(str(p.get('description', ''))))
            
            # Simple color coding for status
            status = str(p.get('status', ''))
            color = None
            if status.startswith("Error"):
                color = Qt.GlobalColor.red
            elif status == "Loaded":
                color = Qt.GlobalColor.darkGreen
            elif status == "No Entry Point":
                color = Qt.GlobalColor.gray
            
            if color:
                self.table.item(row, 0).setForeground(color)

    def update_button_state(self):
        has_selection = (self.table.currentRow() >= 0)
        if hasattr(self, 'btn_remove'):
            self.btn_remove.setEnabled(has_selection)

    def on_reload(self, silent=False):
        # Trigger reload in main manager
        if self.plugin_manager.main_window:
            self.plugin_manager.discover_plugins(self.plugin_manager.main_window)
            self.refresh_plugin_list()
            # Also update main window menu if possible, but that might require a callback or signal
            # For now we assume discover_plugins re-runs autoruns which might duplicate stuff if not careful?
            # Actually discover_plugins clears lists, so re-running is safe logic-wise, 
            # but main_window need to rebuild its menu.
            # We will handle UI rebuild in the main window code by observing or callback.
            
            # For immediate feedback:
            if not silent:
                QMessageBox.information(self, "Reloaded", "Plugins have been reloaded.")
        else:
            self.plugin_manager.discover_plugins()
            self.refresh_plugin_list()

    def on_remove_plugin(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.warning(self, "Warning", "Please select a plugin to remove.")
            return

        # Assuming table row index matches plugins list index (confirmed in refresh_plugin_list)
        if row < len(self.plugin_manager.plugins):
            plugin = self.plugin_manager.plugins[row]
            filepath = plugin.get('filepath')
            
            if filepath and os.path.exists(filepath):
               reply = QMessageBox.question(self, "Remove Plugin", 
                                            f"Are you sure you want to remove '{plugin.get('name', 'Unknown')}'?\n\nFile: {filepath}\nThis cannot be undone.", 
                                            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
               if reply == QMessageBox.StandardButton.Yes:
                   try:
                       os.remove(filepath)
                       self.on_reload(silent=True) # Reload list and plugins
                       QMessageBox.information(self, "Success", f"Removed '{plugin.get('name', 'Unknown')}'.")
                   except Exception as e:
                       QMessageBox.critical(self, "Error", f"Failed to delete file: {e}")
            else:
               QMessageBox.warning(self, "Error", f"Plugin file not found:\n{filepath}")

    def show_plugin_details(self, item):
        row = item.row()
        if row < len(self.plugin_manager.plugins):
            p = self.plugin_manager.plugins[row]
            msg = f"Name: {p.get('name', 'Unknown')}\n" \
                  f"Version: {p.get('version', 'Unknown')}\n" \
                  f"Author: {p.get('author', 'Unknown')}\n" \
                  f"Status: {p.get('status', 'Unknown')}\n" \
                  f"Location: {p.get('filepath', 'Unknown')}\n\n" \
                  f"Description:\n{p.get('description', 'No description available.')}"
            QMessageBox.information(self, "Plugin Details", msg)

    # --- Drag & Drop Support ---
    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        files_installed = []
        errors = []
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if os.path.isfile(file_path) and file_path.endswith('.py'):
                # Extract info and confirm
                info = self.plugin_manager.get_plugin_info_safe(file_path)
                msg = (f"Do you want to install this plugin?\n\n"
                       f"Name: {info['name']}\n"
                       f"Author: {info['author']}\n"
                       f"Version: {info['version']}\n"
                       f"Description: {info['description']}\n\n"
                       f"File: {os.path.basename(file_path)}")
                
                reply = QMessageBox.question(self, "Install Plugin?", msg, 
                                             QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                
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

