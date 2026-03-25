#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
main_window_ui_manager.py
Mixin class separated from main_window.py
"""

import vtk

# PyQt6 Modules
from PyQt6.QtCore import QEvent, Qt, QTimer
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QGraphicsView,
    QMainWindow,
    QMessageBox,
)

try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .custom_interactor_style import CustomInteractorStyle
except ImportError:
    # Fallback to absolute imports for script-style execution
    from modules.custom_interactor_style import CustomInteractorStyle


# --- Classes ---
class MainWindowUiManager:
    """Mixin class separated from main_window.py."""

    def update_status_bar(self, message):
        """Update status bar with worker messages."""
        self.statusBar().showMessage(message)

    def set_mode(self, mode_str):
        prev_mode = getattr(self.scene, "mode", None)
        self.scene.mode = mode_str
        self.view_2d.setMouseTracking(True)
        # Clear ghost when leaving template mode
        if (
            prev_mode
            and prev_mode.startswith("template")
            and not mode_str.startswith("template")
        ):
            self.scene.clear_template_preview()
        elif not mode_str.startswith("template"):
            self.scene.template_preview.hide()

        # Set cursor shape
        if mode_str == "select":
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)
        elif mode_str.startswith(("atom", "bond", "template")):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        elif mode_str.startswith(("charge", "radical")):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)

        if mode_str.startswith("atom"):
            self.scene.current_atom_symbol = mode_str.split("_")[1]
            self.statusBar().showMessage(
                f"Mode: Draw Atom ({self.scene.current_atom_symbol})"
            )
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True)
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0
        elif mode_str.startswith("bond"):
            self.scene.current_atom_symbol = "C"
            parts = mode_str.split("_")
            self.scene.bond_order = int(parts[1])
            self.scene.bond_stereo = int(parts[2]) if len(parts) > 2 else 0
            stereo_text = {0: "", 1: " (Wedge)", 2: " (Dash)"}.get(
                self.scene.bond_stereo, ""
            )
            self.statusBar().showMessage(
                f"Mode: Draw Bond (Order: {self.scene.bond_order}{stereo_text})"
            )
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True)
        elif mode_str.startswith("template"):
            if mode_str.startswith("template_user"):
                # User template mode
                template_name = mode_str.replace("template_user_", "")
                self.statusBar().showMessage(f"Mode: User Template ({template_name})")
            else:
                # Built-in template mode
                self.statusBar().showMessage(
                    f"Mode: {mode_str.split('_')[1].capitalize()} Template"
                )
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "charge_plus":
            self.statusBar().showMessage("Mode: Increase Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "charge_minus":
            self.statusBar().showMessage("Mode: Decrease Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "radical":
            self.statusBar().showMessage("Mode: Toggle Radical (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)

        else:  # Select mode
            self.statusBar().showMessage("Mode: Select")
            self.view_2d.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0

    def set_mode_and_update_toolbar(self, mode_str):
        self.set_mode(mode_str)
        # Map QAction to QToolButton
        toolbar = getattr(self, "toolbar", None)
        action_to_button = {}
        if toolbar:
            for key, action in self.mode_actions.items():
                btn = toolbar.widgetForAction(action)
                if btn:
                    action_to_button[action] = btn

        # Reset all mode buttons
        for key, action in self.mode_actions.items():
            action.setChecked(False)
            btn = action_to_button.get(action)
            if btn:
                btn.setStyleSheet("")

        # Apply style to template buttons
        if mode_str in self.mode_actions:
            action = self.mode_actions[mode_str]
            action.setChecked(True)
            btn = action_to_button.get(action)
            if btn:
                # Blue for templates, clear otherwise
                if mode_str.startswith("template"):
                    btn.setStyleSheet("background-color: #2196F3; color: white;")
                else:
                    btn.setStyleSheet("")

    def activate_select_mode(self):
        self.set_mode("select")
        if "select" in self.mode_actions:
            self.mode_actions["select"].setChecked(True)

    def eventFilter(self, obj, event):
        if obj is self.plotter and event.type() == QEvent.Type.MouseButtonPress:
            self.view_2d.setFocus()
        return super().eventFilter(obj, event)

    def closeEvent(self, event):        # 1. Persist settings
        try:
            modified = getattr(self, "settings_dirty", False) or (hasattr(self, "settings") and hasattr(self, "initial_settings") and self.settings != self.initial_settings)
            if modified and hasattr(self, "save_settings"):
                self.save_settings()
                self.settings_dirty = False
        except (AttributeError, RuntimeError, TypeError, ValueError, OSError):
            # Suppress non-critical settings persistence errors during application shutdown
            pass

        # 2. Handle unsaved changes
        if getattr(self, "has_unsaved_changes", False):
            reply = QMessageBox.question(
                self,
                "Unsaved Changes",
                "You have unsaved changes. Do you want to save them?",
                QMessageBox.StandardButton.Yes
                | QMessageBox.StandardButton.No
                | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Yes,
            )

            if reply == QMessageBox.StandardButton.Yes:
                if hasattr(self, "save_project"):
                    self.save_project()
                if getattr(self, "has_unsaved_changes", False):
                    event.ignore()
                    return
            elif reply == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return

        # 3. Gracefully close child windows and cleanup threads
        try:
            # Close dialogs
            for widget in QApplication.topLevelWidgets():
                if widget != self and isinstance(widget, (QDialog, QMainWindow)):
                    try:
                        widget.close()
                    except (RuntimeError, TypeError):
                        # Suppress errors if a thread is already terminated or non-responsive during bulk teardown.
                        pass
            # Stop calculation threads
            active_threads = list(getattr(self, "_active_calc_threads", []) or [])
            for thr in active_threads:
                try:
                    if hasattr(thr, "quit"): thr.quit()
                    if hasattr(thr, "wait"): thr.wait(200)
                except (RuntimeError, TypeError):
                    # Suppress errors if a thread is already terminated or non-responsive during bulk teardown.
                    pass
        except (AttributeError, RuntimeError, TypeError, ValueError):
            # Suppress non-critical thread/widget cleanup errors during application shutdown
            pass


        event.accept()

    def toggle_3d_edit_mode(self, checked):
        """Toggle 3D Drag mode."""
        if checked:
            # Disable measurement mode when 3D Drag is on
            if self.measurement_mode:
                self.measurement_action.setChecked(False)
                self.toggle_measurement_mode(False)

        self.is_3d_edit_mode = checked
        if checked:
            self.statusBar().showMessage("3D Drag Mode: ON.")
        else:
            self.statusBar().showMessage("3D Drag Mode: OFF.")
        self.view_2d.setFocus()

    def _setup_3d_picker(self):
        self.plotter.picker = vtk.vtkCellPicker()
        self.plotter.picker.SetTolerance(0.025)

        # Create CustomInteractorStyle
        style = CustomInteractorStyle(self)

        # Set interactor style
        self.plotter.interactor.SetInteractorStyle(style)
        self.plotter.interactor.Initialize()
    def dragEnterEvent(self, event):
        """Handle drag enter event."""
        self.handle_drag_enter_event(event)

    def handle_drag_enter_event(self, event):
        """Internal handler for drag enter event (bypasses PyQt type checks in tests)."""
        if not event.mimeData().hasUrls():
            event.ignore()
            return

        for url in event.mimeData().urls():
            if not url.isLocalFile():
                continue
            
            file_lower = url.toLocalFile().lower()
            # 1. Built-in extensions
            if file_lower.endswith((".pmeraw", ".pmeprj", ".mol", ".sdf", ".xyz")):
                event.acceptProposedAction()
                return

            # 2. Plugin drop handlers (accept if any handlers exist)
            plugin_mgr = getattr(self, "plugin_manager", None)
            if plugin_mgr and getattr(plugin_mgr, "drop_handlers", []):
                event.acceptProposedAction()
                return

        event.ignore()

    def dropEvent(self, event):
        """Handle file drop event."""
        self.handle_drop_event(event)

    def handle_drop_event(self, event):
        """Internal handler for file drop event (bypasses PyQt type checks in tests)."""
        urls = event.mimeData().urls()
        file_path = next((u.toLocalFile() for u in urls if u.isLocalFile()), None)
        if not file_path:
            event.ignore()
            return

        # 1. Plugin Handlers
        plugin_mgr = getattr(self, "plugin_manager", None)
        for handler_def in getattr(plugin_mgr, "drop_handlers", []):
            try:
                if handler_def["callback"](file_path):
                    event.acceptProposedAction()
                    return
            except (AttributeError, RuntimeError, TypeError, ValueError) as e:
                # Log but suppress plugin-specific drop handler errors to prevent app-wide crash
                print(f"Error in plugin drop handler: {e}")

        # 2. Built-in Handlers
        file_lower = file_path.lower()
        if file_lower.endswith((".pmeraw", ".pmeprj")):
            self.open_project_file(file_path=file_path)
            QTimer.singleShot(100, self.fit_to_view)
            event.acceptProposedAction()
        elif file_lower.endswith((".mol", ".sdf")):
            # Check if dropped on 3D viewer
            plotter_widget = self.splitter.widget(1)
            drag_point = event.position().toPoint()
            if plotter_widget and plotter_widget.geometry().contains(drag_point):
                self.load_mol_file_for_3d_viewing(file_path=file_path)
            elif hasattr(self, "load_mol_file"):
                self.load_mol_file(file_path=file_path)
            else:
                self.statusBar().showMessage("MOL file import not implemented for 2D editor.")
            
            QTimer.singleShot(100, self.fit_to_view)
            event.acceptProposedAction()
        elif file_lower.endswith(".xyz"):
            self.load_xyz_for_3d_viewing(file_path=file_path)
            QTimer.singleShot(100, self.fit_to_view)
            event.acceptProposedAction()
        else:
            self.statusBar().showMessage(f"Unsupported file type: {file_path}")
            event.ignore()


    def _enable_3d_edit_actions(self, enabled=True):
        """Enable/disable 3D edit actions."""
        actions = [
            "translation_action",
            "move_group_action",
            "alignplane_xy_action",
            "alignplane_xz_action",
            "alignplane_yz_action",
            "align_x_action",
            "align_y_action",
            "align_z_action",
            "bond_length_action",
            "angle_action",
            "dihedral_action",
            "mirror_action",
            "planarize_action",
            "constrained_opt_action",
        ]

        # Also enable/disable menus
        menus = ["align_menu"]

        for action_name in actions:
            if hasattr(self, action_name):
                getattr(self, action_name).setEnabled(enabled)

        for menu_name in menus:
            if hasattr(self, menu_name):
                getattr(self, menu_name).setEnabled(enabled)
    def _enable_3d_features(self, enabled=True):
        """Enable/disable 3D features."""
        # Basic 3D features
        basic_3d_actions = ["optimize_3d_button", "export_button", "analysis_action"]

        for action_name in basic_3d_actions:
            obj = getattr(self, action_name, None)
            if obj is None:
                continue

            try:
                if action_name == "optimize_3d_button":
                    # Optimization is disabled for XYZ-derived or failed-chem-check molecules
                    is_xyz = getattr(self, "is_xyz_derived", False)
                    chem_failed = getattr(self, "chem_check_tried", False) and getattr(self, "chem_check_failed", False)
                    
                    can_optimize = enabled and not (is_xyz or chem_failed)
                    obj.setEnabled(can_optimize)
                else:
                    obj.setEnabled(enabled)
            except (AttributeError, RuntimeError, TypeError, ValueError):
                # Suppress non-critical 3D feature state update errors if widgets are not fully initialized
                pass

        # Always enable these core 3D interactors
        for core_act in ["measurement_action", "edit_3d_action"]:
            obj = getattr(self, core_act, None)
            if obj:
                obj.setEnabled(True)

        self._enable_3d_edit_actions(enabled)


    def _enter_3d_viewer_ui_mode(self):
        """Set UI mode to 3D viewer."""
        self.is_2d_editable = False
        self.cleanup_button.setEnabled(False)
        self.convert_button.setEnabled(False)
        for action in self.tool_group.actions():
            action.setEnabled(False)
        if hasattr(self, "other_atom_action"):
            self.other_atom_action.setEnabled(False)

        self.minimize_2d_panel()

        # Collectively enable 3D-related features
        self._enable_3d_features(True)

    def restore_ui_for_editing(self):
        """Enables all 2D editing UI elements."""
        self.is_2d_editable = True
        self.restore_2d_panel()
        self.cleanup_button.setEnabled(True)
        self.convert_button.setEnabled(True)

        for action in self.tool_group.actions():
            action.setEnabled(True)

        if hasattr(self, "other_atom_action"):
            self.other_atom_action.setEnabled(True)

        # Collectively disable 3D edit functions when returning to 2D mode
        self._enable_3d_edit_actions(False)

    def minimize_2d_panel(self):
        """Minimize (hide) 2D panel."""
        sizes = self.splitter.sizes()
        # Only if not already minimized
        if sizes[0] > 0:
            total_width = sum(sizes)
            self.splitter.setSizes([0, total_width])

    def restore_2d_panel(self):
        """Restore 2D panel."""
        sizes = self.splitter.sizes()

        # Check sizes list before access
        if sizes and sizes[0] == 0:
            self.splitter.setSizes([600, 600])

    def set_panel_layout(self, left_percent, right_percent):
        """Set panel layout ratio."""
        if left_percent + right_percent != 100:
            return

        total_width = self.splitter.width()
        if total_width <= 0:
            total_width = 1200  # Default width

        left_width = int(total_width * left_percent / 100)
        right_width = int(total_width * right_percent / 100)

        self.splitter.setSizes([left_width, right_width])

        # Show feedback
        self.statusBar().showMessage(
            f"Panel layout set to {left_percent}% : {right_percent}%", 2000
        )

    def toggle_2d_panel(self):
        """Toggle 2D panel visibility."""
        sizes = self.splitter.sizes()
        if not sizes:
            return

        if sizes[0] == 0:
            # Restore if hidden
            self.restore_2d_panel()
            self.statusBar().showMessage("2D panel restored", 1500)
        else:
            # Minimize if shown
            self.minimize_2d_panel()
            self.statusBar().showMessage("2D panel minimized", 1500)

    def on_splitter_moved(self, pos, index):
        """Feedback for splitter movement."""
        sizes = self.splitter.sizes()
        if len(sizes) >= 2:
            total = sum(sizes)
            if total > 0:
                left_percent = round(sizes[0] * 100 / total)
                right_percent = round(sizes[1] * 100 / total)

                # Show ratio in tooltip
                if hasattr(self.splitter, "handle"):
                    handle = self.splitter.handle(1)
                    if handle:
                        handle.setToolTip(f"2D: {left_percent}% | 3D: {right_percent}%")

    def setup_splitter_tooltip(self):
        """Set initial splitter tooltip."""
        handle = self.splitter.handle(1)
        if handle:
            handle.setToolTip(
                "Drag to resize panels | Ctrl+1/2/3 for presets | Ctrl+H to toggle 2D panel"
            )
            # Show initial ratio
            self.on_splitter_moved(0, 0)
