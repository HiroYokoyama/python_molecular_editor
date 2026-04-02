#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging  # [REPORT ERROR MISSING ATTRIBUTE]

import vtk

# PyQt6 Modules
from PyQt6.QtCore import QEvent, Qt, QTimer, QObject
from PyQt6.QtWidgets import (
    QApplication,
    QGraphicsView,
    QMainWindow,
    QMessageBox,
    QWidget,
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
    from moleditpy.ui.custom_interactor_style import CustomInteractorStyle


# --- Classes ---
class UIManager(QObject):
    def __init__(self, host=None):
        super().__init__()
        if host is not None:
            self.host = host

    def update_status_bar(self, message):
        """Update status bar with worker messages."""
        self.host.statusBar().showMessage(message)

    def set_mode(self, mode_str):
        if isinstance(mode_str, tuple):
            mode_str = (
                f"bond_{mode_str[0]}_{mode_str[1]}"
                if len(mode_str) == 2
                else str(mode_str[0])
            )

        prev_mode = getattr(self.host.init_manager.scene, "mode", None)
        self.host.init_manager.scene.mode = mode_str
        self.host.init_manager.view_2d.setMouseTracking(True)

        # Trigger immediate scene refresh to show/update template previews
        if hasattr(self.host.init_manager.scene, "refresh_mode_state"):
            self.host.init_manager.scene.refresh_mode_state()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'refresh_mode_state' on object"
            )
        # Clear ghost when leaving template mode
        if (
            prev_mode
            and prev_mode.startswith("template")
            and not mode_str.startswith("template")
        ):
            self.host.init_manager.scene.clear_template_preview()
        elif not mode_str.startswith("template"):
            self.host.init_manager.scene.template_preview.hide()

        # Set cursor shape
        if mode_str == "select":
            self.host.init_manager.view_2d.setCursor(Qt.CursorShape.ArrowCursor)
        elif mode_str.startswith(("atom", "bond", "template")):
            self.host.init_manager.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        elif mode_str.startswith(("charge", "radical")):
            self.host.init_manager.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self.host.init_manager.view_2d.setCursor(Qt.CursorShape.ArrowCursor)

        if mode_str.startswith("atom"):
            self.host.init_manager.scene.current_atom_symbol = mode_str.split("_")[1]
            self.host.statusBar().showMessage(
                f"Mode: Draw Atom ({self.host.init_manager.scene.current_atom_symbol})"
            )
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.host.init_manager.view_2d.setMouseTracking(True)
            self.host.init_manager.scene.bond_order = 1
            self.host.init_manager.scene.bond_stereo = 0
        elif mode_str.startswith("bond"):
            self.host.init_manager.scene.current_atom_symbol = "C"
            parts = mode_str.split("_")
            self.host.init_manager.scene.bond_order = int(parts[1])
            self.host.init_manager.scene.bond_stereo = (
                int(parts[2]) if len(parts) > 2 else 0
            )
            stereo_text = {0: "", 1: " (Wedge)", 2: " (Dash)"}.get(
                self.host.init_manager.scene.bond_stereo, ""
            )
            self.host.statusBar().showMessage(
                f"Mode: Draw Bond (Order: {self.host.init_manager.scene.bond_order}{stereo_text})"
            )
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.host.init_manager.view_2d.setMouseTracking(True)
        elif mode_str.startswith("template"):
            if mode_str.startswith("template_user"):
                # User template mode
                template_name = mode_str.replace("template_user_", "")
                self.host.statusBar().showMessage(
                    f"Mode: User Template ({template_name})"
                )
            else:
                # Built-in template mode
                self.host.statusBar().showMessage(
                    f"Mode: {mode_str.split('_')[1].capitalize()} Template"
                )
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "charge_plus":
            self.host.statusBar().showMessage("Mode: Increase Charge (Click on Atom)")
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "charge_minus":
            self.host.statusBar().showMessage("Mode: Decrease Charge (Click on Atom)")
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == "radical":
            self.host.statusBar().showMessage("Mode: Toggle Radical (Click on Atom)")
            self.host.init_manager.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)

        else:  # Select mode
            self.host.statusBar().showMessage("Mode: Select")
            self.host.init_manager.view_2d.setDragMode(
                QGraphicsView.DragMode.RubberBandDrag
            )
            self.host.init_manager.scene.bond_order = 1
            self.host.init_manager.scene.bond_stereo = 0

    def set_mode_and_update_toolbar(self, mode_str):
        self.set_mode(mode_str)
        # Map QAction to QToolButton
        toolbar = getattr(self.host.init_manager, "toolbar", None)
        action_to_button = {}
        if toolbar:
            for key, action in self.host.init_manager.mode_actions.items():
                btn = toolbar.widgetForAction(action)
                if btn:
                    action_to_button[action] = btn

        # Reset all mode buttons
        for key, action in self.host.init_manager.mode_actions.items():
            action.setChecked(False)
            btn = action_to_button.get(action)
            if btn:
                btn.setStyleSheet("")

        # Apply style to matching mode buttons (exact match or prefix for user templates)
        matched_key = None
        if mode_str in self.host.init_manager.mode_actions:
            matched_key = mode_str
        elif mode_str.startswith("template_user"):
            matched_key = "template_user"

        if matched_key and matched_key in self.host.init_manager.mode_actions:
            action = self.host.init_manager.mode_actions[matched_key]
            action.setChecked(True)
            btn = action_to_button.get(action)
            if btn:
                # Highlight templates with specific color
                if mode_str.startswith("template"):
                    btn.setStyleSheet(
                        "background-color: #2196F3; color: white; border-radius: 4px;"
                    )
                else:
                    btn.setStyleSheet("")

    def activate_select_mode(self):
        self.set_mode("select")
        if "select" in self.host.init_manager.mode_actions:
            self.host.init_manager.mode_actions["select"].setChecked(True)

    def set_atom_from_periodic_table(self, symbol: str) -> None:
        """Helper to set the current mode from periodic table selection."""
        self.set_mode(f"atom_{symbol}")

    def eventFilter(self, obj, event):
        if (
            hasattr(self.host.view_3d_manager, "plotter")
            and obj is self.host.view_3d_manager.plotter
            and event.type() == QEvent.Type.MouseButtonPress
        ):
            self.host.init_manager.view_2d.setFocus()

        if obj is self.host:
            # Handle Drag and Drop via event filter
            if event.type() == QEvent.Type.DragEnter:
                self.handle_drag_enter_event(event)
                return True
            if event.type() == QEvent.Type.Drop:
                self.handle_drop_event(event)
                return True

            # Handle Window Close via event filter
            if event.type() == QEvent.Type.Close:
                if not self.handle_close_event(event):
                    event.ignore()
                    return True  # Stop propagation

        return super().eventFilter(obj, event)

    def handle_close_event(self, event) -> bool:
        """
        Handle application close logic.
        Returns True if close should proceed, False if it should be cancelled.
        """
        # 1. Persist settings
        try:
            modified = getattr(self.host.init_manager, "settings_dirty", False) or (
                self.host.init_manager.settings != self.host.initial_settings
            )
            if modified:
                self.host.init_manager.save_settings()
                self.host.init_manager.settings_dirty = False
        except (AttributeError, RuntimeError, TypeError, ValueError, OSError):
            pass

        # 2. Handle unsaved changes
        if getattr(self.host.state_manager, "has_unsaved_changes", False):
            reply = QMessageBox.question(
                self.host,
                "Unsaved Changes",
                "You have unsaved changes. Do you want to save them?",
                QMessageBox.StandardButton.Yes
                | QMessageBox.StandardButton.No
                | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Yes,
            )

            if reply == QMessageBox.StandardButton.Yes:
                self.host.io_manager.save_project()
                if getattr(self.host.state_manager, "has_unsaved_changes", False):
                    return False
            elif reply == QMessageBox.StandardButton.Cancel:
                return False

        # 3. Gracefully close child windows and cleanup threads
        try:
            for widget in QApplication.topLevelWidgets():
                if (
                    widget is not None
                    and widget != self.host
                    and isinstance(widget, (QWidget, QMainWindow))
                ):
                    try:
                        widget.close()
                    except (RuntimeError, TypeError):
                        pass

            # Stop calculation threads
            active_threads = list(
                getattr(self.host.compute_manager, "_active_calc_threads", []) or []
            )
            for thr in active_threads:
                try:
                    if hasattr(thr, "quit"):
                        thr.quit()
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error("REPORT ERROR: Missing attribute 'quit' on thr")
                    if hasattr(thr, "wait"):
                        thr.wait(200)
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error("REPORT ERROR: Missing attribute 'wait' on thr")
                except (RuntimeError, TypeError):
                    pass
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass

        return True

    def closeEvent(self, event):
        """Deprecated: Use handle_close_event via eventFilter."""
        if self.handle_close_event(event):
            event.accept()
        else:
            event.ignore()

    def toggle_3d_edit_mode(self, checked):
        """Toggle 3D Drag mode."""
        if checked:
            # Disable measurement mode when 3D Drag is on
            if self.host.edit_3d_manager.measurement_mode:
                self.host.init_manager.measurement_action.setChecked(False)
                self.host.edit_3d_manager.toggle_measurement_mode(False)

        self.host.edit_3d_manager.is_3d_edit_mode = checked
        if checked:
            self.host.statusBar().showMessage("3D Drag Mode: ON.")
        else:
            self.host.statusBar().showMessage("3D Drag Mode: OFF.")
        self.host.init_manager.view_2d.setFocus()

    def _setup_3d_picker(self):
        self.host.view_3d_manager.plotter.picker = vtk.vtkCellPicker()
        self.host.view_3d_manager.plotter.picker.SetTolerance(0.025)

        # Create CustomInteractorStyle
        style = CustomInteractorStyle(self.host)

        # Set interactor style
        self.host.view_3d_manager.plotter.interactor.SetInteractorStyle(style)
        self.host.view_3d_manager.plotter.interactor.Initialize()


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
            plugin_mgr = getattr(self.host, "plugin_manager", None)
            if plugin_mgr and getattr(plugin_mgr, "drop_handlers", []):
                event.acceptProposedAction()
                return

        event.ignore()


    def handle_drop_event(self, event):
        """Internal handler for file drop event (bypasses PyQt type checks in tests)."""
        urls = event.mimeData().urls()
        file_path = next((u.toLocalFile() for u in urls if u.isLocalFile()), None)
        if not file_path:
            event.ignore()
            return

        # 1. Plugin Handlers
        plugin_mgr = getattr(self.host, "plugin_manager", None)
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
            self.host.io_manager.open_project_file(file_path=file_path)
            QTimer.singleShot(100, self.host.view_3d_manager.fit_to_view)
            event.acceptProposedAction()
        elif file_lower.endswith((".mol", ".sdf")):
            # Robust check for drop target using childAt
            drag_point = event.position().toPoint()
            target_widget = self.host.childAt(drag_point)

            # Identify if the target widget is the plotter (or one of its children)
            is_on_3d = False
            if hasattr(self.host.view_3d_manager, "plotter"):
                plotter_widget = self.host.init_manager.splitter.widget(1)
                if target_widget == plotter_widget or plotter_widget.isAncestorOf(
                    target_widget
                ):
                    is_on_3d = True
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error("REPORT ERROR: Missing attribute 'plotter' on object")

            if is_on_3d:
                self.host.io_manager.load_mol_file_for_3d_viewing(file_path=file_path)
                # Ensure 3D viewer zooms and renders the dropped molecule
                QTimer.singleShot(
                    100, lambda: self.host.view_3d_manager.plotter.view_isometric()
                )
                QTimer.singleShot(
                    150, lambda: self.host.view_3d_manager.plotter.render()
                )
            else:
                self.host.io_manager.load_mol_file(file_path=file_path)
                QTimer.singleShot(100, self.host.view_3d_manager.fit_to_view)

            event.acceptProposedAction()
        elif file_lower.endswith(".xyz"):
            self.host.io_manager.load_xyz_for_3d_viewing(file_path=file_path)
            QTimer.singleShot(100, self.host.view_3d_manager.fit_to_view)
            event.acceptProposedAction()
        else:
            self.host.statusBar().showMessage(f"Unsupported file type: {file_path}")
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
            if hasattr(self.host, action_name):
                getattr(self.host, action_name).setEnabled(enabled)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    f"REPORT ERROR: Missing attribute {action_name} on self.host"
                )

        for menu_name in menus:
            if hasattr(self.host, menu_name):
                getattr(self.host, menu_name).setEnabled(enabled)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error(
                    f"REPORT ERROR: Missing attribute {menu_name} on self.host"
                )

    def _enable_3d_features(self, enabled=True):
        """Enable/disable 3D features."""
        # Basic 3D features
        basic_3d_actions = ["optimize_3d_button", "export_button", "analysis_action"]

        for action_name in basic_3d_actions:
            # Check both host and init_manager for these UI components
            obj = getattr(self.host, action_name, None)
            if obj is None:
                obj = getattr(self.host.init_manager, action_name, None)

            if obj is None:
                continue

            try:
                if action_name == "optimize_3d_button":
                    # Optimization is disabled for XYZ-derived or failed-chem-check molecules
                    is_xyz = getattr(self.host, "is_xyz_derived", False)
                    chem_failed = getattr(
                        self.host, "chem_check_tried", False
                    ) and getattr(self.host, "chem_check_failed", False)

                    can_optimize = enabled and not (is_xyz or chem_failed)
                    obj.setEnabled(can_optimize)
                else:
                    obj.setEnabled(enabled)
            except (AttributeError, RuntimeError, TypeError, ValueError):
                # Suppress non-critical 3D feature state update errors if widgets are not fully initialized
                pass

        # Always enable these core 3D interactors
        for core_act in ["measurement_action", "edit_3d_action"]:
            obj = getattr(self.host, core_act, None)
            if obj:
                obj.setEnabled(enabled)

        self._enable_3d_edit_actions(enabled)

    def _enter_3d_viewer_ui_mode(self):
        """Set UI mode to 3D viewer."""
        self.host.ui_manager.is_2d_editable = False
        self.host.init_manager.cleanup_button.setEnabled(False)
        self.host.init_manager.convert_button.setEnabled(False)
        for action in self.host.init_manager.tool_group.actions():
            action.setEnabled(False)
        if hasattr(self.host.init_manager, "other_atom_action"):
            self.host.init_manager.other_atom_action.setEnabled(False)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'other_atom_action' on object"
            )

        self.host.ui_manager.minimize_2d_panel()

        # Collectively enable 3D-related features
        self.host.ui_manager._enable_3d_features(True)

    def restore_ui_for_editing(self):
        """Enables all 2D editing UI elements."""
        self.host.ui_manager.is_2d_editable = True
        self.host.ui_manager.restore_2d_panel()
        self.host.init_manager.cleanup_button.setEnabled(True)
        self.host.init_manager.convert_button.setEnabled(True)

        for action in self.host.init_manager.tool_group.actions():
            action.setEnabled(True)

        if hasattr(self.host.init_manager, "other_atom_action"):
            self.host.init_manager.other_atom_action.setEnabled(True)
        else:
            logging.error(
                "REPORT ERROR: Missing attribute 'other_atom_action' on self.host.init_manager"
            )

        # Collectively disable 3D edit functions when returning to 2D mode
        self.host.ui_manager._enable_3d_edit_actions(False)

    def minimize_2d_panel(self):
        """Minimize (hide) 2D panel."""
        sizes = self.host.init_manager.splitter.sizes()
        # Only if not already minimized
        if sizes[0] > 0:
            total_width = sum(sizes)
            self.host.init_manager.splitter.setSizes([0, total_width])

    def restore_2d_panel(self):
        """Restore 2D panel."""
        sizes = self.host.init_manager.splitter.sizes()

        # Check sizes list before access
        if sizes and sizes[0] == 0:
            self.host.init_manager.splitter.setSizes([600, 600])

    def set_panel_layout(self, left_percent, right_percent):
        """Set panel layout ratio."""
        if left_percent + right_percent != 100:
            return

        total_width = self.host.init_manager.splitter.width()
        if total_width <= 0:
            total_width = 1200  # Default width

        left_width = int(total_width * left_percent / 100)
        right_width = int(total_width * right_percent / 100)

        self.host.init_manager.splitter.setSizes([left_width, right_width])

        # Show feedback
        self.host.statusBar().showMessage(
            f"Panel layout set to {left_percent}% : {right_percent}%", 2000
        )

    def toggle_2d_panel(self):
        """Toggle 2D panel visibility."""
        sizes = self.host.init_manager.splitter.sizes()
        if not sizes:
            return

        if sizes[0] == 0:
            # Restore if hidden
            self.host.ui_manager.restore_2d_panel()
            self.host.statusBar().showMessage("2D panel restored", 1500)
        else:
            # Minimize if shown
            self.host.ui_manager.minimize_2d_panel()
            self.host.statusBar().showMessage("2D panel minimized", 1500)

    def _get_live_splitter(self):
        """Return splitter only when wrapper and C++ object are still alive."""
        init_manager = getattr(self.host, "init_manager", None)
        splitter = getattr(init_manager, "splitter", None)
        if splitter is None:
            return None

        if _sip_isdeleted is not None:
            try:
                if _sip_isdeleted(splitter):
                    return None
            except (AttributeError, RuntimeError, TypeError):
                return None

        return splitter

    def on_splitter_moved(self, pos, index):
        """Feedback for splitter movement."""
        splitter = self._get_live_splitter()
        if splitter is None:
            return

        try:
            sizes = splitter.sizes()
        except (AttributeError, RuntimeError, TypeError):
            return

        if len(sizes) >= 2:
            total = sum(sizes)
            if total > 0:
                left_percent = round(sizes[0] * 100 / total)
                right_percent = round(sizes[1] * 100 / total)

                # Show ratio in tooltip
                if hasattr(splitter, "handle"):
                    try:
                        handle = splitter.handle(1)
                    except (AttributeError, RuntimeError, TypeError):
                        return
                    if handle:
                        try:
                            handle.setToolTip(
                                f"2D: {left_percent}% | 3D: {right_percent}%"
                            )
                        except (AttributeError, RuntimeError, TypeError):
                            return
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error("REPORT ERROR: Missing attribute 'handle' on object")

    def setup_splitter_tooltip(self):
        """Set initial splitter tooltip."""
        splitter = self._get_live_splitter()
        if splitter is None:
            return

        try:
            handle = splitter.handle(1)
        except (AttributeError, RuntimeError, TypeError):
            return

        if handle:
            try:
                handle.setToolTip(
                    "Drag to resize panels | Ctrl+1/2/3 for presets | Ctrl+H to toggle 2D panel"
                )
            except (AttributeError, RuntimeError, TypeError):
                return
            # Show initial ratio
            self.on_splitter_moved(0, 0)


# Backward-compat aliases
MainWindowUiManager = UIManager
