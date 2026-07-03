#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

import logging
import os
from typing import Any, Callable, Dict, List, TYPE_CHECKING

from PyQt6.QtGui import QAction, QIcon, QKeySequence
from PyQt6.QtWidgets import QFileDialog, QMenu, QMessageBox

if TYPE_CHECKING:
    from .main_window_init import MainInitManager


class PluginMenuManager:
    """Handles all plugin-related menu and toolbar UI management.

    Extracted from MainInitManager to keep UI lifecycle logic separate from
    application-level initialization concerns.  All methods that discover,
    register, or rebuild plugin UI elements live here.

    The manager holds a reference to its owning MainInitManager (``_im``) so it
    can access shared state (host, buttons, toolbars) without duplicating them.
    """

    def __init__(self, init_manager: MainInitManager) -> None:
        self._im = init_manager

    def _make_safe_callback(self, callback: Callable, plugin_name: str) -> Callable:
        """Wrap a plugin callback so exceptions don't propagate into Qt's signal machinery."""

        def _safe(*args: Any, **kwargs: Any) -> None:
            try:
                callback()
            except Exception as exc:
                logging.exception("Plugin callback error (%s)", plugin_name)
                try:
                    QMessageBox.critical(
                        self._im.host,
                        "Plugin Error",
                        f"An error occurred in plugin '{plugin_name}':\n{exc}",
                    )
                except Exception:
                    logging.debug("Suppressed non-critical error", exc_info=True)

        return _safe

    # ------------------------------------------------------------------
    # Public API — called by MainInitManager and PluginManager
    # ------------------------------------------------------------------

    def update_plugin_menu(self, plugin_menu: QMenu) -> None:
        """Discover plugins and populate the Plugin menu from scratch."""
        if not self._im.host.plugin_manager:
            return

        self._clear_all_plugin_actions(plugin_menu)

        manage_action = QAction("Plugin Manager...", self._im.host)
        manage_action.triggered.connect(lambda: self._show_plugin_manager(plugin_menu))
        plugin_menu.addAction(manage_action)
        plugin_menu.addSeparator()

        plugins = self._im.host.plugin_manager.discover_plugins(self._im.host)

        self.update_style_menu_with_plugins()
        self.add_registered_plugin_actions()
        self.add_plugin_toolbar_actions()
        self._add_legacy_plugin_actions(plugin_menu, plugins)
        self.integrate_plugin_export_actions()
        self.integrate_plugin_file_openers()
        self.integrate_plugin_analysis_tools()
        self.integrate_plugin_optimization_methods()

    def rebuild_plugin_menus(self) -> None:
        """Fully rebuild all plugin-managed UI after an install/uninstall.

        Unlike ``update_plugin_menu``, this does not re-discover — it assumes
        ``PluginManager.discover_plugins`` has already been called and the
        registries are fresh.  It cleans the existing tagged actions and
        re-populates every integration point: menus, toolbars, export,
        file-openers, analysis tools, and 3D styles.
        """
        PLUGIN_ACTION_TAG = "plugin_managed"

        def _clean_menu(menu: QMenu) -> None:
            for action in list(menu.actions()):
                submenu = action.menu()
                if submenu is not None:
                    _clean_menu(submenu)
                    if not any(not a.isSeparator() for a in submenu.actions()):
                        menu.removeAction(action)
                elif action.data() == PLUGIN_ACTION_TAG:
                    menu.removeAction(action)

        try:
            for top_action in list(self._im.host.menuBar().actions()):
                top_menu = top_action.menu()
                if top_menu is not None:
                    _clean_menu(top_menu)
        except Exception as e:
            logging.warning("Plugin rebuild: menu cleanup error: %s", e)

        self._im.plugin_menubar_separator_added = False

        for method, label in [
            (self.add_registered_plugin_actions, "menu actions"),
            (self.add_plugin_toolbar_actions, "toolbar actions"),
            (self.integrate_plugin_export_actions, "export actions"),
            (self.integrate_plugin_file_openers, "file openers"),
            (self.integrate_plugin_analysis_tools, "analysis tools"),
            (self.update_style_menu_with_plugins, "style menu"),
            (self.integrate_plugin_optimization_methods, "optimization methods"),
        ]:
            try:
                method()
            except Exception as e:
                logging.warning("Plugin rebuild: %s error: %s", label, e)

    def add_registered_plugin_actions(self) -> None:
        """Add menu actions explicitly registered by V3/V4 plugins."""
        PLUGIN_ACTION_TAG = "plugin_managed"
        if not self._im.host.plugin_manager.menu_actions:
            return

        for action_def in self._im.host.plugin_manager.menu_actions:
            path = action_def["path"]
            callback = action_def["callback"]
            text = action_def["text"]

            parts = path.split("/")
            top_level_title = parts[0]
            current_menu = next(
                (
                    a.menu()
                    for a in self._im.host.menuBar().actions()
                    if a.menu() and a.text().replace("&", "") == top_level_title
                ),
                None,
            )

            if not current_menu:
                if not self._im.plugin_menubar_separator_added:
                    self._im.host.menuBar().addSeparator()
                    self._im.plugin_menubar_separator_added = True
                current_menu = self._im.host.menuBar().addMenu(top_level_title)

            for part in parts[1:-1]:
                sub = next(
                    (
                        a.menu()
                        for a in current_menu.actions()
                        if a.menu() and a.text().replace("&", "") == part
                    ),
                    None,
                )
                current_menu = sub if sub else current_menu.addMenu(part)

            actions = current_menu.actions()
            if (
                actions
                and not actions[-1].isSeparator()
                and actions[-1].data() != PLUGIN_ACTION_TAG
            ):
                sep = current_menu.addSeparator()
                sep.setData(PLUGIN_ACTION_TAG)

            action = QAction(text or parts[-1], self._im.host)
            action.triggered.connect(
                self._make_safe_callback(callback, action_def.get("plugin", "Plugin"))
            )
            if action_def.get("shortcut"):
                action.setShortcut(QKeySequence(action_def["shortcut"]))
            action.setData(PLUGIN_ACTION_TAG)
            current_menu.addAction(action)

    def add_plugin_toolbar_actions(self) -> None:
        """Populate the plugin toolbar from registered toolbar actions."""
        if not hasattr(self._im, "plugin_toolbar"):
            return

        self._im.plugin_toolbar.clear()  # type: ignore[union-attr]
        if self._im.host.plugin_manager.toolbar_actions:
            self._im.plugin_toolbar.show()  # type: ignore[union-attr]
            for action_def in self._im.host.plugin_manager.toolbar_actions:
                action = QAction(action_def["text"], self._im.host)
                action.triggered.connect(
                    self._make_safe_callback(
                        action_def["callback"], action_def.get("plugin", "Plugin")
                    )
                )
                if action_def["icon"] and os.path.exists(action_def["icon"]):
                    action.setIcon(QIcon(action_def["icon"]))
                if action_def["tooltip"]:
                    action.setToolTip(action_def["tooltip"])
                self._im.plugin_toolbar.addAction(action)  # type: ignore[union-attr]
        else:
            self._im.plugin_toolbar.hide()  # type: ignore[union-attr]

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _show_plugin_manager(self, plugin_menu: QMenu) -> None:
        """Open the Plugin Manager dialog and refresh the menu on close."""
        if not self._im.host.plugin_manager:
            QMessageBox.information(
                self._im.host, "Safe Mode", "Plugins are disabled (safe mode)."
            )
            return
        from ..plugins.plugin_manager_window import PluginManagerWindow

        dlg = PluginManagerWindow(self._im.host.plugin_manager, self._im.host)
        dlg.exec()
        self.update_plugin_menu(plugin_menu)

    def _clear_all_plugin_actions(self, plugin_menu: QMenu) -> None:
        """Remove all tagged plugin actions from every menu and the export button."""
        PLUGIN_ACTION_TAG = "plugin_managed"

        def clear_menu(menu: Any) -> None:
            if not menu:
                return
            for act in list(menu.actions()):
                if act.data() == PLUGIN_ACTION_TAG:
                    menu.removeAction(act)
                elif act.menu():
                    clear_menu(act.menu())

        plugin_menu.clear()
        for top_action in self._im.host.menuBar().actions():
            if top_action.menu():
                clear_menu(top_action.menu())

        export_button = getattr(self._im, "export_button", None)
        if export_button and export_button.menu():
            clear_menu(export_button.menu())

    def update_style_menu_with_plugins(self) -> None:
        """Append custom 3D styles registered by plugins to the style menu."""
        style_button = getattr(self._im, "style_button", None)
        if not style_button or not style_button.menu():
            return

        style_menu = style_button.menu()
        style_group = next(
            (a.actionGroup() for a in style_menu.actions() if a.actionGroup()), None
        )

        if style_group and self._im.host.plugin_manager.custom_3d_styles:
            if style_menu.actions() and not style_menu.actions()[-1].isSeparator():
                style_menu.addSeparator()

            for style_name in self._im.host.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self._im.host)
                    action.setCheckable(True)
                    action.triggered.connect(
                        lambda checked=False,
                        s=style_name: self._im.host.view_3d_manager.set_3d_style(s)
                    )
                    style_menu.addAction(action)
                    style_group.addAction(action)

    def _add_legacy_plugin_actions(
        self, plugin_menu: QMenu, plugins: List[Any]
    ) -> None:
        """Add run()-based (legacy) plugins as menu entries under the Plugin menu."""
        if not plugins:
            no_plugin = QAction("(No plugins found)", self._im.host)
            no_plugin.setEnabled(False)
            plugin_menu.addAction(no_plugin)
            return

        categorized: Dict[str, Any] = {}
        root = []
        for p in plugins:
            if hasattr(p["module"], "run"):
                cat = p.get("category", p.get("rel_folder", "")).strip()
                if cat:
                    categorized.setdefault(cat, []).append(p)
                else:
                    root.append(p)

        for cat in sorted(categorized.keys()):
            current_parent: Any = plugin_menu
            for part in cat.split(os.sep):
                sub = next(
                    (
                        a.menu()
                        for a in current_parent.actions()
                        if a.menu() and a.text().replace("&", "") == part
                    ),
                    None,
                )
                current_parent = (
                    sub if sub is not None else current_parent.addMenu(part)
                )

            for p in sorted(categorized[cat], key=lambda x: x["name"]):
                a = QAction(p["name"], self._im.host)
                a.triggered.connect(
                    lambda checked,
                    mod=p["module"]: self._im.host.plugin_manager.run_plugin(
                        mod, self._im.host
                    )
                )
                current_parent.addAction(a)

        for p in sorted(root, key=lambda x: x["name"]):
            a = QAction(p["name"], self._im.host)
            a.triggered.connect(
                lambda checked,
                mod=p["module"]: self._im.host.plugin_manager.run_plugin(
                    mod, self._im.host
                )
            )
            plugin_menu.addAction(a)

    def integrate_plugin_export_actions(self) -> None:
        """Inject plugin export actions into the File > Export menu and export button."""
        if not self._im.host.plugin_manager.export_actions:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        main_export_menu = None
        for top_action in self._im.host.menuBar().actions():
            if top_action.menu() and top_action.text().replace("&", "") == "File":
                main_export_menu = next(
                    (
                        a.menu()
                        for a in top_action.menu().actions()
                        if a.menu() and a.text().replace("&", "") == "Export"
                    ),
                    None,
                )
                if main_export_menu:
                    break

        targets = []
        export_button = getattr(self._im, "export_button", None)
        if export_button and export_button.menu():
            targets.append(export_button.menu())
        if main_export_menu:
            targets.append(main_export_menu)

        for menu in targets:
            sep = menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for exp in self._im.host.plugin_manager.export_actions:
                a = QAction(exp["label"], self._im.host)
                a.triggered.connect(
                    self._make_safe_callback(
                        exp["callback"], exp.get("plugin", "Plugin")
                    )
                )
                a.setData(PLUGIN_ACTION_TAG)
                menu.addAction(a)

    def integrate_plugin_optimization_methods(self) -> None:
        """Inject plugin optimization methods into the 3D Optimization Settings menu."""
        if not hasattr(self._im.host.plugin_manager, "optimization_methods"):
            return

        for key, entry in self._im.host.plugin_manager.optimization_methods.items():
            self._im.add_optimization_method(entry["label"], key)
            if key in self._im.opt3d_actions:
                self._im.opt3d_actions[key].setData("plugin_managed")

    def integrate_plugin_file_openers(self) -> None:
        """Inject plugin file-opener entries into the File > Import menu."""
        if not self._im.host.plugin_manager.file_openers:
            return

        import_menu = getattr(self._im, "import_menu", None) or getattr(
            self._im.host, "import_menu", None
        )
        if import_menu is None:
            return

        PLUGIN_ACTION_TAG = "plugin_managed"
        sep = import_menu.addSeparator()
        sep.setData(PLUGIN_ACTION_TAG)

        plugin_map: Dict[str, Any] = {}
        for ext, openers in self._im.host.plugin_manager.file_openers.items():
            for info in openers:
                plugin_map.setdefault(info.get("plugin", "Plugin"), {})[ext] = info[
                    "callback"
                ]

        for p_name, ext_map in sorted(plugin_map.items()):
            exts = sorted(ext_map.keys())
            ext_str = "/".join(exts)
            if len(ext_str) > 30:
                cutoff = ext_str.rfind("/", 0, 30)
                ext_str = (ext_str[:cutoff] if cutoff != -1 else ext_str[:30]) + "/..."

            filter_str = (
                f"{p_name} Files ({' '.join(['*' + e for e in exts])});;All Files (*)"
            )

            def make_cb(m: Any, f: Any, n: Any) -> Any:
                def _cb() -> None:
                    fpath, _ = QFileDialog.getOpenFileName(
                        self._im.host, f"Import {n} Files", "", f
                    )
                    if fpath:
                        ext = os.path.splitext(fpath)[1].lower()
                        if ext in m:
                            try:
                                m[ext](fpath)
                            except Exception as exc:
                                logging.exception("Plugin file opener error (%s)", n)
                                QMessageBox.critical(
                                    self._im.host,
                                    "Plugin Error",
                                    f"An error occurred in plugin '{n}':\n{exc}",
                                )
                                return
                            self._im.current_file_path = fpath
                            self._im.host.state_manager.update_window_title()

                return _cb

            a = QAction(f"Import {ext_str} ({p_name})...", self._im.host)
            a.triggered.connect(make_cb(ext_map, filter_str, p_name))
            a.setData(PLUGIN_ACTION_TAG)
            import_menu.addAction(a)

    def integrate_plugin_analysis_tools(self) -> None:
        """Inject plugin analysis tools into the Analysis menu."""
        analysis_menu = next(
            (
                a.menu()
                for a in self._im.host.menuBar().actions()
                if a.text().replace("&", "") == "Analysis"
            ),
            None,
        )
        if analysis_menu and self._im.host.plugin_manager.analysis_tools:
            PLUGIN_ACTION_TAG = "plugin_managed"
            sep = analysis_menu.addSeparator()
            sep.setData(PLUGIN_ACTION_TAG)
            for tool in self._im.host.plugin_manager.analysis_tools:
                a = QAction(
                    f"{tool['label']} ({tool.get('plugin', 'Plugin')})", self._im.host
                )
                a.triggered.connect(
                    self._make_safe_callback(
                        tool["callback"], tool.get("plugin", "Plugin")
                    )
                )
                a.setData(PLUGIN_ACTION_TAG)
                analysis_menu.addAction(a)

    def _clear_plugin_ui_elements(self, plugin_menu: QMenu) -> None:
        """Remove only the legacy-tagged plugin actions from the plugin menu."""
        for action in plugin_menu.actions():
            if action.data() == "plugin_action":
                plugin_menu.removeAction(action)

        self._im.plugin_toolbar.clear()  # type: ignore[union-attr]
        self._im.plugin_toolbar.hide()  # type: ignore[union-attr]
