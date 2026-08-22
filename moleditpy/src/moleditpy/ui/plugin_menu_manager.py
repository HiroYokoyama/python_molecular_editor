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

    # QAction.data() marker stamped on every plugin-owned menu entry/separator
    # so the rebuild pass can find and remove exactly those on the next reload.
    _PLUGIN_ACTION_TAG = "plugin_managed"

    # Title of the menu that hosts "Plugin/<path>" registrations.
    _PLUGIN_MENU_TITLE = "Plugin"

    # Value of add_menu_action(pin=...) requesting the header slot: placed with
    # "Plugin Manager..." above the first divider, for plugins that manage the
    # plugin system itself rather than acting on a molecule.
    _PIN_HEADER = "header"

    def __init__(self, init_manager: MainInitManager) -> None:
        self._im = init_manager

    def _make_safe_callback(self, callback: Callable, plugin_name: str) -> Callable:
        """Wrap a plugin callback so exceptions don't propagate into Qt's signal machinery."""

        def _safe(*args: Any, **kwargs: Any) -> None:
            try:
                callback()
            except Exception:
                logging.exception("Plugin callback error (%s)", plugin_name)

        return _safe

    @staticmethod
    def _remove_action(menu: QMenu, action: QAction) -> None:
        """Remove *action* from *menu* and from its QActionGroup, if any.

        Checkable plugin actions (e.g. custom 3D styles) are added to an
        exclusive QActionGroup. ``menu.removeAction`` alone leaves them in the
        group — and, being parented to the host, alive — so each rebuild would
        leak a stale member into the group and could drop the active item's
        check state. Detaching from the group first keeps the group in sync.
        """
        group = action.actionGroup()
        if group is not None:
            group.removeAction(action)
        menu.removeAction(action)

    def _get_plugin_target_menus(self) -> List[QMenu]:
        """Return a list of all menus that might contain plugin actions."""
        menus = []
        for top_action in list(self._im.host.menuBar().actions()):
            if top_action.menu():
                menus.append(top_action.menu())

        export_button = getattr(self._im, "export_button", None)
        if export_button and export_button.menu():
            menus.append(export_button.menu())

        style_button = getattr(self._im, "style_button", None)
        if style_button and style_button.menu():
            menus.append(style_button.menu())

        import_menu = getattr(self._im, "import_menu", None) or getattr(
            self._im.host, "import_menu", None
        )
        if import_menu:
            menus.append(import_menu)

        return menus

    # ------------------------------------------------------------------
    # Public API — called by MainInitManager and PluginManager
    # ------------------------------------------------------------------

    def update_plugin_menu(self, plugin_menu: QMenu) -> None:
        """Discover plugins and populate the Plugin menu from scratch."""
        if not self._im.host.plugin_manager:
            return

        self._clear_all_plugin_actions()
        self._reset_plugin_menu(plugin_menu)

        plugins = self._im.host.plugin_manager.discover_plugins(self._im.host)

        self.update_style_menu_with_plugins()
        self.add_registered_plugin_actions()
        self.add_plugin_toolbar_actions()
        self._add_legacy_plugin_actions(plugin_menu, plugins)
        self.integrate_plugin_export_actions()
        self.integrate_plugin_file_openers()
        self.integrate_plugin_analysis_tools()
        self.integrate_plugin_optimization_methods()

    def _reset_plugin_menu(self, plugin_menu: QMenu) -> None:
        """Empty the Plugin menu and re-add its fixed header."""
        plugin_menu.clear()
        manage_action = QAction("Plugin Manager...", self._im.host)
        manage_action.triggered.connect(lambda: self._show_plugin_manager(plugin_menu))
        plugin_menu.addAction(manage_action)
        plugin_menu.addSeparator()

    def rebuild_plugin_menus(self) -> None:
        """Fully rebuild all plugin-managed UI after an install/uninstall.

        Unlike ``update_plugin_menu``, this does not re-discover — it assumes
        ``PluginManager.discover_plugins`` has already been called and the
        registries are fresh.  It cleans the existing tagged actions and
        re-populates every integration point: menus, toolbars, export,
        file-openers, analysis tools, legacy run() entries, and 3D styles.

        The Plugin menu is emptied and rebuilt rather than topped up, so a
        reload lands the context-injected entries above the folder-derived ones
        exactly as a fresh launch does. Appending to the surviving (untagged)
        folder entries instead used to flip the two groups on every reload.
        """

        def _clean_menu(menu: QMenu) -> None:
            for action in list(menu.actions()):
                submenu = action.menu()
                if submenu is not None:
                    _clean_menu(submenu)
                    # Only drop submenus the plugin system created. Native
                    # submenus that happen to be empty are not ours to remove.
                    if action.data() == self._PLUGIN_ACTION_TAG and not any(
                        not a.isSeparator() for a in submenu.actions()
                    ):
                        menu.removeAction(action)
                elif action.data() == self._PLUGIN_ACTION_TAG:
                    self._remove_action(menu, action)

        try:
            for menu in self._get_plugin_target_menus():
                _clean_menu(menu)
        except Exception:
            logging.warning("Plugin rebuild: menu cleanup error", exc_info=True)

        self._clear_plugin_menubar_entries()

        plugin_menu = getattr(self._im, "plugin_menu", None)
        if plugin_menu is not None:
            self._reset_plugin_menu(plugin_menu)

        def _rebuild_legacy_actions() -> None:
            if plugin_menu is not None:
                self._add_legacy_plugin_actions(
                    plugin_menu, self._im.host.plugin_manager.plugins
                )

        for method, label in [
            (self.add_registered_plugin_actions, "menu actions"),
            (_rebuild_legacy_actions, "legacy plugin actions"),
            (self.add_plugin_toolbar_actions, "toolbar actions"),
            (self.integrate_plugin_export_actions, "export actions"),
            (self.integrate_plugin_file_openers, "file openers"),
            (self.integrate_plugin_analysis_tools, "analysis tools"),
            (self.update_style_menu_with_plugins, "style menu"),
            (self.integrate_plugin_optimization_methods, "optimization methods"),
        ]:
            try:
                method()
            except Exception:
                logging.warning("Plugin rebuild: %s error", label, exc_info=True)

    def _clear_plugin_menubar_entries(self) -> None:
        """Drop emptied plugin-created top-level menus and the menubar divider.

        Both clear paths need this. ``_get_plugin_target_menus`` walks *into*
        the menus on the bar but never the bar itself, so without this pass an
        uninstalled plugin's own top-level menu survives as an empty stub and
        the divider in front of it is never reclaimed. Resetting the flag lets
        the following populate pass re-add that divider.
        """
        try:
            menu_bar = self._im.host.menuBar()
            for top_action in list(menu_bar.actions()):
                if top_action.data() != self._PLUGIN_ACTION_TAG:
                    continue
                submenu = top_action.menu()
                if submenu is None or not any(
                    not a.isSeparator() for a in submenu.actions()
                ):
                    menu_bar.removeAction(top_action)
        except Exception:
            logging.warning("Plugin menubar cleanup error", exc_info=True)

        self._im.plugin_menubar_separator_added = False

    def _ensure_plugin_separator(self, menu: QMenu) -> None:
        """Add a tagged divider before the first plugin entry appended to *menu*.

        The divider must land where the plugin content actually starts. For a
        nested path like ``Tools/Sub/Item`` that is the *parent* menu, just
        before ``Sub`` is created -- checking the freshly made (and therefore
        empty) leaf instead placed no divider at all, and the next sibling
        entry then saw an untagged submenu as the last action and inserted the
        divider after the plugin submenu, splitting the plugin group instead of
        separating it from the native entries.
        """
        actions = menu.actions()
        if (
            actions
            and not actions[-1].isSeparator()
            and actions[-1].data() != self._PLUGIN_ACTION_TAG
        ):
            sep = menu.addSeparator()
            if sep is not None:
                sep.setData(self._PLUGIN_ACTION_TAG)

    @classmethod
    def _wants_header_slot(
        cls, top_level_title: str, parts: List[str], action_def: Dict[str, Any]
    ) -> bool:
        """Return True for an entry that asked to sit beside the Plugin Manager.

        A plugin that manages the plugin system itself — the Plugin Installer —
        reads as part of the menu's header pair rather than as one more
        installed plugin, and says so with ``pin="header"``. The app therefore
        needs to know no plugin by name.
        """
        return (
            top_level_title == cls._PLUGIN_MENU_TITLE
            and len(parts) == 2
            and action_def.get("pin") == cls._PIN_HEADER
        )

    @staticmethod
    def _trim_trailing_separator(menu: QMenu) -> None:
        """Drop a divider left dangling at the end of *menu*."""
        actions = menu.actions()
        if actions and actions[-1].isSeparator():
            menu.removeAction(actions[-1])

    @staticmethod
    def _add_beside_plugin_manager(menu: QMenu, action: QAction) -> None:
        """Add *action* to the menu's header group, above the first divider."""
        anchor = next((a for a in menu.actions() if a.isSeparator()), None)
        if anchor is None:
            menu.addAction(action)
        else:
            menu.insertAction(anchor, action)

    def add_registered_plugin_actions(self) -> None:
        """Add menu actions explicitly registered by V3/V4 plugins."""
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
                    sep = self._im.host.menuBar().addSeparator()
                    sep.setData(self._PLUGIN_ACTION_TAG)
                    self._im.plugin_menubar_separator_added = True
                current_menu = self._im.host.menuBar().addMenu(top_level_title)
                current_menu.menuAction().setData(self._PLUGIN_ACTION_TAG)

            for part in parts[1:-1]:
                sub = next(
                    (
                        a.menu()
                        for a in current_menu.actions()
                        if a.menu() and a.text().replace("&", "") == part
                    ),
                    None,
                )
                if sub is not None:
                    # Descending into an existing submenu appends nothing to
                    # the parent, so no divider is due here.
                    current_menu = sub
                    continue
                self._ensure_plugin_separator(current_menu)
                current_menu = current_menu.addMenu(part)
                current_menu.menuAction().setData(self._PLUGIN_ACTION_TAG)

            label = text or parts[-1]
            pinned = self._wants_header_slot(top_level_title, parts, action_def)
            if not pinned:
                self._ensure_plugin_separator(current_menu)

            action = QAction(label, self._im.host)
            action.triggered.connect(
                self._make_safe_callback(callback, action_def.get("plugin", "Plugin"))
            )
            if action_def.get("shortcut"):
                action.setShortcut(QKeySequence(action_def["shortcut"]))
            action.setData(self._PLUGIN_ACTION_TAG)
            if pinned:
                self._add_beside_plugin_manager(current_menu, action)
            else:
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

    def _clear_all_plugin_actions(self) -> None:
        """Remove all tagged plugin actions from every menu and the export button.

        The Plugin menu itself is not emptied here — ``_reset_plugin_menu``
        owns that, so both the update and the rebuild path clear it the same
        way.
        """

        def clear_menu(menu: Any) -> None:
            if not menu:
                return
            for act in list(menu.actions()):
                if act.data() == self._PLUGIN_ACTION_TAG:
                    self._remove_action(menu, act)
                elif act.menu():
                    clear_menu(act.menu())

        for menu in self._get_plugin_target_menus():
            clear_menu(menu)
        self._clear_plugin_menubar_entries()

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
                sep = style_menu.addSeparator()
                sep.setData(self._PLUGIN_ACTION_TAG)

            for style_name in self._im.host.plugin_manager.custom_3d_styles:
                if not any(a.text() == style_name for a in style_menu.actions()):
                    action = QAction(style_name, self._im.host)
                    action.setCheckable(True)
                    action.triggered.connect(
                        lambda checked=False,
                        s=style_name: self._im.host.view_3d_manager.set_3d_style(s)
                    )
                    action.setData(self._PLUGIN_ACTION_TAG)
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
                cat = str(p.get("category") or "").strip()
                if cat:
                    categorized.setdefault(cat, []).append(p)
                else:
                    root.append(p)

        # This is the last step of Plugin-menu construction on both the update
        # and the rebuild path, so it also drops the header divider when there
        # turns out to be nothing under it -- an install of only the Plugin
        # Installer would otherwise end on a dangling divider.
        if not categorized and not root:
            self._trim_trailing_separator(plugin_menu)
            return

        # Divide the folder-derived entries below from the context-injected
        # ones ("Plugin/XXX" paths) that add_registered_plugin_actions already
        # placed above: two different organizing schemes should not read as one
        # list. Deferred until we know there is content, so an app with only
        # context-injected plugins gets no dangling trailing divider.
        existing = plugin_menu.actions()
        if existing and not existing[-1].isSeparator():
            plugin_menu.addSeparator()

        for cat in sorted(categorized.keys()):
            current_parent: Any = plugin_menu
            # A folder category arrives with os.sep; a root plugin's own
            # PLUGIN_CATEGORY is hand-written and usually uses "/".
            for part in cat.replace("/", os.sep).split(os.sep):
                if not part:
                    continue
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
            sep.setData(self._PLUGIN_ACTION_TAG)
            for exp in self._im.host.plugin_manager.export_actions:
                a = QAction(exp["label"], self._im.host)
                a.triggered.connect(
                    self._make_safe_callback(
                        exp["callback"], exp.get("plugin", "Plugin")
                    )
                )
                a.setData(self._PLUGIN_ACTION_TAG)
                menu.addAction(a)

    def integrate_plugin_optimization_methods(self) -> None:
        """Inject plugin optimization methods into the 3D Optimization Settings menu.

        Idempotent: previously-added plugin methods are purged first, so a menu
        rebuild (install/uninstall/reload) re-adds only the currently registered
        ones. Without this, _clean_menu removes the menu action while
        opt3d_actions still holds it, and add_optimization_method's dedupe guard
        would then skip re-adding a still-installed method.
        """
        im = self._im
        methods = getattr(im.host.plugin_manager, "optimization_methods", None)
        if not isinstance(methods, dict):
            return

        for key, action in list(im.opt3d_actions.items()):
            if action.data() == self._PLUGIN_ACTION_TAG:
                im.opt_group.removeAction(action)
                im.optimization_menu.removeAction(action)
                del im.opt3d_actions[key]
                im.opt3d_method_labels.pop(key, None)

        for key, entry in methods.items():
            im.add_optimization_method(entry["label"], key)
            if key in im.opt3d_actions:
                im.opt3d_actions[key].setData(self._PLUGIN_ACTION_TAG)

    def integrate_plugin_file_openers(self) -> None:
        """Inject plugin file-opener entries into the File > Import menu."""
        if not self._im.host.plugin_manager.file_openers:
            return

        import_menu = getattr(self._im, "import_menu", None) or getattr(
            self._im.host, "import_menu", None
        )
        if import_menu is None:
            return

        sep = import_menu.addSeparator()
        sep.setData(self._PLUGIN_ACTION_TAG)

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
                            except Exception:
                                logging.exception("Plugin file opener error (%s)", n)
                                return
                            self._im.current_file_path = fpath
                            self._im.host.state_manager.update_window_title()

                return _cb

            a = QAction(f"Import {ext_str} ({p_name})...", self._im.host)
            a.triggered.connect(make_cb(ext_map, filter_str, p_name))
            a.setData(self._PLUGIN_ACTION_TAG)
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
            sep = analysis_menu.addSeparator()
            sep.setData(self._PLUGIN_ACTION_TAG)
            for tool in self._im.host.plugin_manager.analysis_tools:
                a = QAction(
                    f"{tool['label']} ({tool.get('plugin', 'Plugin')})", self._im.host
                )
                a.triggered.connect(
                    self._make_safe_callback(
                        tool["callback"], tool.get("plugin", "Plugin")
                    )
                )
                a.setData(self._PLUGIN_ACTION_TAG)
                analysis_menu.addAction(a)

    def _clear_plugin_ui_elements(self, plugin_menu: QMenu) -> None:
        """Remove only the legacy-tagged plugin actions from the plugin menu."""
        for action in plugin_menu.actions():
            if action.data() == "plugin_action":
                plugin_menu.removeAction(action)

        self._im.plugin_toolbar.clear()  # type: ignore[union-attr]
        self._im.plugin_toolbar.hide()  # type: ignore[union-attr]
