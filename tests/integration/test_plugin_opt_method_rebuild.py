"""Regression test: plugin optimization methods survive a menu rebuild.

On rebuild_plugin_menus, _clean_menu removes the plugin-managed optimization
action from the menu but leaves it in opt3d_actions. Before the fix,
add_optimization_method's dedupe guard then skipped re-adding a still-installed
method, so plugin optimizers vanished from the 3D Optimization Settings menu
after any plugin install/uninstall/reload.
"""

from types import SimpleNamespace

from PyQt6.QtCore import QObject
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QMenu

from moleditpy.ui.main_window_init import MainInitManager
from moleditpy.ui.plugin_menu_manager import PluginMenuManager


def _make_im(methods):
    im = MainInitManager.__new__(MainInitManager)
    im.host = QObject()
    im.host.plugin_manager = SimpleNamespace(optimization_methods=methods)
    im.settings = {"optimization_method": "MMFF_RDKIT"}
    im.optimization_menu = QMenu()
    im.opt3d_separator = im.optimization_menu.addSeparator()
    im.opt_group = QActionGroup(im.optimization_menu)
    im.opt_group.setExclusive(True)
    im.opt3d_method_labels = {}

    mmff = QAction("MMFF94s (RDKit)")
    mmff.setCheckable(True)
    im.opt_group.addAction(mmff)
    im.opt3d_actions = {"MMFF_RDKIT": mmff}
    return im


def _menu_labels(im):
    return {a.text() for a in im.optimization_menu.actions() if not a.isSeparator()}


def test_plugin_method_readded_after_rebuild(app):
    methods = {
        "QUICK UFF": {"plugin": "P", "callback": lambda m: True, "label": "Quick UFF"}
    }
    im = _make_im(methods)
    pmm = PluginMenuManager(im)

    pmm.integrate_plugin_optimization_methods()
    assert "Quick UFF" in _menu_labels(im)

    # Simulate _clean_menu during rebuild: the tagged action is pulled from the
    # menu but not from opt3d_actions.
    im.optimization_menu.removeAction(im.opt3d_actions["QUICK UFF"])
    assert "Quick UFF" not in _menu_labels(im)

    pmm.integrate_plugin_optimization_methods()
    assert "Quick UFF" in _menu_labels(im)
    assert "QUICK UFF" in im.opt3d_actions


def test_uninstalled_plugin_method_purged_on_rebuild(app):
    methods = {
        "QUICK UFF": {"plugin": "P", "callback": lambda m: True, "label": "Quick UFF"}
    }
    im = _make_im(methods)
    pmm = PluginMenuManager(im)
    pmm.integrate_plugin_optimization_methods()

    # Plugin uninstalled: registry empties before the rebuild.
    methods.clear()
    pmm.integrate_plugin_optimization_methods()

    assert "QUICK UFF" not in im.opt3d_actions
    assert "Quick UFF" not in _menu_labels(im)
