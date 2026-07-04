"""Regression tests for restoring a persisted plugin optimization method.

Plugins register their optimization methods *after* the settings menu is built,
so a saved default naming a plugin method could not be restored at menu-init and
silently fell back to MMFF_RDKIT. add_optimization_method must reconcile the
newly added method with the saved default.
"""

from PyQt6.QtCore import QObject
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QMenu

from moleditpy.ui.main_window_init import MainInitManager


def _make_init_manager(saved_opt):
    """Build a MainInitManager with only the state add_optimization_method needs,
    simulating the post-menu-init state where MMFF_RDKIT is the checked fallback.
    """
    im = MainInitManager.__new__(MainInitManager)
    # Real QObject so QAction accepts it as parent (PyQt6 validates parent type);
    # the triggered lambda is only stored, never invoked here.
    im.host = QObject()
    im.settings = {"optimization_method": saved_opt}
    im.optimization_menu = QMenu()
    im.opt3d_separator = im.optimization_menu.addSeparator()
    im.opt_group = QActionGroup(im.optimization_menu)
    im.opt_group.setExclusive(True)
    im.opt3d_method_labels = {}

    mmff = QAction("MMFF94s (RDKit)")
    mmff.setCheckable(True)
    mmff.setChecked(True)
    im.opt_group.addAction(mmff)
    im.opt3d_actions = {"MMFF_RDKIT": mmff}
    im.optimization_method = "MMFF_RDKIT"
    return im


def test_saved_plugin_method_restored_on_registration(app):
    im = _make_init_manager(saved_opt="QUICK UFF")

    im.add_optimization_method("Quick UFF", "QUICK UFF")

    assert im.opt3d_actions["QUICK UFF"].isChecked()
    assert not im.opt3d_actions["MMFF_RDKIT"].isChecked()
    assert im.optimization_method == "QUICK UFF"


def test_non_default_plugin_method_left_unchecked(app):
    im = _make_init_manager(saved_opt="MMFF_RDKIT")

    im.add_optimization_method("Quick UFF", "QUICK UFF")

    assert not im.opt3d_actions["QUICK UFF"].isChecked()
    assert im.opt3d_actions["MMFF_RDKIT"].isChecked()
    assert im.optimization_method == "MMFF_RDKIT"
