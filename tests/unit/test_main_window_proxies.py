"""
Unit tests for MainWindow (ui/main_window.py).

Covers:
  - All manager attributes are assigned during __init__
  - _is_restoring_state initial value
  - start_calculation signal exists
  - current_mol getter / setter proxy → view_3d_manager
  - plotter property proxy → view_3d_manager
  - data property proxy → state_manager
  - scene property proxy → init_manager
  - draw_molecule_3d proxy: sets current_mol and delegates to view_3d_manager
"""

from unittest.mock import patch, MagicMock


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_window():
    """Instantiate MainWindow with all heavy managers patched."""
    patches = [
        "moleditpy.ui.main_window.ExportManager",
        "moleditpy.ui.main_window.View3DManager",
        "moleditpy.ui.main_window.Edit3DManager",
        "moleditpy.ui.main_window.EditActionsManager",
        "moleditpy.ui.main_window.ComputeManager",
        "moleditpy.ui.main_window.DialogManager",
        "moleditpy.ui.main_window.IOManager",
        "moleditpy.ui.main_window.StateManager",
        "moleditpy.ui.main_window.StringImporterManager",
        "moleditpy.ui.main_window.UIManager",
        "moleditpy.ui.main_window.MainInitManager",
    ]
    with patch.multiple("moleditpy.ui.main_window", **{
        p.split(".")[-1]: MagicMock() for p in patches
    }):
        from moleditpy.ui.main_window import MainWindow
        # Patch all at module level
        with patch("moleditpy.ui.main_window.ExportManager"), \
             patch("moleditpy.ui.main_window.View3DManager"), \
             patch("moleditpy.ui.main_window.Edit3DManager"), \
             patch("moleditpy.ui.main_window.EditActionsManager"), \
             patch("moleditpy.ui.main_window.ComputeManager"), \
             patch("moleditpy.ui.main_window.DialogManager"), \
             patch("moleditpy.ui.main_window.IOManager"), \
             patch("moleditpy.ui.main_window.StateManager"), \
             patch("moleditpy.ui.main_window.StringImporterManager"), \
             patch("moleditpy.ui.main_window.UIManager"), \
             patch("moleditpy.ui.main_window.MainInitManager"):
            mw = MainWindow()
    return mw


# ---------------------------------------------------------------------------
# Instantiation
# ---------------------------------------------------------------------------

def test_mainwindow_all_managers_assigned(app):
    mw = _make_window()
    assert hasattr(mw, "export_manager")
    assert hasattr(mw, "view_3d_manager")
    assert hasattr(mw, "edit_3d_manager")
    assert hasattr(mw, "edit_actions_manager")
    assert hasattr(mw, "compute_manager")
    assert hasattr(mw, "dialog_manager")
    assert hasattr(mw, "io_manager")
    assert hasattr(mw, "state_manager")
    assert hasattr(mw, "string_importer_manager")
    assert hasattr(mw, "ui_manager")
    assert hasattr(mw, "init_manager")


def test_mainwindow_is_restoring_state_default(app):
    mw = _make_window()
    assert mw._is_restoring_state is False


def test_mainwindow_start_calculation_signal_exists(app):
    from moleditpy.ui.main_window import MainWindow
    assert hasattr(MainWindow, "start_calculation")


# ---------------------------------------------------------------------------
# Proxy properties
# ---------------------------------------------------------------------------

def test_current_mol_getter_delegates_to_view_3d_manager(app):
    mw = _make_window()
    mock_mol = MagicMock()
    mw.view_3d_manager.current_mol = mock_mol
    assert mw.current_mol is mock_mol


def test_current_mol_setter_delegates_to_view_3d_manager(app):
    mw = _make_window()
    mock_mol = MagicMock()
    mw.current_mol = mock_mol
    assert mw.view_3d_manager.current_mol is mock_mol


def test_plotter_property_delegates_to_view_3d_manager(app):
    mw = _make_window()
    mock_plotter = MagicMock()
    mw.view_3d_manager.plotter = mock_plotter
    assert mw.plotter is mock_plotter


def test_data_property_delegates_to_state_manager(app):
    mw = _make_window()
    mock_data = MagicMock()
    mw.state_manager.data = mock_data
    assert mw.data is mock_data


def test_scene_property_delegates_to_init_manager(app):
    mw = _make_window()
    mock_scene = MagicMock()
    mw.init_manager.scene = mock_scene
    assert mw.scene is mock_scene


# ---------------------------------------------------------------------------
# draw_molecule_3d proxy
# ---------------------------------------------------------------------------

def test_draw_molecule_3d_sets_current_mol(app):
    mw = _make_window()
    mock_mol = MagicMock()
    mw.draw_molecule_3d(mock_mol)
    assert mw.view_3d_manager.current_mol is mock_mol


def test_draw_molecule_3d_delegates_to_view_3d_manager(app):
    mw = _make_window()
    mock_mol = MagicMock()
    mw.draw_molecule_3d(mock_mol)
    mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(mock_mol)


def test_draw_molecule_3d_none_mol(app):
    mw = _make_window()
    mw.draw_molecule_3d(None)
    assert mw.view_3d_manager.current_mol is None
    mw.view_3d_manager.draw_molecule_3d.assert_called_once_with(None)
