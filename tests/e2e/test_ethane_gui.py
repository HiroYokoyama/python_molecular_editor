# -*- coding: utf-8 -*-
"""E2E GUI tests: real MainWindow, headless, platform-appropriate package.

These tests launch the actual MainWindow (with VTK/PyVista mocked for headless
operation) and interact with it the same way a user would:
  1. Draw ethane in the 2D scene
  2. Click the Convert 2D→3D button
  3. Assert the resulting 3D mol is chemically correct

The conversion uses real RDKit embedding + MMFF optimisation — no mocked results.
"""

import sys
import platform
import pytest
import numpy as np
from PyQt6.QtCore import QPointF, Qt

_OS_LABEL = platform.system()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _draw_ethane(window):
    """Draw ethane (2 C atoms, 1 bond) in the 2D scene programmatically."""
    scene = window.init_manager.scene
    window.ui_manager.set_mode("atom_C")
    c1 = scene.create_atom("C", QPointF(0, 0))
    c2 = scene.create_atom("C", QPointF(50, 0))
    scene.create_bond(scene.atom_items[c1], scene.atom_items[c2])
    return c1, c2


def _assert_ethane_3d(mol):
    """Check the 3D mol has 2 carbons with a chemically reasonable C-C distance."""
    assert mol is not None, f"[{_OS_LABEL}] 3D mol is None after conversion"
    assert mol.GetNumConformers() > 0, f"[{_OS_LABEL}] no 3D conformer"

    carbons = [a for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    assert len(carbons) == 2, f"[{_OS_LABEL}] expected 2 carbons, got {len(carbons)}"

    conf = mol.GetConformer()
    p0 = np.array(conf.GetAtomPosition(carbons[0].GetIdx()))
    p1 = np.array(conf.GetAtomPosition(carbons[1].GetIdx()))
    cc_dist = float(np.linalg.norm(p1 - p0))

    assert 1.3 < cc_dist < 1.7, (
        f"[{_OS_LABEL}] C-C distance {cc_dist:.3f} Å outside chemically valid range [1.3, 1.7]"
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.gui
def test_mainwindow_launches(window):
    """MainWindow initialises without error on this platform."""
    assert window is not None
    assert window.state_manager is not None
    assert window.view_3d_manager is not None
    assert window.init_manager.scene is not None


@pytest.mark.gui
def test_ethane_drawn_in_scene(window):
    """2D ethane can be drawn into the scene on this platform."""
    c1, c2 = _draw_ethane(window)

    assert len(window.state_manager.data.atoms) == 2
    assert len(window.state_manager.data.bonds) == 1
    symbols = [v["symbol"] for v in window.state_manager.data.atoms.values()]
    assert symbols.count("C") == 2


@pytest.mark.gui
def test_convert_button_enabled_after_draw(window):
    """Convert 2D→3D button is enabled once atoms are present."""
    _draw_ethane(window)
    assert window.init_manager.convert_button.isEnabled()


@pytest.mark.gui
def test_ethane_2d_to_3d_via_button(window, qtbot):
    """Clicking Convert produces a real 3D mol with correct C-C geometry."""
    _draw_ethane(window)

    qtbot.mouseClick(window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None, timeout=5000
    )

    _assert_ethane_3d(window.view_3d_manager.current_mol)


@pytest.mark.gui
def test_3d_buttons_enabled_after_conversion(window, qtbot):
    """Export and Optimize buttons become active after a successful conversion."""
    _draw_ethane(window)
    qtbot.mouseClick(window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None, timeout=5000
    )

    assert window.init_manager.optimize_3d_button.isEnabled()
    assert window.init_manager.export_button.isEnabled()
    assert window.init_manager.analysis_action.isEnabled()


@pytest.mark.gui
def test_clear_resets_3d_mol(window, qtbot):
    """Clear All resets current_mol back to None."""
    _draw_ethane(window)
    qtbot.mouseClick(window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None, timeout=5000
    )

    window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    window.state_manager.has_unsaved_changes = False

    assert len(window.state_manager.data.atoms) == 0
    assert window.view_3d_manager.current_mol is None
