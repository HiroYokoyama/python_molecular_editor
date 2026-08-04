# -*- coding: utf-8 -*-
"""2D->3D conversion of benzene, driven through the running application.

The real ``CalculationWorker`` QThread does the work here (nothing is
monkeypatched to run synchronously), and the result is rendered into the real
PyVista widget, so these tests cover the whole button-to-geometry path.
"""

import importlib
import math

import numpy as np
import pytest
from PyQt6.QtCore import Qt

from benzene import assert_is_2d_benzene, draw_benzene_by_hand, place_benzene_template
from conftest import PKG

pytestmark = pytest.mark.full_gui

CONVERT_TIMEOUT_MS = 120000

# Aromatic C-C is 1.39 A; MMFF-relaxed benzene stays well inside this window.
CC_MIN, CC_MAX = 1.32, 1.46
CH_MIN, CH_MAX = 0.95, 1.15


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _convert_via_button(window, qtbot):
    """Click Convert 2D->3D and wait for the worker thread to deliver a mol."""
    assert window.init_manager.convert_button.isEnabled(), (
        "Convert button is disabled although the canvas holds atoms"
    )
    qtbot.mouseClick(window.init_manager.convert_button, Qt.MouseButton.LeftButton)
    qtbot.waitUntil(
        lambda: window.view_3d_manager.current_mol is not None,
        timeout=CONVERT_TIMEOUT_MS,
    )
    return window.view_3d_manager.current_mol


def _ring_carbons(mol):
    return [a for a in mol.GetAtoms() if a.GetSymbol() == "C"]


def _pos(mol, idx):
    return np.array(mol.GetConformer().GetAtomPosition(idx))


def assert_is_3d_benzene(mol, label="convert"):
    """Full geometric check: composition, bond lengths, angles, planarity."""
    assert mol is not None, f"[{label}] conversion produced no molecule"
    assert mol.GetNumConformers() > 0, f"[{label}] molecule has no 3D conformer"

    carbons = _ring_carbons(mol)
    hydrogens = [a for a in mol.GetAtoms() if a.GetSymbol() == "H"]
    assert len(carbons) == 6, f"[{label}] expected 6 C, got {len(carbons)}"
    assert len(hydrogens) == 6, f"[{label}] expected 6 H, got {len(hydrogens)}"

    # Ring C-C bond lengths
    cc = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetSymbol() == "C" and b.GetSymbol() == "C":
            cc.append(
                float(np.linalg.norm(_pos(mol, a.GetIdx()) - _pos(mol, b.GetIdx())))
            )
    assert len(cc) == 6, f"[{label}] expected 6 C-C bonds, got {len(cc)}"
    assert all(CC_MIN < d < CC_MAX for d in cc), (
        f"[{label}] C-C lengths outside [{CC_MIN}, {CC_MAX}]: "
        f"{[round(d, 3) for d in cc]}"
    )

    # C-H bond lengths
    ch = []
    for bond in mol.GetBonds():
        syms = {bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()}
        if syms == {"C", "H"}:
            ch.append(
                float(
                    np.linalg.norm(
                        _pos(mol, bond.GetBeginAtomIdx())
                        - _pos(mol, bond.GetEndAtomIdx())
                    )
                )
            )
    assert len(ch) == 6, f"[{label}] expected 6 C-H bonds, got {len(ch)}"
    assert all(CH_MIN < d < CH_MAX for d in ch), (
        f"[{label}] C-H lengths outside [{CH_MIN}, {CH_MAX}]: "
        f"{[round(d, 3) for d in ch]}"
    )

    # Planarity: every carbon within 0.1 A of the best-fit plane of the ring.
    coords = np.array([_pos(mol, a.GetIdx()) for a in carbons])
    centred = coords - coords.mean(axis=0)
    normal = np.linalg.svd(centred)[2][2]
    deviation = np.abs(centred @ normal)
    assert deviation.max() < 0.1, (
        f"[{label}] ring is not planar; max out-of-plane deviation "
        f"{deviation.max():.3f} A"
    )

    # Internal ring angles must sit near the 120 deg of an sp2 hexagon.
    order = _ring_order(mol, carbons)
    for i in range(6):
        p_prev = _pos(mol, order[i - 1])
        p = _pos(mol, order[i])
        p_next = _pos(mol, order[(i + 1) % 6])
        v1, v2 = p_prev - p, p_next - p
        cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        angle = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
        assert 112.0 < angle < 128.0, (
            f"[{label}] ring angle at atom {order[i]} is {angle:.1f} deg, not sp2-like"
        )


def _ring_order(mol, carbons):
    """Walk the carbon ring so consecutive indices are actually bonded."""
    idxs = {a.GetIdx() for a in carbons}
    neighbours = {
        i: [
            n.GetIdx()
            for n in mol.GetAtomWithIdx(i).GetNeighbors()
            if n.GetIdx() in idxs
        ]
        for i in idxs
    }
    start = next(iter(idxs))
    order = [start, neighbours[start][0]]
    while len(order) < 6:
        nxt = [n for n in neighbours[order[-1]] if n != order[-2]]
        assert nxt, "carbon ring is not a closed cycle"
        order.append(nxt[0])
    return order


# ---------------------------------------------------------------------------
# Conversion through the GUI
# ---------------------------------------------------------------------------


def test_convert_button_enabled_after_drawing_benzene(full_window, qtbot):
    place_benzene_template(full_window, qtbot)
    assert full_window.init_manager.convert_button.isEnabled()


def test_benzene_2d_to_3d_via_button(full_window, qtbot):
    """The headline path: draw benzene, press the button, get real geometry."""
    place_benzene_template(full_window, qtbot)
    assert_is_2d_benzene(full_window.state_manager.data)

    mol = _convert_via_button(full_window, qtbot)
    assert_is_3d_benzene(mol, "button")


def test_converted_benzene_is_still_aromatic(full_window, qtbot):
    from rdkit import Chem

    draw_benzene_by_hand(full_window)
    mol = _convert_via_button(full_window, qtbot)

    smiles = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))
    assert smiles == "c1ccccc1", f"3D molecule is no longer benzene: {smiles}"


def test_3d_actions_enable_after_conversion(full_window, qtbot):
    draw_benzene_by_hand(full_window)
    _convert_via_button(full_window, qtbot)

    assert full_window.init_manager.optimize_3d_button.isEnabled()
    assert full_window.init_manager.export_button.isEnabled()
    assert full_window.init_manager.analysis_action.isEnabled()


def test_conversion_populates_the_real_vtk_scene(full_window, qtbot):
    """Geometry actually reaches the renderer -- actors exist and render."""
    draw_benzene_by_hand(full_window)
    _convert_via_button(full_window, qtbot)

    plotter = full_window.view_3d_manager.plotter
    actors = plotter.renderer.actors
    assert actors, "no actors in the 3D renderer after conversion"
    plotter.render()


def test_optimize_3d_keeps_benzene_valid(full_window, qtbot):
    """Running the 3D optimiser on the result must not wreck the ring."""
    draw_benzene_by_hand(full_window)
    first = _convert_via_button(full_window, qtbot)
    before = np.array(first.GetConformer().GetPositions())

    qtbot.mouseClick(
        full_window.init_manager.optimize_3d_button, Qt.MouseButton.LeftButton
    )
    qtbot.waitUntil(
        lambda: full_window.view_3d_manager.current_mol is not None
        and not np.array_equal(
            np.array(
                full_window.view_3d_manager.current_mol.GetConformer().GetPositions()
            ),
            before,
        ),
        timeout=CONVERT_TIMEOUT_MS,
    )
    assert_is_3d_benzene(full_window.view_3d_manager.current_mol, "optimize")


def test_clear_all_drops_the_3d_molecule(full_window, qtbot):
    draw_benzene_by_hand(full_window)
    _convert_via_button(full_window, qtbot)

    full_window.edit_actions_manager.clear_2d_editor(push_to_undo=False)
    full_window.state_manager.has_unsaved_changes = False

    assert len(full_window.state_manager.data.atoms) == 0
    assert full_window.view_3d_manager.current_mol is None


# ---------------------------------------------------------------------------
# Conversion back-ends
# ---------------------------------------------------------------------------


def _run_worker_sync(mol_block, mode):
    CalculationWorker = importlib.import_module(
        f"{PKG}.ui.calculation_worker"
    ).CalculationWorker
    results = []
    worker = CalculationWorker()
    worker.finished.connect(results.append)
    worker.run_calculation(mol_block, {"conversion_mode": mode, "do_optimize": True})
    return results[0][1] if results and results[0] else None


@pytest.mark.parametrize("mode", ["rdkit", "direct", "fallback"])
def test_benzene_conversion_backends(full_window, mode):
    """Every always-available conversion back-end yields the same molecule."""
    draw_benzene_by_hand(full_window)
    mol_block = full_window.state_manager.data.to_mol_block()
    assert_is_3d_benzene(_run_worker_sync(mol_block, mode), mode)


def test_benzene_conversion_obabel_when_available(full_window):
    """Open Babel back-end, where the platform actually ships it (Linux)."""
    pkg = importlib.import_module(PKG)
    if not getattr(pkg, "OBABEL_AVAILABLE", False):
        pytest.skip("Open Babel is not available on this platform")

    draw_benzene_by_hand(full_window)
    mol_block = full_window.state_manager.data.to_mol_block()
    assert_is_3d_benzene(_run_worker_sync(mol_block, "obabel"), "obabel")


# ---------------------------------------------------------------------------
# Round-trip out of the app
# ---------------------------------------------------------------------------


def test_converted_benzene_exports_to_xyz(full_window, qtbot, tmp_path):
    """The 3D result written as XYZ reads back as the same six-carbon ring."""
    draw_benzene_by_hand(full_window)
    mol = _convert_via_button(full_window, qtbot)

    from rdkit import Chem

    target = tmp_path / "benzene.xyz"
    Chem.MolToXYZFile(mol, str(target))

    lines = target.read_text(encoding="utf-8").splitlines()
    assert int(lines[0].strip()) == 12, "XYZ atom count is not 12"
    symbols = [ln.split()[0] for ln in lines[2:14]]
    assert symbols.count("C") == 6 and symbols.count("H") == 6
