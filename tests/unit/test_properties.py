import pytest
from rdkit import Chem
from moleditpy.modules.analysis_window import AnalysisWindow
from unittest.mock import MagicMock
from PyQt6.QtWidgets import QApplication, QWidget


def test_analysis_window_regular_mol(qtbot):
    """Verify AnalysisWindow calculates and displays properties for regular molecules."""
    # Simple Ethanol
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)

    # Mock parent using a real QWidget to satisfy PyQt6 type check
    parent = QWidget()
    # We can still add attributes to it if needed
    parent.statusBar = MagicMock()
    parent.statusBar.return_value = MagicMock()

    try:
        window = AnalysisWindow(mol, parent=parent)
    except Exception as e:
        print(f"DEBUG: AnalysisWindow init failed with: {e}")
        raise e

    qtbot.addWidget(window)

    # Extract values from QLineEdit widgets
    # Formula should be C2H6O
    formula_found = False
    for i in range(window.layout().count()):
        item = window.layout().itemAt(i)
        if hasattr(item, "layout") and item.layout():
            grid = item.layout()
            for r in range(grid.rowCount()):
                label_item = grid.itemAtPosition(r, 0)
                if label_item and label_item.widget():
                    label = label_item.widget().text()
                    if "Formula" in label:
                        value = grid.itemAtPosition(r, 1).widget().text()
                        assert "C2H6O" in value
                        formula_found = True
                        break
    assert formula_found


def test_analysis_window_xyz_derived(qtbot):
    """Verify AnalysisWindow uses manual logic for XYZ-derived structures."""
    mol = Chem.MolFromSmiles("CCO")  # Placeholder
    # Manually attach XYZ metadata
    mol._xyz_atom_data = [
        ("C", 0, 0, 0),
        ("C", 1.5, 0, 0),
        ("O", 2.0, 1.0, 0),
        ("H", 1.5, -0.5, 0),  # Partial hydrogen
    ]

    parent = QWidget()
    window = AnalysisWindow(mol, parent=parent, is_xyz_derived=True)
    qtbot.addWidget(window)

    # Formula should be C2HO (based on manual atom counts, not RDKit)
    formula_val = ""
    smiles_present = False
    for i in range(window.layout().count()):
        item = window.layout().itemAt(i)
        if hasattr(item, "layout") and item.layout():
            grid = item.layout()
            for r in range(grid.rowCount()):
                label_item = grid.itemAtPosition(r, 0)
                if label_item and label_item.widget():
                    label = label_item.widget().text()
                    if "Formula" in label:
                        formula_val = grid.itemAtPosition(r, 1).widget().text()
                    if "SMILES" in label:
                        smiles_present = True
        elif item.widget() and "SMILES and structure-dependent" in item.widget().text():
            # Check for the warning message
            pass

    assert "C2HO" in formula_val
    assert not smiles_present  # SMILES should be withheld for XYZ
