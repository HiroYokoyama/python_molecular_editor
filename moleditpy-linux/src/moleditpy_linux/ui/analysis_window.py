#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QGridLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rd_inchi


class AnalysisWindow(QDialog):
    def __init__(self, mol, parent=None, is_xyz_derived=False):
        super().__init__(parent)
        self.mol = mol
        self.is_xyz_derived = is_xyz_derived  # Flag indicating if derived from XYZ
        self.setWindowTitle("Molecule Analysis")
        self.setMinimumWidth(400)
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)
        grid_layout = QGridLayout()

        # --- Calculate molecular properties ---
        try:
            # Import RDKit modules

            if self.is_xyz_derived:
                # XYZ-derived: Calculate directly from original atomic information
                # (Avoids inaccuracies from bond estimation)

                # Retrieve original atomic data from XYZ file
                if hasattr(self.mol, "_xyz_atom_data"):
                    xyz_atoms = self.mol._xyz_atom_data
                else:
                    # Fallback: Retrieve from RDKit object
                    xyz_atoms = [
                        (atom.GetSymbol(), 0, 0, 0) for atom in self.mol.GetAtoms()
                    ]

                # Summarize atom counts and element types
                atom_counts = {}
                total_atoms = len(xyz_atoms)
                num_heavy_atoms = 0

                for symbol, x, y, z in xyz_atoms:
                    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                    if symbol != "H":  # Non-hydrogen
                        num_heavy_atoms += 1

                # Construct molecular formula manually (following Hill system)
                element_order = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
                formula_parts = []

                # Add elements in defined order
                remaining_counts = atom_counts.copy()
                for element in element_order:
                    if element in remaining_counts:
                        count = remaining_counts[element]
                        if count == 1:
                            formula_parts.append(element)
                        else:
                            formula_parts.append(f"{element}{count}")
                        del remaining_counts[element]

                # Add remaining elements alphabetically
                for element in sorted(remaining_counts.keys()):
                    count = remaining_counts[element]
                    if count == 1:
                        formula_parts.append(element)
                    else:
                        formula_parts.append(f"{element}{count}")

                mol_formula = "".join(formula_parts)

                # Get molecular weight and exact mass from RDKit

                mol_wt = 0.0
                exact_mw = 0.0
                pt = Chem.GetPeriodicTable()

                for symbol, count in atom_counts.items():
                    try:
                        # Get atomic weights and isotope masses from RDKit Periodic Table
                        atomic_num = pt.GetAtomicNumber(symbol)
                        atomic_weight = pt.GetAtomicWeight(atomic_num)
                        exact_mass = pt.GetMostCommonIsotopeMass(atomic_num)

                        mol_wt += atomic_weight * count
                        exact_mw += exact_mass * count
                    except (ValueError, RuntimeError):
                        # Skip unrecognized elements
                        print(
                            f"Warning: Unknown element {symbol}, skipping in mass calculation"
                        )
                        continue

                # Compile properties for display (derived from XYZ data)
                properties = {
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Total Atoms:": str(total_atoms),
                }

                # Add a note about XYZ limitations
                note_label = QLabel(
                    "<i>Note: SMILES and structure-dependent properties are not available for XYZ-derived structures due to potential bond estimation inaccuracies.</i>"
                )
                note_label.setWordWrap(True)
                main_layout.addWidget(note_label)

            else:
                # Regular molecule (from MOL file or 2D editor): Calculate all properties

                # Create temporary H-depleted molecule for SMILES generation
                mol_for_smiles = Chem.RemoveHs(self.mol)
                # Generate canonical SMILES from H-depleted molecule
                smiles = Chem.MolToSmiles(mol_for_smiles, isomericSmiles=True)

                # Calculate various properties
                mol_formula = rdMolDescriptors.CalcMolFormula(self.mol)
                mol_wt = Descriptors.MolWt(self.mol)
                exact_mw = Descriptors.ExactMolWt(self.mol)
                num_heavy_atoms = self.mol.GetNumHeavyAtoms()
                num_rings = rdMolDescriptors.CalcNumRings(self.mol)
                log_p = Descriptors.MolLogP(self.mol)
                tpsa = Descriptors.TPSA(self.mol)
                num_h_donors = rdMolDescriptors.CalcNumHBD(self.mol)
                num_h_acceptors = rdMolDescriptors.CalcNumHBA(self.mol)

                # Generate InChI
                try:
                    inchi = Chem.MolToInchi(self.mol)
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    inchi = "N/A"

                # Generate InChIKey (with fallback)
                try:
                    # Prefer Chem.MolToInchiKey when available
                    inchi_key = None
                    try:
                        inchi_key = Chem.MolToInchiKey(self.mol)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        # Fallback to rdkit.Chem.inchi if present
                        try:
                            inchi_key = rd_inchi.MolToInchiKey(self.mol)
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            inchi_key = None

                    if not inchi_key:
                        inchi_key = "N/A"
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    inchi_key = "N/A"

                # Compile properties for display
                properties = {
                    "SMILES:": smiles,
                    "InChI:": inchi,
                    "InChIKey:": inchi_key,
                    "Molecular Formula:": mol_formula,
                    "Molecular Weight:": f"{mol_wt:.4f}",
                    "Exact Mass:": f"{exact_mw:.4f}",
                    "Heavy Atoms:": str(num_heavy_atoms),
                    "Ring Count:": str(num_rings),
                    "LogP (o/w):": f"{log_p:.3f}",
                    "TPSA (Å²):": f"{tpsa:.2f}",
                    "H-Bond Donors:": str(num_h_donors),
                    "H-Bond Acceptors:": str(num_h_acceptors),
                }
        except (AttributeError, RuntimeError, ValueError) as e:
            main_layout.addWidget(QLabel(f"Error calculating properties: {e}"))
            return

        # --- Display results in UI ---
        row = 0
        for label_text, value_text in properties.items():
            label = QLabel(f"<b>{label_text}</b>")
            value = QLineEdit(value_text)
            value.setReadOnly(True)

            copy_btn = QPushButton("Copy")
            copy_btn.clicked.connect(
                lambda _, v=value: self.copy_to_clipboard(v.text())
            )

            grid_layout.addWidget(label, row, 0)
            grid_layout.addWidget(value, row, 1)
            grid_layout.addWidget(copy_btn, row, 2)
            row += 1

        main_layout.addLayout(grid_layout)

        # --- OK Button ---
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        main_layout.addWidget(ok_button, 0, Qt.AlignmentFlag.AlignCenter)

        self.setLayout(main_layout)

    def copy_to_clipboard(self, text):
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        if self.parent() and hasattr(self.parent(), "statusBar"):
            self.parent().statusBar().showMessage(
                f"Copied '{text}' to clipboard.", 2000
            )
