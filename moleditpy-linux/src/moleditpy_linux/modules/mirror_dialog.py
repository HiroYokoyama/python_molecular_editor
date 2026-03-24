#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import (
    QButtonGroup,
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
)
from rdkit import Chem
from rdkit.Geometry import Point3D


class MirrorDialog(QDialog):  # pragma: no cover
    """Dialog to create a mirror image of the molecule."""

    def __init__(self, mol, main_window, parent=None):
        super().__init__(parent)
        self.mol = mol
        self.main_window = main_window
        self.plane_group = None
        self.xy_radio = None
        self.xz_radio = None
        self.yz_radio = None
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Mirror Molecule")
        self.setMinimumSize(300, 200)

        layout = QVBoxLayout(self)

        # Instructional text
        info_label = QLabel("Select the mirror plane to create molecular mirror image:")
        layout.addWidget(info_label)

        # Radio buttons for mirror plane selection
        self.plane_group = QButtonGroup(self)

        self.xy_radio = QRadioButton("XY plane (Z = 0)")
        self.xz_radio = QRadioButton("XZ plane (Y = 0)")
        self.yz_radio = QRadioButton("YZ plane (X = 0)")

        self.xy_radio.setChecked(True)  # Default selection

        self.plane_group.addButton(self.xy_radio, 0)
        self.plane_group.addButton(self.xz_radio, 1)
        self.plane_group.addButton(self.yz_radio, 2)

        layout.addWidget(self.xy_radio)
        layout.addWidget(self.xz_radio)
        layout.addWidget(self.yz_radio)

        layout.addSpacing(20)

        # Buttons
        button_layout = QHBoxLayout()

        apply_button = QPushButton("Apply Mirror")
        apply_button.clicked.connect(self.apply_mirror)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)

        button_layout.addWidget(apply_button)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

    def apply_mirror(self):
        """Apply mirror transformation across the selected plane."""
        if not self.mol or self.mol.GetNumConformers() == 0:
            QMessageBox.warning(self, "Error", "No 3D coordinates available.")
            return

        # Get the selected plane
        plane_id = self.plane_group.checkedId()

        try:
            conf = self.mol.GetConformer()

            # Transform coordinates for each atom
            for atom_idx in range(self.mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(atom_idx)

                new_pos = [pos.x, pos.y, pos.z]
                if plane_id == 0:  # XY plane (mirror across Z-axis)
                    new_pos = [pos.x, pos.y, -pos.z]
                elif plane_id == 1:  # XZ plane (mirror across Y-axis)
                    new_pos = [pos.x, -pos.y, pos.z]
                elif plane_id == 2:  # YZ plane (mirror across X-axis)
                    new_pos = [-pos.x, pos.y, pos.z]

                # Set new coordinates
                conf.SetAtomPosition(
                    atom_idx, Point3D(new_pos[0], new_pos[1], new_pos[2])
                )

            # Force recalculation of chiral tags after mirror transform (required before 3D rendering)
            try:
                if self.mol.GetNumConformers() > 0:
                    # Clear existing chiral tags
                    for atom in self.mol.GetAtoms():
                        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                    # Calculate new chiral tags from 3D coordinates
                    Chem.AssignAtomChiralTagsFromStructure(self.mol, confId=0)
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                print(f"Error updating chiral tags: {e}")

            # Update 3D view (which also draws 3D chiral labels)
            self.main_window.draw_molecule_3d(self.mol)

            # Update 2D chiral labels
            self.main_window.update_chiral_labels()

            self.main_window.push_undo_state()

            plane_names = ["XY", "XZ", "YZ"]
            self.main_window.statusBar().showMessage(
                f"Molecule mirrored across {plane_names[plane_id]} plane."
            )

        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            QMessageBox.critical(
                self, "Error", f"Failed to apply mirror transformation: {str(e)}"
            )
