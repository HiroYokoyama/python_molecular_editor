#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

"""
main_window_dialog_manager.py
Module separated from MainWindow (main_window.py)
Functional class: MainWindowDialogManager
"""

import json
import os

# PyQt6 Modules
from PyQt6.QtCore import QDateTime
from PyQt6.QtWidgets import QInputDialog, QMessageBox

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except (AttributeError, RuntimeError, TypeError):
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .about_dialog import AboutDialog
    from .align_plane_dialog import AlignPlaneDialog
    from .alignment_dialog import AlignmentDialog
    from .analysis_window import AnalysisWindow
    from .angle_dialog import AngleDialog
    from .bond_length_dialog import BondLengthDialog
    from .constants import VERSION
    from .constrained_optimization_dialog import ConstrainedOptimizationDialog
    from .dihedral_dialog import DihedralDialog
    from .mirror_dialog import MirrorDialog
    from .move_group_dialog import MoveGroupDialog
    from .periodic_table_dialog import PeriodicTableDialog
    from .planarize_dialog import PlanarizeDialog
    from .translation_dialog import TranslationDialog
    from .user_template_dialog import UserTemplateDialog
except ImportError:
    # Fallback to absolute imports for script-style execution
    from modules.about_dialog import AboutDialog
    from modules.align_plane_dialog import AlignPlaneDialog
    from modules.alignment_dialog import AlignmentDialog
    from modules.analysis_window import AnalysisWindow
    from modules.angle_dialog import AngleDialog
    from modules.bond_length_dialog import BondLengthDialog
    from modules.constants import VERSION
    from modules.constrained_optimization_dialog import ConstrainedOptimizationDialog
    from modules.dihedral_dialog import DihedralDialog
    from modules.mirror_dialog import MirrorDialog
    from modules.move_group_dialog import MoveGroupDialog
    from modules.periodic_table_dialog import PeriodicTableDialog
    from modules.planarize_dialog import PlanarizeDialog
    from modules.translation_dialog import TranslationDialog
    from modules.user_template_dialog import UserTemplateDialog


# --- Class Definition ---
class MainWindowDialogManager:
    """Mixin class separated from main_window.py"""

    _cls = None

    def show_about_dialog(self):
        """Show the custom About dialog with Easter egg functionality"""
        dialog = AboutDialog(self, self)
        dialog.exec()

    def open_periodic_table_dialog(self):
        dialog = PeriodicTableDialog(self)
        dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action = self.tool_group.checkedAction()
        if checked_action:
            self.tool_group.setExclusive(False)
            checked_action.setChecked(False)
            self.tool_group.setExclusive(True)
        dialog.exec()

    def open_analysis_window(self):
        if self.current_mol:
            dialog = AnalysisWindow(
                self.current_mol, self, is_xyz_derived=self.is_xyz_derived
            )
            dialog.exec()
        else:
            self.statusBar().showMessage(
                "Please generate a 3D structure first to show analysis."
            )

    def open_template_dialog(self):
        """Open the template dialog"""
        dialog = UserTemplateDialog(self, self)
        dialog.exec()

    def open_template_dialog_and_activate(self):
        """Open the template dialog and activate the selected template for use in the main window"""
        # Check for existing dialog
        _template_dialog = getattr(self, "_template_dialog", None)
        if _template_dialog and not _template_dialog.isHidden():
            # Bring existing dialog to front
            _template_dialog.raise_()
            _template_dialog.activateWindow()
            return

        # Create new dialog
        self._template_dialog = UserTemplateDialog(self, self)
        self._template_dialog.show()

        # Activate if a template is selected after dialog is closed
        def on_dialog_finished():
            if (
                hasattr(self._template_dialog, "selected_template")
                and self._template_dialog.selected_template
            ):
                template_name = self._template_dialog.selected_template.get(
                    "name", "user_template"
                )
                mode_name = f"template_user_{template_name}"

                # Store template data for the scene to use
                self.scene.user_template_data = self._template_dialog.selected_template
                self.set_mode(mode_name)

                # Update status
                self.statusBar().showMessage(f"Template mode: {template_name}")

        self._template_dialog.finished.connect(on_dialog_finished)

    def save_2d_as_template(self):
        """Save current 2D structure as a template"""
        if not self.data.atoms:
            QMessageBox.warning(self, "Warning", "No structure to save as template.")
            return

        # Get template name
        name, ok = QInputDialog.getText(self, "Save Template", "Enter template name:")
        if not ok or not name.strip():
            return

        name = name.strip()

        try:
            # Template directory
            template_dir = os.path.join(self.settings_dir, "user-templates")
            if not os.path.exists(template_dir):
                os.makedirs(template_dir)

            # Convert current structure to template format
            atoms_data = []
            bonds_data = []

            # Convert atoms
            for atom_id, atom_info in self.data.atoms.items():
                pos = atom_info["pos"]
                atoms_data.append(
                    {
                        "id": atom_id,
                        "symbol": atom_info["symbol"],
                        "x": pos.x(),
                        "y": pos.y(),
                        "charge": atom_info.get("charge", 0),
                        "radical": atom_info.get("radical", 0),
                    }
                )

            # Convert bonds
            for (atom1_id, atom2_id), bond_info in self.data.bonds.items():
                bonds_data.append(
                    {
                        "atom1": atom1_id,
                        "atom2": atom2_id,
                        "order": bond_info["order"],
                        "stereo": bond_info.get("stereo", 0),
                    }
                )

            # Create template data
            template_data = {
                "format": "PME Template",
                "version": "1.0",
                "application": "MoleditPy",
                "application_version": VERSION,
                "name": name,
                "created": str(QDateTime.currentDateTime().toString()),
                "atoms": atoms_data,
                "bonds": bonds_data,
            }

            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(template_dir, filename)

            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self,
                    "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return

            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)

            # Mark as saved (no unsaved changes for this operation)
            self.has_unsaved_changes = False
            self.update_window_title()

            QMessageBox.information(
                self, "Success", f"Template '{name}' saved successfully."
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.critical(self, "Error", f"Failed to save template: {str(e)}")

    def open_translation_dialog(self):
        """Open the translation dialog"""
        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = TranslationDialog(self.current_mol, self, parent=self)
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage("Translation applied.")
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(
            lambda: self.remove_dialog_from_list(dialog)
        )  # Remove from list when dialog is closed

    def open_move_group_dialog(self):
        """Open Move Group dialog"""
        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = MoveGroupDialog(self.current_mol, self, parent=self)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage("Group transformation applied.")
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))

    def open_align_plane_dialog(self, plane):
        """Open align dialog"""
        # Get pre-selected atoms (before disabling measurement mode)
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = AlignPlaneDialog(
            self.current_mol, self, plane, preselected_atoms, parent=self
        )
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage(
                f"Atoms alignd to {plane.upper()} plane."
            )
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(
            lambda: self.remove_dialog_from_list(dialog)
        )  # Remove from list when dialog is closed

    def open_planarize_dialog(self, plane=None):
        """Open dialog to project selected atoms to the best-fit plane"""
        # Get pre-selected atoms
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = PlanarizeDialog(self.current_mol, self, preselected_atoms, parent=self)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage(
                "Selection planarized to best-fit plane."
            )
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))

    def open_alignment_dialog(self, axis):
        """Open alignment dialog"""
        # Get pre-selected atoms
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = AlignmentDialog(
            self.current_mol, self, axis, preselected_atoms, parent=self
        )
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage(
                f"Atoms aligned to {axis.upper()}-axis."
            )
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(
            lambda: self.remove_dialog_from_list(dialog)
        )  # Remove from list when dialog is closed

    def open_bond_length_dialog(self):
        """Open bond length adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = BondLengthDialog(
            self.current_mol, self, preselected_atoms, parent=self
        )
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage("Bond length adjusted.")
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(
            lambda: self.remove_dialog_from_list(dialog)
        )  # Remove from list when dialog is closed

    def open_angle_dialog(self):
        """Open angle adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = AngleDialog(self.current_mol, self, preselected_atoms, parent=self)
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Angle adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(
            lambda: self.remove_dialog_from_list(dialog)
        )  # Remove from list when dialog is closed

    def open_dihedral_dialog(self):
        """Open dihedral angle adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = []
        if hasattr(self, "selected_atoms_3d") and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif (
            hasattr(self, "selected_atoms_for_measurement")
            and self.selected_atoms_for_measurement
        ):
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = DihedralDialog(self.current_mol, self, preselected_atoms, parent=self)
        self.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.statusBar().showMessage("Dihedral angle adjusted.")
        )
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))

    def open_mirror_dialog(self):
        """Open mirror function dialog"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = MirrorDialog(self.current_mol, self)
        dialog.exec()  # Display as modal dialog

    def open_constrained_optimization_dialog(self):
        """Open constrained optimization dialog"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return

        # Disable measurement mode
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = ConstrainedOptimizationDialog(self.current_mol, self, parent=self)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))


MainWindowDialogManager._cls = MainWindowDialogManager
