#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""
import logging  # [REPORT ERROR MISSING ATTRIBUTE]

import json
import os
from PyQt6.QtWidgets import QInputDialog, QMessageBox

try:
    from . import sip_isdeleted_safe
except ImportError:
    pass

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .about_dialog import AboutDialog
    from .align_plane_dialog import AlignPlaneDialog
    from .alignment_dialog import AlignmentDialog
    from .analysis_window import AnalysisWindow
    from .angle_dialog import AngleDialog
    from .bond_length_dialog import BondLengthDialog
    from .constrained_optimization_dialog import ConstrainedOptimizationDialog
    from .dihedral_dialog import DihedralDialog
    from .mirror_dialog import MirrorDialog
    from .move_group_dialog import MoveGroupDialog
    from .periodic_table_dialog import PeriodicTableDialog
    from .planarize_dialog import PlanarizeDialog
    from .settings_dialog import SettingsDialog
    from .color_settings_dialog import ColorSettingsDialog
    from .translation_dialog import TranslationDialog
    from .user_template_dialog import UserTemplateDialog
except ImportError:
    # Fallback to absolute imports for script-style execution
    from moleditpy.ui.about_dialog import AboutDialog
    from moleditpy.ui.align_plane_dialog import AlignPlaneDialog
    from moleditpy.ui.alignment_dialog import AlignmentDialog
    from moleditpy.ui.analysis_window import AnalysisWindow
    from moleditpy.ui.angle_dialog import AngleDialog
    from moleditpy.ui.bond_length_dialog import BondLengthDialog
    from moleditpy.ui.constrained_optimization_dialog import (
        ConstrainedOptimizationDialog,
    )
    from moleditpy.ui.dihedral_dialog import DihedralDialog
    from moleditpy.ui.mirror_dialog import MirrorDialog
    from moleditpy.ui.move_group_dialog import MoveGroupDialog
    from moleditpy.ui.periodic_table_dialog import PeriodicTableDialog
    from moleditpy.ui.planarize_dialog import PlanarizeDialog
    from moleditpy.ui.settings_dialog import SettingsDialog
    from moleditpy.ui.color_settings_dialog import ColorSettingsDialog
    from moleditpy.ui.translation_dialog import TranslationDialog
    from moleditpy.ui.user_template_dialog import UserTemplateDialog

# Import VERSION from constants
try:
    from .constants import VERSION
except ImportError:
    from moleditpy.utils.constants import VERSION

class DialogManager:
    """Independent manager for UI dialogs, ported from MainWindowDialogManager mixin."""

    def __init__(self, host):
        self.host = host

    def _get_preselected_atoms_3d(self):
        """Helper to collect preselected atoms from measurement mode (3D Select)."""
        preselected_atoms = []
        if hasattr(self.host, "edit_3d_manager"):
            if self.host.edit_3d_manager.selected_atoms_for_measurement:
                preselected_atoms = list(self.host.edit_3d_manager.selected_atoms_for_measurement)
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(f"REPORT ERROR: Missing attribute 'edit_3d_manager' on self.host")
        return preselected_atoms

    def show_about_dialog(self):
        """Show the custom About dialog with Easter egg functionality"""
        dialog = AboutDialog(self.host, self.host)
        dialog.exec()

    def open_periodic_table_dialog(self):
        dialog = PeriodicTableDialog(self.host)
        dialog.element_selected.connect(self.host.ui_manager.set_atom_from_periodic_table)
        checked_action = self.host.init_manager.tool_group.checkedAction()
        if checked_action:
            self.host.init_manager.tool_group.setExclusive(False)
            checked_action.setChecked(False)
            self.host.init_manager.tool_group.setExclusive(True)
        dialog.exec()

    def open_analysis_window(self):
        if self.host.view_3d_manager.current_mol:
            dialog = AnalysisWindow(
                self.host.view_3d_manager.current_mol,
                self.host,
                is_xyz_derived=self.host.is_xyz_derived,
            )
            dialog.exec()
        else:
            self.host.statusBar().showMessage(
                "Please generate a 3D structure first to show analysis."
            )

    def open_template_dialog(self):
        """Open the template dialog"""
        dialog = UserTemplateDialog(self.host, self.host)
        dialog.exec()

    def open_template_dialog_and_activate(self):
        """Open the template dialog and activate the selected template for use in the main window"""
        # Check for existing dialog
        _template_dialog = getattr(self.host, "_template_dialog", None)
        if _template_dialog and not _template_dialog.isHidden():
            # Bring existing dialog to front
            _template_dialog.raise_()
            _template_dialog.activateWindow()
            return

        # Create new dialog
        self.host._template_dialog = UserTemplateDialog(self.host, self.host)
        self.host._template_dialog.show()

        # Activate if a template is selected after dialog is closed
        def on_dialog_finished():
            if (
                hasattr(self.host._template_dialog, "selected_template")
                and self.host._template_dialog.selected_template
            ):
                template_name = self.host._template_dialog.selected_template.get(
                    "name", "user_template"
                )
                mode_name = f"template_user_{template_name}"

                # Store template data for the scene to use
                self.host.init_manager.scene.user_template_data = (
                    self.host._template_dialog.selected_template
                )
                self.host.ui_manager.set_mode(mode_name)

                # Update status
                self.host.statusBar().showMessage(f"Template mode: {template_name}")

        self.host._template_dialog.finished.connect(on_dialog_finished)

    def save_2d_as_template(self):
        """Save current 2D structure as a template"""
        if not self.host.state_manager.data.atoms:
            QMessageBox.warning(
                self.host, "Warning", "No structure to save as template."
            )
            return

        # Get template name
        name, ok = QInputDialog.getText(
            self.host, "Save Template", "Enter template name:"
        )
        if not ok or not name.strip():
            return

        name = name.strip()

        try:
            # Template directory
            template_dir = os.path.join(self.host.settings_dir, "user-templates")
            if not os.path.exists(template_dir):
                os.makedirs(template_dir)
            # Convert current structure to template format using core method
            template_data = self.host.state_manager.data.to_template_dict(
                name, application_version=VERSION
            )

            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(template_dir, filename)

            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self.host,
                    "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return

            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)

            # Mark as saved (no unsaved changes for this operation)
            self.host.state_manager.has_unsaved_changes = False
            self.host.state_manager.update_window_title()

            QMessageBox.information(
                self.host, "Success", f"Template '{name}' saved successfully."
            )

        except (AttributeError, RuntimeError, ValueError) as e:
            QMessageBox.critical(
                self.host, "Error", f"Failed to save template: {str(e)}"
            )

    def open_translation_dialog(self):
        """Open the translation dialog"""
        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        # Get preselected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        dialog = TranslationDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)  # Keep reference
        dialog.show()  # Use show for modeless display
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage("Translation applied.")
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_move_group_dialog(self):
        """Open Move Group dialog"""
        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        # Get preselected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        dialog = MoveGroupDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage("Group transformation applied.")
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_align_plane_dialog(self, plane):
        """Open align dialog"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = AlignPlaneDialog(
            self.host.view_3d_manager.current_mol, self.host, plane, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage(
                f"Atoms aligned to {plane.upper()} plane."
            )
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_planarize_dialog(self, plane=None):
        """Open dialog to project selected atoms to the best-fit plane"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = PlanarizeDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage(
                "Selection planarized to best-fit plane."
            )
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_alignment_dialog(self, axis):
        """Open alignment dialog"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = AlignmentDialog(
            self.host.view_3d_manager.current_mol, self.host, axis, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage(
                f"Atoms aligned to {axis.upper()}-axis."
            )
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_bond_length_dialog(self):
        """Open bond length adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = BondLengthDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage("Bond length adjusted.")
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_angle_dialog(self):
        """Open angle adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = AngleDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage("Angle adjusted.")
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_dihedral_dialog(self):
        """Open dihedral angle adjustment dialog"""
        # Get pre-selected atoms
        preselected_atoms = self._get_preselected_atoms_3d()

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = DihedralDialog(
            self.host.view_3d_manager.current_mol, self.host, preselected_atoms, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(
            lambda: self.host.statusBar().showMessage("Dihedral angle adjusted.")
        )
        dialog.accepted.connect(self.host.edit_actions_manager.push_undo_state)
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))

    def open_mirror_dialog(self):
        """Open mirror function dialog"""
        if not self.host.view_3d_manager.current_mol:
            self.host.statusBar().showMessage("No 3D molecule loaded.")
            return

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = MirrorDialog(self.host.view_3d_manager.current_mol, self.host)
        dialog.exec()

    def open_settings_dialog(self):
        """Open the application settings dialog."""
        dialog = SettingsDialog(self.host.init_manager.settings, parent=self.host)
        dialog.exec()

    def open_color_settings_dialog(self):
        """Open the CPK color settings dialog."""
        dialog = ColorSettingsDialog(self.host.init_manager.settings, parent=self.host)
        dialog.exec()

    def open_constrained_optimization_dialog(self):
        """Open constrained optimization dialog"""
        if not self.host.view_3d_manager.current_mol:
            self.host.statusBar().showMessage("No 3D molecule loaded.")
            return

        # Disable measurement mode
        if self.host.edit_3d_manager.measurement_mode:
            self.host.init_manager.measurement_action.setChecked(False)
            self.host.edit_3d_manager.toggle_measurement_mode(False)

        dialog = ConstrainedOptimizationDialog(
            self.host.view_3d_manager.current_mol, self.host, parent=self.host
        )
        self.host.edit_3d_manager.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.finished.connect(lambda: self.host.edit_3d_manager.remove_dialog_from_list(dialog))
