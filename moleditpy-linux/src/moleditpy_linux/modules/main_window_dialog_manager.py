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
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowDialogManager
"""


import os
import json 


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QMessageBox, 
    QInputDialog
)



from PyQt6.QtCore import (
    QDateTime
)


# Use centralized Open Babel availability from package-level __init__
# Use per-package modules availability (local __init__).
try:
    from . import OBABEL_AVAILABLE
except Exception:
    from modules import OBABEL_AVAILABLE
# Only import pybel on demand — `moleditpy` itself doesn't expose `pybel`.
if OBABEL_AVAILABLE:
    try:
        from openbabel import pybel
    except Exception:
        # If import fails here, disable OBABEL locally; avoid raising
        pybel = None
        OBABEL_AVAILABLE = False
        print("Warning: openbabel.pybel not available. Open Babel fallback and OBabel-based options will be disabled.")
else:
    pybel = None
    
# Optional SIP helper: on some PyQt6 builds sip.isdeleted is available and
# allows safely detecting C++ wrapper objects that have been deleted. Import
# it once at module import time and expose a small, robust wrapper so callers
# can avoid re-importing sip repeatedly and so we centralize exception
# handling (this reduces crash risk during teardown and deletion operations).
try:
    import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, 'isdeleted', None)
except Exception:
    _sip = None
    _sip_isdeleted = None

try:
    # package relative imports (preferred when running as `python -m moleditpy`)
    from .constants import VERSION
    from .user_template_dialog import UserTemplateDialog
    from .about_dialog import AboutDialog
    from .translation_dialog import TranslationDialog
    from .mirror_dialog import MirrorDialog
    from .move_group_dialog import MoveGroupDialog
    from .align_plane_dialog import AlignPlaneDialog
    from .planarize_dialog import PlanarizeDialog
    from .alignment_dialog import AlignmentDialog
    from .periodic_table_dialog import PeriodicTableDialog
    from .analysis_window import AnalysisWindow
    from .bond_length_dialog import BondLengthDialog
    from .angle_dialog import AngleDialog
    from .dihedral_dialog import DihedralDialog
    from .constrained_optimization_dialog import ConstrainedOptimizationDialog
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.constants import VERSION
    from modules.user_template_dialog import UserTemplateDialog
    from modules.about_dialog import AboutDialog
    from modules.translation_dialog import TranslationDialog
    from modules.mirror_dialog import MirrorDialog
    from modules.move_group_dialog import MoveGroupDialog
    from modules.align_plane_dialog import AlignPlaneDialog
    from modules.planarize_dialog import PlanarizeDialog
    from modules.alignment_dialog import AlignmentDialog
    from modules.periodic_table_dialog import PeriodicTableDialog
    from modules.analysis_window import AnalysisWindow
    from modules.bond_length_dialog import BondLengthDialog
    from modules.angle_dialog import AngleDialog
    from modules.dihedral_dialog import DihedralDialog
    from modules.constrained_optimization_dialog import ConstrainedOptimizationDialog


# --- クラス定義 ---
class MainWindowDialogManager(object):
    """ main_window.py から分離された機能クラス """


    def show_about_dialog(self):
        """Show the custom About dialog with Easter egg functionality"""
        dialog = AboutDialog(self, self)
        dialog.exec()



    def open_periodic_table_dialog(self):
        dialog=PeriodicTableDialog(self); dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action=self.tool_group.checkedAction()
        if checked_action: self.tool_group.setExclusive(False); checked_action.setChecked(False); self.tool_group.setExclusive(True)
        dialog.exec()



    def open_analysis_window(self):
        if self.current_mol:
            dialog = AnalysisWindow(self.current_mol, self, is_xyz_derived=self.is_xyz_derived)
            dialog.exec()
        else:
            self.statusBar().showMessage("Please generate a 3D structure first to show analysis.")



    def open_template_dialog(self):
        """テンプレートダイアログを開く"""
        dialog = UserTemplateDialog(self, self)
        dialog.exec()
    


    def open_template_dialog_and_activate(self):
        """テンプレートダイアログを開き、テンプレートがメイン画面で使用できるようにする"""
        # 既存のダイアログがあるかチェック
        if hasattr(self, '_template_dialog') and self._template_dialog and not self._template_dialog.isHidden():
            # 既存のダイアログを前面に表示
            self._template_dialog.raise_()
            self._template_dialog.activateWindow()
            return
        
        # 新しいダイアログを作成
        self._template_dialog = UserTemplateDialog(self, self)
        self._template_dialog.show()  # モードレスで表示
        
        # ダイアログが閉じられた後、テンプレートが選択されていればアクティブ化
        def on_dialog_finished():
            if hasattr(self._template_dialog, 'selected_template') and self._template_dialog.selected_template:
                template_name = self._template_dialog.selected_template.get('name', 'user_template')
                mode_name = f"template_user_{template_name}"
                
                # Store template data for the scene to use
                self.scene.user_template_data = self._template_dialog.selected_template
                self.set_mode(mode_name)
                
                # Update status
                self.statusBar().showMessage(f"Template mode: {template_name}")
        
        self._template_dialog.finished.connect(on_dialog_finished)
    


    def save_2d_as_template(self):
        """現在の2D構造をテンプレートとして保存"""
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
            template_dir = os.path.join(self.settings_dir, 'user-templates')
            if not os.path.exists(template_dir):
                os.makedirs(template_dir)
            
            # Convert current structure to template format
            atoms_data = []
            bonds_data = []
            
            # Convert atoms
            for atom_id, atom_info in self.data.atoms.items():
                pos = atom_info['pos']
                atoms_data.append({
                    'id': atom_id,
                    'symbol': atom_info['symbol'],
                    'x': pos.x(),
                    'y': pos.y(),
                    'charge': atom_info.get('charge', 0),
                    'radical': atom_info.get('radical', 0)
                })
            
            # Convert bonds
            for (atom1_id, atom2_id), bond_info in self.data.bonds.items():
                bonds_data.append({
                    'atom1': atom1_id,
                    'atom2': atom2_id,
                    'order': bond_info['order'],
                    'stereo': bond_info.get('stereo', 0)
                })
            
            # Create template data
            template_data = {
                'format': "PME Template",
                'version': "1.0",
                'application': "MoleditPy",
                'application_version': VERSION,
                'name': name,
                'created': str(QDateTime.currentDateTime().toString()),
                'atoms': atoms_data,
                'bonds': bonds_data
            }
            
            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(template_dir, filename)
            
            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self, "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return
            
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)
            
            # Mark as saved (no unsaved changes for this operation)
            self.has_unsaved_changes = False
            self.update_window_title()
            
            QMessageBox.information(self, "Success", f"Template '{name}' saved successfully.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save template: {str(e)}")



    def open_translation_dialog(self):
        """平行移動ダイアログを開く"""
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = TranslationDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Translation applied."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    


    def open_move_group_dialog(self):
        """Move Groupダイアログを開く"""
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = MoveGroupDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Group transformation applied."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    


    def open_align_plane_dialog(self, plane):
        """alignダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AlignPlaneDialog(self.current_mol, self, plane, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage(f"Atoms alignd to {plane.upper()} plane."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
        


    def open_planarize_dialog(self, plane=None):
        """選択原子群を最適平面へ投影するダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)

        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)

        dialog = PlanarizeDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)
        dialog.show()
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Selection planarized to best-fit plane."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    


    def open_alignment_dialog(self, axis):
        """アライメントダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AlignmentDialog(self.current_mol, self, axis, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage(f"Atoms aligned to {axis.upper()}-axis."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    


    def open_bond_length_dialog(self):
        """結合長変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = BondLengthDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Bond length adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    


    def open_angle_dialog(self):
        """角度変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = AngleDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Angle adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    


    def open_dihedral_dialog(self):
        """二面角変換ダイアログを開く"""
        # 事前選択された原子を取得（測定モード無効化前に）
        preselected_atoms = []
        if hasattr(self, 'selected_atoms_3d') and self.selected_atoms_3d:
            preselected_atoms = list(self.selected_atoms_3d)
        elif hasattr(self, 'selected_atoms_for_measurement') and self.selected_atoms_for_measurement:
            preselected_atoms = list(self.selected_atoms_for_measurement)
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = DihedralDialog(self.current_mol, self, preselected_atoms)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # execではなくshowを使用してモードレス表示
        dialog.accepted.connect(lambda: self.statusBar().showMessage("Dihedral angle adjusted."))
        dialog.accepted.connect(self.push_undo_state)
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))  # ダイアログが閉じられた時にリストから削除
    


    def open_mirror_dialog(self):
        """ミラー機能ダイアログを開く"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = MirrorDialog(self.current_mol, self)
        dialog.exec()  # モーダルダイアログとして表示



    def open_constrained_optimization_dialog(self):
        """制約付き最適化ダイアログを開く"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule loaded.")
            return
        
        # 測定モードを無効化
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)
        
        dialog = ConstrainedOptimizationDialog(self.current_mol, self)
        self.active_3d_dialogs.append(dialog)  # 参照を保持
        dialog.show()  # モードレス表示
        dialog.finished.connect(lambda: self.remove_dialog_from_list(dialog))
    
