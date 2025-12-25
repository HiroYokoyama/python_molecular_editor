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
    QDialog, QVBoxLayout, QLabel, QWidget, QGridLayout, QScrollArea,
    QHBoxLayout, QPushButton, QGraphicsScene, QInputDialog, QMessageBox
)
from PyQt6.QtGui import QPainter, QFont, QColor, QPen, QBrush
from PyQt6.QtCore import Qt, QPointF, QRectF, QTimer, QDateTime, QLineF
from .template_preview_view import TemplatePreviewView
try:
    from .constants import VERSION, CPK_COLORS
except Exception:
    from modules.constants import VERSION, CPK_COLORS
import os
import json
import logging

class UserTemplateDialog(QDialog):
    """ユーザーテンプレート管理ダイアログ"""
    
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.user_templates = []
        self.selected_template = None
        self.init_ui()
        self.load_user_templates()
    
    def init_ui(self):
        self.setWindowTitle("User Templates")
        self.setModal(False)  # モードレスに変更
        self.resize(800, 600)
        
        # ウィンドウを右上に配置
        if self.parent():
            parent_geometry = self.parent().geometry()
            x = parent_geometry.right() - self.width() - 20
            y = parent_geometry.top() + 50
            self.move(x, y)
        
        layout = QVBoxLayout(self)
        
        # Instructions
        instruction_label = QLabel("Create and manage your custom molecular templates. Click a template to use it in the editor.")
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)
        
        # Template grid
        self.template_widget = QWidget()
        self.template_layout = QGridLayout(self.template_widget)
        self.template_layout.setSpacing(10)
        
        scroll_area = QScrollArea()
        scroll_area.setWidget(self.template_widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumHeight(400)
        layout.addWidget(scroll_area)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.save_current_button = QPushButton("Save Current 2D as Template")
        self.save_current_button.clicked.connect(self.save_current_as_template)
        button_layout.addWidget(self.save_current_button)
        
        button_layout.addStretch()
        
        self.delete_button = QPushButton("Delete Selected")
        self.delete_button.clicked.connect(self.delete_selected_template)
        self.delete_button.setEnabled(False)
        button_layout.addWidget(self.delete_button)
        
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        button_layout.addWidget(close_button)
        
        layout.addLayout(button_layout)

    def closeEvent(self, event):
        """ダイアログクローズ時にモードをリセット"""
        self.cleanup_template_mode()
        super().closeEvent(event)

    def cleanup_template_mode(self):
        """テンプレートモードを終了し、atom_C(炭素描画)モードに戻す (Defensive implementation)"""
        # 1. Reset Dialog State
        self.selected_template = None
        if hasattr(self, 'delete_button'):
            self.delete_button.setEnabled(False)

        # 2. Reset Main Window Mode (UI/Toolbar)
        target_mode = 'atom_C'
        try:
            if hasattr(self.main_window, 'set_mode_and_update_toolbar'):
                 self.main_window.set_mode_and_update_toolbar(target_mode)
            elif hasattr(self.main_window, 'set_mode'):
                 self.main_window.set_mode(target_mode)
            
            # Fallback: set attribute directly if methods fail/don't exist
            if hasattr(self.main_window, 'mode'):
                self.main_window.mode = target_mode
        except Exception as e:
            logging.error(f"Error resetting main window mode: {e}")

        # 3. Reset Scene State (The Source of Truth)
        try:
            if hasattr(self.main_window, 'scene') and self.main_window.scene:
                 scene = self.main_window.scene
                 
                 # A. FORCE MODE
                 scene.mode = target_mode
                 scene.current_atom_symbol = 'C'
                 
                 # B. Clear Data
                 if hasattr(scene, 'user_template_data'):
                     scene.user_template_data = None
                 if hasattr(scene, 'template_context'):
                     scene.template_context = {}
                 
                 # C. Clear/Hide Preview Item
                 if hasattr(scene, 'clear_template_preview'):
                     scene.clear_template_preview()
                 
                 if hasattr(scene, 'template_preview') and scene.template_preview:
                     scene.template_preview.hide()
                     
                 # D. Reset Cursor & View
                 if scene.views():
                     view = scene.views()[0]
                     view.setCursor(Qt.CursorShape.CrossCursor)
                     view.viewport().update()
                 
                 scene.update()
        except Exception as e:
            logging.error(f"Error cleaning up scene state: {e}")
    
    def resizeEvent(self, event):
        """ダイアログリサイズ時にテンプレートプレビューを再フィット"""
        super().resizeEvent(event)
        # Delay the refit to ensure proper widget sizing
        QTimer.singleShot(100, self.refit_all_previews)
    
    def refit_all_previews(self):
        """すべてのテンプレートプレビューを再フィット"""
        try:
            for i in range(self.template_layout.count()):
                item = self.template_layout.itemAt(i)
                if item and item.widget():
                    widget = item.widget()
                    # Find the TemplatePreviewView within this widget
                    for child in widget.findChildren(TemplatePreviewView):
                        if hasattr(child, 'redraw_with_current_size'):
                            # Use redraw for better scaling adaptation
                            child.redraw_with_current_size()
                        elif hasattr(child, 'refit_view'):
                            child.refit_view()
        except Exception as e:
            logging.warning(f"Warning: Failed to refit template previews: {e}")
    
    def showEvent(self, event):
        """ダイアログ表示時にプレビューを適切にフィット"""
        super().showEvent(event)
        # Ensure all previews are properly fitted when dialog becomes visible
        QTimer.singleShot(300, self.refit_all_previews)
    
    def get_template_directory(self):
        """テンプレートディレクトリのパスを取得"""
        template_dir = os.path.join(self.main_window.settings_dir, 'user-templates')
        if not os.path.exists(template_dir):
            os.makedirs(template_dir)
        return template_dir
    
    def load_user_templates(self):
        """ユーザーテンプレートを読み込み"""
        template_dir = self.get_template_directory()
        self.user_templates.clear()
        
        try:
            for filename in os.listdir(template_dir):
                if filename.endswith('.pmetmplt'):
                    filepath = os.path.join(template_dir, filename)
                    template_data = self.load_template_file(filepath)
                    if template_data:
                        template_data['filename'] = filename
                        template_data['filepath'] = filepath
                        self.user_templates.append(template_data)
        except Exception as e:
            logging.error(f"Error loading user templates: {e}")
        
        self.update_template_grid()
    
    def load_template_file(self, filepath):
        """テンプレートファイルを読み込み"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            logging.error(f"Error loading template file {filepath}: {e}")
            return None
    
    def save_template_file(self, filepath, template_data):
        """テンプレートファイルを保存"""
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(template_data, f, indent=2, ensure_ascii=False)
            return True
        except Exception as e:
            logging.error(f"Error saving template file {filepath}: {e}")
            return False
    
    def update_template_grid(self):
        """テンプレートグリッドを更新"""
        # Clear existing widgets
        for i in reversed(range(self.template_layout.count())):
            self.template_layout.itemAt(i).widget().setParent(None)
        
        # Add template previews (left-to-right, top-to-bottom ordering)
        cols = 4
        for i, template in enumerate(self.user_templates):
            row = i // cols
            col = i % cols  # Left-to-right ordering
            
            preview_widget = self.create_template_preview(template)
            self.template_layout.addWidget(preview_widget, row, col)
        
        # Ensure all previews are properly fitted after grid update
        QTimer.singleShot(200, self.refit_all_previews)
    
    def create_template_preview(self, template_data):
        """テンプレートプレビューウィジェットを作成"""
        widget = QWidget()
        widget.setFixedSize(180, 200)
        widget.setStyleSheet("""
            QWidget {
                border: 2px solid #ccc;
                border-radius: 8px;
                background-color: white;
            }
            QWidget:hover {
                border-color: #007acc;
                background-color: #f0f8ff;
            }
        """)
        
        layout = QVBoxLayout(widget)
        
        # Preview graphics - use custom view class for better resize handling
        preview_scene = QGraphicsScene()
        preview_view = TemplatePreviewView(preview_scene)
        preview_view.setFixedSize(160, 140)
        preview_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        preview_view.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        preview_view.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        
        # Set template data for dynamic redrawing
        preview_view.set_template_data(template_data, self)
        
        # Draw template structure with view size for proper scaling
        view_size = (preview_view.width(), preview_view.height())
        self.draw_template_preview(preview_scene, template_data, view_size)
        
        # Improved fitting approach with better error handling
        bounding_rect = preview_scene.itemsBoundingRect()
        if not bounding_rect.isEmpty() and bounding_rect.width() > 0 and bounding_rect.height() > 0:
            # Calculate appropriate padding based on content size
            content_size = max(bounding_rect.width(), bounding_rect.height())
            padding = max(20, content_size * 0.2)  # At least 20 units or 20% of content
            
            padded_rect = bounding_rect.adjusted(-padding, -padding, padding, padding)
            preview_scene.setSceneRect(padded_rect)
            
            # Store original scene rect for proper fitting on resize
            preview_view.original_scene_rect = padded_rect
            
            # Use QTimer to ensure fitInView happens after widget is fully initialized
            QTimer.singleShot(0, lambda: self.fit_preview_view_safely(preview_view, padded_rect))
        else:
            # Default view for empty or invalid content
            default_rect = QRectF(-50, -50, 100, 100)
            preview_scene.setSceneRect(default_rect)
            preview_view.original_scene_rect = default_rect
            QTimer.singleShot(0, lambda: self.fit_preview_view_safely(preview_view, default_rect))
        
        layout.addWidget(preview_view)
        
        # Template name
        name = template_data.get('name', 'Unnamed Template')
        name_label = QLabel(name)
        name_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        name_label.setWordWrap(True)
        layout.addWidget(name_label)
        
        # Mouse events
        widget.mousePressEvent = lambda event: self.select_template(template_data, widget)
        widget.mouseDoubleClickEvent = lambda event: self.use_template(template_data)
        
        return widget
    
    def fit_preview_view_safely(self, view, rect):
        """プレビュービューを安全にフィット"""
        try:
            if view and not rect.isEmpty():
                view.fitInView(rect, Qt.AspectRatioMode.KeepAspectRatio)
        except Exception as e:
            logging.warning(f"Warning: Failed to fit preview view: {e}")
    
    def draw_template_preview(self, scene, template_data, view_size=None):
        """テンプレートプレビューを描画 - fitInView縮小率に基づく動的スケーリング"""
        atoms = template_data.get('atoms', [])
        bonds = template_data.get('bonds', [])
        
        if not atoms:
            # Add placeholder text when no atoms
            text = scene.addText("No structure", QFont('Arial', 12))
            text.setDefaultTextColor(QColor('gray'))
            return
        
        # Calculate molecular dimensions
        positions = [QPointF(atom['x'], atom['y']) for atom in atoms]
        min_x = min(pos.x() for pos in positions)
        max_x = max(pos.x() for pos in positions)
        min_y = min(pos.y() for pos in positions)
        max_y = max(pos.y() for pos in positions)
        
        mol_width = max_x - min_x
        mol_height = max_y - min_y
        mol_size = max(mol_width, mol_height)
        
        # Calculate fit scale factor (how much fitInView will shrink the content)
        if view_size is None:
            view_size = (160, 140)  # Default preview view size
        
        view_width, view_height = view_size
        
        if mol_size > 0 and mol_width > 0 and mol_height > 0:
            # Calculate the padding that will be added
            padding = max(20, mol_size * 0.2)
            padded_width = mol_width + 2 * padding
            padded_height = mol_height + 2 * padding
            
            # Calculate how much fitInView will scale down the content
            # fitInView fits the padded rectangle into the view while maintaining aspect ratio
            fit_scale_x = view_width / padded_width
            fit_scale_y = view_height / padded_height
            fit_scale = min(fit_scale_x, fit_scale_y)  # fitInView uses the smaller scale
            
            # Compensate for the fit scaling to maintain visual thickness
            # When fit_scale is small (content heavily shrunk), we need thicker lines/fonts
            if fit_scale > 0:
                scale_factor = max(0.4, min(4.0, 1.0 / fit_scale))
            else:
                scale_factor = 4.0
            
            # Debug info (can be removed in production)
            # logging.debug(f"Mol size: {mol_size:.1f}, Fit scale: {fit_scale:.3f}, Scale factor: {scale_factor:.2f}")
        else:
            scale_factor = 1.0
        
        # Base sizes that look good at 1:1 scale
        base_bond_width = 1.8
        base_font_size = 11
        base_ellipse_width = 18
        base_ellipse_height = 14
        base_double_bond_offset = 3.5
        base_triple_bond_offset = 2.5
        
        # Apply inverse fit scaling to maintain visual consistency
        bond_width = max(1.0, min(8.0, base_bond_width * scale_factor))
        font_size = max(8, min(24, int(base_font_size * scale_factor)))
        ellipse_width = max(10, min(40, base_ellipse_width * scale_factor))
        ellipse_height = max(8, min(30, base_ellipse_height * scale_factor))
        double_bond_offset = max(2.0, min(10.0, base_double_bond_offset * scale_factor))
        triple_bond_offset = max(1.5, min(8.0, base_triple_bond_offset * scale_factor))
        
        # Create atom ID to index mapping for bond drawing
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get('id', i)  # Use id if available, otherwise use index
            atom_id_to_index[atom_id] = i
        
        # Draw bonds first using original coordinates with dynamic sizing
        for bond in bonds:
            atom1_id, atom2_id = bond['atom1'], bond['atom2']
            
            # Get atom indices from IDs
            atom1_idx = atom_id_to_index.get(atom1_id)
            atom2_idx = atom_id_to_index.get(atom2_id)
            
            if atom1_idx is not None and atom2_idx is not None and atom1_idx < len(atoms) and atom2_idx < len(atoms):
                pos1 = QPointF(atoms[atom1_idx]['x'], atoms[atom1_idx]['y'])
                pos2 = QPointF(atoms[atom2_idx]['x'], atoms[atom2_idx]['y'])
                
                # Draw bonds with proper order - dynamic thickness
                bond_order = bond.get('order', 1)
                pen = QPen(QColor('black'), bond_width)
                
                if bond_order == 2:
                    # Double bond - draw two parallel lines
                    line = QLineF(pos1, pos2)
                    if line.length() > 0:
                        normal = line.normalVector()
                        normal.setLength(double_bond_offset)
                        
                        line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                        line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                        
                        scene.addLine(line1, pen)
                        scene.addLine(line2, pen)
                    else:
                        scene.addLine(line, pen)
                elif bond_order == 3:
                    # Triple bond - draw three parallel lines
                    line = QLineF(pos1, pos2)
                    if line.length() > 0:
                        normal = line.normalVector()
                        normal.setLength(triple_bond_offset)
                        
                        # Center line
                        scene.addLine(line, pen)
                        # Side lines
                        line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                        line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                        
                        scene.addLine(line1, pen)
                        scene.addLine(line2, pen)
                    else:
                        scene.addLine(line, pen)
                else:
                    # Single bond
                    scene.addLine(QLineF(pos1, pos2), pen)
        
        # Draw only non-carbon atom labels with dynamic sizing
        for i, atom in enumerate(atoms):
            try:
                pos = QPointF(atom['x'], atom['y'])
                symbol = atom.get('symbol', 'C')
                
                # Draw atoms - white ellipse background to hide bonds, then CPK colored text
                if symbol != 'C':
                    # All non-carbon atoms including hydrogen: white background ellipse + CPK colored text
                    color = CPK_COLORS.get(symbol, CPK_COLORS.get('DEFAULT', QColor('#FF1493')))
                    
                    # Add white background ellipse to hide bonds - dynamic size
                    pen = QPen(Qt.GlobalColor.white, 0)  # No border
                    brush = QBrush(Qt.GlobalColor.white)
                    ellipse_x = pos.x() - ellipse_width/2
                    ellipse_y = pos.y() - ellipse_height/2
                    ellipse = scene.addEllipse(ellipse_x, ellipse_y, ellipse_width, ellipse_height, pen, brush)
                    
                    # Add CPK colored text label on top - dynamic font size
                    font = QFont("Arial", font_size, QFont.Weight.Bold)
                    text = scene.addText(symbol, font)
                    text.setDefaultTextColor(color)  # CPK colored text
                    text_rect = text.boundingRect()
                    text.setPos(pos.x() - text_rect.width()/2, pos.y() - text_rect.height()/2)
                    
            except Exception:
                continue
    
    def select_template(self, template_data, widget):
        """テンプレートを選択してテンプレートモードに切り替え"""
        # Clear previous selection styling
        for i in range(self.template_layout.count()):
            item = self.template_layout.itemAt(i)
            if item and item.widget():
                item.widget().setStyleSheet("""
                    QWidget {
                        border: 2px solid #ccc;
                        border-radius: 8px;
                        background-color: white;
                    }
                    QWidget:hover {
                        border-color: #007acc;
                        background-color: #f0f8ff;
                    }
                """)
        
        # Highlight selected widget - only border, no background change
        widget.setStyleSheet("""
            QWidget {
                border: 3px solid #007acc;
                border-radius: 8px;
                background-color: white;
            }
        """)
        
        self.selected_template = template_data
        self.delete_button.setEnabled(True)

        # Automatically switch to template mode when template is selected
        template_name = template_data.get('name', 'user_template')
        mode_name = f"template_user_{template_name}"

        # Store template data for the scene to use
        try:
            self.main_window.scene.user_template_data = template_data
        except Exception:
            # Best-effort: ignore if scene or attribute missing
            pass

        # Force the main window into the template mode.
        # Clear or uncheck any existing mode actions if present to avoid staying in another mode.
        try:
            # Uncheck all mode actions first (if a dict of QAction exists)
            if hasattr(self.main_window, 'mode_actions') and isinstance(self.main_window.mode_actions, dict):
                for act in self.main_window.mode_actions.values():
                    try:
                        act.setChecked(False)
                    except Exception:
                        continue

            # If main_window has a set_mode method, call it. Otherwise, try to set a mode attribute.
            if hasattr(self.main_window, 'set_mode') and callable(self.main_window.set_mode):
                self.main_window.set_mode(mode_name)
            else:
                # Fallback: set an attribute and try to update UI
                setattr(self.main_window, 'mode', mode_name)

            # Update UI
            try:
                self.main_window.statusBar().showMessage(f"Template mode: {template_name}")
            except Exception:
                # ignore status bar failures
                pass

            # If there is a matching QAction in mode_actions, check it
            try:
                if hasattr(self.main_window, 'mode_actions') and f"template_user_{template_name}" in self.main_window.mode_actions:
                    self.main_window.mode_actions[f"template_user_{template_name}"].setChecked(True)
            except Exception:
                pass
        except Exception as e:
            logging.warning(f"Warning: Failed to switch main window to template mode: {e}")
    
    def use_template(self, template_data):
        """テンプレートを使用（エディタに適用）"""
        try:
            # Switch to template mode
            template_name = template_data.get('name', 'user_template')
            mode_name = f"template_user_{template_name}"
            
            # Store template data for the scene to use
            try:
                self.main_window.scene.user_template_data = template_data
            except Exception:
                pass

            # Force the main window into the template mode (same approach as select_template)
            try:
                if hasattr(self.main_window, 'mode_actions') and isinstance(self.main_window.mode_actions, dict):
                    for act in self.main_window.mode_actions.values():
                        try:
                            act.setChecked(False)
                        except Exception:
                            continue

                if hasattr(self.main_window, 'set_mode') and callable(self.main_window.set_mode):
                    self.main_window.set_mode(mode_name)
                else:
                    setattr(self.main_window, 'mode', mode_name)

                try:
                    self.main_window.statusBar().showMessage(f"Template mode: {template_name}")
                except Exception:
                    pass

                # Mark selected and keep dialog open
                self.selected_template = template_data
            except Exception as e:
                logging.warning(f"Warning: Failed to switch main window to template mode: {e}")

            # Don't close dialog - keep it open for easy template switching
            # self.accept()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply template: {str(e)}")
    
    def save_current_as_template(self):
        """現在の2D構造をテンプレートとして保存"""
        if not self.main_window.data.atoms:
            QMessageBox.warning(self, "Warning", "No structure to save as template.")
            return
        
        # Get template name
        name, ok = QInputDialog.getText(self, "Save Template", "Enter template name:")
        if not ok or not name.strip():
            return
        
        name = name.strip()
        
        try:
            # Convert current structure to template format
            template_data = self.convert_structure_to_template(name)
            
            # Save to file
            filename = f"{name.replace(' ', '_')}.pmetmplt"
            filepath = os.path.join(self.get_template_directory(), filename)
            
            if os.path.exists(filepath):
                reply = QMessageBox.question(
                    self, "Overwrite Template",
                    f"Template '{name}' already exists. Overwrite?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return
            
            if self.save_template_file(filepath, template_data):
                # Mark main window as saved
                self.main_window.has_unsaved_changes = False
                self.main_window.update_window_title()
                
                QMessageBox.information(self, "Success", f"Template '{name}' saved successfully.")
                self.load_user_templates()  # Refresh the display
            else:
                QMessageBox.critical(self, "Error", "Failed to save template.")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save template: {str(e)}")
    
    def convert_structure_to_template(self, name):
        """現在の構造をテンプレート形式に変換"""
        atoms_data = []
        bonds_data = []
        
        # Convert atoms
        for atom_id, atom_info in self.main_window.data.atoms.items():
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
        for (atom1_id, atom2_id), bond_info in self.main_window.data.bonds.items():
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

        return template_data
    
    def delete_selected_template(self):
        """選択されたテンプレートを削除"""
        if not self.selected_template:
            return
        
        name = self.selected_template.get('name', 'Unknown')
        reply = QMessageBox.question(
            self, "Delete Template",
            f"Are you sure you want to delete template '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            try:
                filepath = self.selected_template['filepath']
                os.remove(filepath)
                QMessageBox.information(self, "Success", f"Template '{name}' deleted successfully.")
                self.load_user_templates()  # Refresh the display
                self.selected_template = None
                self.delete_button.setEnabled(False)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to delete template: {str(e)}")
