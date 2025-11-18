#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_window_ui_manager.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowUiManager
"""


import vtk


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QGraphicsView, QDialog, QMessageBox
)



from PyQt6.QtCore import (
    Qt, QEvent, 
    QTimer
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
    from .custom_interactor_style import CustomInteractorStyle
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.custom_interactor_style import CustomInteractorStyle


# --- クラス定義 ---
class MainWindowUiManager(object):
    """ main_window.py から分離された機能クラス """

    def __init__(self, main_window):
        """ クラスの初期化 """
        self.mw = main_window


    def update_status_bar(self, message):
        """ワーカースレッドからのメッセージでステータスバーを更新するスロット"""
        self.statusBar().showMessage(message)



    def set_mode(self, mode_str):
        prev_mode = getattr(self.scene, 'mode', None)
        self.scene.mode = mode_str
        self.view_2d.setMouseTracking(True)
        # テンプレートモードから離れる場合はゴーストを消す
        if prev_mode and prev_mode.startswith('template') and not mode_str.startswith('template'):
            self.scene.clear_template_preview()
        elif not mode_str.startswith('template'):
            self.scene.template_preview.hide()

        # カーソル形状の設定
        if mode_str == 'select':
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)
        elif mode_str.startswith(('atom', 'bond', 'template')):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        elif mode_str.startswith(('charge', 'radical')):
            self.view_2d.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self.view_2d.setCursor(Qt.CursorShape.ArrowCursor)

        if mode_str.startswith('atom'): 
            self.scene.current_atom_symbol = mode_str.split('_')[1]
            self.statusBar().showMessage(f"Mode: Draw Atom ({self.scene.current_atom_symbol})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True) 
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0
        elif mode_str.startswith('bond'):
            self.scene.current_atom_symbol = 'C'
            parts = mode_str.split('_')
            self.scene.bond_order = int(parts[1])
            self.scene.bond_stereo = int(parts[2]) if len(parts) > 2 else 0
            stereo_text = {0: "", 1: " (Wedge)", 2: " (Dash)"}.get(self.scene.bond_stereo, "")
            self.statusBar().showMessage(f"Mode: Draw Bond (Order: {self.scene.bond_order}{stereo_text})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
            self.view_2d.setMouseTracking(True)
        elif mode_str.startswith('template'):
            if mode_str.startswith('template_user'):
                # User template mode
                template_name = mode_str.replace('template_user_', '')
                self.statusBar().showMessage(f"Mode: User Template ({template_name})")
            else:
                # Built-in template mode
                self.statusBar().showMessage(f"Mode: {mode_str.split('_')[1].capitalize()} Template")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'charge_plus':
            self.statusBar().showMessage("Mode: Increase Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'charge_minus':
            self.statusBar().showMessage("Mode: Decrease Charge (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif mode_str == 'radical':
            self.statusBar().showMessage("Mode: Toggle Radical (Click on Atom)")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)

        else: # Select mode
            self.statusBar().showMessage("Mode: Select")
            self.view_2d.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
            self.scene.bond_order = 1
            self.scene.bond_stereo = 0



    def set_mode_and_update_toolbar(self, mode_str):
        self.set_mode(mode_str)
        # QAction→QToolButtonのマッピングを取得
        toolbar = getattr(self, 'toolbar', None)
        action_to_button = {}
        if toolbar:
            for key, action in self.mode_actions.items():
                btn = toolbar.widgetForAction(action)
                if btn:
                    action_to_button[action] = btn

        # すべてのモードボタンの選択解除＆色リセット
        for key, action in self.mode_actions.items():
            action.setChecked(False)
            btn = action_to_button.get(action)
            if btn:
                btn.setStyleSheet("")

        # テンプレート系（User含む）は全て同じスタイル適用
        if mode_str in self.mode_actions:
            action = self.mode_actions[mode_str]
            action.setChecked(True)
            btn = action_to_button.get(action)
            if btn:
                # テンプレート系は青、それ以外はクリア
                if mode_str.startswith('template'):
                    btn.setStyleSheet("background-color: #2196F3; color: white;")
                else:
                    btn.setStyleSheet("")



    def activate_select_mode(self):
        self.set_mode('select')
        if 'select' in self.mode_actions:
            self.mode_actions['select'].setChecked(True)




    def eventFilter(self, obj, event):
        if obj is self.plotter and event.type() == QEvent.Type.MouseButtonPress:
            self.view_2d.setFocus()
        return super().eventFilter(obj, event)



    def closeEvent(self, event):
        # Persist settings on exit only when explicitly modified (deferred save)
        try:
            if getattr(self, 'settings_dirty', False) or self.settings != self.initial_settings:
                self.save_settings()
                self.settings_dirty = False
        except Exception:
            pass
        
        # 未保存の変更がある場合の処理
        if self.has_unsaved_changes:
            reply = QMessageBox.question(
                self, "Unsaved Changes",
                "You have unsaved changes. Do you want to save them?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Yes
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # 保存処理
                self.save_project()
                
                # 保存がキャンセルされた場合は終了もキャンセル
                if self.has_unsaved_changes:
                    event.ignore()
                    return
                    
            elif reply == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
            # No の場合はそのまま終了処理へ
        
        # 開いているすべてのダイアログウィンドウを閉じる
        try:
            for widget in QApplication.topLevelWidgets():
                if widget != self and isinstance(widget, (QDialog, QMainWindow)):
                    try:
                        widget.close()
                    except Exception:
                        pass
        except Exception:
            pass
        
        # 終了処理
        if self.scene and self.scene.template_preview:
            self.scene.template_preview.hide()

        # Clean up any active per-run calculation threads we spawned.
        try:
            for thr in list(getattr(self, '_active_calc_threads', []) or []):
                try:
                    thr.quit()
                except Exception:
                    pass
                try:
                    thr.wait(200)
                except Exception:
                    pass
        except Exception:
            pass
        
        event.accept()



    def toggle_3d_edit_mode(self, checked):
        """「3D Drag」ボタンの状態に応じて編集モードを切り替える"""
        if checked:
            # 3D Editモードをオンにする時は、Measurementモードを無効化
            if self.measurement_mode:
                self.measurement_action.setChecked(False)
                self.toggle_measurement_mode(False)
        
        self.is_3d_edit_mode = checked
        if checked:
            self.statusBar().showMessage("3D Drag Mode: ON.")
        else:
            self.statusBar().showMessage("3D Drag Mode: OFF.")
        self.view_2d.setFocus()



    def _setup_3d_picker(self):
        self.plotter.picker = vtk.vtkCellPicker()
        self.plotter.picker.SetTolerance(0.025)

        # 新しいカスタムスタイル（原子移動用）のインスタンスを作成
        style = CustomInteractorStyle(self)
        
        # 調査の結果、'style' プロパティへの代入が正しい設定方法と判明
        self.plotter.interactor.SetInteractorStyle(style)
        self.plotter.interactor.Initialize()



    def dragEnterEvent(self, event):
        """ウィンドウ全体で .pmeraw、.pmeprj、.mol、.sdf、.xyz ファイルのドラッグを受け入れる"""
        # Accept if any dragged local file has a supported extension
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            for url in urls:
                try:
                    if url.isLocalFile():
                        file_path = url.toLocalFile()
                        if file_path.lower().endswith(('.pmeraw', '.pmeprj', '.mol', '.sdf', '.xyz')):
                            event.acceptProposedAction()
                            return
                except Exception:
                    continue
        event.ignore()



    def dropEvent(self, event):
        """ファイルがウィンドウ上でドロップされたときに呼び出される"""
        urls = event.mimeData().urls()
        # Find the first local file from the dropped URLs
        file_path = None
        if urls:
            for url in urls:
                try:
                    if url.isLocalFile():
                        file_path = url.toLocalFile()
                        break
                except Exception:
                    continue

        if file_path:
            # ドロップ位置を取得
            drop_pos = event.position().toPoint()
            # 拡張子に応じて適切な読み込みメソッドを呼び出す
            if file_path.lower().endswith((".pmeraw", ".pmeprj")):
                self.open_project_file(file_path=file_path)
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            elif file_path.lower().endswith((".mol", ".sdf")):
                plotter_widget = self.splitter.widget(1)  # 3Dビューアーウィジェット
                plotter_rect = plotter_widget.geometry()
                if plotter_rect.contains(drop_pos):
                    self.load_mol_file_for_3d_viewing(file_path=file_path)
                else:
                    if hasattr(self, "load_mol_file"):
                        self.load_mol_file(file_path=file_path)
                    else:
                        self.statusBar().showMessage("MOL file import not implemented for 2D editor.")
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            elif file_path.lower().endswith(".xyz"):
                self.load_xyz_for_3d_viewing(file_path=file_path)
                QTimer.singleShot(100, self.fit_to_view)  # 遅延でFit
                event.acceptProposedAction()
            else:
                self.statusBar().showMessage(f"Unsupported file type: {file_path}")
                event.ignore()
        else:
            event.ignore()



    def _enable_3d_edit_actions(self, enabled=True):
        """3D編集機能のアクションを統一的に有効/無効化する"""
        actions = [
            'translation_action',
            'move_group_action',
            'alignplane_xy_action',
            'alignplane_xz_action',
            'alignplane_yz_action',
            'align_x_action',
            'align_y_action', 
            'align_z_action',
            'bond_length_action',
            'angle_action',
            'dihedral_action',
            'mirror_action',
            'planarize_action',
            'constrained_opt_action'
        ]
        
        # メニューとサブメニューも有効/無効化
        menus = [
            'align_menu'
        ]
        
        for action_name in actions:
            if hasattr(self, action_name):
                getattr(self, action_name).setEnabled(enabled)
        
        for menu_name in menus:
            if hasattr(self, menu_name):
                getattr(self, menu_name).setEnabled(enabled)



    def _enable_3d_features(self, enabled=True):
        """3D関連機能を統一的に有効/無効化する"""
        # 基本的な3D機能（3D SelectとEditは除外して常に有効にする）
        basic_3d_actions = [
            'optimize_3d_button',
            'export_button', 
            'analysis_action'
        ]
        
        for action_name in basic_3d_actions:
            if hasattr(self, action_name):
                # If enabling globally but chemical sanitization failed earlier, keep Optimize 3D disabled
                # Keep Optimize disabled when any of these conditions are true:
                # - we're globally disabling 3D features (enabled==False)
                # - the current molecule was created via the "skip chemistry checks" XYZ path
                # - a prior chemistry check was attempted and failed
                if action_name == 'optimize_3d_button':
                    try:
                        # If we're disabling all 3D features, ensure Optimize is disabled
                        if not enabled:
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # If the current molecule was marked as XYZ-derived (skip path), keep Optimize disabled
                        if getattr(self, 'is_xyz_derived', False):
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # If a chemistry check was tried and failed, keep Optimize disabled
                        if getattr(self, 'chem_check_tried', False) and getattr(self, 'chem_check_failed', False):
                            getattr(self, action_name).setEnabled(False)
                            continue

                        # Otherwise enable/disable according to the requested global flag
                        getattr(self, action_name).setEnabled(bool(enabled))
                    except Exception:
                        pass
                else:
                    try:
                        getattr(self, action_name).setEnabled(enabled)
                    except Exception:
                        pass
        
        # 3D Selectボタンは常に有効にする
        if hasattr(self, 'measurement_action'):
            self.measurement_action.setEnabled(True)
        
        # 3D Dragボタンも常に有効にする
        if hasattr(self, 'edit_3d_action'):
            self.edit_3d_action.setEnabled(True)
        
        # 3D編集機能も含める
        if enabled:
            self._enable_3d_edit_actions(True)
        else:
            self._enable_3d_edit_actions(False)



    def _enter_3d_viewer_ui_mode(self):
        """3DビューアモードのUI状態に設定する"""
        self.is_2d_editable = False
        self.cleanup_button.setEnabled(False)
        self.convert_button.setEnabled(False)
        for action in self.tool_group.actions():
            action.setEnabled(False)
        if hasattr(self, 'other_atom_action'):
            self.other_atom_action.setEnabled(False)
        
        self.minimize_2d_panel()

        # 3D関連機能を統一的に有効化
        self._enable_3d_features(True)



    def restore_ui_for_editing(self):
        """Enables all 2D editing UI elements."""
        self.is_2d_editable = True
        self.restore_2d_panel()
        self.cleanup_button.setEnabled(True)
        self.convert_button.setEnabled(True)

        for action in self.tool_group.actions():
            action.setEnabled(True)
        
        if hasattr(self, 'other_atom_action'):
            self.other_atom_action.setEnabled(True)
            
        # 2Dモードに戻る時は3D編集機能を統一的に無効化
        self._enable_3d_edit_actions(False)



    def minimize_2d_panel(self):
        """2Dパネルを最小化（非表示に）する"""
        sizes = self.splitter.sizes()
        # すでに最小化されていなければ実行
        if sizes[0] > 0:
            total_width = sum(sizes)
            self.splitter.setSizes([0, total_width])



    def restore_2d_panel(self):
        """最小化された2Dパネルを元のサイズに戻す"""
        sizes = self.splitter.sizes()
        
        # sizesリストが空でないことを確認してからアクセスする
        if sizes and sizes[0] == 0:
            self.splitter.setSizes([600, 600])



    def set_panel_layout(self, left_percent, right_percent):
        """パネルレイアウトを指定した比率に設定する"""
        if left_percent + right_percent != 100:
            return
        
        total_width = self.splitter.width()
        if total_width <= 0:
            total_width = 1200  # デフォルト幅
        
        left_width = int(total_width * left_percent / 100)
        right_width = int(total_width * right_percent / 100)
        
        self.splitter.setSizes([left_width, right_width])
        
        # ユーザーにフィードバック表示
        self.statusBar().showMessage(
            f"Panel layout set to {left_percent}% : {right_percent}%", 
            2000
        )



    def toggle_2d_panel(self):
        """2Dパネルの表示/非表示を切り替える"""
        sizes = self.splitter.sizes()
        if not sizes:
            return
            
        if sizes[0] == 0:
            # 2Dパネルが非表示の場合は表示
            self.restore_2d_panel()
            self.statusBar().showMessage("2D panel restored", 1500)
        else:
            # 2Dパネルが表示されている場合は非表示
            self.minimize_2d_panel()
            self.statusBar().showMessage("2D panel minimized", 1500)



    def on_splitter_moved(self, pos, index):
        """スプリッターが移動された時のフィードバック表示"""
        sizes = self.splitter.sizes()
        if len(sizes) >= 2:
            total = sum(sizes)
            if total > 0:
                left_percent = round(sizes[0] * 100 / total)
                right_percent = round(sizes[1] * 100 / total)
                
                # 現在の比率をツールチップで表示
                if hasattr(self.splitter, 'handle'):
                    handle = self.splitter.handle(1)
                    if handle:
                        handle.setToolTip(f"2D: {left_percent}% | 3D: {right_percent}%")



    def setup_splitter_tooltip(self):
        """スプリッターハンドルの初期ツールチップを設定"""
        handle = self.splitter.handle(1)
        if handle:
            handle.setToolTip("Drag to resize panels | Ctrl+1/2/3 for presets | Ctrl+H to toggle 2D panel")
            # 初期サイズ比率も表示
            self.on_splitter_moved(0, 0)

            
