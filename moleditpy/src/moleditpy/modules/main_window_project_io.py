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
main_window_project_io.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowProjectIo
"""


import pickle
import os
import json 
import traceback


# RDKit imports (explicit to satisfy flake8 and used features)
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QFileDialog, QMessageBox
)



from PyQt6.QtCore import (
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
    pass
except Exception:
    # Fallback to absolute imports for script-style execution
    pass


# --- クラス定義 ---
class MainWindowProjectIo(object):
    """ main_window.py から分離された機能クラス """


    def save_project(self):
        """上書き保存（Ctrl+S）- デフォルトでPMEPRJ形式"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
        # 非ネイティブ形式（.mol, .sdf, .xyz など）は上書き保存せず、必ず「名前を付けて保存」にする
        native_exts = ['.pmeprj', '.pmeraw']
        if self.current_file_path and any(self.current_file_path.lower().endswith(ext) for ext in native_exts):
            # 既存のPMEPRJ/PMERAWファイルの場合は上書き保存
            try:
                if self.current_file_path.lower().endswith('.pmeraw'):
                    # 既存のPMERAWファイルの場合はPMERAW形式で保存
                    save_data = self.get_current_state()
                    with open(self.current_file_path, 'wb') as f: 
                        pickle.dump(save_data, f)
                else:
                    # PMEPRJ形式で保存
                    json_data = self.create_json_data()
                    with open(self.current_file_path, 'w', encoding='utf-8') as f: 
                        json.dump(json_data, f, indent=2, ensure_ascii=False)
                
                # 保存成功時に状態をリセット
                self.has_unsaved_changes = False
                self.update_window_title()
                
                self.statusBar().showMessage(f"Project saved to {self.current_file_path}")
                
            except (OSError, IOError) as e:
                self.statusBar().showMessage(f"File I/O error: {e}")
            except (pickle.PicklingError, TypeError, ValueError) as e:
                self.statusBar().showMessage(f"Data serialization error: {e}")
            except Exception as e: 
                self.statusBar().showMessage(f"Error saving project file: {e}")
                
                traceback.print_exc()
        else:
            # MOL/SDF/XYZなどは上書き保存せず、必ず「名前を付けて保存」にする
            self.save_project_as()



    def save_project_as(self):
        """名前を付けて保存（Ctrl+Shift+S）- デフォルトでPMEPRJ形式"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            # Determine a sensible default filename based on current file (strip extension)
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # Prefer the directory of the currently opened file as default
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save Project As", default_path, 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeprj'): 
                file_path += '.pmeprj'
            
            # JSONデータを保存
            json_data = self.create_json_data()
            with open(file_path, 'w', encoding='utf-8') as f: 
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Replace current file with the newly saved file so subsequent saves go to this path
            self.current_file_path = file_path
            self.update_window_title()
            # Mark this state as the last saved state for undo tracking
            try:
                self._saved_state = copy.deepcopy(self.get_current_state())
            except Exception:
                pass
            
            self.statusBar().showMessage(f"Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving project file: {e}")
            
            traceback.print_exc()



    def save_raw_data(self):
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            save_data = self.get_current_state()
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(self, "Save Project File", default_path, "Project Files (*.pmeraw);;All Files (*)")
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeraw'): 
                file_path += '.pmeraw'
                
            with open(file_path, 'wb') as f: 
                pickle.dump(save_data, f)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Update current file to the newly saved raw file
            self.current_file_path = file_path
            self.update_window_title()
            try:
                self._saved_state = copy.deepcopy(self.get_current_state())
            except Exception:
                pass
            
            self.statusBar().showMessage(f"Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.PicklingError as e:
            self.statusBar().showMessage(f"Data serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving project file: {e}")
            
            traceback.print_exc()




    def load_raw_data(self, file_path=None):
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)")
            if not file_path: 
                return
        
        try:
            with open(file_path, 'rb') as f: 
                loaded_data = pickle.load(f)
            self.restore_ui_for_editing()
            self.set_state_from_data(loaded_data)
            
            # ファイル読み込み時に状態をリセット
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()
            try:
                self._saved_state= copy.deepcopy(self.et_current_state())
            except Exception:
                pass
         
            self.statusBar.showMessage(f"Project loaded from {file_path}")
            
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except pickle.UnpicklingError as e:
            self.statusBar().showMessage(f"Invalid project file format: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading project file: {e}")
            
            traceback.print_exc()



    def save_as_json(self):
        """PMEJSONファイル形式で保存 (3D MOL情報含む)"""
        if not self.data.atoms and not self.current_mol: 
            self.statusBar().showMessage("Error: Nothing to save.")
            return
            
        try:
            # default filename based on current file
            default_name = "untitled"
            try:
                if self.current_file_path:
                    base = os.path.basename(self.current_file_path)
                    default_name = os.path.splitext(base)[0]
            except Exception:
                default_name = "untitled"

            # prefer same directory as current file when available
            default_path = default_name
            try:
                if self.current_file_path:
                    default_path = os.path.join(os.path.dirname(self.current_file_path), default_name)
            except Exception:
                default_path = default_name

            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save as PME Project", default_path, 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path:
                return
                
            if not file_path.lower().endswith('.pmeprj'): 
                file_path += '.pmeprj'
            
            # JSONデータを作成
            json_data = self.create_json_data()
            
            # JSON形式で保存（美しい整形付き）
            with open(file_path, 'w', encoding='utf-8') as f: 
                json.dump(json_data, f, indent=2, ensure_ascii=False)
            
            # 保存成功時に状態をリセット
            self.has_unsaved_changes = False
            # Replace current file with the newly saved PME Project
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"PME Project saved to {file_path}")
            
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except (TypeError, ValueError) as e:
            self.statusBar().showMessage(f"JSON serialization error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error saving PME Project file: {e}")
            
            traceback.print_exc()



    def load_json_data(self, file_path=None):
        """PME Projectファイル形式を読み込み"""
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open PME Project File", "", 
                "PME Project Files (*.pmeprj);;All Files (*)", 
            )
            if not file_path: 
                return
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f: 
                json_data = json.load(f)
            
            # フォーマット検証
            if json_data.get("format") != "PME Project":
                QMessageBox.warning(
                    self, "Invalid Format", 
                    "This file is not a valid PME Project format."
                )
                return
            
            # バージョン確認
            file_version = json_data.get("version", "1.0")
            if file_version != "1.0":
                QMessageBox.information(
                    self, "Version Notice", 
                    f"This file was created with PME Project version {file_version}.\n"
                    "Loading will be attempted but some features may not work correctly."
                )
            
            self.restore_ui_for_editing()
            self.load_from_json_data(json_data)
            # ファイル読み込み時に状態をリセット
            self.reset_undo_stack()
            self.has_unsaved_changes = False
            self.current_file_path = file_path
            self.update_window_title()
            
            self.statusBar().showMessage(f"PME Project loaded from {file_path}")
            
            QTimer.singleShot(0, self.fit_to_view)
            
        except FileNotFoundError:
            self.statusBar().showMessage(f"File not found: {file_path}")
        except json.JSONDecodeError as e:
            self.statusBar().showMessage(f"Invalid JSON format: {e}")
        except (OSError, IOError) as e:
            self.statusBar().showMessage(f"File I/O error: {e}")
        except Exception as e: 
            self.statusBar().showMessage(f"Error loading PME Project file: {e}")
            
            traceback.print_exc()



    def open_project_file(self, file_path=None):
        """プロジェクトファイルを開く（.pmeprjと.pmerawの両方に対応）"""
        # Check for unsaved changes before opening a new project file.
        # Previously this function opened .pmeprj/.pmeraw without prompting the
        # user to save current unsaved work. Ensure we honor the global
        # unsaved-change check like other loaders (SMILES/MOL/etc.).
        if not self.check_unsaved_changes():
            return
        if not file_path:
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Open Project File", "", 
                "PME Project Files (*.pmeprj);;PME Raw Files (*.pmeraw);;All Files (*)", 
            )
            if not file_path: 
                return
        
        # 拡張子に応じて適切な読み込み関数を呼び出し
        if file_path.lower().endswith('.pmeprj'):
            self.load_json_data(file_path)
        elif file_path.lower().endswith('.pmeraw'):
            self.load_raw_data(file_path)
        else:
            # 拡張子不明の場合はJSONとして試行
            try:
                self.load_json_data(file_path)
            except Exception:
                try:
                    self.load_raw_data(file_path)
                except Exception:
                    self.statusBar().showMessage("Error: Unable to determine file format.")

