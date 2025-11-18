#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main_window_compute.py
MainWindow (main_window.py) から分離されたモジュール
機能クラス: MainWindowCompute
"""




# RDKit imports (explicit to satisfy flake8 and used features)
from rdkit import Chem
from rdkit.Chem import AllChem
try:
    pass
except Exception:
    pass

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMenu
)

from PyQt6.QtGui import (
    QColor, QAction
)


from PyQt6.QtCore import (
    QThread, QTimer
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
    from .calculation_worker import CalculationWorker
except Exception:
    # Fallback to absolute imports for script-style execution
    from modules.calculation_worker import CalculationWorker


# --- クラス定義 ---
class MainWindowCompute(object):
    """ main_window.py から分離された機能クラス """

    def __init__(self, main_window):
        """ クラスの初期化 """
        self.mw = main_window


    def set_optimization_method(self, method_name):
        """Set preferred 3D optimization method and persist to settings.

        Supported values: 'GAFF', 'MMFF'
        """
        # Normalize input and validate
        if not method_name:
            return
        method = str(method_name).strip().upper()
        valid_methods = (
            'MMFF_RDKIT', 'MMFF94_RDKIT', 'UFF_RDKIT',
            'UFF_OBABEL', 'GAFF_OBABEL', 'MMFF94_OBABEL', 'GHEMICAL_OBABEL'
        )
        if method not in valid_methods:
            # Unknown method: ignore but notify
            self.statusBar().showMessage(f"Unknown 3D optimization method: {method_name}")
            return

        # Update internal state (store canonical uppercase key)
        self.optimization_method = method

        # Persist to settings
        try:
                self.settings['optimization_method'] = self.optimization_method
                try:
                    self.settings_dirty = True
                except Exception:
                    pass
        except Exception:
            pass

        # Update menu checked state if actions mapping exists
        try:
            if hasattr(self, 'opt3d_actions') and self.opt3d_actions:
                for k, act in self.opt3d_actions.items():
                    try:
                        # keys in opt3d_actions may be mixed-case; compare uppercased
                        act.setChecked(k.upper() == method)
                    except Exception:
                        pass
        except Exception:
            pass

        # Also show user-friendly label if available
        try:
            label = self.opt3d_method_labels.get(self.optimization_method, self.optimization_method)
        except Exception:
            label = self.optimization_method
        self.statusBar().showMessage(f"3D optimization method set to: {label}")



    def show_convert_menu(self, pos):
        """右クリックで表示する一時的な3D変換メニュー。
        選択したモードは一時フラグとして保持され、その後の変換で使用されます（永続化しません）。
        """
        try:
            menu = QMenu(self)
            conv_options = [
                ("RDKit -> Open Babel (fallback)", 'fallback'),
                ("RDKit only", 'rdkit'),
                ("Open Babel only", 'obabel'),
                ("Direct (use 2D coords + add H)", 'direct')
            ]
            for label, key in conv_options:
                a = QAction(label, self)
                # If Open Babel is not available, disable actions that depend on it
                if key in ('obabel', 'fallback') and not globals().get('OBABEL_AVAILABLE', False):
                    a.setEnabled(False)
                a.triggered.connect(lambda checked=False, k=key: self._trigger_conversion_with_temp_mode(k))
                menu.addAction(a)

            # Show menu at button position
            menu.exec_(self.convert_button.mapToGlobal(pos))
        except Exception as e:
            print(f"Error showing convert menu: {e}")




    def _trigger_conversion_with_temp_mode(self, mode_key):
        try:
            # store temporary override and invoke conversion
            self._temp_conv_mode = mode_key
            # Call the normal conversion entry point (it will consume the temp)
            QTimer.singleShot(0, self.trigger_conversion)
        except Exception as e:
            print(f"Failed to start conversion with temp mode {mode_key}: {e}")




    def show_optimize_menu(self, pos):
        """右クリックで表示する一時的な3D最適化メニュー。
        選択したメソッドは一時フラグとして保持され、その後の最適化で使用されます（永続化しません）。
        """
        try:
            menu = QMenu(self)
            opt_list = [
                ("MMFF94s", 'MMFF_RDKIT'),
                ("MMFF94", 'MMFF94_RDKIT'),
                ("UFF", 'UFF_RDKIT')
            ]
            for label, key in opt_list:
                a = QAction(label, self)
                # If opt3d_actions exist, reflect their enabled state
                try:
                    if hasattr(self, 'opt3d_actions') and key in self.opt3d_actions:
                        a.setEnabled(self.opt3d_actions[key].isEnabled())
                except Exception:
                    pass
                a.triggered.connect(lambda checked=False, k=key: self._trigger_optimize_with_temp_method(k))
                menu.addAction(a)

            menu.exec_(self.optimize_3d_button.mapToGlobal(pos))
        except Exception as e:
            print(f"Error showing optimize menu: {e}")




    def _trigger_optimize_with_temp_method(self, method_key):
        try:
            # store temporary override and invoke optimization
            self._temp_optimization_method = method_key
            # Run optimize on next event loop turn so UI updates first
            QTimer.singleShot(0, self.optimize_3d_structure)
        except Exception as e:
            print(f"Failed to start optimization with temp method {method_key}: {e}")



    def trigger_conversion(self):
        # Reset last successful optimization method at start of new conversion
        self.last_successful_optimization_method = None
        
        # 3D変換時に既存の3D制約をクリア
        self.constraints_3d = []

        # 2Dエディタに原子が存在しない場合は3Dビューをクリア
        if not self.data.atoms:
            self.plotter.clear()
            self.current_mol = None
            self.analysis_action.setEnabled(False)
            self.statusBar().showMessage("3D view cleared.")
            self.view_2d.setFocus() 
            return

        # 描画モード変更時に測定モードと3D編集モードをリセット
        if self.measurement_mode:
            self.measurement_action.setChecked(False)
            self.toggle_measurement_mode(False)  # 測定モードを無効化
        if self.is_3d_edit_mode:
            self.edit_3d_action.setChecked(False)
            self.toggle_3d_edit_mode(False)  # 3D編集モードを無効化

        mol = self.data.to_rdkit_mol(use_2d_stereo=False)

        # 分子オブジェクトが作成できない場合でも化学的問題をチェック
        if not mol or mol.GetNumAtoms() == 0:
            # RDKitでの変換に失敗した場合は、独自の化学的問題チェックを実行
            self.check_chemistry_problems_fallback()
            return

        # 原子プロパティを保存（ワーカープロセスで失われるため）
        self.original_atom_properties = {}
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            try:
                original_id = atom.GetIntProp("_original_atom_id")
                self.original_atom_properties[i] = original_id
            except KeyError:
                pass

        problems = Chem.DetectChemistryProblems(mol)
        if problems:
            # 化学的問題が見つかった場合は既存のフラグをクリアしてから新しい問題を表示
            self.scene.clear_all_problem_flags()
            self.statusBar().showMessage(f"Error: {len(problems)} chemistry problem(s) found.")
            # 既存の選択状態をクリア
            self.scene.clearSelection() 

            # 問題のある原子に赤枠フラグを立てる
            for prob in problems:
                atom_idx = prob.GetAtomIdx()
                rdkit_atom = mol.GetAtomWithIdx(atom_idx)
                # エディタ側での原子IDの取得と存在確認
                if rdkit_atom.HasProp("_original_atom_id"):
                    original_id = rdkit_atom.GetIntProp("_original_atom_id")
                    if original_id in self.data.atoms and self.data.atoms[original_id]['item']:
                        item = self.data.atoms[original_id]['item']
                        item.has_problem = True 
                        item.update()

            self.view_2d.setFocus()
            return

        # 化学的問題がない場合のみフラグをクリアして3D変換を実行
        self.scene.clear_all_problem_flags()

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus() 
            return

        # 複数分子の処理に対応
        num_frags = len(Chem.GetMolFrags(mol))
        if num_frags > 1:
            self.statusBar().showMessage(f"Converting {num_frags} molecules to 3D with collision detection...")
        else:
            self.statusBar().showMessage("Calculating 3D structure...")
            
        # CRITICAL FIX: Use the 2D editor's MOL block instead of RDKit's to preserve
        # wedge/dash stereo information that is stored in the 2D editor data.
        # RDKit's MolToMolBlock() doesn't preserve this information.
        mol_block = self.data.to_mol_block()
        if not mol_block:
            mol_block = Chem.MolToMolBlock(mol, includeStereo=True)
        
        # Additional E/Z stereo enhancement: add M CFG lines for explicit E/Z bonds
        mol_lines = mol_block.split('\n')
        
        # Find bonds with explicit E/Z labels from our data and map to RDKit bond indices
        ez_bond_info = {}
        for (id1, id2), bond_data in self.data.bonds.items():
            if bond_data.get('stereo') in [3, 4]:  # E/Z labels
                # Find corresponding atoms in RDKit molecule by _original_atom_id property
                rdkit_idx1 = None
                rdkit_idx2 = None
                for atom in mol.GetAtoms():
                    if atom.HasProp("_original_atom_id"):
                        orig_id = atom.GetIntProp("_original_atom_id")
                        if orig_id == id1:
                            rdkit_idx1 = atom.GetIdx()
                        elif orig_id == id2:
                            rdkit_idx2 = atom.GetIdx()
                
                if rdkit_idx1 is not None and rdkit_idx2 is not None:
                    rdkit_bond = mol.GetBondBetweenAtoms(rdkit_idx1, rdkit_idx2)
                    if rdkit_bond and rdkit_bond.GetBondType() == Chem.BondType.DOUBLE:
                        ez_bond_info[rdkit_bond.GetIdx()] = bond_data['stereo']
        
        # Add M  CFG lines for E/Z stereo if needed
        if ez_bond_info:
            insert_idx = len(mol_lines) - 1  # Before M  END
            for bond_idx, stereo_type in ez_bond_info.items():
                cfg_value = 1 if stereo_type == 3 else 2  # 1=Z, 2=E in MOL format
                cfg_line = f"M  CFG  1 {bond_idx + 1:3d}   {cfg_value}"
                mol_lines.insert(insert_idx, cfg_line)
                insert_idx += 1
            mol_block = '\n'.join(mol_lines)
        
        # Assign a unique ID for this conversion run so it can be halted/validated
        try:
            run_id = int(self.next_conversion_id)
        except Exception:
            run_id = 1
        try:
            self.next_conversion_id = run_id + 1
        except Exception:
            self.next_conversion_id = getattr(self, 'next_conversion_id', 1) + 1

        # Record this run as active. Use a set to track all active worker ids
        # so a Halt request can target every running conversion.
        try:
            self.active_worker_ids.add(run_id)
        except Exception:
            # Ensure attribute exists in case of weird states
            self.active_worker_ids = set([run_id])

        # Change the convert button to a Halt button so user can cancel
        try:
            # keep it enabled so the user can click Halt
            self.convert_button.setText("Halt conversion")
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.clicked.connect(self.halt_conversion)
        except Exception:
            pass

        # Keep cleanup disabled while conversion is in progress
        self.cleanup_button.setEnabled(False)
        # Disable 3D features during calculation
        self._enable_3d_features(False)
        self.statusBar().showMessage("Calculating 3D structure...")
        self.plotter.clear() 
        bg_color_hex = self.settings.get('background_color', '#919191')
        bg_qcolor = QColor(bg_color_hex)
        
        if bg_qcolor.isValid():
            luminance = bg_qcolor.toHsl().lightness()
            text_color = 'black' if luminance > 128 else 'white'
        else:
            text_color = 'white'
        
        text_actor = self.plotter.add_text(
            "Calculating...",
            position='lower_right',
            font_size=15,
            color=text_color,
            name='calculating_text'
        )
        # Keep a reference so we can reliably remove the text actor later
        try:
            self._calculating_text_actor = text_actor
        except Exception:
            # Best-effort: if storing fails, ignore — cleanup will still attempt renderer removal
            pass
        text_actor.GetTextProperty().SetOpacity(1)
        self.plotter.render()
        # Emit skip flag so the worker can ignore sanitization errors if user requested
        # Determine conversion_mode from settings (default: 'fallback').
        # If the user invoked conversion via the right-click menu, a temporary
        # override may be set on self._temp_conv_mode and should be used once.
        conv_mode = getattr(self, '_temp_conv_mode', None)
        if conv_mode:
            try:
                del self._temp_conv_mode
            except Exception:
                try:
                    delattr(self, '_temp_conv_mode')
                except Exception:
                    pass
        else:
            conv_mode = self.settings.get('3d_conversion_mode', 'fallback')

        # Allow a temporary optimization method override as well (used when
        # Optimize 3D is invoked via right-click menu). Do not persist here.
        opt_method = getattr(self, '_temp_optimization_method', None) or self.optimization_method
        if hasattr(self, '_temp_optimization_method'):
            try:
                del self._temp_optimization_method
            except Exception:
                try:
                    delattr(self, '_temp_optimization_method')
                except Exception:
                    pass

        options = {'conversion_mode': conv_mode, 'optimization_method': opt_method}
        # Attach the run id so the worker and main thread can correlate
        try:
            # Attach the concrete run id rather than the single waiting id
            options['worker_id'] = run_id
        except Exception:
            pass

        # Create a fresh CalculationWorker + QThread for this run so multiple
        # conversions can execute in parallel. The worker will be cleaned up
        # automatically after it finishes/errors.
        try:
            thread = QThread()
            worker = CalculationWorker()
            # Share the halt_ids set so user can request cancellation
            try:
                worker.halt_ids = self.halt_ids
            except Exception:
                pass

            worker.moveToThread(thread)

            # Forward status signals to main window handlers
            try:
                worker.status_update.connect(self.update_status_bar)
            except Exception:
                pass

            # When the worker finishes, call existing handler and then clean up
            def _on_worker_finished(result, w=worker, t=thread):
                try:
                    # deliver result to existing handler
                    self.on_calculation_finished(result)
                finally:
                    # Clean up signal connections to avoid stale references
                    # worker used its own start_work signal; no shared-signal
                    # disconnect necessary here.
                    # Remove thread from active threads list
                    try:
                        self._active_calc_threads.remove(t)
                    except Exception:
                        pass
                    try:
                        # ask thread to quit; it will finish as worker returns
                        t.quit()
                    except Exception:
                        pass
                    try:
                        # ensure thread object is deleted when finished
                        t.finished.connect(t.deleteLater)
                    except Exception:
                        pass
                    try:
                        # schedule worker deletion
                        w.deleteLater()
                    except Exception:
                        pass

            # When the worker errors (or halts), call existing handler and then clean up
            def _on_worker_error(error_msg, w=worker, t=thread):
                try:
                    # deliver error to existing handler
                    self.on_calculation_error(error_msg)
                finally:
                    # Clean up signal connections to avoid stale references
                    # worker used its own start_work signal; no shared-signal
                    # disconnect necessary here.
                    # Remove thread from active threads list
                    try:
                        self._active_calc_threads.remove(t)
                    except Exception:
                        pass
                    try:
                        # ask thread to quit; it will finish as worker returns
                        t.quit()
                    except Exception:
                        pass
                    try:
                        # ensure thread object is deleted when finished
                        t.finished.connect(t.deleteLater)
                    except Exception:
                        pass
                    try:
                        # schedule worker deletion
                        w.deleteLater()
                    except Exception:
                        pass

            try:
                worker.error.connect(_on_worker_error)
            except Exception:
                pass

            try:
                worker.finished.connect(_on_worker_finished)
            except Exception:
                pass

            # Start the thread
            thread.start()

            # Start the worker calculation via the worker's own start_work signal
            # (queued to the worker thread). Capture variables into lambda defaults
            # to avoid late-binding issues.
            QTimer.singleShot(10, lambda w=worker, m=mol_block, o=options: w.start_work.emit(m, o))

            # Track the thread so it isn't immediately garbage-collected (diagnostics)
            try:
                self._active_calc_threads.append(thread)
            except Exception:
                pass
        except Exception as e:
            # Fall back: if thread/worker creation failed, create a local
            # worker and start it (runs in main thread). This preserves
            # functionality without relying on the shared MainWindow signal.
            try:
                fallback_worker = CalculationWorker()
                QTimer.singleShot(10, lambda w=fallback_worker, m=mol_block, o=options: w.start_work.emit(m, o))
            except Exception:
                # surface the original error via existing UI path
                self.on_calculation_error(str(e))

        # 状態をUndo履歴に保存
        self.push_undo_state()
        self.update_chiral_labels()
        
        self.view_2d.setFocus()



    def halt_conversion(self):
        """User requested to halt the in-progress conversion.

        This will mark the current waiting_worker_id as halted (added to halt_ids),
        clear the waiting_worker_id, and immediately restore the UI (button text
        and handlers). The worker thread will observe halt_ids and should stop.
        """
        try:
            # Halt all currently-active workers by adding their ids to halt_ids
            wids_to_halt = set(getattr(self, 'active_worker_ids', set()))
            if wids_to_halt:
                try:
                    self.halt_ids.update(wids_to_halt)
                except Exception:
                    pass

            # Clear the active set immediately so UI reflects cancellation
            try:
                if hasattr(self, 'active_worker_ids'):
                    self.active_worker_ids.clear()
            except Exception:
                pass

            # Restore UI immediately
            try:
                try:
                    self.convert_button.clicked.disconnect()
                except Exception:
                    pass
                self.convert_button.setText("Convert 2D to 3D")
                self.convert_button.clicked.connect(self.trigger_conversion)
                self.convert_button.setEnabled(True)
            except Exception:
                pass

            try:
                self.cleanup_button.setEnabled(True)
            except Exception:
                pass

            # Remove any calculating text actor if present
            try:
                actor = getattr(self, '_calculating_text_actor', None)
                if actor is not None:
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
            except Exception:
                pass

            # Give immediate feedback
            self.statusBar().showMessage("3D conversion halted. Waiting for the thread to finish")
        except Exception:
            pass



    def check_chemistry_problems_fallback(self):
        """RDKit変換が失敗した場合の化学的問題チェック（独自実装）"""
        try:
            # 既存のフラグをクリア
            self.scene.clear_all_problem_flags()
            
            # 簡易的な化学的問題チェック
            problem_atoms = []
            
            for atom_id, atom_data in self.data.atoms.items():
                atom_item = atom_data.get('item')
                if not atom_item:
                    continue
                
                symbol = atom_data['symbol']
                charge = atom_data.get('charge', 0)
                
                # 結合数を計算
                bond_count = 0
                for (id1, id2), bond_data in self.data.bonds.items():
                    if id1 == atom_id or id2 == atom_id:
                        bond_count += bond_data.get('order', 1)
                
                # 基本的な価数チェック
                is_problematic = False
                if symbol == 'C' and bond_count > 4:
                    is_problematic = True
                elif symbol == 'N' and bond_count > 3 and charge == 0:
                    is_problematic = True
                elif symbol == 'O' and bond_count > 2 and charge == 0:
                    is_problematic = True
                elif symbol == 'H' and bond_count > 1:
                    is_problematic = True
                elif symbol in ['F', 'Cl', 'Br', 'I'] and bond_count > 1 and charge == 0:
                    is_problematic = True
                
                if is_problematic:
                    problem_atoms.append(atom_item)
            
            if problem_atoms:
                # 問題のある原子に赤枠を設定
                for atom_item in problem_atoms:
                    atom_item.has_problem = True
                    atom_item.update()
                
                self.statusBar().showMessage(f"Error: {len(problem_atoms)} chemistry problem(s) found (valence issues).")
            else:
                self.statusBar().showMessage("Error: Invalid chemical structure (RDKit conversion failed).")
            
            self.scene.clearSelection()
            self.view_2d.setFocus()
            
        except Exception as e:
            print(f"Error in fallback chemistry check: {e}")
            self.statusBar().showMessage("Error: Invalid chemical structure.")
            self.view_2d.setFocus()



    def optimize_3d_structure(self):
        """現在の3D分子構造を力場で最適化する"""
        if not self.current_mol:
            self.statusBar().showMessage("No 3D molecule to optimize.")
            return

        # If a prior chemical/sanitization check was attempted and failed, do not run optimization
        if getattr(self, 'chem_check_tried', False) and getattr(self, 'chem_check_failed', False):
            self.statusBar().showMessage("3D optimization disabled: molecule failed chemical sanitization.")
            # Ensure the Optimize 3D button is disabled to reflect this
            if hasattr(self, 'optimize_3d_button'):
                try:
                    self.optimize_3d_button.setEnabled(False)
                except Exception:
                    pass
            return

        self.statusBar().showMessage("Optimizing 3D structure...")
        QApplication.processEvents() # UIの更新を確実に行う

        try:
            # Allow a temporary optimization method override (right-click menu)
            method = getattr(self, '_temp_optimization_method', None) or getattr(self, 'optimization_method', 'MMFF_RDKIT')
            # Clear temporary override if present
            if hasattr(self, '_temp_optimization_method'):
                try:
                    del self._temp_optimization_method
                except Exception:
                    try:
                        delattr(self, '_temp_optimization_method')
                    except Exception:
                        pass
            method = method.upper() if method else 'MMFF_RDKIT'
            # 事前チェック：コンフォーマがあるか
            if self.current_mol.GetNumConformers() == 0:
                self.statusBar().showMessage("No conformer found: cannot optimize. Embed molecule first.")
                return
            if method in ('MMFF_RDKIT', 'MMFF94_RDKIT'):
                try:
                    # Choose concrete mmffVariant string
                    mmff_variant = "MMFF94s" if method == 'MMFF_RDKIT' else "MMFF94"
                    res = AllChem.MMFFOptimizeMolecule(self.current_mol, maxIters=4000, mmffVariant=mmff_variant)
                    if res != 0:
                        # 非収束や何らかの問題が起きた可能性 -> ForceField API で詳細に試す
                        try:
                            mmff_props = AllChem.MMFFGetMoleculeProperties(self.current_mol)
                            ff = AllChem.MMFFGetMoleculeForceField(self.current_mol, mmff_props, confId=0)
                            ff_ret = ff.Minimize(maxIts=4000)
                            if ff_ret != 0:
                                self.statusBar().showMessage(f"{mmff_variant} minimize returned non-zero status: {ff_ret}")
                                return
                        except Exception as e:
                            self.statusBar().showMessage(f"{mmff_variant} parameterization/minimize failed: {e}")
                            return
                except Exception as e:
                    self.statusBar().showMessage(f"{mmff_variant} (RDKit) optimization error: {e}")
                    return
            elif method == 'UFF_RDKIT':
                try:
                    res = AllChem.UFFOptimizeMolecule(self.current_mol, maxIters=4000)
                    if res != 0:
                        try:
                            ff = AllChem.UFFGetMoleculeForceField(self.current_mol, confId=0)
                            ff_ret = ff.Minimize(maxIts=4000)
                            if ff_ret != 0:
                                self.statusBar().showMessage(f"UFF minimize returned non-zero status: {ff_ret}")
                                return
                        except Exception as e:
                            self.statusBar().showMessage(f"UFF parameterization/minimize failed: {e}")
                            return
                except Exception as e:
                    self.statusBar().showMessage(f"UFF (RDKit) optimization error: {e}")
                    return
            else:
                self.statusBar().showMessage("Selected optimization method is not available. Use MMFF94 (RDKit) or UFF (RDKit).")
                return
        except Exception as e:
            self.statusBar().showMessage(f"3D optimization error: {e}")
        
        # 最適化後の構造で3Dビューを再描画
        try:
            # Remember which concrete optimizer variant succeeded so it
            # can be saved with the project. Normalize internal flags to
            # a human-friendly label: MMFF94s, MMFF94, or UFF.
            try:
                norm_method = None
                m = method.upper() if method else None
                if m in ('MMFF_RDKIT', 'MMFF94_RDKIT'):
                    # The code above uses mmffVariant="MMFF94s" when
                    # method == 'MMFF_RDKIT' and "MMFF94" otherwise.
                    norm_method = 'MMFF94s' if m == 'MMFF_RDKIT' else 'MMFF94'
                elif m == 'UFF_RDKIT' or m == 'UFF':
                    norm_method = 'UFF'
                else:
                    norm_method = getattr(self, 'optimization_method', None)

                # store for later serialization
                if norm_method:
                    self.last_successful_optimization_method = norm_method
            except Exception:
                pass
            # 3D最適化後は3D座標から立体化学を再計算（2回目以降は3D優先）
            if self.current_mol.GetNumConformers() > 0:
                Chem.AssignAtomChiralTagsFromStructure(self.current_mol, confId=0)
            self.update_chiral_labels() # キラル中心のラベルも更新
        except Exception:
            pass
            
        self.draw_molecule_3d(self.current_mol)
        
        # Show which method was used in the status bar (prefer human-readable label).
        # Prefer the actual method used during this run (last_successful_optimization_method
        # set earlier), then any temporary/local override used for this call (method),
        # and finally the persisted preference (self.optimization_method).
        try:
            used_method = (
                getattr(self, 'last_successful_optimization_method', None)
                or locals().get('method', None)
                or getattr(self, 'optimization_method', None)
            )
            used_label = None
            if used_method:
                # opt3d_method_labels keys are stored upper-case; normalize for lookup
                used_label = (getattr(self, 'opt3d_method_labels', {}) or {}).get(str(used_method).upper(), used_method)
        except Exception:
            used_label = None

        if used_label:
            self.statusBar().showMessage(f"3D structure optimization successful. Method: {used_label}")
        else:
            self.statusBar().showMessage("3D structure optimization successful.")
        self.push_undo_state() # Undo履歴に保存
        self.view_2d.setFocus()



    def on_calculation_finished(self, result):
        # Accept either (worker_id, mol) tuple or legacy single mol arg
        worker_id = None
        mol = None
        try:
            if isinstance(result, tuple) and len(result) == 2:
                worker_id, mol = result
            else:
                mol = result
        except Exception:
            mol = result

        # If this finished result is from a stale/halting run, discard it
        try:
            if worker_id is not None:
                # If this worker_id is not in the active set, it's stale/halting
                if worker_id not in getattr(self, 'active_worker_ids', set()):
                    # Cleanup calculating UI and ignore
                    try:
                        actor = getattr(self, '_calculating_text_actor', None)
                        if actor is not None:
                            if hasattr(self.plotter, 'remove_actor'):
                                try:
                                    self.plotter.remove_actor(actor)
                                except Exception:
                                    pass
                            else:
                                if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                    try:
                                        self.plotter.renderer.RemoveActor(actor)
                                    except Exception:
                                        pass
                            try:
                                delattr(self, '_calculating_text_actor')
                            except Exception:
                                try:
                                    del self._calculating_text_actor
                                except Exception:
                                    pass
                    except Exception:
                        pass
                    # Ensure Convert button is restored
                    try:
                        try:
                            self.convert_button.clicked.disconnect()
                        except Exception:
                            pass
                        self.convert_button.setText("Convert 2D to 3D")
                        self.convert_button.clicked.connect(self.trigger_conversion)
                        self.convert_button.setEnabled(True)
                    except Exception:
                        pass
                    try:
                        self.cleanup_button.setEnabled(True)
                    except Exception:
                        pass
                    self.statusBar().showMessage("Ignored result from stale conversion.")
                    return
        except Exception:
            pass

        # Remove the finished worker id from the active set and any halt set
        try:
            if worker_id is not None:
                try:
                    self.active_worker_ids.discard(worker_id)
                except Exception:
                    pass
            # Also remove id from halt set if present
            if worker_id is not None:
                try:
                    if worker_id in getattr(self, 'halt_ids', set()):
                        try:
                            self.halt_ids.discard(worker_id)
                        except Exception:
                            pass
                except Exception:
                    pass
        except Exception:
            pass

        self.dragged_atom_info = None
        self.current_mol = mol
        self.is_xyz_derived = False  # 2Dから生成した3D構造はXYZ由来ではない
        # Record the optimization method used for this conversion if available.
        try:
            opt_method = None
            try:
                # Worker or molecule may have attached a prop with the used method
                if hasattr(mol, 'HasProp') and mol is not None:
                    try:
                        if mol.HasProp('_pme_optimization_method'):
                            opt_method = mol.GetProp('_pme_optimization_method')
                    except Exception:
                        # not all Mol objects support HasProp/GetProp safely
                        pass
            except Exception:
                pass
            if not opt_method:
                opt_method = getattr(self, 'optimization_method', None)
            # normalize common forms
            if opt_method:
                om = str(opt_method).upper()
                if 'MMFF94S' in om or 'MMFF_RDKIT' in om:
                    self.last_successful_optimization_method = 'MMFF94s'
                elif 'MMFF94' in om:
                    self.last_successful_optimization_method = 'MMFF94'
                elif 'UFF' in om:
                    self.last_successful_optimization_method = 'UFF'
                else:
                    # store raw value otherwise
                    self.last_successful_optimization_method = opt_method
        except Exception:
            # non-fatal
            pass
        
        # 原子プロパティを復元（ワーカープロセスで失われたため）
        if hasattr(self, 'original_atom_properties'):
            for i, original_id in self.original_atom_properties.items():
                if i < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(i)
                    atom.SetIntProp("_original_atom_id", original_id)
        
        # 原子IDマッピングを作成
        self.create_atom_id_mapping()
        
        # キラル中心を初回変換時は2Dの立体情報を考慮して設定
        try:
            if mol.GetNumConformers() > 0:
                # 初回変換では、2Dで設定したwedge/dashボンドの立体情報を保持
                # 立体化学の割り当てを行うが、既存の2D立体情報を尊重
                Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
            
            self.update_chiral_labels()
        except Exception:
            # 念のためエラーを握り潰して UI を壊さない
            pass

        self.draw_molecule_3d(mol)
        
        # 複数分子の場合、衝突検出と配置調整を実行
        try:
            frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
            if len(frags) > 1:
                self.statusBar().showMessage(f"Detecting collisions among {len(frags)} molecules...")
                QApplication.processEvents()
                self.adjust_molecule_positions_to_avoid_collisions(mol, frags)
                self.draw_molecule_3d(mol)
                self.update_chiral_labels()
                self.statusBar().showMessage(f"{len(frags)} molecules converted with collision avoidance.")
        except Exception as e:
            print(f"Warning: Collision detection failed: {e}")
            # 衝突検出に失敗してもエラーにはしない

        # Ensure any 'Calculating...' text is removed and the plotter is refreshed
        try:
            actor = getattr(self, '_calculating_text_actor', None)
            if actor is not None:
                try:
                    # Prefer plotter API if available
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            # Some pyvista versions use renderer.RemoveActor
                            if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                try:
                                    self.plotter.renderer.RemoveActor(actor)
                                except Exception:
                                    pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                finally:
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
            # Re-render to ensure the UI updates immediately
            try:
                self.plotter.render()
            except Exception:
                pass
        except Exception:
            pass

        #self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
        # Restore Convert button text/handler in case it was changed to Halt
        try:
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.setText("Convert 2D to 3D")
            self.convert_button.clicked.connect(self.trigger_conversion)
        except Exception:
            pass
        self.push_undo_state()
        self.view_2d.setFocus()
        self.cleanup_button.setEnabled(True)
        
        # 3D関連機能を統一的に有効化
        self._enable_3d_features(True)
            
        self.plotter.reset_camera()
        
        # 3D原子情報ホバー表示を再設定
        self.setup_3d_hover()
        
        # メニューテキストと状態を更新
        self.update_atom_id_menu_text()
        self.update_atom_id_menu_state()



    def create_atom_id_mapping(self):
        """2D原子IDから3D RDKit原子インデックスへのマッピングを作成する（RDKitの原子プロパティ使用）"""
        if not self.current_mol:
            return
            
        self.atom_id_to_rdkit_idx_map = {}
        
        # RDKitの原子プロパティから直接マッピングを作成
        for i in range(self.current_mol.GetNumAtoms()):
            rdkit_atom = self.current_mol.GetAtomWithIdx(i)
            try:
                original_atom_id = rdkit_atom.GetIntProp("_original_atom_id")
                self.atom_id_to_rdkit_idx_map[original_atom_id] = i
            except KeyError:
                # プロパティが設定されていない場合（外部ファイル読み込み時など）
                continue



    def on_calculation_error(self, result):
        """ワーカースレッドからのエラー（またはHalt）を処理する"""
        worker_id = None
        error_message = ""
        try:
            if isinstance(result, tuple) and len(result) == 2:
                worker_id, error_message = result
            else:
                error_message = str(result)
        except Exception:
            error_message = str(result)

        # If this error is from a stale/previous worker (not in active set), ignore it.
        if worker_id is not None and worker_id not in getattr(self, 'active_worker_ids', set()):
            # Stale/late error from a previously-halted worker; ignore to avoid clobbering newer runs
            print(f"Ignored stale error from worker {worker_id}: {error_message}")
            return

        # Clear temporary plotter content and remove calculating text if present
        try:
            self.plotter.clear()
        except Exception:
            pass

        # Also attempt to explicitly remove the calculating text actor if it was stored
        try:
            actor = getattr(self, '_calculating_text_actor', None)
            if actor is not None:
                try:
                    if hasattr(self.plotter, 'remove_actor'):
                        try:
                            self.plotter.remove_actor(actor)
                        except Exception:
                            if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                                try:
                                    self.plotter.renderer.RemoveActor(actor)
                                except Exception:
                                    pass
                    else:
                        if hasattr(self.plotter, 'renderer') and self.plotter.renderer:
                            try:
                                self.plotter.renderer.RemoveActor(actor)
                            except Exception:
                                pass
                finally:
                    try:
                        delattr(self, '_calculating_text_actor')
                    except Exception:
                        try:
                            del self._calculating_text_actor
                        except Exception:
                            pass
        except Exception:
            pass

        self.dragged_atom_info = None
        # Remove this worker id from active set (error belongs to this worker)
        try:
            if worker_id is not None:
                try:
                    self.active_worker_ids.discard(worker_id)
                except Exception:
                    pass
        except Exception:
            pass

        # If this error was caused by an intentional halt and the main thread
        # already cleared waiting_worker_id earlier for other reasons, suppress the error noise.
        try:
            low = (error_message or '').lower()
            # If a halt message and there are no active workers left, the user
            # already saw the halt message — suppress duplicate noise.
            if 'halt' in low and not getattr(self, 'active_worker_ids', set()):
                return
        except Exception:
            pass

        self.statusBar().showMessage(f"Error: {error_message}")
        
        try:
            self.cleanup_button.setEnabled(True)
        except Exception:
            pass
        try:
            # Restore Convert button text/handler
            try:
                self.convert_button.clicked.disconnect()
            except Exception:
                pass
            self.convert_button.setText("Convert 2D to 3D")
            self.convert_button.clicked.connect(self.trigger_conversion)
            self.convert_button.setEnabled(True)
        except Exception:
            pass

        # On calculation error we should NOT enable 3D-only features.
        # Explicitly disable Optimize and Export so the user can't try to operate
        # on an invalid or missing 3D molecule.
        try:
            if hasattr(self, 'optimize_3d_button'):
                self.optimize_3d_button.setEnabled(False)
        except Exception:
            pass
        try:
            if hasattr(self, 'export_button'):
                self.export_button.setEnabled(False)
        except Exception:
            pass

        # Keep 3D feature buttons disabled to avoid inconsistent UI state
        try:
            self._enable_3d_features(False)
        except Exception:
            pass

        # Keep 3D edit actions disabled (no molecule to edit)
        try:
            self._enable_3d_edit_actions(False)
        except Exception:
            pass
        # Some menu items are explicitly disabled on error
        try:
            if hasattr(self, 'analysis_action'):
                self.analysis_action.setEnabled(False)
        except Exception:
            pass
        try:
            if hasattr(self, 'edit_3d_action'):
                self.edit_3d_action.setEnabled(False)
        except Exception:
            pass

        # Force a UI refresh
        try:
            self.plotter.render()
        except Exception:
            pass

        # Ensure focus returns to 2D editor
        self.view_2d.setFocus()

