import numpy as np

from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera

from PyQt6.QtWidgets import QApplication

from PyQt6.QtCore import (
    Qt
)

try:
    from .constants import pt
except Exception:
    from modules.constants import pt
try:
    from .move_group_dialog import MoveGroupDialog
except Exception:
    from modules.move_group_dialog import MoveGroupDialog

class CustomInteractorStyle(vtkInteractorStyleTrackballCamera):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        # カスタム状態を管理するフラグを一つに絞ります
        self._is_dragging_atom = False
        # undoスタックのためのフラグ
        self.is_dragging = False
        # 回転操作を検出するためのフラグ
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None

        self.AddObserver("LeftButtonPressEvent", self.on_left_button_down)
        self.AddObserver("RightButtonPressEvent", self.on_right_button_down)
        self.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.AddObserver("LeftButtonReleaseEvent", self.on_left_button_up)
        self.AddObserver("RightButtonReleaseEvent", self.on_right_button_up)

    def on_left_button_down(self, obj, event):
        """
        クリック時の処理を振り分けます。
        原子を掴めた場合のみカスタム動作に入り、それ以外は親クラス（カメラ回転）に任せます。
        """
        mw = self.main_window
        
        # 前回のドラッグ状態をクリア（トリプルクリック/ダブルクリック対策）
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_moved_during_drag = False
        self._mouse_press_pos = None
        
        # Move Groupダイアログが開いている場合の処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and move_group_dialog.group_atoms:
            # グループが選択されている場合、グループドラッグ処理
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)
                
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5
                        
                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)
            
            # グループ内の原子がクリックされた場合
            if clicked_atom_idx is not None:
                if clicked_atom_idx in move_group_dialog.group_atoms:
                    # 既存グループ内の原子 - ドラッグ準備
                    move_group_dialog._is_dragging_group_vtk = True
                    move_group_dialog._drag_atom_idx = clicked_atom_idx
                    move_group_dialog._drag_start_pos = click_pos
                    move_group_dialog._mouse_moved = False
                    # 初期位置を保存
                    move_group_dialog._initial_positions = {}
                    conf = mw.current_mol.GetConformer()
                    for atom_idx in move_group_dialog.group_atoms:
                        pos = conf.GetAtomPosition(atom_idx)
                        move_group_dialog._initial_positions[atom_idx] = np.array([pos.x, pos.y, pos.z])
                    mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                    return  # カメラ回転を無効化
                else:
                    # グループ外の原子をクリック - BFS/DFSで連結成分を探索
                    visited = set()
                    queue = [clicked_atom_idx]
                    visited.add(clicked_atom_idx)
                    
                    while queue:
                        current_idx = queue.pop(0)
                        for bond_idx in range(mw.current_mol.GetNumBonds()):
                            bond = mw.current_mol.GetBondWithIdx(bond_idx)
                            begin_idx = bond.GetBeginAtomIdx()
                            end_idx = bond.GetEndAtomIdx()
                            
                            if begin_idx == current_idx and end_idx not in visited:
                                visited.add(end_idx)
                                queue.append(end_idx)
                            elif end_idx == current_idx and begin_idx not in visited:
                                visited.add(begin_idx)
                                queue.append(begin_idx)
                    
                    # Ctrlキーが押されている場合のみ複数グループ選択
                    is_ctrl_pressed = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier)
                    
                    if is_ctrl_pressed:
                        # Ctrl + クリック: 追加または解除
                        if visited.issubset(move_group_dialog.group_atoms):
                            # すでに選択されている - 解除
                            move_group_dialog.group_atoms -= visited
                        else:
                            # 新しいグループを追加
                            move_group_dialog.group_atoms |= visited
                    else:
                        # 通常のクリック: 既存の選択を置き換え
                        move_group_dialog.group_atoms = visited.copy()
                    
                    move_group_dialog.selected_atoms.add(clicked_atom_idx)
                    move_group_dialog.show_atom_labels()
                    move_group_dialog.update_display()
                    return
            else:
                # 原子以外をクリック - 全選択を解除
                move_group_dialog.group_atoms.clear()
                move_group_dialog.selected_atoms.clear()
                move_group_dialog.clear_atom_labels()
                move_group_dialog.update_display()
                # カメラ回転を許可
                super(CustomInteractorStyle, self).OnLeftButtonDown()
                return
        
        is_temp_mode = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.AltModifier)
        is_edit_active = mw.is_3d_edit_mode or is_temp_mode
        
        # Ctrl+クリックで原子選択（3D編集用）
        is_ctrl_click = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier)

        # 測定モードが有効な場合の処理
        if mw.measurement_mode and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            self._mouse_press_pos = click_pos  # マウスプレス位置を記録
            self._mouse_moved_during_drag = False  # 移動フラグをリセット
            
            picker = mw.plotter.picker
            
            # 通常のピック処理を実行
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            # 原子がクリックされた場合のみ特別処理
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)

                # 範囲チェックを追加
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # クリック閾値チェック
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            mw.handle_measurement_atom_selection(int(closest_atom_idx))
                            return  # 原子選択処理完了、カメラ回転は無効
            
            # 測定モードで原子以外をクリックした場合は計測選択をクリア
            # ただし、これは通常のカメラ回転も許可する
            self._is_dragging_atom = False
            super().OnLeftButtonDown()
            return
        
        # Ctrl+クリックの原子選択機能は無効化（Move Group機能で代替）
        # if is_ctrl_click and mw.current_mol:
        #     ... (無効化)

        # 3D分子(mw.current_mol)が存在する場合のみ、原子の選択処理を実行
        if is_edit_active and mw.current_mol:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)

            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)

                # 範囲チェックを追加
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    # RDKitのMolオブジェクトから原子を安全に取得
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5

                        if distances[closest_atom_idx] < click_threshold:
                            # 原子を掴むことに成功した場合
                            self._is_dragging_atom = True
                        self.is_dragging = False 
                        mw.dragged_atom_info = {'id': int(closest_atom_idx)}
                        mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                        return  # 親クラスのカメラ回転を呼ばない

        self._is_dragging_atom = False
        super().OnLeftButtonDown()

    def on_right_button_down(self, obj, event):
        """
        右クリック時の処理。Move Groupダイアログが開いている場合はグループ回転を開始。
        """
        mw = self.main_window
        
        # Move Groupダイアログが開いているか確認
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and move_group_dialog.group_atoms:
            # グループが選択されている場合、回転ドラッグを開始
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            clicked_atom_idx = None
            if picker.GetActor() is mw.atom_actor:
                picked_position = np.array(picker.GetPickPosition())
                distances = np.linalg.norm(mw.atom_positions_3d - picked_position, axis=1)
                closest_atom_idx = np.argmin(distances)
                
                if 0 <= closest_atom_idx < mw.current_mol.GetNumAtoms():
                    atom = mw.current_mol.GetAtomWithIdx(int(closest_atom_idx))
                    if atom:
                        atomic_num = atom.GetAtomicNum()
                        vdw_radius = pt.GetRvdw(atomic_num)
                        click_threshold = vdw_radius * 1.5
                        
                        if distances[closest_atom_idx] < click_threshold:
                            clicked_atom_idx = int(closest_atom_idx)
            
            # グループ内の原子がクリックされた場合、回転ドラッグを開始
            if clicked_atom_idx is not None and clicked_atom_idx in move_group_dialog.group_atoms:
                move_group_dialog._is_rotating_group_vtk = True
                move_group_dialog._rotation_start_pos = click_pos
                move_group_dialog._rotation_mouse_moved = False
                move_group_dialog._rotation_atom_idx = clicked_atom_idx  # 掴んだ原子を記録
                
                # 初期位置と重心を保存
                move_group_dialog._initial_positions = {}
                conf = mw.current_mol.GetConformer()
                centroid = np.zeros(3)
                for atom_idx in move_group_dialog.group_atoms:
                    pos = conf.GetAtomPosition(atom_idx)
                    pos_array = np.array([pos.x, pos.y, pos.z])
                    move_group_dialog._initial_positions[atom_idx] = pos_array
                    centroid += pos_array
                centroid /= len(move_group_dialog.group_atoms)
                move_group_dialog._group_centroid = centroid
                
                mw.plotter.setCursor(Qt.CursorShape.ClosedHandCursor)
                return  # カメラ回転を無効化
        
        # 通常の右クリック処理
        super().OnRightButtonDown()

    def on_mouse_move(self, obj, event):
        """
        マウス移動時の処理。原子ドラッグ中か、それ以外（カメラ回転＋ホバー）かをハンドリングします。
        """
        mw = self.main_window
        
        # Move Groupダイアログのドラッグ処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and getattr(move_group_dialog, '_is_dragging_group_vtk', False):
            # グループをドラッグ中 - 移動距離を記録するのみ
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()
            
            dx = current_pos[0] - move_group_dialog._drag_start_pos[0]
            dy = current_pos[1] - move_group_dialog._drag_start_pos[1]
            
            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._mouse_moved = True
            
            return  # カメラ回転を無効化
        
        # グループ回転中の処理
        if move_group_dialog and getattr(move_group_dialog, '_is_rotating_group_vtk', False):
            interactor = self.GetInteractor()
            current_pos = interactor.GetEventPosition()
            
            dx = current_pos[0] - move_group_dialog._rotation_start_pos[0]
            dy = current_pos[1] - move_group_dialog._rotation_start_pos[1]
            
            if abs(dx) > 2 or abs(dy) > 2:
                move_group_dialog._rotation_mouse_moved = True
            
            return  # カメラ回転を無効化
        
        interactor = self.GetInteractor()

        # マウス移動があったことを記録
        if self._mouse_press_pos is not None:
            current_pos = interactor.GetEventPosition()
            if abs(current_pos[0] - self._mouse_press_pos[0]) > 3 or abs(current_pos[1] - self._mouse_press_pos[1]) > 3:
                self._mouse_moved_during_drag = True

        if self._is_dragging_atom and mw.dragged_atom_info is not None:
            # カスタムの原子ドラッグ処理
            self.is_dragging = True
            atom_id = mw.dragged_atom_info['id']
            # We intentionally do NOT update visible coordinates or the
            # authoritative atom position during mouse-move while dragging.
            # The UX requirement here is that atoms need not visibly move
            # while the mouse is being dragged. Compute and apply the final
            # world-coordinate only once on mouse release (on_left_button_up).
            # Keep minimal state: mark that a drag occurred (is_dragging)
            # and allow the release handler to compute the final position.
            # This avoids duplicate updates and simplifies event ordering.
        else:
            # カメラ回転処理を親クラスに任せます
            super().OnMouseMove()

            # その後、カーソルの表示を更新します
            is_edit_active = mw.is_3d_edit_mode or interactor.GetAltKey()
            if is_edit_active:
                # 編集がアクティブな場合のみ、原子のホバーチェックを行う
                atom_under_cursor = False
                click_pos = interactor.GetEventPosition()
                picker = mw.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
                if picker.GetActor() is mw.atom_actor:
                    atom_under_cursor = True

                if atom_under_cursor:
                    mw.plotter.setCursor(Qt.CursorShape.OpenHandCursor)
                else:
                    mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            else:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)

    def on_left_button_up(self, obj, event):
        """
        クリック終了時の処理。状態をリセットします。
        """
        mw = self.main_window
        
        # Move Groupダイアログのドラッグ終了処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        # ダブルクリック/トリプルクリックで状態が混乱するのを防ぐ（Move Group用）
        if move_group_dialog:
            if getattr(move_group_dialog, '_is_dragging_group_vtk', False) and not getattr(move_group_dialog, '_mouse_moved', False):
                # ドラッグしていない状態で複数クリックされた場合は状態をリセット
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, '_initial_positions'):
                    delattr(move_group_dialog, '_initial_positions')
        
        if move_group_dialog and getattr(move_group_dialog, '_is_dragging_group_vtk', False):
            if getattr(move_group_dialog, '_mouse_moved', False):
                # ドラッグが実行された - リリース時に座標を更新
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()
                    
                    # ドラッグ原子の初期位置
                    drag_atom_initial_pos = move_group_dialog._initial_positions[move_group_dialog._drag_atom_idx]
                    
                    # スクリーン座標からワールド座標への変換
                    renderer.SetWorldPoint(drag_atom_initial_pos[0], drag_atom_initial_pos[1], drag_atom_initial_pos[2], 1.0)
                    renderer.WorldToDisplay()
                    display_coords = renderer.GetDisplayPoint()
                    
                    new_display_pos = (current_pos[0], current_pos[1], display_coords[2])
                    renderer.SetDisplayPoint(new_display_pos[0], new_display_pos[1], new_display_pos[2])
                    renderer.DisplayToWorld()
                    new_world_coords = renderer.GetWorldPoint()
                    
                    # 移動ベクトル
                    translation_vector = np.array([
                        new_world_coords[0] - drag_atom_initial_pos[0],
                        new_world_coords[1] - drag_atom_initial_pos[1],
                        new_world_coords[2] - drag_atom_initial_pos[2]
                    ])
                    
                    # グループ全体を移動
                    for atom_idx in move_group_dialog.group_atoms:
                        initial_pos = move_group_dialog._initial_positions[atom_idx]
                        new_pos = initial_pos + translation_vector
                        conf.SetAtomPosition(atom_idx, new_pos.tolist())
                        mw.atom_positions_3d[atom_idx] = new_pos
                    
                    # 3D表示を更新
                    mw.draw_molecule_3d(mw.current_mol)
                    mw.update_chiral_labels()
                    move_group_dialog.show_atom_labels()
                    mw.push_undo_state()
                except Exception as e:
                    print(f"Error finalizing group drag: {e}")
            else:
                # ドラッグがなかった = クリックのみ → トグル処理
                if hasattr(move_group_dialog, '_drag_atom_idx'):
                    clicked_atom = move_group_dialog._drag_atom_idx
                    try:
                        move_group_dialog.on_atom_picked(clicked_atom)
                    except Exception as e:
                        print(f"Error in toggle: {e}")
            
            # 状態をリセット（完全なクリーンアップ）
            move_group_dialog._is_dragging_group_vtk = False
            move_group_dialog._drag_start_pos = None
            move_group_dialog._mouse_moved = False
            if hasattr(move_group_dialog, '_initial_positions'):
                delattr(move_group_dialog, '_initial_positions')
            if hasattr(move_group_dialog, '_drag_atom_idx'):
                delattr(move_group_dialog, '_drag_atom_idx')
            
            # CustomInteractorStyleの状態もクリア
            self._is_dragging_atom = False
            self.is_dragging = False
            self._mouse_moved_during_drag = False
            self._mouse_press_pos = None
            
            try:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except Exception:
                pass
            return

        # 計測モードで、マウスが動いていない場合（つまりクリック）の処理
        if mw.measurement_mode and not self._mouse_moved_during_drag and self._mouse_press_pos is not None:
            click_pos = self.GetInteractor().GetEventPosition()
            picker = mw.plotter.picker
            picker.Pick(click_pos[0], click_pos[1], 0, mw.plotter.renderer)
            
            # 原子がクリックされていない場合は測定選択をクリア
            if picker.GetActor() is not mw.atom_actor:
                mw.clear_measurement_selection()

        if self._is_dragging_atom:
            # カスタムドラッグの後始末
            if self.is_dragging:
                if mw.current_mol and mw.current_mol.GetNumConformers() > 0:
                    try:
                        # Before applying conformer updates, compute the final
                        # world coordinates for the dragged atom based on the
                        # release pointer position. During the drag we did not
                        # update mw.atom_positions_3d (to keep the visuals
                        # static). Now compute the final position for the
                        # dragged atom and store it into mw.atom_positions_3d
                        # so the conformer update loop below will pick it up.
                        atom_id = None
                        try:
                            atom_id = mw.dragged_atom_info.get('id') if mw.dragged_atom_info else None
                        except Exception:
                            atom_id = None

                        if atom_id is not None:
                            try:
                                interactor = self.GetInteractor()
                                renderer = mw.plotter.renderer
                                current_display_pos = interactor.GetEventPosition()
                                conf = mw.current_mol.GetConformer()
                                # Use the atom's current 3D position to obtain a
                                # display-space depth (z) value, then replace the
                                # x/y with the pointer position to project back to
                                # world coordinates at that depth.
                                pos_3d = conf.GetAtomPosition(atom_id)
                                renderer.SetWorldPoint(pos_3d.x, pos_3d.y, pos_3d.z, 1.0)
                                renderer.WorldToDisplay()
                                display_coords = renderer.GetDisplayPoint()
                                new_display_pos = (current_display_pos[0], current_display_pos[1], display_coords[2])
                                renderer.SetDisplayPoint(new_display_pos[0], new_display_pos[1], new_display_pos[2])
                                renderer.DisplayToWorld()
                                new_world_coords_tuple = renderer.GetWorldPoint()
                                new_world_coords = list(new_world_coords_tuple)[:3]
                                # Ensure the container supports assignment
                                try:
                                    mw.atom_positions_3d[atom_id] = new_world_coords
                                except Exception:
                                    # If atom_positions_3d is immutable or shaped
                                    # differently, attempt a safe conversion.
                                    try:
                                        ap = list(mw.atom_positions_3d)
                                        ap[atom_id] = new_world_coords
                                        mw.atom_positions_3d = ap
                                    except Exception:
                                        pass
                            except Exception:
                                # If final-position computation fails, continue
                                # and apply whatever state is available.
                                pass

                        # Apply the (now updated) positions to the RDKit conformer
                        # exactly once. This ensures the conformer is
                        # authoritative and avoids double-moves.
                        conf = mw.current_mol.GetConformer()
                        for i in range(mw.current_mol.GetNumAtoms()):
                            try:
                                pos = mw.atom_positions_3d[i]
                                conf.SetAtomPosition(i, pos.tolist())
                            except Exception:
                                # Skip individual failures but continue applying
                                # other atom positions.
                                pass
                    except Exception:
                        # If applying positions fails, continue to redraw from
                        # whatever authoritative state is available.
                        pass

                    # Redraw once and push undo state
                    try:
                        mw.draw_molecule_3d(mw.current_mol)
                    except Exception:
                        pass
                    mw.push_undo_state()
            mw.dragged_atom_info = None
            # Refresh overlays and labels that depend on atom_positions_3d. Do
            # not overwrite mw.atom_positions_3d here — it already reflects the
            # positions the user dragged to. Only update dependent displays.
            try:
                mw.update_3d_selection_display()
            except Exception:
                pass
            try:
                mw.update_measurement_labels_display()
            except Exception:
                pass
            try:
                mw.update_2d_measurement_labels()
            except Exception:
                pass
            try:
                mw.show_all_atom_info()
            except Exception:
                pass
            except Exception:
                # Do not allow a failure here to interrupt release flow
                pass
        else:
            # カメラ回転の後始末を親クラスに任せます
            super().OnLeftButtonUp()

        # 状態をリセット（完全なクリーンアップ）
        self._is_dragging_atom = False
        self.is_dragging = False
        self._mouse_press_pos = None
        self._mouse_moved_during_drag = False
        
        # Move Group関連の状態もクリア
        try:
            move_group_dialog = None
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
            
            if move_group_dialog:
                move_group_dialog._is_dragging_group_vtk = False
                move_group_dialog._drag_start_pos = None
                move_group_dialog._mouse_moved = False
                if hasattr(move_group_dialog, '_initial_positions'):
                    delattr(move_group_dialog, '_initial_positions')
                if hasattr(move_group_dialog, '_drag_atom_idx'):
                    delattr(move_group_dialog, '_drag_atom_idx')
        except Exception:
            pass
        
        # ピックリセットは測定モードで実際に問題が発生した場合のみ行う
        # （通常のドラッグ回転では行わない）
        
        # ボタンを離した後のカーソル表示を最新の状態に更新
        self.on_mouse_move(obj, event)

        # 2Dビューにフォーカスを戻し、ショートカットキーなどが使えるようにする
        if mw and mw.view_2d:
            mw.view_2d.setFocus()

    def on_right_button_up(self, obj, event):
        """
        右クリック終了時の処理。グループ回転を確定。
        """
        mw = self.main_window
        
        # Move Groupダイアログの回転終了処理
        move_group_dialog = None
        try:
            for widget in QApplication.topLevelWidgets():
                if isinstance(widget, MoveGroupDialog) and widget.isVisible():
                    move_group_dialog = widget
                    break
        except Exception:
            pass
        
        if move_group_dialog and getattr(move_group_dialog, '_is_rotating_group_vtk', False):
            # 回転モードで右クリックリリース - 選択を保持
            if getattr(move_group_dialog, '_rotation_mouse_moved', False):
                # 回転が実行された - リリース時に回転を適用
                try:
                    interactor = self.GetInteractor()
                    renderer = mw.plotter.renderer
                    current_pos = interactor.GetEventPosition()
                    conf = mw.current_mol.GetConformer()
                    centroid = move_group_dialog._group_centroid
                    
                    # 掴んだ原子の初期位置
                    if not hasattr(move_group_dialog, '_rotation_atom_idx'):
                        # 最初に掴んだ原子のインデックスを保存
                        move_group_dialog._rotation_atom_idx = next(iter(move_group_dialog.group_atoms))
                    
                    grabbed_atom_idx = move_group_dialog._rotation_atom_idx
                    grabbed_initial_pos = move_group_dialog._initial_positions[grabbed_atom_idx]
                    
                    # 開始位置のスクリーン座標を取得
                    renderer.SetWorldPoint(grabbed_initial_pos[0], grabbed_initial_pos[1], grabbed_initial_pos[2], 1.0)
                    renderer.WorldToDisplay()
                    start_display = renderer.GetDisplayPoint()
                    
                    # 現在のマウス位置をワールド座標に変換（同じ深度で）
                    renderer.SetDisplayPoint(current_pos[0], current_pos[1], start_display[2])
                    renderer.DisplayToWorld()
                    target_world = renderer.GetWorldPoint()
                    target_pos = np.array([target_world[0], target_world[1], target_world[2]])
                    
                    # 重心から見た、掴んだ原子の初期ベクトルと目標ベクトル
                    v1 = grabbed_initial_pos - centroid
                    v2 = target_pos - centroid
                    
                    # ベクトルを正規化
                    v1_norm = np.linalg.norm(v1)
                    v2_norm = np.linalg.norm(v2)
                    
                    if v1_norm > 1e-6 and v2_norm > 1e-6:
                        v1_normalized = v1 / v1_norm
                        v2_normalized = v2 / v2_norm
                        
                        # 回転軸（外積）
                        rotation_axis = np.cross(v1_normalized, v2_normalized)
                        axis_norm = np.linalg.norm(rotation_axis)
                        
                        if axis_norm > 1e-6:
                            rotation_axis = rotation_axis / axis_norm
                            
                            # 回転角（内積）
                            cos_angle = np.clip(np.dot(v1_normalized, v2_normalized), -1.0, 1.0)
                            angle = np.arccos(cos_angle)
                            
                            # Rodriguesの回転公式で回転行列を作成
                            K = np.array([
                                [0, -rotation_axis[2], rotation_axis[1]],
                                [rotation_axis[2], 0, -rotation_axis[0]],
                                [-rotation_axis[1], rotation_axis[0], 0]
                            ])
                            
                            rot_matrix = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
                            
                            # グループ全体を重心周りに回転
                            for atom_idx in move_group_dialog.group_atoms:
                                initial_pos = move_group_dialog._initial_positions[atom_idx]
                                # 重心からの相対座標
                                relative_pos = initial_pos - centroid
                                # 回転を適用
                                rotated_pos = rot_matrix @ relative_pos
                                # 絶対座標に戻す
                                new_pos = rotated_pos + centroid
                                
                                conf.SetAtomPosition(atom_idx, new_pos.tolist())
                                mw.atom_positions_3d[atom_idx] = new_pos
                            
                            # 3D表示を更新
                            mw.draw_molecule_3d(mw.current_mol)
                            mw.update_chiral_labels()
                            move_group_dialog.show_atom_labels()
                            mw.push_undo_state()
                except Exception as e:
                    print(f"Error finalizing group rotation: {e}")
            
            # 状態をリセット
            move_group_dialog._is_rotating_group_vtk = False
            move_group_dialog._rotation_start_pos = None
            move_group_dialog._rotation_mouse_moved = False
            if hasattr(move_group_dialog, '_initial_positions'):
                delattr(move_group_dialog, '_initial_positions')
            if hasattr(move_group_dialog, '_group_centroid'):
                delattr(move_group_dialog, '_group_centroid')
            if hasattr(move_group_dialog, '_rotation_atom_idx'):
                delattr(move_group_dialog, '_rotation_atom_idx')
            
            try:
                mw.plotter.setCursor(Qt.CursorShape.ArrowCursor)
            except Exception:
                pass
            return
        
        # 通常の右クリックリリース処理
        super().OnRightButtonUp()
