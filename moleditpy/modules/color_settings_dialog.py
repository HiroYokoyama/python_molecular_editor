from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QGridLayout, QPushButton, QHBoxLayout, QLabel,
    QApplication, QColorDialog
)
from PyQt6.QtGui import QColor
try:
    from .constants import CPK_COLORS, DEFAULT_CPK_COLORS
except Exception:
    from modules.constants import CPK_COLORS, DEFAULT_CPK_COLORS

class ColorSettingsDialog(QDialog):
    """Dialog to customize CPK element colors.

    - Click an element to pick a new color for the element (CPK colors).
    - Reset All button to restore defaults for everything.
    """
    def __init__(self, current_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CPK Colors")
        self.parent_window = parent
        self.current_settings = current_settings or {}

        self.changed_cpk = {}  # symbol -> hex
        self._reset_all_flag = False

        layout = QVBoxLayout(self)

        # Color picking for CPK is available in the periodic table and CPK dialog

        # Periodic table grid (buttons like PeriodicTableDialog)
        grid = QGridLayout()
        self.element_buttons = {}
        elements = [
            ('H',1,1), ('He',1,18),
            ('Li',2,1), ('Be',2,2), ('B',2,13), ('C',2,14), ('N',2,15), ('O',2,16), ('F',2,17), ('Ne',2,18),
            ('Na',3,1), ('Mg',3,2), ('Al',3,13), ('Si',3,14), ('P',3,15), ('S',3,16), ('Cl',3,17), ('Ar',3,18),
            ('K',4,1), ('Ca',4,2), ('Sc',4,3), ('Ti',4,4), ('V',4,5), ('Cr',4,6), ('Mn',4,7), ('Fe',4,8),
            ('Co',4,9), ('Ni',4,10), ('Cu',4,11), ('Zn',4,12), ('Ga',4,13), ('Ge',4,14), ('As',4,15), ('Se',4,16),
            ('Br',4,17), ('Kr',4,18),
            ('Rb',5,1), ('Sr',5,2), ('Y',5,3), ('Zr',5,4), ('Nb',5,5), ('Mo',5,6), ('Tc',5,7), ('Ru',5,8),
            ('Rh',5,9), ('Pd',5,10), ('Ag',5,11), ('Cd',5,12), ('In',5,13), ('Sn',5,14), ('Sb',5,15), ('Te',5,16),
            ('I',5,17), ('Xe',5,18),
            ('Cs',6,1), ('Ba',6,2), ('Hf',6,4), ('Ta',6,5), ('W',6,6), ('Re',6,7), ('Os',6,8),
            ('Ir',6,9), ('Pt',6,10), ('Au',6,11), ('Hg',6,12), ('Tl',6,13), ('Pb',6,14), ('Bi',6,15), ('Po',6,16),
            ('At',6,17), ('Rn',6,18),
            ('Fr',7,1), ('Ra',7,2), ('Rf',7,4), ('Db',7,5), ('Sg',7,6), ('Bh',7,7), ('Hs',7,8),
            ('Mt',7,9), ('Ds',7,10), ('Rg',7,11), ('Cn',7,12), ('Nh',7,13), ('Fl',7,14), ('Mc',7,15), ('Lv',7,16),
            ('Ts',7,17), ('Og',7,18),
            ('La',8,3), ('Ce',8,4), ('Pr',8,5), ('Nd',8,6), ('Pm',8,7), ('Sm',8,8), ('Eu',8,9), ('Gd',8,10), ('Tb',8,11),
            ('Dy',8,12), ('Ho',8,13), ('Er',8,14), ('Tm',8,15), ('Yb',8,16), ('Lu',8,17),
            ('Ac',9,3), ('Th',9,4), ('Pa',9,5), ('U',9,6), ('Np',9,7), ('Pu',9,8), ('Am',9,9), ('Cm',9,10), ('Bk',9,11),
            ('Cf',9,12), ('Es',9,13), ('Fm',9,14), ('Md',9,15), ('No',9,16), ('Lr',9,17),
        ]

        for symbol, row, col in elements:
            b = QPushButton(symbol)
            b.setFixedSize(40, 40)
            # Choose override color (if present) else default CPK color
            override = self.current_settings.get('cpk_colors', {}).get(symbol)
            if override:
                q_color = QColor(override)
            else:
                q_color = CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT'])

            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            b.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
            b.clicked.connect(self.on_element_clicked)
            grid.addWidget(b, row, col)
            self.element_buttons[symbol] = b

        layout.addLayout(grid)

        # Ball & Stick bond color (3D) picker - placed near the periodic table for CPK settings
        self.changed_bs_color = None
        try:
            bs_h = QHBoxLayout()
            bs_label = QLabel("Ball & Stick bond color:")
            self.bs_button = QPushButton()
            self.bs_button.setFixedSize(36, 24)
            # initialize from current settings (if provided)
            try:
                cur_bs = self.current_settings.get('ball_stick_bond_color') if self.current_settings else None
            except Exception:
                cur_bs = None
            if not cur_bs and self.parent_window and hasattr(self.parent_window, 'settings'):
                cur_bs = self.parent_window.settings.get('ball_stick_bond_color', '#7F7F7F')
            try:
                self.bs_button.setStyleSheet(f"background-color: {cur_bs}; border: 1px solid #888;")
                self.bs_button.setToolTip(cur_bs)
            except Exception:
                pass
            self.bs_button.clicked.connect(self.pick_bs_bond_color)
            bs_h.addWidget(bs_label)
            bs_h.addWidget(self.bs_button)
            bs_h.addStretch(1)
            layout.addLayout(bs_h)
        except Exception:
            pass

        # Reset button and action buttons
        h = QHBoxLayout()
        reset_button = QPushButton("Reset All")
        reset_button.clicked.connect(self.reset_all)
        h.addWidget(reset_button)
        h.addStretch(1)
        apply_button = QPushButton("Apply")
        apply_button.clicked.connect(self.apply_changes)
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)

        h.addWidget(apply_button); h.addWidget(ok_button); h.addWidget(cancel_button)
        layout.addLayout(h)

        # initialize
        # No 2D bond color control here

    # 2D bond color picker removed — 2D bond color is fixed and not configurable here

    def on_element_clicked(self):
        b = self.sender()
        symbol = b.text()
        # get current color (override if exists else default)
        cur = self.current_settings.get('cpk_colors', {}).get(symbol)
        if not cur:
            cur = CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT']).name()
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            self.changed_cpk[symbol] = color.name()
            # Update button appearance
            brightness = (color.red() * 299 + color.green() * 587 + color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            b.setStyleSheet(f"background-color: {color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")

    def reset_all(self):
        # Clear overrides
        self.changed_cpk = {}
        self._reset_all_flag = True
        
        # 1. B&S結合色もリセット対象（デフォルト値）に設定
        try:
            self.changed_bs_color = self.parent_window.default_settings.get('ball_stick_bond_color', '#7F7F7F') if hasattr(self.parent_window, 'default_settings') else '#7F7F7F'
        except Exception:
            self.changed_bs_color = '#7F7F7F'
            
        # 2. ダイアログ内のCPKボタンの表示をデフォルトに戻す
        for s, btn in self.element_buttons.items():
            q_color = DEFAULT_CPK_COLORS.get(s, DEFAULT_CPK_COLORS['DEFAULT'])
            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = 'white' if brightness < 128 else 'black'
            btn.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
        
        # 3. 3Dプレビューを更新する L.3337〜L.3386 の try...finally ブロックは削除

        # 4. ダイアログ内のB&S結合色ボタンの表示をデフォルトに戻す
        try:
            if hasattr(self, 'bs_button'):
                # self.changed_bs_color に設定したデフォルト値を反映
                hexv = self.changed_bs_color
                self.bs_button.setStyleSheet(f"background-color: {hexv}; border: 1px solid #888;")
                self.bs_button.setToolTip(hexv)
        except Exception:
            pass

    def apply_changes(self):
        # Persist only changed keys
        if self.parent_window:
            if self._reset_all_flag:
                # Remove any cpk overrides
                try:
                    if 'cpk_colors' in self.parent_window.settings:
                        del self.parent_window.settings['cpk_colors']
                except Exception:
                    pass
            if self.changed_cpk:
                # Merge with existing overrides
                cdict = self.parent_window.settings.get('cpk_colors', {}).copy()
                cdict.update(self.changed_cpk)
                self.parent_window.settings['cpk_colors'] = cdict
                self.parent_window.settings_dirty = True
            # After changing settings, update global CPK color map and refresh views
            try:
                self.parent_window.update_cpk_colors_from_settings()
            except Exception:
                pass
            try:
                self.parent_window.apply_3d_settings(redraw=False)
            except Exception:
                pass
            try:
                if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                    self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
            except Exception:
                pass
            # update 2D scene objects
            try:
                if hasattr(self.parent_window, 'scene'):
                    for it in self.parent_window.scene.items():
                        try:
                            if hasattr(it, 'update_style'):
                                it.update_style()
                        except Exception:
                            pass
            except Exception:
                pass
            # update periodic table button styles in the dialog to reflect any overrides
            try:
                for s, btn in self.element_buttons.items():
                    try:
                        q_color = QColor(self.parent_window.settings.get('cpk_colors', {}).get(s, CPK_COLORS.get(s, CPK_COLORS['DEFAULT']).name()))
                        brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
                        text_color = 'white' if brightness < 128 else 'black'
                        btn.setStyleSheet(f"background-color: {q_color.name()}; color: {text_color}; border: 1px solid #555; font-weight: bold;")
                    except Exception:
                        pass
            except Exception:
                pass
            # Refresh any open SettingsDialog instances so the ball & stick color preview updates
            try:
                # Avoid circular import at module level; import SettingsDialog on demand
                try:
                    from .settings_dialog import SettingsDialog
                except Exception:
                    from modules.settings_dialog import SettingsDialog

                for w in QApplication.topLevelWidgets():
                    try:
                        if isinstance(w, SettingsDialog):
                            try:
                                w.update_ui_from_settings(self.parent_window.settings)
                            except Exception:
                                pass
                    except Exception:
                        pass
            except Exception:
                pass
            # Persist changed Ball & Stick color if the user changed it from the CPK dialog
            if getattr(self, 'changed_bs_color', None):
                try:
                    self.parent_window.settings['ball_stick_bond_color'] = self.changed_bs_color
                    try:
                        self.parent_window.settings_dirty = True
                    except Exception:
                        pass
                    # After changing ball-stick color, ensure 3D view updates
                    try:
                        self.parent_window.apply_3d_settings()
                    except Exception:
                        pass
                    try:
                        if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                            self.parent_window.draw_molecule_3d(self.parent_window.current_mol)
                    except Exception:
                        pass
                except Exception:
                    pass
            # Removed 2D bond color control — nothing to persist here
            elif self._reset_all_flag:
                # Reset Ball & Stick 3D bond color to default
                try:
                    # Use a stable default instead of relying on parent_window.default_settings
                    self.parent_window.settings['ball_stick_bond_color'] = '#7F7F7F'
                    try:
                        self.parent_window.settings_dirty = True
                    except Exception:
                        pass
                except Exception:
                    pass
                self.parent_window.update_cpk_colors_from_settings()
                self.parent_window.apply_3d_settings()
                if hasattr(self.parent_window, 'current_mol') and self.parent_window.current_mol:
                    self.parent_window.draw_molecule_3d(self.parent_window.current_mol)

            # update 2D scene
            try:
                if hasattr(self.parent_window, 'scene'):
                    for it in self.parent_window.scene.items():
                        try:
                            # AtomItem.update_style uses CPK_COLORS map
                            if hasattr(it, 'update_style'):
                                it.update_style()
                        except Exception:
                            pass
            except Exception:
                pass

    def accept(self):
        self.apply_changes()
        super().accept()

    def pick_bs_bond_color(self):
        """Pick Ball & Stick 3D bond color from the CPK dialog and update preview immediately."""
        try:
            cur = getattr(self, 'changed_bs_color', None) or (self.current_settings.get('ball_stick_bond_color') if self.current_settings else None)
        except Exception:
            cur = None
        if not cur and self.parent_window and hasattr(self.parent_window, 'settings'):
            cur = self.parent_window.settings.get('ball_stick_bond_color', '#7F7F7F')
        color = QColorDialog.getColor(QColor(cur), self)
        if color.isValid():
            hexv = color.name()
            self.changed_bs_color = hexv
            try:
                self.bs_button.setStyleSheet(f"background-color: {hexv}; border: 1px solid #888;")
                self.bs_button.setToolTip(hexv)
            except Exception:
                pass
