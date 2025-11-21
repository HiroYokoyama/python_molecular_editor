#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import QDialog, QGridLayout, QPushButton
from PyQt6.QtGui import QColor
try:
    from .constants import CPK_COLORS
except Exception:
    from modules.constants import CPK_COLORS

from PyQt6.QtCore import pyqtSignal

class PeriodicTableDialog(QDialog):
    element_selected = pyqtSignal(str)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select an Element")
        layout = QGridLayout(self)
        self.setLayout(layout)

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
            # Lanthanides (placed on a separate row)
            ('La',8,3), ('Ce',8,4), ('Pr',8,5), ('Nd',8,6), ('Pm',8,7), ('Sm',8,8), ('Eu',8,9), ('Gd',8,10), ('Tb',8,11),
            ('Dy',8,12), ('Ho',8,13), ('Er',8,14), ('Tm',8,15), ('Yb',8,16), ('Lu',8,17),
            # Actinides (separate row)
            ('Ac',9,3), ('Th',9,4), ('Pa',9,5), ('U',9,6), ('Np',9,7), ('Pu',9,8), ('Am',9,9), ('Cm',9,10), ('Bk',9,11),
            ('Cf',9,12), ('Es',9,13), ('Fm',9,14), ('Md',9,15), ('No',9,16), ('Lr',9,17),
        ]
        for symbol, row, col in elements:
            b = QPushButton(symbol)
            b.setFixedSize(40,40)

            # Prefer saved user override (from parent.settings), otherwise use CPK_COLORS
            try:
                overrides = parent.settings.get('cpk_colors', {}) if parent and hasattr(parent, 'settings') else {}
                override = overrides.get(symbol)
            except Exception:
                override = None
            q_color = QColor(override) if override else CPK_COLORS.get(symbol, CPK_COLORS['DEFAULT'])

            # 背景色の輝度を計算して、文字色を黒か白に決定
            # 輝度 = (R*299 + G*587 + B*114) / 1000
            brightness = (q_color.red() * 299 + q_color.green() * 587 + q_color.blue() * 114) / 1000
            text_color = "white" if brightness < 128 else "black"

            # ボタンのスタイルシートを設定
            b.setStyleSheet(
                f"background-color: {q_color.name()};"
                f"color: {text_color};"
                "border: 1px solid #555;"
                "font-weight: bold;"
            )

            b.clicked.connect(self.on_button_clicked)
            layout.addWidget(b, row, col)

    def on_button_clicked(self):
        b=self.sender()
        self.element_selected.emit(b.text())
        self.accept()
