#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import QWidget

class SettingsTabBase(QWidget):
    """Base class for all settings tabs."""
    def __init__(self, default_settings, parent=None):
        super().__init__(parent)
        self.default_settings = default_settings

    def update_ui(self, settings_dict):
        """Update UI based on settings dictionary. Must be implemented by subclass."""
        raise NotImplementedError

    def get_settings(self):
        """Get settings values from the current UI. Must be implemented by subclass."""
        raise NotImplementedError

    def reset_to_defaults(self):
        """Reset only the settings of the current tab to defaults."""
        self.update_ui(self.default_settings)
