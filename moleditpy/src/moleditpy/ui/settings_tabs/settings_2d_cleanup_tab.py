#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from collections.abc import Mapping
from typing import Any, Optional

from PyQt6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QWidget,
)
from .settings_tab_base import SettingsTabBase


class Settings2DCleanupTab(SettingsTabBase):
    """Settings tab for the 'Clean Up 2D Structure' layout options."""

    def __init__(
        self, default_settings: Mapping[str, Any], parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(default_settings, parent)
        self.prefer_coordgen_2d_checkbox: Any = None
        self.cleanup_canonical_orientation_2d_checkbox: Any = None
        self.cleanup_use_ring_templates_2d_checkbox: Any = None
        self.cleanup_straighten_bonds_2d_checkbox: Any = None
        self.cleanup_avoid_clashes_2d_checkbox: Any = None
        self._setup_ui()

    def _setup_ui(self) -> None:
        form_layout = QFormLayout(self)

        self.prefer_coordgen_2d_checkbox = QCheckBox()
        self.prefer_coordgen_2d_checkbox.setToolTip(
            "If checked, 'Clean Up 2D Structure' uses RDKit's CoordGen "
            "algorithm, which typically produces nicer ring layouts than "
            "the default generator."
        )
        form_layout.addRow(
            "Prefer CoordGen for Cleanup:", self.prefer_coordgen_2d_checkbox
        )

        self.cleanup_canonical_orientation_2d_checkbox = QCheckBox()
        self.cleanup_canonical_orientation_2d_checkbox.setToolTip(
            "If checked, 'Clean Up 2D Structure' rotates the result into "
            "RDKit's canonical orientation."
        )
        form_layout.addRow(
            "Canonical Orientation:", self.cleanup_canonical_orientation_2d_checkbox
        )

        self.cleanup_use_ring_templates_2d_checkbox = QCheckBox()
        self._ring_templates_tooltip = (
            "If checked, 'Clean Up 2D Structure' uses template shapes for "
            "complex fused ring systems."
        )
        self.cleanup_use_ring_templates_2d_checkbox.setToolTip(
            self._ring_templates_tooltip
        )
        form_layout.addRow(
            "Use Ring Templates:", self.cleanup_use_ring_templates_2d_checkbox
        )

        self.cleanup_straighten_bonds_2d_checkbox = QCheckBox()
        self.cleanup_straighten_bonds_2d_checkbox.setToolTip(
            "If checked, 'Clean Up 2D Structure' rotates the result so most "
            "bonds align to 30°/90° angles."
        )
        form_layout.addRow(
            "Straighten Bonds:", self.cleanup_straighten_bonds_2d_checkbox
        )

        self.cleanup_avoid_clashes_2d_checkbox = QCheckBox()
        self._avoid_clashes_tooltip = (
            "If checked, 'Clean Up 2D Structure' samples extra layouts to "
            "reduce atom/bond overlap. Can be slower on large molecules."
        )
        self.cleanup_avoid_clashes_2d_checkbox.setToolTip(self._avoid_clashes_tooltip)
        form_layout.addRow(
            "Avoid Clashes:", self.cleanup_avoid_clashes_2d_checkbox
        )

        # Ring templates and clash-avoidance sampling are both ignored by
        # RDKit's CoordGen algorithm, so they don't apply when it's active.
        self.prefer_coordgen_2d_checkbox.toggled.connect(
            self._update_coordgen_dependent_options
        )

    def _update_coordgen_dependent_options(self, prefer_coordgen: bool) -> None:
        self.cleanup_use_ring_templates_2d_checkbox.setEnabled(not prefer_coordgen)
        self.cleanup_avoid_clashes_2d_checkbox.setEnabled(not prefer_coordgen)
        not_applicable = " (not applicable when Prefer CoordGen is checked)"
        self.cleanup_use_ring_templates_2d_checkbox.setToolTip(
            self._ring_templates_tooltip
            + (not_applicable if prefer_coordgen else "")
        )
        self.cleanup_avoid_clashes_2d_checkbox.setToolTip(
            self._avoid_clashes_tooltip + (not_applicable if prefer_coordgen else "")
        )

    def update_ui(self, settings_dict: Mapping[str, Any]) -> None:
        prefer_coordgen = settings_dict.get("prefer_coordgen_2d", False)
        self.prefer_coordgen_2d_checkbox.setChecked(prefer_coordgen)
        self._update_coordgen_dependent_options(prefer_coordgen)
        self.cleanup_canonical_orientation_2d_checkbox.setChecked(
            settings_dict.get("cleanup_canonical_orientation_2d", True)
        )
        self.cleanup_use_ring_templates_2d_checkbox.setChecked(
            settings_dict.get("cleanup_use_ring_templates_2d", False)
        )
        self.cleanup_straighten_bonds_2d_checkbox.setChecked(
            settings_dict.get("cleanup_straighten_bonds_2d", False)
        )
        self.cleanup_avoid_clashes_2d_checkbox.setChecked(
            settings_dict.get("cleanup_avoid_clashes_2d", False)
        )

    def get_settings(self) -> dict[str, Any]:
        return {
            "prefer_coordgen_2d": self.prefer_coordgen_2d_checkbox.isChecked(),
            "cleanup_canonical_orientation_2d": (
                self.cleanup_canonical_orientation_2d_checkbox.isChecked()
            ),
            "cleanup_use_ring_templates_2d": (
                self.cleanup_use_ring_templates_2d_checkbox.isChecked()
            ),
            "cleanup_straighten_bonds_2d": (
                self.cleanup_straighten_bonds_2d_checkbox.isChecked()
            ),
            "cleanup_avoid_clashes_2d": (
                self.cleanup_avoid_clashes_2d_checkbox.isChecked()
            ),
        }
