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
        # (checkbox, base tooltip) pairs for options CoordGen silently ignores.
        self._coordgen_ignored: list[tuple[Any, str]] = []
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

        canonical_tooltip = (
            "If checked, 'Clean Up 2D Structure' rotates the result into "
            "RDKit's canonical orientation. Note that cleanup always "
            "regenerates coordinates, so unchecking this keeps the "
            "generator's raw orientation, not the pre-cleanup one."
        )
        self.cleanup_canonical_orientation_2d_checkbox = QCheckBox()
        self.cleanup_canonical_orientation_2d_checkbox.setToolTip(canonical_tooltip)
        form_layout.addRow(
            "Canonical Orientation:", self.cleanup_canonical_orientation_2d_checkbox
        )

        ring_templates_tooltip = (
            "If checked, 'Clean Up 2D Structure' uses RDKit's built-in "
            "template shapes for strained polycyclic ring systems. Only a "
            "small set of ring systems has a template, so most structures "
            "are unaffected."
        )
        self.cleanup_use_ring_templates_2d_checkbox = QCheckBox()
        self.cleanup_use_ring_templates_2d_checkbox.setToolTip(ring_templates_tooltip)
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

        avoid_clashes_tooltip = (
            "If checked, 'Clean Up 2D Structure' samples extra layouts to "
            "reduce atom/bond overlap. Mainly affects flexible acyclic "
            "chains, and can be slower on large molecules."
        )
        self.cleanup_avoid_clashes_2d_checkbox = QCheckBox()
        self.cleanup_avoid_clashes_2d_checkbox.setToolTip(avoid_clashes_tooltip)
        form_layout.addRow("Avoid Clashes:", self.cleanup_avoid_clashes_2d_checkbox)

        # RDKit's CoordGen algorithm silently ignores canonOrient,
        # useRingTemplates and the clash-avoidance sampling parameters, so
        # these are grayed out while it is active. StraightenDepiction runs
        # afterwards on the finished depiction, so it still applies.
        self._coordgen_ignored = [
            (self.cleanup_canonical_orientation_2d_checkbox, canonical_tooltip),
            (self.cleanup_use_ring_templates_2d_checkbox, ring_templates_tooltip),
            (self.cleanup_avoid_clashes_2d_checkbox, avoid_clashes_tooltip),
        ]
        self.prefer_coordgen_2d_checkbox.toggled.connect(
            self._update_coordgen_dependent_options
        )

    def _update_coordgen_dependent_options(self, prefer_coordgen: bool) -> None:
        suffix = " (not applicable when Prefer CoordGen is checked)"
        for checkbox, base_tooltip in self._coordgen_ignored:
            checkbox.setEnabled(not prefer_coordgen)
            checkbox.setToolTip(base_tooltip + (suffix if prefer_coordgen else ""))

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
