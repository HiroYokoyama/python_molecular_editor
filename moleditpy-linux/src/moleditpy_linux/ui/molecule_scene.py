#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations
import logging
from typing import Any, List, Optional

from PyQt6.QtCore import QLineF, Qt, QPointF
from PyQt6.QtGui import QPen
from PyQt6.QtWidgets import (
    QApplication,
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsScene,
)

try:
    from .atom_item import AtomItem
    from .bond_item import BondItem
    from .template_preview_item import TemplatePreviewItem
except ImportError:
    from moleditpy_linux.ui.atom_item import AtomItem
    from moleditpy_linux.ui.bond_item import BondItem
    from moleditpy_linux.ui.template_preview_item import TemplatePreviewItem

try:
    from ..utils.constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE
except ImportError:
    pass

try:
    from PyQt6 import sip as _sip  # type: ignore

    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    from ..utils.sip_isdeleted_safe import sip_isdeleted_safe
except ImportError:
    from moleditpy_linux.utils.sip_isdeleted_safe import sip_isdeleted_safe

try:
    from .molecular_scene_handler import TemplateMixin, KeyboardMixin, SceneQueryMixin
except ImportError:
    from moleditpy_linux.ui.molecular_scene_handler import (
        TemplateMixin,
        KeyboardMixin,
        SceneQueryMixin,
    )


class MoleculeScene(TemplateMixin, KeyboardMixin, SceneQueryMixin, QGraphicsScene):
    def __init__(self, data: Any, window: Any) -> None:
        super().__init__()
        self.data, self.window = data, window
        self.mode: str = "select"
        self.current_atom_symbol: str = "C"
        self.bond_order: int = 1
        self.bond_stereo: int = 0
        self.start_atom: Optional[AtomItem] = None
        self.temp_line: Optional[QGraphicsLineItem] = None
        self.start_pos: Optional[QPointF] = None
        self.press_pos: Optional[QPointF] = None
        self.mouse_moved_since_press: bool = False
        self.data_changed_in_event: bool = False
        self.hovered_item: Optional[QGraphicsItem] = None
        # ... (rest of __init__)

        self.key_to_symbol_map = {
            Qt.Key.Key_C: "C",
            Qt.Key.Key_N: "N",
            Qt.Key.Key_O: "O",
            Qt.Key.Key_S: "S",
            Qt.Key.Key_F: "F",
            Qt.Key.Key_B: "B",
            Qt.Key.Key_I: "I",
            Qt.Key.Key_H: "H",
            Qt.Key.Key_P: "P",
        }
        self.key_to_symbol_map_shift = {
            Qt.Key.Key_C: "Cl",
            Qt.Key.Key_B: "Br",
            Qt.Key.Key_S: "Si",
        }

        self.key_to_bond_mode_map = {
            Qt.Key.Key_1: (1, 0),
            Qt.Key.Key_2: (2, 0),
            Qt.Key.Key_3: (3, 0),
            Qt.Key.Key_4: (1, 1),  # Wedge
            Qt.Key.Key_5: (1, 2),  # Dash
        }

        self.reinitialize_items()

    def get_setting(self, key: str, default: Any = None) -> Any:
        """Safe gateway to access MainWindow settings without deep traversal from items."""
        if (
            hasattr(self, "window")
            and self.window
            and hasattr(self.window, "init_manager")
        ):
            return self.window.init_manager.settings.get(key, default)
        return default

    def update_connected_bonds(self, atoms: List[AtomItem]) -> None:
        """Update the positions of all bonds connected to the specified atom list."""
        bonds_to_update = set()
        for atom in atoms:
            if hasattr(atom, "bonds"):
                bonds_to_update.update(atom.bonds)
            else:  # [REPORT ERROR MISSING ATTRIBUTE]
                logging.error("REPORT ERROR: Missing attribute 'bonds' on atom")

        for bond in bonds_to_update:
            if not sip_isdeleted_safe(bond):
                if hasattr(bond, "update_position"):
                    try:
                        bond.update_position()
                    except (RuntimeError, ValueError, TypeError) as e:
                        logging.debug(f"Failed to update bond position for {bond}: {e}")
                else:  # [REPORT ERROR MISSING ATTRIBUTE]
                    logging.error(
                        "REPORT ERROR: Missing attribute 'update_position' on bond"
                    )

    def update_all_items(self) -> None:
        """Force redraw of all items."""
        if hasattr(self.data, "update_ring_info_2d"):
            self.data.update_ring_info_2d()
        else:  # [REPORT ERROR MISSING ATTRIBUTE]
            logging.error(
                "REPORT ERROR: Missing attribute 'update_ring_info_2d' on self.data"
            )

        for item in self.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.update()
        if self.views():
            self.views()[0].viewport().update()

    def reinitialize_items(self) -> None:
        self.template_preview = TemplatePreviewItem()
        self.addItem(self.template_preview)
        self.template_preview.hide()
        self.template_preview_points = []
        self.template_context = {}
        self._deleted_items = []

        app = QApplication.instance()
        if app is not None and hasattr(app, "aboutToQuit"):
            try:
                app.aboutToQuit.connect(self.purge_deleted_items)
            except (RuntimeError, ValueError, TypeError) as e:
                # Non-fatal during setup; app instance may be invalid or signal already connected
                logging.debug(f"Could not connect aboutToQuit in MoleculeScene: {e}")

    def clear_all_problem_flags(self):
        """Reset the has_problem flag for all AtomItems and redraw them."""
        needs_update = False
        for atom_data in self.data.atoms.values():
            item = atom_data.get("item")
            # hasattr is a safety check
            if item and hasattr(item, "has_problem") and item.has_problem:
                item.has_problem = False
                item.update()
                needs_update = True
        return needs_update

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False

        # Record initial positions by safely checking for deleted objects
        self.initial_positions_in_event = {}
        for item in self.items():
            if isinstance(item, AtomItem) and not sip_isdeleted_safe(item):
                try:
                    self.initial_positions_in_event[item] = item.pos()
                except (RuntimeError, AttributeError):
                    # Skip if the object is inaccessible
                    continue

        if not self.window.ui_manager.is_2d_editable:
            return

        if event.button() == Qt.MouseButton.RightButton:
            item = self.itemAt(event.scenePos(), self.views()[0].transform())
            if not isinstance(item, (AtomItem, BondItem)):
                return  # Do nothing if something other than the target is clicked
            data_changed = False
            # If the user has a rectangular multi-selection and the clicked item
            # is part of that selection, delete all selected items (atoms/bonds).
            try:
                # Use getattr safely for selectedItems if scene state is transitioning
                raw_selected = getattr(self, "selectedItems", lambda: [])()
                selected_items = [
                    it
                    for it in raw_selected
                    if isinstance(it, (AtomItem, BondItem))
                    and not sip_isdeleted_safe(it)
                ]
            except Exception as e:
                # Fallback to empty selection if the scene state is inconsistent during event processing
                logging.debug(
                    f"Failed to retrieve selected items in mousePressEvent: {e}"
                )
                selected_items = []

            if (
                len(selected_items) > 1
                and item in selected_items
                and not self.mode.startswith(("template", "charge", "radical"))
            ):
                # Delete the entire rectangular selection
                data_changed = self.delete_items(set(selected_items))
                if data_changed:
                    self.update_all_items()
                    self.window.edit_actions_manager.push_undo_state()
                self.press_pos = None
                event.accept()
                return
            # --- E/Z mode specific processing ---
            if self.mode == "bond_2_5":
                if isinstance(item, BondItem):
                    try:
                        # Clear E/Z label (revert to normal)
                        if hasattr(item, "stereo") and item.stereo in [3, 4]:
                            item.set_stereo(0)
                            # Also update the data model
                            for (id1, id2), bdata in self.data.bonds.items():
                                if bdata.get("item") is item:
                                    bdata["stereo"] = 0
                                    break
                            self.window.edit_actions_manager.push_undo_state()
                            data_changed = False  # Already added to undo stack, so skip redundant pushes later
                    except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                        logging.error(
                            f"Error in E/Z stereo toggle (mousePressEvent): {e}",
                            exc_info=True,
                        )
                        if hasattr(self.window, "statusBar"):
                            sb = self.window.statusBar()
                            if sb:
                                sb.showMessage(f"Error clearing E/Z label: {e}", 5000)
                        else:  # [REPORT ERROR MISSING ATTRIBUTE]
                            logging.error(
                                "REPORT ERROR: Missing attribute 'statusBar' on self.window"
                            )
                        self.update_all_items()  # Redraw even on error to maintain consistency
                # AtomItem does nothing
            # --- Normal processing ---
            elif isinstance(item, AtomItem):
                # Set radical to 0 if in radical mode
                if self.mode == "radical" and item.radical != 0:
                    item.prepareGeometryChange()
                    item.radical = 0
                    self.data.atoms[item.atom_id]["radical"] = 0
                    item.update_style()
                    data_changed = True
                # Set charge to 0 if in charge mode
                elif self.mode in ["charge_plus", "charge_minus"] and item.charge != 0:
                    item.prepareGeometryChange()
                    item.charge = 0
                    self.data.atoms[item.atom_id]["charge"] = 0
                    item.update_style()
                    data_changed = True
                # Delete atom in modes other than those above (excluding template, charge, radical)
                elif not self.mode.startswith(("template", "charge", "radical")):
                    data_changed = self.delete_items({item})
            elif isinstance(item, BondItem):
                # Delete bond in modes other than template, charge, radical
                if not self.mode.startswith(("template", "charge", "radical")):
                    data_changed = self.delete_items({item})

            if data_changed:
                self.update_all_items()
                self.window.edit_actions_manager.push_undo_state()
            self.press_pos = None
            event.accept()
            return  # Complete right-click processing and prevent proceeding to left-click logic

        if self.mode.startswith("template"):
            self.clearSelection()  # In template mode, no selection is made; only the click position is recorded.
            return

        # Prevent selection processing when in Z,E mode
        if self.mode in ["bond_2_5"]:
            self.clearSelection()
            event.accept()
            return

        if getattr(self, "mode", "") != "select":
            self.clearSelection()
            event.accept()

        item = self.itemAt(self.press_pos, self.views()[0].transform())

        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != "select":
                self.clearSelection()
                self.temp_line = QGraphicsLineItem(
                    QLineF(self.start_atom.pos(), self.press_pos)
                )
                self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
                self.addItem(self.temp_line)
            else:
                super().mousePressEvent(event)
        elif item is None and (
            self.mode.startswith("atom") or self.mode.startswith("bond")
        ):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos))
            self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
            self.addItem(self.temp_line)
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.window.ui_manager.is_2d_editable:
            return

        if self.mode.startswith("template"):
            self.update_template_preview(event.scenePos())

        if not self.mouse_moved_since_press and self.press_pos:
            if (
                event.scenePos() - self.press_pos
            ).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True

        if self.temp_line and not self.mode.startswith("template"):
            start_point = self.start_atom.pos() if self.start_atom else self.start_pos
            if not start_point:
                super().mouseMoveEvent(event)
                return

            current_pos = event.scenePos()
            end_point = current_pos

            target_atom = None
            for item in self.items(current_pos):
                if isinstance(item, AtomItem):
                    target_atom = item
                    break

            is_valid_snap_target = target_atom is not None and (
                self.start_atom is None or target_atom is not self.start_atom
            )

            if is_valid_snap_target:
                end_point = target_atom.pos()

            self.temp_line.setLine(QLineF(start_point, end_point))
        else:
            # Even in template mode, hover events are propagated here
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if not self.window.ui_manager.is_2d_editable:
            return

        end_pos = event.scenePos()
        is_click = (
            self.press_pos
            and (end_pos - self.press_pos).manhattanLength()
            < QApplication.startDragDistance()
        )

        if self.temp_line:
            if not sip_isdeleted_safe(self.temp_line):
                try:
                    scene_func = getattr(self.temp_line, "scene", None)
                    if scene_func:
                        sc = scene_func()
                        if sc and hasattr(self, "removeItem"):
                            self.removeItem(self.temp_line)
                except (RuntimeError, ValueError, TypeError, AttributeError) as e:
                    logging.debug(f"Error removing temp_line in mouseReleaseEvent: {e}")

            self.temp_line = None

        if self.mode.startswith("template") and is_click:
            if self.template_context and self.template_context.get("points"):
                context = self.template_context
                # Check if this is a user template
                if self.mode.startswith("template_user"):
                    self.add_user_template_fragment(context)
                else:
                    self.add_molecule_fragment(
                        context["points"],
                        context["bonds_info"],
                        existing_items=context.get("items", []),
                    )
                self.data_changed_in_event = True
                # Complete event processing here to prevent selection of underlying items
                self.start_atom = None
                self.press_pos = None
                if self.data_changed_in_event:
                    self.update_all_items()
                    self.window.edit_actions_manager.push_undo_state()
                return

        released_item = self.itemAt(end_pos, self.views()[0].transform())

        # 1. Handle special modes (delete/radical/charge)
        if (self.mode == "delete") and is_click and released_item is not None:
            # Safe deletion via unified handler
            if self.delete_items({released_item}):
                self.window.edit_actions_manager.push_undo_state()
            self.press_pos = None
            return

        elif (
            (self.mode == "radical")
            and is_click
            and isinstance(released_item, AtomItem)
        ):
            atom = released_item
            atom.prepareGeometryChange()
            # Toggle radical state (0 -> 1 -> 2 -> 0)
            atom.radical = (atom.radical + 1) % 3
            if atom.atom_id in self.data.atoms:
                self.data.atoms[atom.atom_id]["radical"] = atom.radical
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom = None
            self.start_pos = None
            self.press_pos = None
            if self.data_changed_in_event:
                self.window.edit_actions_manager.push_undo_state()
            return
        elif (
            (self.mode == "charge_plus" or self.mode == "charge_minus")
            and is_click
            and isinstance(released_item, AtomItem)
        ):
            atom = released_item
            atom.prepareGeometryChange()
            delta = 1 if self.mode == "charge_plus" else -1
            atom.charge += delta
            if atom.atom_id in self.data.atoms:
                self.data.atoms[atom.atom_id]["charge"] = atom.charge
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom = None
            self.start_pos = None
            self.press_pos = None
            if self.data_changed_in_event:
                self.window.edit_actions_manager.push_undo_state()
            return

        elif (
            self.mode.startswith("bond")
            and is_click
            and isinstance(released_item, BondItem)
        ):
            b = released_item
            if self.mode == "bond_2_5":
                try:
                    if hasattr(b, "order") and b.order == 2:
                        current_stereo = getattr(b, "stereo", 0)
                        if current_stereo not in [3, 4]:
                            new_stereo = 3  # None -> Z
                        elif current_stereo == 3:
                            new_stereo = 4  # Z -> E
                        else:  # current_stereo == 4
                            new_stereo = 0  # E -> None

                        if hasattr(self, "update_bond_stereo"):
                            self.update_bond_stereo(b, new_stereo)
                            self.update_all_items()  # Force redraw
                            self.window.edit_actions_manager.push_undo_state()  # Push to undo stack here
                        else:  # [REPORT ERROR MISSING ATTRIBUTE]
                            logging.error(
                                "REPORT ERROR: Missing attribute 'update_bond_stereo' on self"
                            )
                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    logging.error(
                        f"Error in E/Z stereo toggle (mouseReleaseEvent): {e}",
                        exc_info=True,
                    )
                    if hasattr(self.window, "statusBar"):
                        sb = self.window.statusBar()
                        if sb:
                            sb.showMessage(
                                f"Error changing E/Z stereochemistry: {e}", 5000
                            )
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error(
                            "REPORT ERROR: Missing attribute 'statusBar' on self.window"
                        )
                    self.update_all_items()  # Redraw even on error to maintain consistency
                return  # Do not proceed further
            elif (
                self.bond_stereo != 0
                and b.order == self.bond_order
                and b.stereo == self.bond_stereo
            ):
                # Invert bond direction
                old_id1, old_id2 = b.atom1.atom_id, b.atom2.atom_id
                # 1. Remove old bond from data
                self.data.remove_bond(old_id1, old_id2)
                # 2. Add bond in reverse direction
                new_key, _ = self.data.add_bond(
                    old_id2, old_id1, self.bond_order, self.bond_stereo
                )
                # 3. Swap atom references and link to new data in BondItem
                b.atom1, b.atom2 = b.atom2, b.atom1
                self.data.bonds[new_key]["item"] = b
                # 4. Update visual state
                b.update_position()
            else:
                # Remove existing bond once
                self.data.remove_bond(b.atom1.atom_id, b.atom2.atom_id)
                # Recreate bond in the direction recorded by BondItem (b.atom1 -> b.atom2)
                # This ensures the corrected add_bond is called and saved with proper directionality.
                new_key, _ = self.data.add_bond(
                    b.atom1.atom_id, b.atom2.atom_id, self.bond_order, self.bond_stereo
                )
                # Update BondItem visual state and data reference
                b.prepareGeometryChange()
                b.order = self.bond_order
                b.stereo = self.bond_stereo
                self.data.bonds[new_key]["item"] = b
                b.update()
            self.clearSelection()
            self.data_changed_in_event = True
        # 3. Create new atom/bond (Allowed in atom_* and all bond_* modes)
        elif self.start_atom and (
            self.mode.startswith("atom") or self.mode.startswith("bond")
        ):
            line = QLineF(self.start_atom.pos(), end_pos)
            end_item = self.itemAt(end_pos, self.views()[0].transform())
            # Determine bond style to use
            # In atom modes, set bond_order/stereo to None so create_bond uses defaults (1, 0)
            # In bond_* modes, use current settings (self.bond_order/stereo)
            order_to_use = self.bond_order if self.mode.startswith("bond") else None
            stereo_to_use = self.bond_stereo if self.mode.startswith("bond") else None
            if is_click:
                # Short click: Update existing atom symbol (atom mode only)
                if (
                    self.mode.startswith("atom")
                    and self.start_atom.symbol != self.current_atom_symbol
                ):
                    self.start_atom.symbol = self.current_atom_symbol
                    self.data.atoms[self.start_atom.atom_id]["symbol"] = (
                        self.current_atom_symbol
                    )
                    self.start_atom.update_style()
                    self.data_changed_in_event = True
            else:
                # Drag: new bond or bond to existing atom
                if isinstance(end_item, AtomItem) and self.start_atom != end_item:
                    self.create_bond(
                        self.start_atom,
                        end_item,
                        bond_order=order_to_use,
                        bond_stereo=stereo_to_use,
                    )
                else:
                    new_id = self.create_atom(self.current_atom_symbol, end_pos)
                    new_item = self.data.atoms[new_id]["item"]
                    self.create_bond(
                        self.start_atom,
                        new_item,
                        bond_order=order_to_use,
                        bond_stereo=stereo_to_use,
                    )
                self.data_changed_in_event = True
        # 4. Create new from empty area (allowed in atom_* and all bond_* modes)
        elif self.start_pos and (
            self.mode.startswith("atom") or self.mode.startswith("bond")
        ):
            line = QLineF(self.start_pos, end_pos)
            # Determine bond type to use
            order_to_use = self.bond_order if self.mode.startswith("bond") else None
            stereo_to_use = self.bond_stereo if self.mode.startswith("bond") else None
            if line.length() < 10:
                self.create_atom(self.current_atom_symbol, end_pos)
                self.data_changed_in_event = True
            else:
                end_item = self.itemAt(end_pos, self.views()[0].transform())
                if isinstance(end_item, AtomItem):
                    start_id = self.create_atom(
                        self.current_atom_symbol, self.start_pos
                    )
                    start_item = self.data.atoms[start_id]["item"]
                    self.create_bond(
                        start_item,
                        end_item,
                        bond_order=order_to_use,
                        bond_stereo=stereo_to_use,
                    )
                else:
                    start_id = self.create_atom(
                        self.current_atom_symbol, self.start_pos
                    )
                    end_id = self.create_atom(self.current_atom_symbol, end_pos)
                    self.create_bond(
                        self.data.atoms[start_id]["item"],
                        self.data.atoms[end_id]["item"],
                        bond_order=order_to_use,
                        bond_stereo=stereo_to_use,
                    )
                self.data_changed_in_event = True
        # 5. Other processing (Select mode, etc.)
        else:
            super().mouseReleaseEvent(event)

        # Safely check for moved objects
        moved_atoms = []
        initial_positions = getattr(self, "initial_positions_in_event", {})
        for item, old_pos in initial_positions.items():
            try:
                # Check if object is valid, in scene, and position changed
                if item.scene() and item.pos() != old_pos:
                    moved_atoms.append(item)
            except RuntimeError:
                # Skip if object is deleted
                continue
        if moved_atoms:
            self.data_changed_in_event = True
            bonds_to_update = set()
            for atom in moved_atoms:
                try:
                    self.data.set_atom_pos(atom.atom_id, atom.pos())
                    bonds_to_update.update(atom.bonds)
                except RuntimeError:
                    # Skip if object is deleted
                    continue
            for bond in bonds_to_update:
                bond.update_position()
            # Update measurement label positions after atom move
            self.window.edit_3d_manager.update_2d_measurement_labels()
            if self.views():
                self.views()[0].viewport().update()

        if getattr(self, "data_changed_in_event", False):
            self.update_all_items()

        self.start_atom = None
        self.start_pos = None
        self.press_pos = None
        self.temp_line = None
        # Clear template context but NOT the template data itself to allow multiple placements
        self.template_context = {}
        if getattr(self, "data_changed_in_event", False):
            self.window.edit_actions_manager.push_undo_state()

    def mouseDoubleClickEvent(self, event):
        """Handle double click events."""
        item = self.itemAt(event.scenePos(), self.views()[0].transform())

        if self.mode in ["charge_plus", "charge_minus", "radical"] and isinstance(
            item, AtomItem
        ):
            if self.mode == "radical":
                item.prepareGeometryChange()
                item.radical = (item.radical + 1) % 3
                self.data.atoms[item.atom_id]["radical"] = item.radical
                item.update_style()
            else:
                item.prepareGeometryChange()
                delta = 1 if self.mode == "charge_plus" else -1
                item.charge += delta
                self.data.atoms[item.atom_id]["charge"] = item.charge
                item.update_style()

            self.update_all_items()
            self.window.edit_actions_manager.push_undo_state()

            event.accept()
            return

        # Select-mode: double-click should select the clicked atom/bond and
        # only the atoms/bonds connected to it (the connected component).
        if self.mode == "select" and isinstance(item, (AtomItem, BondItem)):
            try:
                start_atoms = set()
                if isinstance(item, AtomItem):
                    start_atoms.add(item)
                else:
                    # BondItem: start from both ends if available
                    a1 = getattr(item, "atom1", None)
                    a2 = getattr(item, "atom2", None)
                    if a1 is not None:
                        start_atoms.add(a1)
                    if a2 is not None:
                        start_atoms.add(a2)

                # BFS/DFS over atoms via bond references (defensive checks)
                atoms_to_visit = list(start_atoms)
                connected_atoms = set()
                connected_bonds = set()

                while atoms_to_visit:
                    a = atoms_to_visit.pop()
                    if a is None or sip_isdeleted_safe(a):
                        continue
                    if a in connected_atoms:
                        continue
                    connected_atoms.add(a)
                    # iterate bonds attached to atom
                    atom_bonds = getattr(a, "bonds", [])
                    if atom_bonds:
                        for b in atom_bonds:
                            if b is None or sip_isdeleted_safe(b):
                                continue
                            connected_bonds.add(b)
                            # find the other atom at the bond
                            other = None
                            try:
                                b_atom1 = getattr(b, "atom1", None)
                                b_atom2 = getattr(b, "atom2", None)
                                if b_atom1 is a:
                                    other = b_atom2
                                else:
                                    other = b_atom1
                            except (
                                AttributeError,
                                RuntimeError,
                                ValueError,
                                TypeError,
                            ):
                                other = None
                            if (
                                other is not None
                                and not sip_isdeleted_safe(other)
                                and other not in connected_atoms
                            ):
                                atoms_to_visit.append(other)

                # Apply selection: clear previous and select only these
                self.clearSelection()

                for a in connected_atoms:
                    if not sip_isdeleted_safe(a) and hasattr(a, "setSelected"):
                        try:
                            a.setSelected(True)
                        except (
                            RuntimeError,
                            ValueError,
                            TypeError,
                            AttributeError,
                        ) as e:
                            logging.debug(f"Failed to select atom {a}: {e}")

                for b in connected_bonds:
                    if not sip_isdeleted_safe(b) and hasattr(b, "setSelected"):
                        try:
                            b.setSelected(True)
                        except (
                            RuntimeError,
                            ValueError,
                            TypeError,
                            AttributeError,
                        ) as e:
                            logging.debug(f"Failed to select bond {b}: {e}")
                event.accept()
                return
            except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                # On any unexpected error, fall back to default handling
                logging.error(f"Error in connected component selection: {e}")

        elif self.mode in ["bond_2_5"]:
            event.accept()
            return

        super().mouseDoubleClickEvent(event)

    def purge_deleted_items(self):
        """Purge and release any held deleted-wrapper references during shutdown."""
        if not getattr(self, "_deleted_items", None):
            return

        for obj in list(self._deleted_items):
            if not sip_isdeleted_safe(obj):
                try:
                    if hasattr(obj, "hide"):
                        obj.hide()
                    else:  # [REPORT ERROR MISSING ATTRIBUTE]
                        logging.error("REPORT ERROR: Missing attribute 'hide' on obj")
                    if hasattr(obj, "bonds") and obj.bonds is not None:
                        if hasattr(obj.bonds, "clear"):
                            obj.bonds.clear()
                        else:  # [REPORT ERROR MISSING ATTRIBUTE]
                            logging.error(
                                "REPORT ERROR: Missing attribute 'clear' on object"
                            )
                except (AttributeError, RuntimeError, ValueError, TypeError) as e:
                    logging.debug(f"Error purging item {obj} in MoleculeScene: {e}")

        try:
            self._deleted_items.clear()
        except (AttributeError, RuntimeError, ValueError, TypeError) as e:
            logging.debug(f"Error clearing _deleted_items list: {e}")
            self._deleted_items = []

    def leaveEvent(self, event):
        self.template_preview.hide()

    def refresh_mode_state(self):
        """Immediately update scene state and previews based on the current mouse position."""
        import PyQt6.QtGui

        global_pos = PyQt6.QtGui.QCursor.pos()

        # Find the active view for this scene
        for view in self.views():
            if view.isVisible():
                # Map global cursor position to scene coordinates
                local_pos = view.mapFromGlobal(global_pos)
                scene_pos = view.mapToScene(local_pos)

                # If the mouse is within the viewport, trigger the preview update
                if view.viewport().rect().contains(local_pos):
                    if hasattr(
                        self, "update_template_preview"
                    ) and self.mode.startswith("template"):
                        self.update_template_preview(scene_pos)
                    return

    def set_hovered_item(self, item):
        """Record currently hovered item"""
        self.hovered_item = item
