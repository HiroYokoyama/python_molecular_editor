#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

import logging
import math
from PyQt6.QtCore import QLineF, QPointF, QRectF, Qt
from PyQt6.QtGui import QCursor, QPen
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
    from modules.atom_item import AtomItem
    from modules.bond_item import BondItem
    from modules.template_preview_item import TemplatePreviewItem

try:
    from .constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE
except ImportError:
    from modules.constants import DEFAULT_BOND_LENGTH, SNAP_DISTANCE, SUM_TOLERANCE

try:
    from PyQt6 import sip as _sip  # type: ignore
    _sip_isdeleted = getattr(_sip, "isdeleted", None)
except ImportError:
    _sip = None
    _sip_isdeleted = None

try:
    from . import sip_isdeleted_safe
except ImportError:
    from modules import sip_isdeleted_safe


class MoleculeScene(QGraphicsScene):
    def clear_template_preview(self):
        """Remove all ghost lines for template preview."""
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(
                item, "_is_template_preview", False
            ):
                if sip_isdeleted_safe(item):
                    continue
                sc = getattr(item, "scene", lambda: None)()
                if sc is None:
                    continue
                try:
                    self.removeItem(item)
                except (RuntimeError, ValueError, TypeError):
                    # Best-effort: ignore removal errors during teardown if underlying C++ object is already gone
                    pass
        self.template_context = {}
        if hasattr(self, "template_preview"):
            self.template_preview.hide()

    def __init__(self, data, window):
        super().__init__()
        self.data, self.window = data, window
        self.mode, self.current_atom_symbol = "select", "C"
        self.bond_order, self.bond_stereo = 1, 0
        self.start_atom, self.temp_line, self.start_pos = None, None, None
        self.press_pos = None
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.hovered_item = None

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

    def get_setting(self, key, default=None):
        """Safe gateway to access MainWindow settings without deep traversal from items."""
        if hasattr(self, "window") and self.window and hasattr(self.window, "settings"):
            return self.window.settings.get(key, default)
        return default

    def update_connected_bonds(self, atoms):
        """Update the positions of all bonds connected to the specified atom list."""
        bonds_to_update = set()
        for atom in atoms:
            if hasattr(atom, "bonds"):
                bonds_to_update.update(atom.bonds)

        for bond in bonds_to_update:
            if not sip_isdeleted_safe(bond):
                try:
                    bond.update_position()
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    continue

    def update_all_items(self):
        """Force redraw of all items."""
        if hasattr(self.data, "update_ring_info_2d"):
            self.data.update_ring_info_2d()

        for item in self.items():
            if isinstance(item, (AtomItem, BondItem)):
                item.update()
        if self.views():
            self.views()[0].viewport().update()

    def reinitialize_items(self):
        self.template_preview = TemplatePreviewItem()
        self.addItem(self.template_preview)
        self.template_preview.hide()
        self.template_preview_points = []
        self.template_context = {}
        self._deleted_items = []

        app = QApplication.instance()
        if app is not None:
            try:
                app.aboutToQuit.connect(self.purge_deleted_items)
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Non-fatal during setup; app instance may be invalid or signal already connected
                pass

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
            if isinstance(item, AtomItem):
                try:
                    self.initial_positions_in_event[item] = item.pos()
                except RuntimeError:
                    # Skip if the object has been deleted (common in rapid UI updates)
                    continue

        if not self.window.is_2d_editable:
            return

        if event.button() == Qt.MouseButton.RightButton:
            item = self.itemAt(event.scenePos(), self.views()[0].transform())
            if not isinstance(item, (AtomItem, BondItem)):
                return  # Do nothing if something other than the target is clicked
            data_changed = False
            # If the user has a rectangular multi-selection and the clicked item
            # is part of that selection, delete all selected items (atoms/bonds).
            try:
                selected_items = [
                    it
                    for it in self.selectedItems()
                    if isinstance(it, (AtomItem, BondItem))
                ]
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Fallback to empty selection if the scene state is inconsistent during event processing
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
                    self.window.push_undo_state()
                self.press_pos = None
                event.accept()
                return
            # --- E/Z mode specific processing ---
            if self.mode == "bond_2_5":
                if isinstance(item, BondItem):
                    try:
                        # Clear E/Z label (revert to normal)
                        if item.stereo in [3, 4]:
                            item.set_stereo(0)
                            # Also update the data model
                            for (id1, id2), bdata in self.data.bonds.items():
                                if bdata.get("item") is item:
                                    bdata["stereo"] = 0
                                    break
                            self.window.push_undo_state()
                            data_changed = False  # Already added to undo stack, so skip redundant pushes later
                    except (AttributeError, RuntimeError, ValueError) as e:
                        logging.error(f"Error in E/Z stereo toggle: {e}", exc_info=True)
                        if hasattr(self.window, "statusBar"):
                            self.window.statusBar().showMessage(
                                f"Error clearing E/Z label: {e}", 5000
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
                self.window.push_undo_state()
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
        if not self.window.is_2d_editable:
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
        if not self.window.is_2d_editable:
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
                    sc = getattr(self.temp_line, "scene", lambda: None)()
                    if sc:
                        self.removeItem(self.temp_line)
                except (RuntimeError, ValueError, TypeError):
                    pass  # Suppress removal errors during teardown
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
                self.start_pos = None
                self.press_pos = None
                if self.data_changed_in_event:
                    self.update_all_items()
                    self.window.push_undo_state()
                return

        released_item = self.itemAt(end_pos, self.views()[0].transform())

        # 1. Handle special modes (radical/charge)
        if (
            (self.mode == "radical")
            and is_click
            and isinstance(released_item, AtomItem)
        ):
            atom = released_item
            atom.prepareGeometryChange()
            # Toggle radical state (0 -> 1 -> 2 -> 0)
            atom.radical = (atom.radical + 1) % 3
            self.data.atoms[atom.atom_id]["radical"] = atom.radical
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom = None
            self.start_pos = None
            self.press_pos = None
            if self.data_changed_in_event:
                self.window.push_undo_state()
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
            self.data.atoms[atom.atom_id]["charge"] = atom.charge
            atom.update_style()
            self.data_changed_in_event = True
            self.start_atom = None
            self.start_pos = None
            self.press_pos = None
            if self.data_changed_in_event:
                self.window.push_undo_state()
            return

        elif (
            self.mode.startswith("bond")
            and is_click
            and isinstance(released_item, BondItem)
        ):
            b = released_item
            if self.mode == "bond_2_5":
                try:
                    if b.order == 2:
                        current_stereo = b.stereo
                        if current_stereo not in [3, 4]:
                            new_stereo = 3  # None -> Z
                        elif current_stereo == 3:
                            new_stereo = 4  # Z -> E
                        else:  # current_stereo == 4
                            new_stereo = 0  # E -> None
                        self.update_bond_stereo(b, new_stereo)
                        self.update_all_items()  # Force redraw
                        self.window.push_undo_state()  # Push to undo stack here
                except (AttributeError, RuntimeError, ValueError) as e:
                    logging.error(f"Error in E/Z stereo toggle: {e}", exc_info=True)
                    if hasattr(self.window, "statusBar"):
                        self.window.statusBar().showMessage(
                            f"Error changing E/Z stereochemistry: {e}", 5000
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

        # Safely check for deleted objects
        moved_atoms = []
        for item, old_pos in self.initial_positions_in_event.items():
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
                    self.data.atoms[atom.atom_id]["pos"] = atom.pos()
                    bonds_to_update.update(atom.bonds)
                except RuntimeError:
                    # Skip if object is deleted
                    continue
            for bond in bonds_to_update:
                bond.update_position()
            # Update measurement label positions after atom move
            self.window.update_2d_measurement_labels()
            if self.views():
                self.views()[0].viewport().update()

        if self.data_changed_in_event:
            self.update_all_items()

        self.start_atom = None
        self.start_pos = None
        self.press_pos = None
        self.temp_line = None
        self.template_context = {}
        # Clear user template data when switching modes
        if hasattr(self, "user_template_data"):
            self.user_template_data = None
        if self.data_changed_in_event:
            self.window.push_undo_state()

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
            self.window.push_undo_state()

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
                    if a is None:
                        continue
                    if a in connected_atoms:
                        continue
                    connected_atoms.add(a)
                    # iterate bonds attached to atom
                    for b in getattr(a, "bonds", []) or []:
                        if b is None:
                            continue
                        connected_bonds.add(b)
                        # find the other atom at the bond
                        other = None
                        try:
                            if getattr(b, "atom1", None) is a:
                                other = getattr(b, "atom2", None)
                            else:
                                other = getattr(b, "atom1", None)
                        except (AttributeError, RuntimeError, ValueError, TypeError):
                            other = None
                        if other is not None and other not in connected_atoms:
                            atoms_to_visit.append(other)

                # Apply selection: clear previous and select only these
                self.clearSelection()

                for a in connected_atoms:
                    if not sip_isdeleted_safe(a):
                        try:
                            a.setSelected(True)
                        except (RuntimeError, ValueError, TypeError):
                            pass  # Ignore invalid item states during component search

                for b in connected_bonds:
                    if not sip_isdeleted_safe(b):
                        try:
                            b.setSelected(True)
                        except (RuntimeError, ValueError, TypeError):
                            pass  # Ignore invalid bond states during component search
                event.accept()
                return
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # On any unexpected error, fall back to default handling
                pass  # Suppress errors during connection search

        elif self.mode in ["bond_2_5"]:
            event.accept()
            return

        super().mouseDoubleClickEvent(event)

    def create_atom(self, symbol, pos, charge=0, radical=0):
        atom_id = self.data.add_atom(symbol, pos, charge=charge, radical=radical)
        atom_item = AtomItem(atom_id, symbol, pos, charge=charge, radical=radical)
        self.data.atoms[atom_id]["item"] = atom_item
        self.addItem(atom_item)
        return atom_id

    def create_bond(self, start_atom, end_atom, bond_order=None, bond_stereo=None):
        if start_atom is None or end_atom is None:
            logging.error("Error: Cannot create bond with None atoms")
            return

        try:
            exist_b = self.find_bond_between(start_atom, end_atom)
            if exist_b:
                return

            # Use the order specified in the argument if provided, otherwise use current defaults
            order_to_use = self.bond_order if bond_order is None else bond_order
            stereo_to_use = self.bond_stereo if bond_stereo is None else bond_stereo

            key, status = self.data.add_bond(
                start_atom.atom_id, end_atom.atom_id, order_to_use, stereo_to_use
            )
            if status == "created":
                bond_item = BondItem(start_atom, end_atom, order_to_use, stereo_to_use)
                self.data.bonds[key]["item"] = bond_item
                if hasattr(start_atom, "bonds"):
                    start_atom.bonds.append(bond_item)
                if hasattr(end_atom, "bonds"):
                    end_atom.bonds.append(bond_item)
                self.addItem(bond_item)

            if hasattr(start_atom, "update_style"):
                start_atom.update_style()
            if hasattr(end_atom, "update_style"):
                end_atom.update_style()

        except (AttributeError, RuntimeError, ValueError) as e:
            logging.error(f"Error creating bond: {e}", exc_info=True)
            self.update_all_items()

    def add_molecule_fragment(
        self, points, bonds_info, existing_items=None, symbol="C"
    ):
        """Add a molecular fragment (e.g., benzene template) to the scene.
        
        - Enforce policy of not changing existing bond orders.
        - Benzene template rotation is decided based on fused bond orders 
          so that exactly two new double bonds are created.
        """

        num_points = len(points)
        atom_items = [None] * num_points

        is_benzene_template = num_points == 6 and any(o == 2 for _, _, o in bonds_info)

        def coords(p):
            if hasattr(p, "x") and hasattr(p, "y"):
                return (p.x(), p.y())
            try:
                return (p[0], p[1])
            except (AttributeError, RuntimeError, ValueError, TypeError):
                raise ValueError("point has no x/y")

        def dist_pts(a, b):
            ax, ay = coords(a)
            bx, by = coords(b)
            return math.hypot(ax - bx, ay - by)

        # --- 1) Map already clicked existing_items to template vertices ---
        existing_items = existing_items or []
        used_indices = set()
        ref_lengths = [
            dist_pts(points[i], points[j])
            for i, j, _ in bonds_info
            if i < num_points and j < num_points
        ]
        avg_len = (sum(ref_lengths) / len(ref_lengths)) if ref_lengths else 20.0
        map_threshold = max(0.5 * avg_len, 8.0)

        for ex_item in existing_items:
            try:
                ex_pos = ex_item.pos()
                best_idx, best_d = -1, float("inf")
                for i, p in enumerate(points):
                    if i in used_indices:
                        continue
                    d = dist_pts(p, ex_pos)
                    if best_d is None or d < best_d:
                        best_d, best_idx = d, i
                if best_idx != -1 and best_d <= max(map_threshold, 1.5 * avg_len):
                    atom_items[best_idx] = ex_item
                    used_indices.add(best_idx)
            except (AttributeError, RuntimeError, ValueError, TypeError):  
                pass  # Suppress point distance mapping errors during fragment addition

        # --- 2) Enumerate existing atoms in the scene from self.data.atoms and map them ---
        mapped_atoms = {it for it in atom_items if it is not None}
        for i, p in enumerate(points):
            if atom_items[i] is not None:
                continue

            nearby = None
            best_d = float("inf")

            for atom_data in self.data.atoms.values():
                a_item = atom_data.get("item")
                if not a_item or a_item in mapped_atoms:
                    continue
                try:
                    d = dist_pts(p, a_item.pos())
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    continue
                if d < best_d:
                    best_d, nearby = d, a_item

            if nearby and best_d <= map_threshold:
                atom_items[i] = nearby
                mapped_atoms.add(nearby)

        # --- 3) Create missing vertices ---
        for i, p in enumerate(points):
            if atom_items[i] is None:
                atom_id = self.create_atom(symbol, p)
                atom_items[i] = self.data.atoms[atom_id]["item"]

        # --- 4) Determine bond array for the template (benzene rotation alignment processing) ---
        template_bonds_to_use = list(bonds_info)
        is_6ring = num_points == 6 and len(bonds_info) == 6
        template_has_double = any(o == 2 for (_, _, o) in bonds_info)

        if is_6ring and template_has_double:
            existing_orders = {}  # key: index of bonds_info, value: existing bond order
            for k, (i_idx, j_idx, _) in enumerate(bonds_info):
                if i_idx < len(atom_items) and j_idx < len(atom_items):
                    a, b = atom_items[i_idx], atom_items[j_idx]
                    if a is None or b is None:
                        continue
                    eb = self.find_bond_between(a, b)
                    if eb:
                        existing_orders[k] = getattr(eb, "order", 1)

            if existing_orders:
                orig_orders = [o for (_, _, o) in bonds_info]
                best_rot = 0
                max_score = -999  # Score means "goodness of fit"

                # --- Conditional branching based on the number of fused edges ---
                if len(existing_orders) >= 2:
                    for rot in range(num_points):
                        match_double_count = 0
                        match_bonus = 0
                        mismatch_penalty = 0

                        # Safety check: ensure fused region neighbors are single bonds to avoid valence overflow
                        safe_connection_score = 0

                        # Identify fused region boundaries (assuming continuity)
                        fused_indices = sorted(list(existing_orders.keys()))

                        # Check neighboring edges outside the fused region
                        for k in existing_orders:
                            # Check left neighbor
                            prev_idx = (k - 1 + rot) % num_points
                            # Check right neighbor
                            next_idx = (k + 1 + rot) % num_points

                        # Template bond order array
                        current_template_orders = [
                            orig_orders[(i + rot) % num_points]
                            for i in range(num_points)
                        ]

                        # To identify both ends of the fused region,
                        # collect template-side indices corresponding to "fusing k"
                        used_template_indices = set(
                            (k + rot) % num_points for k in existing_orders
                        )

                        # Score boost if template leg neighbors are single bonds
                        for t_idx in used_template_indices:
                            # Left template neighbor
                            adj_l = (t_idx - 1) % num_points
                            # Right template neighbor
                            adj_r = (t_idx + 1) % num_points

                            # Unused neighbors are connection "legs"
                            if adj_l not in used_template_indices:
                                if orig_orders[adj_l] == 1:
                                    safe_connection_score += 5000

                            if adj_r not in used_template_indices:
                                if orig_orders[adj_r] == 1:
                                    safe_connection_score += 5000

                        # Score calculation
                        for k, exist_order in existing_orders.items():
                            template_ord = orig_orders[(k + rot) % num_points]
                            if template_ord == exist_order:
                                match_bonus += 100
                                if exist_order == 2:
                                    match_double_count += 1
                            else:
                                # Modest penalty for mismatch if connection is otherwise safe
                                mismatch_penalty += 50

                        # Final score: connection safety as tip priority
                        current_score = (
                            safe_connection_score
                            + (match_double_count * 1000)
                            + match_bonus
                            - mismatch_penalty
                        )

                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot

                elif len(existing_orders) == 1:
                    # 1-edge fuse
                    k_fuse = next(iter(existing_orders.keys()))
                    exist_order = existing_orders[k_fuse]

                    for rot in range(num_points):
                        current_score = 0
                        rotated_template_order = orig_orders[
                            (k_fuse + rot) % num_points
                        ]

                        # 1. Leg order matching

                        # Pattern A: Alternating (opposite of existing)
                        if (exist_order == 1 and rotated_template_order == 2) or (
                            exist_order == 2 and rotated_template_order == 1
                        ):
                            current_score += 100

                        # Handle double bond superposition for conjugation
                        elif exist_order == 2 and rotated_template_order == 2:
                            current_score += 100

                        # 2. Verify alternating bond arrangement in adjacent edges
                        m_adj1 = (k_fuse - 1 + rot) % num_points
                        m_adj2 = (k_fuse + 1 + rot) % num_points
                        neighbor_order_1 = orig_orders[m_adj1]
                        neighbor_order_2 = orig_orders[m_adj2]

                        if exist_order == 1:
                            # If leg is single, neighbor should be double
                            if neighbor_order_1 == 2:
                                current_score += 50
                            if neighbor_order_2 == 2:
                                current_score += 50

                        elif exist_order == 2:
                            # If leg is double, neighbor should be single
                            if neighbor_order_1 == 1:
                                current_score += 50
                            if neighbor_order_2 == 1:
                                current_score += 50

                        # 3. Tie-break using consistency of non-contacting edges
                        for k, e_order in existing_orders.items():
                            if k != k_fuse:
                                r_t_order = orig_orders[(k + rot) % num_points]
                                if r_t_order == e_order:
                                    current_score += 10

                        if current_score > max_score:
                            max_score = current_score
                            best_rot = rot

                # Reflect final rotation
                new_tb = []
                for m in range(num_points):
                    i_idx, j_idx, _ = bonds_info[m]
                    new_order = orig_orders[(m + best_rot) % num_points]
                    new_tb.append((i_idx, j_idx, new_order))
                template_bonds_to_use = new_tb

        # --- 5) Bond Creation/Update ---
        for id1_idx, id2_idx, order in template_bonds_to_use:
            if id1_idx < len(atom_items) and id2_idx < len(atom_items):
                a_item, b_item = atom_items[id1_idx], atom_items[id2_idx]
                if not a_item or not b_item or a_item is b_item:
                    continue

                id1, id2 = a_item.atom_id, b_item.atom_id
                if id1 > id2:
                    id1, id2 = id2, id1

                exist_b = self.find_bond_between(a_item, b_item)

                if exist_b:
                    # Keep existing bonds by default
                    should_overwrite = False

                    # Check if benzene template and fused bond is single
                    if is_benzene_template and exist_b.order == 1:
                        # Check if fused single bond is part of conjugation
                        # (atoms at both ends should not have other double bonds)
                        atom1 = exist_b.atom1
                        atom2 = exist_b.atom2

                        # Check if atom1 has other double bonds
                        atom1_has_other_double_bond = any(
                            b.order == 2 for b in atom1.bonds if b is not exist_b
                        )

                        # Check atom2 for double bonds
                        atom2_has_other_double_bond = any(
                            b.order == 2 for b in atom2.bonds if b is not exist_b
                        )

                        # Overwrite if isolated single bond
                        if (
                            not atom1_has_other_double_bond
                            and not atom2_has_other_double_bond
                        ):
                            should_overwrite = True

                    if should_overwrite:
                        # Update bond order if conditions met
                        exist_b.order = order
                        exist_b.stereo = 0
                        self.data.bonds[(id1, id2)]["order"] = order
                        self.data.bonds[(id1, id2)]["stereo"] = 0
                        exist_b.update()
                    else:
                        # Keep existing bond if conditions not met
                        continue
                else:
                    # Create new bond
                    self.create_bond(a_item, b_item, bond_order=order, bond_stereo=0)

        # --- 6) Redraw ---
        for at in atom_items:
            try:
                if at:
                    at.update_style()
            except (AttributeError, RuntimeError, ValueError, TypeError):  
                pass  # Suppress style update errors during fragment addition

        return atom_items

    def update_template_preview(self, pos):  
        mode_parts = self.mode.split("_")

        # Check if this is a user template
        if len(mode_parts) >= 3 and mode_parts[1] == "user":
            self.update_user_template_preview(pos)
            return

        is_aromatic = False
        if mode_parts[1] == "benzene":
            n = 6
            is_aromatic = True
        else:
            try:
                n = int(mode_parts[1])
            except ValueError:
                return

        items_under = self.items(pos)  # top-most first
        item = None
        for it in items_under:
            if isinstance(it, (AtomItem, BondItem)):
                item = it
                break

        points, bonds_info = [], []
        l = DEFAULT_BOND_LENGTH
        self.template_context = {}

        if isinstance(item, AtomItem):
            p0 = item.pos()
            continuous_angle = math.atan2(pos.y() - p0.y(), pos.x() - p0.x())
            snap_angle_rad = math.radians(15)
            snapped_angle = round(continuous_angle / snap_angle_rad) * snap_angle_rad
            p1 = p0 + QPointF(l * math.cos(snapped_angle), l * math.sin(snapped_angle))
            points = self._calculate_polygon_from_edge(p0, p1, n)
            self.template_context["items"] = [item]

        elif isinstance(item, BondItem):
            # Snap to bond
            p0, p1 = item.atom1.pos(), item.atom2.pos()
            points = self._calculate_polygon_from_edge(
                p0, p1, n, cursor_pos=pos, use_existing_length=True
            )
            self.template_context["items"] = [item.atom1, item.atom2]

        else:
            angle_step = 2 * math.pi / n
            start_angle = -math.pi / 2 if n % 2 != 0 else -math.pi / 2 - angle_step / 2
            points = [
                pos
                + QPointF(
                    l * math.cos(start_angle + i * angle_step),
                    l * math.sin(start_angle + i * angle_step),
                )
                for i in range(n)
            ]

        if points:
            if is_aromatic:
                bonds_info = [
                    (i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)
                ]
            else:
                bonds_info = [(i, (i + 1) % n, 1) for i in range(n)]

            self.template_context["points"] = points
            self.template_context["bonds_info"] = bonds_info

            self.template_preview.set_geometry(points, is_aromatic)

            self.template_preview.show()
            if self.views():
                self.views()[0].viewport().update()
        else:
            self.template_preview.hide()
            if self.views():
                self.views()[0].viewport().update()

    def _calculate_polygon_from_edge(
        self, p0, p1, n, cursor_pos=None, use_existing_length=False
    ):
        if n < 3:
            return []
        v_edge = p1 - p0
        edge_length = math.sqrt(v_edge.x() ** 2 + v_edge.y() ** 2)
        if edge_length == 0:
            return []

        target_length = edge_length if use_existing_length else DEFAULT_BOND_LENGTH

        v_edge = (v_edge / edge_length) * target_length

        if not use_existing_length:
            p1 = p0 + v_edge

        points = [p0, p1]

        interior_angle = (n - 2) * math.pi / n
        rotation_angle = math.pi - interior_angle

        if cursor_pos:
            # Note: v_edge is normalized, direction is same
            v_cursor = cursor_pos - p0
            cross_product_z = (p1 - p0).x() * v_cursor.y() - (
                p1 - p0
            ).y() * v_cursor.x()
            if cross_product_z < 0:
                rotation_angle = -rotation_angle

        cos_a, sin_a = math.cos(rotation_angle), math.sin(rotation_angle)

        current_p, current_v = p1, v_edge
        for _ in range(n - 2):
            new_vx = current_v.x() * cos_a - current_v.y() * sin_a
            new_vy = current_v.x() * sin_a + current_v.y() * cos_a
            current_v = QPointF(new_vx, new_vy)
            current_p = current_p + current_v
            points.append(current_p)
        return points

    def delete_items(self, items_to_delete):
        """Safely delete specified items (atoms/bonds) in order"""
        if not items_to_delete:
            return False

        # 1. Sanitize input: only keep live QGraphicsItem wrappers
        sanitized = set()
        for it in items_to_delete:
            if it is not None and not sip_isdeleted_safe(it):
                if isinstance(it, (AtomItem, BondItem, QGraphicsItem)):
                    sanitized.add(it)
        items_to_delete = sanitized

        try:
            atoms_to_delete = {it for it in items_to_delete if isinstance(it, AtomItem)}
            bonds_to_delete = {it for it in items_to_delete if isinstance(it, BondItem)}

            # Include bonds attached to atoms being deleted
            for atom in list(atoms_to_delete):
                if hasattr(atom, "bonds") and atom.bonds:
                    for b in list(atom.bonds):
                        bonds_to_delete.add(b)

            # Determine surviving atoms whose bond lists need pruning
            atoms_to_update = set()
            for bond in list(bonds_to_delete):
                a1 = getattr(bond, "atom1", None)
                a2 = getattr(bond, "atom2", None)
                if a1 and a1 not in atoms_to_delete:
                    atoms_to_update.add(a1)
                if a2 and a2 not in atoms_to_delete:
                    atoms_to_update.add(a2)

            # 2. Update surviving atoms' bond lists
            for atom in list(atoms_to_update):
                if sip_isdeleted_safe(atom):
                    continue
                if hasattr(atom, "bonds") and atom.bonds:
                    try:
                        # Filter out bonds being deleted and SIP-stale wrappers
                        atom.bonds[:] = [
                            b for b in atom.bonds 
                            if not sip_isdeleted_safe(b) and b not in bonds_to_delete
                        ]
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass  # Suppress SIP-stale bond removal errors
                
                if hasattr(atom, "update_style"):
                    atom.update_style()

            # 3. Remove from data model
            for bond in list(bonds_to_delete):
                a1 = getattr(bond, "atom1", None)
                a2 = getattr(bond, "atom2", None)
                if a1 and a2 and hasattr(self, "data"):
                    try:
                        # Try both directions
                        if not self.data.remove_bond(a1.atom_id, a2.atom_id):
                            self.data.remove_bond(a2.atom_id, a1.atom_id)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass  # Suppress bond model data removal errors

            for atom in list(atoms_to_delete):
                if hasattr(atom, "atom_id") and hasattr(self, "data"):
                    try:
                        self.data.remove_atom(atom.atom_id)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass  # Suppress data model/graphic removal errors

            try:
                self._ih_update_counter = getattr(self, "_ih_update_counter", 0) + 1
            except (AttributeError, RuntimeError, ValueError, TypeError):
                # Suppress non-critical internal update counter increment errors.
                # This counter is only used for UI throttling and is non-critical for data integrity.
                pass
            # 4. Remove graphic items from the scene
            current_scene_items = set(self.items())
            
            # Helper to safely remove and hide items
            def safe_remove_and_hide(item_set):
                for item in list(item_set):
                    if sip_isdeleted_safe(item):
                        continue
                    if item in current_scene_items:
                        try:
                            self.removeItem(item)
                        except (RuntimeError, ValueError, TypeError):
                            pass  # Suppress graphic removal errors (stale pointers)
                    
                    try:
                        item.hide()
                        if not hasattr(self, "_deleted_items") or self._deleted_items is None:
                            self._deleted_items = []
                        self._deleted_items.append(item)
                    except (AttributeError, RuntimeError, ValueError, TypeError):
                        pass  # Suppress data model/graphic removal errors

            safe_remove_and_hide(bonds_to_delete)
            safe_remove_and_hide(atoms_to_delete)

            for atom in list(atoms_to_update):
                if hasattr(atom, "update_style"):
                    atom.update_style()

            return True

        except (AttributeError, RuntimeError, ValueError) as e:
            print(f"Error during delete_items operation: {e}")
            self.update_all_items()
            return False

    def purge_deleted_items(self):
        """Purge and release any held deleted-wrapper references during shutdown."""
        if not getattr(self, "_deleted_items", None):
            return

        for obj in list(self._deleted_items):
            if not sip_isdeleted_safe(obj):
                try:
                    obj.hide()
                    if hasattr(obj, "bonds") and obj.bonds is not None:
                        obj.bonds.clear()
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    pass  # Suppress errors during template point mapping

        try:
            self._deleted_items.clear()
        except (AttributeError, RuntimeError, ValueError, TypeError):
            self._deleted_items = []

    def add_user_template_fragment(self, context):
        """Place user template fragment"""
        points = context.get("points", [])
        bonds_info = context.get("bonds_info", [])
        atoms_data = context.get("atoms_data", [])
        attachment_atom = context.get("attachment_atom")

        if not points or not atoms_data:
            return

        # Create atoms
        atom_id_map = {}  # template id -> scene atom id

        for i, (pos, atom_data) in enumerate(zip(points, atoms_data)):
            # Skip first atom if attaching to existing atom
            if i == 0 and attachment_atom:
                atom_id_map[atom_data["id"]] = attachment_atom.atom_id
                continue

            symbol = atom_data.get("symbol", "C")
            charge = atom_data.get("charge", 0)
            radical = atom_data.get("radical", 0)

            atom_id = self.data.add_atom(symbol, pos, charge, radical)
            atom_id_map[atom_data["id"]] = atom_id

            # Create visual atom item
            atom_item = AtomItem(atom_id, symbol, pos, charge, radical)
            self.data.atoms[atom_id]["item"] = atom_item
            self.addItem(atom_item)

        # Create bonds (bonds_info is always id-based)
        # Create index-to-id conversion table first
        index_to_id = [atom_data.get("id", i) for i, atom_data in enumerate(atoms_data)]
        for bond_info in bonds_info:
            if isinstance(bond_info, (list, tuple)) and len(bond_info) >= 2:
                # Convert 1st/2nd elements of bonds_info to id if they are indices
                atom1_idx = bond_info[0]
                atom2_idx = bond_info[1]
                order = bond_info[2] if len(bond_info) > 2 else 1
                stereo = bond_info[3] if len(bond_info) > 3 else 0

                # Index-to-id conversion (leave as-is if already id)
                if isinstance(atom1_idx, int) and atom1_idx < len(index_to_id):
                    template_atom1_id = index_to_id[atom1_idx]
                else:
                    template_atom1_id = atom1_idx
                if isinstance(atom2_idx, int) and atom2_idx < len(index_to_id):
                    template_atom2_id = index_to_id[atom2_idx]
                else:
                    template_atom2_id = atom2_idx

                atom1_id = atom_id_map.get(template_atom1_id)
                atom2_id = atom_id_map.get(template_atom2_id)

                if atom1_id is not None and atom2_id is not None:
                    # Skip if bond already exists
                    existing_bond = None
                    if (atom1_id, atom2_id) in self.data.bonds:
                        existing_bond = (atom1_id, atom2_id)
                    elif (atom2_id, atom1_id) in self.data.bonds:
                        existing_bond = (atom2_id, atom1_id)

                    if not existing_bond:
                        bond_key, _ = self.data.add_bond(
                            atom1_id, atom2_id, order, stereo
                        )
                        # Create visual bond item
                        atom1_item = self.data.atoms[atom1_id]["item"]
                        atom2_item = self.data.atoms[atom2_id]["item"]
                        if atom1_item and atom2_item:
                            bond_item = BondItem(atom1_item, atom2_item, order, stereo)
                            self.data.bonds[bond_key]["item"] = bond_item
                            self.addItem(bond_item)
                            atom1_item.bonds.append(bond_item)
                            atom2_item.bonds.append(bond_item)

        # Update atom visuals
        for atom_id in atom_id_map.values():
            if atom_id in self.data.atoms and self.data.atoms[atom_id]["item"]:
                self.data.atoms[atom_id]["item"].update_style()

    def update_user_template_preview(self, pos):  
        """Update user template preview"""
        # Robust preview: avoid self.data.atoms for preview-only atoms
        if not hasattr(self, "user_template_data") or not self.user_template_data:
            return

        template_data = self.user_template_data
        atoms = template_data.get("atoms", [])
        bonds = template_data.get("bonds", [])

        if not atoms:
            return

        # Find attachment point (first atom or clicked item)
        items_under = self.items(pos)
        attachment_atom = None
        for item in items_under:
            if isinstance(item, AtomItem):
                attachment_atom = item
                break

        # Calculate template positions
        points = []
        # Find template bounds for centering
        center_x = 0.0
        center_y = 0.0
        if atoms:
            min_x = min(atom["x"] for atom in atoms)
            max_x = max(atom["x"] for atom in atoms)
            min_y = min(atom["y"] for atom in atoms)
            max_y = max(atom["y"] for atom in atoms)
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
        # Position template
        if attachment_atom:
            # Attach to existing atom
            attach_pos = attachment_atom.pos()
            offset_x = attach_pos.x() - atoms[0]["x"]
            offset_y = attach_pos.y() - atoms[0]["y"]
        else:
            # Center at cursor position
            offset_x = pos.x() - center_x
            offset_y = pos.y() - center_y
        # Calculate atom positions
        for atom in atoms:
            new_pos = QPointF(atom["x"] + offset_x, atom["y"] + offset_y)
            points.append(new_pos)
        # Create atom ID to index mapping (for preview only)
        atom_id_to_index = {}
        for i, atom in enumerate(atoms):
            atom_id = atom.get("id", i)
            atom_id_to_index[atom_id] = i
        # Generate bonds_info from template bonds
        bonds_info = []
        for bond in bonds:
            atom1_idx = atom_id_to_index.get(bond["atom1"])
            atom2_idx = atom_id_to_index.get(bond["atom2"])
            if atom1_idx is not None and atom2_idx is not None:
                order = bond.get("order", 1)
                stereo = bond.get("stereo", 0)
                bonds_info.append((atom1_idx, atom2_idx, order, stereo))
        # Preview: draw lines from points and bonds_info
        # Save placement context
        self.template_context = {
            "points": points,
            "bonds_info": bonds_info,
            "atoms_data": atoms,
            "attachment_atom": attachment_atom,
        }
        # Clear legacy preview items
        for item in list(self.items()):
            if isinstance(item, QGraphicsLineItem) and getattr(
                item, "_is_template_preview", False
            ):
                self.removeItem(item)

        # Draw preview using TemplatePreviewItem
        self.template_preview.set_user_template_geometry(points, bonds_info, atoms)
        self.template_preview.show()
        if self.views():
            self.views()[0].viewport().update()

    def leaveEvent(self, event):
        self.template_preview.hide()

    def set_hovered_item(self, item):
        """Record currently hovered item"""
        self.hovered_item = item

    def keyPressEvent(self, event):  
        view = self.views()[0]
        cursor_pos = view.mapToScene(view.mapFromGlobal(QCursor.pos()))
        item_at_cursor = self.itemAt(cursor_pos, view.transform())
        key = event.key()
        modifiers = event.modifiers()

        if not self.window.is_2d_editable:
            return

        if key == Qt.Key.Key_4:
            # Case 1: Cursor over atom/bond (one-shot placement)
            if isinstance(item_at_cursor, (AtomItem, BondItem)):
                # Set benzene template parameters
                n, is_aromatic = 6, True
                points, bonds_info, existing_items = [], [], []

                # Calculate placement like update_template_preview
                if isinstance(item_at_cursor, AtomItem):
                    p0 = item_at_cursor.pos()
                    l = DEFAULT_BOND_LENGTH
                    direction = QLineF(p0, cursor_pos).unitVector()
                    p1 = (
                        p0 + direction.p2() * l
                        if direction.length() > 0
                        else p0 + QPointF(l, 0)
                    )
                    points = self._calculate_polygon_from_edge(p0, p1, n)
                    existing_items = [item_at_cursor]

                elif isinstance(item_at_cursor, BondItem):
                    p0, p1 = item_at_cursor.atom1.pos(), item_at_cursor.atom2.pos()
                    points = self._calculate_polygon_from_edge(
                        p0, p1, n, cursor_pos=cursor_pos, use_existing_length=True
                    )
                    existing_items = [item_at_cursor.atom1, item_at_cursor.atom2]

                if points:
                    bonds_info = [
                        (i, (i + 1) % n, 2 if i % 2 == 0 else 1) for i in range(n)
                    ]

                    # Add fragment at calculated position
                    self.add_molecule_fragment(
                        points, bonds_info, existing_items=existing_items
                    )
                    self.update_all_items()
                    self.window.push_undo_state()

            # Case 2: Cursor over empty space (mode switch)
            else:
                self.window.set_mode_and_update_toolbar("template_benzene")

            event.accept()
            return

        # --- 0a. Change radical (.) ---
        if key == Qt.Key.Key_Period:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                for atom in target_atoms:
                    # Toggle radical state (0 -> 1 -> 2 -> 0)
                    atom.prepareGeometryChange()
                    atom.radical = (atom.radical + 1) % 3
                    self.data.atoms[atom.atom_id]["radical"] = atom.radical
                    atom.update_style()
                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 0b. Change charge (+/-) ---
        if key == Qt.Key.Key_Plus or key == Qt.Key.Key_Minus:
            target_atoms = []
            selected = self.selectedItems()
            if selected:
                target_atoms = [item for item in selected if isinstance(item, AtomItem)]
            elif isinstance(item_at_cursor, AtomItem):
                target_atoms = [item_at_cursor]

            if target_atoms:
                delta = 1 if key == Qt.Key.Key_Plus else -1
                for atom in target_atoms:
                    atom.prepareGeometryChange()
                    atom.charge += delta
                    self.data.atoms[atom.atom_id]["charge"] = atom.charge
                    atom.update_style()
                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 1. Atom operations (change symbol) ---
        if isinstance(item_at_cursor, AtomItem):
            new_symbol = None
            if (
                modifiers == Qt.KeyboardModifier.NoModifier
                and key in self.key_to_symbol_map
            ):
                new_symbol = self.key_to_symbol_map[key]
            elif (
                modifiers == Qt.KeyboardModifier.ShiftModifier
                and key in self.key_to_symbol_map_shift
            ):
                new_symbol = self.key_to_symbol_map_shift[key]

            if new_symbol and item_at_cursor.symbol != new_symbol:
                item_at_cursor.prepareGeometryChange()

                item_at_cursor.symbol = new_symbol
                self.data.atoms[item_at_cursor.atom_id]["symbol"] = new_symbol
                item_at_cursor.update_style()

                atoms_to_update = {item_at_cursor}
                for bond in item_at_cursor.bonds:
                    bond.update()
                    other_atom = (
                        bond.atom1 if bond.atom2 is item_at_cursor else bond.atom2
                    )
                    atoms_to_update.add(other_atom)

                for atom in atoms_to_update:
                    atom.update_style()

                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 2. Bond operations (change order/stereo) ---
        target_bonds = []
        if isinstance(item_at_cursor, BondItem):
            target_bonds = [item_at_cursor]
        else:
            target_bonds = [
                it for it in self.selectedItems() if isinstance(it, BondItem)
            ]

        if target_bonds:
            any_bond_changed = False
            for bond in target_bonds:
                # 1. Identify current bond key in data model
                id1, id2 = bond.atom1.atom_id, bond.atom2.atom_id
                current_key = None
                if (id1, id2) in self.data.bonds:
                    current_key = (id1, id2)
                elif (id2, id1) in self.data.bonds:
                    current_key = (id2, id1)

                if not current_key:
                    continue

                # 2. Save previous state
                old_order, old_stereo = bond.order, bond.stereo

                # 3. Update BondItem properties based on key
                if key == Qt.Key.Key_W:
                    if bond.stereo == 1:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1
                        bond.stereo = 1

                elif key == Qt.Key.Key_D:
                    if bond.stereo == 2:
                        bond_data = self.data.bonds.pop(current_key)
                        new_key = (current_key[1], current_key[0])
                        self.data.bonds[new_key] = bond_data
                        bond.atom1, bond.atom2 = bond.atom2, bond.atom1
                        bond.update_position()
                        was_reversed = True
                    else:
                        bond.order = 1
                        bond.stereo = 2

                elif key == Qt.Key.Key_1 and (bond.order != 1 or bond.stereo != 0):
                    bond.order = 1
                    bond.stereo = 0
                elif key == Qt.Key.Key_2 and (bond.order != 2 or bond.stereo != 0):
                    bond.order = 2
                    bond.stereo = 0
                elif key == Qt.Key.Key_3 and bond.order != 3:
                    bond.order = 3
                    bond.stereo = 0

                # 4. Update data model if changed
                if old_order != bond.order or old_stereo != bond.stereo:
                    any_bond_changed = True

                    # 5. Remove data with old key
                    bond_data = self.data.bonds.pop(current_key)
                    bond_data["order"] = bond.order
                    bond_data["stereo"] = bond.stereo

                    # 6. Determine new key and re-register
                    new_key_id1, new_key_id2 = bond.atom1.atom_id, bond.atom2.atom_id
                    if bond.stereo == 0:
                        if new_key_id1 > new_key_id2:
                            new_key_id1, new_key_id2 = new_key_id2, new_key_id1

                    new_key = (new_key_id1, new_key_id2)
                    self.data.bonds[new_key] = bond_data

                    bond.update()

            if any_bond_changed:
                self.update_all_items()
                self.window.push_undo_state()

            if key in [
                Qt.Key.Key_1,
                Qt.Key.Key_2,
                Qt.Key.Key_3,
                Qt.Key.Key_W,
                Qt.Key.Key_D,
            ]:
                event.accept()
                return

        if isinstance(self.hovered_item, BondItem) and self.hovered_item.order == 2:
            if event.key() == Qt.Key.Key_Z:
                self.update_bond_stereo(self.hovered_item, 3)  # Z-isomer
                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return
            elif event.key() == Qt.Key.Key_E:
                self.update_bond_stereo(self.hovered_item, 4)  # E-isomer
                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 3. Atom operations (add atom) ---
        if key in [Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3]:
            target_order = 1
            if key == Qt.Key.Key_2:
                target_order = 2
            elif key == Qt.Key.Key_3:
                target_order = 3

            start_atom = None
            if isinstance(item_at_cursor, AtomItem):
                start_atom = item_at_cursor
            else:
                selected_atoms = [
                    item for item in self.selectedItems() if isinstance(item, AtomItem)
                ]
                if len(selected_atoms) == 1:
                    start_atom = selected_atoms[0]

            if start_atom:
                start_pos = start_atom.pos()
                l = DEFAULT_BOND_LENGTH
                new_pos_offset = QPointF(0, -l)  # Default offset (up)

                # Get non-H neighbors
                neighbor_positions = []
                for bond in start_atom.bonds:
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    if (
                        other_atom.symbol != "H"
                    ):  # Ignore H
                        neighbor_positions.append(other_atom.pos())

                num_non_H_neighbors = len(neighbor_positions)

                if num_non_H_neighbors == 0:
                    # Zero bonds: default direction
                    new_pos_offset = QPointF(0, -l)

                elif num_non_H_neighbors == 1:
                    # One bond: ~120/60 degree angle
                    bond = start_atom.bonds[0]
                    other_atom = bond.atom1 if bond.atom2 is start_atom else bond.atom2
                    existing_bond_vector = start_pos - other_atom.pos()

                    # Rotate 60° clockwise from existing bond
                    angle_rad = math.radians(60)
                    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
                    vx, vy = existing_bond_vector.x(), existing_bond_vector.y()
                    new_vx, new_vy = vx * cos_a - vy * sin_a, vx * sin_a + vy * cos_a
                    rotated_vector = QPointF(new_vx, new_vy)
                    line = QLineF(QPointF(0, 0), rotated_vector)
                    line.setLength(l)
                    new_pos_offset = line.p2()

                elif num_non_H_neighbors == 3:
                    bond_vectors_sum = QPointF(0, 0)
                    for pos in neighbor_positions:
                        # Vector from start_pos to neighbor_pos
                        vec = pos - start_pos
                        # Convert to unit vector
                        line_to_other = QLineF(QPointF(0, 0), vec)
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2()

                    # SUM_TOLERANCE is now a module-level constant
                    if bond_vectors_sum.manhattanLength() > SUM_TOLERANCE:
                        new_direction_line = QLineF(QPointF(0, 0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        new_pos_offset = QPointF(l * 0.7071, -l * 0.7071)

                else:  # 2, 4+ bonds: skeleton continuation or over-bonding
                    bond_vectors_sum = QPointF(0, 0)
                    for bond in start_atom.bonds:
                        other_atom = (
                            bond.atom1 if bond.atom2 is start_atom else bond.atom2
                        )
                        line_to_other = QLineF(start_pos, other_atom.pos())
                        if line_to_other.length() > 0:
                            line_to_other.setLength(1.0)
                            bond_vectors_sum += line_to_other.p2() - line_to_other.p1()

                    if bond_vectors_sum.manhattanLength() > 0.01:
                        new_direction_line = QLineF(QPointF(0, 0), -bond_vectors_sum)
                        new_direction_line.setLength(l)
                        new_pos_offset = new_direction_line.p2()
                    else:
                        # Default (up) if sum is zero
                        new_pos_offset = QPointF(0, -l)

                # SNAP_DISTANCE is a module-level constant
                target_pos = start_pos + new_pos_offset

                # Find nearby atom
                near_atom = self.find_atom_near(target_pos, tol=SNAP_DISTANCE)

                if near_atom and near_atom is not start_atom:
                    # Bond if exists
                    self.create_bond(
                        start_atom, near_atom, bond_order=target_order, bond_stereo=0
                    )
                else:
                    # Create new atom and bond
                    new_atom_id = self.create_atom("C", target_pos)
                    new_atom_item = self.data.atoms[new_atom_id]["item"]
                    self.create_bond(
                        start_atom,
                        new_atom_item,
                        bond_order=target_order,
                        bond_stereo=0,
                    )

                self.clearSelection()
                self.update_all_items()
                self.window.push_undo_state()
                event.accept()
                return

        # --- 4. Global operations (delete, mode switch) ---
        if key == Qt.Key.Key_Delete or key == Qt.Key.Key_Backspace:
            if self.temp_line:
                try:
                    if not sip_isdeleted_safe(self.temp_line):
                        try:
                            if (
                                getattr(self.temp_line, "scene", None)
                                and self.temp_line.scene()
                            ):
                                self.removeItem(self.temp_line)
                        except (AttributeError, RuntimeError, ValueError, TypeError):  
                            pass  # Suppress bond visual state sync errors
                except (AttributeError, RuntimeError, ValueError, TypeError):
                    try:
                        self.removeItem(self.temp_line)
                    except (AttributeError, RuntimeError, ValueError, TypeError):  
                        pass  # Suppress atom visual state sync errors

                self.temp_line = None
                self.start_atom = None
                self.start_pos = None
                self.initial_positions_in_event = {}
                event.accept()
                return

            items_to_process = set(self.selectedItems())
            # Include item under cursor in deletion
            if item_at_cursor and isinstance(item_at_cursor, (AtomItem, BondItem)):
                items_to_process.add(item_at_cursor)

            if self.delete_items(items_to_process):
                self.update_all_items()
                self.window.push_undo_state()
                self.window.statusBar().showMessage("Deleted selected items.")

            # Clear scene if no atoms left
            if not self.data.atoms:
                # 1. Remove all graphics items
                self.clear()

                # 2. Re-initialize required items
                self.reinitialize_items()

                # 3. Reset temporary states
                self.temp_line = None
                self.start_atom = None
                self.start_pos = None
                self.initial_positions_in_event = {}

                # Event handled
                event.accept()
                return

            # Force redraw
            if self.views():
                self.views()[0].viewport().update()
                QApplication.processEvents()

                event.accept()
                return

        if key == Qt.Key.Key_Space:
            if self.mode != "select":
                self.window.activate_select_mode()
            else:
                self.window.select_all()
            event.accept()
            return

        # Global drawing mode switch
        mode_to_set = None

        # 1. Switch to atom mode
        symbol_for_mode_change = None
        if (
            modifiers == Qt.KeyboardModifier.NoModifier
            and key in self.key_to_symbol_map
        ):
            symbol_for_mode_change = self.key_to_symbol_map[key]
        elif (
            modifiers == Qt.KeyboardModifier.ShiftModifier
            and key in self.key_to_symbol_map_shift
        ):
            symbol_for_mode_change = self.key_to_symbol_map_shift[key]

        if symbol_for_mode_change:
            mode_to_set = f"atom_{symbol_for_mode_change}"

        # 2. Switch to bond mode
        elif (
            modifiers == Qt.KeyboardModifier.NoModifier
            and key in self.key_to_bond_mode_map
        ):
            mode_to_set = self.key_to_bond_mode_map[key]

        # Execute mode change
        if mode_to_set:
            if hasattr(self.window, "set_mode_and_update_toolbar"):
                self.window.set_mode_and_update_toolbar(mode_to_set)
                event.accept()
                return

        # Fallback
        super().keyPressEvent(event)

    def find_atom_near(self, pos, tol=14.0):
        # Create a small search rectangle around the position
        search_rect = QRectF(pos.x() - tol, pos.y() - tol, 2 * tol, 2 * tol)
        nearby_items = self.items(search_rect)

        for it in nearby_items:
            if isinstance(it, AtomItem):
                # Check the precise distance only for candidate items
                if QLineF(it.pos(), pos).length() <= tol:
                    return it
        return None

    def find_bond_between(self, atom1, atom2):
        for b in atom1.bonds:
            if (b.atom1 is atom1 and b.atom2 is atom2) or (
                b.atom1 is atom2 and b.atom2 is atom1
            ):
                return b
        return None

    def update_bond_stereo(self, bond_item, new_stereo):
        """Update bond stereochemistry"""
        if bond_item is None:
            return

        try:
            if bond_item.order != 2 or bond_item.stereo == new_stereo:
                return

            a1 = getattr(bond_item, "atom1", None)
            a2 = getattr(bond_item, "atom2", None)
            if not a1 or not a2 or not hasattr(a1, "atom_id") or not hasattr(a2, "atom_id"):
                return

            id1, id2 = a1.atom_id, a2.atom_id
            key_to_update = (id1, id2)
            
            if key_to_update not in self.data.bonds:
                key_to_update = (id2, id1)
                if key_to_update not in self.data.bonds:
                    if hasattr(self.window, "statusBar"):
                        self.window.statusBar().showMessage(
                            f"Warning: Bond {id1}-{id2} not found in model.", 3000
                        )
                    return

            # Update data model and visual representation
            self.data.bonds[key_to_update]["stereo"] = new_stereo
            bond_item.set_stereo(new_stereo)
            self.data_changed_in_event = True

        except (AttributeError, RuntimeError, ValueError) as e:
            pass  # Suppress final coordinate adjustment errors
            if hasattr(self.window, "statusBar"):
                self.window.statusBar().showMessage(f"Error: {e}", 5000)
            self.update_all_items()
