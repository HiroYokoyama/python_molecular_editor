#!/usr/bin/env python3
"""Fix reaction_sketcher/patcher.py for V3 manager architecture."""

import re

PATCHER = r'tmp/plugins/Utility/reaction_sketcher/patcher.py'

with open(PATCHER, 'r', encoding='utf-8') as f:
    src = f.read()

orig = src  # keep for diff

def rep(old, new, src=None):
    global content
    if src is None:
        src = content
    if old not in src:
        print(f"  [MISS] {old[:60]!r}")
        return src
    n = src.count(old)
    result = src.replace(old, new)
    print(f"  [OK x{n}] {old[:60]!r}")
    return result

content = src

# ── 1. View2D class resolution (uses main_window directly) ──────────────────
content = rep(
    "    if hasattr(main_window, 'view_2d') and main_window.view_2d:\n        View2D = main_window.view_2d.__class__",
    "    if hasattr(main_window, 'init_manager') and main_window.init_manager.view_2d:\n        View2D = main_window.init_manager.view_2d.__class__"
)

# ── 2. data class resolution (uses main_window directly) ────────────────────
content = rep(
    "    # Fallback to instance inspection if available (Safest)\n    if hasattr(main_window, 'data'):",
    "    # Fallback to instance inspection if available (Safest)\n    if hasattr(main_window, 'state_manager'):"
)
content = rep(
    "        if AtomItem is None and main_window.data.atoms:",
    "        if AtomItem is None and main_window.state_manager.data.atoms:"
)
content = rep(
    "            for d in main_window.data.atoms.values():",
    "            for d in main_window.state_manager.data.atoms.values():"
)
content = rep(
    "        if BondItem is None and main_window.data.bonds:",
    "        if BondItem is None and main_window.state_manager.data.bonds:"
)
content = rep(
    "            for d in main_window.data.bonds.values():",
    "            for d in main_window.state_manager.data.bonds.values():"
)

# ── 3. patched_set_mode: self -> self.host for _reaction_mode_manager ────────
content = rep(
    "        rmm = getattr(self, '_reaction_mode_manager', None)",
    "        rmm = getattr(self.host, '_reaction_mode_manager', None)"
)

# ── 4. MainWindow delegator lambdas: old manager name ───────────────────────
content = rep(
    "lambda self: self.main_window_edit_actions.select_all()",
    "lambda self: self.edit_actions_manager.select_all()"
)
content = rep(
    "lambda self: self.main_window_edit_actions.delete_selection()",
    "lambda self: self.edit_actions_manager.delete_selection()"
)
content = rep(
    "lambda self: self.main_window_edit_actions.copy_2d_image_to_clipboard()",
    "lambda self: self.edit_actions_manager.copy_2d_image_to_clipboard()"
)

# ── 5. window.push_undo_state in patched_delete_items (MoleculeScene) ────────
content = rep(
    '                if hasattr(window, "push_undo_state"):\n                    window.push_undo_state()',
    '                if hasattr(window, "edit_actions_manager"):\n                    window.edit_actions_manager.push_undo_state()'
)

# ── 6. patched_push_undo_state (on StateManager): self.current_mol ───────────
content = rep(
    "            'mol_3d': self.current_mol.ToBinary() if self.current_mol else None,",
    "            'mol_3d': self.host.view_3d_manager.current_mol.ToBinary() if self.host.view_3d_manager.current_mol else None,"
)

# ── 7. patched_get_current_state (on StateManager): self.scene ───────────────
content = rep(
    "        for item in self.scene.items():\n            if hasattr(item, \"create_json_data\"):\n                 rs_items_data.append(item.create_json_data())",
    "        for item in self.host.scene.items():\n            if hasattr(item, \"create_json_data\"):\n                 rs_items_data.append(item.create_json_data())"
)

# ── 8. patched_set_state_from_data (on StateManager): self.scene, load_handler_core ──
content = rep(
    "        for item in list(self.scene.items()):\n            if hasattr(item, \"create_json_data\"):\n                 self.scene.removeItem(item)",
    "        for item in list(self.host.scene.items()):\n            if hasattr(item, \"create_json_data\"):\n                 self.host.scene.removeItem(item)"
)
content = rep(
    "            load_handler_core(self, state_data['rs_items'])",
    "            load_handler_core(self.host, state_data['rs_items'])"
)

# ── 9. patched_copy_selection (on EditActionsManager) ────────────────────────
content = rep(
    "            selected_items = [i for i in self.scene.selectedItems() if not sip_isdeleted_safe(i)]\n            selected_atoms = [item for item in selected_items if hasattr(item, 'atom_id')]\n            selected_rs_items_raw",
    "            selected_items = [i for i in self.host.scene.selectedItems() if not sip_isdeleted_safe(i)]\n            selected_atoms = [item for item in selected_items if hasattr(item, 'atom_id')]\n            selected_rs_items_raw"
)
content = rep(
    "            for (id1, id2), bond_data in self.data.bonds.items():",
    "            for (id1, id2), bond_data in self.host.state_manager.data.bonds.items():"
)
content = rep(
    '            self.statusBar().showMessage(f"Copied selection ({len(fragment_atoms)} atoms, {len(fragment_rs_items)} reaction items).")',
    '            self.host.statusBar().showMessage(f"Copied selection ({len(fragment_atoms)} atoms, {len(fragment_rs_items)} reaction items).")'
)
content = rep(
    '            self.statusBar().showMessage(f"Error during patched copy: {e}")',
    '            self.host.statusBar().showMessage(f"Error during patched copy: {e}")'
)

# ── 10. patched_paste_from_clipboard (on EditActionsManager) ─────────────────
content = rep(
    "            paste_center_pos = self.view_2d.mapToScene(self.view_2d.mapFromGlobal(QCursor.pos()))\n            self.scene.clearSelection()",
    "            paste_center_pos = self.host.init_manager.view_2d.mapToScene(self.host.init_manager.view_2d.mapFromGlobal(QCursor.pos()))\n            self.host.scene.clearSelection()"
)
content = rep(
    "                new_id = self.scene.create_atom(",
    "                new_id = self.host.scene.create_atom("
)
content = rep(
    "                item = self.data.atoms[new_id]['item']",
    "                item = self.host.state_manager.data.atoms[new_id]['item']"
)
content = rep(
    "                self.scene.create_bond(",
    "                self.host.scene.create_bond("
)
content = rep(
    "                load_handler_core(self, rs_items_data)",
    "                load_handler_core(self.host, rs_items_data)"
)
content = rep(
    '            self.statusBar().showMessage("Pasted selection.")',
    '            self.host.statusBar().showMessage("Pasted selection.")'
)
content = rep(
    "            if hasattr(self, 'activate_select_mode'):\n                self.activate_select_mode()",
    "            if hasattr(self.host, 'ui_manager'):\n                self.host.ui_manager.activate_select_mode()"
)
content = rep(
    '            self.statusBar().showMessage(f"Error during patched paste: {e}")',
    '            self.host.statusBar().showMessage(f"Error during patched paste: {e}")'
)

# ── 11. patched_delete_selection (on EditActionsManager) ─────────────────────
content = rep(
    "        items = self.scene.selectedItems()\n        if not items: return\n        # Delegate to patched scene.delete_items which handles separation\n        self.scene.delete_items(items)",
    "        items = self.host.scene.selectedItems()\n        if not items: return\n        # Delegate to patched scene.delete_items which handles separation\n        self.host.scene.delete_items(items)"
)

# ── 12. patched_select_all (on EditActionsManager) ───────────────────────────
content = rep(
    "        for item in self.scene.items():\n            if sip_isdeleted_safe(item):\n                continue\n            if hasattr(item, \"create_json_data\") or isinstance(item, (AtomItem, BondItem)):\n                item.setSelected(True)",
    "        for item in self.host.scene.items():\n            if sip_isdeleted_safe(item):\n                continue\n            if hasattr(item, \"create_json_data\") or isinstance(item, (AtomItem, BondItem)):\n                item.setSelected(True)"
)

# ── 13. patched_rotate_molecule_2d (on EditActionsManager) ───────────────────
content = rep(
    "            selected_items = [i for i in self.scene.selectedItems() if not sip_isdeleted_safe(i)]\n            \n            # Identify targets",
    "            selected_items = [i for i in self.host.scene.selectedItems() if not sip_isdeleted_safe(i)]\n            \n            # Identify targets"
)
content = rep(
    "                target_atoms = [data['item'] for data in self.data.atoms.values() if data.get('item')]",
    "                target_atoms = [data['item'] for data in self.host.state_manager.data.atoms.values() if data.get('item')]"
)
content = rep(
    "                for item in self.scene.items():\n                    if sip_isdeleted_safe(item):\n                        continue\n                    if hasattr(item, \"rotate_around\"):\n                        target_reaction_items.append(item)",
    "                for item in self.host.scene.items():\n                    if sip_isdeleted_safe(item):\n                        continue\n                    if hasattr(item, \"rotate_around\"):\n                        target_reaction_items.append(item)"
)
content = rep(
    '                self.statusBar().showMessage("No items to rotate.")',
    '                self.host.statusBar().showMessage("No items to rotate.")'
)
content = rep(
    "            self.scene.update_connected_bonds(target_atoms)\n            \n            self.push_undo_state()\n            self.statusBar().showMessage(f\"Rotated {len(target_atoms) + len(target_reaction_items)} items by {angle_degrees} degrees.\")\n            self.scene.update()\n            self.scene.update_all_items()",
    "            self.host.scene.update_connected_bonds(target_atoms)\n            \n            self.push_undo_state()\n            self.host.statusBar().showMessage(f\"Rotated {len(target_atoms) + len(target_reaction_items)} items by {angle_degrees} degrees.\")\n            self.host.scene.update()\n            self.host.scene.update_all_items()"
)
content = rep(
    '            self.statusBar().showMessage(f"Error rotating: {e}")',
    '            self.host.statusBar().showMessage(f"Error rotating: {e}")'
)

# ── 14. patched_clean_up_2d_structure (on EditActionsManager) ────────────────
content = rep(
    '            self.statusBar().showMessage("Error: RDKit is required for structure optimization.")',
    '            self.host.statusBar().showMessage("Error: RDKit is required for structure optimization.")'
)
content = rep(
    '        self.statusBar().showMessage("Optimizing 2D structure (CoG Preserved)...")\n        self.scene.clear_all_problem_flags()\n        if not self.data.atoms:\n            self.statusBar().showMessage("Error: No atoms to optimize.")',
    '        self.host.statusBar().showMessage("Optimizing 2D structure (CoG Preserved)...")\n        self.host.scene.clear_all_problem_flags()\n        if not self.host.state_manager.data.atoms:\n            self.host.statusBar().showMessage("Error: No atoms to optimize.")'
)
content = rep(
    "            mol = self.data.to_rdkit_mol()",
    "            mol = self.host.state_manager.data.to_rdkit_mol()"
)
content = rep(
    "                self.check_chemistry_problems_fallback()",
    "                self.host.compute_manager.check_chemistry_problems_fallback()"
)
content = rep(
    "            selected_items = self.scene.selectedItems()\n            target_atom_ids = set()",
    "            selected_items = self.host.scene.selectedItems()\n            target_atom_ids = set()"
)
content = rep(
    "                    if aid in self.data.atoms:\n                        pos = self.data.atoms[aid]['item'].pos()",
    "                    if aid in self.host.state_manager.data.atoms:\n                        pos = self.host.state_manager.data.atoms[aid]['item'].pos()"
)
content = rep(
    "                    if aid in self.data.atoms:\n                        item = self.data.atoms[aid]['item']\n                        rd_pos = conf.GetAtomPosition(idx)\n                        sx",
    "                    if aid in self.host.state_manager.data.atoms:\n                        item = self.host.state_manager.data.atoms[aid]['item']\n                        rd_pos = conf.GetAtomPosition(idx)\n                        sx"
)
content = rep(
    "                        self.data.atoms[aid]['pos'] = new_pos",
    "                        self.host.state_manager.data.atoms[aid]['pos'] = new_pos"
)
content = rep(
    "            for bond_data in self.data.bonds.values():",
    "            for bond_data in self.host.state_manager.data.bonds.values():"
)
content = rep(
    "            self.update_2d_measurement_labels()\n            self.scene.update()",
    "            self.host.edit_3d_manager.update_2d_measurement_labels()\n            self.host.scene.update()"
)
content = rep(
    '            self.statusBar().showMessage(msg)\n            self.push_undo_state()',
    '            self.host.statusBar().showMessage(msg)\n            self.push_undo_state()'
)
content = rep(
    '            self.statusBar().showMessage(f"Error during CoG optimization: {e}")',
    '            self.host.statusBar().showMessage(f"Error during CoG optimization: {e}")'
)
content = rep(
    "        if hasattr(self, 'view_2d') and self.view_2d:\n                self.view_2d.setFocus()",
    "        if hasattr(self.host, 'init_manager') and self.host.init_manager.view_2d:\n                self.host.init_manager.view_2d.setFocus()"
)

# ── 15. Export methods (on ExportManager): self.scene, self.data, etc. ────────
# _render_2d_to_image
content = rep(
    "            all_visible_items = [i for i in self.scene.items() if i.isVisible()]",
    "            all_visible_items = [i for i in self.host.scene.items() if i.isVisible()]"
)
content = rep(
    "            original_background = self.scene.backgroundBrush()\n            if is_transparent:\n                # Set truly transparent WHITE background",
    "            original_background = self.host.scene.backgroundBrush()\n            if is_transparent:\n                # Set truly transparent WHITE background"
)
content = rep(
    "                self.scene.setBackgroundBrush(QBrush(QColor(255, 255, 255, 0)))\n            \n            # Clear selection and focus to avoid highlights/cursors in export\n            selected_items = self.scene.selectedItems()\n            self.scene.clearSelection()\n            original_focus = self.scene.focusItem()\n            self.scene.setFocusItem(None)",
    "                self.host.scene.setBackgroundBrush(QBrush(QColor(255, 255, 255, 0)))\n            \n            # Clear selection and focus to avoid highlights/cursors in export\n            selected_items = self.host.scene.selectedItems()\n            self.host.scene.clearSelection()\n            original_focus = self.host.scene.focusItem()\n            self.host.scene.setFocusItem(None)"
)
content = rep(
    "                    self.scene.render(painter, QRectF(0, 0, w, h), rect_to_render)",
    "                    self.host.scene.render(painter, QRectF(0, 0, w, h), rect_to_render)"
)
content = rep(
    "                self.scene.setBackgroundBrush(original_background)\n                for item in selected_items:\n                    item.setSelected(True)\n                return None, None\n\n            self.scene.setBackgroundBrush(original_background)",
    "                self.host.scene.setBackgroundBrush(original_background)\n                for item in selected_items:\n                    item.setSelected(True)\n                return None, None\n\n            self.host.scene.setBackgroundBrush(original_background)"
)

# patched_export_2d_png
content = rep(
    "            if not self.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.scene.items()):\n                self.statusBar().showMessage(\"Nothing to export.\")\n                return\n\n            default_name = \"untitled-2d\"\n            try:\n                if self.current_file_path:",
    "            if not self.host.state_manager.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.host.scene.items()):\n                self.host.statusBar().showMessage(\"Nothing to export.\")\n                return\n\n            default_name = \"untitled-2d\"\n            try:\n                if self.host.init_manager.current_file_path:"
)
content = rep(
    "                    default_name = os.path.splitext(os.path.basename(self.current_file_path))[0] + \"-2d\"\n            except: pass\n\n            filePath, _ = QFileDialog.getSaveFileName(self, \"Export 2D as PNG\"",
    "                    default_name = os.path.splitext(os.path.basename(self.host.init_manager.current_file_path))[0] + \"-2d\"\n            except: pass\n\n            filePath, _ = QFileDialog.getSaveFileName(self, \"Export 2D as PNG\""
)
content = rep(
    '                self.statusBar().showMessage(f"2D view exported to {filePath}")\n            else:\n                self.statusBar().showMessage("Failed to save image.")',
    '                self.host.statusBar().showMessage(f"2D view exported to {filePath}")\n            else:\n                self.host.statusBar().showMessage("Failed to save image.")'
)

# patched_copy_to_clipboard (ExportManager)
content = rep(
    "            if not self.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.scene.items()):\n                self.statusBar().showMessage(\"Nothing to copy.\")\n                return\n\n            # Default to Transparent",
    "            if not self.host.state_manager.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.host.scene.items()):\n                self.host.statusBar().showMessage(\"Nothing to copy.\")\n                return\n\n            # Default to Transparent"
)
content = rep(
    '                self.statusBar().showMessage("Copied 2D view to clipboard (Transparent)")\n            else:\n                self.statusBar().showMessage("Failed to copy image.")',
    '                self.host.statusBar().showMessage("Copied 2D view to clipboard (Transparent)")\n            else:\n                self.host.statusBar().showMessage("Failed to copy image.")'
)

# patched_export_2d_svg
content = rep(
    "            if not self.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.scene.items()):\n                self.statusBar().showMessage(\"Nothing to export.\")\n                return\n\n            default_name = \"untitled-2d\"\n            try:\n                if self.current_file_path:\n                    default_name = os.path.splitext(os.path.basename(self.current_file_path))[0] + \"-2d\"",
    "            if not self.host.state_manager.data.atoms and not any(hasattr(i, \"create_json_data\") for i in self.host.scene.items()):\n                self.host.statusBar().showMessage(\"Nothing to export.\")\n                return\n\n            default_name = \"untitled-2d\"\n            try:\n                if self.host.init_manager.current_file_path:\n                    default_name = os.path.splitext(os.path.basename(self.host.init_manager.current_file_path))[0] + \"-2d\""
)
content = rep(
    "            selected_items = [i for i in self.scene.selectedItems() if i.isVisible()]\n            items_to_export = selected_items if selected_items else list(self.scene.items())",
    "            selected_items = [i for i in self.host.scene.selectedItems() if i.isVisible()]\n            items_to_export = selected_items if selected_items else list(self.host.scene.items())"
)
content = rep(
    '                self.statusBar().showMessage("Error: Could not determine molecule bounds for export.")',
    '                self.host.statusBar().showMessage("Error: Could not determine molecule bounds for export.")'
)
content = rep(
    "            original_background = self.scene.backgroundBrush()\n            if reply == QMessageBox.StandardButton.Yes:\n                # Use strictly transparent WHITE background\n                self.scene.setBackgroundBrush(QBrush(QColor(255, 255, 255, 0)))\n            else:\n                # Use white background if user chose No transparency\n                self.scene.setBackgroundBrush(QBrush(Qt.GlobalColor.white))",
    "            original_background = self.host.scene.backgroundBrush()\n            if reply == QMessageBox.StandardButton.Yes:\n                # Use strictly transparent WHITE background\n                self.host.scene.setBackgroundBrush(QBrush(QColor(255, 255, 255, 0)))\n            else:\n                # Use white background if user chose No transparency\n                self.host.scene.setBackgroundBrush(QBrush(Qt.GlobalColor.white))"
)
content = rep(
    "            # Clear selection to avoid highlights in export\n            selected_items = self.scene.selectedItems()\n            self.scene.clearSelection()\n            self.original_focus = self.scene.focusItem()\n            self.scene.setFocusItem(None)",
    "            # Clear selection to avoid highlights in export\n            selected_items = self.host.scene.selectedItems()\n            self.host.scene.clearSelection()\n            self.original_focus = self.host.scene.focusItem()\n            self.host.scene.setFocusItem(None)"
)
content = rep(
    "                self.scene.setBackgroundBrush(original_background)\n                for item in selected_items:\n                    item.setSelected(True)\n                self.statusBar().showMessage(\"Failed to start SVG painter. Check file access.\")",
    "                self.host.scene.setBackgroundBrush(original_background)\n                for item in selected_items:\n                    item.setSelected(True)\n                self.host.statusBar().showMessage(\"Failed to start SVG painter. Check file access.\")"
)
content = rep(
    "                self.scene.render(painter, target_rect, rect_to_render)\n            finally:\n                painter.end()\n                self.scene.setBackgroundBrush(original_background)",
    "                self.host.scene.render(painter, target_rect, rect_to_render)\n            finally:\n                painter.end()\n                self.host.scene.setBackgroundBrush(original_background)"
)
content = rep(
    '            self.statusBar().showMessage(f"2D view exported to {filePath}")\n\n        patch_core(MainWindowExport',
    '            self.host.statusBar().showMessage(f"2D view exported to {filePath}")\n\n        patch_core(MainWindowExport'
)

# patched_copy_svg_to_clipboard
content = rep(
    '            if not self.data.atoms and not any(hasattr(i, "create_json_data") for i in self.scene.items()):\n                self.statusBar().showMessage("Nothing to copy.")',
    '            if not self.host.state_manager.data.atoms and not any(hasattr(i, "create_json_data") for i in self.host.scene.items()):\n                self.host.statusBar().showMessage("Nothing to copy.")'
)
content = rep(
    '                self.statusBar().showMessage("Error: Could not determine molecule bounds for copy.")',
    '                self.host.statusBar().showMessage("Error: Could not determine molecule bounds for copy.")'
)

# patched_copy_2d_image_to_clipboard (on EditActionsManager)
content = rep(
    '                self.statusBar().showMessage("2D Image copied to clipboard (Transparent).", 2000)',
    '                self.host.statusBar().showMessage("2D Image copied to clipboard (Transparent).", 2000)'
)

# ── 16. SVG copy methods: remaining self.scene references ─────────────────────
# These occur in patched_copy_svg_to_clipboard which is also patched onto ExportManager
# selected_items = [i for i in self.scene.selectedItems() if i.isVisible()]
# items_to_export = selected_items if selected_items else list(self.scene.items())
# These may match multiple times — use replace_all=True style
for old, new in [
    ("for i in self.scene.selectedItems()", "for i in self.host.scene.selectedItems()"),
    ("list(self.scene.items())", "list(self.host.scene.items())"),
]:
    count = content.count(old)
    if count:
        content = content.replace(old, new)
        print(f"  [OK x{count}] {old[:60]!r}")
    else:
        print(f"  [MISS] {old[:60]!r}")

# Any remaining self.scene references in export block
for old, new in [
    ("self.scene.backgroundBrush()", "self.host.scene.backgroundBrush()"),
    ("self.scene.setBackgroundBrush(", "self.host.scene.setBackgroundBrush("),
    ("self.scene.clearSelection()", "self.host.scene.clearSelection()"),
    ("self.scene.focusItem()", "self.host.scene.focusItem()"),
    ("self.scene.setFocusItem(None)", "self.host.scene.setFocusItem(None)"),
    ("self.scene.render(", "self.host.scene.render("),
]:
    count = content.count(old)
    if count:
        content = content.replace(old, new)
        print(f"  [OK x{count}] {old[:60]!r} (cleanup pass)")
    # don't warn on 0 since may already be fixed

# ── 17. Remaining self.statusBar() in export block ───────────────────────────
# After all above replacements there should be none left, but catch any stragglers
# Write result FIRST
with open(PATCHER, 'w', encoding='utf-8') as f:
    f.write(content)

remaining_status = content.count("self.statusBar()")
if remaining_status:
    print(f"  [WARN] {remaining_status} self.statusBar() still present")

remaining_scene = content.count("self.scene.")
if remaining_scene:
    print(f"  [WARN] {remaining_scene} self.scene. still present")

changed = content.count('\n') - orig.count('\n')
print(f"Done. Lines delta: {changed:+d}")
