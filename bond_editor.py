import numpy as np
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QPushButton, QMessageBox
)
from PyQt6.QtCore import Qt, QEvent
from rdkit import Chem

PLUGIN_NAME = "Bond Editor"
PLUGIN_VERSION = "2026.03.03"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Edit molecular topology by clicking atoms in the 3D view. Supports adding, deleting, and changing bond orders."

# --- Constants ---
try:
    _pt = Chem.GetPeriodicTable()
except Exception:
    _pt = None


class BondEditorDialog(QDialog):
    """
    A dialog for manually overriding the molecular topology.
    Useful for transition states or hypervalent coordination spheres.
    Click two atoms in the 3D view to select them, then choose a bond order.
    """
    def __init__(self, context, parent=None):
        super().__init__(parent)
        self.context = context
        self.mw = context.get_main_window()
        self.mol = self.mw.current_mol
        self.atom1_idx = None
        self.atom2_idx = None
        self.picking_enabled = False
        self.selection_labels = []
        self._mouse_press_pos = None
        self._mouse_moved = False
        self.setWindowTitle("Bond Editor")
        self.setWindowFlags(Qt.WindowType.Window)

        layout = QVBoxLayout(self)

        # Instructions
        instruction_label = QLabel(
            "Click two atoms in the 3D view to select a bond, then choose the bond order."
        )
        instruction_label.setWordWrap(True)
        layout.addWidget(instruction_label)

        # Selected atoms display
        self.selection_label = QLabel("No atoms selected")
        layout.addWidget(self.selection_label)

        # Bond Order
        order_layout = QHBoxLayout()
        order_layout.addWidget(QLabel("Bond Order:"))
        self.order_combo = QComboBox()
        self.order_combo.addItems(["Single", "Double", "Triple", "Delete (0)"])
        order_layout.addWidget(self.order_combo)
        layout.addLayout(order_layout)

        # Buttons
        btn_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_selection)
        btn_layout.addWidget(self.clear_button)

        btn_layout.addStretch()

        self.apply_button = QPushButton("Apply")
        self.apply_button.clicked.connect(self.apply_changes)
        self.apply_button.setEnabled(False)
        btn_layout.addWidget(self.apply_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        btn_layout.addWidget(close_button)

        layout.addLayout(btn_layout)

        # Pre-populate from existing 3D selection
        self._populate_from_selection()

        # Enable picking
        if self.mol:
            self.enable_picking()

    # ------------------------------------------------------------------
    # Picking (self-contained, no mixin dependency)
    # ------------------------------------------------------------------

    def eventFilter(self, obj, event):
        """Capture mouse clicks on the 3D view for atom picking."""
        if (
            obj == self.mw.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonPress
            and event.button() == Qt.MouseButton.LeftButton
        ):
            self._mouse_press_pos = event.pos()
            self._mouse_moved = False

            try:
                interactor = self.mw.plotter.interactor
                click_pos = interactor.GetEventPosition()
                picker = self.mw.plotter.picker
                picker.Pick(click_pos[0], click_pos[1], 0, self.mw.plotter.renderer)

                if picker.GetActor() is self.mw.atom_actor:
                    picked_position = np.array(picker.GetPickPosition())
                    distances = np.linalg.norm(
                        self.mw.atom_positions_3d - picked_position, axis=1
                    )
                    closest_atom_idx = np.argmin(distances)

                    if 0 <= closest_atom_idx < self.mol.GetNumAtoms():
                        atom = self.mol.GetAtomWithIdx(int(closest_atom_idx))
                        if atom:
                            try:
                                atomic_num = atom.GetAtomicNum()
                                vdw_radius = _pt.GetRvdw(atomic_num) if _pt else 1.5
                                if vdw_radius < 0.1:
                                    vdw_radius = 1.5
                            except Exception:
                                vdw_radius = 1.5
                            click_threshold = vdw_radius * 1.5

                            if distances[closest_atom_idx] < click_threshold:
                                try:
                                    self.mw._picking_consumed = True
                                except Exception:
                                    pass
                                self.on_atom_picked(int(closest_atom_idx))
                                self._mouse_press_pos = None
                                return True

                return False
            except Exception:
                return False

        elif (
            obj == self.mw.plotter.interactor
            and event.type() == QEvent.Type.MouseMove
        ):
            if self._mouse_press_pos is not None:
                diff = event.pos() - self._mouse_press_pos
                if diff.manhattanLength() > 3:
                    self._mouse_moved = True

        elif (
            obj == self.mw.plotter.interactor
            and event.type() == QEvent.Type.MouseButtonRelease
            and event.button() == Qt.MouseButton.LeftButton
        ):
            if self._mouse_press_pos is not None:
                if not self._mouse_moved:
                    self.clear_selection()
                self._mouse_press_pos = None
                self._mouse_moved = False

        return super().eventFilter(obj, event)

    def enable_picking(self):
        self.mw.plotter.interactor.installEventFilter(self)
        self.picking_enabled = True
        try:
            self.mw._picking_consumed = False
        except Exception:
            pass

    def disable_picking(self):
        if self.picking_enabled:
            self.mw.plotter.interactor.removeEventFilter(self)
            self.picking_enabled = False
        try:
            if hasattr(self.mw, "_picking_consumed"):
                self.mw._picking_consumed = False
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Label management
    # ------------------------------------------------------------------

    def clear_atom_labels(self):
        for label_actor in self.selection_labels:
            try:
                self.mw.plotter.remove_actor(label_actor)
            except Exception:
                pass
        self.selection_labels = []

    def show_atom_labels(self):
        self.clear_atom_labels()
        if self.mol is None:
            return

        atoms_and_labels = []
        if self.atom1_idx is not None and self.atom1_idx < self.mol.GetNumAtoms():
            atoms_and_labels.append((self.atom1_idx, "1"))
        if self.atom2_idx is not None and self.atom2_idx < self.mol.GetNumAtoms():
            atoms_and_labels.append((self.atom2_idx, "2"))

        for atom_idx, label_text in atoms_and_labels:
            pos = self.mw.atom_positions_3d[atom_idx]
            label_actor = self.mw.plotter.add_point_labels(
                [pos],
                [label_text],
                point_size=20,
                font_size=12,
                text_color="yellow",
                always_visible=True,
            )
            self.selection_labels.append(label_actor)

    # ------------------------------------------------------------------
    # Dialog logic
    # ------------------------------------------------------------------

    def _populate_from_selection(self):
        """Pre-populate from existing 3D selection if available."""
        if hasattr(self.mw, "selected_atoms_3d") and len(self.mw.selected_atoms_3d) == 2:
            sel_list = list(self.mw.selected_atoms_3d)
            if self.mol:
                self.atom1_idx = sel_list[0]
                self.atom2_idx = sel_list[1]
                self.update_display()
                self.show_atom_labels()

    def on_atom_picked(self, atom_idx):
        if self.mol is None:
            return

        if self.atom1_idx is None:
            self.atom1_idx = atom_idx
        elif self.atom2_idx is None:
            self.atom2_idx = atom_idx
        else:
            self.atom1_idx = atom_idx
            self.atom2_idx = None

        self.show_atom_labels()
        self.update_display()

    def update_display(self):
        labels = []
        if self.atom1_idx is not None and self.atom1_idx < self.mol.GetNumAtoms():
            try:
                a1 = self.mol.GetAtomWithIdx(self.atom1_idx)
                sym = a1.GetSymbol()
                try:
                    orig_id = a1.GetIntProp('_original_atom_id')
                    labels.append(f"Atom 1: {sym} (ID {orig_id + 1})")
                except Exception:
                    labels.append(f"Atom 1: {sym} (idx {self.atom1_idx})")
            except Exception:
                labels.append("Atom 1: ?")

        if self.atom2_idx is not None and self.atom2_idx < self.mol.GetNumAtoms():
            try:
                a2 = self.mol.GetAtomWithIdx(self.atom2_idx)
                sym = a2.GetSymbol()
                try:
                    orig_id = a2.GetIntProp('_original_atom_id')
                    labels.append(f"Atom 2: {sym} (ID {orig_id + 1})")
                except Exception:
                    labels.append(f"Atom 2: {sym} (idx {self.atom2_idx})")
            except Exception:
                labels.append("Atom 2: ?")

        if labels:
            self.selection_label.setText("Selected: " + ", ".join(labels))
        else:
            self.selection_label.setText("No atoms selected")

        self.apply_button.setEnabled(
            self.atom1_idx is not None and self.atom2_idx is not None
        )

    def clear_selection(self):
        self.atom1_idx = None
        self.atom2_idx = None
        self.update_display()
        self.clear_atom_labels()

    def apply_changes(self):
        if self.atom1_idx is None or self.atom2_idx is None:
            QMessageBox.warning(self, "Invalid Input", "Please select two atoms first.")
            return

        if self.atom1_idx == self.atom2_idx:
            QMessageBox.warning(self, "Invalid Input", "Atom 1 and Atom 2 cannot be the same.")
            return

        # Get original 2D atom IDs for the data model
        try:
            atom1 = self.mol.GetAtomWithIdx(self.atom1_idx)
            atom2 = self.mol.GetAtomWithIdx(self.atom2_idx)
            id1 = atom1.GetIntProp("_original_atom_id")
            id2 = atom2.GetIntProp("_original_atom_id")
        except Exception:
            QMessageBox.warning(
                self, "Error",
                "Could not resolve original atom IDs. "
                "Hydrogen atoms added during 3D conversion cannot be used for bond editing."
            )
            return

        order_idx = self.order_combo.currentIndex()
        order_map = [1.0, 2.0, 3.0, 0.0]
        order = order_map[order_idx]

        try:
            if hasattr(self.mw, "apply_bond_editor_changes"):
                self.mw.apply_bond_editor_changes(id1, id2, order)
            else:
                # Direct fallback
                self.mw.data.force_bond(id1, id2, order)
                self.mw.scene.reinitialize_items()
                self.mw.scene.update()
                try:
                    self.mw.trigger_conversion()
                except Exception as e2:
                    self.mw.statusBar().showMessage(f"Failed to generate 3D: {e2}")
                self.mw.push_undo_state()
                self.mw.statusBar().showMessage("Applied manual bond edit.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to apply bond edit: {e}")

    def reject(self):
        self.disable_picking()
        self.clear_atom_labels()
        super().reject()
        try:
            if self.mw.current_mol:
                self.mw.draw_molecule_3d(self.mw.current_mol)
        except Exception:
            pass


def initialize(context):
    mw = context.get_main_window()

    def show_bond_editor():
        if not hasattr(mw, '_bond_editor_dialog') or mw._bond_editor_dialog is None:
            mw._bond_editor_dialog = BondEditorDialog(context, parent=mw)

        # Refresh mol reference
        mw._bond_editor_dialog.mol = mw.current_mol

        mw._bond_editor_dialog.show()
        mw._bond_editor_dialog.raise_()
        mw._bond_editor_dialog.activateWindow()

    context.add_menu_action("3D Edit/Edit Bond...", show_bond_editor)


def run(mw):
    """Legacy support."""
    class LegacyContext:
        def __init__(self, mw):
            self._mw = mw
        def get_main_window(self):
            return self._mw
        @property
        def current_molecule(self):
            return self._mw.current_mol
        @current_molecule.setter
        def current_molecule(self, mol):
            self._mw.current_mol = mol
            self._mw.draw_molecule_3d(mol)

    ctx = LegacyContext(mw)
    if not hasattr(mw, '_bond_editor_dialog') or mw._bond_editor_dialog is None:
        mw._bond_editor_dialog = BondEditorDialog(ctx, parent=mw)

    mw._bond_editor_dialog.mol = mw.current_mol
    mw._bond_editor_dialog.show()
    mw._bond_editor_dialog.raise_()
    mw._bond_editor_dialog.activateWindow()
