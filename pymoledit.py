import sys
import numpy as np
import pickle
import copy 

# PyQt6 Modules
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QPushButton, QSplitter, QGraphicsView, QGraphicsScene, QGraphicsItem,
    QToolBar, QStatusBar, QGraphicsTextItem, QGraphicsLineItem, QDialog, QGridLayout,
    QFileDialog
)
from PyQt6.QtGui import (
    QPen, QBrush, QColor, QPainter, QAction, QActionGroup, QFont, QPolygonF,
    QPainterPath, QFontMetrics, QKeySequence
)
from PyQt6.QtCore import Qt, QPointF, QRectF, QLineF, QObject, QThread, pyqtSignal

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# PyVista
import pyvista as pv
from pyvistaqt import QtInteractor

# --- Constants ---
ATOM_RADIUS = 18
BOND_OFFSET = 3.5
CPK_COLORS = {
    'H': QColor("white"), 'C': QColor("#333333"), 'N': QColor("#0000FF"),
    'O': QColor("#FF0000"), 'F': QColor("#00FF00"), 'S': QColor("#FFC000"),
    'Cl': QColor("#00FF00"), 'Br': QColor("#A52A2A"), 'I': QColor("#9400D3"),
    'P': QColor("#FFA500"), 'Si': QColor("#DAA520"), 'B': QColor("#FA8072"),
}
CPK_COLORS_PV = {
    'H': [0.9, 0.9, 0.9], 'C': [0.2, 0.2, 0.2], 'N': [0.0, 0.0, 1.0],
    'O': [1.0, 0.0, 0.0], 'F': [0.0, 1.0, 0.0], 'S': [1.0, 0.8, 0.0],
    'Cl': [0.0, 1.0, 0.0], 'Br': [0.6, 0.2, 0.0], 'I': [0.4, 0.0, 0.6],
    'P': [1.0, 0.65, 0.0], 'Si': [0.85, 0.65, 0.125], 'B': [0.98, 0.5, 0.45],
}
pt = Chem.GetPeriodicTable()
VDW_RADII = {pt.GetElementSymbol(i): pt.GetRvdw(i) * 0.3 for i in range(1, 119)}


# --- Data Model ---
class MolecularData:
    def __init__(self):
        self.atoms = {}; self.bonds = {}; self._next_atom_id = 0
    def add_atom(self, symbol, pos):
        atom_id = self._next_atom_id
        self.atoms[atom_id] = {'symbol': symbol, 'pos': pos, 'item': None}
        self._next_atom_id += 1; return atom_id
    def add_bond(self, id1, id2, order=1):
        if id1 > id2: id1, id2 = id2, id1
        if (id1, id2) in self.bonds:
            self.bonds[(id1, id2)]['order'] = order; return (id1, id2), 'updated'
        else:
            self.bonds[(id1, id2)] = {'order': order, 'item': None}; return (id1, id2), 'created'
    def remove_atom(self, atom_id):
        if atom_id in self.atoms:
            del self.atoms[atom_id]
            bonds_to_remove = [key for key in self.bonds if atom_id in key]
            for key in bonds_to_remove: del self.bonds[key]
    def remove_bond(self, id1, id2):
        if id1 > id2: id1, id2 = id2, id1
        if (id1, id2) in self.bonds: del self.bonds[(id1, id2)]
    def to_mol_block(self):
        if not self.atoms: return None
        atom_map = {old_id: new_id for new_id, old_id in enumerate(self.atoms.keys())}
        num_atoms, num_bonds = len(self.atoms), len(self.bonds)
        mol_block = "\n  PyQtEditor\n\n"
        mol_block += f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n"
        for old_id, atom in self.atoms.items():
            x, y, z, symbol = atom['item'].pos().x(), -atom['item'].pos().y(), 0.0, atom['symbol']
            mol_block += f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n"
        for (id1, id2), bond in self.bonds.items():
            idx1, idx2, order = atom_map[id1] + 1, atom_map[id2] + 1, bond['order']
            mol_block += f"{idx1:3d}{idx2:3d}{order:3d}  0  0  0  0\n"
        mol_block += "M  END\n"
        return mol_block


# --- Custom 2D Graphics Items ---
class AtomItem(QGraphicsItem):
    def __init__(self, atom_id, symbol, pos):
        super().__init__()
        self.atom_id, self.symbol, self.bonds = atom_id, symbol, []
        self.setPos(pos)
        self.setFlags(QGraphicsItem.GraphicsItemFlag.ItemIsMovable | QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(1); self.font = QFont("Arial", 16, QFont.Weight.Bold); self.update_style()
    def boundingRect(self): return QRectF(-ATOM_RADIUS, -ATOM_RADIUS, ATOM_RADIUS*2, ATOM_RADIUS*2)
    def paint(self, painter, option, widget):
        if self.is_visible:
            painter.setFont(self.font); fm = painter.fontMetrics()
            text_rect = fm.boundingRect(self.symbol); text_rect.moveCenter(QPointF(0, 0).toPoint())
            if self.scene():
                bg_brush = self.scene().backgroundBrush(); bg_rect = text_rect.adjusted(-3, -3, 3, 3)
                painter.setBrush(bg_brush); painter.setPen(Qt.PenStyle.NoPen); painter.drawRect(bg_rect)
            if self.symbol == 'H':
                path = QPainterPath(); path.addText(text_rect.left(), text_rect.bottom(), self.font, self.symbol)
                painter.setPen(QPen(Qt.GlobalColor.black, 1)); painter.setBrush(QBrush(CPK_COLORS.get(self.symbol)))
                painter.drawPath(path)
            else:
                painter.setPen(QPen(CPK_COLORS.get(self.symbol, QColor("pink"))))
                painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, self.symbol)
        if self.isSelected():
            painter.setBrush(Qt.BrushStyle.NoBrush); painter.setPen(QPen(QColor(0, 100, 255), 3))
            painter.drawRect(self.boundingRect())
    def update_style(self):
        self.is_visible = not (self.symbol == 'C' and len(self.bonds) > 0); self.update()
    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            for bond in self.bonds: bond.update_position()
        return super().itemChange(change, value)

class BondItem(QGraphicsItem):
    def __init__(self, atom1_item, atom2_item, order=1):
        super().__init__()
        self.atom1, self.atom2, self.order = atom1_item, atom2_item, order
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable); self.pen = QPen(Qt.GlobalColor.black, 2)
    def boundingRect(self):
        extra = (self.order - 1) * BOND_OFFSET + 5; return self.shape().controlPointRect().adjusted(-extra, -extra, extra, extra)
    def shape(self):
        path = QPainterPath(); line = QLineF(self.mapFromItem(self.atom1,0,0), self.mapFromItem(self.atom2,0,0))
        length = line.length()
        if length == 0: perp_dx, perp_dy = 0, 5
        else: dx, dy = line.dx(), line.dy(); perp_dx, perp_dy = -dy/length*5, dx/length*5
        offset = QPointF(perp_dx, perp_dy)
        p = QPolygonF([line.p1()-offset, line.p1()+offset, line.p2()+offset, line.p2()-offset]); path.addPolygon(p)
        return path
    def paint(self, painter, option, widget):
        painter.setPen(self.pen); line = QLineF(self.mapFromItem(self.atom1,0,0), self.mapFromItem(self.atom2,0,0))
        if self.isSelected(): painter.setPen(QPen(QColor("blue"), 3))
        if self.order == 1: painter.drawLine(line)
        else:
            v = line.unitVector().normalVector(); offset = QPointF(v.dx(), v.dy()) * BOND_OFFSET
            if self.order == 2: painter.drawLine(line.translated(offset)); painter.drawLine(line.translated(-offset))
            elif self.order == 3: painter.drawLine(line); painter.drawLine(line.translated(offset)); painter.drawLine(line.translated(-offset))
    def update_position(self): self.prepareGeometryChange(); self.update()


# --- 2D Editor Scene Class ---
class MoleculeScene(QGraphicsScene):
    def __init__(self, data, window):
        super().__init__()
        self.data, self.window = data, window
        self.mode, self.current_atom_symbol, self.bond_order = 'select', 'C', 1
        self.start_atom, self.temp_line, self.start_pos = None, None, None; self.press_pos = None
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False

    def mousePressEvent(self, event):
        self.press_pos = event.scenePos()
        self.mouse_moved_since_press = False
        self.data_changed_in_event = False
        self.initial_positions_in_event = {item: item.pos() for item in self.items() if isinstance(item, AtomItem)}
        
        item = self.itemAt(self.press_pos, self.views()[0].transform())
        if isinstance(item, AtomItem):
            self.start_atom = item
            if self.mode != 'select':
                self.temp_line=QGraphicsLineItem(QLineF(self.start_atom.pos(),self.press_pos)); self.temp_line.setPen(QPen(Qt.GlobalColor.red,2,Qt.PenStyle.DotLine)); self.addItem(self.temp_line)
            else: super().mousePressEvent(event)
        elif item is None and self.mode.startswith('atom'):
            self.start_pos = self.press_pos
            self.temp_line = QGraphicsLineItem(QLineF(self.start_pos, self.press_pos))
            self.temp_line.setPen(QPen(Qt.GlobalColor.red, 2, Qt.PenStyle.DotLine))
            self.addItem(self.temp_line)
        else: super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.mouse_moved_since_press and self.press_pos:
            if (event.scenePos() - self.press_pos).manhattanLength() > QApplication.startDragDistance():
                self.mouse_moved_since_press = True
        
        if self.temp_line:
            start_point = self.start_atom.pos() if self.start_atom else self.start_pos
            if start_point:
                self.temp_line.setLine(QLineF(start_point, event.scenePos()))
        else:
            super().mouseMoveEvent(event)


    def mouseReleaseEvent(self, event):
        end_pos = event.scenePos()
        is_click = not self.mouse_moved_since_press
        released_item = self.itemAt(end_pos, self.views()[0].transform())
        
        if self.mode.startswith('atom') and is_click and isinstance(released_item, BondItem):
            b=released_item; new_order=(b.order%3)+1; id1,id2=b.atom1.atom_id,b.atom2.atom_id
            if id1>id2: id1,id2=id2,id1
            self.data.bonds[(id1,id2)]['order']=new_order; b.order=new_order; b.update()
            self.data_changed_in_event = True
            if self.temp_line: self.removeItem(self.temp_line); self.temp_line=None
        elif self.start_atom and self.temp_line: # Drag started from an existing atom
            self.removeItem(self.temp_line); self.temp_line = None
            line = QLineF(self.start_atom.pos(), end_pos); end_item = self.itemAt(end_pos, self.views()[0].transform())
            if line.length() < 10:
                if self.mode.startswith('atom'):
                    if self.start_atom.symbol != self.current_atom_symbol:
                        self.start_atom.symbol=self.current_atom_symbol; self.data.atoms[self.start_atom.atom_id]['symbol']=self.current_atom_symbol; self.start_atom.update_style()
                        self.data_changed_in_event = True
                elif self.mode == 'select': self.start_atom.setSelected(True)
            else:
                if isinstance(end_item, AtomItem) and self.start_atom!=end_item: self.create_bond(self.start_atom, end_item)
                else:
                    new_id = self.create_atom(self.current_atom_symbol, end_pos); new_item = self.data.atoms[new_id]['item']
                    self.create_bond(self.start_atom, new_item)
                self.data_changed_in_event = True
        elif self.start_pos: # Drag started from an empty space
            if self.temp_line:
                self.removeItem(self.temp_line)
                self.temp_line = None
            
            line = QLineF(self.start_pos, end_pos)
            if line.length() < 10: # A click
                self.create_atom(self.current_atom_symbol, end_pos)
                self.data_changed_in_event = True
            else: # A drag
                start_id = self.create_atom(self.current_atom_symbol, self.start_pos)
                end_id = self.create_atom(self.current_atom_symbol, end_pos)
                start_item = self.data.atoms[start_id]['item']
                end_item = self.data.atoms[end_id]['item']
                
                original_order = self.bond_order
                self.bond_order = 1
                self.create_bond(start_item, end_item)
                self.bond_order = original_order
                self.data_changed_in_event = True
        else:
            super().mouseReleaseEvent(event)
        
        for item, old_pos in self.initial_positions_in_event.items():
            if item in self.items() and item.pos() != old_pos:
                self.data_changed_in_event = True
                break

        self.start_atom=None; self.start_pos = None; self.press_pos = None
        
        if self.data_changed_in_event:
            self.window.push_undo_state()
            
    def create_atom(self, symbol, pos):
        atom_id = self.data.add_atom(symbol, pos); atom_item = AtomItem(atom_id, symbol, pos)
        self.data.atoms[atom_id]['item'] = atom_item; self.addItem(atom_item); return atom_id

    def create_bond(self, start_atom, end_atom):
        key, status = self.data.add_bond(start_atom.atom_id, end_atom.atom_id, self.bond_order)
        if status == 'created':
            bond_item=BondItem(start_atom, end_atom, self.bond_order); self.data.bonds[key]['item']=bond_item
            start_atom.bonds.append(bond_item); end_atom.bonds.append(bond_item); self.addItem(bond_item)
        else: self.data.bonds[key]['item'].order=self.bond_order; self.data.bonds[key]['item'].update()
        start_atom.update_style(); end_atom.update_style()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Delete or event.key() == Qt.Key.Key_Backspace:
            selected_items = self.selectedItems()
            if not selected_items:
                return

            atoms_to_delete = {item for item in selected_items if isinstance(item, AtomItem)}
            bonds_to_delete = {item for item in selected_items if isinstance(item, BondItem)}

            for atom in atoms_to_delete:
                for bond in atom.bonds[:]:
                    bonds_to_delete.add(bond)
            
            atoms_to_update = set()

            for bond in bonds_to_delete:
                atom1, atom2 = bond.atom1, bond.atom2
                
                if atom1 not in atoms_to_delete:
                    atoms_to_update.add(atom1)
                if atom2 not in atoms_to_delete:
                    atoms_to_update.add(atom2)

                if bond in atom1.bonds:
                    atom1.bonds.remove(bond)
                if bond in atom2.bonds:
                    atom2.bonds.remove(bond)

                self.data.remove_bond(atom1.atom_id, atom2.atom_id)
                self.removeItem(bond)

            for atom in atoms_to_delete:
                self.data.remove_atom(atom.atom_id)
                self.removeItem(atom)

            for atom in atoms_to_update:
                atom.update_style()
            
            self.window.push_undo_state()
        else:
            super().keyPressEvent(event)


# --- Worker Thread for Calculations ---
class CalculationWorker(QObject):
    finished=pyqtSignal(object); error=pyqtSignal(str)
    def run_calculation(self, mol_block):
        try:
            if not mol_block: raise ValueError("No atoms to convert.")
            mol=Chem.MolFromMolBlock(mol_block, removeHs=False)
            if mol is None: raise ValueError("Failed to create molecule from MOL block.")
            mol=Chem.AddHs(mol); AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.MMFFOptimizeMolecule(mol)
            self.finished.emit(mol)
        except Exception as e: self.error.emit(str(e))


# --- Periodic Table Dialog ---
class PeriodicTableDialog(QDialog):
    element_selected = pyqtSignal(str)
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select an Element"); layout = QGridLayout(self); self.setLayout(layout)
        elements=[('H',1,1),('He',1,18),('Li',2,1),('Be',2,2),('B',2,13),('C',2,14),('N',2,15),('O',2,16),('F',2,17),('Ne',2,18),('Na',3,1),('Mg',3,2),('Al',3,13),('Si',3,14),('P',3,15),('S',3,16),('Cl',3,17),('Ar',3,18),('K',4,1),('Ca',4,2),('Sc',4,3),('Ti',4,4),('V',4,5),('Cr',4,6),('Mn',4,7),('Fe',4,8),('Co',4,9),('Ni',4,10),('Cu',4,11),('Zn',4,12),('Ga',4,13),('Ge',4,14),('As',4,15),('Se',4,16),('Br',4,17),('Kr',4,18),('Rb',5,1),('Sr',5,2),('Y',5,3),('Zr',5,4),('Nb',5,5),('Mo',5,6),('Tc',5,7),('Ru',5,8),('Rh',5,9),('Pd',5,10),('Ag',5,11),('Cd',5,12),('In',5,13),('Sn',5,14),('Sb',5,15),('Te',5,16),('I',5,17),('Xe',5,18)]
        for symbol, row, col in elements:
            b=QPushButton(symbol); b.setFixedSize(40,40); b.clicked.connect(self.on_button_clicked); layout.addWidget(b, row, col)
    def on_button_clicked(self):
        b=self.sender(); self.element_selected.emit(b.text()); self.accept()


# --- Main Window ---
class MainWindow(QMainWindow):
    start_calculation = pyqtSignal(str)
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Python Molecular Editor"); self.setGeometry(100, 100, 1600, 900)
        self.data = MolecularData(); self.current_mol = None
        self.undo_stack = []
        self.redo_stack = []
        self.init_ui()
        self.init_worker_thread()
        self.reset_undo_stack()

    def init_ui(self):
        self.init_menu_bar()

        splitter=QSplitter(Qt.Orientation.Horizontal)
        self.setCentralWidget(splitter)

        left_pane=QWidget()
        left_layout=QVBoxLayout(left_pane)
        left_layout.setContentsMargins(0,0,0,0)

        self.scene=MoleculeScene(self.data,self)
        self.scene.setSceneRect(-800,-800,1600,1600)
        self.scene.setBackgroundBrush(QColor("#FFFFFF"))

        self.view_2d=QGraphicsView(self.scene)
        self.view_2d.setRenderHint(QPainter.RenderHint.Antialiasing)
        left_layout.addWidget(self.view_2d)

        self.cleanup_button=QPushButton("Optimize 2D")
        self.cleanup_button.clicked.connect(self.clean_up_2d_structure)
        left_layout.addWidget(self.cleanup_button)
        splitter.addWidget(left_pane)

        right_pane=QWidget()
        right_layout=QVBoxLayout(right_pane)
        self.plotter=QtInteractor(right_pane)
        right_layout.addWidget(self.plotter)

        self.convert_button=QPushButton("Convert to 3D")
        self.convert_button.clicked.connect(self.trigger_conversion)
        right_layout.addWidget(self.convert_button)
        splitter.addWidget(right_pane)
        splitter.setSizes([800, 800])

        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)
        self.tool_group = QActionGroup(self)
        self.tool_group.setExclusive(True)

        actions = {"Select":'select', "C":'atom_C', "N":'atom_N', "O":'atom_O', 
                   "S":'atom_S', "F":'atom_F', "Cl":'atom_Cl', "Br":'atom_Br', 
                   "I":'atom_I', "H":'atom_H', "Other...":'atom_other'}
        for text, mode in actions.items():
            if text == "C":
                toolbar.addSeparator()
            action=QAction(text,self,checkable=(mode!='atom_other'))
            if mode=='atom_other':
                action.triggered.connect(self.open_periodic_table_dialog)
            else:
                action.triggered.connect(lambda c,m=mode: self.set_mode(m))
            toolbar.addAction(action)
            if mode!='atom_other':
                self.tool_group.addAction(action)
            if text=="C":
                action.setChecked(True)
                self.set_mode('atom_C')

    def init_menu_bar(self):
        menu_bar = self.menuBar()
        
        file_menu = menu_bar.addMenu("&File")
        load_mol_action = QAction("Open MOL/SDF...", self); load_mol_action.triggered.connect(self.load_mol_file)
        file_menu.addAction(load_mol_action)
        file_menu.addSeparator()
        save_mol_action = QAction("Save 2D as MOL...", self); save_mol_action.triggered.connect(self.save_as_mol)
        file_menu.addAction(save_mol_action)
        save_xyz_action = QAction("Save 3D as XYZ...", self); save_xyz_action.triggered.connect(self.save_as_xyz)
        file_menu.addAction(save_xyz_action)
        file_menu.addSeparator()
        save_raw_action = QAction("Save Project...", self); save_raw_action.triggered.connect(self.save_raw_data)
        file_menu.addAction(save_raw_action)
        load_raw_action = QAction("Open Project...", self); load_raw_action.triggered.connect(self.load_raw_data)
        file_menu.addAction(load_raw_action)
        
        edit_menu = menu_bar.addMenu("&Edit")
        self.undo_action = QAction("Undo", self); self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.triggered.connect(self.undo); edit_menu.addAction(self.undo_action)
        
        self.redo_action = QAction("Redo", self); self.redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.redo_action.triggered.connect(self.redo); edit_menu.addAction(self.redo_action)
        
        edit_menu.addSeparator()
        
        select_all_action = QAction("Select All", self); select_all_action.setShortcut(QKeySequence.StandardKey.SelectAll)
        select_all_action.triggered.connect(self.select_all); edit_menu.addAction(select_all_action)
        
        clear_all_action = QAction("Clear All", self)
        clear_all_action.triggered.connect(self.clear_all); edit_menu.addAction(clear_all_action)

    # --- Undo/Redo and State Management ---
    def get_current_state(self):
        atoms = {atom_id: {'symbol': data['symbol'], 'pos': (data['item'].pos().x(), data['item'].pos().y())}
                 for atom_id, data in self.data.atoms.items()}
        bonds = {key: {'order': data['order']} for key, data in self.data.bonds.items()}
        state = {'atoms': atoms, 'bonds': bonds, '_next_atom_id': self.data._next_atom_id}
        
        if self.current_mol:
            state['mol_3d'] = self.current_mol.ToBinary()
            
        return state

    def set_state_from_data(self, state_data):
        self.clear_2d_editor(push_to_undo=False)
        
        loaded_data = state_data
        raw_atoms = loaded_data.get('atoms', {})
        raw_bonds = loaded_data.get('bonds', {})

        for atom_id, data in raw_atoms.items():
            pos = QPointF(data['pos'][0], data['pos'][1])
            atom_item = AtomItem(int(atom_id), data['symbol'], pos)
            self.data.atoms[int(atom_id)] = {'symbol': data['symbol'], 'pos': pos, 'item': atom_item}
            self.scene.addItem(atom_item)
        
        self.data._next_atom_id = loaded_data.get('_next_atom_id', max(self.data.atoms.keys()) + 1 if self.data.atoms else 0)

        for (id1, id2), data in raw_bonds.items():
            if id1 in self.data.atoms and id2 in self.data.atoms:
                atom1_item = self.data.atoms[id1]['item']
                atom2_item = self.data.atoms[id2]['item']
                self.scene.create_bond(atom1_item, atom2_item)
                key = (id1, id2) if id1 < id2 else (id2, id1)
                self.data.bonds[key]['order'] = data['order']
                if self.data.bonds[key]['item']:
                    self.data.bonds[key]['item'].order = data['order']
                    self.data.bonds[key]['item'].update()
        
        for atom_data in self.data.atoms.values():
            if atom_data['item']:
                atom_data['item'].update_style()
        self.scene.update()

        if 'mol_3d' in state_data:
            try:
                self.current_mol = Chem.Mol(state_data['mol_3d'])
                self.draw_molecule_3d(self.current_mol)
            except Exception as e:
                self.statusBar().showMessage(f"Could not load 3D model from project: {e}")
                self.current_mol = None
        else:
            self.current_mol = None
            self.plotter.clear()

    def push_undo_state(self):
        state = self.get_current_state()
        if self.undo_stack and self.undo_stack[-1] == state:
            return
        self.undo_stack.append(state)
        self.redo_stack.clear()
        self.update_undo_redo_actions()

    def reset_undo_stack(self):
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.push_undo_state()

    def undo(self):
        if len(self.undo_stack) > 1:
            self.redo_stack.append(self.undo_stack.pop())
            state = self.undo_stack[-1]
            self.set_state_from_data(state)
        self.update_undo_redo_actions()

    def redo(self):
        if self.redo_stack:
            state = self.redo_stack.pop()
            self.undo_stack.append(state)
            self.set_state_from_data(state)
        self.update_undo_redo_actions()
        
    def update_undo_redo_actions(self):
        self.undo_action.setEnabled(len(self.undo_stack) > 1)
        self.redo_action.setEnabled(len(self.redo_stack) > 0)

    # --- Edit Menu Actions ---
    def select_all(self):
        for item in self.scene.items():
            item.setSelected(True)

    def clear_all(self):
        if not self.data.atoms: 
            if self.current_mol is None: return
        self.clear_2d_editor(push_to_undo=True)
        self.current_mol = None
        self.plotter.clear()
        self.push_undo_state()
        
    def clear_2d_editor(self, push_to_undo=True):
        self.data = MolecularData()
        self.scene.data = self.data
        self.scene.clear()
        if push_to_undo:
            self.push_undo_state()

    # --- File I/O Methods ---
    def load_mol_file(self):
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Open MOL/SDF File", "", "Chemical Files (*.mol *.sdf);;All Files (*)", options=options)
        if not file_path: return
        try:
            mol = Chem.MolFromMolFile(file_path, removeHs=False)
            if mol is None: raise ValueError("Failed to read molecule from file.")
            self.clear_2d_editor(push_to_undo=False)
            self.current_mol = None
            self.plotter.clear()
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer(); SCALE_FACTOR = 50.0
            rdkit_idx_to_my_id = {}
            positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            center_x = sum(p.x for p in positions)/len(positions) if positions else 0
            center_y = sum(p.y for p in positions)/len(positions) if positions else 0
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i); pos = conf.GetAtomPosition(i)
                scene_x=(pos.x-center_x)*SCALE_FACTOR; scene_y=-(pos.y-center_y)*SCALE_FACTOR
                atom_id = self.scene.create_atom(atom.GetSymbol(), QPointF(scene_x, scene_y))
                rdkit_idx_to_my_id[i] = atom_id
            for bond in mol.GetBonds():
                b_idx,e_idx,b_type=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx(),bond.GetBondTypeAsDouble()
                a1_id, a2_id = rdkit_idx_to_my_id[b_idx], rdkit_idx_to_my_id[e_idx]
                a1_item,a2_item=self.data.atoms[a1_id]['item'],self.data.atoms[a2_id]['item']
                orig_order=self.scene.bond_order; self.scene.bond_order=int(b_type)
                self.scene.create_bond(a1_item, a2_item); self.scene.bond_order=orig_order
            self.statusBar().showMessage(f"Successfully loaded {file_path}")
            self.reset_undo_stack()
        except Exception as e: self.statusBar().showMessage(f"Error loading file: {e}")

    def save_raw_data(self):
        if not self.data.atoms: self.statusBar().showMessage("Error: Nothing to save."); return
        save_data = self.get_current_state()
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
        if file_path:
            if not file_path.lower().endswith('.pmeraw'): file_path += '.pmeraw'
            try:
                with open(file_path, 'wb') as f: pickle.dump(save_data, f)
                self.statusBar().showMessage(f"Project saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving project file: {e}")

    def load_raw_data(self):
        options = QFileDialog.Option.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Project File", "", "Project Files (*.pmeraw);;All Files (*)", options=options)
        if not file_path: return
        try:
            with open(file_path, 'rb') as f: loaded_data = pickle.load(f)
            self.set_state_from_data(loaded_data)
            self.statusBar().showMessage(f"Project loaded from {file_path}")
            self.reset_undo_stack()
        except Exception as e: self.statusBar().showMessage(f"Error loading project file: {e}")

    def save_as_mol(self):
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No 2D data to save."); return
        options=QFileDialog.Option.DontUseNativeDialog
        file_path,_=QFileDialog.getSaveFileName(self,"Save 2D MOL File","","MOL Files (*.mol);;All Files (*)",options=options)
        if file_path:
            if not file_path.lower().endswith('.mol'): file_path += '.mol'
            try:
                with open(file_path,'w') as f: f.write(mol_block)
                self.statusBar().showMessage(f"2D data saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")

    def save_as_xyz(self):
        if not self.current_mol: self.statusBar().showMessage("Error: Please generate a 3D structure first."); return
        options=QFileDialog.Option.DontUseNativeDialog
        file_path,_=QFileDialog.getSaveFileName(self,"Save 3D XYZ File","","XYZ Files (*.xyz);;All Files (*)",options=options)
        if file_path:
            if not file_path.lower().endswith('.xyz'): file_path += '.xyz'
            try:
                conf=self.current_mol.GetConformer(); num_atoms=self.current_mol.GetNumAtoms()
                xyz_lines=[str(num_atoms)]; smiles=Chem.MolToSmiles(Chem.RemoveHs(self.current_mol))
                xyz_lines.append(f"Generated by Python Molecular Editor. SMILES: {smiles}")
                for i in range(num_atoms):
                    pos=conf.GetAtomPosition(i); symbol=self.current_mol.GetAtomWithIdx(i).GetSymbol()
                    xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
                with open(file_path,'w') as f: f.write("\n".join(xyz_lines))
                self.statusBar().showMessage(f"Successfully saved to {file_path}")
            except Exception as e: self.statusBar().showMessage(f"Error saving file: {e}")

    # --- Worker and Mode Management ---
    def init_worker_thread(self):
        self.thread=QThread();self.worker=CalculationWorker();self.worker.moveToThread(self.thread)
        self.start_calculation.connect(self.worker.run_calculation)
        self.worker.finished.connect(self.on_calculation_finished); self.worker.error.connect(self.on_calculation_error)
        self.thread.start()
    
    def open_periodic_table_dialog(self):
        dialog=PeriodicTableDialog(self); dialog.element_selected.connect(self.set_atom_from_periodic_table)
        checked_action=self.tool_group.checkedAction()
        if checked_action: self.tool_group.setExclusive(False); checked_action.setChecked(False); self.tool_group.setExclusive(True)
        dialog.exec()

    def set_atom_from_periodic_table(self, symbol): self.set_mode(f'atom_{symbol}')

    def set_mode(self, mode_str):
        parts=mode_str.split('_')
        self.scene.mode=parts[0]
        if self.scene.mode=='atom': 
            self.scene.current_atom_symbol=parts[1]
            self.statusBar().showMessage(f"Mode: Draw Atom ({parts[1]})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        elif self.scene.mode=='bond': 
            self.scene.bond_order=int(parts[1])
            self.statusBar().showMessage(f"Mode: Draw Bond (Order: {parts[1]})")
            self.view_2d.setDragMode(QGraphicsView.DragMode.NoDrag)
        else: # Select mode
            self.scene.bond_order=1
            self.statusBar().showMessage("Mode: Select")
            self.view_2d.setDragMode(QGraphicsView.DragMode.RubberBandDrag)

    # --- 2D and 3D Conversion ---
    def clean_up_2d_structure(self):
        self.statusBar().showMessage("Optimizing 2D structure...")
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No atoms to optimize."); return
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if mol is None: self.statusBar().showMessage("Error: Failed to create molecule for optimization."); return
        try:
            AllChem.Compute2DCoords(mol); conf=mol.GetConformer(); original_ids=list(self.data.atoms.keys()); SCALE=50.0
            positions=[conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            if not positions: self.statusBar().showMessage("Optimization complete."); return
            cx=sum(p.x for p in positions)/len(positions); cy=sum(p.y for p in positions)/len(positions)
            for i in range(mol.GetNumAtoms()):
                item=self.data.atoms[original_ids[i]]['item']; new_pos=conf.GetAtomPosition(i)
                sx=(new_pos.x-cx)*SCALE; sy=-(new_pos.y-cy)*SCALE; item.setPos(sx,sy)
            for bond_data in self.data.bonds.values():
                if bond_data.get('item'): bond_data['item'].update_position()
            self.statusBar().showMessage("2D structure optimization successful.")
            self.push_undo_state()
        except Exception as e: self.statusBar().showMessage(f"Error during 2D optimization: {e}")

    def trigger_conversion(self):
        mol_block = self.data.to_mol_block()
        if not mol_block: self.statusBar().showMessage("Error: No atoms to convert."); return
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
        if mol is None: self.statusBar().showMessage("Error: Invalid chemical structure."); return
        if len(Chem.GetMolFrags(mol)) > 1: self.statusBar().showMessage("Error: 3D conversion not supported for multiple molecules."); return
        self.convert_button.setEnabled(False)
        self.statusBar().showMessage("Calculating 3D structure...")
        self.start_calculation.emit(mol_block)
        
    def on_calculation_finished(self, mol):
        self.current_mol=mol
        self.draw_molecule_3d(mol)
        self.statusBar().showMessage("3D conversion successful.")
        self.convert_button.setEnabled(True)
        self.push_undo_state()
        
    def on_calculation_error(self, error_message):
        self.statusBar().showMessage(f"Error: {error_message}")
        self.convert_button.setEnabled(True)
        
    def draw_molecule_3d(self, mol):
        self.plotter.clear(); conf = mol.GetConformer()
        pos=np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        sym=[a.GetSymbol() for a in mol.GetAtoms()]
        rad=np.array([VDW_RADII.get(s,0.4) for s in sym]); col=np.array([CPK_COLORS_PV.get(s,[0.5,0.5,0.5]) for s in sym])
        poly=pv.PolyData(pos); poly['colors']=col; poly['radii']=rad
        glyphs=poly.glyph(scale='radii',geom=pv.Sphere(radius=1.0),orient=False)
        
        # 変更点: 塗りつぶし表示に戻し、エッジを薄いグレーで半透明(opacity=0.3)にする
        self.plotter.add_mesh(glyphs, scalars='colors', rgb=True, smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3, line_width=0.1)
        
        for bond in mol.GetBonds():
            sp=np.array(conf.GetAtomPosition(bond.GetBeginAtomIdx())); ep=np.array(conf.GetAtomPosition(bond.GetEndAtomIdx()))
            bt=bond.GetBondType(); c=(sp+ep)/2; d=ep-sp; h=np.linalg.norm(d)
            if h==0: continue
            color=[0.5, 0.5, 0.5] # Gray
            if bt==Chem.rdchem.BondType.SINGLE or bt==Chem.rdchem.BondType.AROMATIC:
                cyl=pv.Cylinder(center=c,direction=d,radius=0.1,height=h)
                # 変更点: エッジを薄いグレーで半透明(opacity=0.3)にする
                self.plotter.add_mesh(cyl, color=color, smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
            else:
                v1=d/h; v_arb=np.array([0,0,1])
                if np.allclose(np.abs(np.dot(v1,v_arb)),1.0): v_arb=np.array([0,1,0])
                off_dir=np.cross(v1,v_arb); off_dir/=np.linalg.norm(off_dir)
                if bt==Chem.rdchem.BondType.DOUBLE:
                    r=0.09; s=0.15; c1=c+off_dir*(s/2); c2=c-off_dir*(s/2)
                    cyl1=pv.Cylinder(center=c1,direction=d,radius=r,height=h); cyl2=pv.Cylinder(center=c2,direction=d,radius=r,height=h)
                    # 変更点: エッジを薄いグレーで半透明(opacity=0.3)にする
                    self.plotter.add_mesh(cyl1,color=color,smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
                    self.plotter.add_mesh(cyl2,color=color,smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
                elif bt==Chem.rdchem.BondType.TRIPLE:
                    r=0.08; s=0.18; cc=pv.Cylinder(center=c,direction=d,radius=r,height=h)
                    c1=pv.Cylinder(center=c+off_dir*s,direction=d,radius=r,height=h); c2=pv.Cylinder(center=c-off_dir*s,direction=d,radius=r,height=h)
                    # 変更点: エッジを薄いグレーで半透明(opacity=0.3)にする
                    self.plotter.add_mesh(cc,color=color,smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
                    self.plotter.add_mesh(c1,color=color,smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
                    self.plotter.add_mesh(c2,color=color,smooth_shading=True, show_edges=True, edge_color='gray', edge_opacity=0.3)
        self.plotter.reset_camera()
        
    def closeEvent(self, event):
        self.thread.quit(); self.thread.wait(); super().closeEvent(event)

# --- Application Execution ---
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
    