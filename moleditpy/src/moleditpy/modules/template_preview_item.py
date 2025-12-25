#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from PyQt6.QtWidgets import (
    QGraphicsItem
)

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QFont, QPolygonF
)


from PyQt6.QtCore import (
    Qt, QPointF, QRectF, QLineF
)

try:
    from .constants import CPK_COLORS
except Exception:
    from modules.constants import CPK_COLORS

class TemplatePreviewItem(QGraphicsItem):
    def __init__(self):
        super().__init__()
        self.setZValue(2)
        self.pen = QPen(QColor(80, 80, 80, 180), 2)
        self.polygon = QPolygonF()
        self.is_aromatic = False
        self.user_template_points = []
        self.user_template_bonds = []
        self.user_template_atoms = []
        self.is_user_template = False

    def set_geometry(self, points, is_aromatic=False):
        self.prepareGeometryChange()
        self.polygon = QPolygonF(points)
        self.is_aromatic = is_aromatic
        self.is_user_template = False
        self.update()
    
    def set_user_template_geometry(self, points, bonds_info, atoms_data):
        self.prepareGeometryChange()
        self.user_template_points = points
        self.user_template_bonds = bonds_info
        self.user_template_atoms = atoms_data
        self.is_user_template = True
        self.is_aromatic = False
        self.polygon = QPolygonF()
        self.update()

    def boundingRect(self):
        if self.is_user_template and self.user_template_points:
            # Calculate bounding rect for user template
            min_x = min(p.x() for p in self.user_template_points)
            max_x = max(p.x() for p in self.user_template_points)
            min_y = min(p.y() for p in self.user_template_points)
            max_y = max(p.y() for p in self.user_template_points)
            return QRectF(min_x - 20, min_y - 20, max_x - min_x + 40, max_y - min_y + 40)
        return self.polygon.boundingRect().adjusted(-5, -5, 5, 5)

    def paint(self, painter, option, widget):
        if self.is_user_template:
            self.paint_user_template(painter)
        else:
            self.paint_regular_template(painter)
    
    def paint_regular_template(self, painter):
        painter.setPen(self.pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        if not self.polygon.isEmpty():
            painter.drawPolygon(self.polygon)
            if self.is_aromatic:
                center = self.polygon.boundingRect().center()
                radius = QLineF(center, self.polygon.first()).length() * 0.6
                painter.drawEllipse(center, radius, radius)
    
    def paint_user_template(self, painter):
        if not self.user_template_points:
            return
        
        # Draw bonds first with better visibility
        # Draw bonds first with better visibility
        # Use gray (ghost) color for template preview to distinguish from real bonds
        bond_pen = QPen(QColor(80, 80, 80, 180), 2.5)
        painter.setPen(bond_pen)
        
        for bond_info in self.user_template_bonds:
            if len(bond_info) >= 3:
                atom1_idx, atom2_idx, order = bond_info[:3]
            else:
                atom1_idx, atom2_idx = bond_info[:2]
                order = 1
                
            if atom1_idx < len(self.user_template_points) and atom2_idx < len(self.user_template_points):
                pos1 = self.user_template_points[atom1_idx]
                pos2 = self.user_template_points[atom2_idx]
                
                if order == 2:
                    # Double bond - draw two parallel lines
                    line = QLineF(pos1, pos2)
                    normal = line.normalVector()
                    normal.setLength(4)
                    
                    line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                    line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                    
                    painter.drawLine(line1)
                    painter.drawLine(line2)
                elif order == 3:
                    # Triple bond - draw three parallel lines
                    line = QLineF(pos1, pos2)
                    normal = line.normalVector()
                    normal.setLength(6)
                    
                    painter.drawLine(line)
                    line1 = QLineF(pos1 + normal.p2() - normal.p1(), pos2 + normal.p2() - normal.p1())
                    line2 = QLineF(pos1 - normal.p2() + normal.p1(), pos2 - normal.p2() + normal.p1())
                    
                    painter.drawLine(line1)
                    painter.drawLine(line2)
                else:
                    # Single bond
                    painter.drawLine(QLineF(pos1, pos2))
        
        # Draw atoms - white ellipse background to hide bonds, then CPK colored text
        for i, pos in enumerate(self.user_template_points):
            if i < len(self.user_template_atoms):
                atom_data = self.user_template_atoms[i]
                symbol = atom_data.get('symbol', 'C')
                
                # Draw all non-carbon atoms including hydrogen with white background ellipse + CPK colored text
                if symbol != 'C':
                    # Get CPK color for text
                    color = CPK_COLORS.get(symbol, CPK_COLORS.get('DEFAULT', QColor('#FF1493')))
                    
                    # Draw white background ellipse to hide bonds
                    painter.setPen(QPen(Qt.GlobalColor.white, 0))  # No border
                    painter.setBrush(QBrush(Qt.GlobalColor.white))
                    painter.drawEllipse(int(pos.x() - 12), int(pos.y() - 8), 24, 16)
                    
                    # Draw CPK colored text on top
                    painter.setPen(QPen(color))
                    font = QFont("Arial", 12, QFont.Weight.Bold)  # Larger font
                    painter.setFont(font)
                    metrics = painter.fontMetrics()
                    text_rect = metrics.boundingRect(symbol)
                    text_pos = QPointF(pos.x() - text_rect.width()/2, pos.y() + text_rect.height()/3)
                    painter.drawText(text_pos, symbol)
