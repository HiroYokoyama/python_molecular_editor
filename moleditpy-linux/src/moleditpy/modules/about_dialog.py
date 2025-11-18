from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton, QHBoxLayout
)
from PyQt6.QtGui import QPixmap, QPainter, QPen, QCursor
from PyQt6.QtCore import Qt
import os
try:
    from .constants import VERSION
except Exception:
    from modules.constants import VERSION

class AboutDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("About MoleditPy")
        self.setFixedSize(250, 300)
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Create a clickable image label
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        # Load the original icon image
        icon_path = os.path.join(os.path.dirname(__file__), 'assets', 'icon.png')
        if os.path.exists(icon_path):
            original_pixmap = QPixmap(icon_path)
            # Scale to 2x size (160x160)
            pixmap = original_pixmap.scaled(160, 160, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
        else:
            # Fallback: create a simple placeholder if icon.png not found
            pixmap = QPixmap(160, 160)
            pixmap.fill(Qt.GlobalColor.lightGray)
            painter = QPainter(pixmap)
            painter.setPen(QPen(Qt.GlobalColor.black, 2))
            painter.drawText(pixmap.rect(), Qt.AlignmentFlag.AlignCenter, "MoleditPy")
            painter.end()
        
        self.image_label.setPixmap(pixmap)
        try:
            self.image_label.setCursor(QCursor(Qt.CursorShape.ArrowCursor))
        except Exception:
            pass

        self.image_label.mousePressEvent = self.image_mouse_press_event
        
        layout.addWidget(self.image_label)
        
        # Add text information
        info_text = f"MoleditPy for Linux Ver. {VERSION}\nAuthor: Hiromichi Yokoyama\nLicense: Apache-2.0"
        info_label = QLabel(info_text)
        info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(info_label)
        
        # Add OK button
        ok_button = QPushButton("OK")
        ok_button.setFixedSize(80, 30)  # 小さいサイズに固定
        ok_button.clicked.connect(self.accept)
        
        # Center the button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(ok_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)
    
    def image_clicked(self, event):
        """Easter egg: Clear all and load bipyrimidine from SMILES"""
        # Clear the current scene
        self.main_window.clear_all()

        bipyrimidine_smiles = "C1=CN=C(N=C1)C2=NC=CC=N2"
        self.main_window.load_from_smiles(bipyrimidine_smiles)

        # Close the dialog
        self.accept()

    def image_mouse_press_event(self, event):
        """Handle mouse press on the image: trigger easter egg only for right-click."""
        try:
            if event.button() == Qt.MouseButton.RightButton:
                self.image_clicked(event)
            else:
                event.ignore()
        except Exception:
            try:
                event.ignore()
            except Exception:
                pass
