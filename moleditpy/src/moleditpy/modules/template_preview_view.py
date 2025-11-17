from PyQt6.QtWidgets import QGraphicsView

from PyQt6.QtCore import Qt, QTimer

class TemplatePreviewView(QGraphicsView):
    """テンプレートプレビュー用のカスタムビュークラス"""
    
    def __init__(self, scene):
        super().__init__(scene)
        self.original_scene_rect = None
        self.template_data = None  # Store template data for dynamic redrawing
        self.parent_dialog = None  # Reference to parent dialog for redraw access
    
    def set_template_data(self, template_data, parent_dialog):
        """テンプレートデータと親ダイアログの参照を設定"""
        self.template_data = template_data
        self.parent_dialog = parent_dialog
    
    def resizeEvent(self, event):
        """リサイズイベントを処理してプレビューを再フィット"""
        super().resizeEvent(event)
        if self.original_scene_rect and not self.original_scene_rect.isEmpty():
            # Delay the fitInView call to ensure proper widget sizing
            QTimer.singleShot(10, self.refit_view)
    
    def refit_view(self):
        """ビューを再フィット"""
        try:
            if self.original_scene_rect and not self.original_scene_rect.isEmpty():
                self.fitInView(self.original_scene_rect, Qt.AspectRatioMode.KeepAspectRatio)
        except Exception as e:
            print(f"Warning: Failed to refit template preview: {e}")
    
    def showEvent(self, event):
        """表示イベントを処理"""
        super().showEvent(event)
        # Ensure proper fitting when widget becomes visible
        if self.original_scene_rect:
            QTimer.singleShot(50, self.refit_view)
    
    def redraw_with_current_size(self):
        """現在のサイズに合わせてテンプレートを再描画"""
        if self.template_data and self.parent_dialog:
            try:
                # Clear current scene
                self.scene().clear()
                
                # Redraw with current view size for proper fit-based scaling
                view_size = (self.width(), self.height())
                self.parent_dialog.draw_template_preview(self.scene(), self.template_data, view_size)
                
                # Refit the view
                bounding_rect = self.scene().itemsBoundingRect()
                if not bounding_rect.isEmpty() and bounding_rect.width() > 0 and bounding_rect.height() > 0:
                    content_size = max(bounding_rect.width(), bounding_rect.height())
                    padding = max(20, content_size * 0.2)
                    padded_rect = bounding_rect.adjusted(-padding, -padding, padding, padding)
                    self.scene().setSceneRect(padded_rect)
                    self.original_scene_rect = padded_rect
                    QTimer.singleShot(10, lambda: self.fitInView(padded_rect, Qt.AspectRatioMode.KeepAspectRatio))
            except Exception as e:
                print(f"Warning: Failed to redraw template preview: {e}")
