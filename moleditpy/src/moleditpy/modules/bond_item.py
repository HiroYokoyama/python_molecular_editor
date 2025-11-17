from PyQt6.QtWidgets import QGraphicsItem, QGraphicsScene

from PyQt6.QtGui import (
    QPen, QBrush, QColor, QFont, QPolygonF,
    QPainterPath, QPainterPathStroker, QFontMetricsF
)

from PyQt6.QtCore import (
    Qt, QPointF, QRectF, QLineF
)

try:
    from .constants import (
        EZ_LABEL_BOX_SIZE, EZ_LABEL_TEXT_OUTLINE, EZ_LABEL_MARGIN,
        BOND_OFFSET, FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH, DESIRED_BOND_PIXEL_WIDTH,
    )
except Exception:
    from modules.constants import (
        EZ_LABEL_BOX_SIZE, EZ_LABEL_TEXT_OUTLINE, EZ_LABEL_MARGIN,
        BOND_OFFSET, FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD,
        HOVER_PEN_WIDTH, DESIRED_BOND_PIXEL_WIDTH,
    )

class BondItem(QGraphicsItem):

    def get_ez_label_rect(self):
        """E/Zラベルの描画範囲（シーン座標）を返す。ラベルが無い場合はNone。"""
        if self.order != 2 or self.stereo not in [3, 4]:
            return None
        line = self.get_line_in_local_coords()
        center = line.center()
        label_width = EZ_LABEL_BOX_SIZE
        label_height = EZ_LABEL_BOX_SIZE
        label_rect = QRectF(center.x() - label_width/2, center.y() - label_height/2, label_width, label_height)
        # シーン座標に変換
        return self.mapToScene(label_rect).boundingRect()
    def set_stereo(self, new_stereo):
        try:
            # ラベルを消す場合は、消す前のboundingRectをscene().invalidateで強制的に無効化
            if new_stereo == 0 and self.stereo in [3, 4] and self.scene():
                rect = self.mapToScene(self.boundingRect()).boundingRect()
                self.scene().invalidate(rect, QGraphicsScene.SceneLayer.BackgroundLayer | QGraphicsScene.SceneLayer.ForegroundLayer)
            
            self.prepareGeometryChange()
            self.stereo = new_stereo
            self.update()
            
            if self.scene() and self.scene().views():
                try:
                    self.scene().views()[0].viewport().update()
                except (IndexError, RuntimeError):
                    # Handle case where views are being destroyed
                    pass
                    
        except Exception as e:
            print(f"Error in BondItem.set_stereo: {e}")
            # Continue without crashing
            self.stereo = new_stereo

    def set_order(self, new_order):
        self.prepareGeometryChange()
        self.order = new_order
        self.update()
        if self.scene() and self.scene().views():
            self.scene().views()[0].viewport().update()
    def __init__(self, atom1_item, atom2_item, order=1, stereo=0):
        super().__init__()
        # Validate input parameters
        if atom1_item is None or atom2_item is None:
            raise ValueError("BondItem requires non-None atom items")
        self.atom1, self.atom2, self.order, self.stereo = atom1_item, atom2_item, order, stereo
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.pen = QPen(Qt.GlobalColor.black, 2)
        self.setZValue(0)
        self.update_position()
        self.setAcceptHoverEvents(True)
        self.hovered = False


    def get_line_in_local_coords(self):
        if self.atom1 is None or self.atom2 is None:
            return QLineF(0, 0, 0, 0)
        try:
            p2 = self.mapFromItem(self.atom2, 0, 0)
            return QLineF(QPointF(0, 0), p2)
        except (RuntimeError, TypeError):
            # Handle case where atoms are deleted from scene
            return QLineF(0, 0, 0, 0)

    def boundingRect(self):
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            line = QLineF(0, 0, 0, 0)
        bond_offset = globals().get('BOND_OFFSET', 2)
        extra = (getattr(self, 'order', 1) - 1) * bond_offset + 20
        rect = QRectF(line.p1(), line.p2()).normalized().adjusted(-extra, -extra, extra, extra)

        # E/Zラベルの描画範囲も考慮して拡張（QFontMetricsFで正確に）
        if self.order == 2 and self.stereo in [3, 4]:
            font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD)
            font.setItalic(True)
            text = "Z" if self.stereo == 3 else "E"
            fm = QFontMetricsF(font)
            text_rect = fm.boundingRect(text)
            outline = EZ_LABEL_TEXT_OUTLINE  # 輪郭の太さ分
            margin = EZ_LABEL_MARGIN   # 追加余白
            center = line.center()
            label_rect = QRectF(center.x() - text_rect.width()/2 - outline - margin,
                                center.y() - text_rect.height()/2 - outline - margin,
                                text_rect.width() + 2*outline + 2*margin,
                                text_rect.height() + 2*outline + 2*margin)
            rect = rect.united(label_rect)
        return rect

    def shape(self):
        path = QPainterPath()
        try:
            line = self.get_line_in_local_coords()
        except Exception:
            return path 
        if line.length() == 0:
            return path

        scene = self.scene()
        if not scene or not scene.views():
            return super().shape()

        view = scene.views()[0]
        scale = view.transform().m11()

        scene_width = DESIRED_BOND_PIXEL_WIDTH / scale

        stroker = QPainterPathStroker()
        stroker.setWidth(scene_width)
        stroker.setCapStyle(Qt.PenCapStyle.RoundCap)  
        stroker.setJoinStyle(Qt.PenJoinStyle.RoundJoin) 

        center_line_path = QPainterPath(line.p1())
        center_line_path.lineTo(line.p2())
        
        return stroker.createStroke(center_line_path)

    def paint(self, painter, option, widget):
        if self.atom1 is None or self.atom2 is None:
            return
        line = self.get_line_in_local_coords()
        if line.length() == 0: return

        # --- 1. 選択状態に応じてペンとブラシを準備 ---
        if self.isSelected():
            selection_color = QColor("blue")
            painter.setPen(QPen(selection_color, 3))
            painter.setBrush(QBrush(selection_color))
        else:
            # Allow bond color override from app settings (2D color)
            try:
                sc = self.scene()
                if sc is not None and hasattr(sc, 'window') and sc.window is not None:
                    bond_hex = sc.window.settings.get('bond_color', '#222222')
                    bond_color = QColor(bond_hex)
                    painter.setPen(QPen(bond_color, 2))
                else:
                    painter.setPen(self.pen)
            except Exception:
                painter.setPen(self.pen)
            painter.setBrush(QBrush(Qt.GlobalColor.black))

        # --- 立体化学 (Wedge/Dash) の描画 ---
        if self.order == 1 and self.stereo in [1, 2]:
            vec = line.unitVector()
            normal = vec.normalVector()
            p1 = line.p1() + vec.p2() * 5
            p2 = line.p2() - vec.p2() * 5

            if self.stereo == 1: # Wedge (くさび形)
                offset = QPointF(normal.dx(), normal.dy()) * 6.0
                poly = QPolygonF([p1, p2 + offset, p2 - offset])
                painter.drawPolygon(poly)
            
            elif self.stereo == 2: # Dash (破線)
                painter.save()
                if not self.isSelected():
                    pen = painter.pen()
                    pen.setWidthF(2.5) 
                    painter.setPen(pen)
                
                num_dashes = 8
                for i in range(num_dashes + 1):
                    t = i / num_dashes
                    start_pt = p1 * (1 - t) + p2 * t
                    width = 12.0 * t
                    offset = QPointF(normal.dx(), normal.dy()) * width / 2.0
                    painter.drawLine(start_pt - offset, start_pt + offset)
                painter.restore()
        
        # --- 通常の結合 (単/二重/三重) の描画 ---
        else:
            if self.order == 1:
                painter.drawLine(line)
            else:
                v = line.unitVector().normalVector()
                offset = QPointF(v.dx(), v.dy()) * BOND_OFFSET

                if self.order == 2:
                    # -------------------- ここから差し替え --------------------)
                    line1 = line.translated(offset)
                    line2 = line.translated(-offset)
                    painter.drawLine(line1)
                    painter.drawLine(line2)

                    # E/Z ラベルの描画処理
                    if self.stereo in [3, 4]:
                        painter.save() # 現在の描画設定を保存

                        # --- ラベルの設定 ---
                        font = QFont(FONT_FAMILY, FONT_SIZE_LARGE, FONT_WEIGHT_BOLD)
                        font.setItalic(True)
                        text_color = QColor("gray")
                        # 輪郭の色を背景色と同じにする（scene()がNoneのときは安全なフォールバックを使う）
                        outline_color = None
                        try:
                            sc = self.scene()
                            if sc is not None:
                                outline_color = sc.backgroundBrush().color()
                        except Exception:
                            outline_color = None
                        if outline_color is None:
                            # デフォルトでは白背景を想定して黒系の輪郭が見やすい
                            outline_color = QColor(255, 255, 255)

                        # --- 描画パスの作成 ---
                        text = "Z" if self.stereo == 3 else "E"
                        path = QPainterPath()
                        
                        # テキストが正確に中央に来るように位置を計算
                        fm = QFontMetricsF(font)
                        text_rect = fm.boundingRect(text)
                        text_rect.moveCenter(line.center())
                        path.addText(text_rect.topLeft(), font, text)

                        # --- 輪郭の描画 ---
                        stroker = QPainterPathStroker()
                        stroker.setWidth(EZ_LABEL_TEXT_OUTLINE) # 輪郭の太さ
                        outline_path = stroker.createStroke(path)
                        
                        painter.setBrush(outline_color)
                        painter.setPen(Qt.PenStyle.NoPen)
                        painter.drawPath(outline_path)

                        # --- 文字本体の描画 ---
                        painter.setBrush(text_color)
                        painter.setPen(text_color)
                        painter.drawPath(path)

                        painter.restore() # 描画設定を元に戻す

                elif self.order == 3:
                    painter.drawLine(line)
                    painter.drawLine(line.translated(offset))
                    painter.drawLine(line.translated(-offset))

        # --- 2. ホバー時のエフェクトを上から重ねて描画 ---
        if (not self.isSelected()) and getattr(self, 'hovered', False):
            try:
                # ホバー時のハイライトを太めの半透明な線で描画
                hover_pen = QPen(QColor(144, 238, 144, 180), HOVER_PEN_WIDTH) # LightGreen, 半透明
                hover_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                painter.setPen(hover_pen)
                painter.drawLine(line) 
            except Exception:
                pass



    def update_position(self):
        try:
            self.prepareGeometryChange()
            if self.atom1:
                self.setPos(self.atom1.pos())
            self.update()
        except Exception as e:
            print(f"Error updating bond position: {e}")
            # Continue without crashing


    def hoverEnterEvent(self, event):
        scene = self.scene()
        mode = getattr(scene, 'mode', '')
        self.hovered = True
        self.update()
        if self.scene():
            self.scene().set_hovered_item(self)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        if self.hovered:
            self.hovered = False
            self.update()
        if self.scene():
            self.scene().set_hovered_item(None)
        super().hoverLeaveEvent(event)
